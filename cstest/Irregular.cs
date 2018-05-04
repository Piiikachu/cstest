using System;
using MPI_Request = System.Int32;
using bigint = System.Int64;
using System.IO;
using System.Text;

namespace cstest
{
    public class Irregular
    {
        public static int[] proc_recv_copy;

        public void create_procs(int n, int[] proclist, int sort = 0)
        {
            int i, m;

            // setup for collective comm
            // work1 = 1 for procs I send to, set self to 0
            // work2 = 1 for all procs, used for ReduceScatter
            // nsend = # of procs I send messages to, not including self

            for (i = 0; i < nprocs; i++)
            {
                work1[i] = 0;
                work2[i] = 1;
            }
            for (i = 0; i < n; i++) work1[proclist[i]] = 1;

            nsend = n;
            if (work1[me]!=0)
            {
                work1[me] = 0;
                nsend--;
            }
            // nrecv = # of procs I receive messages from, not including self
            // options for performing ReduceScatter operation
            // some are more efficient on some machines at big sizes

            //sparta.mpi.MPI_Reduce_scatter(work1, &nrecv, work2, MPI_INT, MPI_SUM, world);

            // proc_send = procs I send to
            // to balance pattern of send messages:
            //   each proc starts with iproc > me, continues until iproc = me
            // reset work1 to store which send message each proc corresponds to
            //   used by augmen_data()
            for (i = 0; i < nprocs; i++) work1[i] = 0;
            for (i = 0; i < n; i++) work1[proclist[i]] = 1;
            work1[me] = 0;

            int iproc = me;
            int isend = 0;
            for (i = 0; i < nprocs; i++)
            {
                iproc++;
                if (iproc == nprocs) iproc = 0;
                if (iproc == me) continue;
                if (work1[iproc]!=0)
                {
                    proc_send[isend] = iproc;
                    work1[iproc] = isend;
                    isend++;
                }
            }

            // tell receivers I send to them

            m = 0;
            for (i = 0; i < nsend; i++)
                sparta.mpi.MPI_Send(ref m, 0, MPI.MPI_INT, proc_send[i], 0, sparta.world);

            // receive incoming messages
            // proc_recv = procs I recv from

            for (i = 0; i < nrecv; i++)
            {
                sparta.mpi.MPI_Recv(ref m, 0, MPI.MPI_INT, MPI.MPI_ANY_SOURCE, 0, sparta.world,ref status[0]);
                proc_recv[i] = status[0].MPI_SOURCE;
            }

            // sort proc_recv by proc ID if requested
            // useful for debugging to insure reproducible ordering of received datums

            if (sort!=0)
            {
                int[] order = new int[nrecv];
                int[] proc_recv_ordered = new int[nrecv];

                for (i = 0; i < nrecv; i++) order[i] = i;
                proc_recv_copy = proc_recv;
                //qsort(order, nrecv, sizeof(int), compare_standalone);
                Array.Sort(order);
                int j;
                for (i = 0; i < nrecv; i++)
                {
                    j = order[i];
                    proc_recv_ordered[i] = proc_recv[j];
                }

                //memcpy(proc_recv, proc_recv_ordered, nrecv * sizeof(int));
                //delete[] order;
                //delete[] proc_recv_ordered;
            }

            // proc2recv[I] = which recv the Ith proc ID is
            // will only be accessed by procs I actually receive from

            for (i = 0; i < nrecv; i++) proc2recv[proc_recv[i]] = i;

            // barrier to insure all MPI_ANY_SOURCE messages are received
            // else another proc could proceed to augment_data() and send to me

            sparta.mpi.MPI_Barrier(sparta.world);

        }
        public int create_data_uniform(int n, int[] proclist, int sort = 0)
        {
            int i, m;

            // setup for collective comm
            // work1 = # of datums I send to each proc, set self to 0
            // work2 = 1 for all procs, used for ReduceScatter

            for (i = 0; i < nprocs; i++)
            {
                work1[i] = 0;
                work2[i] = 1;
            }
            for (i = 0; i < n; i++) work1[proclist[i]] = 1;
            work1[me] = 0;

            // nrecv = # of procs I receive messages from, not including self
            // options for performing ReduceScatter operation
            // some are more efficient on some machines at big sizes

            // work1 = # of datums I send to each proc, including self
            // nsend = # of procs I send messages to, not including self

            for (i = 0; i < nprocs; i++) work1[i] = 0;
            for (i = 0; i < n; i++) work1[proclist[i]]++;

            nsend = 0;
            for (i = 0; i < nprocs; i++)
                if (work1[i]!=0) nsend++;
            if (work1[me]!=0) nsend--;

            // reallocate send and self index lists if necessary
            // could use n-work1[me] for length of index_send to be more precise

            if (n > indexmax)
            {
                indexmax = n;
                index_send = new int[indexmax];
            }

            if (work1[me] > indexselfmax)
            {
                indexselfmax = work1[me];
                index_self = new int[indexselfmax];
            }
            // proc_send = procs I send to
            // num_send = # of datums I send to each proc
            // num_self = # of datums I copy to self
            // to balance pattern of send messages:
            //   each proc starts with iproc > me, continues until iproc = me
            // reset work1 to store which send message each proc corresponds to

            int iproc = me;
            int isend = 0;
            for (i = 0; i < nprocs; i++)
            {
                iproc++;
                if (iproc == nprocs) iproc = 0;
                if (iproc == me)
                {
                    num_self = work1[iproc];
                    work1[iproc] = 0;
                }
                else if (work1[iproc]!=0)
                {
                    proc_send[isend] = iproc;
                    num_send[isend] = work1[iproc];
                    work1[iproc] = isend;
                    isend++;
                }
            }
            // work2 = offsets into index_send for each proc I send to
            // m = ptr into index_self
            // index_send = list of which datums to send to each proc
            //   1st N1 values are datum indices for 1st proc,
            //   next N2 values are datum indices for 2nd proc, etc
            // index_self = list of which datums to copy to self

            work2[0] = 0;
            for (i = 1; i < nsend; i++) work2[i] = work2[i - 1] + num_send[i - 1];

            m = 0;
            for (i = 0; i < n; i++)
            {
                iproc = proclist[i];
                if (iproc == me) index_self[m++] = i;
                else
                {
                    isend = work1[iproc];
                    index_send[work2[isend]++] = i;
                }
            }
            // tell receivers how many datums I send them
            // sendmax = largest # of datums I send in a single message

            sendmax = 0;
            for (i = 0; i < nsend; i++)
            {
                sparta.mpi.MPI_Send(ref num_send[i], 1, MPI.MPI_INT, proc_send[i], 0, sparta.world);
                sendmax =Math.Max(sendmax, num_send[i]);
            }
            // receive incoming messages
            // proc_recv = procs I recv from
            // num_recv = # of datums each proc sends me
            // nrecvdatum = total # of datums I recv

            nrecvdatum = 0;
            for (i = 0; i < nrecv; i++)
            {
                sparta.mpi.MPI_Recv(ref num_recv[i], 1, MPI.MPI_INT, MPI.MPI_ANY_SOURCE, 0, sparta.world,ref status[0]);
                proc_recv[i] = status[0].MPI_SOURCE;
                nrecvdatum += num_recv[i];
            }
            nrecvdatum += num_self;
            // sort proc_recv and num_recv by proc ID if requested
            // useful for debugging to insure reproducible ordering of received datums

            if (sort!=0)
            {
                int[] order = new int[nrecv];
                int[] proc_recv_ordered = new int[nrecv];
                int[] num_recv_ordered = new int[nrecv];

                for (i = 0; i < nrecv; i++) order[i] = i;
                proc_recv_copy = proc_recv;
                //qsort(order, nrecv, sizeof(int), compare_standalone);
                Array.Sort(order);
                int j;
                for (i = 0; i < nrecv; i++)
                {
                    j = order[i];
                    proc_recv_ordered[i] = proc_recv[j];
                    num_recv_ordered[i] = num_recv[j];
                }

                //memcpy(proc_recv, proc_recv_ordered, nrecv * sizeof(int));
                //memcpy(num_recv, num_recv_ordered, nrecv * sizeof(int));
                proc_recv = proc_recv_ordered;
                num_recv = num_recv_ordered;
            }
            // proc2recv[I] = which recv the Ith proc ID is
            // will only be accessed by procs I actually receive from

            for (i = 0; i < nrecv; i++) proc2recv[proc_recv[i]] = i;

            // barrier to insure all MPI_ANY_SOURCE messages are received
            // else another proc could proceed to exchange_data() and send to me

            sparta.mpi.MPI_Barrier(sparta.world);

            // return # of datums I will receive

            return nrecvdatum;
        }
        public int create_data_variable(int n, int[] proclist, int[] sizes,
                                    out int recvsize, int sort = 0)
        {
            int i;
            int nrecvdatum = create_data_uniform(n, proclist, sort);

            if (size_send==null)
            {
                size_send = new int[nprocs];
                size_recv = new int[nprocs];
            }

            if (n> offsetmax)
            {
                offsetmax = n;
                offset_send = new int[offsetmax];
            }

            int offset = 0;
            for (i = 0; i < n; i++)
            {
                offset_send[i] = offset;
                offset += sizes[i];
            }
            // work1 = # of bytes to send to each proc, including self

            for (i = 0; i < nprocs; i++) work1[i] = 0;
            for (i = 0; i < n; i++) work1[proclist[i]] += sizes[i];

            // size_send = # of bytes I send to each proc
            // size_self = # of bytes I copy to self
            // to balance pattern of send messages:
            //   each proc starts with iproc > me, continues until iproc = me

            int iproc = me;
            int isend = 0;
            for (i = 0; i < nprocs; i++)
            {
                iproc++;
                if (iproc == nprocs) iproc = 0;
                if (iproc == me) size_self = work1[iproc];
                else if (work1[iproc] > 0) size_send[isend++] = work1[iproc];
            }

            // tell receivers how many bytes I send them
            // sendmaxbytes = largest # of bytes I send in a single message

            sendmaxbytes = 0;
            for (i = 0; i < nsend; i++)
            {
                sparta.mpi.MPI_Send(ref size_send[i], 1, MPI.MPI_INT, proc_send[i], 1, sparta.world);
                sendmaxbytes = Math.Max(sendmaxbytes, size_send[i]);
            }


            // receive incoming messages
            // num_recv = # of datums each proc sends me
            // nrecvdatum = total # of datums I recv

            int nbytes=0;
            bigint brecvsize = 0;
            for (i = 0; i < nrecv; i++)
            {
                sparta.mpi.MPI_Recv(ref nbytes, 1, MPI.MPI_INT, MPI.MPI_ANY_SOURCE, 1, sparta.world,ref status[0]);
                size_recv[proc2recv[status[0].MPI_SOURCE]] = nbytes;
                brecvsize += nbytes;
            }
            brecvsize += size_self;
            
            if (brecvsize > Run.MAXSMALLINT) 
                sparta.error.one("Irregular comm recv buffer exceeds 2 GB");
            recvsize = (int)brecvsize;

            // return # of datums I will receive

            return nrecvdatum;

        }
        public int augment_data_uniform(int n, int[] proclist)
        {
            int i, m=0, iproc, isend;

            // tally count of messages to each proc in num_send and num_self

            num_self = 0;
            for (i = 0; i < nsend; i++) work2[proc_send[i]] = 0;
            work2[me] = 0;
            for (i = 0; i < n; i++) work2[proclist[i]]++;
            for (i = 0; i < nsend; i++) num_send[i] = work2[proc_send[i]];
            num_self = work2[me];

            // reallocate send and self index lists if necessary
            // could use n-num_self for length of index_send to be more precise

            if (n > indexmax)
            {
                indexmax = n;
                //memory->destroy(index_send);
                //memory->create(index_send, indexmax, "irregular:index_send");
                index_send = new int[indexmax];
            }

            if (num_self > indexselfmax)
            {
                indexselfmax = num_self;
                //memory->destroy(index_self);
                //memory->create(index_self, indexselfmax, "irregular:index_self");
                index_self = new int[indexselfmax];
            }

            // work2 = offsets into index_send for each proc I send to
            // m = ptr into index_self
            // index_send = list of which datums to send to each proc
            //   1st N1 values are datum indices for 1st proc,
            //   next N2 values are datum indices for 2nd proc, etc
            // index_self = list of which datums to copy to self

            work2[0] = 0;
            for (i = 1; i < nsend; i++) work2[i] = work2[i - 1] + num_send[i - 1];

            if (num_self!=0)
            {
                m = 0;
                for (i = 0; i < n; i++)
                {
                    iproc = proclist[i];
                    if (iproc == me) index_self[m++] = i;
                    else
                    {
                        isend = work1[iproc];
                        index_send[work2[isend]++] = i;
                    }
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    isend = work1[proclist[i]];
                    index_send[work2[isend]++] = i;
                }
            }

            // tell receivers how many datums I send them
            // sendmax = largest # of datums I send in a single message

            sendmax = 0;
            for (i = 0; i < nsend; i++)
            {
                sparta.mpi.MPI_Send(ref num_send[i], 1, MPI.MPI_INT, proc_send[i], 0, sparta.world);
                sendmax = Math.Max(sendmax, num_send[i]);
            }

            // receive incoming messages
            // num_recv = # of datums each proc sends me
            // nrecvdatum = total # of datums I recv

            nrecvdatum = 0;
            for (i = 0; i < nrecv; i++)
            {
                sparta.mpi.MPI_Recv(ref m, 1, MPI.MPI_INT, MPI.MPI_ANY_SOURCE, 0, sparta.world,ref status[0]);
                iproc = status[0].MPI_SOURCE;
                num_recv[proc2recv[iproc]] = m;
                nrecvdatum += m;
            }
            nrecvdatum += num_self;

            // barrier to insure all MPI_ANY_SOURCE messages are received
            // else another proc could proceed to exchange_data() and send to me

            sparta.mpi.MPI_Barrier(sparta.world);

            // return # of datums I will receive

            return nrecvdatum;
        }
        public void exchange_uniform(string sendbuf, int nbytes,string recvbuf)
        {
            int i, m, n, offset, count;

            // post all receives, starting after self copies

            offset = num_self * nbytes;
            for (int irecv = 0; irecv < nrecv; irecv++)
            {

                sparta.mpi.MPI_Irecv(ref recvbuf, num_recv[irecv] * nbytes, MPI.MPI_CHAR,
                      proc_recv[irecv], 0, sparta.world,ref  request[irecv]);
                offset += num_recv[irecv] * nbytes;
            }

            // reallocate buf for largest send if necessary

            if (sendmax * nbytes > bufmax)
            {
                bufmax = sendmax * nbytes;
                
            }

            // send each message
            // pack buf with list of datums
            // m = index of datum in sendbuf

            n = 0;
            for (int isend = 0; isend < nsend; isend++)
            {
                count = num_send[isend];
                for (i = 0; i < count; i++)
                {
                    m = index_send[n++];
                    //memcpy(&buf[i * nbytes], &sendbuf[m * nbytes], nbytes);
                    buf = Encoding.UTF8.GetBytes( sendbuf);


                }
                sparta.mpi.MPI_Send(ref buf, count * nbytes, MPI.MPI_CHAR, proc_send[isend], 0, sparta.world);
            }

            // copy datums to self, put at beginning of recvbuf

            for (i = 0; i < num_self; i++)
            {
                m = index_self[i];
                //memcpy(&recvbuf[i * nbytes], &sendbuf[m * nbytes], nbytes);
                recvbuf = sendbuf;
            }

            // wait on all incoming messages

            if (nrecv != 0) Console.WriteLine("irregular.exchange_uniform()->MPI_waitall");//sparta.mpi.MPI_Waitall(nrecv, request, status);
        }
        public void exchange_variable(byte[] sendbuf, int[] nbytes, byte[] recvbuf)
        {
            int i, m, n, offset, count;

            // post all receives, starting after self copies

            offset = size_self;
            for (int irecv = 0; irecv < nrecv; irecv++)
            {
                sparta.mpi.MPI_Irecv(ref recvbuf[offset], size_recv[irecv], MPI.MPI_CHAR,
                      proc_recv[irecv], 0, sparta.world, ref request[irecv]);
                offset += size_recv[irecv];
            }
            // reallocate buf for largest send if necessary

            if (sendmaxbytes > bufmax)
            {
                //memory->destroy(buf);
                bufmax = sendmaxbytes;
                buf = new byte[bufmax];
                //memory->create(buf, bufmax, "irregular:buf");
            }
            // send each message
            // pack buf with list of datums
            // m = index of datum in sendbuf
            // offset_send[m] = starting loc of datum in sendbuf

            n = 0;
            for (int isend = 0; isend < nsend; isend++)
            {
                offset = 0;
                count = num_send[isend];
                for (i = 0; i < count; i++)
                {
                    m = index_send[n++];
                    //memcpy(&buf[offset], &sendbuf[offset_send[m]], nbytes[m]);
                    Array.Copy(sendbuf, offset_send[m], buf, offset, buf.Length - offset);
                    offset += nbytes[m];
                }
                sparta.mpi.MPI_Send(ref buf, size_send[isend], MPI.MPI_CHAR, proc_send[isend], 0, sparta.world);
            }

            // copy datums to self, put at beginning of recvbuf

            offset = 0;
            for (i = 0; i < num_self; i++)
            {
                m = index_self[i];
                //memcpy(&recvbuf[offset], &sendbuf[offset_send[m]], nbytes[m]);
                Array.Copy(sendbuf, offset_send[m], buf, offset, buf.Length - offset);
                offset += nbytes[m];
            }

            // wait on all incoming messages

            if (nrecv != 0) Console.WriteLine("irregular.exchange_variable()->MPI_Waitall");//sparta.mpi.MPI_Waitall(nrecv, request, status);
        }
        //public void reverse(int, int*);


        protected int me, nprocs;

        // plan for irregular communication of datums
        // same for uniform or variable sized datums

        protected int nsend;                 // # of messages to send, no self
        protected int nrecv;                 // # of messages to recv, no self
        protected int sendmax;               // # of datums in largest send message
        protected int nrecvdatum;            // total # of datums I recv
        protected int num_self;              // # of datums to copy to self
        protected int indexmax;              // current size of index_send
        protected int indexselfmax;          // current size of index_self
        protected int bufmax;                // current size of buf in bytes
        protected int[] proc_send;            // list of procs to send to
        protected int[] num_send;             // # of datums to send to each proc
        protected int[] proc_recv;            // list of procs to recv from
        protected int[] num_recv;             // # of datums to recv from each proc
        protected int[] index_send;           // list of datum indices to send to each proc
        protected int[] index_self;           // list of datum indices to copy to self
        protected int[] proc2recv;            // mapping from proc IDs to recv list
        protected int[] work1,work2;         // work vectors
        protected MPI_Request[] request;      // MPI requests for posted recvs
        protected MPI._MPI_Status[] status;        // MPI statuses for WaitAll
        protected byte[] buf;                 // buffer for largest single send message
        protected int copymode;              // 1 if copy of class (prevents deallocation of
                                   //   base class when child copy is destroyed)

        // only defined for variable sized datums

        protected int sendmaxbytes;          // # of bytes in largest send message
        protected int size_self;             // # of bytes in datums copied to self
        protected int offsetmax;             // current size of offset_send
        protected int[] size_send;            // # of bytes of send to each proc
        protected int[] size_recv;            // # of bytes to recv from each proc
        protected int[] offset_send;          // list of byte offsets for each send datum
        private SPARTA sparta;
        public Irregular(SPARTA sparta)
        {
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);
            sparta.mpi.MPI_Comm_size(sparta.world, ref nprocs);

            // allocate fixed-length and work vectors for plan
            proc_send = new int[nprocs];
            num_send = new int[nprocs];
            proc_recv = new int[nprocs];
            num_recv = new int[nprocs];
            proc2recv = new int[nprocs];





            request = new MPI_Request[nprocs];
            status = new MPI._MPI_Status[nprocs];

            size_send = null;
            size_recv = null;

            work1 = new int[nprocs];
            work2 = new int[nprocs];

            indexmax = 0;
            index_send = null;
            indexselfmax = 0;
            index_self = null;
            offsetmax = 0;
            offset_send = null;
            bufmax = 0;
            buf = null;

            copymode = 0;
        }
    }
}