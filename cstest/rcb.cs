using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MPI_Op = System.Int32;
using MPI_Comm = System.Int32;
using MPI_Request = System.Int32;
using MPI_Datatype = System.Int32;
using System.Runtime.InteropServices;

namespace cstest
{
    class RCB
    {
        public const double MYHUGE = 1.0e30;
        public const double TINY = 1.0e-6;

        // set this to bigger number after debugging

        public const int DELTA = 0;
        int me, nprocs;

        // point to balance on

        public struct Dot
        {
            public double[] x;          // coord of point
            public double wt;            // weight of point
            public int proc;             // owning proc
            public int index;            // index on owning proc
        }

        // tree of RCB cuts

        public struct Tree
        {
            public double cut;         // position of cut
            public int dim;            // dimension = 0/1/2 of cut
        }

        // inversion message

        public struct Invert
        {
            public int rindex;         // index on receiving proc
            public int sproc;          // sending proc
            public int sindex;         // index on sending proc
        };

        Dot[] dots;        // dots on this proc
        int ndot;         // # of dots on this proc
        int maxdot;       // allocated size of dots
        int ndotorig;

        int nlist;
        int maxlist;
        int[] dotlist;
        int[] dotmark;

        int maxbuf;
        Dot[] buf;

        int maxrecv, maxsend;

        Irregular irregular;

        MPI_Op box_op, med_op;
        MPI_Datatype box_type, med_type;

        int reuse;        // 1/0 to use/not use previous cuts
        int dottop;       // dots >= this index are new
        BBox rcbbox;      // bounding box of final RCB sub-domain
        Tree[] tree;       // tree of RCB cuts, used by reuse()
        int[] counters=new int[7];  // diagnostic counts
                          // 0 = # of median iterations
                          // 1 = # of points sent
                          // 2 = # of points received
                          // 3 = most points this proc ever owns
                          // 4 = most point memory this proc ever allocs
                          // 5 = # of times a previous cut is re-used
                          // 6 = # of reallocs of point vector

        // set by compute()

        public int noriginal;              // # of dots I own before balancing
        public int nfinal;                 // # of dots I own after balancing
        public int nkeep;                  // how many dots of noriginal I still own
                                     // will be first nkept of nfinal list
        public int[] recvproc;              // proc IDs of nfinal dots
        public int[] recvindex;             // index of nfinal dots on owning procs
                                     // based on input list for compute()
        public double[] lo,hi;             // final bounding box of my RCB sub-domain
         
             // set by invert()
         
        public int[] sendproc;              // proc to send each of my noriginal dots to
        public int[] sendindex;             // index of dot in receiver's nfinal list

        private SPARTA sparta;
        public RCB(SPARTA sparta)
        {
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);
            sparta.mpi.MPI_Comm_size(sparta.world, ref nprocs);

            ndot = maxdot = 0;
            dots = null;

            nlist = maxlist = 0;
            dotlist = dotmark = null;

            maxbuf = 0;
            buf = null;

            maxrecv = maxsend = 0;
            recvproc = recvindex = sendproc = sendindex = null;

            //tree = (Tree*)memory->smalloc(nprocs * sizeof(Tree), "RCB:tree");
            tree = new Tree[nprocs];
            irregular = null;

            // create MPI data and function types for box and median AllReduce ops

            sparta.mpi.MPI_Type_contiguous(6, MPI.MPI_DOUBLE, ref box_type);
            sparta.mpi.MPI_Type_commit(ref box_type);
            sparta.mpi.MPI_Type_contiguous(Marshal.SizeOf(typeof(Median)), MPI.MPI_CHAR, ref med_type);
            sparta.mpi.MPI_Type_commit(ref med_type);

            sparta.mpi.MPI_Op_create(box_merge, 1, ref box_op);
            sparta.mpi.MPI_Op_create(median_merge, 1, ref med_op);

            reuse = 0;
        }
        public void compute(int n, double[,] x, double[] wt, int flip = 0)
        {
            int i, j, k;
            int keep, outgoing, incoming, incoming2=0;
            int dim;
            bool markactive;
            int indexlo=0, indexhi=0;
            int first_iteration, breakflag;
            double wttot=0, wtlo, wthi, wtsum, wtok, wtupto=0, wtmax;
            double targetlo, targethi;
            double valuemin, valuemax, valuehalf;
            double tolerance=1;
            MPI_Comm comm=1, comm_half=0;
            MPI_Request request=0, request2=0;
            MPI._MPI_Status status=new MPI._MPI_Status();
            Median med, medme;
            med = new Median();
            // create list of my Dots

            ndot = nkeep = noriginal = n;
            if (ndot > maxdot)
            {
                maxdot = ndot;
                dots = new Dot[ndot];

            }
            if (flip == 0)
            {
                for (i = 0; i < ndot; i++)
                {
                    dots[i].x[0] = x[i, 0];
                    dots[i].x[1] = x[i, 1];
                    dots[i].x[2] = x[i, 2];
                    dots[i].proc = me;
                    dots[i].index = i;
                }
            }
            else
            {
                for (i = 0; i < ndot; i++)
                {
                    dots[i].x[0] = -x[i, 0];
                    dots[i].x[1] = -x[i, 1];
                    dots[i].x[2] = -x[i, 2];
                    dots[i].proc = me;
                    dots[i].index = i;
                }
            }

            if (wt != null)
                for (i = 0; i < ndot; i++) dots[i].wt = wt[i];
            else
                for (i = 0; i < ndot; i++) dots[i].wt = 1.0;

            // shrink-wrap initial bounding box around dots

            BBox boxtmp;
            boxtmp.lo = new double[3];
            boxtmp.hi = new double[3];
            boxtmp.lo[0] = boxtmp.lo[1] = boxtmp.lo[2] = MYHUGE;
            boxtmp.hi[0] = boxtmp.hi[1] = boxtmp.hi[2] = -MYHUGE;

            for (i = 0; i < ndot; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    if (dots[i].x[j] < boxtmp.lo[j])
                        boxtmp.lo[j] = dots[i].x[j];
                    if (dots[i].x[j] > boxtmp.hi[j])
                        boxtmp.hi[j] = dots[i].x[j];
                }
            }
            sparta.mpi.MPI_Allreduce(ref boxtmp, ref rcbbox, 1, box_type, box_op, sparta.world);
            // initialize counters

            counters[0] = 0;
            counters[1] = 0;
            counters[2] = 0;
            counters[3] = ndot;
            counters[4] = maxdot;
            counters[5] = 0;
            counters[6] = 0;

            // create communicator for use in recursion

            sparta.mpi.MPI_Comm_dup(sparta.world, ref comm);

            // recurse until partition is a single proc = me
            // proclower,procupper = lower,upper procs in partition
            // procmid = 1st proc in upper half of partition

            int procpartner, procpartner2=0;

            int procmid;
            int proclower = 0;
            int procupper = nprocs - 1;

            while (proclower != procupper)
            {
                // if odd # of procs, lower partition gets extra one

                procmid = proclower + (procupper - proclower) / 2 + 1;

                // determine communication partner(s)
                // readnumber = # of proc partners to read from

                if (me < procmid)
                    procpartner = me + (procmid - proclower);
                else
                    procpartner = me - (procmid - proclower);

                int readnumber = 1;
                if (procpartner > procupper)
                {
                    readnumber = 0;
                    procpartner--;
                }
                if (me == procupper && procpartner != procmid - 1)
                {
                    readnumber = 2;
                    procpartner2 = procpartner + 1;
                }

                // wttot = summed weight of entire partition
                // search tolerance = largest single weight (plus epsilon
                // targetlo = desired weight in lower half of partition
                // targethi = desired weight in upper half of partition

                wtmax = wtsum = 0.0;
                for (i = 0; i < ndot; i++)
                {
                    wtsum += dots[i].wt;
                    if (dots[i].wt > wtmax) wtmax = dots[i].wt;
                }

                sparta.mpi.MPI_Allreduce(ref wtsum, ref wttot, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, comm);
                sparta.mpi.MPI_Allreduce(ref wtmax, ref tolerance, 1, MPI.MPI_DOUBLE, MPI.MPI_MAX, comm);

                tolerance *= 1.0 + TINY;
                targetlo = wttot * (procmid - proclower) / (procupper + 1 - proclower);
                targethi = wttot - targetlo;

                // dim = dimension to bisect on
                dim = 0;
                if (rcbbox.hi[1] - rcbbox.lo[1] > rcbbox.hi[0] - rcbbox.lo[0])
                    dim = 1;
                if (dim == 0 && rcbbox.hi[2] - rcbbox.lo[2] >
                    rcbbox.hi[0] - rcbbox.lo[0])
                    dim = 2;
                if (dim == 1 && rcbbox.hi[2] - rcbbox.lo[2] >
                rcbbox.hi[1] - rcbbox.lo[1])
                    dim = 2;

                // create active list and mark array for dots
                // initialize active list to all dots
                if (ndot > maxlist)
                {


                    maxlist = maxdot;
                    dotlist = new int[maxlist];
                    dotmark = new int[maxlist];

                }

                nlist = ndot;
                for (i = 0; i < nlist; i++) dotlist[i] = i;
                // median iteration
                // zoom in on bisector until correct # of dots in each half of partition
                // as each iteration of median-loop begins, require:
                //   all non-active dots are marked with 0/1 in dotmark
                //   valuemin <= every active dot <= valuemax
                //   wtlo, wthi = total wt of non-active dots
                // when leave median-loop, require only:
                //   valuehalf = correct cut position
                //   all dots <= valuehalf are marked with 0 in dotmark
                //   all dots >= valuehalf are marked with 1 in dotmark
                // markactive = which side of cut is active = 0/1
                // indexlo,indexhi = indices of dot closest to median

                wtlo = wthi = 0.0;
                valuemin = rcbbox.lo[dim];
                valuemax = rcbbox.hi[dim];
                first_iteration = 1;

                while (true)
                {
                    // choose bisector value
                    // use old value on 1st iteration if old cut dimension is the same
                    // on 2nd option: could push valuehalf towards geometric center 
                    //   with "1.0-factor" to force overshoot
                    if (first_iteration != 0 && reuse != 0 && dim == tree[procmid].dim)
                    {
                        counters[5]++;
                        valuehalf = tree[procmid].cut;
                        if (valuehalf < valuemin || valuehalf > valuemax)
                            valuehalf = 0.5 * (valuemin + valuemax);
                    }
                    else if (wt != null)
                    {
                        valuehalf = valuemin + (targetlo - wtlo) /
                          (wttot - wtlo - wthi) * (valuemax - valuemin);
                    }
                    else
                    {
                        valuehalf = 0.5 * (valuemin + valuemax);
                    }

                    first_iteration = 0;

                    // initialize local median data structure

                    medme.totallo = medme.totalhi = 0.0;
                    medme.valuelo = -MYHUGE;
                    medme.valuehi = MYHUGE;
                    medme.wtlo = medme.wthi = 0.0;
                    medme.countlo = medme.counthi = 0;
                    medme.proclo = medme.prochi = me;

                    // mark all active dots on one side or other of bisector
                    // also set all fields in median data struct
                    // save indices of closest dots on either side

                    for (j = 0; j < nlist; j++)
                    {
                        i = dotlist[j];
                        if (dots[i].x[dim] <= valuehalf)
                        {            // in lower part
                            medme.totallo += dots[i].wt;
                            dotmark[i] = 0;
                            if (dots[i].x[dim] > medme.valuelo)
                            {       // my closest dot
                                medme.valuelo = dots[i].x[dim];
                                medme.wtlo = dots[i].wt;
                                medme.countlo = 1;
                                indexlo = i;
                            }
                            else if (dots[i].x[dim] == medme.valuelo)
                            {   // tied for closest
                                medme.wtlo += dots[i].wt;
                                medme.countlo++;
                            }
                        }
                        else
                        {
                            medme.totalhi += dots[i].wt;
                            dotmark[i] = 1;
                            if (dots[i].x[dim] < medme.valuehi)
                            {       // my closest dot
                                medme.valuehi = dots[i].x[dim];
                                medme.wthi = dots[i].wt;
                                medme.counthi = 1;
                                indexhi = i;
                            }
                            else if (dots[i].x[dim] == medme.valuehi)
                            {   // tied for closest
                                medme.wthi += dots[i].wt;
                                medme.counthi++;
                            }
                        }
                    }
                    // combine median data struct across current subset of procs

                    counters[0]++;
                    sparta.mpi.MPI_Allreduce(ref medme, ref med, 1, med_type, med_op, comm);

                    // test median guess for convergence
                    // move additional dots that are next to cut across it

                    if (wtlo + med.totallo < targetlo)    // lower half TOO SMALL
                    {
                        wtlo += med.totallo;
                        valuehalf = med.valuehi;
                        if (med.counthi == 1)
                        {                  // only one dot to move
                            if (wtlo + med.wthi < targetlo)
                            {  // move it, keep iterating
                                if (me == med.prochi) dotmark[indexhi] = 0;
                            }
                            else
                            {                                 // only move if beneficial
                                if (wtlo + med.wthi - targetlo < targetlo - wtlo)
                                    if (me == med.prochi) dotmark[indexhi] = 0;
                                break;                               // all done
                            }
                        }
                        else
                        {                                   // multiple dots to move
                            breakflag = 0;
                            wtok = 0.0;
                            if (medme.valuehi == med.valuehi) wtok = medme.wthi;
                            if (wtlo + med.wthi >= targetlo)
                            {                // all done
                                sparta.mpi.MPI_Scan(ref wtok, ref wtupto, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, comm);
                                wtmax = targetlo - wtlo;
                                if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
                                breakflag = 1;
                            }                                      // wtok = most I can move
                            for (j = 0, wtsum = 0.0; j < nlist && wtsum < wtok; j++)
                            {
                                i = dotlist[j];
                                if (dots[i].x[dim] == med.valuehi)
                                { // only move if better
                                    if (wtsum + dots[i].wt - wtok < wtok - wtsum)
                                        dotmark[i] = 0;
                                    wtsum += dots[i].wt;
                                }
                            }
                            if (breakflag!=0) break;                   // done if moved enough
                        }
                        wtlo += med.wthi;
                        if (targetlo - wtlo <= tolerance) break;  // close enough

                        valuemin = med.valuehi;                   // iterate again
                        markactive = true;
                    }
                    else if (wthi + med.totalhi < targethi)     // upper half TOO SMALL
                    {
                        wthi += med.totalhi;
                        valuehalf = med.valuelo;

                        if (med.countlo == 1)
                        {                  // only one dot to move
                            if (wthi + med.wtlo < targethi)
                            {  // move it, keep iterating
                                if (me == med.proclo) dotmark[indexlo] = 1;
                            }
                            else
                            {                                 // only move if beneficial
                                if (wthi + med.wtlo - targethi < targethi - wthi)
                                    if (me == med.proclo) dotmark[indexlo] = 1;
                                break;                               // all done
                            }
                        }
                        else
                        {                                   // multiple dots to move
                            breakflag = 0;
                            wtok = 0.0;
                            if (medme.valuelo == med.valuelo) wtok = medme.wtlo;
                            if (wthi + med.wtlo >= targethi)
                            {                // all done
                                sparta.mpi.MPI_Scan(ref wtok, ref wtupto, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, comm);
                                wtmax = targethi - wthi;
                                if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
                                breakflag = 1;
                            }                                      // wtok = most I can move
                            for (j = 0, wtsum = 0.0; j < nlist && wtsum < wtok; j++)
                            {
                                i = dotlist[j];
                                if (dots[i].x[dim] == med.valuelo)
                                { // only move if better
                                    if (wtsum + dots[i].wt - wtok < wtok - wtsum)
                                        dotmark[i] = 1;
                                    wtsum += dots[i].wt;
                                }
                            }
                            if (breakflag!=0) break;
                        }
                        wthi += med.wtlo;
                        if (targethi - wthi <= tolerance) break;  // close enough

                        valuemax = med.valuelo;                   // iterate again
                        markactive = false;
                    }
                    else                  // Goldilocks result: both partitions just right
                        break;
                    // shrink the active list

                    k = 0;
                    for (j = 0; j < nlist; j++)
                    {
                        i = dotlist[j];
                        if (dotmark[i] ==0&&markactive==false) dotlist[k++] = i;
                        if (dotmark[i] != 0 && markactive == true) dotlist[k++] = i;
                    }
                    nlist = k;
                }
                // found median
                // store cut info in tree only if I am procmid

                if (me == procmid)
                {
                    tree[me].dim = dim;
                    tree[me].cut = valuehalf;
                }
                // use cut to shrink RCB bounding box

                if (me < procmid) rcbbox.hi[dim] = valuehalf;
                else rcbbox.lo[dim] = valuehalf;

                // outgoing = number of dots to ship to partner
                // nkeep = number of dots that have never migrated

                markactive = (me < procpartner);
                for (i = 0, keep = 0, outgoing = 0; i < ndot; i++)
                    if (dotmark[i] == 0 && markactive == false) outgoing++;
                    else if (dotmark[i] != 0 && markactive == true) outgoing++;
                    else if (i < nkeep) keep++;
                nkeep = keep;
                // alert partner how many dots I'll send, read how many I'll recv

                sparta.mpi.MPI_Send(ref outgoing, 1, MPI.MPI_INT, procpartner, 0, sparta.world);
                incoming = 0;
                if (readnumber!=0)
                {
                    sparta.mpi.MPI_Recv(ref incoming, 1, MPI.MPI_INT, procpartner, 0, sparta.world, ref status);
                    if (readnumber == 2)
                    {
                        sparta.mpi.MPI_Recv(ref incoming2, 1, MPI.MPI_INT, procpartner2, 0, sparta.world, ref status);
                        incoming += incoming2;
                    }
                }

                // check if need to alloc more space

                int ndotnew = ndot - outgoing + incoming;
                if (ndotnew > maxdot)
                {
                    while (maxdot < ndotnew) maxdot += DELTA;
                    dots = new Dot[maxdot];
                    counters[6]++;
                }
                counters[1] += outgoing;
                counters[2] += incoming;
                if (ndotnew > counters[3]) counters[3] = ndotnew;
                if (maxdot > counters[4]) counters[4] = maxdot;

                // malloc comm send buffer

                if (outgoing > maxbuf)
                {
                    //memory->sfree(buf);
                    maxbuf = outgoing;
                    buf = new Dot[maxbuf];
                }
                // fill buffer with dots that are marked for sending
                // pack down the unmarked ones

                keep = outgoing = 0;
                for (i = 0; i < ndot; i++)
                {
                    if ((dotmark[i] ==0&& markactive==false)|| (dotmark[i] != 0 && markactive == true))
                        buf[outgoing++] = dots[i];
                    else
                        dots[keep++] = dots[i];
                }// post receives for dots

                if (readnumber > 0)
                {
                    sparta.mpi.MPI_Irecv(ref dots[keep], incoming * Marshal.SizeOf(typeof(Dot)), MPI.MPI_CHAR,
                              procpartner, 1, sparta.world, ref request);
                    if (readnumber == 2)
                    {
                        keep += incoming - incoming2;
                        sparta.mpi.MPI_Irecv(ref dots[keep], incoming2 * Marshal.SizeOf(typeof(Dot)), MPI.MPI_CHAR,
                                      procpartner2, 1, sparta.world, ref request2);
                    }
                }

                // handshake before sending dots to insure recvs have been posted

                if (readnumber > 0)
                {
                    //sparta.mpi.MPI_Send(null, 0, MPI.MPI_INT, procpartner, 0, sparta.world);
                    //if (readnumber == 2) sparta.mpi.MPI_Send(null, 0, MPI.MPI_INT, procpartner2, 0, sparta.world);
                    Console.WriteLine("MPI Stub WARNING: Should not send message to self\n");


                }
                //MPI_Recv(NULL, 0, MPI_INT, procpartner, 0, world, &status);
                Console.WriteLine("MPI Stub WARNING: Should not recv message from self\n");
                // send dots to partner

                //MPI_Rsend(buf, outgoing * sizeof(Dot), MPI_CHAR, procpartner, 1, world);
                Console.WriteLine("MPI Stub WARNING: Should not recv message from self\n");
                // wait until all dots are received

                if (readnumber > 0)
                {
                    //MPI_Wait(&request, &status);
                    //if (readnumber == 2) MPI_Wait(&request2, &status);
                    Console.WriteLine("MPI Stub WARNING: Should not recv message from self\n");
                }

                ndot = ndotnew;

                // cut partition in half, create new communicators of 1/2 size

                int split;
                if (me < procmid)
                {
                    procupper = procmid - 1;
                    split = 0;
                }
                else
                {
                    proclower = procmid;
                    split = 1;
                }

                sparta.mpi.MPI_Comm_split(comm, split, me, ref comm_half);
                sparta.mpi.MPI_Comm_free(ref comm);
                comm = comm_half;
            }
            // clean up

            sparta.mpi.MPI_Comm_free(ref comm);

            // set public variables with results of rebalance

            nfinal = ndot;

            if (nfinal > maxrecv)
            {
                maxrecv = nfinal;
                recvproc = new int[maxrecv];
                recvindex = new int[maxrecv];
                
            }

            for (i = 0; i < nfinal; i++)
            {
                recvproc[i] = dots[i].proc;
                recvindex[i] = dots[i].index;
            }

            lo = rcbbox.lo;
            hi = rcbbox.hi;
        }
        public void invert()
        {
            if (irregular==null) irregular = new Irregular(sparta);

            // nsend = # of dots to request from other procs

            int nsend = nfinal - nkeep;

            int[] proclist;
            proclist = new int[nsend];

            Invert[] sinvert = new Invert[nsend];

            int m = 0;
            for (int i = nkeep; i < nfinal; i++)
            {
                proclist[m] = recvproc[i];
                sinvert[m].rindex = recvindex[i];
                sinvert[m].sproc = me;
                sinvert[m].sindex = i;
                m++;
            }
            // perform inversion via irregular comm
            // nrecv = # of my dots to send to other procs

            int nrecv = irregular.create_data_uniform(nsend, proclist,sparta.comm.commsortflag);
            Invert[] rinvert = new Invert[nrecv]; 
            irregular.exchange_uniform(sinvert.ToString(), Marshal.SizeOf(typeof(Invert)),rinvert.ToString());

            // set public variables from requests to send my dots
        }
        //public void check();
        //public void stats(int);


        void box_merge(IntPtr inn, IntPtr inout,ref int a,ref MPI_Datatype b)
        {
            //BBox box1 = (BBox) inn;
            //BBox box2 = (BBox)inout;

            //for (int i = 0; i < 3; i++)
            //{
            //    if (box1->lo[i] < box2->lo[i])
            //        box2->lo[i] = box1->lo[i];
            //    if (box1->hi[i] > box2->hi[i])
            //        box2->hi[i] = box1->hi[i];
            //}
            Console.WriteLine("RCB.box_merge()->don't know how to fix it");
        }

        void median_merge(IntPtr inn, IntPtr inout, ref int a, ref MPI_Datatype b)

        {
            //RCB::Median* med1 = (RCB::Median*) in;
            //RCB::Median* med2 = (RCB::Median*)inout;

            //med2->totallo += med1->totallo;
            //if (med1->valuelo > med2->valuelo)
            //{
            //    med2->valuelo = med1->valuelo;
            //    med2->wtlo = med1->wtlo;
            //    med2->countlo = med1->countlo;
            //    med2->proclo = med1->proclo;
            //}
            //else if (med1->valuelo == med2->valuelo)
            //{
            //    med2->wtlo += med1->wtlo;
            //    med2->countlo += med1->countlo;
            //    if (med1->proclo < med2->proclo) med2->proclo = med1->proclo;
            //}

            //med2->totalhi += med1->totalhi;
            //if (med1->valuehi < med2->valuehi)
            //{
            //    med2->valuehi = med1->valuehi;
            //    med2->wthi = med1->wthi;
            //    med2->counthi = med1->counthi;
            //    med2->prochi = med1->prochi;
            //}
            //else if (med1->valuehi == med2->valuehi)
            //{
            //    med2->wthi += med1->wthi;
            //    med2->counthi += med1->counthi;
            //    if (med1->prochi < med2->prochi) med2->prochi = med1->prochi;
            //}
            Console.WriteLine("RCB.median_merge()->don't know how to fix it, too");
        }

        // RCB cut info

        public struct Median
        {
            public double totallo, totalhi;   // weight in each half of active partition
            public double valuelo, valuehi;   // position of dot(s) nearest to cut
            public double wtlo, wthi;         // total weight of dot(s) at that position
            public int countlo, counthi;      // # of dots at that position
            public int proclo, prochi;       // unique proc who owns a nearest dot
        }

        // bounding box

        public struct BBox
        {
            public double[] lo;
            public double[] hi;
        }
    }
}
