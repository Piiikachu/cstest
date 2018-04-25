using System;
using System.Runtime.InteropServices;
using System.Text;
using bigint = System.Int64;
using MPI_Request = System.Int32;
namespace cstest
{
    public delegate void callbackHandler(int nsize,StringBuilder buf);
    public class Comm
    {
        enum Enum1{ PKEEP, PINSERT, PDONE, PDISCARD, PENTRY, PEXIT, PSURF };   // several files
        public int me, nprocs;                    // proc info
        public bigint ncomm;                     // dummy statistic for now
         
        public int commsortflag;                 // 1 to force sort in all irregular comms
                                           //   useful for debugging to insure
                                           //   reproducible ordering of recv datums
        public int commpartstyle;                // 1 for neighbor, 0 for all
                                                 //   changes how irregular comm for
                                                 //   particles is performed


        public void init()
        {

        }
        public void reset_neighbors()
        {
            neighflag = 0;
            if (commpartstyle !=0|| sparta.grid.clumped==0) return;
            neighflag = 1;

            if (neighlist == null)
            {
                neighlist = new int[nprocs];
                //memory->create(neighlist, nprocs, "comm:neighlist");
            }

            for (int i = 0; i < nprocs; i++) neighlist[i] = 0;

            Grid.ChildCell[] cells = sparta.grid.cells;
            int nglocal = sparta.grid.nlocal;
            int ntotal = nglocal + sparta.grid.nghost;

            for (int icell = nglocal; icell < ntotal; icell++)
                neighlist[cells[icell].proc] = 1;
            neighlist[me] = 0;

            nneigh = 0;
            for (int i = 0; i < nprocs; i++)
                if (neighlist[i]!=0) neighlist[nneigh++] = i;

            iparticle.create_procs(nneigh, neighlist, commsortflag);
        }
        public int migrate_particles(int nmigrate, int[] plist)
        {
            //int i, j;

            //Grid.ChildCell[] cells = sparta.grid.cells;
            //Particle.OnePart[] particles = sparta.particle.particles;

            //int ncustom = sparta.particle.ncustom;
            //int nbytes_particle = Marshal.SizeOf(typeof(Particle.OnePart));
            //int nbytes_custom = sparta.particle.sizeof_custom();
            //int nbytes = nbytes_particle + nbytes_custom;

            //// grow pproc and sbuf if necessary

            //if (nmigrate > maxpproc)
            //{
            //    maxpproc = nmigrate;
            //    //memory->destroy(pproc);
            //    //memory->create(pproc, maxpproc, "comm:pproc");
            //    pproc = new int[maxpproc];

            //}
            //if (nmigrate * nbytes > maxsendbuf)
            //{
            //    maxsendbuf = nmigrate * nbytes;
            //    //memory->destroy(sbuf);
            //    //memory->create(sbuf, maxsendbuf, "comm:sbuf");
            //    sbuf = new char[maxsendbuf];
            //}

            //// fill proclist with procs to send to
            //// pack sbuf with particles to migrate
            //// if flag == PDISCARD, particle is deleted but not sent
            //// change icell of migrated particle to owning cell on receiving proc
            //// nsend = particles that actually migrate
            //// if no custom attributes, pack particles directly via memcpy()
            //// else pack_custom() performs packing into sbuf

            //int nsend = 0;
            //int offset = 0;

            //if (ncustom==0)
            //{
            //    for (i = 0; i < nmigrate; i++)
            //    {
            //        j = plist[i];
            //        if (particles[j].flag == (int)Enum1.PDISCARD) continue;
            //        pproc[nsend++] = cells[particles[j].icell].proc;
            //        particles[j].icell = cells[particles[j].icell].ilocal;
            //        memcpy(&sbuf[offset], &particles[j], nbytes_particle);
            //        offset += nbytes_particle;
            //    }
            //}
            //else
            //{
            //    for (i = 0; i < nmigrate; i++)
            //    {
            //        j = plist[i];
            //        if (particles[j].flag == (int)Enum1.PDISCARD) continue;
            //        pproc[nsend++] = cells[particles[j].icell].proc;
            //        particles[j].icell = cells[particles[j].icell].ilocal;
            //        memcpy(&sbuf[offset], &particles[j], nbytes_particle);
            //        offset += nbytes_particle;
            //        sparta.particle.pack_custom(j, &sbuf[offset]);
            //        offset += nbytes_custom;
            //    }
            //}

            //// compress my list of particles

            //sparta.particle.compress_migrate(nmigrate, plist);
            //int ncompress = sparta.particle.nlocal;

            //// create or augment irregular communication plan
            //// nrecv = # of incoming particles

            //int nrecv;
            //if (neighflag!=0)
            //    nrecv = sparta.particle.augment_data_uniform(nsend, pproc);
            //else
            //    nrecv = sparta.particle.create_data_uniform(nsend, pproc, commsortflag);

            //// extend particle list if necessary

            //sparta.particle.grow(nrecv);

            //// perform irregular communication
            //// if no custom attributes, append recv particles directly to particle list
            //// else receive into rbuf, unpack particles one by one via unpack_custom()

            //if (ncustom==0)
            //    sparta.particle.
            //      exchange_uniform(sbuf, nbytes,
            //                       (char*)&sparta.particle.particles[sparta.particle.nlocal]);

            //else
            //{
            //    if (nrecv * nbytes > maxrecvbuf)
            //    {
            //        maxrecvbuf = nrecv * nbytes;
            //        memory->destroy(rbuf);
            //        memory->create(rbuf, maxrecvbuf, "comm:rbuf");
            //    }

            //    sparta.particle.exchange_uniform(sbuf, nbytes, rbuf);

            //    offset = 0;
            //    int nlocal = sparta.particle.nlocal;
            //    for (i = 0; i < nrecv; i++)
            //    {
            //        memcpy(&sparta.particle.particles[nlocal], &rbuf[offset], nbytes_particle);
            //        offset += nbytes_particle;
            //        sparta.particle.unpack_custom(&rbuf[offset], nlocal);
            //        offset += nbytes_custom;
            //        nlocal++;
            //    }
            //}

            //sparta.particle.nlocal += nrecv;
            //ncomm += nsend;
            //return ncompress;
            Console.WriteLine("comm.migrate_particles");
        }
        public virtual void migrate_cells(int nmigrate)
        {
            int i, n;

            Grid.ChildCell[] cells = sparta.grid.cells;
            int nglocal = sparta.grid.nlocal;

            // grow proc and size lists if needed

            if (nmigrate > maxgproc)
            {
                maxgproc = nmigrate;
                //memory->destroy(gproc);
                //memory->destroy(gsize);
                gproc = new int[maxgproc];
                gsize = new int[maxgproc];
                //memory->create(gproc, maxgproc, "comm:gproc");
                //memory->create(gsize, maxgproc, "comm:gsize");
            }

            // fill proclist with procs to send to
            // compute byte count needed to pack cells

            int nsend = 0;
            bigint boffset = 0;
            StringBuilder sb = new StringBuilder();
            for (int icell = 0; icell < nglocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                if (cells[icell].proc == me) continue;
                gproc[nsend] = cells[icell].proc;
                n = sparta.grid.pack_one(icell,1, 1, 0,ref sb);
                gsize[nsend++] = n;
                boffset += n;
            }

            if (boffset > Run.MAXSMALLINT)
                sparta.error.one("Migrate cells send buffer exceeds 2 GB");
            int offset = (int)boffset;

            // reallocate sbuf as needed

            if (offset > maxsendbuf)
            {
                //memory->destroy(sbuf);
                maxsendbuf = offset;
                sbuf = new char[maxsendbuf];
                //memory->create(sbuf, maxsendbuf, "comm:sbuf");
                //memset(sbuf, 0, maxsendbuf);
            }
            sbuf = new char[maxsendbuf];
            rbuf = new char[maxsendbuf];
            // pack cell info into sbuf
            // only called for unsplit and split cells I no longer own

            offset = 0;
            for (int icell = 0; icell < nglocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                if (cells[icell].proc == me) continue;
                char[] tmp = new char[maxsendbuf - offset];
                Array.Copy(sbuf, offset, tmp, 0, maxsendbuf - offset);
                StringBuilder sb1 = new StringBuilder($"{sbuf}"); 
                offset += sparta.grid.pack_one(icell, 1, 1, 1, ref sb1);
            }

            // compress my list of owned grid cells to remove migrated cells
            // compress particle list to remove particles in migrating cells
            // unset particle sorted since compress_rebalance() does
            //   and may receive particles which will then be unsorted

            if (nmigrate!=0)
            {
                sparta.grid.compress();
                sparta.particle.compress_rebalance();
            }
            else sparta.particle.sorted = 0;

            // create irregular communication plan with variable size datums
            // nrecv = # of incoming grid cells
            // recvsize = total byte size of incoming grid + particle info
            // DEBUG: append a sort=1 arg so that messages from other procs
            //        are received in repeatable order, thus grid cells stay in order

            if (igrid==null) igrid = new Irregular(sparta);
            int recvsize;
            int nrecv =
              igrid.create_data_variable(nmigrate, gproc, gsize,out recvsize, commsortflag);

            // reallocate rbuf as needed

            if (recvsize > maxrecvbuf)
            {
               // memory->destroy(rbuf);
                maxrecvbuf = recvsize;
                rbuf = new char[maxrecvbuf];
                //memory->create(rbuf, maxrecvbuf, "comm:rbuf");
                //memset(rbuf, 0, maxrecvbuf);
            }

            // perform irregular communication
            byte[] bytes = new byte[sbuf.Length * sizeof(char)];
            bytes = Encoding.UTF8.GetBytes(sbuf.ToString());
            byte[] byter = new byte[rbuf.Length * sizeof(char)];
            byter = Encoding.UTF8.GetBytes(rbuf.ToString());
            igrid.exchange_variable(bytes, gsize, byter);

            // unpack received grid cells with their particles

            offset = 0;
            for (i = 0; i < nrecv; i++)
            {
                offset += sparta.grid.unpack_one(new StringBuilder(), 1, 1);
            }
            Console.Write("comm.migrate_cells()->unpackone\n");
            for (int a = offset; a < rbuf.Length-1-offset; a++)
            {
                Console.Write("{0}",rbuf[a]);
            }
        }
        //public int send_cells_adapt(int, int*, char*, char**);
        //public int irregular_uniform(int, int*, char*, int, char**);
        public void ring(int n, int nper, StringBuilder inbuf, int messtag,
                callbackHandler callback, StringBuilder outbuf, int self)
        {
            MPI_Request request=0;
            MPI._MPI_Status status;

            int nbytes = n * nper;
            int maxbytes=0;
            sparta.mpi.MPI_Allreduce(ref nbytes, ref maxbytes, 1, MPI.MPI_INT, MPI.MPI_MAX, sparta.world);

            StringBuilder buf,bufcopy;
            bufcopy = new StringBuilder();
            buf = new StringBuilder();
            buf.Append(inbuf);

            int next = me + 1;
            int prev = me - 1;
            if (next == nprocs) next = 0;
            if (prev < 0) prev = nprocs - 1;

            for (int loop = 0; loop < nprocs; loop++)
            {
                if (me != next)
                {
                    sparta.mpi.MPI_Irecv(ref bufcopy, maxbytes, MPI.MPI_CHAR, prev, messtag, sparta.world, ref request);
                    sparta.mpi.MPI_Send(ref buf, nbytes, MPI.MPI_CHAR, next, messtag, sparta.world);
                    //sparta.mpi.MPI_Wait(&request, &status);
                    //sparta.mpi.MPI_Get_count(&status, MPI.MPI_CHAR, &nbytes);
                    //memcpy(buf, bufcopy, nbytes);
                    buf.Append(bufcopy);
                }
                if (self!=0 || loop != nprocs - 1) callback(nbytes / nper, buf);
            }

            if (outbuf!=null)
            {
                //memcpy(outbuf, buf, nbytes);
                outbuf.Append(buf);
            }

            //memory->destroy(buf);
            //memory->destroy(bufcopy);
        }


        protected Irregular iparticle,igrid,iuniform;
        protected char[] sbuf,rbuf;
        protected int maxsendbuf, maxrecvbuf;
        protected int[] pproc,gproc,gsize;
        protected int maxpproc, maxgproc;
         
        protected int neighflag;                    // 1 if nearest-neighbor particle comm
        protected int nneigh;                       // # of procs I own ghost cells of
        protected int[] neighlist;                   // list of ghost procs
         
        protected int copymode;                     // 1 if copy of class (prevents deallocation of
                                                    //  base class when child copy is destroyed)

        private SPARTA sparta;

        public Comm(SPARTA sparta)
        {
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);
            sparta.mpi.MPI_Comm_size(sparta.world, ref nprocs);

            ncomm = 0;
            commsortflag = 0;
            commpartstyle = 1;

            neighflag = 0;
            neighlist = null;

            iparticle = new Irregular(sparta);
            igrid = null;
            iuniform = null;

            pproc = null;
            maxpproc = 0;
            gproc = gsize = null;
            maxgproc = 0;
            sbuf = rbuf = null;
            maxsendbuf = maxrecvbuf = 0;

            copymode = 0;
        }
    }
}