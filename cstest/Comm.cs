using bigint = System.Int64;
namespace cstest
{
    public class Comm
    {
        
        public int me, nprocs;                    // proc info
        public bigint ncomm;                     // dummy statistic for now
         
        public int commsortflag;                 // 1 to force sort in all irregular comms
                                           //   useful for debugging to insure
                                           //   reproducible ordering of recv datums
        public int commpartstyle;                // 1 for neighbor, 0 for all
                                          //   changes how irregular comm for
                                          //   particles is performed


        //void init() { }
        //void reset_neighbors();
        //int migrate_particles(int, int*);
        //virtual void migrate_cells(int);
        //int send_cells_adapt(int, int*, char*, char**);
        //int irregular_uniform(int, int*, char*, int, char**);
        //void ring(int, int, void*, int, void (*)(int, char*),
        //          void*, int self = 1);

        
        protected Irregular iparticle,igrid,iuniform;
        protected string sbuf,rbuf;
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