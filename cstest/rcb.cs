using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MPI_Op = System.Int32;
using MPI_Datatype = System.Int32;
using System.Runtime.InteropServices;

namespace cstest
{
    class RCB
    {
        int me, nprocs;

        // point to balance on

        struct Dot
        {
            double[] x;          // coord of point
            double wt;            // weight of point
            int proc;             // owning proc
            int index;            // index on owning proc
        };

        // tree of RCB cuts

        struct Tree
        {
            double cut;         // position of cut
            int dim;            // dimension = 0/1/2 of cut
        };

        // inversion message

        struct Invert
        {
            int rindex;         // index on receiving proc
            int sproc;          // sending proc
            int sindex;         // index on sending proc
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
        //void compute(int, double**, double*, int flip = 0);
        //void invert();
        //void check();
        //void stats(int);


        void box_merge(IntPtr inn, IntPtr inout,ref int a,ref MPI_Datatype b)
        {
            BBox box1 = (BBox) inn;
            BBox box2 = (BBox)inout;

            for (int i = 0; i < 3; i++)
            {
                if (box1->lo[i] < box2->lo[i])
                    box2->lo[i] = box1->lo[i];
                if (box1->hi[i] > box2->hi[i])
                    box2->hi[i] = box1->hi[i];
            }
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
