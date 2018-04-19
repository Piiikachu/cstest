using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class Cut3d
    {
        enum Enum1{ UNKNOWN, OUTSIDE, INSIDE, OVERLAP };     // several files
        enum Enum2{ CTRI, CTRIFACE, FACEPGON, FACE };
        enum Enum3{ EXTERIOR, INTERIOR, BORDER };
        enum Enum4{ ENTRY, EXIT, TWO, CORNER };              // same as Cut2d

        public const double EPSCELL = 1.0e-10;   // tolerance for pushing surf pts to cell surface

        // cell ID for 2d or 3d cell

        //#define VERBOSE
        public const int VERBOSE_ID = 4873;
        //#define VERBOSE_ID 27810406321L

        public int npushmax;           // # of push options to try
        public int[] npushcell=new int[4];       // tally of cells that required surf point push
        private SPARTA sparta;
        public Cut3d(SPARTA sparta)
        {
            this.sparta = sparta;
            cut2d = new Cut2d(sparta, 0);
            for (int i = 0; i <= cut2d.npushmax; i++)
            {
                cut2d.npushcell[i] = 0;
            }
            path1 = new double[12, 3];
            path2 = new double[12, 3];
            //memory->create(path1, 12, 3, "cut3d:path1");
            //memory->create(path2, 12, 3, "cut3d:path2");

            npushmax = 2;    // if increase this, increase push vec size in cut3d.h

            pushlo_vec[0] = -1.0;
            pushhi_vec[0] = 1.0;
            pushvalue_vec[0] = 0.0;
            pushlo_vec[1] = -1.0;
            pushhi_vec[1] = 1.0;
            pushvalue_vec[1] = 1.0;

            /*
            pushlo_vec[0] = -1.0;
            pushhi_vec[0] = 1.0;
            pushvalue_vec[0] = 1.0;
            pushlo_vec[1] = -1.0;
            pushhi_vec[1] = 1.0;
            pushvalue_vec[1] = 0.0;
            */

            if (sparta.surf.pushflag == 0) npushmax = 0;
            if (sparta.surf.pushflag == 2) npushmax = 3;
            if (sparta.surf.pushflag == 2)
            {
                pushlo_vec[2] = sparta.surf.pushlo;
                pushhi_vec[2] = sparta.surf.pushhi;
                pushvalue_vec[2] = sparta.surf.pushvalue;
            }

            for (int i = 0; i <= npushmax; i++) npushcell[i] = 0;

            // DEBUG
            //totcell = totsurf = totvert = totedge = 0;
        }
        //public int surf2grid(cellint, double[], double[], int[], int);
        //public int surf2grid_list(cellint, double[], double[], int, int[], int[], int);
        public int split(int id_caller, double[] lo_caller, double[] hi_caller,
                 int nsurf_caller, int[] surfs_caller,
                 out double[] vols_caller, int[] surfmap,
                 int[] corners, out int xsub, double[] xsplit)
       {
            Console.WriteLine("cut3d.split");
            vols_caller = null;
            xsub = 0;
            return 0;
//            id = id_caller;
//            lo = lo_caller;
//            hi = hi_caller;
//            nsurf = nsurf_caller;
//            surfs = surfs_caller;

//            // perform cut/split
//            // first attempt is with no pushing of surface points
//            // if fails, then try again with push options
//            // if all push options fail, then print error message

//            int nsplit, errflag;
//            pushflag = 0;


//            // debug 
//            //if (id == VERBOSE_ID) npushmax = 0;

//            while (true)
//            {
//                errflag = add_tris();
//                if (errflag != 0)
//                {
//                    if (push_increment() != 0) continue;
//                    break;
//                }

////# ifdef VERBOSE
////                if (id == VERBOSE_ID) print_bpg("BPG after added tris");
////#endif

//                int grazeflag = clip_tris();

//                // DEBUG
//                //totcell++;
//                //totsurf += nsurf;
//                //totvert += verts.n;
//                //totedge += edges.n;

////# ifdef VERBOSE
////                if (id == VERBOSE_ID) print_bpg("BPG after clipped tris");
////#endif

//                // all triangles just touched cell surface
//                // mark corner points based on grazeflag and in/out tri orientation
//                // return vol = 0.0 for UNKNOWN/INSIDE, full cell vol for OUTSIDE
//                // vol is changed in Grid::set_inout() if OVERLAP cell corners are marked

//                if (empty != 0)
//                {
//                    if (pushflag != 0) npushcell[pushflag]++;

//                    int mark = UNKNOWN;
//                    if (grazeflag || inout == INSIDE) mark = INSIDE;
//                    else if (inout == OUTSIDE) mark = OUTSIDE;
//                    corners[0] = corners[1] = corners[2] = corners[3] =
//                      corners[4] = corners[5] = corners[6] = corners[7] = mark;


//                    /*
//                    double ctr[3];
//                    ctr[0] = 0.5*(lo[0]+hi[0]);
//                    ctr[1] = 0.5*(lo[1]+hi[1]);
//                    ctr[2] = 0.5*(lo[2]+hi[2]);
//                    int check = 0;
//                    if (mark == INSIDE && 
//                        (fabs(ctr[0]) > 1.0 || fabs(ctr[1]) > 1.0 || fabs(ctr[2]) > 1.0))
//                      check = 1;
//                    if (mark == OUTSIDE && 
//                        (fabs(ctr[0]) < 1.0 && fabs(ctr[1]) < 1.0 && fabs(ctr[2]) < 1.0))
//                      check = 1;
//                    if (mark == UNKNOWN) check = 1;
//                    if (check) {
//                      printf("BAD MARKING %d %g %g %g: mark %d: "
//                             "nsurf %d pushflag %d grazeflag %d inout %d\n",
//                             id,ctr[0],ctr[1],ctr[2],mark,nsurf,pushflag,grazeflag,inout);
//                    }
//                    */


//                    double vol = 0.0;
//                    if (mark == OUTSIDE) vol = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);

//                    vols.grow(1);
//                    vols[0] = vol;
//                    vols_caller = &vols[0];
//                    return 1;
//                }

//                ctri_volume();
//                errflag = edge2face();
//                if (errflag != 0)
//                {
//                    if (push_increment() != 0) continue;
//                    break;
//                }

//                double[] lo2d=new double[2], hi2d = new double[2];

//                for (int iface = 0; iface < 6; iface++)
//                {



//                    // debug 
//                    //if (id == VERBOSE_ID) printf("FACE %d\n",iface);



//                    if (facelist[iface].n)
//                    {
//                        face_from_cell(iface, lo2d, hi2d);
//                        edge2clines(iface);
//                        errflag = cut2d->split_face(id, iface, lo2d, hi2d);
//                        if (errflag != 0) break;
//                        errflag = add_face_pgons(iface);
//                        if (errflag != 0) break;
//                    }
//                    else
//                    {
//                        face_from_cell(iface, lo2d, hi2d);
//                        errflag = add_face(iface, lo2d, hi2d);
//                        if (errflag != 0) break;
//                    }
//                }

//                if (errflag != 0)
//                {
//                    if (push_increment()) continue;
//                    break;
//                }

//                remove_faces();

//# ifdef VERBOSE
//                if (id == VERBOSE_ID) print_bpg("BPG after faces");
//#endif

//                errflag = check();
//                if (errflag)
//                {
//                    if (push_increment()) continue;
//                    break;
//                }

//                walk();

//# ifdef VERBOSE
//                if (id == VERBOSE_ID) print_loops();
//#endif

//                errflag = loop2ph();
//                if (errflag != 0)
//                {
//                    if (push_increment()) continue;
//                    break;
//                }

//                nsplit = phs.n;
//                if (nsplit > 1)
//                {
//                    create_surfmap(surfmap);
//                    errflag = split_point(surfmap, xsplit, xsub);
//                }
//                if (errflag != 0)
//                {
//                    if (push_increment()) continue;
//                    break;
//                }

//                // successful cut/split
//                // set corners = OUTSIDE if corner pt is in list of edge points
//                // else set corners = INSIDE

//                int icorner;
//                double[] p1,p2;

//                corners[0] = corners[1] = corners[2] = corners[3] =
//                  corners[4] = corners[5] = corners[6] = corners[7] = INSIDE;

//                int nedge = edges.n;
//                for (int iedge = 0; iedge < nedge; iedge++)
//                {
//                    if (!edges[iedge].active) continue;
//                    p1 = edges[iedge].p1;
//                    p2 = edges[iedge].p2;
//                    icorner = corner(p1);
//                    if (icorner >= 0) corners[icorner] = OUTSIDE;
//                    icorner = corner(p2);
//                    if (icorner >= 0) corners[icorner] = OUTSIDE;
//                }

//                // store volumes in vector so can return ptr to it

//                vols.grow(nsplit);
//                for (int i = 0; i < nsplit; i++) vols[i] = phs[i].volume;
//                vols_caller = &vols[0];

//                break;
//            }

//            // could not perform cut/split -> fatal error
//            // print info about cell and final error message
//            // 2-letter prefix is which method encountered error
//            // NOTE: store errflag_original for no-push error to print this message?

//            if (errflag!=0)
//            {
//                failed_cell();

//                if (errflag == 1)
//                    sparta.error.one("FE: Found edge in same direction");
//                if (errflag == 2)
//                    sparta.error.one("EF: Singlet BPG edge not on cell face");
//                if (errflag == 3)
//                    sparta.error.one("EF: BPG edge on more than 2 faces");
//                if (errflag == 4)
//                    sparta.error.one("LP: No positive volumes in cell");
//                if (errflag == 5)
//                    sparta.error.one("LP: More than one positive volume with a negative volume");
//                if (errflag == 6)
//                    sparta.error.one("LP: Single volume is negative, inverse donut");
//                if (errflag == 7)
//                    sparta.error.one("SP: Could not find split point in split cell");

//                if (errflag == 11)
//                    sparta.error.one("CH: Vertex has less than 3 edges");
//                if (errflag == 12)
//                    sparta.error.one("CH: Vertex contains invalid edge");
//                if (errflag == 13)
//                    sparta.error.one("CH: Vertex contains edge that doesn't point to it");
//                if (errflag == 14)
//                    sparta.error.one("CH: Vertex contains duplicate edge");
//                if (errflag == 15)
//                    sparta.error.one("CH: Vertex pointers to last edge are invalid");
//                if (errflag == 16)
//                    sparta.error.one("CH: Edge not part of 2 vertices");
//                if (errflag == 17)
//                    sparta.error.one("CH: Edge part of same vertex twice");
//                if (errflag == 18)
//                    sparta.error.one("CH: Edge part of invalid vertex");
//                if (errflag == 19)
//                    sparta.error.one("CH: Edge part of invalid vertex");

//                if (errflag == 21)
//                    sparta.error.one("WB: Point appears first in more than one CLINE");
//                if (errflag == 22)
//                    sparta.error.one("WB: Point appears last in more than one CLINE");
//                if (errflag == 23)
//                    sparta.error.one("WB: Singlet CLINES point not on cell border");
//                if (errflag == 24)
//                    sparta.error.one("LP: No positive areas in cell");
//                if (errflag == 25)
//                    sparta.error.one("LP: More than one positive area with a negative area");
//                if (errflag == 26)
//                    sparta.error.one("LP: Single area is negative, inverse donut");
//            }

//            if (pushflag!=0) npushcell[pushflag]++;

//            return nsplit;
        }
        //public int clip_external(double[], double[], double[],
        //                  double[], double[], double[]);


        private int id;            // ID of cell being worked on
        private double[] lo,hi;        // opposite corner pts of cell
        private int nsurf;             // # of surf elements in cell
        private int[] surfs;            // indices of surf elements in cell
         
        private int pushflag;          // 0 for no push, else push surf points near cell surf
        private double pushlo, pushhi;  // lo/hi ranges to push on
        private double pushvalue;      // new position to push to
        private double[] pushlo_vec=new double[3], pushhi_vec = new double[3]
            , pushvalue_vec = new double[3];  // push values to try
        private int inout;             // orientation of triangles that just touch cell
         
        private double[,] path1,path2;

  // DEBUG
  //int totcell,totsurf,totvert,totedge;

         List<double> vols;    // vols of each flow polyhedra found

        int empty;

        struct Vertex
        {
            int active;      // 1/0 if active or not
            int style;       // CTRI or CTRIFACE or FACEPGON or FACE
            int label;       // index in list of tris that intersect this cell
                             //   for CTRI or CTRIFACE
                             // face index (0-5) for FACEPGON or FACE
            int next;        // index of next vertex when walking a loop
            int nedge;       // # of edges in this vertex
            int first;       // first edge in vertex
            int dirfirst;    // dir of first edge in vertex
            int last;        // last edge in vertex
            int dirlast;     // dir of last edge in vertex
            double volume;   // volume of vertex projected against lower z face of cell
            double[] norm;    // ptr to norm of tri, NULL for other styles
        };

        struct Edge
        {
            double[] p1, p2;  // 2 points in edge
            int active;          // 1/0 if active or not
            int style;           // CTRI or CTRIFACE or FACEPGON or FACE
            int clipped;         // 1/0 if already clipped during face iteration
            int nvert;           // flag for verts containing this edge
                                 // 0 = no verts
                                 // 1 = just 1 vert in forward dir
                                 // 2 = just 1 vert in reverse dir
                                 // 3 = 2 verts in both dirs
                                 // all vecs are [0] in forward dir, [1] in reverse dir
            int[] verts;        // index of vertices containing this edge, -1 if not
            int[] next;         // index of next edge for each vertex, -1 for end
            int[] dirnext;      // next edge for each vertex is forward/reverse (0,1)
            int[] prev;         // index of prev edge for each vertex, -1 for start
            int[] dirprev;      // prev edge for each vertex is forward/reverse (0,1)
            //public Edge()
            //{
            //    p1 = new double[3];
            //    p2 = new double[3];
            //    int[] verts[2];
            //    int[] next[2];
            //    int[] dirnext[2];
            //    int[] prev[2];
            //    int[] dirprev[2];
            //}
        };

        struct Loop
        {
            double volume;        // volume of loop
            int flag;             // INTERIOR (if all CTRI vertices) or BORDER
            int n;                // # of vertices in loop
            int first;            // index of first vertex in loop
            int next;             // index of next loop in same PH, -1 if last loop
        };

        struct PH
        {
            double volume;
            int n;
            int first;
        };

        List<Vertex> verts;       // list of vertices in BPG
        List<Edge> edges;         // list of edges in BPG
        List<Loop> loops;         // list of loops of vertices = polyhedra
        List<PH> phs;             // list of polyhedrons = one or more loops

        List<int> facelist=new List<int>(6);    // list of edges on each cell face
        List<int> used;           // 0/1 flag for each vertex when walking loops
        List<int> stack;          // list of vertices to check when walking loops

        Cut2d cut2d;

        //int clip(double[], double[], double[]);
        //int add_tris();
        //int clip_tris();
        //void ctri_volume();
        //int edge2face();
        //void edge2clines(int);
        //int add_face_pgons(int);
        //int add_face(int, double[], double[]);
        //void remove_faces();
        //int check();
        //void walk();
        //int loop2ph();
        //void create_surfmap(int[]);
        //int split_point(int[], double[], int &);

        //void edge_insert(int, int, int, int, int, int, int);
        //void edge_remove(Edge*);
        //void edge_remove(Edge*, int);
        //void vertex_remove(Vertex*);
        //int grazing(Vertex*);
        //int which_faces(double[], double[], int[]);
        //void face_from_cell(int, double[], double[]);
        //void compress2d(int, double[], double[]);
        //void expand2d(int, double, double[], double[]);

        //int findedge(double[], double[], int, int &);
        //void between(double[], double[], int, double, double[]);
        //int samepoint(double[], double[]);
        //int corner(double[]);
        //int ptflag(double[]);
        //int push_increment();
        //void push(double[]);

        //void failed_cell();
        //void print_bpg(const char*);
        //void print_loops();
    }
}
