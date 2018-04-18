using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class Cut3d
    {
        
        public int npushmax;           // # of push options to try
        public int[] npushcell=new int[4];       // tally of cells that required surf point push
         
        public Cut3d(SPARTA sparta)
        {
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
