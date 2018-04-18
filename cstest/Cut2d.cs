using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class Cut2d
    {
        
        public int npushmax;          // # of push options to try
        public int[] npushcell=new int[4];      // tally of cells that required surf point push

        public struct Cline
        {
            public double[] x, y;   // coords of end points of line clipped to cell
            public int line;           // index in list of lines that intersect this cell
        };

        public struct Point
        {
            public double[] x;        // coords of point
            public int type;           // type of pt = ENTRY,EXIT,TWO,CORNER
                                       // ENTRY/EXIT = only part of one Cline,
                                       //   could also be a geometric corner pt
                                       // TWO = part of two Clines
                                       // CORNER = part of no Cline, is a geometric corner pt
            public int next;           // index of next point when walking a loop
                                       // set for ENTRY/TWO pts between ENTRY and EXIT
                                       // set for EXIT/CORNER points around cell perimeter,
                                       //   though may not be walked
            public int line;           // original line (as stored by Cline) the pt starts,
                                       //   only set for ENTRY and TWO pts
            public int corner;         // 1,2,3,4 if x is a corner point, else 0
                                       // could be ENTRY,EXIT,CORNER pt, but not a TWO pt
            public int cprev, cnext;    // indices of pts in linked list around cell perimeter
            public int side;           // which side of cell (0,1,2,3) pt is on
                                       // only for ENTRY/EXIT/CORNER pts to make linked list
            public double value;       // coord along the side
                                // only for ENTRY/EXIT/CORNER pts to make linked list
        };

        public struct Loop
        {
            public double area;        // area of loop
            public int active;         // 1/0 if active or not
            public int flag;           // INTERIOR (if all TWO points) or BORDER
            public int n;              // # of points in loop
            public int first;          // index of first point in loop
            public int next;           // index of next loop in same PG, -1 if last loop
        };

        public struct PG
        {
            public double area;        // summed area (over loops) of PG
            public int n;              // # of loops in PG
            public int first;          // index of first loop in PG
        };

        public List<Cline> clines;  // list of Clines
        public List<Point> points;  // list of Points = Weiler/Atherton data structure
        public List<Loop> loops;    // list of loops in Points
        public List<PG> pgs;        // list of polygons = one or more loops
        
        public Cut2d(SPARTA sparta, int caller_axisymmetric)
        {
            axisymmetric = caller_axisymmetric;

            npushmax = 2;    // if increase this, increase push vec size in cut2d.h

            pushlo_vec[0] = -1.0;
            pushhi_vec[0] = 1.0;
            pushvalue_vec[0] = 0.0;
            pushlo_vec[1] = -1.0;
            pushhi_vec[1] = 1.0;
            pushvalue_vec[1] = 1.0;

            if (sparta.surf.pushflag == 0) npushmax = 0;
            if (sparta.surf.pushflag == 2) npushmax = 3;
            if (sparta.surf.pushflag == 2)
            {
                pushlo_vec[2] = sparta.surf.pushlo;
                pushhi_vec[2] = sparta.surf.pushhi;
                pushvalue_vec[2] = sparta.surf.pushvalue;
            }

            for (int i = 0; i <= npushmax; i++) npushcell[i] = 0;
        }
        //public int surf2grid(cellint, double*, double*, int*, int);
        //public int surf2grid_list(cellint, double*, double*, int, int*, int*, int);
        //public int split(cellint, double*, double*, int, int*,
        //           double*&, int*, int*, int &, double*);
        //public int split_face(int, int, double*, double*);
        //public int clip_external(double*, double*, double*, double*, double*);

        
        private int axisymmetric;
        private int id;            // ID of cell being worked on
        private double[] lo,hi;        // opposite corner pts of cell
        private int nsurf;             // # of surf elements in cell
        private int[] surfs;            // indices of surf elements in cell
         
        private int pushflag;          // 0 for no push, else push surf points near cell surf
        private double pushlo, pushhi;  // lo/hi ranges to push on
        private double pushvalue;      // new position to push to
        private double[] pushlo_vec=new double[3], pushhi_vec = new double[3], pushvalue_vec = new double[3];  // push values to try
        private int inout;             // orientation of lines that just touch cell
         
        private List<double> areas;   // areas of each flow polygon found
        private List<int> used;       // 0/1 flag for each point when walking loops
         
        //private int build_clines();
        //private int weiler_build();
        //private void weiler_loops();
        //private int loop2pg();
        //private void create_surfmap(int*);
        //private int split_point(int*, double*, int &);
         
        //private int cliptest(double*, double*);
        //private void clip(double*, double*, double*, double*);
         
        //private int ptflag(double*);
        //private int push_increment();
        //private void push(double*);
        //private int sameedge(double*, double*);
        //private int whichside(double*);
         
        //private void failed_cell();
        //private void print_clines();
        //private void print_points();
        //private void print_loops();
    }
}
