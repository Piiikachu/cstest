
using System.Collections.Generic;
using bigint = System.Int64;
namespace cstest
{
    public class Surf
    {
        enum Enum1{ TALLYAUTO, TALLYREDUCE, TALLYLOCAL };         // same as Update
        enum Enum2{ REGION_ALL, REGION_ONE, REGION_CENTER };      // same as Grid
        enum Enum3{ TYPE, MOLECULE, ID };
        enum Enum4{ LT, LE, GT, GE, EQ, NEQ, BETWEEN };

        public const int DELTA = 4;
        public const double EPSSQ = 1.0e-12;
        public const double BIG = 1.0e20;
        public const int MAXGROUP = 32;

        public int exist;                // 1 if any surfaces are defined, else 0
        public int surf_collision_check; // flag for whether init() check is required
                                    // for assign of collision models to surfs
          
        public double[] bblo, bbhi=new double[3];   // bounding box around surfs
        public int tally_comm;           // style of comm for surf tallies
          
        public int nreact_one;           // surface reactions in current step
        public bigint nreact_running;    // running count of surface reactions
          
        public int ngroup;               // # of defined groups
        public string[] gnames;            // name of each group
        public int[] bitmask;             // one-bit mask for each group
        public int[] inversemask;         // inverse mask for each group

        public struct Point
        {
            double[] x;
        };

        public struct Line
        {
            int type, mask;          // type and mask of the element
            int isc, isr;            // index of surface collision and reaction models
                                     // -1 if unassigned
            int p1, p2;              // indices of points in line segment
                                     // rhand rule: Z x (p2-p1) = outward normal
            double[] norm;         // outward normal to line segment
        };

        public struct Tri
        {
            int type, mask;          // type and mask of the element
            int isc, isr;            // index of surface collision and reaction models
                                     // -1 if unassigned
            int p1, p2, p3;           // indices of points in triangle
                                      // rhand rule: (p2-p1) x (p3-p1) = outward normal
            double[] norm;         // outward normal to triangle
        };

        public Point[] pts;               // global list of points
        public Line[] lines;              // global list of lines
        public Tri[] tris;                // global list of tris
        public int npoint, nline, ntri;    // number of each
          
        public int[] mysurfs;             // indices of surf elements I own
        public int nlocal;               // # of surf elements I own
          
        public int nsc, nsr;              // # of surface collision and reaction models
        public List<SurfCollide> sc;   // list of surface collision models
        public List<SurfReact> sr;     // list of surface reaction models
          
        public int pushflag;             // set to 1 to push surf pts near grid cell faces
        public double pushlo, pushhi;     // lo/hi ranges to push on
        public double pushvalue;         // new position to push to

        public Surf(SPARTA sparta)
        {
            exist = 0;
            surf_collision_check = 1;

            gnames = new string[MAXGROUP] ;
            bitmask = new int[MAXGROUP];
            inversemask = new int[MAXGROUP];

            for (int i = 0; i < MAXGROUP; i++) bitmask[i] = 1 << i;
            for (int i = 0; i < MAXGROUP; i++) inversemask[i] = bitmask[i] ^ ~0;

            ngroup = 1;
            int n = "all".Length + 1;
            gnames[0] = string.Copy("all");

            npoint = nline = ntri = 0;
            //pts = null;
            //lines = null;
            //tris = null;
            pushflag = 1;

            nlocal = 0;
            mysurfs = null;

            nsc = maxsc = 0;
            sc = null;

            nsr = maxsr = 0;
            sr = null;

            tally_comm = (int)Enum1.TALLYAUTO;
        }


        //void modify_params(int, char**);
        //void init();
        //int nelement();
        //void setup_surf();

        //void compute_line_normal(int, int);
        //void compute_tri_normal(int, int);
        //void quad_corner_point(int, double*, double*, double*);
        //void hex_corner_point(int, double*, double*, double*);
        //double line_size(int);
        //double axi_line_size(int);
        //double tri_size(int, double &);

        //void check_watertight_2d(int, int);
        //void check_watertight_3d(int, int);
        //void check_point_inside(int, int);

        //void add_collide(int, char**);
        //int find_collide(const char*);
        //void add_react(int, char**);
        //int find_react(const char*);

        //void group(int, char**);
        //int add_group(const char*);
        //int find_group(const char*);

        //void collate_vector(int, int*, double*, int, double*);
        //void collate_array(int, int, int*, double**, double**);

        //void write_restart(FILE*);
        //void read_restart(FILE*);
        //bigint memory_usage();

        //private:
        private int maxsc;                // max # of models in sc
        private int maxsr;                // max # of models in sr

        //private void collate_vector_allreduce(int, int*, double*, int, double*);
        //private void collate_vector_irregular(int, int*, double*, int, double*);
        //private void collate_array_allreduce(int, int, int*, double**, double**);
        //private void collate_array_irregular(int, int, int*, double**, double**);
        
    }
}