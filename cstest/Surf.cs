
using System;
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
            public double[] x;
        };

        public struct Line
        {
            public int type, mask;          // type and mask of the element
            public int isc, isr;            // index of surface collision and reaction models
                                      // -1 if unassigned
            public int p1, p2;              // indices of points in line segment
                                      // rhand rule: Z x (p2-p1) = outward normal
            public double[] norm;         // outward normal to line segment
        };

        public struct Tri
        {
            public int type, mask;          // type and mask of the element
            public int isc, isr;            // index of surface collision and reaction models
                                      // -1 if unassigned
            public int p1, p2, p3;           // indices of points in triangle
                                       // rhand rule: (p2-p1) x (p3-p1) = outward normal
            public double[] norm;         // outward normal to triangle
        };

        public List<Point> pts;               // global list of points
        public List<Line> lines;              // global list of lines
        public List<Tri> tris;                // global list of tris
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
            this.sparta = sparta;
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


        public void modify_params(int narg, string[] args)
        {
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);
            if (narg < 2) sparta.error.all("Illegal surf_modify command");
            int igroup = find_group(arg[0]);
            if (igroup < 0) sparta.error.all("Surf_modify surface group is not defined");
            int groupbit = bitmask[igroup];

            int iarg = 1;
            while (iarg < narg)
            {
                switch (arg[iarg])
                {
                    case "collide":
                        if (iarg + 2 > narg) sparta.error.all("Illegal surf_modify command");

                        int isc = find_collide(arg[iarg + 1]);
                        if (isc < 0) sparta.error.all("Could not find surf_modify sc-ID");

                        // set surf collision model for each surf in surface group

                        if (sparta.domain.dimension == 2)
                        {
                            for (int i = 0; i < nline; i++)
                                if (lines[i].mask!=0 & groupbit!=0)
                                {
                                    Line line = new Line();
                                    line = lines[i];
                                    line.isc = isc;
                                    lines[i] = line;
                                }
                        }
                        if (sparta.domain.dimension == 3)
                        {
                            for (int i = 0; i < ntri; i++)
                                if (tris[i].mask!=0 & groupbit!=0)
                                {
                                    Tri tri = new Tri();
                                    tri = tris[i];
                                    tri.isc = isc;
                                    tris[i] = tri;
                                }
                        }

                        iarg += 2;
                        break;
                    case "react":
                        if (iarg + 2 > narg) sparta.error.all("Illegal surf_modify command");

                        int isr;
                        if (string.Equals(arg[iarg + 1], "none")) isr = -1;
                        else
                        {
                            isr = find_react(arg[iarg + 1]);
                            if (isr < 0) sparta.error.all("Could not find surf_modify sr-ID");
                        }

                        // set surf reaction model for each surf in surface group

                        if (sparta.domain.dimension == 2)
                        {
                            for (int i = 0; i < nline; i++)
                                if (lines[i].mask!=0 & groupbit!=0)
                                {
                                    Line line = new Line();
                                    line = lines[i];
                                    line.isr = isr;
                                    lines[i] = line;
                                    
                                }
                        }
                        if (sparta.domain.dimension == 3)
                        {
                            for (int i = 0; i < ntri; i++)
                                if (tris[i].mask != 0 & groupbit != 0)
                                {
                                    Tri tri = new Tri();
                                    tri = tris[i];
                                    tri.isr = isr;
                                    tris[i] = tri;
                                }
                        }

                        iarg += 2;
                        break;
                    default:
                        sparta.error.all("Illegal surf_modify command");
                        break;
                }
            }
        }
        public void init()
        {
            // check that every element is assigned to a surf collision model
            // skip if caller turned off the check, e.g. BalanceGrid

            if (surf_collision_check!=0)
            {
                int flag = 0;
                if (sparta.domain.dimension == 2)
                {
                    for (int i = 0; i < nline; i++)
                        if (lines[i].isc < 0) flag++;
                }
                if (sparta.domain.dimension == 3)
                {
                    for (int i = 0; i < ntri; i++)
                        if (tris[i].isc < 0) flag++;
                }
                if (flag != 0)
                {
                    string str=string.Format( "{0} surface elements not assigned to a collision model", flag);
                    sparta.error.all(str);
                }
            }

            // initialize surf collision and reaction models

            for (int i = 0; i < nsc; i++) sc[i].init();
            for (int i = 0; i < nsr; i++) sr[i].init();
        }
        //public int nelement();
        //public void setup_surf();

        //public void compute_line_normal(int, int);
        //public void compute_tri_normal(int, int);
        //public void quad_corner_point(int, double*, double*, double*);
        //public void hex_corner_point(int, double*, double*, double*);
        //public double line_size(int);
        //public double axi_line_size(int);
        //public double tri_size(int, double &);

        //public void check_watertight_2d(int, int);
        //public void check_watertight_3d(int, int);
        //public void check_point_inside(int, int);
        private SPARTA sparta;

        public void add_collide(int narg, string[] args)
        {
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);
            if (narg < 2) sparta.error.all("Illegal surf_collide command");

            // error check

            for (int i = 0; i < nsc; i++)
                if (string.Equals(arg[0], sc[i].id))
                    sparta.error.all("Reuse of surf_collide ID");

            // extend SurfCollide list if necessary

            if (nsc == maxsc)
            {
                maxsc += DELTA;
                sc = new List<SurfCollide>(maxsc);
            }
            // create new SurfCollide class

            if (sparta.suffix_enable!=0)
            {
                if (sparta.suffix!=null)
                {
                    string estyle = string.Format("{0}{1}", arg[1], sparta.suffix);
                }
            }
            switch (arg[1])
            {
                case "diffuse":
                    sc.Add(new SurfCollideDiffuse(sparta, narg, arg));
                    break;
                default:
                    Console.WriteLine("Unrecognized surf_collide style");
                    break;
            }

            nsc++;

        }
        public int find_collide(string id)
        {
            int isc;
            for (isc = 0; isc < nsc; isc++)
                if (id.Equals(sc[isc].id)) break;
            if (isc == nsc) return -1;
            return isc;
        }
        //public void add_react(int, char**);
        public int find_react(string id)
        {
            int isr;
            for (isr = 0; isr < nsr; isr++)
                if (string.Equals(id, sr[isr].id)) break;
            if (isr == nsr) return -1;
            return isr;
        }

        //public void group(int, char**);
        //public int add_group(const char*);
        public int find_group(string id)
        {
            int igroup;
            for (igroup = 0; igroup < ngroup; igroup++)
                if (string.Equals(id, gnames[igroup])) break;
            if (igroup == ngroup) return -1;
            return igroup;
        }

        //public void collate_vector(int, int*, double*, int, double*);
        //public void collate_array(int, int, int*, double**, double**);

        //public void write_restart(FILE*);
        //public void read_restart(FILE*);
        //public bigint memory_usage();

        //private:
        private int maxsc;                // max # of models in sc
        private int maxsr;                // max # of models in sr

        //private void collate_vector_allreduce(int, int*, double*, int, double*);
        //private void collate_vector_irregular(int, int*, double*, int, double*);
        //private void collate_array_allreduce(int, int, int*, double**, double**);
        //private void collate_array_irregular(int, int, int*, double**, double**);
        
    }
}