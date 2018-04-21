
using System;
using System.Collections.Generic;
using System.Linq;
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
          
        public double[] bblo=new double[3], bbhi=new double[3];   // bounding box around surfs
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
        public int nelement()
        {
            if (sparta.domain.dimension == 2) return nline;
            return ntri;
        }
        public void setup_surf()
        {
            int me = sparta.comm.me;
            int nprocs = sparta.comm.nprocs;

            int n = nelement();

            // assign every Pth surf element to this proc

            nlocal = n / nprocs;
            if (me < n % nprocs) nlocal++;

            //memory->destroy(mysurfs);
            //memory->create(mysurfs, nlocal, "surf:mysurfs");
            mysurfs = new int[nlocal];

            nlocal = 0;
            for (int m = me; m < n; m += nprocs)
                mysurfs[nlocal++] = m;

            // set bounding box of all surfs based on pts
            // for 2d, set zlo,zhi to box bounds

            int i;
            for (i = 0; i < 3; i++)
            {
                bblo[i] = BIG;
                bbhi[i] = -BIG;
            }

            double[] x;
            for (int ipt = 0; ipt < npoint; ipt++)
            {
                x = pts[ipt].x;
                for (i = 0; i < 3; i++)
                {
                    bblo[i] =Math.Min(bblo[i], x[i]);
                    bbhi[i] =Math.Max(bbhi[i], x[i]);
                }
            }

            if (sparta.domain.dimension == 2)
            {
                bblo[2] = sparta.domain.boxlo[2];
                bbhi[2] = sparta.domain.boxhi[2];
            }
        }

        public void compute_line_normal(int nstart, int n)
        {
            int p1, p2;
            double[] z=new double[3], delta = new double[3], norm = new double[3];

            z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;

            int m = nstart;
            for (int i = 0; i < n; i++)
            {
                p1 = lines[m].p1;
                p2 = lines[m].p2;
                MathExtra.sub3(pts[p2].x, pts[p1].x, delta);
                MathExtra.cross3(z, delta, norm);
                MathExtra.norm3(norm);
                norm[2] = 0.0;
                Line line = lines[m];
                line.norm= norm;
                lines[m] = line; 
                        
                m++;
            }
        }
        public void compute_tri_normal(int nstart, int n)
        {
            int p1, p2, p3;
            double[] delta12=new double[3], delta13 = new double[3];

            int m = nstart;
            for (int i = 0; i < n; i++)
            {
                p1 = tris[m].p1;
                p2 = tris[m].p2;
                p3 = tris[m].p3;
                MathExtra.sub3(pts[p2].x, pts[p1].x, delta12);
                MathExtra.sub3(pts[p3].x, pts[p1].x, delta13);
                MathExtra.cross3(delta12, delta13, tris[m].norm);
                MathExtra.norm3(tris[m].norm);
                m++;
            }
        }
        //public void quad_corner_point(int, double*, double*, double*);
        //public void hex_corner_point(int, double*, double*, double*);
        public double line_size(int m)
        {
            double[] delta=new double[3];
            MathExtra.sub3(pts[lines[m].p2].x, pts[lines[m].p1].x, delta);
            return MathExtra.len3(delta);
        }
        //public double axi_line_size(int);
        public double tri_size(int m,out double len)
        {
            double[] delta12=new double[3], delta13 = new double[3], 
                delta23 = new double[3], cross = new double[3];

            MathExtra.sub3(pts[tris[m].p2].x, pts[tris[m].p1].x, delta12);
            MathExtra.sub3(pts[tris[m].p3].x, pts[tris[m].p1].x, delta13);
            MathExtra.sub3(pts[tris[m].p3].x, pts[tris[m].p2].x, delta23);
            len =Math.Min(MathExtra.len3(delta12), MathExtra.len3(delta13));
            len =Math.Min(len, MathExtra.len3(delta23));

            MathExtra.cross3(delta12, delta13, cross);
            double area = 0.5 * MathExtra.len3(cross);
            return area;
        }

        public void check_watertight_2d(int npoint_old, int nline_old)
        {
            int p1, p2;

            int npoint_new = npoint - npoint_old;
            int nline_new = nline - nline_old;

            // count[I] = # of lines that vertex I is part of

            int[] count;
            count = new int[npoint_new];
            //memory->create(count, npoint_new, "readsurf:count");
            for (int i = 0; i < npoint_new; i++) count[i] = 0;

            int ndup = 0;
            int m = nline_old;
            for (int i = 0; i < nline_new; i++)
            {
                p1 = lines[m].p1 - npoint_old;
                p2 = lines[m].p2 - npoint_old;
                count[p1]++;
                count[p2]++;
                if (count[p1] > 2) ndup++;
                if (count[p2] > 2) ndup++;
                m++;
            }

            if (ndup!=0)
            {
                string str=string.Format("Surface check failed with {0} duplicate points", ndup);
                sparta.error.all( str);
            }

            // check that all counts are 2
            // allow for exception if point on box surface

            double[] boxlo = sparta.domain.boxlo;
            double[] boxhi = sparta.domain.boxhi;

            int nbad = 0;
            for (int i = 0; i < npoint_new; i++)
            {
                if (count[i] == 0) nbad++;
                else if (count[i] == 1)
                {
                    if (Geometry.point_on_hex(pts[i + npoint_old].x, boxlo, boxhi)==0) nbad++;
                }
            }

            if (nbad!=0)
            {
                string str=string.Format("Surface check failed with {0} unmatched points", nbad);
                sparta.error.all( str);
            }

            // clean up

           // memory->destroy(count);
        }
        public void check_watertight_3d(int a, int ntri_old)
        {
            int ntri_new = ntri - ntri_old;

            // hash directed edges of all triangles
            // key = directed edge, value = # of times it appears in any triangle
            // NOTE: could prealloc hash to correct size here


            Dictionary<bigint, int> hash = new Dictionary<bigint, int>();

            // insert each edge into hash
            // should appear once in each direction
            // error if any duplicate edges

            bigint p1, p2, p3, key;

            int ndup = 0;
            int m = ntri_old;
            for (int i = 0; i < ntri_new; i++)
            {
                p1 = tris[m].p1;
                p2 = tris[m].p2;
                p3 = tris[m].p3;
                key = (p1 << 32) | p2;
                if (hash[key] != hash[hash.Keys.Last()]) ndup++;
                else hash[key] = 0;
                key = (p2 << 32) | p3;
                if (hash[key] != hash[hash.Keys.Last()]) ndup++;
                else hash[key] = 0;
                key = (p3 << 32) | p1;
                if (hash[key] != hash[hash.Keys.Last()]) ndup++;
                else hash[key] = 0;
                m++;
            }

            if (ndup!=0)
            {
                string str=string.Format( "Surface check failed with {0} duplicate edges", ndup);
                sparta.error.all( str);
            }

            // check that each edge has an inverted match
            // allow for exception if edge on box surface

            double[] boxlo = sparta.domain.boxlo;
            double[] boxhi = sparta.domain.boxhi;

            int nbad = 0;
            foreach (var it in hash)
            {
                p1 = it.Key >> 32;
                p2 = it.Key & Run.MAXSMALLINT;
                key = (p2 << 32) | p1;
                if (hash[key] == hash[hash.Keys.Last()])
                {
                    if (Geometry.point_on_hex(pts[(int)p1].x, boxlo, boxhi)==0 ||
                        Geometry.point_on_hex(pts[(int)p2].x, boxlo, boxhi)==0) nbad++;
                }

            }
            

            if (nbad!=0)
            {
                string str=string.Format("Surface check failed with {0} unmatched edges", nbad);
                sparta.error.all( str);
            }

        }
        public void check_point_inside(int npoint_old, int npoint_new)
        {
            double[] x;

            double[] boxlo = sparta.domain.boxlo;
            double[] boxhi = sparta.domain.boxhi;

            int m = npoint_old;
            int nbad = 0;
            for (int i = 0; i < npoint_new; i++)
            {
                x = pts[m].x;
                if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
                x[1] < boxlo[1] || x[1] > boxhi[1] ||
                x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
                m++;
            }

            if (nbad!=0)
            {
                string str=string.Format("{0} surface points are not inside simulation box",nbad);
                sparta.error.all(str);
            }
        }
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
        public int add_group(string id)
        {
            if (ngroup == MAXGROUP)
                sparta.error.all("Cannot have more than 32 surface groups");

            int n = id.Length + 1;
            gnames[ngroup] = string.Copy(id);

            for (int i = 0; i < n - 1; i++)
                if (!char.IsLetterOrDigit(id[i]) && id[i] != '_')
                    sparta.error.all("Group ID must be alphanumeric or underscore characters");

            ngroup++;
            return ngroup - 1;
        }
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