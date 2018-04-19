using System;
using System.Collections.Generic;
using System.IO;

namespace cstest
{
    public class Domain
    {
        enum Enum1:int{ XLO, XHI, YLO, YHI, ZLO, ZHI, INTERIOR };         // several files
        enum Enum2:int{ PERIODIC, OUTFLOW, REFLECT, SURFACE, AXISYM };  // several files

        public const int DELTAREGION = 4;


        private int[] surf_collide = new int[6];              // index of SurfCollide model
        private int[] surf_react = new int[6];                // index of SurfReact model
                                                              // for each bflag = SURFACE boundary


        public int box_exist;                    // 0 = not yet created, 1 = exists
        public int dimension;                    // 2,3
        public int axisymmetric;                 // 1 for yes, 0 for no, only allowed in 2d
        public int boundary_collision_check;  // flag for whether init() check is required
                                              // for assign of collision models to boundaries

        public double[] boxlo = new double[3];
        public double[] boxhi = new double[3];         // box global bounds
        public double xprd, yprd, zprd;            // global box dimensions
        public double[] prd = new double[3];                    // array form of dimensions
          
        public int[] bflag = new int[6];                     // boundary flags
        public double[,] norm = new double[6, 3];                // boundary normals
          
        public int surfreactany;                 // 1 if any boundary has surf reactions
          
        public int copy, copymode;                // 1 if copy of class (prevents deallocation of
                                             //  base class when child copy is destroyed)
          
        public int nregion;                      // # of defined Regions
        public int maxregion;                    // max # regions can hold
        public List<Region> regions;          // list of defined Regions

        private SPARTA sparta;
        public Domain(SPARTA sparta)
        {
            this.sparta = sparta;
            box_exist = 0;
            dimension = 3;
            axisymmetric = 0;
            boundary_collision_check = 1;

            for (int i = 0; i < 6; i++) bflag[i] = (int)Enum2.PERIODIC;
            for (int i = 0; i < 6; i++) surf_collide[i] = surf_react[i] = -1;

            // surface normals of 6 box faces pointed inward towards particles

            norm[(int)Enum1.XLO,0] = 1.0; norm[(int)Enum1.XLO,1] = 0.0; norm[(int)Enum1.XLO,2] = 0.0;
            norm[(int)Enum1.XHI,0] = -1.0; norm[(int)Enum1.XHI,1] = 0.0; norm[(int)Enum1.XHI,2] = 0.0;
            norm[(int)Enum1.YLO,0] = 0.0; norm[(int)Enum1.YLO,1] = 1.0; norm[(int)Enum1.YLO,2] = 0.0;
            norm[(int)Enum1.YHI,0] = 0.0; norm[(int)Enum1.YHI,1] = -1.0; norm[(int)Enum1.YHI,2] = 0.0;
            norm[(int)Enum1.ZLO,0] = 0.0; norm[(int)Enum1.ZLO,1] = 0.0; norm[(int)Enum1.ZLO,2] = 1.0;
            norm[(int)Enum1.ZHI,0] = 0.0; norm[(int)Enum1.ZHI,1] = 0.0; norm[(int)Enum1.ZHI,2] = -1.0;

            nregion = maxregion = 0;
            regions = null;
            copy = copymode = 0;

        }
        public void init()
        {
            if (axisymmetric!=0 && dimension != 2)
                sparta.error.all("Axi-symmetry only allowed for 2d simulation");
            if (dimension == 2 && (bflag[(int)Enum1.ZLO] != (int)Enum2.PERIODIC || bflag[(int)Enum1.ZHI] != (int)Enum2.PERIODIC))
                sparta.error.all("Z dimension must be periodic for 2d simulation");

            // check that every SURFACE boundary is assigned to a surf collision model
            // skip if caller turned off the check, e.g. BalanceGrid

            if (boundary_collision_check != 0)
            {
                for (int i = 0; i < 6; i++)
                    if (bflag[i] == (int)Enum2.SURFACE && surf_collide[i] < 0)
                        sparta.error.all("Box boundary not assigned a surf_collide ID");
            }

            int cutflag = 0;
            if (bflag[0] == (int)Enum2.PERIODIC && sparta.grid.cutoff > xprd) cutflag = 1;
            if (bflag[2] == (int)Enum2.PERIODIC && sparta.grid.cutoff > yprd) cutflag = 1;
            if (dimension == 3 && bflag[4] == (int)Enum2.PERIODIC && sparta.grid.cutoff > zprd)
                cutflag = 1;
            if (cutflag != 0) sparta.error.all("Grid cutoff is longer than box length in a periodic dimension");

            // surfreactany = 1 if any face has surface reactions assigned to it

            surfreactany = 0;
            for (int i = 0; i < 6; i++)
                if (surf_react[i] >= 0) surfreactany = 1;
        }
        public void set_initial_box()
        {
            if (boxlo[0] >= boxhi[0] || boxlo[1] >= boxhi[1] || boxlo[2] >= boxhi[2])
                sparta.error.all("Box bounds are invalid");
        }
        public void set_global_box()
        {
            prd[0] = xprd = boxhi[0] - boxlo[0];
            prd[1] = yprd = boxhi[1] - boxlo[1];
            prd[2] = zprd = boxhi[2] - boxlo[2];
        }
        public void set_boundary(int narg, string[] args)
        {
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);

            if (box_exist!=0)
                sparta.error.all( "Boundary command after simulation box is defined");

            if (narg != 3) sparta.error.all( "Illegal boundary command");

            char c;
            int m = 0;
            for (int idim = 0; idim < 3; idim++)
                for (int iside = 0; iside < 2; iside++)
                {
                    if (iside == 0) c = arg[idim][0];
                    else if (iside == 1 && arg[idim].Length == 1) c = arg[idim][0];
                    else c = arg[idim][1];

                    if (c == 'o') bflag[m] = (int)Enum2.OUTFLOW;
                    else if (c == 'p') bflag[m] = (int)Enum2.PERIODIC;
                    else if (c == 'r') bflag[m] = (int)Enum2.REFLECT;
                    else if (c == 's') bflag[m] = (int)Enum2.SURFACE;
                    else if (c == 'a') bflag[m] = (int)Enum2.AXISYM;
                    else sparta.error.all( "Illegal boundary command");

                    surf_collide[m] = surf_react[m] = -1;
                    m++;
                }

            if (dimension == 2 && (bflag[(int)Enum1.ZLO] != (int)Enum2.PERIODIC || bflag[(int)Enum1.ZHI] != (int)Enum2.PERIODIC))
                sparta.error.all( "Z dimension must be periodic for 2d simulation");

            if (bflag[(int)Enum1.XLO] == (int)Enum2.AXISYM || bflag[(int)Enum1.XHI] == (int)Enum2.AXISYM ||
                bflag[(int)Enum1.YHI] == (int)Enum2.AXISYM || bflag[(int)Enum1.ZLO] == (int)Enum2.AXISYM || bflag[(int)Enum1.ZHI] == (int)Enum2.AXISYM)
                sparta.error.all( "Only ylo boundary can be axi-symmetric");

            if (bflag[(int)Enum1.YLO] == (int)Enum2.AXISYM)
            {
                axisymmetric = 1;
                if (bflag[(int)Enum1.YHI] == (int)Enum2.PERIODIC)
                    sparta.error.all( "Y cannot be periodic for axi-symmetric");
            }

            for (m = 0; m < 6; m += 2)
                if (bflag[m] == (int)Enum2.PERIODIC || bflag[m + 1] == (int)Enum2.PERIODIC)
                {
                    if (bflag[m] != (int)Enum2.PERIODIC || bflag[m + 1] != (int)Enum2.PERIODIC)
                        sparta.error.all( "Both sides of boundary must be periodic");
                }
        }
        //public int periodic(int*);
        //public void boundary_modify(int, char**);
        //public virtual int collide(Particle::OnePart*&, int, int, double*, double &,
        //             Particle::OnePart*&);
        public virtual void uncollide(int face, double[] x)
        {
            switch ((Enum1)face)
            {
                case Enum1.XLO:
                    x[0] = boxlo[0];
                    break;
                case Enum1.XHI:
                    x[0] = boxhi[0];
                    break;
                case Enum1.YLO:
                    x[1] = boxlo[1];
                    break;
                case Enum1.YHI:
                    x[1] = boxhi[1];
                    break;
                case Enum1.ZLO:
                    x[2] = boxlo[2];
                    break;
                case Enum1.ZHI:
                    x[2] = boxhi[2];
                    break;
            }
        }
        //public void add_region(int, char**);
        //public void delete_region(int, char**);
        public int find_region(string name)
        {
            for (int iregion = 0; iregion < nregion; iregion++)
                if (string.Equals(name, regions[iregion].id)) return iregion;
            return -1;
        }
        public void print_box(string str)
        {
            string tmp = string.Format("{0}orthogonal box = ({1:G} {2:G} {3:G}) to ({4:G} {5:G} {6:G})\n",
                        str, boxlo[0], boxlo[1], boxlo[2], boxhi[0], boxhi[1], boxhi[2]);
            if (sparta.comm.me == 0)
            {
                if (sparta.screen!=null)
                {
                    Console.WriteLine(tmp);
                    new StreamWriter(sparta.screen).WriteLine(str);
                }

                if (sparta.logfile != null)
                    new StreamWriter(sparta.logfile).WriteLine(str);
            }
        }
    }
}
