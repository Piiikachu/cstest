using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{

    class FixEmitFace : FixEmit
    {
        enum Enum1{ XLO, XHI, YLO, YHI, ZLO, ZHI, INTERIOR };         // same as Domain
        enum Enum2{ PERIODIC, OUTFLOW, REFLECT, SURFACE, AXISYM };  // same as Domain
        enum Enum3{ UNKNOWN, OUTSIDE, INSIDE, OVERLAP };           // same as Grid
        enum Enum4{ PKEEP, PINSERT, PDONE, PDISCARD, PENTRY, PEXIT, PSURF };   // several files
        enum Enum5{ NCHILD, NPARENT, NUNKNOWN, NPBCHILD, NPBPARENT, NPBUNKNOWN, NBOUND };  // Grid
        enum Enum6{ NOSUBSONIC, PTBOTH, PONLY };

        public const int DELTATASK = 256;
        public const double TEMPLIMIT = 1.0e5;


        //public virtual void init();
        //public virtual void post_compress_grid();

        // one insertion task for a cell and a face

        public struct Task
        {
            public double[] lo;               // lower-left corner of face
            public double[] hi;               // upper-right corner of face
            public double[] normal;           // inward normal from external boundary
            public double area;                // area of cell face
            public double ntarget;             // # of mols to insert for all species
            public double nrho;                // from mixture or adjacent subsonic cell
            public double temp_thermal;        // from mixture or adjacent subsonic cell
            public double temp_rot;            // from mixture or subsonic temp_thermal
            public double temp_vib;            // from mixture or subsonic temp_thermal
            public double[] vstream;          // from mixture or adjacent subsonic cell
            public double[] ntargetsp;          // # of mols to insert for each species,
                                         //   only defined for PERSPECIES
            public double[] vscale;             // vscale for each species,
                                        //   only defined for subsonic_style PONLY

            public int icell;                  // associated cell index, unsplit or split cell
            public int iface;                  // which face of unsplit or split cell
            public int pcell;                  // associated cell index for particles
                                         // unsplit or sub cell (not split cell)
            public int ndim;                   // dim (0,1,2) normal to face
            public int pdim, qdim;              // 2 dims (0,1,2) parallel to face
        };


        public FixEmitFace(SPARTA sparta, int narg, string[] arg) : base(sparta, narg, arg)
        {
            if (narg < 4) sparta.error.all("Illegal fix emit/face command");

            imix = sparta.particle.find_mixture(arg[2]);
            if (imix < 0) sparta.error.all("Fix emit/face mixture ID does not exist");

            // flag specified faces

            faces[(int)Enum1.XLO] = faces[(int)Enum1.XHI] = faces[(int)Enum1.YLO] = faces[(int)Enum1.YHI] =
              faces[(int)Enum1.ZLO] = faces[(int)Enum1.ZHI] = 0;

            int iarg = 3;
            while (iarg < narg)
            {
                if (string.Compare(arg[iarg], "all") == 0)
                {
                    if (sparta.domain.dimension == 3)
                        faces[(int)Enum1.XLO] = faces[(int)Enum1.XHI] = faces[(int)Enum1.YLO] = faces[(int)Enum1.YHI] =
                          faces[(int)Enum1.ZLO] = faces[(int)Enum1.ZHI] = 1;
                    else faces[(int)Enum1.XLO] = faces[(int)Enum1.XHI] = faces[(int)Enum1.YLO] = faces[(int)Enum1.YHI] = 1;
                }
                else if (string.Compare(arg[iarg], "xlo") == 0) faces[(int)Enum1.XLO] = 1;
                else if (string.Compare(arg[iarg], "xhi") == 0) faces[(int)Enum1.XHI] = 1;
                else if (string.Compare(arg[iarg], "ylo") == 0) faces[(int)Enum1.YLO] = 1;
                else if (string.Compare(arg[iarg], "yhi") == 0) faces[(int)Enum1.YHI] = 1;
                else if (string.Compare(arg[iarg], "zlo") == 0) faces[(int)Enum1.ZLO] = 1;
                else if (string.Compare(arg[iarg], "zhi") == 0) faces[(int)Enum1.ZHI] = 1;
                else break;
                iarg++;
            }

            // optional args

            np = 0;
            subsonic = 0;
            subsonic_style = (int)Enum6.NOSUBSONIC;
            subsonic_warning = 0;
            twopass = 0;
            string[] argss = new string[narg - iarg];
            Array.Copy(arg, iarg, argss, 0, narg - iarg);
            options(narg - iarg, argss);

            // error checks

            if (sparta.domain.dimension == 2 && (faces[(int)Enum1.ZLO] !=0|| faces[(int)Enum1.ZHI]!=0))
                sparta.error.all("Cannot use fix emit/face in z dimension for 2d simulation");
            if (sparta.domain.axisymmetric !=0&& faces[(int)Enum1.YLO]!=0)
                sparta.error.all("Cannot use fix emit/face on ylo face for axisymmetric model");
            if (np > 0 && perspecies==0)
                sparta.error.all("Cannot use fix emit/face n > 0 with perspecies yes");
            if (np > 0 && subsonic==0)
                sparta.error.all("Cannot use fix emit/face n > 0 with subsonic");

            // task list and subsonic data structs

            tasks = null;
            ntask = ntaskmax = 0;

            maxactive = 0;
            activecell = null;
        }

        int imix, np, subsonic, subsonic_style, subsonic_warning;
        int[] faces=new int[6];
        int npertask, nthresh, twopass;
        double psubsonic, tsubsonic, nsubsonic;
        double tprefactor, soundspeed_mixture;

        // copies of data from other classes

        int dimension, nspecies;
        double fnum, dt;
        double[] fraction,cummulative;

        List<Surf.Point> pts;
        List<Surf.Line> lines;
        List<Surf.Tri> tris;

        // ntask = # of tasks is stored by parent class
        List<Task> tasks;           // list of particle insertion tasks
        int ntaskmax;          // max # of tasks allocated

        // active grid cells assigned to tasks, used by subsonic sorting

        int maxactive;
        int[] activecell;

        // protected methods

        //protected int create_task(int);
        //protected virtual void perform_task();
        //protected void perform_task_onepass();
        //protected virtual void perform_task_twopass();
         
        //protected int split(int, int);
         
        //protected void subsonic_inflow();
        //protected void subsonic_sort();
        //protected void subsonic_grid();
         
        //protected virtual int pack_task(int, char*, int);
        //protected virtual int unpack_task(char*, int);
        //protected virtual void copy_task(int, int, int, int);
        //protected virtual void grow_task();
        //protected virtual void realloc_nspecies();
         
        //protected int option(int, char**);



    }
}
