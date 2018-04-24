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


        public new  void init()
        {
            // copies of class data before invoking parent init() and count_task()

            dimension =sparta.domain.dimension;
            fnum = sparta.update.fnum;
            dt = sparta.update.dt;

            nspecies = sparta.particle.mixture[imix].nspecies;
            fraction = sparta.particle.mixture[imix].fraction;
            cummulative = sparta.particle.mixture[imix].cummulative;

            pts = sparta.surf.pts;
            lines = sparta.surf.lines;
            tris = sparta.surf.tris;

            // subsonic prefactor

            tprefactor = sparta.update.mvv2e / (3.0 * sparta.update.boltz);

            // mixture soundspeed, used by subsonic PONLY as default cell property

            double avegamma = 0.0;
            double avemass = 0.0;

            for (int m = 0; m < nspecies; m++)
            {
                int ispecies = sparta.particle.mixture[imix].species[m];
                avemass += fraction[m] * sparta.particle.species[ispecies].mass;
                avegamma += fraction[m] * (1.0 + 2.0 /
                                           (3.0 + sparta.particle.species[ispecies].rotdof));
            }

            soundspeed_mixture = Math.Sqrt(avegamma * sparta.update.boltz *
                                      sparta.particle.mixture[imix].temp_thermal / avemass);

            // cannot inflow thru periodic boundary

            for (int i = 0; i < 6; i++)
                if (faces[i]!=0 &&sparta.domain.bflag[i] == (int) Enum2.PERIODIC)
                    sparta.error.all("Cannot use fix emit/face on periodic boundary");

            // cannot have inflow on yhi if axisymmetric

            double[] vstream = sparta.particle.mixture[imix].vstream;

            if (sparta.domain.axisymmetric !=0 && faces[(int)Enum1.YHI]!=0 && vstream[1] != 0.0)
                sparta.error.all("Cannot use fix emit on axisymmetric yhi if streaming velocity has a y-component");

            // warn if any inflow face does not have an inward normal
            //   in direction of streaming velocity

            double[] normal=new double[3];
            int flag = 0;

            for (int i = 0; i < 6; i++)
            {
                if (faces[i]==0) continue;
                normal[0] = normal[1] = normal[2] = 0.0;
                if (i % 2 == 0) normal[i / 2] = 1.0;
                else normal[i / 2] = -1.0;
                double indot = vstream[0] * normal[0] + vstream[1] * normal[1] +
                  vstream[2] * normal[2];
                if (indot < 0.0) flag = 1;
            }

            if (flag!=0 && sparta.comm.me == 0)
                sparta.error.all(
                               "warning One or more fix inflow faces oppose streaming velocity");

            // if used, reallocate ntargetsp and vscale for each task
            // b/c nspecies count of mixture may have changed

            realloc_nspecies();

            // invoke FixEmit::init() to populate task list
            // it calls create_task() for each grid cell

            ntask = 0;
            base.init();

            // if Np > 0, nper = # of insertions per task
            // set nthresh so as to achieve exactly Np insertions
            // tasks > tasks_with_no_extra need to insert 1 extra particle
            // NOTE: currently setting same # of insertions per task
            //       could instead weight by cell face area

            if (np > 0)
            {
                int all=0, nupto=0, tasks_with_no_extra;
                sparta.mpi.MPI_Allreduce(ref ntask, ref all, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
                if (all!=0)
                {
                    npertask = np / all;
                    tasks_with_no_extra = all - (np % all);
                }
                else npertask = tasks_with_no_extra = 0;

                sparta.mpi.MPI_Scan(ref ntask,ref nupto, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
                if (tasks_with_no_extra < nupto - ntask) nthresh = 0;
                else if (tasks_with_no_extra >= nupto) nthresh = ntask;
                else nthresh = tasks_with_no_extra - (nupto - ntask);
            }
        }
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

        private SPARTA sparta;
        public FixEmitFace(SPARTA sparta, int narg, string[] arg) : base(sparta, narg, arg)
        {
            this.sparta = sparta;
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

        protected override int create_task(int icell)
        {
            int i, j, n, flag, isp, extflag;
            Enum1 iface;
            int[] cflags;
            double indot, area, ntargetsp;

            Grid.ChildCell[] cells = sparta.grid.cells;
            Grid.ChildInfo[] cinfo = sparta.grid.cinfo;

            double nrho = sparta.particle.mixture[imix].nrho;
            double[] vstream = sparta.particle.mixture[imix].vstream;
            double[] vscale = sparta.particle.mixture[imix].vscale;

            // corners[i][j] = J corner points of face I of a grid cell
            // works for 2d quads and 3d hexes

            int[,] corners = new int[6, 4]{{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7},
               {0,1,2,3}, {4,5,6,7}};
            int nface_pts = 4;
            if (sparta.domain.dimension == 2) nface_pts = 2;

            // loop over 6 faces of icell

            int ntaskorig = ntask;
            int nmask = cells[icell].nmask;
            for ( i = 0; i < 6; i++)
            {
                if (i == 0) iface = Enum1.XLO;
                else if (i == 1) iface = Enum1.XHI;
                else if (i == 2) iface = Enum1.YLO;
                else if (i == 3) iface = Enum1.YHI;
                else if (i == 4) iface = Enum1.ZLO;
                else  iface = Enum1.ZHI;

                // flag = 1 if insertion happens on iface of cell
                // only if face adjoins global boundary with inflow defined
                // if cell is OVERLAP:
                //   allow if any face corner point is OUTSIDE and none is INSIDE
                //   disallow if any pt of any line/tri in cell touches face

                flag = 0;
                if (faces[(int)iface]!=0 && sparta.grid.neigh_decode(nmask, iface) == (int)Enum5.NBOUND)
                {
                    if (cinfo[icell].type == (int)Enum3.OUTSIDE) flag = 1;
                    else if (cinfo[icell].type == (int)Enum3.OVERLAP)
                    {
                        flag = 1;
                        cflags = cinfo[icell].corner;

                        extflag = 0;
                        for (j = 0; j < nface_pts; j++)
                        {
                            if (cflags[corners[(int)iface,j]] == (int)Enum3.OUTSIDE) extflag = 1;
                            else if (cflags[corners[(int)iface,j]] == (int)Enum3.INSIDE) flag = 0;
                        }
                        if (extflag==0) flag = 0;

                        if (flag!=0 && dimension == 2)
                        {
                            for (j = 0; j < cells[icell].nsurf; j++)
                            {
                                n = cells[icell].csurfs[j];
                                if (Geometry.line_quad_face_touch(pts[lines[n].p1].x,
                                                         pts[lines[n].p2].x,
                                                         (int)iface, cells[icell].lo, cells[icell].hi)!=0)
                                {
                                    flag = 0;
                                    break;
                                }
                            }
                        }
                        else if (flag!=0 && dimension == 3)
                        {
                            for (j = 0; j < cells[icell].nsurf; j++)
                            {
                                n = cells[icell].csurfs[j];
                                if (Geometry.tri_hex_face_touch(pts[tris[n].p1].x,
                                                       pts[tris[n].p2].x,
                                                       pts[tris[n].p3].x,
                                                       (int)iface, cells[icell].lo, cells[icell].hi)!=0)
                                {
                                    flag = 0;
                                    break;
                                }
                            }
                        }
                    }
                }


                // no insertions on this face

                if (flag==0) continue;

                // set cell parameters of task
                // pcell = sub cell for particles if a split cell

                if (ntask == ntaskmax) grow_task();
                Task task = tasks[ntask];
                task.icell = icell;
                task.iface = (int)iface;
                if (cells[icell].nsplit > 1) task.pcell = split(icell, iface);
                else task.pcell = icell;

                // set face-dependent params of task
                task.lo = new double[3];
                task.hi = new double[3];
                task.lo[0] = cells[icell].lo[0];
                task.hi[0] = cells[icell].hi[0];
                task.lo[1] = cells[icell].lo[1];
                task.hi[1] = cells[icell].hi[1];
                task.lo[2] = cells[icell].lo[2];
                task.hi[2] = cells[icell].hi[2];
                if (dimension == 2) task.lo[2] = task.hi[2] = 0.0;
                task.normal[0] = 0.0;
                task.normal[1] = 0.0;
                task.normal[2] = 0.0;

                if (iface == XLO || iface == XHI)
                {
                    task.ndim = 0;
                    task.pdim = 1;
                    task.qdim = 2;
                    if (iface == XLO) tasks[ntask].hi[0] = cells[icell].lo[0];
                    else tasks[ntask].lo[0] = cells[icell].hi[0];
                    if (iface == XLO) tasks[ntask].normal[0] = 1.0;
                    else tasks[ntask].normal[0] = -1.0;
                }
                else if (iface == YLO || iface == YHI)
                {
                    task.ndim = 1;
                    task.pdim = 0;
                    task.qdim = 2;
                    if (iface == YLO) tasks[ntask].hi[1] = cells[icell].lo[1];
                    else tasks[ntask].lo[1] = cells[icell].hi[1];
                    if (iface == YLO) tasks[ntask].normal[1] = 1.0;
                    else tasks[ntask].normal[1] = -1.0;
                }
                else if (iface == ZLO || iface == ZHI)
                {
                    task.ndim = 2;
                    task.pdim = 0;
                    task.qdim = 1;
                    if (iface == ZLO) tasks[ntask].hi[2] = cells[icell].lo[2];
                    else tasks[ntask].lo[2] = cells[icell].hi[2];
                    if (iface == ZLO) tasks[ntask].normal[2] = 1.0;
                    else tasks[ntask].normal[2] = -1.0;
                }

                // indot = dot product of vstream with inward face normal

                indot = vstream[0] * tasks[ntask].normal[0] +
                  vstream[1] * tasks[ntask].normal[1] +
                  vstream[2] * tasks[ntask].normal[2];

                // area = area for insertion
                // depends on dimension and axisymmetry

                if (iface == XLO || iface == XHI)
                {
                    if (dimension == 3)
                        area = (cells[icell].hi[1] - cells[icell].lo[1]) *
                          (cells[icell].hi[2] - cells[icell].lo[2]);
                    else if (domain.axisymmetric)
                        area = (cells[icell].hi[1] * cells[icell].hi[1] -
                                cells[icell].lo[1] * cells[icell].lo[1]) * MY_PI;
                    else area = cells[icell].hi[1] - cells[icell].lo[1];
                }
                else if (iface == YLO || iface == YHI)
                {
                    if (dimension == 3)
                        area = (cells[icell].hi[0] - cells[icell].lo[0]) *
                          (cells[icell].hi[2] - cells[icell].lo[2]);
                    else if (domain.axisymmetric)
                        area = 2.0 * MY_PI * cells[icell].hi[1] *
                          (cells[icell].hi[0] - cells[icell].lo[0]);
                    else area = cells[icell].hi[0] - cells[icell].lo[0];
                }
                else if (iface == ZLO || iface == ZHI)
                {
                    area = (cells[icell].hi[0] - cells[icell].lo[0]) *
                  (cells[icell].hi[1] - cells[icell].lo[1]);
                }
                task.area = area;

                // set ntarget and ntargetsp via mol_inflow()
                // skip task if final ntarget = 0.0, due to large outbound vstream
                // do not skip for subsonic since it resets ntarget every step

                task.ntarget = 0.0;
                for (isp = 0; isp < nspecies; isp++)
                {
                    ntargetsp = mol_inflow(indot, vscale[isp], fraction[isp]);
                    ntargetsp *= nrho * area * dt / fnum;
                    ntargetsp /= cinfo[icell].weight;
                    task.ntarget += ntargetsp;
                    if (perspecies) tasks[ntask].ntargetsp[isp] = ntargetsp;
                }

                if (subsonic==0)
                {
                    if (tasks[ntask].ntarget == 0.0) continue;
                    if (tasks[ntask].ntarget >= MAXSMALLINT)
                        error.one(FLERR,
                                   "Fix emit/face insertion count exceeds 32-bit int");
                }

                // initialize other task values with mixture properties
                // may be overwritten by subsonic methods

                task.nrho = sparta.particle.mixture[imix].nrho;
                task.temp_thermal = sparta.particle.mixture[imix].temp_thermal;
                task.temp_rot = sparta.particle.mixture[imix].temp_rot;
                task.temp_vib = sparta.particle.mixture[imix].temp_vib;
                task.vstream[0] = sparta.particle.mixture[imix].vstream[0];
                task.vstream[1] = sparta.particle.mixture[imix].vstream[1];
                task.vstream[2] = sparta.particle.mixture[imix].vstream[2];

                // increment task counter
                tasks[ntask] = task;
                ntask++;
            }

            




        }
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
        protected virtual void realloc_nspecies()
        {
            if (perspecies!=0)
            {
                for (int i = 0; i < ntask; i++)
                {
                    Task task = tasks[i];
                    task.ntargetsp = new double[nspecies];
                    tasks[i] = task;
                }
            }
            if (subsonic_style == (int) Enum6.PONLY)
            {
                for (int i = 0; i < ntask; i++)
                {
                    Task task = tasks[i];
                    task.vscale = new double[nspecies];
                    tasks[i] = task;
                }
            }
        }

        //protected int option(int, char**);



    }
}
