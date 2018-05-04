using System;
using System.Collections.Generic;
using bigint = System.Int64;
using cellint = System.Int32;

namespace cstest
{
    public delegate void MovePerTube(double dt, double[] x, double[] v);
    public delegate void MovePtr(int a, int b);
    public class Update
    {
        enum Enum1 { XLO, XHI, YLO, YHI, ZLO, ZHI, INTERIOR };         // same as Domain
        enum Enum2 { PERIODIC, OUTFLOW, REFLECT, SURFACE, AXISYM };  // same as Domain
        enum Enum3 { OUTSIDE, INSIDE, ONSURF2OUT, ONSURF2IN };      // several files
        enum Enum4 : int { PKEEP, PINSERT, PDONE, PDISCARD, PENTRY, PEXIT, PSURF };   // several files
        enum Enum5 { NCHILD, NPARENT, NUNKNOWN, NPBCHILD, NPBPARENT, NPBUNKNOWN, NBOUND };  // Grid
        enum Enum6 { TALLYAUTO, TALLYREDUCE, TALLYLOCAL };         // same as Surf

        public const int MAXSTUCK = 20;
        public const double EPSPARAM = 1.0e-7;

        // either set ID or PROC/INDEX, set other to -1

        //#define MOVE_DEBUG 1              // un-comment to debug one particle
        public const int MOVE_DEBUG_ID = 537898343;   // particle ID
        public const int MOVE_DEBUG_PROC = 0;        // owning proc
        public const int MOVE_DEBUG_INDEX = 16617;       // particle index on owning proc
        public const int MOVE_DEBUG_STEP = 16;    // timestep


        public bigint ntimestep;               // current timestep
        public int nsteps;                     // # of steps to run
        public int runflag;                    // 0 for unset, 1 for run
        public bigint firststep, laststep;      // 1st & last step of this run
        public bigint beginstep, endstep;       // 1st and last step of multiple runs
        public int first_update;               // 0 before initial update, 1 after
        public double dt;                      // timestep size

        public string unit_style;      // style of units used throughout simulation
        public double boltz;          // Boltzmann constant (eng/degree K)
        public double mvv2e;          // conversion of mv^2 to energy

        public double fnum;           // ratio of real particles to simulation particles
        public double nrho;           // number density of background gas
        public double[] vstream = new double[3];     // streaming velocity of background gas
        public double temp_thermal;   // thermal temperature of background gas
        public double[] gravity = new double[3];     // acceleration vector of gravity

        public int nmigrate;          // # of particles to migrate to new procs
        public int[] mlist;            // indices of particles to migrate

        // current step counters
        public int niterate;          // iterations of move/comm
        public int ntouch_one;        // particle-cell touches
        public int ncomm_one;         // particles migrating to new procs
        public int nboundary_one;     // particles colliding with global boundary
        public int nexit_one;         // particles exiting outflow boundary
        public int nscheck_one;       // surface elements checked for collisions
        public int nscollide_one;     // particle/surface collisions

        public bigint first_running_step; // timestep running counts start on
        public int niterate_running;      // running count of move/comm interations
        public bigint nmove_running;      // running count of total particle moves
        public bigint ntouch_running;     // running count of current step counters
        public bigint ncomm_running;
        public bigint nboundary_running;
        public bigint nexit_running;
        public bigint nscheck_running;
        public bigint nscollide_running;

        public MovePerTube moveperturb;
        //public MovePtr moveptr;
        int[] moveptr ;

        public int nstuck;                // # of particles stuck on surfs and deleted

        public int reorder_period;        // # of timesteps between particle reordering

        public int copymode;          // 1 if copy of class (prevents deallocation of
                                      //  base class when child copy is destroyed)

        public RanMars ranmaster;   // master random number generator

        double[] rcblo, rcbhi = new double[3];    // debug info from RCB for dump image

        private SPARTA sparta;

        public Update(SPARTA sparta)
        {
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);
            sparta.mpi.MPI_Comm_size(sparta.world, ref nprocs);

            ntimestep = 0;
            firststep = laststep = 0;
            beginstep = endstep = 0;
            runflag = 0;

            unit_style = null;
            set_units("si");

            fnum = 1.0;
            nrho = 1.0;
            vstream[0] = vstream[1] = vstream[2] = 0.0;
            temp_thermal = 273.15;
            gravity[0] = gravity[1] = gravity[2] = 0.0;

            maxmigrate = 0;
            mlist = null;

            nslist_compute = nblist_compute = 0;
            slist_compute = blist_compute = null;
            slist_active = blist_active = null;

            ranmaster = new RanMars(sparta);

            reorder_period = 0;

            copymode = 0;
        }
        public void set_units(string style)
        {
            // physical constants from:
            // http://physics.nist.gov/cuu/Constants/Table/allascii.txt

            if (string.Compare(style, "cgs") == 0)
            {
                boltz = 1.3806488e-16;
                mvv2e = 1.0;
                dt = 1.0;

            }
            else if (string.Compare(style, "si") == 0)
            {
                boltz = 1.3806488e-23;
                mvv2e = 1.0;
                dt = 1.0;

            }
            else sparta.error.all("Illegal units command");

            //delete[] unit_style;
            int n = style.Length + 1;
            unit_style = string.Copy(style);
        }
        public virtual void init()
        {
            // init the Update class if performing a run, else just return
            // only set first_update if a run is being performed

            if (runflag == 0) return;
            first_update = 1;

            // choose the appropriate move method

            if (sparta.domain.dimension == 3)
            {
                if (sparta.surf.exist != 0) moveptr =new int[]{3,1};
                else moveptr = new int[] { 3, 0 };
            }
            else if (sparta.domain.axisymmetric!=0)
            {
                if (sparta.surf.exist != 0) moveptr = new int[] { 1, 1 };
                else moveptr = new int[] { 1, 0 };
            }
            else if (sparta.domain.dimension == 2)
            {
                if (sparta.surf.exist!=0) moveptr = new int[] { 2, 1 };
                else moveptr = new int[] { 2, 0 };
            }

            // check gravity vector

            if (sparta.domain.dimension == 2 && gravity[2] != 0.0)
                sparta.error.all("Gravity in z not allowed for 2d");
            if (sparta.domain.axisymmetric != 0 && gravity[1] != 0.0)
                sparta.error.all("Gravity in y not allowed for axi-symmetric model");

            // moveperturb method is set if particle motion is perturbed

            moveperturb = null;
            if (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0)
            {
                if (sparta.domain.dimension == 3) moveperturb =new MovePerTube( gravity3d);
                if (sparta.domain.dimension == 2) moveperturb =new MovePerTube( gravity2d);
            }

            if (moveperturb != null) perturbflag = 1;
            else perturbflag = 0;
        }
        public virtual void setup()
        {
            // initialize counters in case stats outputs them
            // initialize running stats before each run

            ntouch_one = ncomm_one = 0;
            nboundary_one = nexit_one = 0;
            nscheck_one = nscollide_one = 0;
            sparta.surf.nreact_one = 0;

            first_running_step = sparta.update.ntimestep;
            niterate_running = 0;
            nmove_running = ntouch_running = ncomm_running = 0;
            nboundary_running = nexit_running = 0;
            nscheck_running = nscollide_running = 0;
            sparta.surf.nreact_running = 0;
            nstuck = 0;

            collide_react = collide_react_setup();
            bounce_tally = bounce_setup();

            sparta.modify.setup();
            sparta.output.setup(1);
        }
        public virtual void run(int nsteps)
        {
            int n_start_of_step = sparta.modify.n_start_of_step;
            int n_end_of_step = sparta.modify.n_end_of_step;
            //int dynamic = 0;

            // cellweightflag = 1 if grid-based particle weighting is ON

            int cellweightflag = 0;
            if (sparta.grid.cellweightflag!=0) cellweightflag = 1;

            // loop over timesteps

            for (int i = 0; i < nsteps; i++)
            {

                ntimestep++;
                if (collide_react!=0) collide_react_update();
                if (bounce_tally!=0) bounce_set(ntimestep);

                sparta.timer.stamp();

                // start of step fixes

                if (n_start_of_step!=0)
                {
                    sparta.modify.start_of_step();
                    sparta.timer.stamp((int)Timer.Enum1.TIME_MODIFY);
                }

                //if (dynamic) domain.dynamic();

                // move particles

                if (cellweightflag != 0) sparta.particle.pre_weight();
                //moveptr();
                move(moveptr[0], moveptr[1]);
                sparta.timer.stamp((int)Timer.Enum1.TIME_MOVE);

                // communicate particles

                sparta.comm.migrate_particles(nmigrate, mlist);
                if (cellweightflag != 0) sparta.particle.post_weight();
                sparta.timer.stamp((int)Timer.Enum1.TIME_COMM);

                if (sparta.collide != null)
                {
                    sparta.particle.sort();
                    sparta.timer.stamp((int)Timer.Enum1.TIME_SORT);
                    sparta.collide.collisions();
                    
                    sparta.timer.stamp((int)Timer.Enum1.TIME_COLLIDE);
                }

                // diagnostic fixes

                if (n_end_of_step!=0)
                {
                    sparta.modify.end_of_step();
                    
                    sparta.timer.stamp((int)Timer.Enum1.TIME_MODIFY);
                }

                // all output

                if (ntimestep == sparta.output.next)
                {
                    sparta.output.write(ntimestep);
                    
                    sparta.timer.stamp((int)Timer.Enum1.TIME_OUTPUT);
                }
            }
        }
        public void global(int nnarg, string[] arg)
        {
            int narg = nnarg + 1;
            if (narg < 1) sparta.error.all("Illegal global command");
            int iarg = 1;
            while (iarg < narg)
            {
                if (string.Equals(arg[iarg], "fnum"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    fnum = System.Double.Parse(arg[iarg + 1]);
                    if (fnum <= 0.0) sparta.error.all("Illegal global command");
                    iarg += 2;
                }
                else if (string.Equals(arg[iarg], "nrho"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    nrho = System.Double.Parse(arg[iarg + 1]);
                    if (nrho <= 0.0) sparta.error.all("Illegal global command");
                    iarg += 2;
                }
                else if (string.Equals(arg[iarg], "vstream"))
                {
                    if (iarg + 4 > narg) sparta.error.all("Illegal global command");
                    vstream[0] = System.Double.Parse(arg[iarg + 1]);
                    vstream[1] = System.Double.Parse(arg[iarg + 2]);
                    vstream[2] = System.Double.Parse(arg[iarg + 3]);
                    iarg += 4;
                }
                else if (string.Equals(arg[iarg], "temp"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    temp_thermal = System.Double.Parse(arg[iarg + 1]);
                    if (temp_thermal <= 0.0) sparta.error.all("Illegal global command");
                    iarg += 2;
                }
                else if (string.Equals(arg[iarg], "gravity"))
                {
                    if (iarg + 5 > narg) sparta.error.all("Illegal global command");
                    double gmag = System.Double.Parse(arg[iarg + 1]);
                    gravity[0] = System.Double.Parse(arg[iarg + 2]);
                    gravity[1] = System.Double.Parse(arg[iarg + 3]);
                    gravity[2] = System.Double.Parse(arg[iarg + 4]);
                    if (gmag < 0.0) sparta.error.all("Illegal global command");
                    if (gmag > 0.0 &&
                        gravity[0] == 0.0 && gravity[1] == 0.0 && gravity[2] == 0.0)
                        sparta.error.all("Illegal global command");
                    if (gmag > 0.0) MathExtra.snorm3(gmag, gravity);
                    else gravity[0] = gravity[1] = gravity[2] = 0.0;
                    iarg += 5;

                }
                else if (string.Equals(arg[iarg], "surfmax"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    if (sparta.surf.exist != 0)
                        sparta.error.all(
                                   "Cannot set global surfmax when surfaces already exist");
                    sparta.grid.maxsurfpercell = System.Int32.Parse(arg[iarg + 1]);
                    if (sparta.grid.maxsurfpercell <= 0) sparta.error.all("Illegal global command");
                    // reallocate paged data structs for variable-length surf into
                    sparta.grid.allocate_surf_arrays();
                    iarg += 2;
                }
                else if (string.Equals(arg[iarg], "surfpush"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    if (string.Equals(arg[iarg + 1], "no"))
                    {
                        sparta.surf.pushflag = 0;
                        iarg += 2;
                    }
                    else if (string.Equals(arg[iarg + 1], "yes"))
                    {
                        sparta.surf.pushflag = 1;
                        iarg += 2;
                    }
                    else
                    {
                        if (iarg + 4 > narg) sparta.error.all("Illegal global command");
                        sparta.surf.pushflag = 2;
                        sparta.surf.pushlo = sparta.input.numeric(arg[iarg + 1]);
                        sparta.surf.pushhi = sparta.input.numeric(arg[iarg + 2]);
                        sparta.surf.pushvalue = sparta.input.numeric(arg[iarg + 3]);
                        if (sparta.surf.pushlo > sparta.surf.pushhi)
                            sparta.error.all("Illegal global command");
                        iarg += 4;
                    }
                }
                else if (string.Equals(arg[iarg], "gridcut"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    sparta.grid.cutoff = double.Parse(arg[iarg + 1]);
                    if (sparta.grid.cutoff < 0.0 && sparta.grid.cutoff != -1.0)
                        sparta.error.all("Illegal global command");
                    // force ghost info to be regenerated with new cutoff
                    sparta.grid.remove_ghosts();
                    iarg += 2;

                }
                else if (string.Equals(arg[iarg], "comm/sort"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    if (string.Equals(arg[iarg + 1], "yes")) sparta.comm.commsortflag = 1;
                    else if (string.Equals(arg[iarg + 1], "no")) sparta.comm.commsortflag = 0;
                    else sparta.error.all("Illegal global command");
                    iarg += 2;
                }
                else if (string.Equals(arg[iarg], "comm/style"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    if (string.Equals(arg[iarg + 1], "neigh")) sparta.comm.commpartstyle = 1;
                    else if (string.Equals(arg[iarg + 1], "all")) sparta.comm.commpartstyle = 0;
                    else sparta.error.all("Illegal global command");
                    iarg += 2;
                }
                else if (string.Equals(arg[iarg], "surf/comm"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    if (string.Equals(arg[iarg + 1], "auto")) sparta.surf.tally_comm = (int)Enum6.TALLYAUTO;
                    else if (string.Equals(arg[iarg + 1], "all")) sparta.surf.tally_comm = (int)Enum6.TALLYREDUCE;
                    else if (string.Equals(arg[iarg + 1], "local")) sparta.surf.tally_comm = (int)Enum6.TALLYLOCAL;
                    else sparta.error.all("Illegal global command");
                    iarg += 2;
                }
                else if (string.Equals(arg[iarg], "weight"))
                {
                    if (iarg + 3 > narg) sparta.error.all("Illegal global command");
                    string[] temp = new string[narg];
                    Array.Copy(arg, iarg + 3, temp, 0, narg);
                    if (string.Equals(arg[iarg + 1], "cell"))
                    {
                        sparta.grid.weight(1, temp);
                        Console.WriteLine("update.global().try copy array");
                    }
                    else sparta.error.all("Illegal global command");
                    iarg += 3;
                }
                else if (string.Equals(arg[iarg], "particle/reordor"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    reorder_period = int.Parse(arg[iarg + 1]);
                    if (reorder_period < 0)
                    {
                        sparta.error.all("Illegal global command");
                    }
                    iarg += 2;
                }
                else
                {
                    sparta.error.all("Illegal global command");
                }



            }



        }
        //      public void reset_timestep(int, char**);

        public int split3d(int icell, double[] x)
        {
            int m, cflag, isurf, hitflag, side, minsurfindex = 0;
            double param, minparam;
            double[] xc = new double[3];
            Surf.Tri tri;

            Grid.ChildCell[] cells = sparta.grid.cells;
            Grid.SplitInfo[] sinfo = sparta.grid.sinfo;
            List<Surf.Tri> tris = sparta.surf.tris;
            List<Surf.Point> pts = sparta.surf.pts;

            // check for collisions with lines in cell
            // find 1st surface hit via minparam
            // only consider tris that are mapped via csplits to a split cell
            //   unmapped tris only touch cell surf at xnew
            //   another mapped tri should include same xnew
            // not considered a collision if particles starts on surf, moving out
            // not considered a collision if 2 params are tied and one is INSIDE surf

            int nsurf = cells[icell].nsurf;
            int[] csurfs = cells[icell].csurfs;
            int isplit = cells[icell].isplit;
            int[] csplits = sinfo[isplit].csplits;
            double[] xnew = sinfo[isplit].xsplit;

            cflag = 0;
            minparam = 2.0;
            for (m = 0; m < nsurf; m++)
            {
                if (csplits[m] < 0) continue;
                isurf = csurfs[m];
                tri = tris[isurf];
                hitflag = (Geometry.line_tri_intersect(x, xnew,
                                     pts[tri.p1].x, pts[tri.p2].x, pts[tri.p3].x,
                                     tri.norm, xc, out param, out side)) ? 0 : 1;

                if (hitflag != 0 && side != (int)Enum3.INSIDE && param < minparam)
                {
                    cflag = 1;
                    minparam = param;
                    minsurfindex = m;
                }
            }
            if (cflag == 0) return sinfo[isplit].csubs[sinfo[isplit].xsub];
            int index = csplits[minsurfindex];
            return sinfo[isplit].csubs[index];

        }
        public int split2d(int icell, double[] x)
        {
            int m, cflag, isurf, hitflag, side, minsurfindex = 0;
            double param, minparam;
            double[] xc = new double[3];
            Surf.Line line;

            Grid.ChildCell[] cells = sparta.grid.cells;
            Grid.SplitInfo[] sinfo = sparta.grid.sinfo;
            List<Surf.Line> lines = sparta.surf.lines;
            List<Surf.Point> pts = sparta.surf.pts;

            // check for collisions with lines in cell
            // find 1st surface hit via minparam
            // only consider lines that are mapped via csplits to a split cell
            //   unmapped lines only touch cell surf at xnew
            //   another mapped line should include same xnew
            // not considered a collision if particle starts on surf, moving out
            // not considered a collision if 2 params are tied and one is INSIDE surf

            int nsurf = cells[icell].nsurf;
            int[] csurfs = cells[icell].csurfs;
            int isplit = cells[icell].isplit;
            int[] csplits = sinfo[isplit].csplits;
            double[] xnew = sinfo[isplit].xsplit;

            cflag = 0;
            minparam = 2.0;
            for (m = 0; m < nsurf; m++)
            {
                if (csplits[m] < 0) continue;
                isurf = csurfs[m];
                line = lines[isurf];
                hitflag = (Geometry.line_line_intersect(x, xnew,
                                      pts[line.p1].x, pts[line.p2].x, line.norm,
                                      xc, out param, out side)) ? 0 : 1; ;

                if (hitflag != 0 && side != (int)Enum3.INSIDE && param < minparam)
                {
                    cflag = 1;
                    minparam = param;
                    minsurfindex = m;
                }
            }

            if (cflag == 0) return sinfo[isplit].csubs[sinfo[isplit].xsub];
            int index = csplits[minsurfindex];
            return sinfo[isplit].csubs[index];


        }

        //protected:
        protected int me, nprocs;
        protected int maxmigrate;            // max # of particles in mlist
        RanPark random;     // RNG for particle timestep moves

        protected int collide_react;         // 1 if any SurfCollide or React classes defined
        protected int nsc, nsr;               // copy of Collide/React data in Surf class
        protected List<SurfCollide> sc;
        protected List<SurfReact> sr;

        int bounce_tally;               // 1 if any bounces are ever tallied
        int nslist_compute;             // # of computes that tally surf bounces
        int nblist_compute;             // # of computes that tally boundary bounces
        List<Compute> slist_compute;  // list of all surf bounce Computes
        List<Compute> blist_compute;  // list of all boundary bounce Computes

        int nsurf_tally;         // # of Cmp tallying surf bounce info this step
        int nboundary_tally;     // # of Cmp tallying boundary bounce info this step
        List<Compute> slist_active;   // list of active surf Computes this step
        List<Compute> blist_active;   // list of active boundary Computes this step

        int surf_pre_tally;       // 1 to log particle stats before surf collide
        int boundary_pre_tally;   // 1 to log particle stats before boundary collide

        int collide_react_setup()
        {
            nsc = sparta.surf.nsc;
            sc = sparta.surf.sc;
            nsr = sparta.surf.nsr;
            sr = sparta.surf.sr;

            if (sc!=null || sr != null) return 1;
            return 0;
        }
        void collide_react_update()
        {
            for (int i = 0; i < nsc; i++) sc[i].tally_update();
            for (int i = 0; i < nsr; i++) sr[i].tally_update();
        }

        int bounce_setup()
        {
            slist_compute = blist_compute = null;
            nslist_compute = nblist_compute = 0;
            for (int i = 0; i < sparta.modify.ncompute; i++)
            {
                if (sparta.modify.compute[i].surf_tally_flag!=0) nslist_compute++;
                if (sparta.modify.compute[i].boundary_tally_flag != 0) nblist_compute++;
            }

            if (nslist_compute!=0) slist_compute = new List<Compute>(nslist_compute);
            if (nblist_compute!=0) blist_compute = new List<Compute>(nblist_compute);
            if (nslist_compute!=0) slist_active =  new List<Compute>(nslist_compute);
            if (nblist_compute != 0) blist_active =new List<Compute>(nblist_compute);

            nslist_compute = nblist_compute = 0;
            for (int i = 0; i < sparta.modify.ncompute; i++)
            {
                if (sparta.modify.compute[i].surf_tally_flag!=0)
                    slist_compute[nslist_compute++] = sparta.modify.compute[i];
                if (sparta.modify.compute[i].boundary_tally_flag != 0)
                    blist_compute[nblist_compute++] = sparta.modify.compute[i];
            }

            if (nslist_compute != 0 || nblist_compute != 0) return 1;
            nsurf_tally = nboundary_tally = 0;
            return 0;
        }
        void bounce_set(bigint ntimestep)
        {
            int i;

            nsurf_tally = 0;
            if (nslist_compute!=0)
            {
                for (i = 0; i < nslist_compute; i++)
                    if (slist_compute[i].matchstep(ntimestep)!=0)
                    {
                        slist_active[nsurf_tally++] = slist_compute[i];
                        slist_compute[i].clear();
                    }
            }

            nboundary_tally = 0;
            if (nblist_compute != 0)
            {
                for (i = 0; i < nblist_compute; i++)
                    if (blist_compute[i].matchstep(ntimestep)!=0)
                    {
                        blist_active[nboundary_tally++] = blist_compute[i];
                        blist_compute[i].clear();
                    }
            }
        }
        //      void reset_timestep(bigint);

        //      //int axi_vertical_line(double, double *, double *, double, double, double,
        //      //                     double &);

        //      // remap x and v components into axisymmetric plane
        //      // input x at end of linear move (x = xold + dt*v)
        //      // change x[1] = Math.Sqrt(x[1]^2 + x[2]^2), x[2] = 0.0
        //      // change vy,vz by rotation into axisymmetric plane

        //      inline void axi_remap(double* x, double* v)
        //      {
        //          double ynew = x[1];
        //          double znew = x[2];
        //          x[1] = Math.Sqrt(ynew * ynew + znew * znew);
        //          x[2] = 0.0;
        //          double rn = ynew / x[1];
        //          double wn = znew / x[1];
        //          double vy = v[1];
        //          double vz = v[2];
        //          v[1] = vy * rn + vz * wn;
        //          v[2] = -vy * wn + vz * rn;
        //      };

        //      typedef void (Update::* FnPtr) ();
        //FnPtr moveptr;             // ptr to move method
        void move(int DIM, int SURF)
        {
            bool hitflag = false;
            int m, icell, icell_original, nmask, outface, bflag, itmp;
            int side = 0, minside = 0, minsurf = 0, nsurf = 0, cflag, isurf, exclude, stuck_iterate;
            int pstart = 0, pstop = 0, entryexit, any_entryexit=0;
            Enum4 pflag;
            Enum5 nflag;
            int[] csurfs;
            cellint[] neigh;
            double dtremain = 0, frac, newfrac, param = 0, minparam, rnew, dtsurf = 0, tc, tmp;
            double[] xnew = new double[3], xhold = new double[3], xc = new double[3],
                vc = new double[3], minxc = new double[3], minvc = new double[3];
            double[] x, v, lo, hi;
            Surf.Tri tri;
            Surf.Line line;
            Particle.OnePart iorig=new Particle.OnePart();
            Particle.OnePart[] particles;
            Particle.OnePart? ipart=null;
            Particle.OnePart? jpart=null;

            // for 2d and axisymmetry only
            // xnew,xc passed to geometry routines which use or set z component

            if (DIM < 3) xnew[2] = xc[2] = 0.0;

            // extend migration list if necessary

            int nlocal = sparta.particle.nlocal;
            int maxlocal = sparta.particle.maxlocal;

            if (nlocal > maxmigrate)
            {
                maxmigrate = maxlocal;
                mlist = new int[maxmigrate];
            }

            // counters

            niterate = 0;
            ntouch_one = ncomm_one = 0;
            nboundary_one = nexit_one = 0;
            nscheck_one = nscollide_one = 0;
            sparta.surf.nreact_one = 0;

            // move/migrate iterations

            Grid.ChildCell[] cells = sparta.grid.cells;
            List<Surf.Tri> tris = sparta.surf.tris;
            List<Surf.Line> lines = sparta.surf.lines;
            List<Surf.Point> pts = sparta.surf.pts;
            double dt = this.dt;
            int notfirst = 0;

            while (true)
            {
                niterate++;
                particles = sparta.particle.particles;
                //Console.WriteLine(nlocal);
                nmigrate = 0;
                entryexit = 0;

                if (notfirst == 0)
                {
                    notfirst = 1;
                    pstart = 0;
                    pstop = nlocal;
                }

                for (int i = pstart; i < pstop; i++)
                {
                    pflag = (Enum4)particles[i].flag;

                    // received from another proc and move is done
                    // if first iteration, PDONE is from a previous step,
                    //   set pflag to PKEEP so move the particle on this step
                    // else do nothing

                    if (pflag == Enum4.PDONE)
                    {
                        pflag = Enum4.PKEEP;
                        particles[i].flag = (int)Enum4.PKEEP;
                        if (niterate > 1) continue;
                    }

                    x = particles[i].x;
                    v = particles[i].v;
                    exclude = -1;

                    // apply moveperturb() to PKEEP and PINSERT since are computing xnew
                    // not to PENTRY,PEXIT since are just re-computing xnew of sender
                    // set xnew[2] to linear move for axisymmetry, will be remapped later
                    // let pflag = PEXIT persist to check during axisymmetric cell crossing

                    switch (pflag)
                    {
                        case Enum4.PKEEP:
                            dtremain = dt;
                            xnew[0] = x[0] + dtremain * v[0];
                            xnew[1] = x[1] + dtremain * v[1];
                            if (DIM != 2) xnew[2] = x[2] + dtremain * v[2];
                            if (perturbflag != 0) moveperturb(dtremain, xnew, v);
                            break;
                        case Enum4.PINSERT:
                            dtremain = particles[i].dtremain;
                            xnew[0] = x[0] + dtremain * v[0];
                            xnew[1] = x[1] + dtremain * v[1];
                            if (DIM != 2) xnew[2] = x[2] + dtremain * v[2];
                            if (perturbflag != 0) moveperturb(dtremain, xnew, v);
                            break;
                        case Enum4.PDONE:

                            break;
                        case Enum4.PDISCARD:
                            break;
                        case Enum4.PENTRY:
                            icell = particles[i].icell;
                            if (cells[icell].nsplit > 1)
                            {
                                if (DIM == 3 && SURF != 0) icell = split3d(icell, x);
                                if (DIM < 3 && SURF != 0) icell = split2d(icell, x);
                                particles[i].icell = icell;
                            }
                            dtremain = particles[i].dtremain;
                            xnew[0] = x[0] + dtremain * v[0];
                            xnew[1] = x[1] + dtremain * v[1];
                            if (DIM != 2) xnew[2] = x[2] + dtremain * v[2];
                            break;
                        case Enum4.PEXIT:
                            dtremain = particles[i].dtremain;
                            xnew[0] = x[0] + dtremain * v[0];
                            xnew[1] = x[1] + dtremain * v[1];
                            if (DIM != 2) xnew[2] = x[2] + dtremain * v[2];
                            break;
                        case Enum4.PSURF:
                            //dtremain = particles[i].dtremain;
                            //xnew[0] = x[0] + dtremain * v[0];
                            //xnew[1] = x[1] + dtremain * v[1];
                            //if (DIM != 2) xnew[2] = x[2] + dtremain * v[2];
                            //if (pflag > Enum4.PSURF) exclude = pflag - Enum4.PSURF - 1;
                            //break;
                        default:
                            break;
                    }
                    if ((int)pflag >= (int)Enum4.PSURF)
                    {
                        dtremain = particles[i].dtremain;
                        xnew[0] = x[0] + dtremain * v[0];
                        xnew[1] = x[1] + dtremain * v[1];
                        if (DIM != 2) xnew[2] = x[2] + dtremain * v[2];
                        if (pflag > Enum4.PSURF) exclude = pflag - Enum4.PSURF - 1;

                    }

                    particles[i].flag = (int)Enum4.PKEEP;
                    icell = particles[i].icell;
                    lo = cells[icell].lo;
                    hi = cells[icell].hi;
                    neigh = cells[icell].neigh;
                    nmask = cells[icell].nmask;
                    stuck_iterate = 0;
                    ntouch_one++;

                    // advect one particle from cell to cell and thru surf collides til done

                    //int iterate = 0;

                    while (true)
                    {
                        outface = (int)Enum1.INTERIOR;
                        frac = 1.0;

                        if (xnew[0] < lo[0])
                        {
                            frac = (lo[0] - x[0]) / (xnew[0] - x[0]);
                            outface = (int)Enum1.XLO;
                        }
                        else if (xnew[0] >= hi[0])
                        {
                            frac = (hi[0] - x[0]) / (xnew[0] - x[0]);
                            outface = (int)Enum1.XHI;
                        }

                        if (DIM != 1)
                        {
                            if (xnew[1] < lo[1])
                            {
                                newfrac = (lo[1] - x[1]) / (xnew[1] - x[1]);
                                if (newfrac < frac)
                                {
                                    frac = newfrac;
                                    outface = (int)Enum1.YLO;
                                }
                            }
                            else if (xnew[1] >= hi[1])
                            {
                                newfrac = (hi[1] - x[1]) / (xnew[1] - x[1]);
                                if (newfrac < frac)
                                {
                                    frac = newfrac;
                                    outface = (int)Enum1.YHI;
                                }
                            }
                        }

                        if (DIM == 1)
                        {
                            if (pflag == Enum4.PEXIT && x[1] == lo[1])
                            {
                                frac = 0.0;
                                outface = (int)Enum1.YLO;
                            }
                            else if (Geometry.axi_horizontal_line(dtremain, x, v, lo[1], out itmp, out tc, out tmp))
                            {
                                newfrac = tc / dtremain;
                                if (newfrac < frac)
                                {
                                    frac = newfrac;
                                    outface = (int)Enum1.YLO;
                                }
                            }
                            if (pflag == Enum4.PEXIT && x[1] == hi[1])
                            {
                                frac = 0.0;
                                outface = (int)Enum1.YHI;
                            }
                            else
                            {
                                rnew = Math.Sqrt(xnew[1] * xnew[1] + xnew[2] * xnew[2]);
                                if (rnew >= hi[1])
                                {
                                    if (Geometry.axi_horizontal_line(dtremain, x, v, hi[1], out itmp, out tc, out tmp))
                                    {
                                        newfrac = tc / dtremain;
                                        if (newfrac < frac)
                                        {
                                            frac = newfrac;
                                            outface = (int)Enum1.YHI;
                                        }
                                    }
                                }
                            }

                            pflag = 0;
                        }
                        if (DIM == 3)
                        {
                            if (xnew[2] < lo[2])
                            {
                                newfrac = (lo[2] - x[2]) / (xnew[2] - x[2]);
                                if (newfrac < frac)
                                {
                                    frac = newfrac;
                                    outface = (int)Enum1.ZLO;
                                }
                            }
                            else if (xnew[2] >= hi[2])
                            {
                                newfrac = (hi[2] - x[2]) / (xnew[2] - x[2]);
                                if (newfrac < frac)
                                {
                                    frac = newfrac;
                                    outface = (int)Enum1.ZHI;
                                }
                            }
                        }


                        if (SURF != 0)
                        {
                            nsurf = cells[icell].nsurf;
                            if (nsurf != 0)
                            {


                                // particle crosses cell face, reset xnew exactly on face of cell
                                // so surface check occurs only for particle path within grid cell
                                // xhold = saved xnew so can restore below if no surf collision

                                if (outface != (int)Enum1.INTERIOR)
                                {
                                    xhold[0] = xnew[0];
                                    xhold[1] = xnew[1];
                                    if (DIM != 2) xhold[2] = xnew[2];

                                    xnew[0] = x[0] + frac * (xnew[0] - x[0]);
                                    xnew[1] = x[1] + frac * (xnew[1] - x[1]);
                                    if (DIM != 2) xnew[2] = x[2] + frac * (xnew[2] - x[2]);

                                    if (outface == (int)Enum1.XLO) xnew[0] = lo[0];
                                    else if (outface == (int)Enum1.XHI) xnew[0] = hi[0];
                                    else if (outface == (int)Enum1.YLO) xnew[1] = lo[1];
                                    else if (outface == (int)Enum1.YHI) xnew[1] = hi[1];
                                    else if (outface == (int)Enum1.ZLO) xnew[2] = lo[2];
                                    else if (outface == (int)Enum1.ZHI) xnew[2] = hi[2];
                                }
                                // for axisymmetric, dtsurf = time that particle stays in cell
                                // used as arg to axi_line_intersect()

                                if (DIM == 1)
                                {
                                    if (outface == (int)Enum1.INTERIOR) dtsurf = dtremain;
                                    else dtsurf = dtremain * frac;
                                }
                                // check for collisions with triangles or lines in cell
                                // find 1st surface hit via minparam
                                // skip collisions with previous surf, but not for axisymmetric
                                // not considered collision if 2 params are tied and one INSIDE surf
                                // if collision occurs, perform collision with surface model
                                // reset x,v,xnew,dtremain and continue single particle trajectory

                                cflag = 0;
                                minparam = 2.0;
                                csurfs = cells[icell].csurfs;
                                tri = default(Surf.Tri);
                                line = default(Surf.Line);
                                for (m = 0; m < nsurf; m++)
                                {
                                    isurf = csurfs[m];

                                    if (DIM > 1)
                                    {
                                        if (isurf == exclude) continue;
                                    }
                                    if (DIM == 3)
                                    {
                                        tri = tris[isurf];
                                        hitflag = Geometry.line_tri_intersect(x, xnew,
                                                             pts[tri.p1].x, pts[tri.p2].x,
                                                             pts[tri.p3].x, tri.norm, xc, out param, out side);
                                    }
                                    if (DIM == 2)
                                    {
                                        line = lines[isurf];
                                        hitflag = Geometry.line_line_intersect(x, xnew,
                                                              pts[line.p1].x, pts[line.p2].x,
                                                              line.norm, xc, out param, out side);
                                    }
                                    if (DIM == 1)
                                    {
                                        line = lines[isurf];
                                        hitflag = Geometry.axi_line_intersect(dtsurf, x, v, outface, lo, hi,
                                                             pts[line.p1].x, pts[line.p2].x,
                                                             line.norm, (exclude == isurf) ? 0 : 1,
                                                             xc, vc, out param, out side);
                                    }
                                    if (hitflag && param < minparam && side == (int)Enum3.OUTSIDE)
                                    {
                                        cflag = 1;
                                        minparam = param;
                                        minside = side;
                                        minsurf = isurf;
                                        minxc[0] = xc[0];
                                        minxc[1] = xc[1];
                                        if (DIM == 3) minxc[2] = xc[2];
                                        if (DIM == 1)
                                        {
                                            minvc[1] = vc[1];
                                            minvc[2] = vc[2];
                                        }
                                    }
                                }

                                nscheck_one += nsurf;

                                if (cflag != 0)
                                {
                                    // NOTE: this check is no longer needed?
                                    if (minside == (int)Enum3.INSIDE)
                                    {
                                        string str = string.Format("Particle {0} on proc {1} hit inside of surf {2} on step {3} ",
                                                 i, me, minsurf, sparta.update.ntimestep);
                                        sparta.error.one(str);
                                    }

                                    if (DIM == 3) tri = tris[minsurf];
                                    if (DIM != 3) line = lines[minsurf];

                                    // set x to collision point
                                    // if axisymmetric, set v to remapped velocity at collision pt

                                    x[0] = minxc[0];
                                    x[1] = minxc[1];
                                    if (DIM == 3) x[2] = minxc[2];
                                    if (DIM == 1)
                                    {
                                        v[1] = minvc[1];
                                        v[2] = minvc[2];
                                    }
                                    // perform surface collision using surface collision model
                                    // surface chemistry may destroy particle or create new one
                                    // must update particle's icell to current icell so that
                                    //   if jpart is created, it will be added to correct cell
                                    // if jpart, add new particle to this iteration via pstop++
                                    // tally surface statistics if requested using iorig

                                    Particle.OnePart p = particles[i];
                                    p.icell = icell;
                                    ipart = p;

                                    dtremain *= 1.0 - minparam * frac;

                                    if (nsurf_tally != 0)
                                    {
                                        iorig = particles[i];
                                        //memcpy(&iorig, &particles[i], sizeof(Particle::OnePart));
                                    }

                                    if (DIM == 3)
                                        jpart = sparta.surf.sc[tri.isc].
                                          collide(ref ipart, tri.norm, dtremain, tri.isr);
                                    if (DIM != 3)
                                        jpart = sparta.surf.sc[line.isc].
                                          collide(ref ipart, line.norm, dtremain, line.isr);

                                    if (jpart != null)
                                    {
                                        Particle.OnePart pp = (Particle.OnePart)jpart;
                                        particles = sparta.particle.particles;
                                        x = particles[i].x;
                                        v = particles[i].v;
                                        pp.flag = (int)Enum4.PSURF + 1 + minsurf;
                                        pp.dtremain = dtremain;
                                        pp.weight = particles[i].weight;
                                        pstop++;
                                    }

                                    if (nsurf_tally != 0)
                                        for (m = 0; m < nsurf_tally; m++)
                                            slist_active[m].surf_tally(minsurf, iorig, ipart, jpart);

                                    // nstuck = consective iterations particle is immobile

                                    if (minparam == 0.0) stuck_iterate++;
                                    else stuck_iterate = 0;

                                    // reset post-bounce xnew

                                    xnew[0] = x[0] + dtremain * v[0];
                                    xnew[1] = x[1] + dtremain * v[1];
                                    if (DIM != 2) xnew[2] = x[2] + dtremain * v[2];

                                    exclude = minsurf;
                                    nscollide_one++;
                                    if (ipart == null) particles[i].flag = (int)Enum4.PDISCARD;
                                    else if (stuck_iterate < MAXSTUCK) continue;
                                    else
                                    {
                                        particles[i].flag = (int)Enum4.PDISCARD;
                                        nstuck++;
                                    }
                                }
                                if (outface != (int)Enum1.INTERIOR)
                                {
                                    xnew[0] = xhold[0];
                                    xnew[1] = xhold[1];
                                    if (DIM != 2) xnew[2] = xhold[2];
                                }
                            }
                        }
                        // break from advection loop if discarding particle

                        if (particles[i].flag == (int)Enum4.PDISCARD) break;

                        // no cell crossing and no surface collision
                        // set final particle position to xnew, then break from advection loop
                        // for axisymmetry, must first remap linear xnew and v
                        // if migrating to another proc,
                        //   flag as PDONE so new proc won't move it more on this step

                        if (outface == (int)Enum1.INTERIOR)
                        {
                            if (DIM == 1) Console.WriteLine("axi_remap(xnew, v);");
                            x[0] = xnew[0];
                            x[1] = xnew[1];
                            if (DIM == 3) x[2] = xnew[2];
                            if (cells[icell].proc != me) particles[i].flag = (int)Enum4.PDONE;
                            break;
                        }

                        // particle crosses cell face
                        // decrement dtremain in case particle is passed to another proc
                        // for axisymmetry, must then remap linear x and v
                        // reset particle x to be exactly on cell face
                        // for axisymmetry, must reset xnew for next iteration since v changed

                        dtremain *= 1.0 - frac;
                        exclude = -1;

                        x[0] += frac * (xnew[0] - x[0]);
                        x[1] += frac * (xnew[1] - x[1]);
                        if (DIM != 2) x[2] += frac * (xnew[2] - x[2]);
                        if (DIM == 1) Console.WriteLine("axi_remap(x, v);");

                        if (outface == (int)Enum1.XLO) x[0] = lo[0];
                        else if (outface == (int)Enum1.XHI) x[0] = hi[0];
                        else if (outface == (int)Enum1.YLO) x[1] = lo[1];
                        else if (outface == (int)Enum1.YHI) x[1] = hi[1];
                        else if (outface == (int)Enum1.ZLO) x[2] = lo[2];
                        else if (outface == (int)Enum1.ZHI) x[2] = hi[2];

                        if (DIM == 1)
                        {
                            xnew[0] = x[0] + dtremain * v[0];
                            xnew[1] = x[1] + dtremain * v[1];
                            xnew[2] = x[2] + dtremain * v[2];
                        }
                        // nflag = type of neighbor cell: child, parent, unknown, boundary
                        // if parent, use id_find_child to identify child cell
                        //   result of id_find_child could be unknown:
                        //     particle is hitting face of a ghost child cell which extends
                        //     beyond my ghost halo, cell on other side of face is a parent,
                        //     it's child which the particle is in is entirely beyond my halo
                        // if new cell is child and surfs exist, check if a split cell

                        nflag = (Enum5)sparta.grid.neigh_decode(nmask, outface);
                        icell_original = icell;

                        switch (nflag)
                        {
                            case Enum5.NCHILD:
                                icell = neigh[outface];
                                if (DIM == 3 && SURF!=0)
                                {
                                    if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                                        icell = split3d(icell, x);
                                }
                                if (DIM < 3 && SURF != 0)
                                {
                                    if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                                        icell = split2d(icell, x);
                                }
                                break;
                            case Enum5.NPARENT:
                                icell = sparta.grid.id_find_child(neigh[outface], x);
                                if (icell >= 0)
                                {
                                    if (DIM == 3 && SURF != 0)
                                    {
                                        if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                                            icell = split3d(icell, x);
                                    }
                                    if (DIM < 3 && SURF != 0)
                                    {
                                        if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                                            icell = split2d(icell, x);
                                    }
                                }
                                break;
                            case Enum5.NUNKNOWN:
                                icell = -1;
                                // neighbor cell is global boundary
                                // tally boundary stats if requested using iorig
                                // collide() updates x,v,xnew as needed due to boundary interaction
                                //   may also update dtremain (piston BC)
                                // for axisymmetric, must recalculate xnew since v may have changed
                                // surface chemistry may destroy particle or create new one
                                // if jpart, add new particle to this iteration via pstop++
                                // OUTFLOW: exit with particle flag = PDISCARD
                                // PERIODIC: new cell via same logic as above for child/parent/unknown
                                // other = reflected particle stays in same grid cell
                                break;
                            
                            default:
                                ipart = particles[i];

                                if (nboundary_tally!=0)
                                {
                                    iorig = particles[i];
                                    //memcpy(&iorig, &particles[i], sizeof(Particle::OnePart));
                                }

                                bflag = sparta.domain.collide(ipart, outface, icell, xnew, dtremain, jpart);

                                if (jpart!=null)
                                {
                                    particles = sparta.particle.particles;
                                    x = particles[i].x;
                                    v = particles[i].v;
                                }

                                if (nboundary_tally!=0)
                                    for (m = 0; m < nboundary_tally; m++)
                                        blist_active[m].boundary_tally(outface, bflag, iorig, ipart, jpart);

                                if (DIM == 1)
                                {
                                    xnew[0] = x[0] + dtremain * v[0];
                                    xnew[1] = x[1] + dtremain * v[1];
                                    xnew[2] = x[2] + dtremain * v[2];
                                }

                                if (bflag == (int)Enum2.OUTFLOW)
                                {
                                    particles[i].flag = (int)Enum4.PDISCARD;
                                    nexit_one++;
                                    break;

                                }
                                else if (bflag == (int)Enum2.PERIODIC)
                                {
                                    if (nflag == Enum5.NPBCHILD)
                                    {
                                        icell = neigh[outface];
                                        if (DIM == 3 && SURF!=0)
                                        {
                                            if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                                                icell = split3d(icell, x);
                                        }
                                        if (DIM < 3 && SURF != 0)
                                        {
                                            if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                                                icell = split2d(icell, x);
                                        }
                                    }
                                    else if (nflag == Enum5.NPBPARENT)
                                    {
                                        icell = sparta.grid.id_find_child(neigh[outface], x);
                                        if (icell >= 0)
                                        {
                                            if (DIM == 3 && SURF != 0)
                                            {
                                                if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                                                    icell = split3d(icell, x);
                                            }
                                            if (DIM < 3 && SURF != 0)
                                            {
                                                if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                                                    icell = split2d(icell, x);
                                            }
                                        }
                                        else sparta.domain.uncollide(outface, x);
                                    }
                                    else if (nflag == Enum5.NPBUNKNOWN)
                                    {
                                        icell = -1;
                                        sparta.domain.uncollide(outface, x);
                                    }

                                }
                                else if (bflag == (int)Enum2.SURFACE)
                                {
                                    if (ipart == null)
                                    {
                                        particles[i].flag = (int)Enum4.PDISCARD;
                                        break;
                                    }
                                    else if (jpart!=null)
                                    {
                                        Particle.OnePart p = (Particle.OnePart)jpart;
                                        p.flag = (int)Enum4.PSURF;
                                        p.dtremain = dtremain;
                                        p.weight = particles[i].weight;
                                        jpart = p;
                                        pstop++;
                                    }
                                    nboundary_one++;
                                    ntouch_one--;    // decrement here since will increment below

                                }
                                else
                                {
                                    nboundary_one++;
                                    ntouch_one--;    // decrement here since will increment below
                                }
                                break;
                        }

                        // neighbor cell is unknown
                        // reset icell to original icell which must be a ghost cell
                        // exit with particle flag = PEXIT, so receiver can identify neighbor

                        if (icell < 0)
                        {
                            icell = icell_original;
                            particles[i].flag = (int)Enum4.PEXIT;
                            particles[i].dtremain = dtremain;
                            entryexit = 1;
                            break;
                        }

                        // if nsurf < 0, new cell is EMPTY ghost
                        // exit with particle flag = PENTRY, so receiver can continue move

                        if (cells[icell].nsurf < 0)
                        {
                            particles[i].flag = (int)Enum4.PENTRY;
                            particles[i].dtremain = dtremain;
                            entryexit = 1;
                            break;
                        }

                        // move particle into new grid cell for next stage of move

                        lo = cells[icell].lo;
                        hi = cells[icell].hi;
                        neigh = cells[icell].neigh;
                        nmask = cells[icell].nmask;
                        ntouch_one++;

                    }
                    // END of while loop over advection of single particle

                    // move is complete, or as much as can be done on this proc
                    // update particle's grid cell
                    // if particle flag set, add particle to migrate list
                    // if discarding, migration will delete particle

                    particles[i].icell = icell;

                    if (particles[i].flag != (int)Enum4.PKEEP)
                    {
                        mlist[nmigrate++] = i;
                        if (particles[i].flag != (int)Enum4.PDISCARD)
                        {
                            if (cells[icell].proc == me)
                            {
                                string str=string.Format("Particle {0} on proc {1} being sent to self on step {2}",
                                        i, me, sparta.update.ntimestep);
                                sparta.error.one( str);
                            }
                            ncomm_one++;
                        }
                    }

                }

                // END of pstart/pstop loop advecting all particles

                // if gridcut >= 0.0, check if another iteration of move is required
                // only the case if some particle flag = PENTRY/PEXIT
                //   in which case perform particle migration
                // if not, move is done and final particle comm will occur in run()
                // if iterating, reset pstart/pstop and extend migration list if necessary

                if (sparta.grid.cutoff < 0.0) break;

                sparta.mpi.MPI_Allreduce(ref entryexit, ref any_entryexit, 1, MPI.MPI_INT, MPI.MPI_MAX, sparta.world);
                if (any_entryexit!=0)
                {
                    sparta.timer.stamp((int)Timer.Enum1.TIME_MOVE);
                    pstart = sparta.comm.migrate_particles(nmigrate, mlist);
                    sparta.timer.stamp((int)Timer.Enum1.TIME_COMM);
                    pstop = sparta.particle.nlocal;
                    if (pstop - pstart > maxmigrate)
                    {
                        maxmigrate = pstop - pstart;
                        //memory.destroy(mlist);
                        //memory.create(mlist, maxmigrate, "particle:mlist");
                        mlist = new int[maxmigrate];
                    }
                }
                else break;

                // END of single move/migrate iteration



            }




        }
        //      typedef void (Update::* FnPtr2) (double, double*, double*);
        //FnPtr2 moveperturb;        // ptr to moveperturb method

        //      // variants of moveperturb method
        //      // adjust end-of-move x,v due to perturbation on straight-line advection

        int perturbflag;
        void gravity2d(double dt, double[] x, double[] v)
        {
            double dtsq = 0.5 * dt * dt;
            x[0] += dtsq * gravity[0];
            x[1] += dtsq * gravity[1];
            v[0] += dt * gravity[0];
            v[1] += dt * gravity[1];
        }

        void gravity3d(double dt, double[] x, double[] v)
        {
            double dtsq = 0.5 * dt * dt;
            x[0] += dtsq * gravity[0];
            x[1] += dtsq * gravity[1];
            x[2] += dtsq * gravity[2];
            v[0] += dt * gravity[0];
            v[1] += dt * gravity[1];
            v[2] += dt * gravity[2];
        }
    }
}