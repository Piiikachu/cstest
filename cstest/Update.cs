using System;
using System.Collections.Generic;
using bigint = System.Int64;

namespace cstest
{
    public class Update
    {
        enum Enum1 { XLO, XHI, YLO, YHI, ZLO, ZHI, INTERIOR };         // same as Domain
        enum Enum2 { PERIODIC, OUTFLOW, REFLECT, SURFACE, AXISYM };  // same as Domain
        enum Enum3 { OUTSIDE, INSIDE, ONSURF2OUT, ONSURF2IN };      // several files
        enum Enum4 { PKEEP, PINSERT, PDONE, PDISCARD, PENTRY, PEXIT, PSURF };   // several files
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
            else sparta.error.all( "Illegal units command");

            //delete[] unit_style;
            int n = style.Length + 1;
            unit_style=string.Copy(style);
        }
        //      public virtual void init();
        //      public virtual void setup();
        //      public virtual void run(int);
        public void global(int nnarg,string[] arg)
        {
            int narg = nnarg + 1;
            if (narg < 1) sparta.error.all("Illegal global command");
            int iarg = 1;
            while (iarg<narg)
            {
                if (string.Equals(arg[iarg],"fnum"))
                {
                    if (iarg + 2 > narg) sparta.error.all( "Illegal global command");
                    fnum = System.Double.Parse(arg[iarg + 1]);
                    if (fnum <= 0.0) sparta.error.all( "Illegal global command");
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
                    if (iarg + 4 > narg) sparta.error.all( "Illegal global command");
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
                    if (iarg + 5 > narg) sparta.error.all( "Illegal global command");
                    double gmag = System.Double.Parse(arg[iarg + 1]);
                    gravity[0] = System.Double.Parse(arg[iarg + 2]);
                    gravity[1] = System.Double.Parse(arg[iarg + 3]);
                    gravity[2] = System.Double.Parse(arg[iarg + 4]);
                    if (gmag < 0.0) sparta.error.all( "Illegal global command");
                    if (gmag > 0.0 &&
                        gravity[0] == 0.0 && gravity[1] == 0.0 && gravity[2] == 0.0)
                        sparta.error.all( "Illegal global command");
                    if (gmag > 0.0) MathExtra.snorm3(gmag, gravity);
                    else gravity[0] = gravity[1] = gravity[2] = 0.0;
                    iarg += 5;

                }
                else if (string.Equals(arg[iarg], "surfmax"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    if (sparta.surf.exist!=0)
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
                    if (iarg + 2 > narg) sparta.error.all( "Illegal global command");
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
                        if (iarg + 4 > narg) sparta.error.all( "Illegal global command");
                        sparta.surf.pushflag = 2;
                        sparta.surf.pushlo = sparta.input.numeric( arg[iarg + 1]);
                        sparta.surf.pushhi = sparta.input.numeric( arg[iarg + 2]);
                        sparta.surf.pushvalue = sparta.input.numeric( arg[iarg + 3]);
                        if (sparta.surf.pushlo > sparta.surf.pushhi)
                            sparta.error.all( "Illegal global command");
                        iarg += 4;
                    }
                }
                else if (string.Equals(arg[iarg], "gridcut"))
                {
                    if (iarg + 2 > narg) sparta.error.all( "Illegal global command");
                    sparta.grid.cutoff = double.Parse(arg[iarg + 1]);
                    if (sparta.grid.cutoff < 0.0 && sparta.grid.cutoff != -1.0)
                        sparta.error.all( "Illegal global command");
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
                    if (string.Equals(arg[iarg + 1], "auto")) sparta.surf.tally_comm= (int)Enum6.TALLYAUTO;
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
                        Console.WriteLine("update.global()->try copy array");
                    }
                    else sparta.error.all("Illegal global command");
                    iarg += 3;
                }
                else if (string.Equals(arg[iarg], "particle/reordor"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal global command");
                    reorder_period = int.Parse(arg[iarg + 1]);
                    if (reorder_period<0)
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

        //      public int split3d(int, double*);
        //      public int split2d(int, double*);

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

        //      int collide_react_setup();
        //      void collide_react_update();

        //      int bounce_setup();
        //      void bounce_set(bigint);
        //      void reset_timestep(bigint);

        //      //int axi_vertical_line(double, double *, double *, double, double, double,
        //      //                     double &);

        //      // remap x and v components into axisymmetric plane
        //      // input x at end of linear move (x = xold + dt*v)
        //      // change x[1] = sqrt(x[1]^2 + x[2]^2), x[2] = 0.0
        //      // change vy,vz by rotation into axisymmetric plane

        //      inline void axi_remap(double* x, double* v)
        //      {
        //          double ynew = x[1];
        //          double znew = x[2];
        //          x[1] = sqrt(ynew * ynew + znew * znew);
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
        //      template<int, int> void move();

        //      int perturbflag;
        //      typedef void (Update::* FnPtr2) (double, double*, double*);
        //FnPtr2 moveperturb;        // ptr to moveperturb method

        //      // variants of moveperturb method
        //      // adjust end-of-move x,v due to perturbation on straight-line advection

        //      inline void gravity2d(double dt, double* x, double* v)
        //      {
        //          double dtsq = 0.5 * dt * dt;
        //          x[0] += dtsq * gravity[0];
        //          x[1] += dtsq * gravity[1];
        //          v[0] += dt * gravity[0];
        //          v[1] += dt * gravity[1];
        //      };

        //      inline void gravity3d(double dt, double* x, double* v)
        //      {
        //          double dtsq = 0.5 * dt * dt;
        //          x[0] += dtsq * gravity[0];
        //          x[1] += dtsq * gravity[1];
        //          x[2] += dtsq * gravity[2];
        //          v[0] += dt * gravity[0];
        //          v[1] += dt * gravity[1];
        //          v[2] += dt * gravity[2];
        //      };
    }
}