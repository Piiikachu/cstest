using System.Collections.Generic;
using System.Text;
using bigint = System.Int64;
namespace cstest
{
    public class Modify
    {
        public const int DELTA = 4;

        // mask settings - same as in fix.cpp

        public const int START_OF_STEP = 1;
        public const int END_OF_STEP = 2;

        public int nfix, maxfix;
        public int n_start_of_step, n_end_of_step;
        public int n_pergrid, n_add_particle, n_gas_react, n_surf_react;
         
        public List<Fix> fix;           // list of fixes
        public int[] fmask;                // bit mask for when each fix is applied
         
        public int ncompute, maxcompute;   // list of computes
        public List<Compute> compute;

        public void init()
        {
            int i;

            // create lists of fixes with masks for calling at each stage of run

            list_init(START_OF_STEP,ref n_start_of_step,ref list_start_of_step);
            list_init_end_of_step(END_OF_STEP,ref n_end_of_step, ref list_end_of_step);

            // create other lists of fixes and computes

            list_init_fixes();
            list_init_computes();

            // init each fix

            for (i = 0; i < nfix; i++) fix[i].init();

            // init each compute
            // set invoked_scalar,vector,etc to -1 to force new run to re-compute them
            // add initial timestep to all computes that store invocation times
            //   since any of them may be invoked by initial thermo
            // do not clear out invocation times stored within a compute,
            //   b/c some may be holdovers from previous run, like for ave fixes

            for (i = 0; i < ncompute; i++)
            {
                compute[i].init();
                compute[i].invoked_scalar = -1;
                compute[i].invoked_vector = -1;
                compute[i].invoked_array = -1;
                compute[i].invoked_per_particle = -1;
                compute[i].invoked_per_grid = -1;
                compute[i].invoked_per_surf = -1;
            }
            addstep_compute_all(sparta.update.ntimestep);

        }
        public void setup()
        {
            // setup each fix

            for (int i = 0; i < nfix; i++) fix[i].setup();
        }
        public virtual void start_of_step()
        {
            for (int i = 0; i < n_start_of_step; i++)
                fix[list_start_of_step[i]].start_of_step();
        }
        public virtual void end_of_step()
        {
            for (int i = 0; i < n_end_of_step; i++)
                if (sparta.update.ntimestep % end_of_step_every[i] == 0)
                    fix[list_end_of_step[i]].end_of_step();
        }

        //public virtual void add_grid_one(int, int);
        public virtual int pack_grid_one(int icell,ref StringBuilder buf, int memflag)
        {
            int n=0;
            for (int i = 0; i < n_pergrid; i++)
                n += fix[list_pergrid[i]].pack_grid_one(icell,ref buf, memflag);
            return n;
        }
        //public virtual int unpack_grid_one(int, char*);
        public virtual void compress_grid(int flag)
        {
            if (flag == 0)
                for (int i = 0; i < n_pergrid; i++)
                    fix[list_pergrid[i]].compress_grid();
            else
                for (int i = 0; i < n_pergrid; i++)
                    fix[list_pergrid[i]].post_compress_grid();
        }

        public void add_fix(int narg, string[] arg)
        {

            if (sparta.domain.box_exist == 0)
                sparta.error.all("Fix command before simulation box is defined");
            if (narg < 2) sparta.error.all("Illegal fix command");

            // if fix ID exists:
            //   set newflag = 0 so create new fix in same location in fix list
            //   error if new style does not match old style
            //     since can't replace it (all when-to-invoke ptrs would be invalid)
            //   delete old fix, but do not call update_callback(),
            //     since will replace this fix and thus other fix locs will not change
            //   set ptr to NULL in case new fix scans list of fixes,
            //     e.g. scan will occur in add_callback() if called by new fix
            // if fix ID does not exist:
            //   set newflag = 1 so create new fix
            //   extend fix and fmask lists as necessary

            int ifix, newflag;
            for (ifix = 0; ifix < nfix; ifix++)
                if (string.Equals(arg[0], fix[ifix].id)) break;

            if (ifix < nfix)
            {
                newflag = 0;
                if (!string.Equals(arg[1], fix[ifix].style))
                    sparta.error.all("Replacing a fix, but new style != old style");

                fix[ifix] = null;
            }
            else
            {
                newflag = 1;
                if (nfix == maxfix)
                {
                    maxfix += DELTA;
                    fix = new List<Fix>(maxfix);
                    fmask = new int[maxfix];
                }
            }

            //todo :create the fix

            int found = 0;
            if (sparta.suffix_enable!=0)
            {
                if (sparta.suffix!=null)
                {
                    string estyle = string.Format("{0}/{1}", arg[1], sparta.suffix);
                }
            }

            switch (arg[1])
            {
                case "emit/face":
                    fix.Add(new FixEmitFace(sparta, narg, arg));
                    break;
                default:
                    sparta.error.all("Unrecognized fix style");
                    break;
            }

            fmask[ifix] = fix[ifix].setmask();
            if (newflag!=0) nfix++;

        }
        //public void delete_fix(const char*);
        //public int find_fix(const char*);

        //public void add_compute(int, char**);
        //public void delete_compute(const char*);
        //public int find_compute(const char*);

        public void clearstep_compute()
        {
            for (int icompute = 0; icompute < ncompute; icompute++)
                compute[icompute].invoked_flag = 0;
        }
        public void addstep_compute(bigint newstep)
        {
            for (int icompute = 0; icompute < n_timeflag; icompute++)
                if (compute[list_timeflag[icompute]].invoked_flag!=0)
                    compute[list_timeflag[icompute]].addstep(newstep);
        }
        public void addstep_compute_all(bigint newstep)
        {
            for (int icompute = 0; icompute < n_timeflag; icompute++)
                if (compute[list_timeflag[icompute]].invoked_flag!=0)
                    compute[list_timeflag[icompute]].addstep(newstep);
        }

        public void list_init_fixes()
        {
            n_pergrid = n_add_particle = n_gas_react = n_surf_react = 0;
            for (int i = 0; i < nfix; i++)
            {
                if (fix[i].gridmigrate!=0) n_pergrid++;
                if (fix[i].flag_add_particle != 0) n_add_particle++;
                if (fix[i].flag_gas_react != 0) n_gas_react++;
                if (fix[i].flag_surf_react != 0) n_surf_react++;
            }

            list_pergrid = new int[n_pergrid];
            list_add_particle = new int[n_add_particle];
            list_gas_react = new int[n_gas_react];
            list_surf_react = new int[n_surf_react];

            n_pergrid = n_add_particle = n_gas_react = n_surf_react = 0;
            for (int i = 0; i < nfix; i++)
            {
                if (fix[i].gridmigrate != 0) list_pergrid[n_pergrid++] = i;
                if (fix[i].flag_add_particle != 0) list_add_particle[n_add_particle++] = i;
                if (fix[i].flag_gas_react != 0) list_gas_react[n_gas_react++] = i;
                if (fix[i].flag_surf_react != 0) list_surf_react[n_surf_react++] = i;
            }
        }
        public void list_init_computes()
        {
            n_timeflag = 0;
            for (int i = 0; i < ncompute; i++)
                if (compute[i].timeflag!=0) n_timeflag++;
            list_timeflag = new int[n_timeflag];

            n_timeflag = 0;
            for (int i = 0; i < ncompute; i++)
                if (compute[i].timeflag != 0) list_timeflag[n_timeflag++] = i;
        }

        public virtual void add_particle(int index, double temp_thermal,
                          double temp_rot, double temp_vib, double[] vstream)
        {
            for (int i = 0; i < n_add_particle; i++)
                fix[list_add_particle[i]].add_particle(index, temp_thermal, temp_rot,
                                                        temp_vib, vstream);
        }
        //public virtual void gas_react(int);
        //public virtual void surf_react(Particle::OnePart*, int &, int &);

        public bigint memory_usage()
        {
            bigint bytes = 0;
            for (int i = 0; i < nfix; i++) bytes += (long)fix[i].memory_usage();
            for (int i = 0; i < ncompute; i++) bytes += compute[i].memory_usage();
            return bytes;
        }
        private SPARTA sparta;

        public Modify(SPARTA sparta)
        {
            this.sparta = sparta;
            nfix = maxfix = 0;
            n_start_of_step = n_end_of_step = 0;

            fix = null;
            fmask = null;
            list_start_of_step = list_end_of_step = null;

            end_of_step_every = null;
            list_pergrid = null;
            list_add_particle = null;
            list_gas_react = null;
            list_surf_react = null;
            list_timeflag = null;

            ncompute = maxcompute = 0;
            compute = null;
        }

        protected int[] list_start_of_step,list_end_of_step;
         
        protected int[] end_of_step_every;
         
        protected int[] list_pergrid;         // list of fixes that store per grid cell info
        protected int[] list_add_particle;    // list of fixes with add_particle() method
        protected int[] list_gas_react;       // list of fixes with gas_react() method
        protected int[] list_surf_react;      // list of fixes with surf_react() method
         
        protected int n_timeflag;            // list of computes that store time invocation
        protected int[] list_timeflag;

        protected void list_init(int mask,ref int n,ref int[] list)
        {
            n = 0;
            for (int i = 0; i < nfix; i++) if (fmask[i]!=0 & mask != 0) n++;
            list = new int[n];

            n = 0;
            for (int i = 0; i < nfix; i++) if (fmask[i] != 0 & mask != 0) list[n++] = i;
        }
        protected void list_init_end_of_step(int mask,ref int n,ref int[] list)
        {
            n = 0;
            for (int i = 0; i < nfix; i++) if (fmask[i] != 0 & mask != 0) n++;
            list = new int[n];
            end_of_step_every = new int[n];

            n = 0;
            for (int i = 0; i < nfix; i++)
                if (fmask[i] != 0 & mask != 0)
                {
                    list[n] = i;
                    end_of_step_every[n++] = fix[i].nevery;
                }
        }

    }
}