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

        //public void init();
        //public void setup();
        //public virtual void start_of_step();
        //public virtual void end_of_step();

        //public virtual void add_grid_one(int, int);
        public virtual int pack_grid_one(int icell,ref StringBuilder buf, int memflag)
        {
            int n=0;
            for (int i = 0; i < n_pergrid; i++)
                n += fix[list_pergrid[i]].pack_grid_one(icell,ref buf, memflag);
            return n;
        }
        //public virtual int unpack_grid_one(int, char*);
        //public virtual void compress_grid(int);

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

        //public void clearstep_compute();
        //public void addstep_compute(bigint);
        //public void addstep_compute_all(bigint);

        //public void list_init_fixes();
        //public void list_init_computes();

        //public virtual void add_particle(int, double, double, double, double*);
        //public virtual void gas_react(int);
        //public virtual void surf_react(Particle::OnePart*, int &, int &);

        //public bigint memory_usage();
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
         
        //protected void list_init(int, int &, int[]&);
        //protected void list_init_end_of_step(int, int &, int[]&);

    }
}