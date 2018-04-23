using bigint = System.Int64;
namespace cstest
{
    public class Compute
    {
        public const int DELTA = 4;
        public string id,style;
         
        public double scalar;            // computed global scalar
        public double[] vector;           // computed global vector
        public double[,] array;           // computed global array
        public double[] vector_particle;  // computed per-particle vector
        public double[,] array_particle;  // computed per-particle array
        public double[] vector_grid;      // computed per-grid vector
        public double[,] array_grid;      // computed per-grid array
                 
                 // vec/array_surf are length nslocal = # of owned surf elements
                 // tally vec/array are length nlocal = # of unique surf elements tallied
                 // tally info is accessed by callers via surfinfo()
                 
        public double[] vector_surf;        // computed per-surf vector
        public double[,] array_surf;        // computed per-surf array
        public double[] vector_surf_tally;  // computed per-surf tally vector
        public double[,] array_surf_tally;  // computed per-surf tally array
         
         // NOTE: get rid of these fields?
        public double[,] array_grid_extra;   // extra per-grid array
        public double[,] norm_grid_extra;    // extra per-grid normalizations
         
        public int scalar_flag;          // 0/1 if compute_scalar() function exists
        public int vector_flag;          // 0/1 if compute_vector() function exists
        public int array_flag;           // 0/1 if compute_array() function exists
        public int size_vector;          // length of global vector
        public int size_array_rows;      // rows in global array
        public int size_array_cols;      // columns in global array
         
        public int per_particle_flag;      // 0/1 if compute_per_particle() function exists
        public int size_per_particle_cols; // 0 = vector, N = columns in per-particle array
         
        public int per_grid_flag;          // 0/1 if compute_per_grid() function exists
        public int size_per_grid_cols;     // 0 = vector, N = columns in per-grid array
         
        public int post_process_grid_flag;   // 1 if requires post_processing for output
         
         // NOTE: get rid of this field?
        public int size_per_grid_extra_cols; // 0 = none, N = columns in extra per-grid array
         
        public int per_surf_flag;          // 0/1 if compute_per_surf() function exists
        public int size_per_surf_cols;     // 0 = vector, N = columns in per-surf array
         
        public int surf_tally_flag;        // 1 if compute tallies surface bounce info
        public int boundary_tally_flag;    // 1 if compute tallies boundary bounce info
         
        public int timeflag;       // 1 if Compute stores list of timesteps it's called on
        public int ntime;          // # of entries in time list
        public int maxtime;        // max # of entries time list can hold
        public bigint[] tlist;      // list of timesteps the Compute is called on
         
        public int invoked_flag;       // non-zero if invoked or accessed this step, 0 if not
        public bigint invoked_scalar;  // last timestep on which compute_scalar() was invoked
        public bigint invoked_vector;       // ditto for compute_vector()
        public bigint invoked_array;        // ditto for compute_array()
        public bigint invoked_per_particle; // ditto for compute_per_particle()
        public bigint invoked_per_grid;     // ditto for compute_per_grid()
        public bigint invoked_per_surf;     // ditto for compute_per_surf()

        public Compute(SPARTA sparta, int narg, string[] arg)
        {
            if (narg < 2) sparta.error.all("Illegal compute command");

            // compute ID and style
            // ID must be all alphanumeric chars or underscores

            int n = arg[0].Length + 1;
            id = string.Copy(arg[0]);

            for (int i = 0; i < n - 1; i++)
                if (char.IsLetterOrDigit(id[i]) && id[i] != '_')
                    sparta.error.all("Compute ID must be alphanumeric or underscore characters");

            //n = strlen(arg[1]) + 1;
            style = string.Copy( arg[1]);

            // set child class defaults

            scalar_flag = vector_flag = array_flag = 0;
            per_particle_flag = per_grid_flag = per_surf_flag = 0;
            post_process_grid_flag = 0;
            size_per_grid_extra_cols = 0;
            surf_tally_flag = boundary_tally_flag = 0;

            timeflag = 0;
            ntime = maxtime = 0;
            tlist = null;

            invoked_scalar = invoked_vector = invoked_array = -1;
            invoked_per_particle = invoked_per_grid = invoked_per_surf = -1;

            kokkos_flag = 0;
            copymode = 0;
        }
        public virtual void init() { }
         
        public virtual double compute_scalar() { return 0.0; }
        public virtual void compute_vector() { }
        public virtual void compute_array() { }
        public virtual void compute_per_particle() { }
        public virtual void compute_per_grid() { }
        public virtual void compute_per_surf() { }
        public virtual void clear() { }
        //public virtual void surf_tally(int, Particle::OnePart*,
        //                        Particle::OnePart*, Particle::OnePart*)
        //{ }
        //public virtual void boundary_tally(int, int, Particle::OnePart*,
        //                            Particle::OnePart*, Particle::OnePart*)
        //{ }
        //
        //public virtual int query_tally_grid(int, double[,]&, int*&) { return 0; }
        //public virtual double post_process_grid(int, int, int, double[,], int*,
        //                                 double[], int)
        //{ return 0.0; }
        //


        // Kokkos methods

        public int kokkos_flag;          // 1 if Kokkos-enabled
        public int copy, copymode;        // 1 if copy of class (prevents deallocation of
                                          //  base class when child copy is destroyed)

        // NOTE: get rid of these methods
        //public virtual void post_process_grid_old(void*, void*, int, int, double[], int) { }
        //public virtual void normwhich(int, int &, int &) { }
        //public virtual double[] normptr(int) { return NULL; }

        //public virtual int surfinfo(int*&) { return 0; }

        public void addstep(bigint ntimestep)
        {
            // i = location in list to insert ntimestep

            int i;
            for (i = ntime - 1; i >= 0; i--)
            {
                if (ntimestep == tlist[i]) return;
                if (ntimestep < tlist[i]) break;
            }
            i++;

            // extend list as needed

            if (ntime == maxtime)
            {
                maxtime += DELTA;
                tlist = new bigint[maxtime];
                //memory->grow(tlist, maxtime, "compute:tlist");
            }

            // move remainder of list upward and insert ntimestep

            for (int j = ntime - 1; j >= i; j--) tlist[j + 1] = tlist[j];
            tlist[i] = ntimestep;
            ntime++;
        }
        //public int matchstep(bigint);
        //public void clearstep();

        //public virtual void reallocate() { }
        //public virtual bigint memory_usage();

    }
}