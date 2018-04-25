using System.Text;

namespace cstest
{
    public class Fix
    {
        public string id,style;
        public enum ExecutionSpace { Host, Device };

        public int nevery;                    // how often to call an end_of_step fix
        public int time_depend;               // 1 if requires continuous timestepping
        public int gridmigrate;               // 0/1 if per grid cell info must migrate
        public int flag_add_particle;         // 0/1 if has add_particle() method
        public int flag_gas_react;            // 0/1 if has gas_react() method
        public int flag_surf_react;           // 0/1 if has surf_react() method
         
        public int scalar_flag;               // 0/1 if compute_scalar() function exists
        public int vector_flag;               // 0/1 if compute_vector() function exists
        public int array_flag;                // 0/1 if compute_array() function exists
        public int size_vector;               // length of global vector
        public int size_array_rows;           // rows in global array
        public int size_array_cols;           // columns in global array
        public int global_freq;               // frequency s/v data is available at
         
        public int per_particle_flag;         // 0/1 if per-particle data is stored
        public int size_per_particle_cols;    // 0 = vector, N = cols in per-particle array
        public int per_particle_freq;         // frequency per-particle data is available at
         
        public int per_grid_flag;             // 0/1 if per-grid data is stored
        public int size_per_grid_cols;        // 0 = vector, N = cols in per-grid array
        public int per_grid_freq;             // frequency per-grid data is available at
         
        public int per_surf_flag;             // 0/1 if per-surf data is stored
        public int size_per_surf_cols;        // 0 = vector, N = cols in per-surf array
        public int per_surf_freq;             // frequency per-surf data is available at
         
        public double[] vector_particle;       // computed per-particle vector
        public double[,] array_particle;       // computed per-particle array
        public double[] vector_grid;           // computed per-grid vector
        public double[,] array_grid;           // computed per-grid array
        public double[] vector_surf;           // computed per-surf vector
        public double[,] array_surf;           // computed per-surf array
         
        public int START_OF_STEP, END_OF_STEP;    // mask settings
         
        public int kokkosable;                // 0/1 if Kokkos fix
        public int copymode;                 // 1 if copy of class (prevents deallocation of
                                       //  base class when child copy is destroyed)
        public ExecutionSpace execution_space;
        public uint datamask_read, datamask_modify;
        private SPARTA sparta;

        public Fix(SPARTA sparta, int narg, string[] arg)
        {
            this.sparta = sparta;
            // fix ID and style
            // ID must be all alphanumeric chars or underscores

            int n = arg[0].Length + 1;
            id = string.Copy( arg[0]);

            for (int i = 0; i < n - 1; i++)
                if (!char.IsLetterOrDigit(id[i]) && id[i] != '_')
                    sparta.error.all("Fix ID must be alphanumeric or underscore characters");

            n = arg[1].Length + 1;
            style = string.Copy(arg[1]);

            // set child class defaults

            time_depend = 0;
            gridmigrate = 0;
            flag_add_particle = flag_gas_react = flag_surf_react = 0;

            scalar_flag = vector_flag = array_flag = 0;
            per_particle_flag = per_grid_flag = per_surf_flag = 0;

            // mask settings - same as in modify.cpp

            START_OF_STEP = 1;
            END_OF_STEP = 2;

            kokkosable = 0;
            copymode = 0;

            execution_space = (int)ExecutionSpace.Host;
            datamask_read =SpartaMasks.ALL_MASK;
            datamask_modify =SpartaMasks.ALL_MASK;
        }


        public virtual int setmask()
        {
            return 0;
        }

        public virtual void init()
        {
            System.Console.WriteLine("fix virtual init");
        }
        public virtual void setup() { System.Console.WriteLine("fix virtual setup"); }

        public virtual void start_of_step() { System.Console.WriteLine("fix virtual start_of_step"); }
        //public virtual void end_of_step() { }
        //public virtual void add_particle(int, double, double, double, double[]) { }
        //public virtual void gas_react(int) { }
        //public virtual void surf_react(Particle::OnePart*, int &, int &) { }

        //public virtual void add_grid_one(int, int) { }
        public virtual int pack_grid_one(int icell,ref StringBuilder buf, int memflag) { return 0; }
        //public virtual int unpack_grid_one(int, char*) { return 0; }
        public virtual void compress_grid() { }
        public virtual void post_compress_grid() { }

        //public virtual double compute_scalar() { return 0.0; }
        //public virtual double compute_vector(int) { return 0.0; }
        //public virtual double compute_array(int, int) { return 0.0; }

        public virtual double memory_usage()
        {
            System.Console.WriteLine("Fix virtual memory_usage");
            return 0.0;
        }
    }
}