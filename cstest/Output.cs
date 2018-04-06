using bigint = System.Int64;

namespace cstest
{
    public class Output
    {
        bigint next;                 // next timestep for any kind of output

        bigint next_stats;           // next timestep for stats output
        int stats_every;             // stats output every this many steps
        bigint last_stats;           // last timestep stats was output
        string var_stats;             // variable name for stats frequency
        int ivar_stats;              // variable index for stats frequency
        Stats stats;          // statistical output

        int ndump;                   // # of Dumps defined
        int max_dump;                // max size of Dump list
        bigint next_dump_any;        // next timestep for any Dump
        int[] every_dump;             // output of each Dump every this many steps
        bigint[] next_dump;           // next timestep to do each Dump
        bigint[] last_dump;           // last timestep each snapshot was output
        string[] var_dump;             // variable name for dump frequency
        int[] ivar_dump;              // variable index for dump frequency
        Dump[] dump;           // list of defined Dumps   class Dump **dump

        int restart_flag;            // 1 if any restart files are written
        int restart_flag_single;     // 1 if single restart files are written
        int restart_flag_double;     // 1 if double restart files are written
        bigint next_restart;         // next timestep to write any restart file
        bigint next_restart_single;  // next timestep to write a single restart file
        bigint next_restart_double;  // next timestep to write a double restart file
        int restart_every_single;    // single restart file write freq, 0 if var
        int restart_every_double;    // double restart file write freq, 0 if var
        bigint last_restart;         // last timestep any restart file was output
        int restart_toggle;          // 0 if use restart2a as prefix, 1 if restart2b
        string var_restart_single;    // variable name for single restart freq
        string var_restart_double;    // variable name for double restart freq
        int ivar_restart_single;     // index of var_restart_single
        int ivar_restart_double;     // index of var_restart_double
        string restart1;              // name single restart file
        string restart2a;  // names of double restart files  char *restart2a,*restart2b
        string[] restart2b;

        WriteRestart restart; // class for writing restart files

        public Output(SPARTA sparta)
        {
            stats = new Stats(sparta);
        }
        //public void init();
        //public void setup(int);                   // initial output before run/min
        //public void write(bigint);                // output for current timestep
        //public void write_dump(bigint);           // force output of dump snapshots
        //public void write_restart(bigint);        // force output of a restart file
        //public void reset_timestep(bigint);       // reset next timestep for all output
         
        //public void add_dump(int, char**);       // add a Dump to Dump list
        //public void modify_dump(int, char**);    // modify a Dump
        //public void delete_dump(char*);          // delete a Dump from Dump list
        
        //public void set_stats(int, char**);      // set stats output frequency
        //public void create_stats(int, char**);   // create a Stats style
        //public void create_restart(int, char**); // create Restart and restart files
        
        //public void memory_usage();               // print out memory usage
    }
}