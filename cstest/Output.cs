﻿using System;
using System.IO;
using bigint = System.Int64;

namespace cstest
{
    public class Output
    {
        public const int DELTA = 1;
        public bigint next;                 // next timestep for any kind of output
         
        public bigint next_stats;           // next timestep for stats output
        public int stats_every;             // stats output every this many steps
        public bigint last_stats;           // last timestep stats was output
        public string var_stats;             // variable name for stats frequency
        public int ivar_stats;              // variable index for stats frequency
        public Stats stats;          // statistical output
         
        public int ndump;                   // # of Dumps defined
        public int max_dump;                // max size of Dump list
        public bigint next_dump_any;        // next timestep for any Dump
        public int[] every_dump;             // output of each Dump every this many steps
        public bigint[] next_dump;           // next timestep to do each Dump
        public bigint[] last_dump;           // last timestep each snapshot was output
        public string[] var_dump;             // variable name for dump frequency
        public int[] ivar_dump;              // variable index for dump frequency
        public Dump[] dump;           // list of defined Dumps   class Dump **dump
         
        public int restart_flag;            // 1 if any restart files are written
        public int restart_flag_single;     // 1 if single restart files are written
        public int restart_flag_double;     // 1 if double restart files are written
        public bigint next_restart;         // next timestep to write any restart file
        public bigint next_restart_single;  // next timestep to write a single restart file
        public bigint next_restart_double;  // next timestep to write a double restart file
        public int restart_every_single;    // single restart file write freq, 0 if var
        public int restart_every_double;    // double restart file write freq, 0 if var
        public bigint last_restart;         // last timestep any restart file was output
        public int restart_toggle;          // 0 if use restart2a as prefix, 1 if restart2b
        public string var_restart_single;    // variable name for single restart freq
        public string var_restart_double;    // variable name for double restart freq
        public int ivar_restart_single;     // index of var_restart_single
        public int ivar_restart_double;     // index of var_restart_double
        public string restart1;              // name single restart file
        public string restart2a,restart2b;  // names of double restart files  char *restart2a,*restart2b
         
        public WriteRestart restart; // class for writing restart files

        private SPARTA sparta;

        public Output(SPARTA sparta)
        {
            this.sparta = sparta;
            stats = new Stats(sparta);

            stats_every = 0;
            var_stats = null;

            ndump = 0;
            max_dump = 0;
            every_dump = null;
            next_dump = null;
            last_dump = null;
            var_dump = null;
            ivar_dump = null;
            dump = null;

            restart_flag = restart_flag_single = restart_flag_double = 0;
            restart_every_single = restart_every_double = 0;
            last_restart = -1;
            restart1 = restart2a = restart2b = null;
            var_restart_single = var_restart_double = null;
            restart = null;

        }
        public void init()
        {
            stats.init();
            if (var_stats != null)
            {
                ivar_stats = sparta.input.variable.find(var_stats);
                if (ivar_stats < 0)
                    sparta.error.all("Variable name for stats every does not exist");
                if (sparta.input.variable.equal_style(ivar_stats)==0)
                    sparta.error.all("Variable for stats every is invalid style");
            }

            for (int i = 0; i < ndump; i++) dump[i].init();
            for (int i = 0; i < ndump; i++)
                if (every_dump[i] == 0)
                {
                    ivar_dump[i] = sparta.input.variable.find(var_dump[i]);
                    if (ivar_dump[i] < 0)
                        sparta.error.all("Variable name for dump every does not exist");
                    if (sparta.input.variable.equal_style(ivar_dump[i])==0)
                        sparta.error.all("Variable for dump every is invalid style");
                }

            if (restart_flag_single != 0 && restart_every_single == 0)
            {
                ivar_restart_single = sparta.input.variable.find(var_restart_single);
                if (ivar_restart_single < 0)
                    sparta.error.all("Variable name for restart does not exist");
                if (sparta.input.variable.equal_style(ivar_restart_single)==0)
                    sparta.error.all("Variable for restart is invalid style");
            }
            if (restart_flag_double != 0 && restart_every_double == 0)
            {
                ivar_restart_double = sparta.input.variable.find(var_restart_double);
                if (ivar_restart_double < 0)
                    sparta.error.all("Variable name for restart does not exist");
                if (sparta.input.variable.equal_style(ivar_restart_double)==0)
                    sparta.error.all("Variable for restart is invalid style");
            }
        }
        public void setup(int memflag)                   // initial output before run/min
        {
            bigint ntimestep = sparta.update.ntimestep;

            // perform dump at start of run only if:
            //   current timestep is multiple of every and last dump not >= this step
            //   this is first run after dump created and firstflag is set
            //   note that variable freq will not write unless triggered by firstflag
            // set next_dump to multiple of every or variable value
            // set next_dump_any to smallest next_dump
            // wrap dumps that invoke computes and variable eval with clear/add
            // if dump not written now, use addstep_compute_all() since don't know
            //   what computes the dump write would invoke
            // if no dumps, set next_dump_any to last+1 so will not influence next

            int writeflag;
            if (ndump!=0)
            {
                for (int idump = 0; idump < ndump; idump++)
                {
                    if (dump[idump].clearstep!=0 || every_dump[idump] == 0)
                        sparta.modify.clearstep_compute();
                    writeflag = 0;
                    if (every_dump[idump] != 0 && ntimestep % every_dump[idump] == 0 &&
                        last_dump[idump] != ntimestep) writeflag = 1;
                    if (last_dump[idump] < 0 && dump[idump].first_flag == 1) writeflag = 1;

                    if (writeflag != 0)
                    {
                        dump[idump].write();
                        last_dump[idump] = ntimestep;
                    }
                    if (every_dump[idump] != 0)
                        next_dump[idump] =
                          (ntimestep / every_dump[idump]) * every_dump[idump] + every_dump[idump];
                    else
                    {
                        bigint nextdump = Convert.ToInt64
                          (sparta.input.variable.compute_equal(ivar_dump[idump]));
                        if (nextdump <= ntimestep)
                            sparta.error.all("Dump every variable returned a bad timestep");
                        next_dump[idump] = nextdump;
                    }
                    if (dump[idump].clearstep!=0 || every_dump[idump] == 0)
                    {
                        if (writeflag!=0) sparta.modify.addstep_compute(next_dump[idump]);
                        else sparta.modify.addstep_compute_all(next_dump[idump]);
                    }
                    if (idump!=0) next_dump_any = Math.Min(next_dump_any, next_dump[idump]);
                    else next_dump_any = next_dump[0];
                }
            }
            else next_dump_any = sparta.update.laststep + 1;

            // do not write restart files at start of run
            // set next_restart values to multiple of every or variable value
            // wrap variable eval with clear/add
            // if no restarts, set next_restart to last+1 so will not influence next

            if (restart_flag!=0)
            {
                if (restart_flag_single != 0)
                {
                    if (restart_every_single != 0)
                        next_restart_single =
                          (ntimestep / restart_every_single) * restart_every_single +
                          restart_every_single;
                    else
                    {
                        bigint nextrestart = Convert.ToInt64
                          (sparta.input.variable.compute_equal(ivar_restart_single));
                        if (nextrestart <= ntimestep)
                            sparta.error.all("Restart variable returned a bad timestep");
                        next_restart_single = nextrestart;
                    }
                }
                else next_restart_single = sparta.update.laststep + 1;
                if (restart_flag_double != 0)
                {
                    if (restart_every_double != 0)
                        next_restart_double =
                          (ntimestep / restart_every_double) * restart_every_double +
                          restart_every_double;
                    else
                    {
                        bigint nextrestart = Convert.ToInt64
                          (sparta.input.variable.compute_equal(ivar_restart_double));
                        if (nextrestart <= ntimestep)
                            sparta.error.all("Restart variable returned a bad timestep");
                        next_restart_double = nextrestart;
                    }
                }
                else next_restart_double = sparta.update.laststep + 1;
                next_restart = Math.Min(next_restart_single, next_restart_double);
            }
            else next_restart = sparta.update.laststep + 1;

            // print memory usage unless being called between multiple runs
            if (memflag != 0) memory_usage();

            // set next_stats to multiple of every or variable eval if var defined
            // insure stats output on last step of run
            // stats may invoke computes so wrap with clear/add

            sparta.modify.clearstep_compute();

            stats.header();
            stats.compute(0);
            last_stats = ntimestep;

            if (var_stats!=null)
            {
                next_stats = Convert.ToInt64
                  (sparta.input.variable.compute_equal(ivar_stats));
                if (next_stats <= ntimestep)
                    sparta.error.all("Stats every variable returned a bad timestep");
            }
            else if (stats_every!=0)
            {
                next_stats = (ntimestep / stats_every) * stats_every + stats_every;
                next_stats = Math.Min(next_stats, sparta.update.laststep);
            }
            else next_stats = sparta.update.laststep;

            sparta.modify.addstep_compute(next_stats);

            // next = next timestep any output will be done

            next = Math.Min(next_dump_any, next_restart);
            next = Math.Min(next, next_stats);
        }
        public void write(bigint ntimestep)                // output for current timestep
        {
            // next_dump does not force output on last step of run
            // wrap dumps that invoke computes or eval of variable with clear/add
            // download data from GPU if necessary

            if (next_dump_any == ntimestep)
            {
                for (int idump = 0; idump < ndump; idump++)
                {
                    if (next_dump[idump] == ntimestep)
                    {
                        if (dump[idump].clearstep !=0|| every_dump[idump] == 0)
                            sparta.modify.clearstep_compute();
                        if (last_dump[idump] != ntimestep)
                        {
                            dump[idump].write();
                            last_dump[idump] = ntimestep;
                        }
                        if (every_dump[idump]!=0) next_dump[idump] += every_dump[idump];
                        else
                        {
                            bigint nextdump = Convert.ToInt64
                              (sparta.input.variable.compute_equal(ivar_dump[idump]));
                            if (nextdump <= ntimestep)
                                sparta.error.all("Dump every variable returned a bad timestep");
                            next_dump[idump] = nextdump;
                        }
                        if (dump[idump].clearstep != 0 || every_dump[idump] == 0)
                            sparta.modify.addstep_compute(next_dump[idump]);
                    }
                    if (idump!=0) next_dump_any = Math.Min(next_dump_any, next_dump[idump]);
                    else next_dump_any = next_dump[0];
                }
            }

            // next_restart does not force output on last step of run
            // for toggle = 0, replace "*" with current timestep in restart filename
            // eval of variable may invoke computes so wrap with clear/add

            if (next_restart == ntimestep)
            {
                if (next_restart_single == ntimestep)
                {
                    Console.WriteLine("Output.write()->next_restart_single == ntimestep");
                    //char* file = new char[strlen(restart1) + 16];
                    //char* ptr = strchr(restart1, '*');
                    //*ptr = '\0';
                    //sprintf(file, "%s" BIGINT_FORMAT "%s", restart1, ntimestep, ptr + 1);
                    //*ptr = '*';
                    //if (last_restart != ntimestep) restart.write(file);
                    //delete[] file;
                    //if (restart_every_single!=0) next_restart_single += restart_every_single;
                    //else
                    //{
                    //    sparta.modify.clearstep_compute();
                    //    bigint nextrestart = Convert.ToInt64
                    //      (sparta.input.variable.compute_equal(ivar_restart_single));
                    //    if (nextrestart <= ntimestep)
                    //        sparta.error.all("Restart variable returned a bad timestep");
                    //    next_restart_single = nextrestart;
                    //    sparta.modify.addstep_compute(next_restart_single);
                    //}
                }
                if (next_restart_double == ntimestep)
                {
                    if (last_restart != ntimestep)
                    {
                        if (restart_toggle == 0)
                        {
                            restart.write(restart2a);
                            restart_toggle = 1;
                        }
                        else
                        {
                            restart.write(restart2b);
                            restart_toggle = 0;
                        }
                    }
                    if (restart_every_double!=0) next_restart_double += restart_every_double;
                    else
                    {
                        sparta.modify.clearstep_compute();
                        bigint nextrestart = Convert.ToInt64
                          (sparta.input.variable.compute_equal(ivar_restart_double));
                        if (nextrestart <= ntimestep)
                            sparta.error.all("Restart variable returned a bad timestep");
                        next_restart_double = nextrestart;
                        sparta.modify.addstep_compute(next_restart_double);
                    }
                }
                last_restart = ntimestep;
                next_restart = Math.Min(next_restart_single, next_restart_double);
            }

            // insure next_thermo forces output on last step of run
            // thermo may invoke computes so wrap with clear/add

            if (next_stats == ntimestep)
            {
                sparta.modify.clearstep_compute();
                if (last_stats != ntimestep) stats.compute(1);
                last_stats = ntimestep;
                if (var_stats != null)
                {
                    next_stats = Convert.ToInt64
                      (sparta.input.variable.compute_equal(ivar_stats));
                    if (next_stats <= ntimestep)
                        sparta.error.all("Stats every variable returned a bad timestep");
                }
                else if (stats_every!=0) next_stats += stats_every;
                else next_stats = sparta.update.laststep;
                next_stats = Math.Min(next_stats, sparta.update.laststep);
                sparta.modify.addstep_compute(next_stats);
            }

            // next = next timestep any output will be done

            next = Math.Min(next_dump_any, next_restart);
            next = Math.Min(next, next_stats);
        }
        //public void write_dump(bigint);           // force output of dump snapshots
        //public void write_restart(bigint);        // force output of a restart file
        //public void reset_timestep(bigint);       // reset next timestep for all output

        //public void add_dump(int, char**);       // add a Dump to Dump list
        //public void modify_dump(int, char**);    // modify a Dump
        //public void delete_dump(char*);          // delete a Dump from Dump list

        public void set_stats(int narg, string[] arg)      // set stats output frequency
        {
            if (narg != 1) sparta.error.all("Illegal stats command");

            if (arg[1].Contains("v_"))
            {
                //delete[] var_stats;
                //int n = strlen(&arg[0][2]) + 1;
                //var_stats = new char[n];
                //strcpy(var_stats, &arg[0][2]);
                System.Console.WriteLine("output.set_stats().starts with v_");
            }
            else
            {
                stats_every = int.Parse(arg[1]);
                if (stats_every < 0) sparta.error.all("Illegal stats command");
            }
        }
        public void create_stats(int narg,string[] args)   // create a Stats style
        {
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);
            if (narg < 1) sparta.error.all("Illegal stats_style command");
            stats.set_fields(narg, arg);





        }
        //public void create_restart(int, char**); // create Restart and restart files

        public void memory_usage()               // print out memory usage
        {
            bigint pbytes, gbytes, sbytes, bytes;
            pbytes = sparta.particle.memory_usage();
            gbytes = sparta.grid.memory_usage();
            sbytes = sparta.surf.memory_usage();
            bytes = pbytes + gbytes + sbytes;
            bytes += sparta.modify.memory_usage();

            double scale = 1.0 / 1024.0 / 1024.0;

            bigint ave=0, min=0, max=0;

            sparta.mpi.MPI_Allreduce(ref pbytes, ref ave, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            double pave = scale * ave / sparta.comm.nprocs;
            sparta.mpi.MPI_Allreduce(ref pbytes, ref min, 1, MPI.MPI_LONG_LONG, MPI.MPI_MIN, sparta.world);
            double pmin = scale * min;
            sparta.mpi.MPI_Allreduce(ref pbytes, ref max, 1, MPI.MPI_LONG_LONG, MPI.MPI_MAX, sparta.world);
            double pmax = scale * max;

            sparta.mpi.MPI_Allreduce(ref gbytes, ref ave, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            double gave = scale * ave / sparta.comm.nprocs;
            sparta.mpi.MPI_Allreduce(ref gbytes, ref min, 1, MPI.MPI_LONG_LONG, MPI.MPI_MIN, sparta.world);
            double gmin = scale * min;
            sparta.mpi.MPI_Allreduce(ref gbytes, ref max, 1, MPI.MPI_LONG_LONG, MPI.MPI_MAX, sparta.world);
            double gmax = scale * max;

            sparta.mpi.MPI_Allreduce(ref sbytes, ref ave, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            double save = scale * ave / sparta.comm.nprocs;
            sparta.mpi.MPI_Allreduce(ref sbytes, ref min, 1, MPI.MPI_LONG_LONG, MPI.MPI_MIN, sparta.world);
            double smin = scale * min;
            sparta.mpi.MPI_Allreduce(ref sbytes, ref max, 1, MPI.MPI_LONG_LONG, MPI.MPI_MAX, sparta.world);
            double smax = scale * max;

            sparta.mpi.MPI_Allreduce(ref bytes, ref ave, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            double tave = scale * ave / sparta.comm.nprocs;
            sparta.mpi.MPI_Allreduce(ref bytes, ref min, 1, MPI.MPI_LONG_LONG, MPI.MPI_MIN, sparta.world);
            double tmin = scale * min;
            sparta.mpi.MPI_Allreduce(ref bytes, ref max, 1, MPI.MPI_LONG_LONG, MPI.MPI_MAX, sparta.world);
            double tmax = scale * max;

            if (sparta.comm.me == 0)
            {
                string str1 = string.Format("Memory usage per proc in Mbytes:\n");
                string str2 = string.Format("  particles (ave,min,max) = {0} {1} {2}\n",
                        pave, pmin, pmax);
                string str3 = string.Format("  grid      (ave,min,max) = {0} {1} {2}\n",
                        gave, gmin, gmax);
                string str4 = string.Format("  surf      (ave,min,max) = {0} {1} {2}\n",
                        save, smin, smax);
                string str5 = string.Format("  total     (ave,min,max) = {0} {1} {2}\n",
                        tave, tmin, tmax);
                Console.WriteLine(str1 + str2 + str3 + str4 + str5);
                if (sparta.screen!=null)
                {
                    new StreamWriter(sparta.screen).WriteLine(str1+str2+str3+str4+str5);
                }
                if (sparta.logfile!=null)
                {
                    new StreamWriter(sparta.logfile).WriteLine(str1 + str2 + str3 + str4 + str5);
                }
            }
        }
    }
}