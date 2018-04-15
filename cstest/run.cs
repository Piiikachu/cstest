using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using bigint = System.Int64;

namespace cstest
{
    public class Run
    {
        private SPARTA sparta;

        public const int MAXSMALLINT = 2147483647;
        public const Int64 MAXBIGINT = 9223372036854775807;

        public Run(SPARTA sparta)
        {
            this.sparta = sparta;
        }

        public void command(int narg, string[] args)
        {
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);
            if (narg < 1) sparta.error.all("Illegal run command");

            if (sparta.grid.exist == 0)
                sparta.error.all("Run command before grid is defined");
            if (sparta.grid.exist_ghost == 0)
                sparta.error.all("Run command before grid ghost cells are defined");

            bigint nsteps_input = Int64.Parse(arg[0]);

            // parse optional args

            int uptoflag = 0;
            int startflag = 0;
            int stopflag = 0;
            bigint start = 0, stop = 0;
            int preflag = 1;
            int postflag = 1;
            int nevery = 0;
            int ncommands = 0;
            int first = 0, last = 0;

            int iarg = 1;

            while (iarg < narg)
            {
                switch (arg[iarg])
                {
                    case "upto":
                        if (iarg + 1 > narg) sparta.error.all("Illegal run command");
                        uptoflag = 1;
                        iarg += 1;
                        break;
                    case "start":
                        if (iarg + 2 > narg) sparta.error.all("Illegal run command");
                        startflag = 1;
                        start = Int64.Parse(arg[iarg + 1]);
                        iarg += 2;
                        break;
                    case "stop":
                        if (iarg + 2 > narg) sparta.error.all("Illegal run command");
                        stopflag = 1;
                        stop = Int64.Parse(arg[iarg + 1]);
                        iarg += 2;
                        break;
                    case "pre":
                        if (iarg + 2 > narg) sparta.error.all("Illegal run command");
                        if (string.Equals(arg[iarg + 1], "no")) preflag = 0;
                        else if (string.Equals(arg[iarg + 1], "yes")) preflag = 1;
                        else sparta.error.all("Illegal run command");
                        iarg += 2;
                        break;
                    case "post":
                        if (iarg + 2 > narg) sparta.error.all("Illegal run command");
                        if (string.Equals(arg[iarg + 1], "no")) postflag = 0;
                        else if (string.Equals(arg[iarg + 1], "yes")) postflag = 1;
                        else sparta.error.all("Illegal run command");
                        iarg += 2;

                        // all remaining args are commands
                        // first,last = arg index of first/last commands
                        // set ncommands = 0 if single command and it is NULL
                        break;
                    case "every":
                        if (iarg + 3 > narg) sparta.error.all("Illegal run command");
                        nevery = int.Parse(arg[iarg + 1]);
                        if (nevery <= 0) sparta.error.all("Illegal run command");
                        first = iarg + 2;
                        last = narg - 1;
                        ncommands = last - first + 1;
                        if (ncommands == 1 && string.Equals(arg[first], "NULL")) ncommands = 0;
                        iarg = narg;
                        break;
                    default:
                        sparta.error.all("Illegal run command");
                        break;
                }





            }
            int nsteps;
            if (uptoflag == 0)
            {
                if (nsteps_input < 0 || nsteps_input > MAXSMALLINT)
                    sparta.error.all("Invalid run command N value");
                nsteps = Convert.ToInt32(nsteps_input);
            }
            else
            {
                bigint delta = nsteps_input - sparta.update.ntimestep;
                if (delta < 0 || delta > MAXSMALLINT)
                    sparta.error.all("Invalid run command upto value");
                nsteps = Convert.ToInt32(delta);
            }

            // error check

            if (startflag != 0)
            {
                if (start < 0 || start > MAXBIGINT)
                    sparta.error.all("Invalid run command start/stop value");
                if (start > sparta.update.ntimestep)
                    sparta.error.all("Run command start value is after start of run");
            }
            if (stopflag != 0)
            {
                if (stop < 0 || stop > MAXBIGINT)
                    sparta.error.all("Invalid run command start/stop value");
                if (stop < sparta.update.ntimestep + nsteps)
                    sparta.error.all("Run command stop value is before end of run");
            }

            // if nevery, make copies of arg strings that are commands
            // required because re-parsing commands via input->one() will wipe out args

            string[] commands = null;
            if (nevery != 0 && ncommands > 0)
            {
                commands = new string[ncommands];
                ncommands = 0;
                for (int i = first; i <= last; i++)
                {
                    int n = arg[i].Length + 1;
                    commands[ncommands] = string.Copy(arg[i]);
                    ncommands++;
                }
            }

            // perform a single run
            // use start/stop to set begin/end step
            // if pre or 1st run, do System init/setup,
            //   else just init timer and setup output
            // if post, do full Finish, else just print time

            sparta.update.runflag = 1;

            if (nevery == 0)
            {
                sparta.update.nsteps = nsteps;
                sparta.update.firststep = sparta.update.ntimestep;
                sparta.update.laststep = sparta.update.ntimestep + nsteps;
                if (sparta.update.laststep < 0 || sparta.update.laststep > MAXBIGINT)
                    sparta.error.all("Too many timesteps");

                if (startflag != 0) sparta.update.beginstep = start;
                else sparta.update.beginstep = sparta.update.firststep;
                if (stopflag != 0) sparta.update.endstep = stop;
                else sparta.update.endstep = sparta.update.laststep;

                if (preflag != 0 || sparta.update.first_update == 0)
                {
                    sparta.init();
                    //sparta.update.setup();
                }
                //else sparta.output.setup(0);


            }
        }
    }
}
