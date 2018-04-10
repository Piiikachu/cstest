using System;
using System.Collections.Generic;
using System.IO;

namespace cstest
{
    public class Input
    {
        public const int DELTALINE = 256;
        public const int DELTA = 4;
        public int narg;                    // # of command args
        public string[] arg;                  // parsed args for command
        public Variable variable;    // defined variables

        private int me;                      // proc ID
        private string command;               // ptr to current command
        private int maxarg;                  // max # of args in arg
        private string line,copy,work;      // input line & copy and work string
        private int maxline, maxcopy, maxwork; // max lengths of char strings
        private int echo_screen;             // 0 = no, 1 = yes
        private int echo_log;                // 0 = no, 1 = yes
        private int nfile, maxfile;           // current # and max # of open input files
        private int label_active;            // 0 = no label, 1 = looking for label
        private string labelstr;              // label string being looked for
        private int jump_skip;               // 1 if skipping next jump, 0 otherwise
        private int ifthenelse_flag;         // 1 if executing commands inside an if-then-else

        private List<FileStream> infiles;              // list of open input files

        private SPARTA sparta;


        public Input(SPARTA sparta, string[] args)
        {
            arg = args;
            int argc = args.Length;
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);

            maxline = maxcopy = maxwork = 0;
            line = copy = work = null;
            narg = maxarg = 0;
            arg = null;

            echo_screen = 0;
            echo_log = 1;

            label_active = 0;
            labelstr = null;
            jump_skip = 0;
            ifthenelse_flag = 0;

            if (me == 0)
            {
                nfile = maxfile = 1;
                //infiles = (FILE**)memory.smalloc(sizeof(FILE*), "input:infiles");
                infiles = new List<FileStream>();
                infiles.Add(sparta.infile);
            }
            else infiles = null;

            variable = new Variable(sparta);

            // process command-line args
            // check for args "-var" and "-echo"
            // caller has already checked that sufficient arguments exist

            int iarg = 0;
            while (iarg <argc )
            {
                if (string.Compare(args[iarg], "-var") == 0 || string.Compare(args[iarg], "-v") == 0)
                {
                    int jarg = iarg + 3;
                    while (jarg < argc && args[jarg][0] != '-') jarg++;
                    //variable.set(args[iarg + 1], jarg - iarg - 2, ref args[iarg + 2]);
                    iarg = jarg;
                }
                else if (string.Compare(args[iarg], "-echo") == 0 ||
                         string.Compare(args[iarg], "-e") == 0)
                {
                    narg = 1;
                    string[] tmp = arg;        // trick echo() into using args instead of arg
                   // arg = &args[iarg + 1];
                    echo();
                    arg = tmp;
                    iarg += 2;
                }
                else iarg++;
            }
        }
        public void file()                    // process all input
        {
            System.Console.WriteLine("input.file()");
            int m, n=0;
            using (StreamReader sr=new StreamReader(sparta.infile))
            {
                string line;
                while ((line=sr.ReadLine())!=null)
                {
                    Parseline(line);
                    //System.Console.WriteLine(line);
                }



            }


            //while (true)
            //{

            //    // read a line from input script
            //    // n = length of line including str terminator, 0 if end of file
            //    // if line ends in continuation char '&', concatenate next line

            //    if (me == 0)
            //    {
            //        m = 0;
            //        while (true)
            //        {
            //            if (maxline - m < 2) reallocate(ref line,ref maxline, 0);
            //            if (fgets(&line[m], maxline - m, infile) == null)
            //            {
            //                if (m) n = strlen(line) + 1;
            //                else n = 0;
            //                break;
            //            }
            //            m = strlen(line);
            //            if (line[m - 1] != '\n') continue;

            //            m--;
            //            while (m >= 0 && isspace(line[m])) m--;
            //            if (m < 0 || line[m] != '&')
            //            {
            //                line[m + 1] = '\0';
            //                n = m + 2;
            //                break;
            //            }
            //        }
            //    }

            //    // bcast the line
            //    // if n = 0, end-of-file
            //    // error if label_active is set, since label wasn't encountered
            //    // if original input file, code is done
            //    // else go back to previous input file

            //    sparta.mpi.MPI_Bcast( n, 1,MPI.MPI_INT, 0, sparta.world);
            //    if (n == 0)
            //    {
            //        if (label_active) sparta.error.all( "Label wasn't found in input script");
            //        if (me == 0)
            //        {
            //            if (infile != stdin)
            //            {
            //                fclose(infile);
            //                infile = null;
            //            }
            //            nfile--;
            //        }
            //        MPI_Bcast(&nfile, 1, MPI_INT, 0, world);
            //        if (nfile == 0) break;
            //        if (me == 0) infile = infiles[nfile - 1];
            //        continue;
            //    }

            //    if (n > maxline) reallocate(line, maxline, n);
            //    MPI_Bcast(line, n, MPI_CHAR, 0, world);

            //    // echo the command unless scanning for label

            //    if (me == 0 && label_active == 0)
            //    {
            //        if (echo_screen && screen) fprintf(screen, "%s\n", line);
            //        if (echo_log && logfile) fprintf(logfile, "%s\n", line);
            //    }

            //    // parse the line
            //    // if no command, skip to next line in input script

            //    parse();
            //    if (command == null) continue;

            //    // if scanning for label, skip command unless it's a label command

            //    if (label_active && string.Compare(command, "label") != 0) continue;

            //    // execute the command

            //    if (execute_command())
            //    {
            //        char* str = new char[maxline + 32];
            //        sprintf(str, "Unknown command: %s", line);
            //        sparta.error.all( str);
            //    }
            //}
        }                  
        public void file(string filename)       // process an input script
        {
            System.Console.WriteLine("input.file({0})",filename);
            // error if another nested file still open, should not be possible
            // open new filename and set infile, infiles[0], nfile
            // call to file() will close filename and decrement nfile

            //if (me == 0)
            //{
            //    if (nfile > 1)
            //        error->one(FLERR, "Invalid use of library file() function");

            //    if (infile && infile != stdin) fclose(infile);
            //    infile = fopen(filename, "r");
            //    if (infile == null)
            //    {
            //        char str[128];
            //        sprintf(str, "Cannot open input script %s", filename);
            //        error->one(FLERR, str);
            //    }
            //    infiles[0] = infile;
            //    nfile = 1;
            //}

            //file();
        }
        public string one(string single)       // process a single command
        {
            int n = single.Length + 1;
            if (n > maxline) reallocate(ref line,ref maxline, n);
           line=string.Copy(single);

            // echo the command unless scanning for label

            if (me == 0 && label_active == 0)
            {
                if (echo_screen!=0 && sparta.screen!=null)
                {
                    SPARTA.fprintf(sparta.screen, "{0}\n", line);
                }

                if (echo_log!=0 && sparta.logfile!=null)
                {
                    SPARTA.fprintf(sparta.logfile, "{0}\n", line);
                }
            }

            // parse the line
            // if no command, just return null

            parse();
            if (command == null) return null;

            // if scanning for label, skip command unless it's a label command

            if (label_active!=0 && string.Compare(command, "label") != 0) return null;

            // execute the command and return its name

            if (execute_command())
            {
                string str = string.Format("Unknown command: {0}", line);
                //char* str = new char[maxline + 32];
                //sprintf(str, "Unknown command: %s", line);
                sparta.error.all(str);
            }

            return command;
        }
        //public void substitute(char*&, char*&, int &, int &, int);
        //// substitute for variables in a string
        public int expand_args(int narg, string[] arg, int mode, ref string[] earg)
        {
            return 0;
        }  // expand args due to wildcard

        public double numeric(string str)    // arg checking
        {
            if (string.IsNullOrEmpty(str))
                sparta.error.all("Expected floating point parameter in input script or data file");
            int n = str.Length;
            if (n == 0)
                sparta.error.all("Expected floating point parameter in input script or data file");

            for (int i = 0; i < n; i++)
            {
                if (char.IsDigit(str[i])) continue;
                if (str[i] == '-' || str[i] == '+' || str[i] == '.') continue;
                if (str[i] == 'e' || str[i] == 'E') continue;
                sparta.error.all("Expected floating point parameter in input script or data file");
            }

            return double.Parse(str);
        }
        public int inumeric(string str)
        {
            if (string.IsNullOrWhiteSpace(str))
               sparta.error.all("Expected integer parameter in input script or data file");
            int n = str.Length;
            if (n == 0)
                sparta.error.all("Expected integer parameter in input script or data file");

            for (int i = 0; i < n; i++)
            {
                if (char.IsDigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
                sparta.error.all("Expected integer parameter in input script or data file");
            }

            return int.Parse(str);
        }
        //public bigint bnumeric(const char*, int, char*);
        //public void bounds(char*, int, int &, int &, int nmin = 1);
        //public int count_words(char*);

        ////private:

        /* ----------------------------------------------------------------------
        parse copy of command line by inserting string terminators
        strip comment = all chars from # on
        replace all $ via variable substitution
        command = first word
        narg = # of args
        arg[] = individual args
        treat text between single/double quotes as one arg
        ------------------------------------------------------------------------- */
        private void Parseline(string line)
        {
            if (string.IsNullOrWhiteSpace(line))
            {
                return;
            }
            if (line.StartsWith("#"))
            {
                line.Remove(0);
                return;
            }



            
            string[] tempstr = line.Split();
            List<string> strings=new List<string>();
            foreach (string str in tempstr)
            {
                if (!string.IsNullOrEmpty(str))
                {
                    if (str.StartsWith("#"))
                    {
                        break;
                    }
                    strings.Add(str);
                }
            }
            string[] strArray = strings.ToArray();
            narg = strArray.Length - 1;
            arg = strArray;
            System.Console.Write("\n"+strArray[0]+"...........");
            execute_command(strings[0]);
            


        }

        private void parse()                          // parse an input text line
        {
            // duplicate line into copy string to break into words

            int n = line.Length + 1;
            if (n > maxcopy) copy = string.Copy(line);
            

            // strip any # comment by replacing it with 0
            // do not strip # inside single/double quotes

            //char quote = '\0';
            if (copy.StartsWith("#"))
            {
                copy.Remove(0);
            }


            // perform $ variable substitution (print changes)
            // except if searching for a label since earlier variable may not be defined

            //if (label_active==0) substitute(copy, work, maxcopy, maxwork, 1);

            //// command = 1st arg in copy string
            
            //string next="";
            //command = nextword(copy, ref next);
            //if (command == null) return;

            // point arg[] at each subsequent arg in copy string
            // nextword() inserts string terminators into copy string to delimit args
            // nextword() treats text between single/double quotes as one arg

            //narg = 0;
            //ptr = next;
            //while (ptr)
            //{
            //    if (narg == maxarg)
            //    {
            //        maxarg += DELTA;
            //        arg = (char**)memory->srealloc(arg, maxarg * sizeof(char*), "input:arg");
            //    }
            //    arg[narg] = nextword(ptr, &next);
            //    if (!arg[narg]) break;
            //    narg++;
            //    ptr = next;
            //}
        }
        //private char* nextword(char*, char**);       // find next word in string with quotes
        private void reallocate(ref string str,ref int max, int n)  // reallocate a char string
        {
            if (n!=0)
            {
                while (n > max) max += DELTALINE;
            }
            else max += DELTALINE;

            //str = (char*)memory->srealloc(str, max * sizeof(char), "input:str");
            str = "input.reallocate";
            System.Console.WriteLine(str);
        }
        private bool execute_command(string command)                 // execute a single command
        {
            //todo: why not switch
            bool flag = true;
            
            if (string.Equals(command, "clear")) clear();
            else if (string.Equals(command, "echo")) echo();
            else if (string.Equals(command, "if")) ifthenelse();
            else if (string.Equals(command, "include")) include();
            else if (string.Equals(command, "jump")) jump();
            else if (string.Equals(command, "label")) label();
            else if (string.Equals(command, "log")) log();
            else if (string.Equals(command, "next")) next_command();
            else if (string.Equals(command, "partition")) partition();
            else if (string.Equals(command, "print")) print();
            else if (string.Equals(command, "quit")) quit();
            else if (string.Equals(command, "shell")) shell();
            else if (string.Equals(command, "variable")) variable_command();

            else if (string.Equals(command, "boundary")) boundary();
            else if (string.Equals(command, "bound_modify")) bound_modify();
            else if (string.Equals(command, "balance_grid")) balance_grid();
            else if (string.Equals(command, "create_box")) create_box();
            else if (string.Equals(command, "create_grid")) create_grid();
            else if (string.Equals(command, "collide")) collide_command();
            else if (string.Equals(command, "collide_modify")) collide_modify();
            else if (string.Equals(command, "compute")) compute();
            else if (string.Equals(command, "dimension")) dimension();
            else if (string.Equals(command, "dump")) dump();
            else if (string.Equals(command, "dump_modify")) dump_modify();
            else if (string.Equals(command, "fix")) fix();
            else if (string.Equals(command, "global")) global();
            else if (string.Equals(command, "group")) group();
            else if (string.Equals(command, "package")) package();
            else if (string.Equals(command, "mixture")) mixture();
            else if (string.Equals(command, "react")) react_command();
            else if (string.Equals(command, "react_modify")) react_modify();
            else if (string.Equals(command, "read_surf")) read_surf();
            else if (string.Equals(command, "region")) region();
            else if (string.Equals(command, "reset_timestep")) reset_timestep();
            else if (string.Equals(command, "restart")) restart();
            else if (string.Equals(command, "seed")) seed();
            else if (string.Equals(command, "species")) species();
            else if (string.Equals(command, "stats")) stats();
            else if (string.Equals(command, "stats_modify")) stats_modify();
            else if (string.Equals(command, "stats_style")) stats_style();
            else if (string.Equals(command, "surf_collide")) surf_collide();
            else if (string.Equals(command, "surf_modify")) surf_modify();
            else if (string.Equals(command, "surf_react")) surf_react();
            else if (string.Equals(command, "timestep")) timestep();
            else if (string.Equals(command, "uncompute")) uncompute();
            else if (string.Equals(command, "undump")) undump();
            else if (string.Equals(command, "unfix")) unfix();
            else if (string.Equals(command, "units")) units();

            else
            {
                flag = false;
                sparta.error.universe_all("unknown command1");
                //System.Console.WriteLine("unknown command1");
            }

            if (sparta.suffix_enable!=0&&sparta.suffix!=null)
            {
                string command2 = string.Format("{0}/{1}", command, sparta.suffix);
                Console.WriteLine(command2);




            }


            return !flag;
            

            // return if command was listed above
        }

        private void read_surf()
        {
            new ReadSurf(sparta).command(narg, arg);
            Console.WriteLine("done");
        }

        private void balance_grid()
        {
            new BalanceGrid(sparta).command(narg, arg);
            Console.WriteLine("done");
        }

        private void create_grid()
        {
            new CreateGrid(sparta).command(narg, arg);
            Console.WriteLine("done");
        }

        private void create_box()
        {
            new CreateBox(sparta).command(narg, arg);
            Console.WriteLine("done");

        }

        private bool execute_command()
        {
            return execute_command(command);
        }


        private void clear()            // input script commands
        {
            System.Console.WriteLine("input.clear()");

            if (narg > 0)
            {
                sparta.error.all("Illegal clear command");
            }
            sparta.create();
            sparta.post_create();
        }                
        private void echo()
        {
            
            System.Console.WriteLine("input.echo()");
            if (narg != 1) sparta.error.all( "Illegal echo command");

            if (string.Compare(arg[0], "none") == 0)
            {
                echo_screen = 0;
                echo_log = 0;
            }
            else if (string.Compare(arg[0], "screen") == 0)
            {
                echo_screen = 1;
                echo_log = 0;
            }
            else if (string.Compare(arg[0], "log") == 0)
            {
                echo_screen = 0;
                echo_log = 1;
            }
            else if (string.Compare(arg[0], "both") == 0)
            {
                echo_screen = 1;
                echo_log = 1;
            }
            else sparta.error.all("Illegal echo command");
        }
        private void ifthenelse(){}
        private void include()
        {
            if (narg != 1) sparta.error.all("Illegal include command");

            // do not allow include inside an if command
            // NOTE: this check will fail if a 2nd if command was inside the if command
            //       and came before the include

            if (ifthenelse_flag!=0)
                sparta.error.all("Cannot use include command within an if command");

            if (me == 0)
            {
                if (nfile == maxfile)
                {
                    maxfile++;
                    //infiles = (FILE**)
                    //  memory->srealloc(infiles, maxfile * sizeof(FILE*), "input:infiles");
                    infiles = new List<FileStream>();
                }
                sparta.infile = new FileStream(arg[0],FileMode.Open, FileAccess.Read);
                if (sparta.infile == null)
                {
                    //char str[128];
                    //sprintf(str, "Cannot open input script %s", arg[0]);
                    //error->one(FLERR, str);
                    string str = string.Format("Cannot open input script {0}", arg[0]);
                    sparta.error.one(str);

                }
                infiles[nfile++] = sparta.infile;
            }
        }
        private void jump()
        {
            if (narg < 1 || narg > 2) sparta.error.all( "Illegal jump command");

            if (jump_skip!=0)
            {
                jump_skip = 0;
                return;
            }

            if (me == 0)
            {
                if (string.Compare(arg[0], "SELF") == 0) System.Console.WriteLine("input.jump()->rewind(infile);"); 
                else
                {
                    if (sparta.infile != null) sparta.infile.Close();
                    sparta.infile = new FileStream(arg[0], FileMode.Open,FileAccess.Read);
                    if (sparta.infile == null)
                    {
                        //char str[128];
                        //sprintf(str, "Cannot open input script %s", arg[0]);
                        //error->one(FLERR, str);
                        string str = string.Format("Cannot open input script {0}", arg[0]);
                        sparta.error.one(str);
                    }
                    infiles[nfile - 1] = sparta.infile;
                }
            }

            if (narg == 2)
            {
                label_active = 1;
                if (labelstr != null) System.Console.WriteLine("input.jump()->delete[] labelstr;"); 
                int n = arg[1].Length + 1;
                labelstr = string.Copy(arg[1]);
            }
        }
        private void label(){}
        private void log(){}
        private void next_command(){}
        private void partition(){}
        private void print(){}
        private void quit(){}
        private void shell(){}
        private void variable_command(){}

        //// SPARTA commands

        private void boundary()
        {
            sparta.domain.set_boundary(narg, arg);
            Console.WriteLine("done");
        }
        private void bound_modify()
        {

        }
        private void collide_command()
        {
            string[] argss = new string[narg];
            Array.Copy(arg, 1, argss, 0, narg);

            if (narg < 1) sparta.error.all("Illegal collide command");

            if (sparta.suffix_enable!=0)
            {
                if (sparta.suffix!=null)
                {
                    string estyle = string.Format("{0}/{1}", argss[0], sparta.suffix);
                }
            }

            if (argss[0].Equals("none"))
            {
                sparta.collide = null;
            }

            switch (argss[0])
            {
                case "vss":
                    sparta.collide = new CollideVSS(sparta, narg, argss);
                    break;
                default:
                    break;
            }

            Console.WriteLine("done");


        }
        private void collide_modify()
        {

        }
        private void compute()
        {

        }
        private void dimension()
        {
            //System.Console.WriteLine("input.dimension");
            if (narg != 1) sparta.error.all( "Illegal dimension command");
            if (sparta.domain.box_exist!=0)
                sparta.error.all( "Dimension command after simulation box is defined");
            sparta.domain.dimension = System.Int32.Parse(arg[1]);
            if (sparta.domain.dimension != 2 && sparta.domain.dimension != 3)
                sparta.error.all( "Illegal dimension command");
            Console.WriteLine("done");
        }
        private void dump(){}
        private void dump_modify(){}
        private void fix()
        {
            string[] argss = new string[narg];
            Array.Copy(arg, 1, argss, 0, narg);
            sparta.modify.add_fix(narg, argss);
            Console.WriteLine("done");
        }
        private void global()
        {
            sparta.update.global(narg, arg);
            Console.WriteLine("done");
        }
        private void group(){}
        private void mixture()
        {
            string[] argss = new string[narg];
            Array.Copy(arg, 1, argss, 0, narg);
            sparta.particle.add_mixture(narg, argss);
            Console.WriteLine("done");
        }
        private void package(){}
        private void react_command(){}
        private void react_modify(){}
        private void region(){}
        private void reset_timestep(){}
        private void restart(){}
        private void seed()
        {
            if (narg != 1)
            {
                sparta.error.all( "Illegal seed command");
            }

            int seed =System.Int32.Parse(arg[1]);
            if (seed <= 0)
            {
                sparta.error.all( "Illegal seed command");
            }

            //update->ranmaster->init(seed);
            sparta.update.ranmaster.init(seed);

            Console.WriteLine("done");
            //System.Console.WriteLine("input.seed()");
        }
        private void species()
        {
            sparta.particle.add_species(narg, arg);
            Console.WriteLine("done");

        }
        private void stats()
        {
            sparta.output.set_stats(narg, arg);
        }
        private void stats_modify(){}
        private void stats_style(){}
        private void surf_collide()
        {
            sparta.surf.add_collide(narg, arg);
            Console.WriteLine("done");
        }
        private void surf_modify()
        {
            sparta.surf.modify_params(narg, arg);
            Console.WriteLine("done");
        }
        private void surf_react(){}
        private void timestep()
        {
            if (narg != 1) sparta.error.all("Illegal timestep command");
            double dt = double.Parse(arg[1]);
            if (dt <= 0.0) sparta.error.all("Illegal timestep command");
            sparta.update.dt = dt;
            Console.WriteLine("down");
        }
        private void uncompute(){}
        private void undump(){}
        private void unfix(){}
        private void units(){}
        private void weight(){}
    }
}