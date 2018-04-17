using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using MPI_Comm = System.Int32;
using smallint = System.Int32;
using bigint = System.Int64;
//todo 太多重复代码 需要优化 

namespace cstest
{

    public class SPARTA
    {
        public MPI mpi;
        public Memory memory;          // memory allocation functions
        public Error error;            // error handling
        public Universe universe;      // universe of processors
        public Input input;            // input script processing
                                       // ptrs to top-level SPARTA-specific classes
        public Particle particle;      // particles
        public Update update;          // timestepper
        public Comm comm;              // inter-processor communication
        public Domain domain;          // simulation box
        public Modify modify;          // fixes and computes
        public Grid grid;              // volumetric grid cells
        public Surf surf;              // surface elements
        public Collide collide;        // collisions and chemistry
        public React react;            // chemistry reactions
        public Output output;          // stats/dump/restart
        public Timer timer;            // CPU timing info

        public MPI_Comm world;                // MPI communicator
        public FileStream infile;                  // infile
        public FileStream screen;                  // screen output
        public FileStream logfile;                 // logfile

        public string suffix;                  // suffix to add to input script style names
        public int suffix_enable;             // 1 if suffixes are enabled, 0 if disabled
        string[][] packargs;              // arguments for cmdline package commands
        public int num_package;               // number of cmdline package commands

        public KokkosSPARTA kokkos;    // KOKKOS accelerator class
        public MemoryKokkos memoryKK;  // KOKKOS version of Memory class

        // other top-level SPARTA classes and variables

        public SPARTA(string[] args,MPI mpi, MPI_Comm communicator)
        {
            int narg = args.Length;
            string[] arg = args;
            this.mpi = mpi;
            error = new Error(this);
            universe = new Universe(this, communicator);
            output = null;

            screen = null;
            logfile = null;

            // parse input switches

            int inflag = 0;
            int screenflag = 0;
            int logflag = 0;
            int partscreenflag = 0;
            int partlogflag = 0;
            int kokkosflag = 0;
            int helpflag = 0;

            suffix = null;
            suffix_enable = 0;
            packargs = null;
            num_package = 0;
            int kkfirst, kklast;

            int npack = 0;
            int[] pfirst = null;
            int[] plast = null;

            //int iarg = 1;
            int iarg = 0;



            while (iarg < narg)
            {
                if (string.Compare(arg[iarg], "-partition") == 0 ||
                string.Compare(arg[iarg], "-p") == 0)
                {
                    universe.existflag = 1;
                    if (iarg + 2 > narg)
                        error.universe_all("Invalid command-line argument");
                    iarg++;
                    while (iarg < narg && arg[iarg][0] != '-')
                    {
                        universe.Add_world(arg[iarg]);
                        iarg++;
                    }
                }
                else if (string.Equals(arg[iarg], "-in") ||
                     string.Compare(arg[iarg], "-i") == 0)
                {
                    if (iarg + 2 > narg)
                        error.universe_all("Invalid command-line argument");
                    inflag = iarg + 1;
                    iarg += 2;
                }
                else if (string.Compare(arg[iarg], "-screen") == 0 ||
                     string.Compare(arg[iarg], "-sc") == 0)
                {
                    if (iarg + 2 > narg)
                        error.universe_all("Invalid command-line argument");
                    screenflag = iarg + 1;
                    iarg += 2;
                }
                else if (string.Compare(arg[iarg], "-log") == 0 ||
                     string.Compare(arg[iarg], "-l") == 0)
                {
                    if (iarg + 2 > narg)
                        error.universe_all("Invalid command-line argument");
                    logflag = iarg + 1;
                    iarg += 2;
                }
                else if (string.Compare(arg[iarg], "-var") == 0 ||
                     string.Compare(arg[iarg], "-v") == 0)
                {
                    if (iarg + 3 > narg)
                        error.universe_all("Invalid command-line argument");
                    iarg += 3;
                    while (iarg < narg && arg[iarg][0] != '-') iarg++;
                }
                else if (string.Compare(arg[iarg], "-echo") == 0 ||
                     string.Compare(arg[iarg], "-e") == 0)
                {
                    if (iarg + 2 > narg)
                        error.universe_all("Invalid command-line argument");
                    iarg += 2;
                }
                else if (string.Compare(arg[iarg], "-pscreen") == 0 ||
                     string.Compare(arg[iarg], "-ps") == 0)
                {
                    if (iarg + 2 > narg)
                        error.universe_all("Invalid command-line argument");
                    partscreenflag = iarg + 1;
                    iarg += 2;
                }
                else if (string.Compare(arg[iarg], "-plog") == 0 ||
                     string.Compare(arg[iarg], "-pl") == 0)
                {
                    if (iarg + 2 > narg)
                        error.universe_all("Invalid command-line argument");
                    partlogflag = iarg + 1;
                    iarg += 2;
                }
                else if (string.Compare(arg[iarg], "-kokkos") == 0 ||
                         string.Compare(arg[iarg], "-k") == 0)
                {
                    if (iarg + 2 > narg)
                        error.universe_all("Invalid command-line argument");
                    if (string.Compare(arg[iarg + 1], "on") == 0) kokkosflag = 1;
                    else if (string.Compare(arg[iarg + 1], "off") == 0) kokkosflag = 0;
                    else error.universe_all("Invalid command-line argument");
                    iarg += 2;
                    // delimit any extra args for the Kokkos instantiation
                    kkfirst = iarg;
                    while (iarg < narg && arg[iarg][0] != '-') iarg++;
                    kklast = iarg;
                }
                else if (string.Compare(arg[iarg], "-package") == 0 ||
                         string.Compare(arg[iarg], "-pk") == 0)
                {
                    if (iarg + 2 > narg)
                        error.universe_all("Invalid command-line argument");
                    //memory.grow(pfirst, npack + 1, "sparta:pfirst");
                    //memory.grow(plast, npack + 1, "sparta:plast");

                    // delimit args for package command invocation
                    // any package arg with leading "-" will be followed by numeric digit
                    iarg++;
                    pfirst[npack] = iarg;
                    while (iarg < narg)
                    {
                        if (arg[iarg][0] != '-') iarg++;
                        else if (char.IsDigit(arg[iarg][1])) iarg++;
                        else break;
                    }
                    plast[npack++] = iarg;
                }
                else if (string.Compare(arg[iarg], "-suffix") == 0 ||
                         string.Compare(arg[iarg], "-sf") == 0)
                {
                    if (iarg + 2 > narg)
                        error.universe_all("Invalid command-line argument");
                    //delete[] suffix;
                    int n = arg[iarg + 1].Length + 1;

                    suffix = string.Copy(arg[iarg + 1]);
                    suffix_enable = 1;
                    iarg += 2;

                }
                else if (string.Compare(arg[iarg], "-help") == 0 ||
                     string.Compare(arg[iarg], "-h") == 0)
                {
                    if (iarg + 1 > narg)
                        error.universe_all("Invalid command-line argument");
                    helpflag = 1;
                    iarg += 1;
                }
                else error.universe_all("Invalid command-line argument");
            }
            // if no partition command-line switch, universe is one world with all procs

            if (universe.existflag == 0) universe.Add_world(null);

            // sum of procs in all worlds must equal total # of procs

            if (universe.Consistent() == 0)
                error.universe_all("Processor partitions are inconsistent");

            // universe cannot use stdin for input file

            if (universe.existflag != 0 && inflag == 0)
                error.universe_all("Must use -in switch with multiple partitions");

            // if no partition command-line switch, cannot use -pscreen option

            if (universe.existflag == 0 && partscreenflag != 0)
                error.universe_all("Can only use -pscreen with multiple partitions");

            // if no partition command-line switch, cannot use -plog option

            if (universe.existflag == 0 && partlogflag != 0)
                error.universe_all("Can only use -plog with multiple partitions");

            // set universe screen and logfile

            if (universe.me == 0)
            {
                if (screenflag == 0)
                {
                    Console.WriteLine("universe.uscreen=stdout");
                    universe.uscreen = new FileStream("screen", FileMode.Create, FileAccess.Write);
                }
                //universe.uscreen = stdout;
                else if (string.Compare(arg[screenflag], "none") == 0)
                    universe.uscreen = null;
                else
                {
                    universe.uscreen = new FileStream(arg[screenflag], FileMode.Create, FileAccess.Write);
                    //universe.uscreen = fopen(arg[screenflag], "w");
                    if (universe.uscreen == null)
                        error.universe_one("Cannot open universe screen file");
                }
                if (logflag == 0)
                {
                    universe.ulogfile = new FileStream("log.sparta", FileMode.Create, FileAccess.Write);
                    //universe.ulogfile = fopen("log.sparta", "w");
                    if (universe.ulogfile == null)
                        error.universe_one("Cannot open log.sparta");
                }
                else if (string.Compare(arg[logflag], "none") == 0)
                    universe.ulogfile = null;
                else
                {
                    universe.ulogfile = new FileStream(arg[logflag], FileMode.Create, FileAccess.Write);
                    //universe.ulogfile = fopen(arg[logflag], "w");
                    if (universe.ulogfile == null)
                        error.universe_one("Cannot open universe log file");
                }
            }
            if (universe.me > 0)
            {
                if (screenflag == 0) Console.WriteLine("universe.uscreen = stdout;");
                else universe.uscreen = null;
                universe.ulogfile = null;
            }

            // make universe and single world the same, since no partition switch
            // world inherits settings from universe
            // set world screen, logfile, communicator, infile
            // open input script if from file

            if (universe.existflag == 0)
            {
                screen = universe.uscreen;
                logfile = universe.ulogfile;
                world = universe.uworld;
                infile = null;

                if (universe.me == 0)
                {
                    if (inflag == 0) Console.WriteLine("infile = stdin;");
                    else
                    {
                        infile = new FileStream(arg[inflag], FileMode.Open, FileAccess.Read);
                        //infile = fopen(arg[inflag], "r");
                    }

                    if (infile == null)
                    {
                        string str = string.Format("Cannot open input script {0}", arg[inflag]);
                        //char str[128];
                        //sprintf(str, "Cannot open input script %s", arg[inflag]);
                        error.one(str);
                    }
                }

                if (universe.me == 0)
                {
                    string strstr = string.Format("SPARTA ({0})\n", universe.version);
                    if (screen != null)
                    {
                        new StreamWriter(screen).Write(strstr);
                        Console.WriteLine(strstr);
                    }

                    if (logfile != null)
                    {
                        new StreamWriter(logfile).Write(strstr);
                    }
                }

                // universe is one or more worlds, as setup by partition switch
                // split universe communicator into separate world communicators
                // set world screen, logfile, communicator, infile
                // open input script

            }
            else
            {
                int me = 0;
                mpi.MPI_Comm_split(universe.uworld, universe.iworld, 0, ref world);
                mpi.MPI_Comm_rank(world, ref me);

                if (me == 0)
                    if (partscreenflag == 0)
                        if (screenflag == 0)
                        {
                            string str = string.Format("screen.{0}", universe.iworld);
                            //char str[32];
                            //sprintf(str, "screen.%d", universe.iworld);
                            screen = new FileStream(str, FileMode.Create, FileAccess.Write);
                            //screen = fopen(str, "w");
                            if (screen == null) error.one("Cannot open screen file");
                        }
                        else if (string.Compare(arg[screenflag], "none") == 0)
                            screen = null;
                        else
                        {
                            string str = string.Format("{0}.{1}", arg[screenflag], universe.iworld);
                            //char str[128];
                            //sprintf(str, "%s.%d", arg[screenflag], universe.iworld);
                            screen = new FileStream(str, FileMode.Create, FileAccess.Write);
                            //screen = fopen(str, "w");
                            if (screen == null) error.one("Cannot open screen file");
                        }
                    else if (string.Compare(arg[partscreenflag], "none") == 0)
                        screen = null;
                    else
                    {
                        string str = string.Format("{0}.{1}", arg[partscreenflag], universe.iworld);
                        //char str[128];
                        //sprintf(str, "%s.%d", arg[partscreenflag], universe.iworld);
                        screen = new FileStream(str, FileMode.Create, FileAccess.Write);
                        //screen = fopen(str, "w");
                        if (screen == null) error.one("Cannot open screen file");
                    }
                else screen = null;

                if (me == 0)
                    if (partlogflag == 0)
                        if (logflag == 0)
                        {
                            string str = string.Format("log.sparta.{0}", universe.iworld);
                            logfile = new FileStream(str, FileMode.Create, FileAccess.Write);

                            //char str[32];
                            //sprintf(str, "log.sparta.%d", universe.iworld);
                            //logfile = fopen(str, "w");
                            if (logfile == null) error.one("Cannot open logfile");
                        }
                        else if (string.Compare(arg[logflag], "none") == 0)
                            logfile = null;
                        else
                        {
                            string str = string.Format("{0}.{1}", arg[logflag], universe.iworld);
                            logfile = new FileStream(str, FileMode.Create, FileAccess.Write);
                            //char str[128];
                            //sprintf(str, "%s.%d", arg[logflag], universe.iworld);
                            //logfile = fopen(str, "w");
                            if (logfile == null) error.one("Cannot open logfile");
                        }
                    else if (string.Compare(arg[partlogflag], "none") == 0)
                        logfile = null;
                    else
                    {
                        string str = string.Format("{0}.{1}", arg[partlogflag], universe.iworld);
                        logfile = new FileStream(str, FileMode.Create, FileAccess.Write);
                        //char str[128];
                        //sprintf(str, "%s.%d", arg[partlogflag], universe.iworld);
                        //logfile = fopen(str, "w");
                        if (logfile == null) error.one("Cannot open logfile");
                    }
                else logfile = null;

                if (me == 0)
                {
                    infile = new FileStream(arg[inflag], FileMode.Open, FileAccess.Read);
                    //infile = fopen(arg[inflag], "r");
                    if (infile == null)
                    {
                        string str = string.Format("Cannot open input script {0}", arg[inflag]);
                        //char str[128];
                        //sprintf(str, "Cannot open input script %s", arg[inflag]);
                        error.one(str);
                    }
                }
                else infile = null;

                // screen and logfile messages for universe and world

                if (universe.me == 0)
                {
                    string str1 = string.Format("SPARTA ({0})\n", universe.version);
                    string str2 = string.Format("Running on {0} partitions of processors\n", universe.nworlds);
                    if (universe.uscreen != null)
                    {
                        new StreamWriter(universe.uscreen).Write(str1);
                        new StreamWriter(universe.uscreen).Write(str2);
                        Console.WriteLine(str1);
                        Console.WriteLine(str2);
                        //fprintf(universe.uscreen, "SPARTA (%s)\n", universe.version);
                        //fprintf(universe.uscreen, "Running on %d partitions of processors\n",
                        //    universe.nworlds);
                    }
                    if (universe.ulogfile != null)
                    {
                        new StreamWriter(universe.ulogfile).Write(str1);
                        new StreamWriter(universe.ulogfile).Write(str2);
                        //fprintf(universe.ulogfile, "SPARTA (%s)\n", universe.version);
                        //fprintf(universe.ulogfile, "Running on %d partitions of processors\n",
                        //    universe.nworlds);
                    }
                }

                if (me == 0)
                {

                    string str1 = string.Format("SPARTA ({0})\n", universe.version);
                    string str2 = string.Format("Processor partition = {0}\n", universe.iworld);
                    if (screen != null)
                    {
                        new StreamWriter(screen).Write(str1);
                        new StreamWriter(screen).Write(str2);
                        Console.WriteLine(str1);
                        Console.WriteLine(str2);
                        //fprintf(screen, "SPARTA (%s)\n", universe.version);
                        //fprintf(screen, "Processor partition = %d\n", universe.iworld);
                    }
                    if (logfile != null)
                    {
                        new StreamWriter(logfile).Write(str1);
                        new StreamWriter(logfile).Write(str2);
                        //fprintf(logfile, "SPARTA (%s)\n", universe.version);
                        //fprintf(logfile, "Processor partition = %d\n", universe.iworld);
                    }
                }
            }

            // check datatype settings in spatype.h

            //if (sizeof(smallint) != sizeof(int))
            //    error.all(  "Smallint setting in spatype.h is invalid");
            //if (sizeof(bigint) < sizeof(smallint))
            //    error.all(  "Bigint setting in spatype.h is invalid");

            int mpisize=1;
            mpi.MPI_Type_size(MPI.MPI_LONG_LONG,ref mpisize);
            if (mpisize != sizeof(bigint))
                error.all(  "MPI_SPARTA_BIGINT and bigint in spatype.h are not compatible");

            //if (sizeof(smallint) != 4 || sizeof(bigint) != 8)
            //    error.all(  "Small,big integers are not sized correctly");

            // error check on accelerator packages

            // create Kokkos class if KOKKOS installed, unless explicitly switched off
            // instantiation creates dummy Kokkos class if KOKKOS is not installed
            // add args between kkfirst and kklast to Kokkos instantiation

            kokkos = null;
            if (kokkosflag == 1)
            {
                Console.WriteLine("new kokkosSPARTA");
                //kokkos = new KokkosSPARTA(this, kklast - kkfirst, &arg[kkfirst]);
                //if (!kokkos->kokkos_exists)
                //    error->all(FLERR, "Cannot use -kokkos on without KOKKOS installed");
            }

            // allocate input class now that MPI is fully setup

            input = new Input(this, arg);

            // copy package cmdline arguments

            if (npack > 0)
            {
                num_package = npack;
                Console.WriteLine("package command,not test yet");
                //packargs = new char**[npack];
                for (int i = 0; i < npack; ++i)
                {
                    int n = plast[i] - pfirst[i];
                    //packargs[i] = new char*[n + 1];
                    for (int j = 0; j < n; ++j)
                       // packargs[i][j] = strdup(arg[pfirst[i] + j]);
                    packargs[i][n] = null;
                }
                //memory->destroy(pfirst);
                //memory->destroy(plast);
            }
            // allocate fundamental classes

            create();
            post_create();

            // if helpflag set, print help and exit

            if (helpflag!=0)
            {
                if (universe.me == 0) print_styles();
                error.done();
            }
        }

   //allocate single instance of top-level classes
   //fundamental classes are allocated in constructor
   //some classes have package variants
   
        public void create()
        {
            //if (kokkos) update = new UpdateKokkos(this);
            //else update = new Update(this);

            //if (kokkos) particle = new ParticleKokkos(this);
            //else particle = new Particle(this);

            //if (kokkos) comm = new CommKokkos(this);
            //else comm = new Comm(this);

            //if (kokkos) domain = new DomainKokkos(this);
            //else domain = new Domain(this);

            //if (kokkos) grid = new GridKokkos(this);
            //else grid = new Grid(this);

            //if (kokkos) surf = new SurfKokkos(this);
            //else surf = new Surf(this);

            update = new Update(this);
            particle = new Particle(this);
            comm = new Comm(this);
            domain = new Domain(this);
            grid = new Grid(this);
            surf = new Surf(this);



            collide = null;
            react = null;

            //if (kokkos) modify = new ModifyKokkos(this);
            //else modify = new Modify(this);

            modify = new Modify(this);

            output = new Output(this);
            timer = new Timer(this);
        }
   //     check suffix consistency with installed packages
   //invoke package-specific deafult package commands
   //  only invoke if suffix is set and enabled
   //called from SPARTA constructor and after clear() command
   //  so that package-specific core classes have been instantiated
        public void post_create()
        {
            // default package commands triggered by "-k on"

            //if (kokkos && kokkos->kokkos_exists) input->one("package kokkos");

            // suffix will always be set if suffix_enable = 1
            // check that KOKKOS package classes was instantiated

            if (suffix_enable==0) return;

            //if (string.Equals(suffix, "kk") &&
            //    (kokkos == null || kokkos.kokkos_exists == 0))
            //    error.all( "Using suffix kk without KOKKOS package enabled");

            // invoke any command-line package commands

            if (num_package!=0)
            {
                string str;
                for (int i = 0; i < num_package; i++)
                {
                    str= "package";
                    //for (char** ptr = packargs[i]; *ptr != NULL; ++ptr)
                    //{
                    //    if (strlen(str) + strlen(*ptr) + 2 > 256)
                    //        error->all(FLERR, "Too many -pk arguments in command line");
                    //    strcat(str, " ");
                    //    strcat(str, *ptr);
                    //}
                    foreach (string item in packargs[i])
                    {
                        string.Concat(str, " ", item);
                    }

                    input.one(str);
                }
            }
        }
        public void init()
        {
            update.init();
            particle.init();
            //todo: comm init is empty?
            comm.init();
            domain.init();
            grid.init();
            surf.init();
            if (react != null)
            {
                //react.init();
                Console.WriteLine("react init");
            }
            if (collide != null)
            {
                collide.init();
            }
            modify.init();
            //output.init();
            //timer.init();
        }
        public void destroy()
        {

        }

        void print_styles()
        {
            Console.WriteLine("\nList of style options included in this executable:\n\n");
            Console.WriteLine("Collide styles:");
        }
        public static void fprintf(FileStream file,string str,params object[] param)
        {
            string tempstr=string.Format(str, param);
            new StreamWriter(file).Write(tempstr);
        }
    }
}
