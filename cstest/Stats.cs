using System.Collections.Generic;
using System.IO;
using System.Text;
using bigint = System.Int64;
namespace cstest
{
    public delegate void ComputeHandler();
    public class Stats
    {
        // customize a new keyword by adding to this list:

        // step,elapsed,elaplong,dt,cpu,tpcpu,spcpu,wall,
        // np,ntouch,ncomm,nbound,nexit,nscoll,nscheck,ncoll,nattempt,nreact,nsreact,
        // npave,ntouchave,ncommave,nboundave,nexitave,nscollave,nscheckave,
        // ncollave,nattemptave,nreactave,nsreactave,
        // nparent,nchild,
        // vol,lx,ly,lz,xlo,xhi,ylo,yhi,zlo,zhi

        enum Enum1{ INT, FLOAT, BIGINT };
        enum Enum2{ SCALAR, VECTOR, ARRAY };

        public const int INVOKED_SCALAR = 1;
        public const int INVOKED_VECTOR = 2;
        public const int INVOKED_ARRAY = 4;

        public const int MAXLINE = 8192;               // make this 4x longer than Input::MAXLINE
        public const int DELTA = 8;

        private SPARTA sparta;

        public Stats(SPARTA sparta)
        {
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);

            sparta.mpi.MPI_Barrier(sparta.world);
            wall0 = sparta.mpi.MPI_Wtime();

            //line = new char[MAXLINE];

            keyword = null;
            //vfunc = null;
            vtype = null;

            field2index = null;
            argindex1 = null;
            argindex2 = null;

            // default args

            //string[] arg = new string[3];
            //arg[0] = (string)"step";
            //arg[1] = (string)"cpu";
            //arg[2] = (string)"np";
            string[] arg = new string[] { "step", "cpu", "np" };

            nfield = 3;
            Allocate();
            //set_fields(3, arg);

            //delete[] arg;

            // stats_modify defaults

            flushflag = 0;

            // format strings

            string bigint_format = "%lld";

            format_float_def = "%12.8g";
            format_int_def = "%8d";
            format_bigint_def=string.Format( "%%8%s", bigint_format);

            format_line_user = null;
            format_float_user = null;
            format_int_user = null;
            format_bigint_user = null;





        }

        public void init()
        {
            int i, n;

            // set format string for each field
            // add trailing '/n' to last value

            string format_line = null ;
            if (format_line_user!=null)
            {
                format_line=string.Copy(format_line_user);
            }
            //string ptr, format_line_ptr;
            //for (i = 0; i < nfield; i++)
            //{
            //    format[i][0] = '\0';

            //    if (format_line != null)
            //    {
            //        if (i == 0) format_line_ptr = strtok(format_line, " \0");
            //        else format_line_ptr = strtok(NULL, " \0");
            //    }

            //    if (format_column_user[i]!=null) ptr = format_column_user[i];
            //    else if (vtype[i] == FLOAT)
            //    {
            //        if (format_float_user) ptr = format_float_user;
            //        else if (format_line_user) ptr = format_line_ptr;
            //        else ptr = format_float_def;
            //    }
            //    else if (vtype[i] == INT)
            //    {
            //        if (format_int_user) ptr = format_int_user;
            //        else if (format_line_user) ptr = format_line_ptr;
            //        else ptr = format_int_def;
            //    }
            //    else if (vtype[i] == BIGINT)
            //    {
            //        if (format_bigint_user) ptr = format_bigint_user;
            //        else if (format_line_user) ptr = format_line_ptr;
            //        else ptr = format_bigint_def;
            //    }

            //    n = strlen(format[i]);
            //    sprintf(&format[i][n], "%s ", ptr);
            //}
            //strcat(format[nfield - 1], "\n");

            //delete[] format_line;
        }
        //public void modify_params(int, string[]);
        public void set_fields(int narg,string[] arg)
        {
            //deallocate();

            nfield = narg;
            Allocate();
            nfield = 0;
            foreach (string item in arg)
            {
                switch (item)
                {
                    case "step":
                        addfield("Step", compute_step, (int)Enum1.BIGINT);
                        break;
                    //case "elapsed":
                    //    addfield("Elapsed", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "elaplong":
                    //    addfield("elaplong", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "dt":
                    //    addfield("dt", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    case "cpu":
                        addfield("CPU", compute_cpu, (int)Enum1.FLOAT);
                        break;
                    //case "tpcpu":
                    //    addfield("tpcpu", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "spcpu":
                    //    addfield("spcpu", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "wall":
                    //    addfield("wall", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    case "np":
                        addfield("Np", compute_np, (int)Enum1.BIGINT);
                        break;
                    //case "ntouch":
                    //    addfield("ntouch", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "ncomm":
                    //    addfield("ncomm", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nbound":
                    //    addfield("nbound", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    case "nscoll":
                        addfield("Nscoll", compute_nscoll, (int)Enum1.BIGINT);
                        break;
                    case "nscheck":
                        addfield("Nscheck", compute_nscheck, (int)Enum1.BIGINT);
                        break;
                    case "ncoll":
                        addfield("Ncoll", compute_ncoll, (int)Enum1.BIGINT);
                        break;
                    case "nattempt":
                        addfield("Nattempt", compute_nattempt, (int)Enum1.BIGINT);
                        break;
                    //case "nreact":
                    //    addfield("nreact", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nsreact":
                    //    addfield("nsreact", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "npave":
                    //    addfield("npave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "ntouchave":
                    //    addfield("ntouchave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "ncommave":
                    //    addfield("ncommave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nboundave":
                    //    addfield("nboundave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nexitave":
                    //    addfield("nexitave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nscollave":
                    //    addfield("nscollave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nscheckave":
                    //    addfield("nscheckave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "ncollave":
                    //    addfield("ncollave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nattemptave":
                    //    addfield("nattemptave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nreactave":
                    //    addfield("nreactave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nsreactave":
                    //    addfield("nsreactave", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nparent":
                    //    addfield("nparent", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nchild":
                    //    addfield("nchild", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "nsplit":
                    //    addfield("nsplit", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "vol":
                    //    addfield("vol", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "lx":
                    //    addfield("lx", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "ly":
                    //    addfield("ly", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "lz":
                    //    addfield("lz", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "xlo":
                    //    addfield("xlo", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "xhi":
                    //    addfield("xhi", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "ylo":
                    //    addfield("ylo", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "yhi":
                    //    addfield("yhi", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "zlo":
                    //    addfield("zlo", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    //case "zhi":
                    //    addfield("zhi", compute_elapsed, (int)Enum1.BIGINT);
                    //    break;
                    default:
                        sparta.error.all( "Invalid keyword in stats_style command");
                        break;
                }
            }



        }
        public void header()
        {
            if (line==null)
            {
                line = new StringBuilder();
            }
            for (int i = 0; i < nfield; i++)
            {
                line .Append( string.Format(" {0,10} ",keyword[i]));
                
            }

            

            if (me == 0)
            {
                System.Console.WriteLine(line);
                if (sparta.screen!=null) new StreamWriter(sparta.screen).WriteLine(line);
                if (sparta.logfile!=null) new StreamWriter(sparta.logfile).WriteLine(line);
            }
        }
        public void compute(int flag)
        {
            int i;

            firststep = flag;
            line.Clear();
            // invoke Compute methods needed for stats keywords

            for (i = 0; i < ncompute; i++)
                if (compute_which[i] == (int)Enum2.SCALAR)
                {
                    if ((computes[i].invoked_flag & INVOKED_SCALAR)==0)
                    {
                        computes[i].compute_scalar();
                        computes[i].invoked_flag |= INVOKED_SCALAR;
                    }
                }
                else if (compute_which[i] == (int)Enum2.VECTOR)
                {
                    if ((computes[i].invoked_flag & INVOKED_VECTOR)==0)
                    {
                        computes[i].compute_vector();
                        computes[i].invoked_flag |= INVOKED_VECTOR;
                    }
                }
                else if (compute_which[i] == (int)Enum2.ARRAY)
                {
                    if ((computes[i].invoked_flag & INVOKED_ARRAY)==0)
                    {
                        computes[i].compute_array();
                        computes[i].invoked_flag |= INVOKED_ARRAY;
                    }
                }

            // add each stat value to line with its specific format

            for (ifield = 0; ifield < nfield; ifield++)
            {
                vfunc[ifield]();
                format[ifield] = " {0,10:G6} ";
                if (vtype[ifield] == (int)Enum1.FLOAT)
                {
                    string str = string.Format(format[ifield], dvalue);
                    line.Append(str);
                }
                else if (vtype[ifield] == (int)Enum1.INT)
                {
                    string str = string.Format(format[ifield], ivalue);
                    line.Append(str);
                }
                else if (vtype[ifield] == (int)Enum1.BIGINT)
                {
                    string str = string.Format(format[ifield], bivalue);
                    line.Append(str);
                }
            }

            // print line to screen and logfile

            if (me == 0)
            {
                System.Console.WriteLine(line);
                if (sparta.screen!=null) new StreamWriter(sparta.screen).WriteLine(line);
                if (sparta.logfile!=null)
                {
                    new StreamWriter(sparta.screen).WriteLine(line);
                    if (flushflag!=0) sparta.logfile.Flush();
                }
            }
        }
        //public int evaluate_keyword(string, double*);

        private StringBuilder line;
        private string[] keyword;
        private int[] vtype;

        private int nfield;
        private int me;

        private string[] format;
        private string format_line_user;
        private string format_float_user;
        private string[] format_int_user,format_bigint_user;
        private string[] format_column_user;

        private string format_float_def;
        private string format_int_def;
        private string format_bigint_def;

        private int firststep;
        private int flushflag, lineflag;

        private double wall0;
        private double last_tpcpu, last_spcpu;
        private double last_time;
        private bigint last_step;
        // data used by routines that compute single values
        private int ivalue;            // integer value to print
        private double dvalue;         // double value to print
        private bigint bivalue;        // big integer value to print
        private int ifield;            // which field in thermo output is being computed
        private int[] field2index;      // which compute,fix,variable calcs this field
        private int[] argindex1;        // indices into compute,fix scalar,vector
        private int[] argindex2;

        int ncompute;                // # of Compute objects called by stats
        string[] id_compute;           // their IDs
        int[] compute_which;          // 0/1/2 if should call scalar,vector,array
        List<Compute> computes;    // list of ptrs to the Compute objects

        int nfix;                    // # of Fix objects called by stats
        string[] id_fix;               // their IDs
        List<Fix> fixes;           // list of ptrs to the Fix objects

        int nsurfcollide;            // # of SurfCollide objs called by stats
        string[] id_surf_collide;      // their IDs
        List<SurfCollide> sc;      // list of ptrs to SurfCollide objects

        int nsurfreact;             // # of SurfReact objects called by stats
        string[] id_surf_react;       // their IDs
        List<SurfReact> sr;       // list of ptrs to SurfReact objects

        int nvariable;               // # of variables evaulated by stats
        string[] id_variable;          // list of variable names
        int[] variables;              // list of Variable indices

        



        // private methods

        void Allocate()
        {
            int n = nfield;

            keyword =new string[n];

            vfunc = new ComputeHandler[n];
            vtype = new int[n];

            format = new string[n];
            format_column_user = new string[n];
            for (int i = 0; i < n; i++) format_column_user[i] = null;

            field2index = new int[n];
            argindex1 = new int[n];
            argindex2 = new int[n];

            // memory for computes, fixes, variables

            ncompute = 0;
            id_compute = new string[n];
            compute_which = new int[n];
            computes = new List<Compute>(n);

            nfix = 0;
            id_fix = new string[n];
            fixes = new List<Fix>(n);

            nsurfcollide = 0;
            id_surf_collide = new string[n];
            sc = new List<SurfCollide>(n);

            nsurfreact = 0;
            id_surf_react = new string[n];
            sr = new List<SurfReact>(n);

            nvariable = 0;
            id_variable = new string[n];
            variables = new int[n];
        }
        //void deallocate();

        //      int add_compute(const string, int);
        //int add_fix(const string);
        //      int add_surf_collide(const string);
        //      int add_surf_react(const string);
        //      int add_variable(const string);

        //      typedef void (Stats::* FnPtr) ();
        void addfield(string key, ComputeHandler func, int typeflag)
        {
            keyword[nfield]=string.Copy(key);
            vfunc[nfield] = func;
            vtype[nfield] = typeflag;
            nfield++;

        }
        ComputeHandler[] vfunc;                // list of ptrs to functions

        //void compute_compute();        // functions that compute a single value
        //void compute_fix();            // via calls to Compute,Fix,
        //void compute_surf_collide();   //   SurfCollide,SurfReact,Variable classes
        //void compute_surf_react();
        //void compute_variable();

        // functions that compute a single value
        // customize a new keyword by adding a method prototype

        void compute_step()
        {
            bivalue = sparta.update.ntimestep;
        }
        //void compute_elapsed()
        //{
        //    System.Console.WriteLine("compute elapsed" );
        //}
        //void compute_elaplong();
        //void compute_dt();
        void compute_cpu()
        {
            if (firststep == 0) dvalue = 0.0;
            else dvalue = sparta.timer.elapsed((int)Timer.Enum1.TIME_LOOP);
        }
        //void compute_tpcpu();
        //void compute_spcpu();
        //void compute_wall();

        void compute_np()
        {
            bigint n = sparta.particle.nlocal;
            sparta.mpi.MPI_Allreduce(ref n,ref sparta.particle.nglobal, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            bivalue = sparta.particle.nglobal;
        }
        //void compute_ntouch();
        //void compute_ncomm();
        //void compute_nbound();
        //void compute_nexit();
        void compute_nscoll()
        {
            bigint n = sparta.update.nscollide_one;
            sparta.mpi.MPI_Allreduce(ref n, ref bivalue, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
        }
        void compute_nscheck()
        {
            bigint n = sparta.update.nscheck_one;
            sparta.mpi.MPI_Allreduce(ref n, ref bivalue, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
        }
        void compute_ncoll()
        {
            if (sparta.collide==null) bivalue = 0;
            else
            {
                bigint n = sparta.collide.ncollide_one;
                sparta.mpi.MPI_Allreduce(ref n, ref bivalue, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            }
        }
        void compute_nattempt()
        {
            if (sparta.collide == null) bivalue = 0;
            else
            {
                bigint n = sparta.collide.nattempt_one;
                sparta.mpi.MPI_Allreduce(ref n, ref bivalue, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            }
        }
        //void compute_nreact();
        //void compute_nsreact();

        //void compute_npave();
        //void compute_ntouchave();
        //void compute_ncommave();
        //void compute_nboundave();
        //void compute_nexitave();
        //void compute_nscollave();
        //void compute_nscheckave();
        //void compute_ncollave();
        //void compute_nattemptave();
        //void compute_nreactave();
        //void compute_nsreactave();

        //void compute_nparent();
        //void compute_nchild();
        //void compute_nsplit();

        //void compute_vol();
        //void compute_lx();
        //void compute_ly();
        //void compute_lz();

        //void compute_xlo();
        //void compute_xhi();
        //void compute_ylo();
        //void compute_yhi();
        //void compute_zlo();
        //void compute_zhi();
    }
}