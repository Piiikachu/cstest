using bigint = System.Int64;
namespace cstest
{
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

            line = new char[MAXLINE];

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

        //public void Init()
        //{
        //    int i, n;

        //    // set format string for each field
        //    // add trailing '/n' to last value

        //    string format_line = null;
        //    if (format_line_user!=null)
        //    {
        //        format_line =string.Copy(format_line_user);
        //    }

        //    string ptr,format_line_ptr;
        //    for (i = 0; i < nfield; i++)
        //    {
        //        format[i][0] = '\0';

        //        if (format_line!=null)
        //        {
        //            if (i == 0) format_line_ptr = strtok(format_line, " \0");
        //            else format_line_ptr = strtok(NULL, " \0");



        //        }

        //        if (format_column_user[i]!=null)
        //        {
        //            ptr = format_column_user[i];
        //        }
        //        else if (vtype[i] == Enum1.FLOAT)
        //        {
        //            if (format_float_user) ptr = format_float_user;
        //            else if (format_line_user) ptr = format_line_ptr;
        //            else ptr = format_float_def;
        //        }
        //        else if (vtype[i] == INT)
        //        {
        //            if (format_int_user) ptr = format_int_user;
        //            else if (format_line_user) ptr = format_line_ptr;
        //            else ptr = format_int_def;
        //        }
        //        else if (vtype[i] == BIGINT)
        //        {
        //            if (format_bigint_user) ptr = format_bigint_user;
        //            else if (format_line_user) ptr = format_line_ptr;
        //            else ptr = format_bigint_def;
        //        }

        //        n = strlen(format[i]);
        //        sprintf(&format[i][n], "%s ", ptr);
        //    }
        //    strcat(format[nfield - 1], "\n");

        //    delete[] format_line;

        //    // find current ptr for each SurfCollide and SurfReact ID

        //    int m;

        //    for (int i = 0; i < nsurfcollide; i++)
        //    {
        //        m = surf->find_collide(id_surf_collide[i]);
        //        if (m < 0) error->all(FLERR, "Could not find stats surf collide ID");
        //        sc[i] = surf->sc[m];
        //    }
        //    for (int i = 0; i < nsurfreact; i++)
        //    {
        //        m = surf->find_react(id_surf_react[i]);
        //        if (m < 0) error->all(FLERR, "Could not find stats surf collide ID");
        //        sr[i] = surf->sr[m];
        //    }

        //    // find current ptr for each Compute ID

        //    for (int i = 0; i < ncompute; i++)
        //    {
        //        m = modify->find_compute(id_compute[i]);
        //        if (m < 0) error->all(FLERR, "Could not find stats compute ID");
        //        computes[i] = modify->compute[m];
        //    }

        //    // find current ptr for each Fix ID
        //    // check that fix frequency is acceptable with stats output frequency

        //    for (int i = 0; i < nfix; i++)
        //    {
        //        m = modify->find_fix(id_fix[i]);
        //        if (m < 0) error->all(FLERR, "Could not find stats fix ID");
        //        fixes[i] = modify->fix[m];
        //        if (output->stats_every % fixes[i]->global_freq)
        //            error->all(FLERR, "Stats and fix not computed at compatible times");
        //    }

        //    // find current ptr for each Variable ID

        //    for (int i = 0; i < nvariable; i++)
        //    {
        //        m = input->variable->find(id_variable[i]);
        //        if (m < 0) error->all(FLERR, "Could not find stats variable name");
        //        variables[i] = m;
        //    }
        //}
        //public void modify_params(int, string[]);
        public void set_fields(int narg,string[] arg)
        {

        }
        //public void header();
        //public void compute(int);
        //public int evaluate_keyword(string, double*);

        private char[] line;
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
        Compute[] computes;    // list of ptrs to the Compute objects

        int nfix;                    // # of Fix objects called by stats
        string[] id_fix;               // their IDs
        Fix[] fixes;           // list of ptrs to the Fix objects

        int nsurfcollide;            // # of SurfCollide objs called by stats
        string[] id_surf_collide;      // their IDs
        SurfCollide[] sc;      // list of ptrs to SurfCollide objects

        int nsurfreact;             // # of SurfReact objects called by stats
        string[] id_surf_react;       // their IDs
        SurfReact[] sr;       // list of ptrs to SurfReact objects

        int nvariable;               // # of variables evaulated by stats
        string[] id_variable;          // list of variable names
        int[] variables;              // list of Variable indices





        // private methods

        void Allocate()
        {
            int n = nfield;

            keyword =new string[n];
            //for (int i = 0; i < n; i++) keyword[i] = new string;
            //vfunc = new FnPtr[n];
            vtype = new int[n];

            format = new string[n];
            //for (int i = 0; i < n; i++) format[i] = new char[32];
            format_column_user = new string[n];
            for (int i = 0; i < n; i++) format_column_user[i] = null;

            field2index = new int[n];
            argindex1 = new int[n];
            argindex2 = new int[n];

            // memory for computes, fixes, variables

            ncompute = 0;
            id_compute = new string[n];
            compute_which = new int[n];
            computes = new Compute[n];

            nfix = 0;
            id_fix = new string[n];
            fixes = new Fix[n];

            nsurfcollide = 0;
            id_surf_collide = new string[n];
            sc = new SurfCollide[n];

            nsurfreact = 0;
            id_surf_react = new string[n];
            sr = new SurfReact[n];

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
  //void addfield(const string, FnPtr, int);
  //FnPtr* vfunc;                // list of ptrs to functions

  //      void compute_compute();        // functions that compute a single value
  //      void compute_fix();            // via calls to Compute,Fix,
  //      void compute_surf_collide();   //   SurfCollide,SurfReact,Variable classes
  //      void compute_surf_react();
  //      void compute_variable();

  //      // functions that compute a single value
  //      // customize a new keyword by adding a method prototype

  //      void compute_step();
  //      void compute_elapsed();
  //      void compute_elaplong();
  //      void compute_dt();
  //      void compute_cpu();
  //      void compute_tpcpu();
  //      void compute_spcpu();
  //      void compute_wall();

  //      void compute_np();
  //      void compute_ntouch();
  //      void compute_ncomm();
  //      void compute_nbound();
  //      void compute_nexit();
  //      void compute_nscoll();
  //      void compute_nscheck();
  //      void compute_ncoll();
  //      void compute_nattempt();
  //      void compute_nreact();
  //      void compute_nsreact();

  //      void compute_npave();
  //      void compute_ntouchave();
  //      void compute_ncommave();
  //      void compute_nboundave();
  //      void compute_nexitave();
  //      void compute_nscollave();
  //      void compute_nscheckave();
  //      void compute_ncollave();
  //      void compute_nattemptave();
  //      void compute_nreactave();
  //      void compute_nsreactave();

  //      void compute_nparent();
  //      void compute_nchild();
  //      void compute_nsplit();

  //      void compute_vol();
  //      void compute_lx();
  //      void compute_ly();
  //      void compute_lz();

  //      void compute_xlo();
  //      void compute_xhi();
  //      void compute_ylo();
  //      void compute_yhi();
  //      void compute_zlo();
  //      void compute_zhi();
    }
}