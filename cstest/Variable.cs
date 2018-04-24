using System;
using System.IO;

namespace cstest
{

    public class Variable
    {
        public const int VARDELTA = 4;
        public const int MAXLEVEL = 4;
        public const int MAXLINE = 256;
        public const int CHUNK = 1024;
        public const int VALUELENGTH = 64;

        //public const int MYROUND(a) (( a-floor(a) ) >= .5) ? ceil(a) : floor(a)

        enum Enum1 : int
        {
            INDEX, LOOP, WORLD, UNIVERSE, ULOOP, STRING, GETENV,
            SCALARFILE, FORMAT, EQUAL, PARTICLE, GRID, SURF, INTERNAL
        };
        enum Enum2 { ARG, OP };

        // customize by adding a function
        // if add before OR,
        // also set precedence level in constructor and precedence length in *.h

        enum Precedence : int
        {
            DONE = 0, ADD = 5, SUBTRACT = 5, MULTIPLY = 6, DIVIDE = 6, CARAT = 7, MODULO = 6, UNARY = 8,
            NOT = 8, EQ = 3, NE = 3, LT = 4, LE = 4, GT = 4, GE = 4, AND = 2, OR = 1,
            SQRT, EXP, LN, LOG, ABS, SIN, COS, TAN, ASIN, ACOS, ATAN, ATAN2,
            RANDOM, NORMAL, CEIL, FLOOR, ROUND, RAMP, STAGGER, LOGFREQ, STRIDE,
            VDISPLACE, SWIGGLE, CWIGGLE,
            VALUE, ARRAY, PARTARRAYDOUBLE, PARTARRAYINT, SPECARRAY
        };

        // customize by adding a special function

        enum Enum4 { SUM, XMIN, XMAX, AVE, TRAP, SLOPE };

        public const int INVOKED_SCALAR = 1;
        public const int INVOKED_VECTOR = 2;
        public const int INVOKED_ARRAY = 4;
        public const int INVOKED_PER_PARTICLE = 8;
        public const int INVOKED_PER_GRID = 16;
        public const int INVOKED_PER_SURF = 32;

        public const double BIG = 1.0e20;

        private SPARTA sparta;

        private int me;
        private int nvar;                // # of defined variables
        private int maxvar;              // max # of variables following lists can hold
        private string[] names;            // name of each variable
        private int[] style;              // style of each variable
        private int[] num;                // # of values for each variable
        private int[] which;              // next available value for each variable
        private int[] pad;                // 1 = pad loop/uloop variables with 0s, 0 = no pad
        //private VarReader[] reader;   // variable that reads from file
        private VarReader reader;
        private string[][] data;            // str value of each variable's values
        private double[] dvalue;          // single numeric value for internal variables

        private int[] eval_in_progress;   // flag if evaluation of variable is in progress

        RanPark randomequal;     // RNG for equal-style vars
        RanPark randomparticle;  // RNG for particle-style vars

        //private int[] precedence=new int[17];      // precedence level of math operators
        // set length to include up to OR in enum

        private int treestyle;           // tree being used for particle or grid-style var

        // local copies of compute vector_grid vectors
        private int nvec_storage;        // # of vectors currently stored locally
        private int maxvec_storage;      // max # of vectors in vec_storage
        private double[][] vec_storage;    // list of vector copies
        private int[] maxlen_storage;     // allocated length of each vector

        struct Tree
        {            // parse tree for particle-style variables
            double value;          // single scalar  
            double[] array;         // per-atom or per-type list of doubles
            string carray;          // ptr into data struct with nstride = sizeof(struct)
            int type;              // operation, see enum{} in variable.cpp
            int nstride;           // stride between atoms if array is a 2d array
            int selfalloc;         // 1 if array is allocated here, else 0
            int ivalue1, ivalue2;   // extra values for needed for gmask,rmask,grmask
            //Tree* left,*middle,*right;    // ptrs further down tree
        };
        /* ---------------------------------------------------------------------- */

        public Variable(SPARTA sparta)
        {
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);

            nvar = maxvar = 0;
            names = null;
            style = null;
            num = null;
            which = null;
            pad = null;
            reader = null;
            data = null;
            dvalue = null;

            eval_in_progress = null;

            randomequal = null;
            randomparticle = null;

            //precedence[Enum3.DONE] = 0;
            //precedence[OR] = 1;
            //precedence[AND] = 2;
            //precedence[EQ] = precedence[NE] = 3;
            //precedence[LT] = precedence[LE] = precedence[GT] = precedence[GE] = 4;
            //precedence[ADD] = precedence[SUBTRACT] = 5;
            //precedence[MULTIPLY] = precedence[DIVIDE] = precedence[MODULO] = 6;
            //precedence[CARAT] = 7;
            //precedence[UNARY] = precedence[NOT] = 8;

            // local storage of compute vector_grid values
            // stored when a grid-style variable compute uses post_process_grid_flag

            maxvec_storage = 0;
            vec_storage = null;
            maxlen_storage = null;
        }
        public void Set(string[] args)
        {
            Console.WriteLine("variable.set");
            //int narg = args.Length;
            //if (narg < 2) sparta.error.all(FLERR, "Illegal variable command");

            //int replaceflag = 0;

            //// DELETE
            //// doesn't matter if variable no longer exists

            //if (string.Equals(args[1], "delete"))
            //{
            //    if (narg != 2) sparta.error.all(FLERR, "Illegal variable command");
            //    if (Find(args[0]) >= 0) remove(Find(args[0]));
            //    return;

            //    // INDEX
            //    // num = listed args, which = 1st value, data = copied args

            //}
            //else if (string.Equals(args[1], "index"))
            //{
            //    if (narg < 3) sparta.error.all(FLERR, "Illegal variable command");
            //    if (Find(args[0]) >= 0) return;
            //    if (nvar == maxvar) grow();
            //    style[nvar] = (int)Enum1.INDEX;
            //    num[nvar] = narg - 2;
            //    which[nvar] = 0;
            //    pad[nvar] = 0;
            //    data[nvar] = new string[num[nvar]];
            //    Copy(num[nvar], args[2], data[nvar]);

            //    // LOOP
            //    // 1 arg + pad: num = N, which = 1st value, data = single string
            //    // 2 args + pad: num = N2, which = N1, data = single string

            //}
            //else if (strcmp(args[1], "loop") == 0)
            //{
            //    if (Find(args[0]) >= 0) return;
            //    if (nvar == maxvar) grow();
            //    style[nvar] = LOOP;
            //    int nfirst, nlast;
            //    if (narg == 3 || (narg == 4 && strcmp(args[3], "pad") == 0))
            //    {
            //        nfirst = 1;
            //        nlast = atoi(args[2]);
            //        if (nlast <= 0) error.all(FLERR, "Illegal variable command");
            //        if (narg == 4 && strcmp(args[3], "pad") == 0)
            //        {
            //            char digits[12];
            //            sprintf(digits, "%d", nlast);
            //            pad[nvar] = strlen(digits);
            //        }
            //        else pad[nvar] = 0;
            //    }
            //    else if (narg == 4 || (narg == 5 && strcmp(args[4], "pad") == 0))
            //    {
            //        nfirst = atoi(args[2]);
            //        nlast = atoi(args[3]);
            //        if (nfirst > nlast || nlast < 0)
            //            error.all(FLERR, "Illegal variable command");
            //        if (narg == 5 && strcmp(args[4], "pad") == 0)
            //        {
            //            char digits[12];
            //            sprintf(digits, "%d", nlast);
            //            pad[nvar] = strlen(digits);
            //        }
            //        else pad[nvar] = 0;
            //    }
            //    else error.all(FLERR, "Illegal variable command");
            //    num[nvar] = nlast;
            //    which[nvar] = nfirst - 1;
            //    data[nvar] = new char*[1];
            //    data[nvar][0] = NULL;

            //    // WORLD
            //    // num = listed args, which = partition this proc is in, data = copied args
            //    // error check that num = # of worlds in universe

            //}
            //else if (strcmp(args[1], "world") == 0)
            //{
            //    if (narg < 3) error.all(FLERR, "Illegal variable command");
            //    if (Find(args[0]) >= 0) return;
            //    if (nvar == maxvar) grow();
            //    style[nvar] = WORLD;
            //    num[nvar] = narg - 2;
            //    if (num[nvar] != universe.nworlds)
            //        error.all(FLERR, "World variable count doesn't match # of partitions");
            //    which[nvar] = universe.iworld;
            //    pad[nvar] = 0;
            //    data[nvar] = new char*[num[nvar]];
            //    copy(num[nvar], &args[2], data[nvar]);

            //    // UNIVERSE and ULOOP
            //    // for UNIVERSE: num = listed args, data = copied args
            //    // for ULOOP: num = N, data = single string
            //    // which = partition this proc is in
            //    // universe proc 0 creates lock file
            //    // error check that all other universe/uloop variables are same length

            //}
            //else if (strcmp(args[1], "universe") == 0 || strcmp(args[1], "uloop") == 0)
            //{
            //    if (strcmp(args[1], "universe") == 0)
            //    {
            //        if (narg < 3) error.all(FLERR, "Illegal variable command");
            //        if (Find(args[0]) >= 0) return;
            //        if (nvar == maxvar) grow();
            //        style[nvar] = UNIVERSE;
            //        num[nvar] = narg - 2;
            //        pad[nvar] = 0;
            //        data[nvar] = new char*[num[nvar]];
            //        copy(num[nvar], &args[2], data[nvar]);
            //    }
            //    else if (strcmp(args[1], "uloop") == 0)
            //    {
            //        if (narg < 3 || narg > 4 || (narg == 4 && strcmp(args[3], "pad") != 0))
            //            error.all(FLERR, "Illegal variable command");
            //        if (Find(args[0]) >= 0) return;
            //        if (nvar == maxvar) grow();
            //        style[nvar] = ULOOP;
            //        num[nvar] = atoi(args[2]);
            //        data[nvar] = new char*[1];
            //        data[nvar][0] = NULL;
            //        if (narg == 4)
            //        {
            //            char digits[12];
            //            sprintf(digits, "%d", num[nvar]);
            //            pad[nvar] = strlen(digits);
            //        }
            //        else pad[nvar] = 0;
            //    }

            //}
        }

        //public void set(char*, int, string[]);
        //public int next(int, string[]);
        public int find(string name)
        {
            for (int i = 0; i < nvar; i++)
                if (string.Equals(name, names[i])) return i;
            return -1;
        }

        public int equal_style(int ivar)
        {
            if (style[ivar] == (int)Enum1.EQUAL || style[ivar] == (int)Enum1.INTERNAL) return 1;
            return 0;
        }
        //public int particle_style(int);
        //public int grid_style(int);
        //public int surf_style(int);
        //public int internal_style(int);

        //public char* retrieve(char*);
        public double compute_equal(int ivar)
        {
            if (eval_in_progress[ivar]!=0)
                sparta.error.all("Variable has circular dependency");

            eval_in_progress[ivar] = 1;

            double value=0;
            if (style[ivar] ==(int)Enum1.EQUAL) value = evaluate(data[ivar][0], null);
            else if (style[ivar] == (int)Enum1.INTERNAL) value = dvalue[ivar];

            eval_in_progress[ivar] = 0;
            return value;
        }
        public double compute_equal(string str)
        {
            return evaluate(str, null);
        }
        //public void compute_particle(int, double*, int, int);
        //public void compute_grid(int, double*, int, int);
        //public void compute_surf(int, double*, int, int) { }
        //public void internal_set(int, double);

        //public int int_between_brackets(char*&, int);
        //public double evaluate_boolean(char*);

        ////private:


        //private void remove(int);
        //private void grow()
        //{

        //}
        //private void Copy(int narg, string[] from, string[] to)
        //{
        //    Array.Copy(from, to, narg);

        //    //int n;
        //    //for (int i = 0; i < narg; i++)
        //    //{
        //    //    n = strlen(from[i]) + 1;
        //    //    to[i] = new char[n];
        //    //    strcpy(to[i], from[i]);
        //    //}
        //}
        private double evaluate(string str,Tree[] tree=null)
        {
            Console.WriteLine("variable.evaluate");
            return 0;
            //int op, opprevious;
            //double value1, value2;
            //char onechar;
            //string ptr;

            //double[] argstack=new double[MAXLEVEL];
            //Tree[] treestack=new Tree[MAXLEVEL];
            //int[] opstack= new int[MAXLEVEL];
            //int nargstack = 0;
            //int ntreestack = 0;
            //int nopstack = 0;

            //int i = 0;
            //int expect = (int)Enum2.ARG;

            //while (true)
            //{
            //    onechar = str[i];

            //    // whitespace: just skip

            //    if (char.IsWhiteSpace(onechar)) i++;

            //    // ----------------
            //    // parentheses: recursively evaluate contents of parens
            //    // ----------------
            //    else if (onechar=='(')
            //    {
            //        if (expect==(int)Enum2.OP)
            //        {
            //            sparta.error.all("Invalid syntax in variable formula");
            //        }
            //        expect = (int)Enum2.OP;

            //        i = str.IndexOf(')', i);
            //        string contents = str.Substring(str.IndexOf('('), i);
            //        i++;

            //        // evaluate contents and push on stack

            //        if (tree!=null)
            //        {
            //            Tree newtree;
            //            evaluate(contents, newtree);
            //            treestack[ntreestack++] = newtree;
            //        }
            //        else
            //        {
            //            argstack[nargstack++] = evaluate(contents);
            //        }

            //    }
            //}


        }
        //private double collapse_tree(Tree*);
        //private double eval_tree(Tree*, int);
        //private void free_tree(Tree*);
        //private int find_matching_paren(char*, int, char*&);
        //private int math_function(char*, char*, Tree**, Tree**, int &, double*, int &);
        //private int special_function(char*, char*, Tree**, Tree**, int &, double*, int &);
        //private int is_particle_vector(char*);
        //private void particle_vector(char*, Tree**, Tree**, int &);
        //private int is_constant(char*);
        //private double constant(char*);
        //private char* find_next_comma(char*);
        //private void print_tree(Tree*, int);
        //private double* add_storage(double*);


        private class VarReader
        {
            int me, style;
            FileStream fp;

            private SPARTA sparta;

            public VarReader(SPARTA sparta,string str,string file,int flag)
            {
                this.sparta = sparta;
                me = sparta.comm.me;
                style = flag;

                if (me == 0)
                {
                    fp = new FileStream( file, FileMode.Open,FileAccess.Read);
                    if (fp == null)
                    {
                        string strstr=string.Format( "Cannot open file variable file {0}", file);
                        sparta.error.one(strstr);
                    }
                }
                else fp = null;
            }
            //public int read_scalar(string str)
            //{
            //    int n=0;
            //    string ptr;

            //    // read one string from file

            //    if (me == 0)
            //    {
            //        while (true)
            //        {
            //            if (fgets(str, MAXLINE, fp) == NULL) n = 0;
            //            else n = strlen(str);
            //            if (n == 0) break;                                 // end of file
            //            str[n - 1] = '\0';                                   // strip newline
            //            if ((ptr = strchr(str, '#'))) *ptr = '\0';          // strip comment
            //            if (strtok(str, " \t\n\r\f") == NULL) continue;     // skip if blank
            //            n = strlen(str) + 1;
            //            break;
            //        }
            //    }

            //    sparta.mpi.MPI_Bcast(ref n, 1, MPI.MPI_INT, 0, sparta.world);
            //    if (n == 0) return 1;
            //    sparta.mpi.MPI_Bcast(ref str, n, MPI.MPI_CHAR, 0, sparta.world);
            //    return 0;
            //}
        }
    }
}