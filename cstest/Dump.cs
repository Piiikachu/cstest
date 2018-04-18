using System.IO;
using bigint = System.Int64;

namespace cstest
{
    public class Dump
    {
        public const double BIG = 1.0e20;
        public const int IBIG = 2147483647;
        public const double EPSILON = 1.0e-6;

        public const int ONEFIELD = 32;
        public const int DELTA = 1048576;

        enum Enum1{ INT, DOUBLE, BIGINT, STRING };    // many dump files

        enum Enum2{ PERIODIC, OUTFLOW, REFLECT, SURFACE, AXISYM };  // same as Domain

        public string id;                  // user-defined name of Dump
        public string style;               // style of Dump
         
        public int first_flag;            // 0 if no initial dump, 1 if yes initial dump
        public int clearstep;             // 1 if dump invokes computes, 0 if not
         
        public int comm_forward;          // size of forward communication (0 if none)
        public int comm_reverse;          // size of reverse communication (0 if none)
        private SPARTA sparta;
        public Dump(SPARTA sparta, int narg, string[] arg)
        {
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);
            sparta.mpi.MPI_Comm_size(sparta.world, ref nprocs);

            
            id = string.Copy(arg[0]);

            
            style = string.Copy(arg[1]);

            
            filename = string.Copy( arg[4]);

            first_flag = 0;
            flush_flag = 1;

            format = null;
            format_default = null;

            format_line_user = null;
            format_float_user = null;
            format_int_user = null;
            format_bigint_user = null;

            clearstep = 0;
            append_flag = 0;
            buffer_allow = 0;
            buffer_flag = 0;
            padflag = 0;

            maxbuf = 0;
            buf = null;
            maxsbuf = 0;
            sbuf = null;

            // parse filename for special syntax
            // if contains '%', write one file per proc and replace % with proc-ID
            // if contains '*', write one file per timestep and replace * with timestep
            // check file suffixes
            //   if ends in .bin = binary file
            //   else if ends in .gz = gzipped text file
            //   else ASCII text file

            fp = null;
            singlefile_opened = 0;
            compressed = 0;
            binary = 0;
            multifile = 0;

            multiproc = 0;
            nclusterprocs = nprocs;
            filewriter = 0;
            if (me == 0) filewriter = 1;
            fileproc = 0;
            multiname = null;
            if (filename.Contains("%"))
            {
                int index = filename.IndexOf('%');
                multiproc = 1;
                nclusterprocs = 1;
                filewriter = 1;
                fileproc = me;
                sparta.mpi.MPI_Comm_split(sparta.world, me, 0, ref clustercomm);
                multiname = new char[filename.Length + 16].ToString();


                //sprintf(multiname, "%s%d%s", filename, me, ptr + 1);
                multiname = string.Format("{0}{1}{2}", filename, me, filename[index + 1]);
                
            }
            //char* ptr;
            //if ((ptr = strchr(filename, '%')))
            //{
            //    multiproc = 1;
            //    nclusterprocs = 1;
            //    filewriter = 1;
            //    fileproc = me;
            //    sparta.mpi.MPI_Comm_split(sparta.world, me, 0, ref clustercomm);
            //    multiname = new char[filename.Length + 16].ToString();
            //    *ptr = '\0';

            //    //sprintf(multiname, "%s%d%s", filename, me, ptr + 1);
            //    multiname=string.Format("{0}{1}{2}", filename, me, ptr + 1)
            //    *ptr = '%';
            //}

            //if (strchr(filename, '*')) multifile = 1;
            if (filename.Contains("*"))
            {
                multifile = 1;
            }
            System.Console.WriteLine("dump-> suffix");
            //string suffix = filename + filename.Length - ".bin".Length;
            //if (suffix > filename && strcmp(suffix, ".bin") == 0) binary = 1;
            //suffix = filename + strlen(filename) - strlen(".gz");
            //if (suffix > filename && strcmp(suffix, ".gz") == 0) compressed = 1;
        }
        //public void init();
        //public virtual void write();
        public virtual void reset_grid()
        {

        }
        //public void modify_params(int, string*);
        //public virtual bigint memory_usage();


        protected int me, nprocs;             // proc info

        protected string filename;            // user-specified file
        protected int compressed;            // 1 if dump file is written compressed, 0 no
        protected int binary;                // 1 if dump file is written binary, 0 no
        protected int multifile;             // 0 = one big file, 1 = one file per timestep
        protected int multiproc;             // 0 = proc 0 writes for all, 1 = one file/proc
                                    // else # of procs writing files
        protected int nclusterprocs;         // # of procs in my cluster that write to one file
        protected int filewriter;            // 1 if this proc writes a file, else 0
        protected int fileproc;              // ID of proc in my cluster who writes to file
        protected string multiname;           // filename with % converted to cluster ID
        protected int clustercomm;      // MPI communicator within my cluster of procs

        protected int header_flag;           // 0 = item, 2 = xyz
        protected int flush_flag;            // 0 if no flush, 1 if flush every dump
        protected int append_flag;           // 1 if open file in append mode, 0 if not
        protected int buffer_allow;          // 1 if style allows for buffer_flag, 0 if not
        protected int buffer_flag;           // 1 if buffer output as one big string, 0 if not
        protected int padflag;               // timestep padding in filename
        protected int singlefile_opened;     // 1 = one big file, already opened, else 0
         
        protected char[] boundstr=new char[9];          // encoding of boundary flags
         
        protected string format;              // format string for the file write
        protected string format_default;      // default format string

        protected string format_line_user;    // user-specified format strings
        protected string format_float_user;
        protected string format_int_user;
        protected string format_bigint_user;
        protected string[] format_column_user;
         
        protected FileStream fp;                  // file to write dump to
        protected int size_one;              // # of quantities for one entity
        protected int nme;                   // # of entities in this dump from me
        protected int nsme;                  // # of chars in string output from me
         
        protected double boxxlo, boxxhi;      // local copies of domain values
        protected double boxylo, boxyhi;
        protected double boxzlo, boxzhi;
         
        protected bigint ntotal;             // # of per-atom lines in snapshot
         
        protected int maxbuf;                // size of buf
        protected double[] buf;               // memory for dumped quantities
        protected int maxsbuf;               // size of sbuf
        protected string sbuf;                // memory for atom quantities in string format
         
        protected int[] vtype;                // type of each field (INT, DOUBLE, etc)
        protected string[] vformat;            // format string for each field
        
  //      protected int convert_string(int, double*);

  //      virtual void init_style() = 0;
  //virtual void openfile();
  //      virtual int modify_param(int, string*) { return 0; }
  //      virtual void write_header(bigint) = 0;
  //virtual int count() = 0;
  //virtual void pack() = 0;
  //virtual void write_data(int, double*) = 0;
    }
}