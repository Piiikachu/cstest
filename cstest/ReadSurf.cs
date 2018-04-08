using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    class ReadSurf
    {
        enum Enum1{ NEITHER, BAD, GOOD };
        enum Enum2{ NONE, CHECK, KEEP };
        enum Enum3{ UNKNOWN, OUTSIDE, INSIDE, OVERLAP };           // several files

        public const int MAXLINE = 256;
        public const int CHUNK = 1024;
        public const double EPSILON_NORM = 1.0e-12;
        public const double EPSILON_GRID = 1.0e-3;
        public const double BIG = 1.0e20;
        public const int DELTA = 128;           // must be 2 or greater 

        private SPARTA sparta;

        public ReadSurf(SPARTA sparta)
        {
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);
            //line = new char[MAXLINE];
            //keyword = new char[MAXLINE];
            //buffer = new char[CHUNK * MAXLINE];
        }
        public virtual void command(int narg, string[] arg)
        {
            if (sparta.grid.exist==0)
                sparta.error.all( "Cannot read_surf before grid is defined");

            sparta.surf.exist = 1;
            dim = sparta.domain.dimension;

            if (narg < 1) sparta.error.all( "Illegal read_surf command");

            // read header info

            if (me == 0)
            {
                if (sparta.screen!=null)
                {
                    Console.WriteLine("Reading surf file ...\n");
                    new StreamWriter(sparta.screen).Write("Reading surf file ...\n");
                }
                try
                {
                    fp = new FileStream(arg[0], FileMode.Open, FileAccess.Read);
                }
                catch (Exception e)
                {

                    Console.WriteLine(e.Message.ToString());
                }
                //open(arg[0]);
            }

            sparta.mpi.MPI_Barrier(sparta.world);
            double time1 = sparta.mpi.MPI_Wtime();

            header();

            // extend pts,lines,tris data structures

            pts = sparta.surf.pts;
            lines = sparta.surf.lines;
            tris = sparta.surf.tris;

            npoint_old = sparta.surf.npoint;
            nline_old = sparta.surf.nline;
            ntri_old = sparta.surf.ntri;

            maxpoint = npoint_old + npoint_new;
            maxline = nline_old + nline_new;
            maxtri = ntri_old + ntri_new;

            grow_surf();

            // read and store Points and Lines/Tris sections

            parse_keyword(1);
            if (strcmp(keyword, "Points") != 0)
                sparta.error.all(
                       "Read_surf did not find points section of surf file");
            read_points();

            parse_keyword(0);
            if (dim == 2)
            {
                if (strcmp(keyword, "Lines") != 0)
                    sparta.error.all(
                         "Read_surf did not find lines section of surf file");
                read_lines();
            }
            else
            {
                if (strcmp(keyword, "Triangles") != 0)
                    sparta.error.all(
                           "Read_surf did not find triangles section of surf file");
                read_tris();
            }

            // close file

            if (me == 0)
            {
                if (compressed) pclose(fp);
                else fclose(fp);
            }

            // apply optional keywords for geometric transformations
            // store optional keywords for group and type information
            // store optional keyword for file output
        }


        protected int me;
        protected string line,keyword,buffer;
        protected FileStream fp;
        protected int compressed;
         
        protected int dim, isc;
        protected double[] origin=new double[3];
         
        protected Surf.Point[] pts;
        protected Surf.Line[] lines;
        protected Surf.Tri[] tris;
        protected int npoint_old, nline_old, ntri_old;
        protected int npoint_new, nline_new, ntri_new;
        protected int maxpoint, maxline, maxtri;
         
        protected int[,] edge;
        protected int nedge, maxedge;
        //todo: fix me first
        protected void header()
        {
            int n;
            
            if (me==0)
            {
                string eof = new StreamReader(fp).ReadLine();
                if (eof==null)
                {
                    sparta.error.one("Unexpected end of data file");
                }
            }

            npoint_new = nline_new = ntri_new = 0;

            while (true)
            {
                n = 0;
                if (me==0)
                {
                    string line = new StreamReader(fp).ReadLine();
                    if (line==null)
                    {
                        n = 0;
                    }
                    else
                    {
                        n = line.Length + 1;
                    }
                }

                sparta.mpi.MPI_Bcast(ref n, 1, MPI.MPI_INT, 0, sparta.world);


            }


        }
        //protected void read_points();
        //protected void read_lines();
        //protected void read_tris();

        //protected void translate(double, double, double);
        //protected void scale(double, double, double);
        //protected void rotate(double, double, double, double);
        //protected void invert();
        //protected void clip2d();
        //protected void clip3d();

        //protected void push_points_to_boundary(double);
        //protected void check_neighbor_norm_2d();
        //protected void check_neighbor_norm_3d();
        //protected void check_point_near_surf_2d();
        //protected void check_point_near_surf_3d();

        //protected void point_line_compare(int, Surf::Line*, double, int &, int &);
        //protected void point_tri_compare(int, Surf::Tri*, double, int &, int &, int, int, int);

        //protected int find_edge(int, int);
        //protected void add_edge(int, int, int);

        //protected double shortest_line();
        //protected void smallest_tri(double &, double &);

        //protected void open(char*);
        //protected void parse_keyword(int);
        //protected int count_words(char*);

        //protected virtual void grow_surf();
    }
}
