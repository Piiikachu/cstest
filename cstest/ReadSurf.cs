using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class ReadSurf
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
        string[] arg;
        public virtual void command(int narg, string[] args)
        {
            arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);
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


            parse_keyword(0);


            // close file

            if (me == 0)
            {
                //if (compressed!=0) fp.Close();
                //else fclose(fp);
                fp.Close();
            }

            // apply optional keywords for geometric transformations
            // store optional keywords for group and type information
            // store optional keyword for file output

            origin[0] = origin[1] = origin[2] = 0.0;
            int grouparg = 0;
            int typeadd = 0;
            int partflag =(int) Enum2.NONE;
            int filearg = 0;

            int iarg = 1;

            while (iarg < narg)
            {
                switch (arg[iarg])
                {
                    case "origion":
                        if (iarg + 4 > narg) sparta.error.all("Invalid read_surf command");
                        double ox = double.Parse(arg[iarg + 1]);
                        double oy = double.Parse(arg[iarg + 2]);
                        double oz = double.Parse(arg[iarg + 3]);
                        if (dim == 2 && oz != 0.0)
                            sparta.error.all("Invalid read_surf geometry transformation for 2d simulation");
                        origin[0] = ox;
                        origin[1] = oy;
                        origin[2] = oz;
                        iarg += 4;
                        break;
                    case "trans":
                        if (iarg + 4 > narg) sparta.error.all("Invalid read_surf command");
                        double dx = double.Parse(arg[iarg + 1]);
                        double dy = double.Parse(arg[iarg + 2]);
                        double dz = double.Parse(arg[iarg + 3]);
                        if (dim == 2 && dz != 0.0)
                            sparta.error.all("Invalid read_surf geometry transformation for 2d simulation");
                        origin[0] += dx;
                        origin[1] += dy;
                        origin[2] += dz;
                        translate(dx, dy, dz);
                        iarg += 4;
                        break;
                    case "atrans":
                        if (iarg + 4 > narg) sparta.error.all("Invalid read_surf command");
                        double ax = double.Parse(arg[iarg + 1]);
                        double ay = double.Parse(arg[iarg + 2]);
                        double az = double.Parse(arg[iarg + 3]);
                        if (dim == 2 && az != 0.0)
                            sparta.error.all("Invalid read_surf geometry transformation for 2d simulation");
                        double dx1 = ax - origin[0];
                        double dy1 = ay - origin[1];
                        double dz1 = az - origin[2];
                        origin[0] = ax;
                        origin[1] = ay;
                        origin[2] = az;
                        translate(dx1, dy1, dz1);
                        iarg += 4;
                        break;
                    case "ftrans":
                        if (iarg + 4 > narg) sparta.error.all("Invalid read_surf command");
                        double fx = double.Parse(arg[iarg + 1]);
                        double fy = double.Parse(arg[iarg + 2]);
                        double fz = double.Parse(arg[iarg + 3]);
                        if (dim == 2 && fz != 0.5)
                            sparta.error.all("Invalid read_surf geometry transformation for 2d simulation");
                        double ax1 = sparta.domain.boxlo[0] + fx * sparta.domain.xprd;
                        double ay1 = sparta.domain.boxlo[1] + fy * sparta.domain.yprd;
                        double az1=0.0;
                        if (dim == 3) az = sparta.domain.boxlo[2] + fz * sparta.domain.zprd;
                        double dx2 = ax1 - origin[0];
                        double dy2 = ay1 - origin[1];
                        double dz2 = az1 - origin[2];
                        origin[0] = ax1;
                        origin[1] = ay1;
                        origin[2] = az1;
                        translate(dx2, dy2, dz2);
                        iarg += 4;
                        break;
                    case "scale":
                        if (iarg + 4 > narg) sparta.error.all("Invalid read_surf command");
                        double sx = double.Parse(arg[iarg + 1]);
                        double sy = double.Parse(arg[iarg + 2]);
                        double sz = double.Parse(arg[iarg + 3]);
                        if (dim == 2 && sz != 1.0)
                            sparta.error.all("Invalid read_surf geometry transformation for 2d simulation");
                        scale(sx, sy, sz);
                        iarg += 4;
                        break;
                    case "rotate":
                        if (iarg + 5 > narg) sparta.error.all("Invalid read_surf command");
                        double theta = double.Parse(arg[iarg + 1]);
                        double rx = double.Parse(arg[iarg + 2]);
                        double ry = double.Parse(arg[iarg + 3]);
                        double rz = double.Parse(arg[iarg + 4]);
                        if (dim == 2 && (rx != 0.0 || ry != 0.0 || rz != 1.0))
                            sparta.error.all("Invalid read_surf geometry transformation for 2d simulation");
                        if (rx == 0.0 && ry == 0.0 && rz == 0.0)
                            sparta.error.all("Invalid read_surf geometry transformation for 2d simulation");
                        rotate(theta, rx, ry, rz);
                        iarg += 5;
                        break;
                    case "invert":
                        invert();
                        iarg += 1;
                        break;
                    case "clip":
                        double frac = 0.0;
                        if (iarg + 1 < narg)
                        {
                            char c = arg[iarg + 1][0];
                            if (char.IsDigit(c) || c == '-' || c == '+' || c == '.')
                            {
                                frac = double.Parse(arg[iarg + 1]);
                                if (frac < 0.0 || frac >= 0.5)
                                    sparta.error.all("Invalid read_surf command");
                                iarg++;
                            }
                        }
                        if (frac > 0.0) push_points_to_boundary(frac);
                        if (dim == 2) clip2d();
                        else clip3d();
                        iarg++;
                        break;
                    case "group":
                        if (iarg + 2 > narg) sparta.error.all("Invalid read_surf command");
                        grouparg = iarg + 1;
                        iarg += 2;
                        break;
                    case "typeadd":
                        if (iarg + 2 > narg) sparta.error.all("Invalid read_surf command");
                        typeadd = int.Parse(arg[iarg + 1]);
                        if (typeadd < 0) sparta.error.all("Invalid read_surf command");
                        iarg += 2;
                        break;
                    case "particle":
                        if (iarg + 2 > narg) sparta.error.all("Invalid read_surf command");
                        if (string.Equals(arg[iarg + 1], "none")) partflag = (int)Enum2.NONE;
                        else if (string.Equals(arg[iarg + 1], "check")) partflag = (int)Enum2.CHECK;
                        else if (string.Equals(arg[iarg + 1], "keep")) partflag = (int)Enum2.KEEP;
                        else sparta.error.all("Invalid read_surf command");
                        iarg += 2;
                        break;
                    case "file":
                        if (iarg + 2 > narg) sparta.error.all("Invalid read_surf command");
                        filearg = iarg + 1;
                        iarg += 2;
                        break;
                    default:
                        sparta.error.all("Invalid read_surf command");
                        break;
                }
            }


        }


        protected int me;
        protected string line,keyword,buffer;
        protected FileStream fp;
        protected int compressed;
         
        protected int dim, isc;
        protected double[] origin=new double[3];
         
        protected List<Surf.Point> pts;
        protected List<Surf.Line> lines;
        protected List<Surf.Tri> tris;
        protected int npoint_old, nline_old, ntri_old;
        protected int npoint_new, nline_new, ntri_new;
        protected int maxpoint, maxline, maxtri;
         
        protected int[,] edge;
        protected int nedge, maxedge;
        //todo: fix me first
        protected void header()
        {
            int n=0;
            using (StreamReader sr = new StreamReader(fp))
            {
                if (me == 0)
                {
                    string eof =sr.ReadLine();
                    if (eof == null)
                    {
                        sparta.error.one("Unexpected end of data file");
                    }
                }

                npoint_new = nline_new = ntri_new = 0;

            
                string line;


                while ((line = sr.ReadLine()) != null)
                {
                    n = line.Length + 1;
                    sparta.mpi.MPI_Bcast(ref n, 1, MPI.MPI_INT, 0, sparta.world);
                    sparta.mpi.MPI_Bcast(ref line, 1, MPI.MPI_INT, 0, sparta.world);
                    if (line.StartsWith("#"))
                    {
                        continue;
                    }
                    if (string.IsNullOrWhiteSpace(line))
                    {
                        continue;
                    }
                    if (line.Contains("points"))
                    {
                        npoint_new = ParseInt(line);
                    }
                    else if (line.Contains("lines"))
                    {
                        if (dim == 3)
                        {
                            sparta.error.all("Surf file cannot contain lines for 3d simulation");
                        }
                        nline_new = ParseInt(line);
                    }
                    else if (line.Contains("triangles"))
                    {
                        if (dim == 2)
                        {
                            sparta.error.all("Surf file cannot contain triangles for 2d simulation");
                        }
                        ntri_new = ParseInt(line);
                    }
                    else break;

                }
                 
                
            }
            if (npoint_new == 0) sparta.error.all( "Surf file does not contain points");
            if (dim == 2 && nline_new == 0)
                sparta.error.all( "Surf file does not contain lines");
            if (dim == 3 && ntri_new == 0)
                sparta.error.all( "Surf file does not contain triangles");


        }
        protected int read_points(string line)
        {
            if (string.IsNullOrWhiteSpace(line))
            {
                return 1;
            }
            if (line.Contains("Lines")||line.Contains("Triangles"))
            {
                return 0;
            }
            string[] words = line.Split();
            if (dim == 2 && words.Length != 3)
            {
                sparta.error.all("Incorrect point format in surf file");
            }
            if (dim == 3 && words.Length != 4)
            {
                sparta.error.all("Incorrect point format in surf file");
            }
            double[] x = new double[3];
            x[0] = double.Parse(words[1]);
            x[1] = double.Parse(words[2]);
            if (dim==3)
            {
                x[2] = double.Parse(words[3]);
            }
            else
            {
                x[2] = 0.0;
            }
            Surf.Point pt = new Surf.Point()
            {
                x = x
            };
            pts.Add(pt);

            return 1;

        }
        protected int read_lines(string line)
        {
            if (string.IsNullOrWhiteSpace(line))
            {
                return 1;
            }
            if (line.Contains("Triangles"))
            {
                return 0;
            }
            string[] words = line.Split();
            int type, p1, p2;


            if (words.Length != 3 && words.Length != 4) 
            {
                sparta.error.all("Incorrect line format in surf file");
            }
            int typeflag = 0;
            if (words.Length==4)
            {
                typeflag = 1;
            }
            if (typeflag != 0)
            {
                type = sparta.input.inumeric(words[1]);
            }
            else
            {
                type = 1;
            }
            p1 = int.Parse(words[1]);
            p2 = int.Parse(words[2]);
            if (p1 < 1 || p1 > npoint_new || p2 < 1 || p2 > npoint_new || p1 == p2)
            {
                sparta.error.all("Invalid point index in line");
            }
            Surf.Line aline = new Surf.Line();
            aline.type = type;
            aline.mask = 1;
            aline.isc = aline.isr = -1;
            aline.p1 = p1 - 1 + npoint_old;
            aline.p2 = p2 - 1 + npoint_old;
            lines.Add(aline);
            return 1;
        }

        
        protected int read_tris(string line)
        {
            if (string.IsNullOrWhiteSpace(line))
            {
                return 1;
            }

            string[] words = line.Split();
            int type, p1, p2, p3;

            if (words.Length!=4||words.Length!=5)
            {
                sparta.error.all("Incorrect line format in surf file");
            }
            int typeflag = 0;
            if (words.Length==5)
            {
                typeflag = 1;
            }
            if (typeflag!=0)
            {
                type = sparta.input.inumeric(words[1]);
            }
            else
            {
                type = 1;
            }
            p1 = int.Parse(words[1]);
            p2 = int.Parse(words[2]);
            p3 = int.Parse(words[3]);
            if (p1 < 1 || p1 > npoint_new || p2 < 1 || p2 > npoint_new || p3 < 1 || p3 > npoint_new || p1 == p2 || p2 == p3)
            {
                sparta.error.all("Invalid point index in triangle");
            }
            Surf.Tri tri = new Surf.Tri();
            tri.type = type;
            tri.mask = 1;
            tri.isc = tri.isr = -1;
            tri.p1 = p1 - 1 + npoint_old;
            tri.p2 = p2 - 1 + npoint_old;
            tri.p3 = p3 - 1 + npoint_old;
            tris.Add(tri);
            return 1;
        }

        protected void translate(double dx, double dy, double dz)
        {
            int m = npoint_old;
            for (int i = 0; i < npoint_new; i++)
            {
                pts[m].x[0] += dx;
                pts[m].x[1] += dy;
                pts[m].x[2] += dz;
                m++;
            }
        }
        protected void scale(double sx , double sy, double sz)
        {
            int m = npoint_old;
            for (int i = 0; i < npoint_new; i++)
            {
                pts[m].x[0] = sx * (pts[m].x[0] - origin[0]) + origin[0];
                pts[m].x[1] = sy * (pts[m].x[1] - origin[1]) + origin[1];
                if (dim == 3) pts[m].x[2] = sz * (pts[m].x[2] - origin[2]) + origin[2];
                m++;
            }
        }
        protected void rotate(double theta, double rx, double ry, double rz)
        {
            double[] r = new double[3];
            double[] q = new double[4];
            double[] d = new double[3];
            double[] dnew = new double[3];
            double[,] rotmat = new double[3, 3];

            theta *= MyConst.MY_PI / 180.0;
            r[0] = rx;
            r[1] = ry;
            r[2] = rz;
            MathExtra.norm3(r);
            MathExtra.axisangle_to_quat(r, theta, q);
            MathExtra.quat_to_mat(q, rotmat);

            int m = npoint_old;
            for (int i = 0; i < npoint_new; i++)
            {
                d[0] = pts[m].x[0] - origin[0];
                d[1] = pts[m].x[1] - origin[1];
                d[2] = pts[m].x[2] - origin[2];
                MathExtra.matvec(rotmat, d, dnew);
                pts[m].x[0] = dnew[0] + origin[0];
                pts[m].x[1] = dnew[1] + origin[1];
                if (dim == 3) pts[m].x[2] = dnew[2] + origin[2];
                m++;
            }
            
        }
        protected void invert()
        {
          
            if (dim == 2)
            {
                int m = nline_old;

                for (int i = 0; i < nline_new; i++)
                {
                    Surf.Line aline = new Surf.Line();
                    int tmp1 = lines[m].p1;
                    int tmp2= lines[m].p2;
                    aline = lines[m];
                    aline.p1 = tmp2;
                    aline.p2 = tmp1;
                    lines[m] = aline;
                    m++;
                }
            }

            if (dim == 3)
            {
                int m = ntri_old;
                for (int i = 0; i < ntri_new; i++)
                {
                    Surf.Tri tri = new Surf.Tri();
                    int tmp1 = tris[m].p2;
                    int tmp2 = tris[m].p3;
                    tri = tris[m];
                    tri.p2 = tmp2;
                    tri.p3 = tmp1;
                    tris[m] = tri;
                    m++;
                    
                }
            }
        }
        protected void clip2d()
        {
            Console.WriteLine("ReadSurf.clip2d");
        }
        protected void clip3d()
        {
            Console.WriteLine("ReadSurf.clip3d");
        }

        protected void push_points_to_boundary(double frac)
        {
            double[] x;

            double[] boxlo =sparta.domain.boxlo;
            double[] boxhi =sparta.domain.boxhi;

            double xdelta = frac * (boxhi[0] - boxlo[0]);
            double ydelta = frac * (boxhi[1] - boxlo[1]);
            double zdelta = frac * (boxhi[2] - boxlo[2]);

            int n = npoint_old + npoint_new;
            for (int i = npoint_old; i < n; i++)
            {
                x = pts[i].x;
                if (x[0] >= boxlo[0] && x[0] <= boxhi[0])
                {
                    if (x[0] - boxlo[0] < xdelta) x[0] = boxlo[0];
                    else if (boxhi[0] - x[0] < xdelta) x[0] = boxhi[0];
                }
                if (x[1] >= boxlo[1] && x[1] <= boxhi[1])
                {
                    if (x[1] - boxlo[1] < ydelta) x[1] = boxlo[1];
                    else if (boxhi[1] - x[1] < ydelta) x[1] = boxhi[1];
                }
                if (dim == 2) continue;
                if (x[2] >= boxlo[2] && x[2] <= boxhi[2])
                {
                    if (x[2] - boxlo[2] < zdelta) x[2] = boxlo[2];
                    else if (boxhi[2] - x[2] < zdelta) x[2] = boxhi[2];
                }
            }
        }
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
        //todo: it has its own open method (to open compressed file)
        int ParseInt(string str)
        {
            int n = str.Length;
            char[] tmp = new char[n];
            int j = 0;
            for (int i = 0; i < n; i++)
            {
                if (char.IsDigit(str[i]))
                {
                    tmp[j] = str[i];
                    j++;
                }
                
            }
            return int.Parse(new string(tmp));
        }
        protected void open(string file)
        {

        }
        protected void parse_keyword(int first)
        {
            //int eof = 0;

            // proc 0 reads upto non-blank line plus 1 following line
            // eof is set to 1 if any read hits end-of-file

            if (me == 0)
            {
                try
                {
                    fp = new FileStream(arg[0], FileMode.Open, FileAccess.Read);
                }
                catch (Exception e)
                {

                    Console.WriteLine(e.Message.ToString());
                }
                StreamReader sr = new StreamReader(fp);
                string line;
                int startflag = 0;
                if (first==0)
                {
                    if (dim==2)
                    {
                        while ((line = sr.ReadLine()) != null)
                        {
                            if (line.Contains("Lines"))
                            {
                                keyword = "Lines";
                                startflag = 1;
                                continue;
                            }
                            if (startflag==1)
                            {
                                startflag = read_lines(line);
                            }

                        }
                        if (!keyword.Equals("Lines"))
                        {
                            sparta.error.all("Read_surf did not find lines section of surf file");
                        }
                        if (lines.Count != nline_new)
                        {
                            sparta.error.all("Incorrect line number in surf file");
                        }
                        else
                        {
                            if (me == 0)
                            {
                                if (sparta.screen != null)
                                {
                                    Console.WriteLine("  {0} lines\n", nline_new);
                                    new StreamWriter(sparta.screen).WriteLine("  {0} lines\n", nline_new);
                                    
                                }

                                if (sparta.logfile != null)
                                {
                                    new StreamWriter(sparta.logfile).WriteLine("  {0} lines\n", nline_new);
                                    
                                }
                            }
                        }
                    }
                    else if (dim==3)
                    {
                        while ((line = sr.ReadLine()) != null)
                        {
                            if (line.Contains("Triangles"))
                            {
                                keyword = "Triangles";
                                startflag = 1;
                                continue;
                                
                            }
                            if (startflag == 1)
                            {
                                startflag = read_tris(line);
                            }

                        }
                        if (!keyword.Equals("Triangles"))
                        {
                            sparta.error.all("Read_surf did not find triangles section of surf file");
                        }
                        if (tris.Count != ntri_new)
                        {
                            sparta.error.all("Incorrect triangle number in surf file");
                        }
                        else
                        {
                            if (me == 0)
                            {
                                if (sparta.screen != null)
                                {
                                    Console.WriteLine("  {0} triangles\n", nline_new);
                                    new StreamWriter(sparta.screen).WriteLine("  {0} triangles\n", nline_new);
                                }

                                if (sparta.logfile != null)
                                {
                                    new StreamWriter(sparta.logfile).WriteLine("  {0} triangles\n", nline_new);
                                }
                            }
                        }
                    }
                }
                else if (first==1)
                {
                    while ((line=sr.ReadLine())!=null)
                    {
                        if (line.Contains("Points"))
                        {
                            keyword = "Points";
                            startflag = 1;
                            continue;
                        }
                        if (startflag==1)
                        {
                            startflag=read_points(line);
                        }
                    }
                    if (!keyword.Equals("Points"))
                    {
                        sparta.error.all("Read_surf did not find points section of surf file");
                    }
                    if (pts.Count!=npoint_new)
                    {
                        sparta.error.all("Incorrect point number in surf file");
                    }
                    else
                    {
                        if (me==0)
                        {
                            if (sparta.screen!=null)
                            {
                                Console.WriteLine("  {0} points\n", npoint_new);
                                new StreamWriter(sparta.screen).WriteLine("  {0} points\n", npoint_new);
                            }

                            if (sparta.logfile!=null)
                            {
                                new StreamWriter(sparta.logfile).WriteLine("  {0} points\n", npoint_new);
                            }
                        }
                    }
                }
                
            }

            //// if eof, set keyword empty and return

            //MPI_Bcast(&eof, 1, MPI_INT, 0, world);
            //if (eof)
            //{
            //    keyword[0] = '\0';
            //    return;
            //}

            //// bcast keyword line to all procs

            //int n;
            //if (me == 0) n = strlen(line) + 1;
            //MPI_Bcast(&n, 1, MPI_INT, 0, world);
            //MPI_Bcast(line, n, MPI_CHAR, 0, world);

            //// copy non-whitespace portion of line into keyword

            //int start = strspn(line, " \t\n\r");
            //int stop = strlen(line) - 1;
            //while (line[stop] == ' ' || line[stop] == '\t'
            //   || line[stop] == '\n' || line[stop] == '\r') stop--;
            //line[stop + 1] = '\0';
            //strcpy(keyword, &line[start]);
        }
        //protected int count_words(char*);

        protected virtual void grow_surf()
        {
            pts = new List<Surf.Point>(maxpoint);
            lines = new List<Surf.Line>(maxline);
            tris = new List<Surf.Tri>(maxtri);
        }
    }
}
