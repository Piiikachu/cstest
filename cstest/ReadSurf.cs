using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using bigint = System.Int64;

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
            if (sparta.grid.exist == 0)
                sparta.error.all("Cannot read_surf before grid is defined");

            sparta.surf.exist = 1;
            dim = sparta.domain.dimension;

            if (narg < 1) sparta.error.all("Illegal read_surf command");

            // read header info

            if (me == 0)
            {
                if (sparta.screen != null)
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
            int partflag = (int)Enum2.NONE;
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
                        double az1 = 0.0;
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
            // error test on particles

            if (sparta.particle.exist != 0 && partflag == (int)Enum2.NONE)
                sparta.error.all("Using read_surf particle none when particles exist");

            // if specified, apply group and typeadd keywords
            // these reset per-element mask/type info
            if (grouparg != 0)
            {
                int igroup = sparta.surf.find_group(arg[grouparg]);
                if (igroup < 0) igroup = sparta.surf.add_group(arg[grouparg]);
                int groupbit = sparta.surf.bitmask[igroup];
                if (dim == 2)
                {
                    int n = nline_old + nline_new;
                    for (int i = nline_old; i < n; i++)
                    {
                        Surf.Line line = lines[i];
                        line.mask |= groupbit;
                        //lines[i].mask |= groupbit;
                        lines[i] = line;

                    }
                }
                else
                {
                    int n = ntri_old + ntri_new;
                    for (int i = ntri_old; i < n; i++)
                    {
                        Surf.Tri tri = tris[i];
                        tri.mask |= groupbit;
                        tris[i] = tri;

                    }
                }
            }

            if (typeadd != 0)
            {
                if (dim == 2)
                {
                    int n = nline_old + nline_new;
                    for (int i = nline_old; i < n; i++)
                    {
                        Surf.Line line = lines[i];
                        line.type += typeadd;
                        lines[i] = line;
                    }
                }
                else
                {
                    int n = ntri_old + ntri_new;
                    for (int i = ntri_old; i < n; i++)
                    {
                        Surf.Tri tri = tris[i];
                        tri.type += typeadd;
                        tris[i] = tri;
                    }
                }
            }
            // update Surf data structures

            sparta.surf.pts = pts;
            sparta.surf.lines = lines;
            sparta.surf.tris = tris;

            sparta.surf.npoint = npoint_old + npoint_new;
            sparta.surf.nline = nline_old + nline_new;
            sparta.surf.ntri = ntri_old + ntri_new;

            // extent of surf after geometric transformations
            // compute sizes of smallest surface elements

            double[,] extent = new double[3, 2];
            extent[0, 0] = extent[1, 0] = extent[2, 0] = BIG;
            extent[0, 1] = extent[1, 1] = extent[2, 1] = -BIG;

            int m = npoint_old;
            for (int i = 0; i < npoint_new; i++)
            {
                extent[0, 0] = Math.Min(extent[0, 0], pts[m].x[0]);
                extent[0, 1] = Math.Max(extent[0, 1], pts[m].x[0]);
                extent[1, 0] = Math.Min(extent[1, 0], pts[m].x[1]);
                extent[1, 1] = Math.Max(extent[1, 1], pts[m].x[1]);
                extent[2, 0] = Math.Min(extent[2, 0], pts[m].x[2]);
                extent[2, 1] = Math.Max(extent[2, 1], pts[m].x[2]);
                m++;
            }

            double minlen = 0, minarea = 0;
            if (dim == 2) minlen = shortest_line();
            if (dim == 3) smallest_tri(out minlen, out minarea);

            if (me == 0)
            {

                if (sparta.screen != null)
                {
                    string str1 = string.Format("  {0:G6} {1:G6} xlo xhi\n", extent[0, 0], extent[0, 1]);
                    string str2 = string.Format("  {0:G6} {1:G6} ylo yhi\n", extent[1, 0], extent[1, 1]);
                    string str3 = string.Format("  {0:G6} {1:G6} zlo zhi\n", extent[2, 0], extent[2, 1]);
                    if (dim == 2)
                    {
                        string str4 = string.Format("  {0} min line length\n", minlen);
                        Console.WriteLine(str1 + str2 + str3 + str4);
                        new StreamWriter(sparta.screen).WriteLine(str1 + str2 + str3 + str4);
                    }
                    if (dim == 3)
                    {
                        string str4 = string.Format("  {0} min triangle edge length\n", minlen);
                        string str5 = string.Format("  {0} min triangle area\n", minarea);
                        Console.WriteLine(str1 + str2 + str3 + str4 + str5);
                        new StreamWriter(sparta.screen).WriteLine(str1 + str2 + str3 + str4 + str5);
                    }
                }
                if (sparta.logfile != null)
                {
                    string str1 = string.Format("  {0:G6} {1:G6} xlo xhi\n", extent[0, 0], extent[0, 1]);
                    string str2 = string.Format("  {0:G6} {1:G6} ylo yhi\n", extent[1, 0], extent[1, 1]);
                    string str3 = string.Format("  {0:G6} {1:G6} zlo zhi\n", extent[2, 0], extent[2, 1]);
                    if (dim == 2)
                    {
                        string str4 = string.Format("  {0:G6} min line length\n", minlen);
                        new StreamWriter(sparta.logfile).WriteLine(str1 + str2 + str3 + str4);
                    }
                    if (dim == 3)
                    {
                        string str4 = string.Format("  {0} min triangle edge length\n", minlen);
                        string str5 = string.Format("  {0} min triangle area\n", minarea);
                        new StreamWriter(sparta.logfile).WriteLine(str1 + str2 + str3 + str4 + str5);
                    }
                }
            }
            // compute normals of new lines or triangles

            if (dim == 2) sparta.surf.compute_line_normal(nline_old, nline_new);
            else sparta.surf.compute_tri_normal(ntri_old, ntri_new);

            // error check on new points,lines,tris
            // all points must be inside or on surface of simulation box

            sparta.surf.check_point_inside(npoint_old, npoint_new);

            // write out new surf file if requested
            // do this before assigning surfs to grid cells, in case an error occurs
            if (filearg != 0)
            {
                WriteSurf wf = new WriteSurf(sparta);
                if (sparta.comm.me == 0)
                {
                    FileStream fp = new FileStream(arg[filearg], FileMode.Open, FileAccess.Write);
                    if (fp == null)
                    {
                        string str = string.Format("Cannot open surface file {0}", arg[0]);
                        sparta.error.one(str);
                    }
                    wf.write_file(fp);
                    fp.Close();
                }
                //delete wf;
            }

            // -----------------------
            // map surfs to grid cells
            // -----------------------

            sparta.mpi.MPI_Barrier(sparta.world);
            double time2 = sparta.mpi.MPI_Wtime();

            // sort particles

            if (sparta.particle.exist != 0) sparta.particle.sort();

            // make list of surf elements I own
            // assign surfs to grid cells
            // error checks to flag bad surfs

            sparta.surf.setup_surf();

            sparta.grid.unset_neighbors();
            sparta.grid.remove_ghosts();

            if (sparta.particle.exist != 0 && sparta.grid.nsplitlocal != 0)
            {
                Grid.ChildCell[] cells = sparta.grid.cells;
                int nglocal = sparta.grid.nlocal;
                for (int icell = 0; icell < nglocal; icell++)
                    if (cells[icell].nsplit > 1)
                        sparta.grid.combine_split_cell_particles(icell, 1);
            }

            sparta.grid.clear_surf();

            sparta.mpi.MPI_Barrier(sparta.world);
            double time3 = sparta.mpi.MPI_Wtime();

            // error checks that can be done before surfs are mapped to grid cells

            if (dim == 2)
            {
                sparta.surf.check_watertight_2d(npoint_old, nline_old);
                check_neighbor_norm_2d();
            }
            else
            {
                sparta.surf.check_watertight_3d(npoint_old, ntri_old);
                check_neighbor_norm_3d();
            }

            sparta.mpi.MPI_Barrier(sparta.world);
            double time4 = sparta.mpi.MPI_Wtime();

            // map surfs to grid cells then error check
            // check done on per-grid-cell basis, too expensive to do globally

            sparta.grid.surf2grid(1, 1);

            if (dim == 2) check_point_near_surf_2d();
            else check_point_near_surf_3d();

            sparta.mpi.MPI_Barrier(sparta.world);
            double time5 = sparta.mpi.MPI_Wtime();

            // re-setup grid ghosts and neighbors

            sparta.grid.setup_owned();
            sparta.grid.acquire_ghosts();
            sparta.grid.reset_neighbors();
            sparta.comm.reset_neighbors();

            sparta.mpi.MPI_Barrier(sparta.world);
            double time6 = sparta.mpi.MPI_Wtime();

            // flag cells and corners as OUTSIDE or INSIDE

            sparta.grid.set_inout();
            sparta.grid.type_check();

            // DEBUG
            //sparta.grid.debug();

            sparta.mpi.MPI_Barrier(sparta.world);
            double time7 = sparta.mpi.MPI_Wtime();

            // remove particles in any cell that is now INSIDE or has new surfs
            // reassign particles in split cells to sub cell owner
            // compress particles if any flagged for deletion

            bigint ndeleted = 0;
            if (sparta.particle.exist != 0)
            {
                Grid.ChildCell[] cells = sparta.grid.cells;
                Grid.ChildInfo[] cinfo = sparta.grid.cinfo;
                int nglocal = sparta.grid.nlocal;
                int delflag = 0;

                for (int icell = 0; icell < nglocal; icell++)
                {
                    if (cinfo[icell].type == (int)Enum3.INSIDE)
                    {
                        if (partflag == (int)Enum2.KEEP)
                            sparta.error.one("Particles are inside new surfaces");
                        if (cinfo[icell].count != 0) delflag = 1;
                        sparta.particle.remove_all_from_cell(cinfo[icell].first);
                        cinfo[icell].count = 0;
                        cinfo[icell].first = -1;
                        continue;
                    }
                    if (cells[icell].nsurf != 0 && cells[icell].nsplit >= 1)
                    {
                        int nsurf = cells[icell].nsurf;
                        int[] csurfs = cells[icell].csurfs;
                        int a;
                        if (dim == 2)
                        {
                            for (a = 0; a < nsurf; a++)
                            {
                                if (csurfs[a] >= nline_old) break;
                            }
                        }
                        else
                        {
                            for (a = 0; a < nsurf; a++)
                            {
                                if (csurfs[a] >= ntri_old) break;
                            }
                        }
                        if (a < nsurf && partflag == (int)Enum2.CHECK)
                        {
                            if (cinfo[icell].count != 0) delflag = 1;
                            sparta.particle.remove_all_from_cell(cinfo[icell].first);
                            cinfo[icell].count = 0;
                            cinfo[icell].first = -1;
                        }
                    }
                    if (cells[icell].nsplit > 1)
                        sparta.grid.assign_split_cell_particles(icell);
                }
                int nlocal_old = sparta.particle.nlocal;
                if (delflag != 0) sparta.particle.compress_rebalance();
                bigint delta = nlocal_old - sparta.particle.nlocal;
                sparta.mpi.MPI_Allreduce(ref delta, ref ndeleted, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            }

            sparta.mpi.MPI_Barrier(sparta.world);
            double time8 = sparta.mpi.MPI_Wtime();

            double time_total = time6 - time1;

            if (sparta.comm.me == 0)
            {
                if (sparta.particle.exist != 0)
                {
                    string str1 = string.Format("  {0} deleted particles\n", ndeleted);
                    string str2 = string.Format("  CPU time = {0} secs\n", time_total);
                    string str3 = string.Format("  read/sort/check/surf2grid/ghost/inout/particle percent = \n{0:G6}   {1:G6}   {2:G6}   {3:G6}   {4:G6}   {5:G6}   {6:G6}\n",
                            100.0 * (time2 - time1) / time_total, 100.0 * (time3 - time2) / time_total,
                            100.0 * (time4 - time3) / time_total, 100.0 * (time5 - time4) / time_total,
                            100.0 * (time6 - time5) / time_total, 100.0 * (time7 - time6) / time_total,
                            100.0 * (time8 - time7) / time_total);
                    Console.WriteLine(str1 + str2 + str3);
                    if (sparta.logfile != null)
                    {
                        new StreamWriter(sparta.logfile).WriteLine(str1 + str2 + str3);
                    }
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
        protected void check_neighbor_norm_2d()
        {
            int p1, p2;

            // count[I] = # of lines that vertex I is part of

            int[] count;
            int[,] p2e;
            count = new int[npoint_new];
            p2e = new int[npoint_new, 2];
            //memory.create(count, npoint_new, "readsurf:count");
            //memory.create(p2e, npoint_new, 2, "readsurf:count");
            for (int i = 0; i < npoint_new; i++) count[i] = 0;

            int m = nline_old;
            for (int i = 0; i < nline_new; i++)
            {
                p1 = lines[m].p1 - npoint_old;
                p2 = lines[m].p2 - npoint_old;
                p2e[p1,count[p1]++] = m;
                p2e[p2,count[p2]++] = m;
                m++;
            }

            // check that norms of adjacent lines are not in opposite directions

            double dot;
            double[] norm1,norm2;

            int nerror = 0;
            int nwarn = 0;
            for (int i = 0; i < npoint_new; i++)
            {
                if (count[i] == 1) continue;
                norm1 = lines[p2e[i,0]].norm;
                norm2 = lines[p2e[i,1]].norm;
                dot = MathExtra.dot3(norm1, norm2);
                if (dot <= -1.0) nerror++;
                else if (dot < -1.0 + EPSILON_NORM) nwarn++;
            }

            if (nerror!=0)
            {
               string str=string.Format("Surface check failed with {0} infinitely thin line pairs", nerror);
                sparta.error.all(str);
            }
            if (nwarn!=0)
            {
                string str = string.Format("Surface check found {0} nearly infinitely thin line pairs", nwarn);
                if (me == 0) sparta.error.one( str);
            }

        }
        protected void check_neighbor_norm_3d()
        {
            Dictionary<bigint, int> hash = new Dictionary<bigint, int>();
            // insert each edge into hash with triangle as value

            bigint p1, p2, p3, key;

            int m = ntri_old;
            for (int i = 0; i < ntri_new; i++)
            {
                p1 = tris[m].p1;
                p2 = tris[m].p2;
                p3 = tris[m].p3;
                key = (p1 << 32) | p2;
                hash[key] = m;
                key = (p2 << 32) | p3;
                hash[key] = m;
                key = (p3 << 32) | p1;
                hash[key] = m;
                m++;
            }
            // check that norms of adjacent triangles are not in opposite directions

            double dot;
            double[] norm1,norm2;

            int nerror = 0;
            int nwarn = 0;
            foreach (var it in hash)
            {
                p1 = it.Key >> 32;
                p2 = it.Key & Run.MAXSMALLINT;
                key = (p2 << 32) | p1;
                if (hash[key] == hash[hash.Keys.Last()]) continue;
                norm1 = tris[it.Value].norm;
                norm2 = tris[hash[key]].norm;
                dot = MathExtra.dot3(norm1, norm2);
                if (dot <= -1.0) nerror++;
                else if (dot < -1.0 + EPSILON_NORM) nwarn++;
            }
            

            if (nerror!=0)
            {
                string str=string.Format("Surface check failed with {0} infinitely thin triangle pairs", nerror);
                sparta.error.all(str);
            }
            if (nwarn!=0)
            {
                string str = string.Format("Surface check found {0} nearly infinitely thin triangle pairs", nwarn);
                if (me == 0) sparta.error.one(str);
            }
        }
        protected void check_point_near_surf_2d()
        {
            int i, j, n, p1, p2;
            int[] csurfs;
            double side, epssq;
            double[] lo,hi;
            Surf.Line line;

            Grid.ChildCell[] cells = sparta.grid.cells;
            int nglocal = sparta.grid.nlocal;

            int nerror = 0;
            int nwarn = 0;

            for (int icell = 0; icell < nglocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                n = cells[icell].nsurf;
                if (n == 0) continue;

                lo = cells[icell].lo;
                hi = cells[icell].hi;
                side = Math.Min(hi[0] - lo[0], hi[1] - lo[1]);
                epssq = (EPSILON_GRID * side) * (EPSILON_GRID * side);

                csurfs = cells[icell].csurfs;
                for (i = 0; i < n; i++)
                {
                    line = lines[csurfs[i]];
                    for (j = 0; j < n; j++)
                    {
                        if (i == j) continue;
                        p1 = lines[csurfs[j]].p1;
                        p2 = lines[csurfs[j]].p2;
                        point_line_compare(p1, line, epssq,ref nerror, ref nwarn);
                        point_line_compare(p2, line, epssq, ref nerror, ref nwarn);
                    }
                }
            }

            int all=0;
            sparta.mpi.MPI_Allreduce(ref nerror, ref all, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            if (all!=0)
            {
                string str=string.Format("Surface check failed with %d points on lines", all);
                sparta.error.all(str);
            }

            sparta.mpi.MPI_Allreduce(ref nwarn, ref all, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            if (all!=0)
            {
                string str = string.Format("Surface check found %d points nearly on lines", all);
                if (me == 0) sparta.error.all( str);
            }
        }
        protected void check_point_near_surf_3d()
        {
            Console.WriteLine("readsurf.check_point_near_surf_3d()"); 
        }

        protected void point_line_compare(int i, Surf.Line line, double epssq,
                                  ref int nerror,ref int nwarn)
        {
            if (i == line.p1 || i == line.p2) return;
            double rsq =
              Geometry.distsq_point_line(pts[i].x, pts[line.p1].x, pts[line.p2].x);
            if (rsq == 0.0) nerror++;
            else if (rsq < epssq) nwarn++;
        }
        //protected void point_tri_compare(int, Surf.Tri*, double, int &, int &, int, int, int);

        //protected int find_edge(int, int);
        //protected void add_edge(int, int, int);

        protected double shortest_line()
        {
            double len = BIG;
            int m = nline_old;
            for (int i = 0; i < nline_new; i++)
            {
                len = Math.Min(len, sparta.surf.line_size(m));
                m++;
            }
            return len;
        }
        protected void smallest_tri(out double len,out double area)
        {
            double lenone, areaone;

            len = area = BIG;
            int m = ntri_old;
            for (int i = 0; i < ntri_new; i++)
            {
                areaone = sparta.surf.tri_size(m,out lenone);
                len =  Math.Min(len, lenone);
                area = Math.Min(area, areaone);
                m++;
            }
        }
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
