using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class Cut2d
    {
        enum Enum1{ UNKNOWN, OUTSIDE, INSIDE, OVERLAP };     // several files
        enum Enum2{ EXTERIOR, INTERIOR, BORDER, INTBORD };
        enum Enum3{ ENTRY, EXIT, TWO, CORNER };              // same as Cut3d

        public const double EPSCELL = 1.0e-10;    // tolerance for pushing surf pts to cell surface

        // cell ID for 2d or 3d cell

        //#define VERBOSE
        public const int VERBOSE_ID = 37;
        //#define VERBOSE_ID 27810406321L

        public int npushmax;          // # of push options to try
        public int[] npushcell=new int[4];      // tally of cells that required surf point push

        public struct Cline
        {
            public double[] x, y;   // coords of end points of line clipped to cell
            public int line;           // index in list of lines that intersect this cell
        };

        public struct Point
        {
            public double[] x;        // coords of point
            public int type;           // type of pt = ENTRY,EXIT,TWO,CORNER
                                       // ENTRY/EXIT = only part of one Cline,
                                       //   could also be a geometric corner pt
                                       // TWO = part of two Clines
                                       // CORNER = part of no Cline, is a geometric corner pt
            public int next;           // index of next point when walking a loop
                                       // set for ENTRY/TWO pts between ENTRY and EXIT
                                       // set for EXIT/CORNER points around cell perimeter,
                                       //   though may not be walked
            public int line;           // original line (as stored by Cline) the pt starts,
                                       //   only set for ENTRY and TWO pts
            public int corner;         // 1,2,3,4 if x is a corner point, else 0
                                       // could be ENTRY,EXIT,CORNER pt, but not a TWO pt
            public int cprev, cnext;    // indices of pts in linked list around cell perimeter
            public int side;           // which side of cell (0,1,2,3) pt is on
                                       // only for ENTRY/EXIT/CORNER pts to make linked list
            public double value;       // coord along the side
                                // only for ENTRY/EXIT/CORNER pts to make linked list
        };

        public struct Loop
        {
            public double area;        // area of loop
            public int active;         // 1/0 if active or not
            public int flag;           // INTERIOR (if all TWO points) or BORDER
            public int n;              // # of points in loop
            public int first;          // index of first point in loop
            public int next;           // index of next loop in same PG, -1 if last loop
        };

        public struct PG
        {
            public double area;        // summed area (over loops) of PG
            public int n;              // # of loops in PG
            public int first;          // index of first loop in PG
        };

        public List<Cline> clines;  // list of Clines
        public List<Point> points;  // list of Points = Weiler/Atherton data structure
        public List<Loop> loops;    // list of loops in Points
        public List<PG> pgs;        // list of polygons = one or more loops
        private SPARTA sparta;
        
        public Cut2d(SPARTA sparta, int caller_axisymmetric)
        {
            this.sparta = sparta;

            axisymmetric = caller_axisymmetric;

            npushmax = 2;    // if increase this, increase push vec size in cut2d.h

            pushlo_vec[0] = -1.0;
            pushhi_vec[0] = 1.0;
            pushvalue_vec[0] = 0.0;
            pushlo_vec[1] = -1.0;
            pushhi_vec[1] = 1.0;
            pushvalue_vec[1] = 1.0;

            if (sparta.surf.pushflag == 0) npushmax = 0;
            if (sparta.surf.pushflag == 2) npushmax = 3;
            if (sparta.surf.pushflag == 2)
            {
                pushlo_vec[2] = sparta.surf.pushlo;
                pushhi_vec[2] = sparta.surf.pushhi;
                pushvalue_vec[2] = sparta.surf.pushvalue;
            }

            for (int i = 0; i <= npushmax; i++) npushcell[i] = 0;
        }
        public int surf2grid(int id_caller, double[] lo_caller, double[] hi_caller,
                     int[] surfs_caller, int max)
        {
            id = id_caller;
            lo = lo_caller;
            hi = hi_caller;
            surfs = surfs_caller;

            List<Surf.Point> pts = sparta.surf.pts;
            List<Surf.Line> lines = sparta.surf.lines;
            int nline = sparta.surf.nline;

            double[] x1,x2;

            nsurf = 0;
            for (int m = 0; m < nline; m++)
            {
                x1 = pts[lines[m].p1].x;
                x2 = pts[lines[m].p2].x;

                if (Math.Max(x1[0], x2[0]) < lo[0]) continue;
                if (Math.Min(x1[0], x2[0]) > hi[0]) continue;
                if (Math.Max(x1[1], x2[1]) < lo[1]) continue;
                if (Math.Min(x1[1], x2[1]) > hi[1]) continue;

                if (cliptest(x1, x2)!=0)
                {
                    if (nsurf == max) return -1;
                    surfs[nsurf++] = m;
                }
            }

            return nsurf;
        }
        //public int surf2grid_list(cellint, double[], double[], int, int[], int[], int);
        public int split(int id_caller, double[] lo_caller, double[] hi_caller,
                 int nsurf_caller, int[] surfs_caller,
                 out double[] areas_caller, int[] surfmap,
                 int[] corners,out int xsub, double[] xsplit)
        {
            id = id_caller;
            lo = lo_caller;
            hi = hi_caller;
            nsurf = nsurf_caller;
            surfs = surfs_caller;

            // perform cut/split
            // first attempt is with no pushing of surface points
            // if fails, then try again with push options
            // if all push options fail, then print error message

            int nsplit=0, errflag;
            pushflag = 0;

            while (true)
            {
                int grazeflag = build_clines();

                // all triangles just touched cell surface
                // mark corner points based on grazeflag and in/out line orientation
                // return area = 0.0 for UNKNOWN/INSIDE, full cell area for OUTSIDE
                // vol is changed in Grid::set_inout() if OVERLAP cell corners are marked

                if (clines.Count == 0)
                {
                    if (pushflag!=0) npushcell[pushflag]++;

                    int mark =(int)Enum1.UNKNOWN;
                    if (grazeflag!=0 || inout == (int)Enum1.INSIDE) mark = (int)Enum1.INSIDE;
                    else if (inout == (int)Enum1.OUTSIDE) mark = (int)Enum1.OUTSIDE;
                    corners[0] = corners[1] = corners[2] = corners[3] = mark;

                    double area = 0.0;
                    if (mark == (int)Enum1.OUTSIDE)
                    {
                        if (axisymmetric!=0)
                            area =MyConst.MY_PI * (hi[1] * hi[1] - lo[1] * lo[1]) * (hi[0] - lo[0]);
                        else area = (hi[0] - lo[0]) * (hi[1] - lo[1]);
                    }

                    areas=new List<double>(1);
                    areas[0] = area;
                    areas_caller = areas.ToArray();
                    xsub = 0;
                    return 1;
                }

                // 3 operations can generate errors: weiler_build, loop2pg, split_point
                // value of errflag corresponds to unique error message (listed below)

                errflag = weiler_build();
                if (errflag!=0)
                {
                    if (push_increment()!=0) continue;
                    break;
                }

                weiler_loops();
                errflag = loop2pg();
                if (errflag != 0)
                {
                    if (push_increment() != 0) continue;
                    break;
                }

                nsplit = pgs.Count;
                if (nsplit > 1)
                {
                    create_surfmap(surfmap);
                    xsub = 0;
                    errflag = split_point(surfmap, xsplit,ref xsub);
                }
                if (errflag != 0)
                {
                    if (push_increment() != 0) continue;
                    break;
                }

                // successful cut/split
                // set corners = OUTSIDE if corner pt is in list of points in PGs
                // else set corners = INSIDE

                corners[0] = corners[1] = corners[2] = corners[3] = (int)Enum1.INSIDE;

                int iloop, nloop, mloop, ipt, npt, mpt;

                int npg = pgs.Count;
                for (int ipg = 0; ipg < npg; ipg++)
                {
                    nloop = pgs[ipg].n;
                    mloop = pgs[ipg].first;
                    for (iloop = 0; iloop < nloop; iloop++)
                    {
                        npt = loops[mloop].n;
                        mpt = loops[mloop].first;
                        for (ipt = 0; ipt < npt; ipt++)
                        {
                            if (points[mpt].corner >= 0)
                                corners[points[mpt].corner] = (int)Enum1.OUTSIDE;
                            mpt = points[mpt].next;
                        }
                        mloop = loops[mloop].next;
                    }
                }

                // store areas in vector so can return ptr to it

                areas=new List<double>(nsplit);
                for (int i = 0; i < nsplit; i++) areas[i] = pgs[i].area;
                areas_caller = areas.ToArray();

                break;
            }

            // could not perform cut/split -> fatal error
            // print info about cell and final error message
            // 2-letter prefix is which method encountered error
            // NOTE: store errflag_original for no-push error to print this message?

            if (errflag != 0)
            {
                failed_cell();

                if (errflag == 1)
                    sparta.error.one("WB: Point appears first in more than one CLINE");
                if (errflag == 2)
                    sparta.error.one("WB: Point appears last in more than one CLINE");
                if (errflag == 3)
                    sparta.error.one("WB: Singlet CLINES point not on cell border");
                if (errflag == 4)
                    sparta.error.one("LP: No positive areas in cell");
                if (errflag == 5)
                    sparta.error.one("LP: More than one positive area with a negative area");
                if (errflag == 6)
                    sparta.error.one("LP: Single area is negative, inverse donut");
                if (errflag == 7)
                    sparta.error.one("SP: Could not find split point in split cell");
            }

            if (pushflag != 0) npushcell[pushflag]++;
            areas_caller = null;
            xsub = 0;
            return nsplit;

        }
        //public int split_face(int, int, double[], double[]);
        //public int clip_external(double[], double[], double[], double[], double[]);


        private int axisymmetric;
        private int id;            // ID of cell being worked on
        private double[] lo,hi;        // opposite corner pts of cell
        private int nsurf;             // # of surf elements in cell
        private int[] surfs;            // indices of surf elements in cell
         
        private int pushflag;          // 0 for no push, else push surf points near cell surf
        private double pushlo, pushhi;  // lo/hi ranges to push on
        private double pushvalue;      // new position to push to
        private double[] pushlo_vec=new double[3], pushhi_vec = new double[3], pushvalue_vec = new double[3];  // push values to try
        private int inout;             // orientation of lines that just touch cell
         
        private List<double> areas;   // areas of each flow polygon found
        private List<int> used;       // 0/1 flag for each point when walking loops

        private int build_clines()
        {
            int m;
            double[] p1 = new double[2], p2 = new double[2], cbox = new double[3],
                cmid = new double[3], l2b = new double[3];
            double[] x, y, norm, pp1, pp2;
            Surf.Line line;
            Cline cline;

            List<Surf.Point> pts = sparta.surf.pts;
            List<Surf.Line> lines = sparta.surf.lines;

            clines = new List<Cline>(nsurf);


            int grazeflag = 0;
            int noutside = 0;
            int ninside = 0;
            int n = 0;

            for (int i = 0; i < nsurf; i++)
            {
                m = surfs[i];
                line = lines[m];

                p1 = pts[line.p1].x;
                p2 = pts[line.p2].x;

                // if pushflag is set, push line pts near cell surface

                if (pushflag != 0)
                {
                    push(p1);
                    push(p2);
                }

                cline = clines[n];
                cline.line = i;

                // clip PQ to cell and store as XY in Cline

                x = cline.x;
                y = cline.y;
                clip(p1, p2, x, y);

                // discard clipped line if only one point
                // if all lines are removed, clines will be empty,
                //   which will result in cell corner pts being left UNKNOWN in split()
                // to try and avoid this, tally the in/out orientation of all removed lines
                // "out" orientation means line norm from line ctr points towards cell ctr 
                // "in" orientation means line norm from line ctr points away from cell ctr
                // some lines may not follow this rule, but majority should
                // cbox = cell center pt, cmid = line center pt, l2b = cbox-cmid

                if (x[0] == y[0] && x[1] == y[1])
                {
                    cbox[0] = 0.5 * (lo[0] + hi[0]);
                    cbox[1] = 0.5 * (lo[1] + hi[1]);
                    cbox[2] = 0.0;
                    pp1 = sparta.surf.pts[line.p1].x;
                    pp2 = sparta.surf.pts[line.p2].x;
                    cmid[0] = 0.5 * (pp1[0] + pp2[0]);
                    cmid[1] = 0.5 * (pp1[1] + pp2[1]);
                    cmid[2] = 0.0;
                    MathExtra.sub3(cbox, cmid, l2b);
                    double dot = MathExtra.dot3(line.norm, l2b);
                    if (dot > 0.0) noutside++;
                    if (dot < 0.0) ninside++;

                    continue;
                }

                // discard clipped line if lies along one cell edge
                //   and normal is not into cell
                // leave grazeflag incremented in this case

                grazeflag++;
                if (ptflag(x) == (int)Enum2.BORDER && ptflag(y) == (int)Enum2.BORDER)
                {
                    int edge = sameedge(x, y);
                    if (edge != 0)
                    {
                        norm = line.norm;
                        if (edge == 1 && norm[0] < 0.0) continue;
                        if (edge == 2 && norm[0] > 0.0) continue;
                        if (edge == 3 && norm[1] < 0.0) continue;
                        if (edge == 4 && norm[1] > 0.0) continue;
                        grazeflag--;
                    }
                }

                n++;
            }

            // if empty, also return in/out status based on deleted line orientations
            // only mark corner pts for INSIDE/OUTSIDE if all deleted lines
            //   had same orientation, else leave them UNKNOWN
            // NOTE: marking if majority are INSIDE/OUTSIDE (commented out)
            //   triggered later marking error for cone_test/in.cone problem

            inout = (int)Enum1.UNKNOWN;
            if (n == 0)
            {
                if (ninside!=0 && noutside == 0) inout = (int)Enum1.INSIDE;
                else if (noutside!=0 && ninside == 0) inout = (int)Enum1.OUTSIDE;
                //if (ninside > noutside) inout = INSIDE;
                //else if (noutside > ninside) inout = OUTSIDE;
            }

            return grazeflag;
        }
        private int weiler_build()
        {
            int i, j;
            int firstpt, lastpt, nextpt;
            double[] pt;

            // insure space in points for 2 pts per cline + 4 corner pts

            int nlines = clines.Count;
            points=new List<Point>(2 * nlines + 4);

            // add each cline end pt to points
            // set x,type,next,line for each pt
            // NOTE: could use hash to find existing pts in O(1) time
            //   see Brian Adams code (howto/c++) for doing this on a vec of doubles

            int npt = 0;

            for (i = 0; i < nlines; i++)
            {

                // 1st point in cline

                pt = clines[i].x;
                for (j = 0; j < npt; j++)
                    if (pt[0] == points[j].x[0] && pt[1] == points[j].x[1]) break;

                // pt already exists

                if (j < npt)
                {
                    Point point = points[j];
                    if (point.type == (int)Enum3.ENTRY || point.type == (int)Enum3.TWO) return 1;

                    point.type = (int)Enum3.TWO;
                    point.line = clines[i].line;
                    points[j] = point;
                    firstpt = j;

                }

                // new pt

                else
                {
                    Point point = points[npt];
                    point.x[0] = pt[0];
                    point.x[1] = pt[1];
                    point.type = (int)Enum3.ENTRY;
                    point.line = clines[i].line;
                    points[npt] = point;
                    firstpt = npt;
                    npt++;
                }

                // 2nd point in cline

                pt = clines[i].y;
                for (j = 0; j < npt; j++)
                    if (pt[0] == points[j].x[0] && pt[1] == points[j].x[1]) break;

                // pt already exists

                if (j < npt)
                {
                    Point point = points[j];
                    if (point.type == (int)Enum3.EXIT || point.type == (int)Enum3.TWO) return 2;
                    point.type = (int)Enum3.TWO;
                    points[j] = point;
                    point = points[firstpt];
                    point.next = j;
                    points[firstpt] = point;
                }

                // new pt

                else
                {
                    Point point = points[npt];
                    point.x[0] = pt[0];
                    point.x[1] = pt[1];
                    point.type = (int)Enum3.EXIT;
                    points[npt] = point;
                    point = points[firstpt];
                    point.next = npt;
                    points[firstpt] = point;
                    npt++;
                }
            }

            // error check that every singlet point is on cell border

            for (i = 0; i < npt; i++)
                if (points[i].type != (int)Enum3.TWO && ptflag(points[i].x) != (int)Enum2.BORDER) return 3;

            // add 4 cell CORNER pts to points
            // only if corner pt is not already an ENTRY or EXIT pt
            // if a TWO pt, still add corner pt as CORNER
            // each corner pt is at beginning of side it is on
            // NOTE: could use hash to find existing pts in O(1) time
            // corner flag = 0,1,2,3 for LL,LR,UL,UR, same as in Grid::ChildInfo,
            //   but ordering of corner pts in linked list is LL,LR,UR,UL
            // side = 0,1,2,3 for lower,right,upper,left = traversal order in linked list
            // value = x-coord for lower/upper sides, y-coord for left,right sides

            double[] cpt=new double[2];
            int ipt1, ipt2, ipt3, ipt4;

            for (i = 0; i < npt; i++)
            {
                Point point = points[i];
                point.corner = -1;
                points[i] = point;
            }

            cpt[0] = lo[0]; cpt[1] = lo[1];
            for (j = 0; j < npt; j++)
                if (cpt[0] == points[j].x[0] && cpt[1] == points[j].x[1]) break;
            if (j == npt || points[j].type ==(int)Enum3.TWO)
            {
                Point point = points[npt];

                point.x[0] = cpt[0];
                point.x[1] = cpt[1];
                point.type = (int)Enum3.CORNER;
                points[npt]=point;
                ipt1 = npt++;
            }
            else ipt1 = j;
            Point apoint = points[ipt1];
            apoint.corner = 0;
            apoint.side = 0;
            apoint.value = lo[0];
            points[ipt1] = apoint;

            cpt[0] = hi[0]; cpt[1] = lo[1];
            for (j = 0; j < npt; j++)
                if (cpt[0] == points[j].x[0] && cpt[1] == points[j].x[1]) break;
            if (j == npt || points[j].type == (int)Enum3.TWO)
            {
                Point point = points[npt];
                point.x[0] = cpt[0];
                point.x[1] = cpt[1];
                point.type = (int)Enum3.CORNER;
                points[npt] = point;
                ipt2 = npt++;
            }
            else ipt2 = j;
            apoint = points[ipt2];
            apoint.corner = 1;
            apoint.side = 1;
            apoint.value = lo[1];
            points[ipt2] = apoint;

            cpt[0] = hi[0]; cpt[1] = hi[1];
            for (j = 0; j < npt; j++)
                if (cpt[0] == points[j].x[0] && cpt[1] == points[j].x[1]) break;
            if (j == npt || points[j].type == (int)Enum3.TWO)
            {
                Point point = points[npt];
                point.x[0] = cpt[0];
                point.x[1] = cpt[1];
                point.type = (int)Enum3.CORNER;
                points[npt] = point;
                ipt3 = npt++;
            }
            else ipt3 = j;
            apoint = points[ipt3];

            apoint.corner = 3;
            apoint.side = 2;
            apoint.value = hi[0];
            points[ipt3]=apoint;

            cpt[0] = lo[0]; cpt[1] = hi[1];
            for (j = 0; j < npt; j++)
                if (cpt[0] == points[j].x[0] && cpt[1] == points[j].x[1]) break;
            if (j == npt || points[j].type == (int)Enum3.TWO)
            {
                Point point = points[npt];
                point.x[0] = cpt[0];
                point.x[1] = cpt[1];
                point.type = (int)Enum3.CORNER;
                points[npt] = point;
                ipt4 = npt++;
            }
            else ipt4 = j;
            apoint = points[ipt4];
            apoint.corner = 2;
            apoint.side = 3;
            apoint.value = hi[1];
            points[ipt3]=apoint;

            // n = final # of points


            // create counter-clockwise linked list of 4 pts around cell perimeter

            firstpt = ipt1;
            lastpt = ipt4;
            Point p1= points[ipt1], p2= points[ipt2], p3= points[ipt3],p4= points[ipt4];
            p1.cprev = -1;
            p1.cnext = ipt2;
            p2.cprev = ipt1;
            p2.cnext = ipt3;
            p3.cprev = ipt2;
            p3.cnext = ipt4;
            p3.cprev = ipt3;
            p3.cnext = -1;

            points[ipt1] = p1;
            points[ipt2] = p2;
            points[ipt3] = p3;
            points[ipt4] = p4;


            // add all non-corner ENTRY/EXIT pts to counter-clockwise linked list
            // side = 0,1,2,3 for lower,right,upper,left sides
            // value = coord of pt along the side it is on

            int ipt, iprev=0, side;
            double value;

            for (i = 0; i < npt; i++)
            {
                if (points[i].type == (int)Enum3.TWO || points[i].type == (int)Enum3.CORNER) continue;
                if (points[i].corner >= 0) continue;

                side = whichside(points[i].x);
                if (side % 2!=0) value = points[i].x[1];
                else value = points[i].x[0];

                // interleave Ith point into linked list between firstpt and lastpt
                // insertion location is between iprev and ipt
                // special logic if inserting at end of list

                ipt = firstpt;
                while (ipt >= 0)
                {
                    if (side < points[ipt].side) break;
                    if (side == points[ipt].side)
                    {
                        if (side < 2)
                        {
                            if (value < points[ipt].value) break;
                        }
                        else
                        {
                            if (value > points[ipt].value) break;
                        }
                    }
                    iprev = ipt;
                    ipt = points[ipt].cnext;
                }
                apoint = points[i];
                apoint.side = side;
                apoint.value = value;
                apoint.cprev = iprev;
                apoint.cnext = ipt;
                points[i] = apoint;

                apoint = points[iprev];
                apoint.cnext = i;
                points[iprev] = apoint;

                if (ipt >= 0)
                {
                    Point point = points[ipt];
                    point.cprev = i;
                    points[ipt] = point;
                }
                else lastpt = i;
            }

            // set next field for cell perimeter points in linked list
            // this completes loops for next fields already set between ENTRY/EXIT pts
            // do not reset next for ENTRY pts
            // after loop, explicitly connect lastpt to firstpt

            ipt = firstpt;
            while (ipt >= 0)
            {
                nextpt = points[ipt].cnext;
                if (points[ipt].type != (int)Enum3.ENTRY)
                {
                    Point point = points[ipt];
                    point.next = nextpt;
                    points[ipt] = point;
                }

                ipt = nextpt;
            }
            if (points[lastpt].type != (int)Enum3.ENTRY)
            {
                Point point = points[lastpt];
                point.next = firstpt;
                points[lastpt] = point;
            }

            // successful return

            return 0;
        }
        private void weiler_loops()
        {
            // used = 0/1 flag for whether a point is already part of a loop

            int n = points.Count;
            used=new List<int>(n);
            for (int i = 0; i < n; i++) used[i] = 0;

            // iterate over all pts
            // start a loop at any unused pt
            // walk loop via next field:
            //   mark pts as used
            //   compute area along the way
            // stop when reach initial pt:
            //   closed loop, add to loops data structure
            //   loop of just 4 corner pts is a valid loop
            // stop when reach used pt:
            //   discard loop, just traversed non-loop corner pts

            int ipt, iflag, cflag, ncount, firstpt, nextpt;
            double area;
            double[] x,y;

            int nloop = 0;

            for (int i = 0; i < n; i++)
            {
                if (used[i]!=0) continue;
                area = 0.0;
                iflag = cflag = 1;
                ncount = 0;

                ipt = firstpt = i;
                x = points[ipt].x;

                while (used[ipt]==0)
                {
                    used[ipt] = 1;
                    ncount++;
                    if (points[ipt].type != (int)Enum3.TWO) iflag = 0;
                    if (points[ipt].type != (int)Enum3.CORNER) cflag = 0;
                    nextpt = points[ipt].next;
                    y = points[nextpt].x;
                    if (axisymmetric!=0)
                        area -= MyConst.MY_PI3 * (x[1] * x[1] + x[1] * y[1] + y[1] * y[1]) * (y[0] - x[0]);
                    else area -= (0.5 * (x[1] + y[1]) - lo[1]) * (y[0] - x[0]);
                    x = y;
                    ipt = nextpt;
                    if (ipt == firstpt) break;
                }
                if (ipt != firstpt) continue;

                loops=new List<Loop>(nloop + 1);
                Loop loop = loops[nloop];
                loop.area = area;
                loop.active = 1;
                if (iflag!=0) loop.flag = (int)Enum2.INTERIOR;
                else if (cflag!=0) loop.flag = (int)Enum2.BORDER;
                else loop.flag = (int)Enum2.INTBORD;
                loop.n = ncount;
                loop.first = firstpt;
                loops[nloop] = loop;

                nloop++;
            }

        }
        private int loop2pg()
        {
            int positive = 0;
            int negative = 0;

            int nloop = loops.Count;
            for (int i = 0; i < nloop; i++)
                if (loops[i].area > 0.0) positive++;
                else negative++;
            if (positive == 0) return 4;
            if (positive > 1 && negative!=0) return 5;

            pgs=new List<PG>(positive);

            // if multiple positive loops, mark BORDER loop as inactive if exists
            // don't want entire cell border to be a loop in this case

            if (positive > 1)
            {
                for (int i = 0; i < nloop; i++)
                    if (loops[i].flag == (int)Enum2.BORDER)
                    {
                        Loop loop = loops[i];
                        loop.active = 0;
                        loops[i] = loop;
                        positive--;
                    }
            }

            // positive = 1 means 1 PG with area = sum of all pos/neg loops
            // positive > 1 means each loop is a PG

            if (positive == 1)
            {
                double area = 0.0;
                int prev = -1;
                int count = 0;
                int first=0;

                for (int i = 0; i < nloop; i++)
                {
                    if (loops[i].active==0) continue;
                    area += loops[i].area;
                    count++;
                    if (prev < 0) first = i;
                    else
                    {
                        Loop loop = loops[prev];
                        loop.next = i;
                        loops[prev] = loop;
                    }

                    prev = i;
                }
                Loop aloop = loops[prev];
                aloop.next = -1;
                loops[prev] = aloop;

                if (area < 0.0) return 6;
                PG pg = pgs[0];
                pg.area = area;
                pg.n = count;
                pg.first = first;
                pgs[0] = pg;

            }
            else
            {
                int m = 0;
                for (int i = 0; i < nloop; i++)
                {
                    if (loops[i].active == 0) continue;
                    PG pg = pgs[m];
                    pg.area = loops[i].area;
                    pg.n = 1;
                    pg.first = i;
                    pgs[m] = pg;
                    m++;
                    Loop loop = loops[i];
                    loop.next = -1;
                }
            }

            return 0;
        }
        private void create_surfmap(int[] surfmap)
        {
            for (int i = 0; i < nsurf; i++) surfmap[i] = -1;

            int iloop, nloop, mloop, ipt, npt, mpt;

            int npg = pgs.Count;
            for (int ipg = 0; ipg < npg; ipg++)
            {
                nloop = pgs[ipg].n;
                mloop = pgs[ipg].first;
                for (iloop = 0; iloop < nloop; iloop++)
                {
                    npt = loops[mloop].n;
                    mpt = loops[mloop].first;
                    for (ipt = 0; ipt < npt; ipt++)
                    {
                        if (points[mpt].type == (int)Enum3.TWO || points[mpt].type == (int)Enum3.ENTRY)
                            surfmap[points[mpt].line] = ipg;
                        mpt = points[mpt].next;
                    }
                    mloop = loops[mloop].next;
                }
            }
        }
        private int split_point(int[] surfmap, double[] xsplit, ref int xsub)
        {
            int iline;
            double[] x1,x2;
            double[] a=new double[2], b = new double[2];

            List<Surf.Point> pts = sparta.surf.pts;
            List < Surf.Line> lines = sparta.surf.lines;

            // if end pt of any line with non-negative surfmap is in/on cell, return

            for (int i = 0; i < nsurf; i++)
            {
                if (surfmap[i] < 0) continue;
                iline = surfs[i];
                x1 = pts[lines[iline].p1].x;
                x2 = pts[lines[iline].p2].x;
                if (ptflag(x1) != (int)Enum2.EXTERIOR)
                {
                    xsplit[0] = x1[0]; xsplit[1] = x1[1];
                    xsub = surfmap[i];
                    return 0;
                }
                if (ptflag(x2) != (int)Enum2.EXTERIOR)
                {
                    xsplit[0] = x2[0]; xsplit[1] = x2[1];
                    xsub = surfmap[i];
                    return 0;
                }
            }

            // clip 1st line with non-negative surfmap to cell, and return clip point

            for (int i = 0; i < nsurf; i++)
            {
                if (surfmap[i] < 0) continue;
                iline = surfs[i];
                x1 = pts[lines[iline].p1].x;
                x2 = pts[lines[iline].p2].x;
                clip(x1, x2, a, b);
                xsplit[0] = a[0]; xsplit[1] = a[1];
                xsub = surfmap[i];
                return 0;
            }

            // error return
            return 7;
        }

        private int cliptest(double[] p, double[] q)
        {
            double x, y;

            if (p[0] >= lo[0] && p[0] <= hi[0] &&
                p[1] >= lo[1] && p[1] <= hi[1]) return 1;
            if (q[0] >= lo[0] && q[0] <= hi[0] &&
                q[1] >= lo[1] && q[1] <= hi[1]) return 1;

            double[] a=new double[2], b = new double[2];
            a[0] = p[0]; a[1] = p[1];
            b[0] = q[0]; b[1] = q[1];

            if (a[0] < lo[0] && b[0] < lo[0]) return 0;
            if (a[0] < lo[0] || b[0] < lo[0])
            {
                y = a[1] + (lo[0] - a[0]) / (b[0] - a[0]) * (b[1] - a[1]);
                if (a[0] < lo[0])
                {
                    a[0] = lo[0]; a[1] = y;
                }
                else
                {
                    b[0] = lo[0]; b[1] = y;
                }
            }
            if (a[0] > hi[0] && b[0] > hi[0]) return 0;
            if (a[0] > hi[0] || b[0] > hi[0])
            {
                y = a[1] + (hi[0] - a[0]) / (b[0] - a[0]) * (b[1] - a[1]);
                if (a[0] > hi[0])
                {
                    a[0] = hi[0]; a[1] = y;
                }
                else
                {
                    b[0] = hi[0]; b[1] = y;
                }
            }

            if (a[1] < lo[1] && b[1] < lo[1]) return 0;
            if (a[1] < lo[1] || b[1] < lo[1])
            {
                x = a[0] + (lo[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
                if (a[1] < lo[1])
                {
                    a[0] = x; a[1] = lo[1];
                }
                else
                {
                    b[0] = x; b[1] = lo[1];
                }
            }
            if (a[1] > hi[1] && b[1] > hi[1]) return 0;
            if (a[1] > hi[1] || b[1] > hi[1])
            {
                x = a[0] + (hi[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
                if (a[1] > hi[1])
                {
                    a[0] = x; a[1] = hi[1];
                }
                else
                {
                    b[0] = x; b[1] = hi[1];
                }
            }

            return 1;
        }
        private void clip(double[] p, double[] q, double[] a, double[] b)
        {
            double x, y;

            a[0] = p[0]; a[1] = p[1];
            b[0] = q[0]; b[1] = q[1];

            if (p[0] >= lo[0] && p[0] <= hi[0] &&
                p[1] >= lo[1] && p[1] <= hi[1] &&
                q[0] >= lo[0] && q[0] <= hi[0] &&
                q[1] >= lo[1] && q[1] <= hi[1]) return;

            if (a[0] < lo[0] || b[0] < lo[0])
            {
                y = a[1] + (lo[0] - a[0]) / (b[0] - a[0]) * (b[1] - a[1]);
                if (a[0] < lo[0])
                {
                    a[0] = lo[0]; a[1] = y;
                }
                else
                {
                    b[0] = lo[0]; b[1] = y;
                }
            }
            if (a[0] > hi[0] || b[0] > hi[0])
            {
                y = a[1] + (hi[0] - a[0]) / (b[0] - a[0]) * (b[1] - a[1]);
                if (a[0] > hi[0])
                {
                    a[0] = hi[0]; a[1] = y;
                }
                else
                {
                    b[0] = hi[0]; b[1] = y;
                }
            }
            if (a[1] < lo[1] || b[1] < lo[1])
            {
                x = a[0] + (lo[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
                if (a[1] < lo[1])
                {
                    a[0] = x; a[1] = lo[1];
                }
                else
                {
                    b[0] = x; b[1] = lo[1];
                }
            }
            if (a[1] > hi[1] || b[1] > hi[1])
            {
                x = a[0] + (hi[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
                if (a[1] > hi[1])
                {
                    a[0] = x; a[1] = hi[1];
                }
                else
                {
                    b[0] = x; b[1] = hi[1];
                }
            }
        }

        private int ptflag(double[] pt)
        {
            double x = pt[0];
            double y = pt[1];
            if (x < lo[0] || x > hi[0] || y < lo[1] || y > hi[1]) return (int)Enum2.EXTERIOR;
            if (x > lo[0] && x < hi[0] && y > lo[1] && y < hi[1]) return (int)Enum2.INTERIOR;
            return (int)Enum2.BORDER;
        }
        private int push_increment()
        {
            if (pushflag == npushmax) return 0;
            pushlo = pushlo_vec[pushflag];
            pushhi = pushhi_vec[pushflag]; ;
            pushvalue = pushvalue_vec[pushflag];
            pushflag++;
            return 1;
        }
        private void push(double[] pt)
        {
            double x = pt[0];
            double y = pt[1];
            double epsx = EPSCELL * (hi[0] - lo[0]);
            double epsy = EPSCELL * (hi[1] - lo[1]);

            //if (x < lo[0]-pushhi*epsx || x > hi[0]+pushhi*epsx || 
            //    y < lo[1]-pushhi*epsy || y > hi[1]+pushhi*epsy) return;

            if (x > lo[0] - pushhi * epsx && x < lo[0] - pushlo * epsx) x = lo[0] - pushvalue * epsx;
            if (x > hi[0] + pushlo * epsx && x < hi[0] + pushhi * epsx) x = hi[0] + pushvalue * epsx;
            if (y > lo[1] - pushhi * epsy && y < lo[1] - pushlo * epsy) y = lo[1] - pushvalue * epsy;
            if (y > hi[1] + pushlo * epsy && y < hi[1] + pushhi * epsy) y = hi[1] + pushvalue * epsy;

            double[] boxlo = sparta.domain.boxlo;
            double[] boxhi = sparta.domain.boxhi;
            x = Math.Max(x, boxlo[0]);
            x = Math.Min(x, boxhi[0]);
            y = Math.Max(y, boxlo[1]);
            y = Math.Min(y, boxhi[1]);

            if (x != pt[0] || y != pt[1])
            {
                pt[0] = x;
                pt[1] = y;
            }
        }
        private int sameedge(double[] a, double[] b)
        {
            if (a[0] == lo[0] && b[0] == lo[0]) return 1;
            if (a[0] == hi[0] && b[0] == hi[0]) return 2;
            if (a[1] == lo[1] && b[1] == lo[1]) return 3;
            if (a[1] == hi[1] && b[1] == hi[1]) return 4;
            return 0;
        }
        private int whichside(double[] pt)
        {
            if (pt[0] == lo[0]) return 3;
            if (pt[0] == hi[0]) return 1;
            if (pt[1] == lo[1]) return 0;
            if (pt[1] == hi[1]) return 2;
            return -1;
        }

        private void failed_cell()
        {
            Console.WriteLine("Cut2d failed in cell ID: {0} \n", id);

            int ichild;
            int iparent = sparta.grid.id_find_parent(id, out ichild);
            while (iparent >= 0)
            {
                int nx = sparta.grid.pcells[iparent].nx;
                int ny = sparta.grid.pcells[iparent].ny;
                int ix = (ichild - 1) % nx;
                int iy = ((ichild - 1) / nx) % ny;
                int iz = (ichild - 1) / (nx * ny);
                Console.WriteLine("  parent {0:G} level {1}: NxNyNz {2} {3} {4}: child {5:G} {6} {7} {8} \n",
                       sparta.grid.pcells[iparent].id,
                       sparta.grid.pcells[iparent].level,
                       sparta.grid.pcells[iparent].nx,
                       sparta.grid.pcells[iparent].ny,
                       sparta.grid.pcells[iparent].nz,
                       ichild, ix, iy, iz);
                if (iparent == 0) break;
                iparent = sparta.grid.id_find_parent(sparta.grid.pcells[iparent].id,out ichild);
            }

            Console.WriteLine("  lo corner {0:G} {1:G}\n", lo[0], lo[1]);
            Console.WriteLine("  hi corner {0:G} {1:G}\n", hi[0], hi[1]);
            Console.WriteLine("  # of surfs = {0} out of {1}\n", nsurf, sparta.surf.nline);
            Console.WriteLine("  # of surfs = {0}\n", nsurf);
            Console.WriteLine("  surfs:");
            for (int i = 0; i < nsurf; i++) Console.WriteLine(" {0}", surfs[i]);
            Console.WriteLine("\n");
        }
        //private void print_clines();
        //private void print_points();
        //private void print_loops();
    }
}
