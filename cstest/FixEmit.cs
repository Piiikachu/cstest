using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class FixEmit : Fix
    {
        enum Enum1{ UNKNOWN, OUTSIDE, INSIDE, OVERLAP };           // same as Grid

        public const int DELTAGRID = 1024;
        public const int DELTACELL = 1024;

        public override int setmask()
        {
            int mask = 0;
            mask |= START_OF_STEP;
            return mask;
        }
        public override void init()
        {
            sparta.particle.exist = 1;
            ntotal = 0;

            int dimension = sparta.domain.dimension;
            Grid.ChildCell[] cells = sparta.grid.cells;
            Grid.ChildInfo[] cinfo = sparta.grid.cinfo;
            nglocal = sparta.grid.nlocal;

            if (nglocal > nglocalmax)
            {
                nglocalmax = nglocal;
                c2list = new int[nglocalmax];
                //memory.create(c2list, nglocalmax, "emit:c2list");
            }

            // upsplit, split, sub cells store c2list flag
            // upsplit, split cells can store clist data, but only if have tasks
            // no tasks for a cell inside surface
            // no tasks if cell is entirely outside region bounding box

            int flag, ntaskcell, ntaskfirst;

            nlist = 0;
            ntaskfirst = 0;

            for (int icell = 0; icell < nglocal; icell++)
            {
                c2list[icell] = -1;
                if (cells[icell].nsplit <= 0) continue;
                if (cinfo[icell].type == (int) Enum1.INSIDE) continue;
                if (region!=null && region.bboxflag!=0)
                {
                    flag = 1;
                    if (cells[icell].hi[0] > region.extent_xlo &&
                        cells[icell].lo[0] < region.extent_xhi) flag = 0;
                    if (cells[icell].hi[1] > region.extent_ylo &&
                        cells[icell].lo[1] < region.extent_yhi) flag = 0;
                    if (dimension == 3)
                    {
                        if (cells[icell].hi[2] > region.extent_zlo &&
                            cells[icell].lo[2] < region.extent_zhi) flag = 0;
                    }
                    if (flag!=0) continue;
                }

                ntaskcell = create_task(icell);
                if (ntaskcell)
                {
                    if (nlist == nlistmax)
                    {
                        nlistmax += DELTAGRID;
                        memory.grow(clist, nlistmax, "emit:clist");
                        memory.grow(clistnum, nlistmax, "emit:clistnum");
                        memory.grow(clistfirst, nlistmax, "emit:clistfirst");
                    }
                    c2list[icell] = nlist;
                    clist[nlist] = icell;
                    clistnum[nlist] = ntaskcell;
                    clistfirst[nlist] = ntaskfirst;
                    ntaskfirst += ntaskcell;
                    nlist++;
                }
            }

            active_current = 0;
        }
        //public void start_of_step();
        //public double compute_vector(int);

        //public void add_grid_one(int, int);
        //public int pack_grid_one(int, char*, int);
        //public int unpack_grid_one(int, char*);
        //public void compress_grid();
        //public virtual void post_compress_grid() { }
        private SPARTA sparta;
        public FixEmit(SPARTA sparta, int narg, string[] arg) : base(sparta, narg, arg)
        {
            this.sparta = sparta;
            vector_flag = 1;
            size_vector = 2;
            global_freq = 1;
            gridmigrate = 1;

            // RNG

            int me = sparta.comm.me;
            random = new RanPark(sparta.update.ranmaster.uniform());
            double seed = sparta.update.ranmaster.uniform();
            random.reset(seed, me, 100);

            // local storage of emit data structures

            c2list = null;
            nglocal = nglocalmax = 0;
            clist = clistnum = clistfirst = null;
            nlist = nlistmax = 0;

            // counters common to all emit styles for output from fix

            nsingle = ntotal = 0;
        }

        protected int perspecies;
        protected Region region;
        protected RanPark random;
        protected int nsingle, ntotal;
         
        protected int nglocal;         // copy of cell.nlocal
        protected int nglocalmax;      // max size of c2list
        protected int[] c2list;         // index into clist for each owned cell
                              // -1 if no tasks for the cell
                              // only for unsplit and split cells
         
        protected int nlist;           // # of owned cells with insert tasks
        protected int nlistmax;        // max size of clist,clistnum,clistfirst
        protected int[] clist;          // local indices of cells with insert tasks
        protected int[] clistnum;       // # of insert tasks in each cell with tasks
        protected int[] clistfirst;     // first task ID for each cell with tasks
         
        protected int ntask;           // # of insert tasks in underlying child class
         
        protected int active_current;  // set to 0 if grid cell data struct changes
                                       // triggers rebuild of active cell list in child classes

        protected virtual int create_task(int icell)
        {
            return 0;
        }
        //protected virtual void perform_task() = 0;
        //protected virtual int pack_task(int, char*, int) = 0;
        //protected virtual int unpack_task(char*, int) = 0;
        //protected virtual void copy_task(int, int, int, int) = 0;
           
        //protected void grow_percell(int);
        //protected void grow_list();
        //protected double mol_inflow(double, double, double);
        //protected int subsonic_temperature_check(int, double);
        protected void options(int narg, string[] arg)
        {
            nevery = 1;
            perspecies = 1;
            region = null;

            int iarg = 0;

            while (iarg<narg)
            {
                switch (arg[iarg])
                {
                    case "nevery":
                        if (iarg + 2 > narg) sparta.error.all("Illegal fix emit command");
                        nevery = int.Parse(arg[iarg + 1]);
                        if (nevery <= 0) sparta.error.all("Illegal fix emit command");
                        iarg += 2;
                        break;
                    case "perspecies":
                        if (iarg + 2 > narg) sparta.error.all("Illegal fix emit command");
                        if (string.Equals(arg[iarg + 1], "yes")) perspecies = 1;
                        else if (string.Equals(arg[iarg + 1], "no") ) perspecies = 0;
                        else sparta.error.all("Illegal fix emit command");
                        iarg += 2;
                        break;
                    case "region":
                        if (iarg + 2 > narg) sparta.error.all("Illegal fix emit command");
                        int iregion = sparta.domain.find_region(arg[iarg + 1]);
                        if (iregion < 0)
                            sparta.error.all("Fix emit region does not exist");
                        region = sparta.domain.regions[iregion];
                        iarg += 2;
                        break;
                    default:
                        string[] argss = new string[narg - iarg];
                        Array.Copy(arg, iarg, argss, 0, narg - iarg);
                        iarg += option(narg - iarg, argss);
                        break;
                }



            }


        }
        protected virtual int option(int narg, string[] arg)
        {
            sparta.error.all("Illegal fix emit command");
            return 0;
        }


    }
}
