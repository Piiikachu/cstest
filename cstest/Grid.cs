using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using bigint = System.Int64;
using cellint = System.Int32;
namespace cstest
{
    public partial class Grid
    {
        public const int DELTA = 8192;
        public const double BIG = 1.0e20;
        public const int MAXGROUP = 32;
        public const int MAXSPLITPERCELL = 10;

        // default value, can be overridden by global command

        public const int MAXSURFPERCELL = 100;

        enum Enum1{ XLO, XHI, YLO, YHI, ZLO, ZHI, INTERIOR };         // same as Domain
        enum Enum2{ PERIODIC, OUTFLOW, REFLECT, SURFACE, AXISYM };  // same as Domain
        enum Enum3{ REGION_ALL, REGION_ONE, REGION_CENTER };      // same as Surf

        // cell type = OUTSIDE/INSIDE/OVERLAP if entirely outside/inside surfs
        //   or has any overlap with surfs including grazing or touching
        // corner point = OUTSIDE/INSIDE (not OVERLAP) if outside/inside
        //   if exactly on surf, is still marked OUTSIDE or INSIDE by cut2d/cut3d
        //   corner pts are only set if cell type = OVERLAP

        enum Enum4{ UNKNOWN, OUTSIDE, INSIDE, OVERLAP };           // several files
        enum Enum5{ NCHILD, NPARENT, NUNKNOWN, NPBCHILD, NPBPARENT, NPBUNKNOWN, NBOUND };  // Update
        enum Enum6{ NOWEIGHT, VOLWEIGHT, RADWEIGHT };

        // allocate space for static class variable

        //Grid* Grid::gptr;

        // corners[i][j] = J corner points of face I of a grid cell
        // works for 2d quads and 3d hexes

        int[,] corners=new int[6,4]{{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, 
                     {0,1,2,3}, {4,5,6,7}};


        public int exist;            // 1 if grid is defined
        public int exist_ghost;      // 1 if ghost cells exist
        public int clumped;          // 1 if grid ownership is clumped, due to RCB
                               // if not, some operations are not allowed
         
        public bigint ncell;         // global count of child cells (unsplit+split, no sub)
        public bigint nunsplit;      // global count of unsplit cells
        public int nsplit;           // global count of split cells
        public int nsub;             // global count of split sub cells
        public int maxsurfpercell;   // max surf elements in one child cell
        public int maxlevel;         // max level of any child cell in grid, 0 = root
        public int uniform;          // 1 if all child cells are at same level, else 0
        public int unx, uny, unz;      // if uniform, effective global Nx,Ny,Nz of finest grid
        public double cutoff;        // cutoff for ghost cells, -1.0 = infinite
        public double cell_epsilon;  // half of smallest cellside of any cell in any dim
        public int cellweightflag;   // 0/1+ for no/yes usage of cellwise fnum weighting
         
        public int ngroup;               // # of defined groups
        public string[] gnames;            // name of each group
        public int[] bitmask;             // one-bit mask for each group
        public int[] inversemask;         // inverse mask for each group
         
        public int copy, copymode;    // 1 if copy of class (prevents deallocation of
                                      //  base class when child copy is destroyed)

        // hash for all cell IDs (owned,ghost,parent)

        public Dictionary<cellint, int> hash;

        public int hashfilled;             // 1 if hash is filled with parent/child IDs

        // list data structs

        MyPage<int> csurfs;        // lists of surf indices for
                                    // owned + ghost child cells with surfs
        MyPage<int> csplits;       // lists of sub cell offsets for
                                    // owned + ghost split info
        MyPage<int> csubs;         // lists of sub cell indices for
                                    // owned + ghost split info

        // owned or ghost child cell
        // includes unsplit cells, split cells, sub cells in any order
        // ghost cells are appended to owned

        public struct ChildCell
        {
            public cellint id;               // cell ID in bitwise format
            public int iparent;              // index of parent in pcells
            public int proc;                 // proc that owns this cell
            public int ilocal;               // index of this cell on owning proc
                                       // must be correct for all kinds of ghost cells
             
            public cellint[] neigh;         // info on 6 neighbor cells in cells/pcells
                                      //   that fully overlap face
                                      // order = XLO,XHI,YLO,YHI,ZLO,ZHI
                                      // nmask stores flags for what all 6 represent
                                      // if an index, store index into cells or pcells
                                      // if unknown, store ID of neighbor child cell
                                      // if non-periodic global boundary, ignored
            public int nmask;                // 3 bits for each of 6 values in neigh
                                             // 0 = index of child neighbor
                                             // 1 = index of parent neighbor
                                             // 2 = unknown child neighbor
                                             // 3 = index of PBC child neighbor
                                             // 4 = index of PBC parent neighbor
                                             // 5 = unknown PBC child neighbor
                                             // 6 = non-PBC boundary or ZLO/ZHI in 2d

            public double[] lo, hi;       // opposite corner pts of cell
            public int nsurf;                // # of surf elements in cell
                                             // -1 = empty ghost cell
            public int[] csurfs;              // indices of surf elements in cell
                                              // for sub cells, lo/hi/nsurf/csurfs
                                              //   are same as in split cell containing them

            public int nsplit;               // 1, unsplit cell
                                             // N > 1, split cell with N sub cells
                                             // N <= 0, neg of sub cell index (0 to Nsplit-1)
            public int isplit;               // index into sinfo
                                      // set for split and sub cells, -1 if unsplit
        };

        // info specific to owned child cell
        // includes unsplit cells, split cells, sub cells in same order as cells

        public struct ChildInfo
        {
            public int count;                // # of particles in this cell, 0 if split cell
            public int first;                // index of 1st particle in this cell, -1 if none
             
            public int mask;                 // grid group mask
            public int type;                 // OUTSIDE,INSIDE,OVERLAP,UNKNOWN
            public int[] corner;            // corner flags, 4/8 in 2d/3d
                                      // OUTSIDE,INSIDE,UNKNOWN
                                      // no OVERLAP is used for this anymore I think
                                      // ordered x first, y next, z last
                                      // for sub cells, type/corner
                                      //   are same as in split cell containing them

            public double volume;            // flow volume of cell or sub cell
                                      // entire cell volume for split cell
            public double weight;            // fnum weighting for this cell
        };

        // additional info for owned or ghost split cell
        // ghost split cell info is appended to owned split cell info

        public struct SplitInfo
        {
            public int icell;                // index of split cell in cells this belongs to
            public int xsub;                 // which sub cell (0 to Nsplit-1) xsplit is in
            public double[] xsplit;         // coords of point in split cell
            public int[] csplits;             // sub cell (0 to Nsplit-1) each Nsurf belongs to
            public int[] csubs;               // indices in cells of Nsplit sub cells
        };

        // parent cell
        // global list of parent cells is stored by all procs

        public struct ParentCell
        {
            public cellint id;               // cell ID in bitwise format, 0 = root
            public int mask;                 // grid group mask
            public int level;                // level in hierarchical grid, 0 = root
            public int nbits;                // # of bits to encode my ID, also my siblings
            public int newbits;              // # of additional bits to encode my children
            public int iparent;              // index of parent, -1 if id=root
            public int grandparent;          // 1 if this cell is a grandparent, 0 if not
            public int nx, ny, nz;             // sub grid within cell
            public double[] lo, hi;       // opposite corner pts of cell
        };

        public int nlocal;                 // # of child cells I own (all 3 kinds)
        public int nghost;                 // # of ghost child cells I store (all 3 kinds)
        public int nempty;                 // # of empty ghost cells I store
        public int nunsplitlocal;          // # of unsplit cells I own
        public int nunsplitghost;          // # of ghost unsplit cells I store
        public int nsplitlocal;            // # of split cells I own
        public int nsplitghost;            // # of ghost split cells I store
        public int nsublocal;              // # of sub cells I own
        public int nsubghost;              // # of ghost sub cells I store
        public int nparent;                // # of parent cells
          
        public int maxlocal;               // size of cinfo
          
        public ChildCell[] cells;           // list of owned and ghost child cells
        public ChildInfo[] cinfo;           // extra info for nlocal owned cells
        public SplitInfo[] sinfo;           // extra info for owned and ghost split cells
        public ParentCell[] pcells;         // global list of parent cells
         
         // restart buffers, filled by read_restart
         
        public int nlocal_restart;
        public cellint id_restart;
        public int[] nsplit_restart;

        private SPARTA sparta;

        // methods
        public Grid(SPARTA sparta)
        {
            this.sparta = sparta;
            exist = exist_ghost = clumped = 0;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);

            //gnames = (char**)memory.smalloc(MAXGROUP * sizeof(char*), "grid:gnames");
            //bitmask = (int*)memory.smalloc(MAXGROUP * sizeof(int), "grid:bitmask");
            //inversemask = (int*)memory.smalloc(MAXGROUP * sizeof(int),
            //                                      "grid:inversemask");

            gnames = new string[MAXGROUP];
            bitmask = new int[MAXGROUP];
            inversemask = new int[MAXGROUP];

            for (int i = 0; i < MAXGROUP; i++) bitmask[i] = 1 << i;
            for (int i = 0; i < MAXGROUP; i++) inversemask[i] = bitmask[i] ^ ~0;

            ngroup = 1;
            int n = "all".Length + 1;
            gnames[0] = string.Copy( "all");

            ncell = nunsplit = nsplit = nsub = 0;

            nlocal = nghost = maxlocal = maxcell = 0;
            nsplitlocal = nsplitghost = maxsplit = 0;
            nsublocal = nsubghost = 0;
            nparent = maxparent = 0;

            //cells = null;
            //cinfo = null;
            //sinfo = null;
            //pcells = null;

            maxbits = 8 * sizeof(cellint) - 1;

            maxsurfpercell = MAXSURFPERCELL;
            csurfs = null; csplits = null; csubs = null;
            allocate_surf_arrays();

            neighshift[(int)Enum1.XLO] = 0;
            neighshift[(int)Enum1.XHI] = 3;
            neighshift[(int)Enum1.YLO] = 6;
            neighshift[(int)Enum1.YHI] = 9;
            neighshift[(int)Enum1.ZLO] = 12;
            neighshift[(int)Enum1.ZHI] = 15;

            neighmask[(int)Enum1.XLO] = 7 << neighshift[(int)Enum1.XLO];
            neighmask[(int)Enum1.XHI] = 7 << neighshift[(int)Enum1.XHI];
            neighmask[(int)Enum1.YLO] = 7 << neighshift[(int)Enum1.YLO];
            neighmask[(int)Enum1.YHI] = 7 << neighshift[(int)Enum1.YHI];
            neighmask[(int)Enum1.ZLO] = 7 << neighshift[(int)Enum1.ZLO];
            neighmask[(int)Enum1.ZHI] = 7 << neighshift[(int)Enum1.ZHI];

            cutoff = -1.0;
            cellweightflag =(int)Enum6.NOWEIGHT;

            // allocate hash for cell IDs

            //# ifdef SPARTA_MAP
            //            hash = new std::map<cellint, int>();
            //#elif defined SPARTA_UNORDERED_MAP
            //            hash = new std::unordered_map<cellint, int>();
            //#else
            //            hash = new std::tr1::unordered_map<cellint, int>();
            //#endif

            hash = new Dictionary<cellint, cellint>();

            hashfilled = 0;
            copy = copymode = 0;
        }

        //      public void remove();
        public void init()
        {
            ncustom = sparta.particle.ncustom;
            //nbytes_particle = sizeof(Particle.OnePart);
            nbytes_particle = 96;
            nbytes_custom = sparta.particle.sizeof_custom();
            nbytes_total = nbytes_particle + nbytes_custom;
        }
        public void add_child_cell(cellint id, int iparent, double[] lo, double[] hi)
        {
            grow_cells(1, 1);

            int ncorner;
            if (sparta.domain.dimension == 3) ncorner = 8;
            else ncorner = 4;

            ChildCell c = cells[nlocal];
            c.id = id;
            c.iparent = iparent;
            c.proc = me;
            c.ilocal = nlocal;

            //c.lo = lo;
            //c.hi = hi;
            c.lo = new double[3];
            c.hi = new double[3];
            c.lo[0] = lo[0];
            c.lo[1] = lo[1];
            c.lo[2] = lo[2];
            c.hi[0] = hi[0];
            c.hi[1] = hi[1];
            c.hi[2] = hi[2];
            c.nsurf = 0;
            c.csurfs = null;
            c.nsplit = 1;
            c.isplit = -1;

            ChildInfo ci = cinfo[nlocal];
            ci.count = 0;
            ci.first = -1;
            ci.mask = 1;
            ci.type = (int) Enum4.OUTSIDE;
            int[] corner = new int[ncorner];
            for (int i = 0; i < ncorner; i++)
            {
                corner[i] = (int)Enum4.UNKNOWN;
            }
            ci.corner = corner;
            ci.weight = 1.0;

            if (sparta.domain.dimension == 3)
                ci.volume = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
            else if (sparta.domain.axisymmetric!=0)
                ci.volume = MyConst.MY_PI * (hi[1] * hi[1] - lo[1] * lo[1]) * (hi[0] - lo[0]);
            else
                ci.volume = (hi[0] - lo[0]) * (hi[1] - lo[1]);

            // increment both since are adding an unsplit cell
            cells[nlocal] = c;
            cinfo[nlocal] = ci;
            nunsplitlocal++;
            nlocal++;
        }


        public void add_parent_cell(cellint id, int iparent,
                           int nx, int ny, int nz, double[] lo, double[] hi)
        {
            grow_pcells(1);

            ParentCell p = pcells[nparent];
            p.id = id;
            p.mask = 1;
            if (iparent >= 0)
            {
                p.level = pcells[iparent].level + 1;
                p.nbits = pcells[iparent].nbits + pcells[iparent].newbits;
            }
            else p.level = p.nbits = 0;
            p.newbits = id_bits(nx, ny, nz);
            p.iparent = iparent;
            p.grandparent = 0;                // set by caller

            if (p.nbits + p.newbits > maxbits)
                sparta.error.one("Cell ID has too many bits");

            p.nx = nx;
            p.ny = ny;
            p.nz = nz;
            p.lo = new double[3];
            p.hi = new double[3];
            p.lo[0] = lo[0]; p.lo[1] = lo[1]; p.lo[2] = lo[2];
            p.hi[0] = hi[0]; p.hi[1] = hi[1]; p.hi[2] = hi[2];
            pcells[nparent] = p;

            nparent++;
        }
        public void add_split_cell(int ownflag)
        {
            grow_sinfo(1);
            if (ownflag!=0) nsplitlocal++;
            else nsplitghost++;
        }
        public void add_sub_cell(int icell, int ownflag)
        {
            grow_cells(1, 1);

            int inew;
            if (ownflag!=0) inew = nlocal;
            else inew = nlocal + nghost;

            cells[inew] = cells[icell];
            //memcpy(&cells[inew], &cells[icell], sizeof(ChildCell));
            if (ownflag!=0)
            {
                cinfo[inew] = cinfo[icell];
                //memcpy(&cinfo[inew], &cinfo[icell], sizeof(ChildInfo));
            }

            if (ownflag != 0)
            {
                nsublocal++;
                nlocal++;
            }
            else
            {
                nsubghost++;
                nghost++;
            }
        }
        public void remove_ghosts()
        {
            exist_ghost = 0;
            nghost = nunsplitghost = nsplitghost = nsubghost = 0;
        }
        public void setup_owned()
        {
            // global counts for ncell, nunsplit, nsplit, nsub

            bigint one = nlocal - nsublocal;
            sparta.mpi.MPI_Allreduce(ref one, ref ncell, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            one = nunsplitlocal = nlocal - nsplitlocal - nsublocal;
            sparta.mpi.MPI_Allreduce(ref one, ref nunsplit, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            sparta.mpi.MPI_Allreduce(ref nsplitlocal, ref nsplit, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            sparta.mpi.MPI_Allreduce(ref nsublocal, ref nsub, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            one = nunsplitlocal = nlocal - nsplitlocal - nsublocal;

            // set cell_epsilon to 1/2 the smallest dimension of any grid cell

            int dimension = sparta.domain.dimension;

            double eps = BIG;
            for (int i = 0; i < nlocal; i++)
            {
                if (cells[i].nsplit <= 0) continue;
                eps = Math.Min(eps, cells[i].hi[0] - cells[i].lo[0]);
                eps = Math.Min(eps, cells[i].hi[1] - cells[i].lo[1]);
                if (dimension == 3) eps = Math.Min(eps, cells[i].hi[2] - cells[i].lo[2]);
            }

            sparta.mpi.MPI_Allreduce(ref eps, ref cell_epsilon, 1, MPI.MPI_DOUBLE, MPI.MPI_MIN, sparta.world);
            cell_epsilon *= 0.5;
        }
        public void acquire_ghosts()
        {
            if (cutoff < 0.0) acquire_ghosts_all();
            else if (clumped!=0) acquire_ghosts_near();
            else if (sparta.comm.me == 0)
                sparta.error.one("Could not acquire nearby ghost cells b/c grid partition is not clumped");
        }
        public void rehash()
        {
            hash.Clear();
            for (int icell = 0; icell < nlocal + nghost; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                hash[cells[icell].id] = icell + 1;
            }
            for (int icell = 0; icell < nparent; icell++)
                hash[pcells[icell].id] = -(icell + 1);

            hashfilled = 1;
        }
        public void find_neighbors()
        {
            int icell, index, nmask, boundary, periodic;
            cellint[] neigh;
            cellint id;
            double[] lo,hi;
            double[] mout=new double[3];

            if (exist_ghost==0) return;

            int dim = sparta.domain.dimension;
            int[] bflag = sparta.domain.bflag;
            double[] boxlo = sparta.domain.boxlo;
            double[] boxhi = sparta.domain.boxhi;

            // insure all cell IDs (owned, ghost, parent) are hashed

            rehash();

            // set neigh flags and nmask for each owned and ghost child cell
            // sub cells have same lo/hi as split cell, so their neigh info is the same
            for (icell = 0; icell < nlocal + nghost; icell++)
            {
                lo = cells[icell].lo;
                hi = cells[icell].hi;
                neigh = new cellint[6];
                //neigh = cells[icell].neigh;
                nmask = 0;

                // generate a point cell_epsilon away from face midpoint, respecting PBC
                // id_find_face() walks from root cell to find the parent or child cell
                //   furthest down heirarchy containing pt and entire lo/hi face of icell

                // XLO

                mout[1] = 0.5 * (lo[1] + hi[1]);
                mout[2] = 0.5 * (lo[2] + hi[2]);

                if (lo[0] == boxlo[0]) boundary = 1;
                else boundary = 0;
                if (bflag[(int)Enum1.XLO] == (int)Enum2.PERIODIC) periodic = 1;
                else periodic = 0;

                if (boundary != 0 && periodic == 0)
                {
                    neigh[(int)Enum1.XLO] = 0;
                    nmask = neigh_encode((int)Enum5.NBOUND, nmask, (int)Enum1.XLO);
                }
                else
                {
                    if (boundary != 0) mout[0] = boxhi[0] - cell_epsilon;
                    else mout[0] = lo[0] - cell_epsilon;

                    id = id_find_face(mout, 0, 0, lo, hi);
                    if (hash[id] == hash[hash.Keys.Last()])
                    {
                        neigh[(int)Enum1.XLO] = id;
                        if (boundary != 0) nmask = neigh_encode((int)Enum5.NPBUNKNOWN, nmask, (int)Enum1.XLO);
                        else nmask = neigh_encode((int)Enum5.NUNKNOWN, nmask, (int)Enum1.XLO);
                    }
                    else
                    {
                        index = hash[id];
                        if (index > 0)
                        {
                            neigh[(int)Enum1.XLO] = index - 1;
                            if (boundary != 0) nmask = neigh_encode((int)Enum5.NPBCHILD, nmask, (int)Enum1.XLO);
                            else nmask = neigh_encode((int)Enum5.NCHILD, nmask, (int)Enum1.XLO);
                        }
                        else
                        {
                            neigh[(int)Enum1.XLO] = -index - 1;
                            if (boundary != 0) nmask = neigh_encode((int)Enum5.NPBPARENT, nmask, (int)Enum1.XLO);
                            else nmask = neigh_encode((int)Enum5.NPARENT, nmask, (int)Enum1.XLO);
                        }
                    }
                }

                // XHI
                if (hi[0] == boxhi[0]) boundary = 1;
                else boundary = 0;
                if (bflag[(int)Enum1.XHI] == (int)Enum2.PERIODIC) periodic = 1;
                else periodic = 0;
                if (boundary!=0 && periodic==0)
                {
                    neigh[(int)Enum1.XHI] = 0;
                    nmask = neigh_encode((int)Enum5.NBOUND, nmask, (int)Enum1.XHI);
                }
                else
                {
                    if ( boundary!=0 ) mout[0] = boxlo[0] + cell_epsilon;
                    else mout[0] = hi[0] + cell_epsilon;

                    id = id_find_face(mout, 0, 0, lo, hi);

                    if (hash[id] == hash[hash.Keys.Last()])
                    {
                        neigh[(int)Enum1.XHI] = id;
                        if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBUNKNOWN, nmask, (int)Enum1.XHI);
                        else nmask = neigh_encode((int)Enum5.NUNKNOWN, nmask, (int)Enum1.XHI);
                    }
                    else
                    {
                        index = hash[id];
                        if (index > 0)
                        {
                            neigh[(int)Enum1.XHI] = index - 1;
                            if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBCHILD, nmask, (int)Enum1.XHI);
                            else nmask = neigh_encode((int)Enum5.NCHILD, nmask, (int)Enum1.XHI);
                        }
                        else
                        {
                            neigh[(int)Enum1.XHI] = -index - 1;
                            if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBPARENT, nmask, (int)Enum1.XHI);
                            else nmask = neigh_encode((int)Enum5.NPARENT, nmask, (int)Enum1.XHI);
                        }
                    }
                }

                // YLO

                mout[0] = 0.5 * (lo[0] + hi[0]);
                mout[2] = 0.5 * (lo[2] + hi[2]);

                if (lo[1] == boxlo[1]) boundary = 1;
                else boundary = 0;
                if (bflag[(int)Enum1.YLO] == (int)Enum2.PERIODIC) periodic = 1;
                else periodic = 0;
                if (boundary!=0 && periodic==0)
                {
                    neigh[(int)Enum1.YLO] = 0;
                    nmask = neigh_encode((int)Enum5.NBOUND, nmask, (int)Enum1.YLO);
                }
                else
                {
                    if ( boundary!=0 ) mout[1] = boxhi[1] - cell_epsilon;
                    else mout[1] = lo[1] - cell_epsilon;

                    id = id_find_face(mout, 0, 1, lo, hi);

                    if (hash[id] == hash[hash.Keys.Last()])
                    {
                        neigh[(int)Enum1.YLO] = id;
                        if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBUNKNOWN, nmask, (int)Enum1.YLO);
                        else nmask = neigh_encode((int)Enum5.NUNKNOWN, nmask, (int)Enum1.YLO);
                    }
                    else
                    {
                        index = hash[id];
                        if (index > 0)
                        {
                            neigh[(int)Enum1.YLO] = index - 1;
                            if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBCHILD, nmask, (int)Enum1.YLO);
                            else nmask = neigh_encode((int)Enum5.NCHILD, nmask, (int)Enum1.YLO);
                        }
                        else
                        {
                            neigh[(int)Enum1.YLO] = -index - 1;
                            if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBPARENT, nmask, (int)Enum1.YLO);
                            else nmask = neigh_encode((int)Enum5.NPARENT, nmask, (int)Enum1.YLO);
                        }
                    }
                }

                // YHI

                if (hi[1] == boxhi[1]) boundary = 1;
                else boundary = 0;
                if (bflag[(int)Enum1.YHI] == (int)Enum2.PERIODIC) periodic = 1;
                else periodic = 0;

                if (boundary!=0 && periodic==0)
                {
                    neigh[(int)Enum1.YHI] = 0;
                    nmask = neigh_encode((int)Enum5.NBOUND, nmask, (int)Enum1.YHI);
                }
                else
                {
                    if ( boundary!=0 ) mout[1] = boxlo[1] + cell_epsilon;
                    else mout[1] = hi[1] + cell_epsilon;

                    id = id_find_face(mout, 0, 1, lo, hi);

                    if (hash[id] == hash[hash.Keys.Last()])
                    {
                        neigh[(int)Enum1.YHI] = id;
                        if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBUNKNOWN, nmask, (int)Enum1.YHI);
                        else nmask = neigh_encode((int)Enum5.NUNKNOWN, nmask, (int)Enum1.YHI);
                    }
                    else
                    {
                        index = hash[id];
                        if (index > 0)
                        {
                            neigh[(int)Enum1.YHI] = index - 1;
                            if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBCHILD, nmask, (int)Enum1.YHI);
                            else nmask = neigh_encode((int)Enum5.NCHILD, nmask, (int)Enum1.YHI);
                        }
                        else
                        {
                            neigh[(int)Enum1.YHI] = -index - 1;
                            if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBPARENT, nmask, (int)Enum1.YHI);
                            else nmask = neigh_encode((int)Enum5.NPARENT, nmask, (int)Enum1.YHI);
                        }
                    }
                }


                // ZLO
                // treat boundary as non-periodic if 2d, so is flagged as NBOUND

                mout[0] = 0.5 * (lo[0] + hi[0]);
                mout[1] = 0.5 * (lo[1] + hi[1]);

                if (lo[2] == boxlo[2]) boundary = 1;
                else boundary = 0;
                if (bflag[(int)Enum1.ZLO] == (int)Enum2.PERIODIC && dim == 3) periodic = 1;
                else periodic = 0;

                if (boundary!=0 && periodic==0)
                {
                    neigh[(int)Enum1.ZLO] = 0;
                    nmask = neigh_encode((int)Enum5.NBOUND, nmask, (int)Enum1.ZLO);
                }
                else
                {
                    if ( boundary!=0 ) mout[2] = boxhi[2] - cell_epsilon;
                    else mout[2] = lo[2] - cell_epsilon;

                    id = id_find_face(mout, 0, 2, lo, hi);

                    if (hash[id] == hash[hash.Keys.Last()])
                    {
                        neigh[(int)Enum1.ZLO] = id;
                        if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBUNKNOWN, nmask, (int)Enum1.ZLO);
                        else nmask = neigh_encode((int)Enum5.NUNKNOWN, nmask, (int)Enum1.ZLO);
                    }
                    else
                    {
                        index = hash[id];
                        if (index > 0)
                        {
                            neigh[(int)Enum1.ZLO] = index - 1;
                            if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBCHILD, nmask, (int)Enum1.ZLO);
                            else nmask = neigh_encode((int)Enum5.NCHILD, nmask, (int)Enum1.ZLO);
                        }
                        else
                        {
                            neigh[(int)Enum1.ZLO] = -index - 1;
                            if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBPARENT, nmask, (int)Enum1.ZLO);
                            else nmask = neigh_encode((int)Enum5.NPARENT, nmask, (int)Enum1.ZLO);
                        }
                    }
                }
                // ZHI
                // treat boundary as non-periodic if 2d, so is flagged as NBOUND

                if (hi[2] == boxhi[2]) boundary = 1;
                else boundary = 0;
                if (bflag[(int)Enum1.ZHI] == (int)Enum2.PERIODIC && dim == 3) periodic = 1;
                else periodic = 0;

                if (boundary!=0 && periodic==0)
                {
                    neigh[(int)Enum1.ZHI] = 0;
                    nmask = neigh_encode((int)Enum5.NBOUND, nmask, (int)Enum1.ZHI);
                }
                else
                {
                    if ( boundary!=0 ) mout[2] = boxlo[2] + cell_epsilon;
                    else mout[2] = hi[2] + cell_epsilon;

                    id = id_find_face(mout, 0, 2, lo, hi);

                    if (hash[id] == hash[hash.Keys.Last()])
                    {
                        neigh[(int)Enum1.ZHI] = id;
                        if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBUNKNOWN, nmask, (int)Enum1.ZHI);
                        else nmask = neigh_encode((int)Enum5.NUNKNOWN, nmask, (int)Enum1.ZHI);
                    }
                    else
                    {
                        index = hash[id];
                        if (index > 0)
                        {
                            neigh[(int)Enum1.ZHI] = index - 1;
                            if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBCHILD, nmask, (int)Enum1.ZHI);
                            else nmask = neigh_encode((int)Enum5.NCHILD, nmask, (int)Enum1.ZHI);
                        }
                        else
                        {
                            neigh[(int)Enum1.ZHI] = -index - 1;
                            if ( boundary!=0 ) nmask = neigh_encode((int)Enum5.NPBPARENT, nmask, (int)Enum1.ZHI);
                            else nmask = neigh_encode((int)Enum5.NPARENT, nmask, (int)Enum1.ZHI);
                        }
                    }
                }
                cells[icell].neigh = neigh;
                cells[icell].nmask = nmask;
            }
            // insure no UNKNOWN neighbors for owned cell
            // else cannot move particle to new proc to continue move
            int n1, n2, n3, n4, n5, n6;

            int flag = 0;
            for (icell = 0; icell < nlocal; icell++)
            {
                nmask = cells[icell].nmask;
                n1 = neigh_decode(nmask, (int)Enum1.XLO);
                n2 = neigh_decode(nmask, (int)Enum1.XHI);
                n3 = neigh_decode(nmask, (int)Enum1.YLO);
                n4 = neigh_decode(nmask, (int)Enum1.YHI);
                n5 = neigh_decode(nmask, (int)Enum1.ZLO);
                n6 = neigh_decode(nmask, (int)Enum1.ZHI);
                if (n1 == (int)Enum5.NUNKNOWN || n2 ==   (int)Enum5.NUNKNOWN || n3 ==   (int)Enum5.NUNKNOWN ||
                    n4 == (int)Enum5.NUNKNOWN || n5 ==   (int)Enum5.NUNKNOWN || n6 ==   (int)Enum5.NUNKNOWN) flag++;
                if (n1 == (int)Enum5.NPBUNKNOWN || n2 == (int)Enum5.NPBUNKNOWN || n3 == (int)Enum5.NPBUNKNOWN ||
                    n4 == (int)Enum5.NPBUNKNOWN || n5 == (int)Enum5.NPBUNKNOWN || n6 == (int)Enum5.NPBUNKNOWN) flag++;
            }
            int flagall=0;
            sparta.mpi.MPI_Allreduce(ref flag, ref flagall, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);

            if (flagall!=0)
            {
                string str = string.Format("Owned cells with unknown neighbors = {0}", flagall);
                sparta.error.all(str);
            }
        }
        public void unset_neighbors()
        {
            if (exist_ghost==0) return;

            // no change in neigh[] needed if nflag = NUNKNOWN, NPBUNKNOWN, or NBOUND

            int i, index, nmask, nflag;
            cellint[] neigh;

            for (int icell = 0; icell < nlocal; icell++)
            {
                neigh = cells[icell].neigh;
                nmask = cells[icell].nmask;

                for (i = 0; i < 6; i++)
                {
                    index = neigh[i];
                    nflag = neigh_decode(nmask, i);
                    if (nflag == (int)Enum5.NCHILD || nflag == (int)Enum5.NPBCHILD)
                        neigh[i] = cells[index].id;
                    else if (nflag == (int)Enum5.NPARENT || nflag == (int)Enum5.NPBPARENT)
                        neigh[i] = pcells[index].id;
                }
            }
        }
        public void reset_neighbors()
        {
            if (exist_ghost==0) return;

            // insure all cell IDs (owned, ghost, parent) are hashed

            rehash();

            // set neigh[] and nmask of each owned and ghost child cell
            // hash lookup can reset nmask to CHILD or UNKNOWN
            // no change in neigh[] or nmask needed if nflag = NBOUND

            int i, nmask, nflag;
            cellint[] neigh;

            for (int icell = 0; icell < nlocal + nghost; icell++)
            {
                neigh = cells[icell].neigh;
                nmask = cells[icell].nmask;

                for (i = 0; i < 6; i++)
                {
                    nflag = neigh_decode(nmask, i);
                    if (nflag ==(int)Enum5.NCHILD || nflag == (int)Enum5.NPBCHILD)
                    {
                        if (hash[neigh[i]] == hash[hash.Keys.Last()])
                        {
                            if (nflag == (int)Enum5.NCHILD) nmask = neigh_encode((int)Enum5.NUNKNOWN, nmask, i);
                            else nmask = neigh_encode((int)Enum5.NPBUNKNOWN, nmask, i);
                        }
                        else neigh[i] = hash[neigh[i]] - 1;

                    }
                    else if (nflag == (int)Enum5.NPARENT || nflag == (int)Enum5.NPBPARENT)
                    {
                        neigh[i] = -hash[neigh[i]] - 1;

                    }
                    else if (nflag == (int)Enum5.NUNKNOWN || nflag == (int)Enum5.NPBUNKNOWN)
                    {
                        if (hash[neigh[i]] != hash[hash.Keys.Last()])
                        {
                            neigh[i] = hash[neigh[i]] - 1;
                            if (nflag == (int)Enum5.NUNKNOWN) nmask = neigh_encode((int)Enum5.NCHILD, nmask, i);
                            else nmask = neigh_encode((int)Enum5.NPBCHILD, nmask, i);
                        }
                    }
                }

                cells[icell].nmask = nmask;
            }
        }
        public void set_inout()
        {
            int i, j, m, icell, jcell, pcell, itype, jtype, marktype;
            int nflag, iface, icorner, ctype, ic;
            int iset, nset, nsetnew;
            int nsend, maxsend, nrecv, maxrecv;
            int[] set,setnew,jcorner;
            int[] cflags;
            int[] proclist;
            double[] xcorner=new double[3];
            Connect[] sbuf,rbuf;

            if (exist_ghost==0)
                sparta.error.all("Cannot mark grid cells as inside/outside surfs because ghost cells do not exist");

            int[] faceflip = { (int)Enum1.XHI, (int)Enum1.XLO, (int)Enum1.YHI, (int)Enum1.YLO, (int)Enum1.ZHI, (int)Enum1.ZLO };

            int me = sparta.comm.me;
            int dimension = sparta.domain.dimension;
            int ncorner, nface, nface_pts;
            if (dimension == 3)
            {
                ncorner = 8;
                nface = 6;
                nface_pts = 4;
            }
            else
            {
                ncorner = 4;
                nface = 4;
                nface_pts = 2;
            }

            // create irregular communicator for exchanging off-processor cell types

            Irregular irregular = new Irregular(sparta);

            // create set1 and set2 lists so can swap between them

            int[] set1,set2;
            set1 = new int[nlocal* nface];
            set2 = new int[nlocal * nface];
            //memory->create(set1, nlocal * nface, "grid:set1");
            //memory->create(set2, nlocal * nface, "grid:set2");

            // initial set list = overlapped cells with corner values which are set

            nset = nsetnew = 0;
            set = set1;
            setnew = set2;

            for (icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                if (cinfo[icell].type == (int)Enum4.OVERLAP && cinfo[icell].corner[0] != (int)Enum4.UNKNOWN)
                    set[nset++] = icell;
            }

            // loop until no more cell types are set

            maxsend = maxrecv = 0;
            proclist = null;
            sbuf = rbuf = null;

            while (true)
            {
                nsend = 0;

                // process local set list
                // for unmarked neighbor cells I own, mark them and add to new set list
                // for neighbor cells I would mark but don't own, add to comm list

                while (nset!=0)
                {
                    nsetnew = 0;
                    for (iset = 0; iset < nset; iset++)
                    {
                        icell = set[iset];
                        itype = cinfo[icell].type;
                        cflags = cinfo[icell].corner;

                        // loop over my cell's faces

                        for (iface = 0; iface < nface; iface++)
                        {

                            // if I am an OUTSIDE/INSIDE cell: marktype = my itype
                            // if I am an OVERLAP cell: 
                            //   can only mark neighbor if all corner pts on face are same
                            //   marktype = value of those corner pts = OUTSIDE or INSIDE

                            if (itype != (int)Enum4.OVERLAP) marktype = itype;
                            else
                            {
                                ctype = cflags[corners[iface,0]];
                                for (m = 1; m < nface_pts; m++)
                                    if (cflags[corners[iface,m]] != ctype) break;
                                if (m < nface_pts) continue;
                                if (ctype == (int)Enum4.OUTSIDE) marktype = (int)Enum4.OUTSIDE;
                                else if (ctype == (int)Enum4.INSIDE) marktype = (int)Enum4.INSIDE;
                                else continue;
                            }

                            jcell = cells[icell].neigh[iface];
                            nflag = neigh_decode(cells[icell].nmask, iface);

                            if (nflag == (int)Enum5.NCHILD || nflag == (int)Enum5.NPBCHILD)
                            {

                                // this proc owns neighbor cell
                                // perform this marking logic (here and below):
                                // (1) if jtype = UNKNOWN:
                                //     set jtype = marktype = OUTSIDE/INSIDE
                                //     add jcell to new set
                                // (2) else if jtype = OVERLAP:
                                //     skip if jcell's corners are already marked
                                //     set jcell corner pts = marktype = OUTSIDE/INSIDE
                                //     also set its volume to full cell or zero
                                //     do not add jcell to new set, since is OVERLAP cell
                                //       maybe could add, but not sure need to, 
                                //       and could be marked wrong ?
                                // (3) else jcell = OUTSIDE/INSIDE:
                                //     if jcell mark is different than marktype,
                                //     error b/c markings are inconsistent

                                if (cells[jcell].proc == me)
                                {
                                    jtype = cinfo[jcell].type;

                                    if (jtype == (int)Enum4.UNKNOWN)
                                    {
                                        cinfo[jcell].type = marktype;
                                        setnew[nsetnew++] = jcell;
                                    }
                                    else if (jtype == (int)Enum4.OVERLAP)
                                    {
                                        // don't think this restriction needed (Mar 2017)
                                        //if (itype == OVERLAP) continue;
                                        jcorner = cinfo[jcell].corner;
                                        if (jcorner[0] != (int)Enum4.UNKNOWN) continue;
                                        for (icorner = 0; icorner < ncorner; icorner++)
                                            jcorner[icorner] = marktype;
                                        if (marktype == (int)Enum4.INSIDE) cinfo[jcell].volume = 0.0;
                                        else if (marktype == (int)Enum4.OUTSIDE)
                                        {
                                            double[] lo = cells[jcell].lo;
                                            double[] hi = cells[jcell].hi;
                                            if (dimension == 3)
                                                cinfo[jcell].volume =
                                                  (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
                                            else if (sparta.domain.axisymmetric!=0)
                                                cinfo[jcell].volume =
                                                  MyConst.MY_PI * (hi[1] * hi[1] - lo[1] * lo[1]) * (hi[0] - lo[0]);
                                            else
                                                cinfo[jcell].volume = (hi[0] - lo[0]) * (hi[1] - lo[1]);
                                        }
                                        // do not add to setnew, see comment above
                                        //setnew[nsetnew++] = jcell;
                                    }
                                    else
                                    {
                                        if (jtype != marktype)
                                        {
                                            Console.WriteLine("ICELL1 {0} id {1} iface {2} jcell {3} id {4} marktype {5} jtype {6}\n",
                                                   icell, cells[icell].id, iface, jcell, cells[jcell].id,
                                                   marktype, jtype);
                                            sparta.error.one("Cell type mis-match when marking on self");
                                        }
                                    }

                                    // this proc does not own neighbor cell
                                    // pack sbuf with info to send

                                }
                                else
                                {
                                    if (nsend == maxsend)
                                    {
                                        maxsend += DELTA;
                                        //memory->grow(proclist, maxsend, "grid:proclist");
                                        proclist = new int[maxsend];
                                        sbuf = new Connect[maxsend];
                                          //memory->srealloc(sbuf, maxsend * sizeof(Connect), "grid:sbuf");
                                    }
                                    proclist[nsend] = cells[jcell].proc;
                                    Connect cbuf = sbuf[nsend];
                                    cbuf.itype = itype;
                                    cbuf.marktype = marktype;
                                    cbuf.jcell = cells[jcell].ilocal;
                                    sbuf[nsend] = cbuf;
                                    nsend++;
                                }

                            }
                            else if (nflag == (int)Enum5.NPARENT || nflag == (int)Enum5.NPBPARENT)
                            {

                                // neighbor is a parent cell, not a child cell
                                // for each corner points of icell face:
                                //   find jcell = child cell owner of the face corner pt
                                //   if I own the child cell, mark it in same manner as above

                                pcell = jcell;
                                for (m = 0; m < nface_pts; m++)
                                {
                                    ic = corners[iface,m];
                                    if (ic % 2!=0) xcorner[0] = cells[icell].hi[0];
                                    else xcorner[0] = cells[icell].lo[0];
                                    if ((ic / 2) % 2 != 0) xcorner[1] = cells[icell].hi[1];
                                    else xcorner[1] = cells[icell].lo[1];
                                    if (ic / 4 != 0) xcorner[2] = cells[icell].hi[2];
                                    else xcorner[2] = cells[icell].lo[2];

                                    if (nflag == (int)Enum5.NPBPARENT)
                                        sparta.domain.uncollide(faceflip[iface], xcorner);

                                    jcell = id_find_child(pcell, xcorner);
                                    if (jcell < 0) sparta.error.one("Parent cell child missing");

                                    // this proc owns neighbor cell
                                    // perform same marking logic as above

                                    if (cells[jcell].proc == me)
                                    {
                                        jtype = cinfo[jcell].type;

                                        if (jtype == (int)Enum4.UNKNOWN)
                                        {
                                            cinfo[jcell].type = marktype;
                                            setnew[nsetnew++] = jcell;
                                        }
                                        else if (jtype == (int)Enum4.OVERLAP)
                                        {
                                            // don't think this restriction needed (Mar 2017)
                                            //if (itype == OVERLAP) continue;
                                            jcorner = cinfo[jcell].corner;
                                            if (jcorner[0] != (int)Enum4.UNKNOWN) continue;
                                            for (icorner = 0; icorner < ncorner; icorner++)
                                                jcorner[icorner] = marktype;
                                            if (marktype == (int)Enum4.INSIDE) cinfo[jcell].volume = 0.0;
                                            else if (marktype == (int)Enum4.OUTSIDE)
                                            {
                                                double[] lo = cells[jcell].lo;
                                                double[] hi = cells[jcell].hi;
                                                if (dimension == 3)
                                                    cinfo[jcell].volume =
                                                      (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
                                                else if (sparta.domain.axisymmetric!=0)
                                                    cinfo[jcell].volume =
                                                      MyConst.MY_PI * (hi[1] * hi[1] - lo[1] * lo[1]) * (hi[0] - lo[0]);
                                                else
                                                    cinfo[jcell].volume = (hi[0] - lo[0]) * (hi[1] - lo[1]);
                                            }
                                            // do not add to setnew, see comment above
                                            //setnew[nsetnew++] = jcell;
                                        }
                                        else
                                        {
                                            if (jtype != marktype)
                                            {
                                                Console.WriteLine("ICELL2 {0} id {1} iface {2} jcell {3} id {4} marktype {5} jtype {6}\n",
                                                       icell, cells[icell].id, iface, jcell, cells[jcell].id,
                                                       marktype, jtype);
                                                sparta.error.one(
                                                           "Cell type mis-match when marking on self");
                                            }
                                        }

                                        // this proc does not own neighbor cell
                                        // pack sbuf with info to send

                                    }
                                    else
                                    {
                                        if (nsend == maxsend)
                                        {
                                            maxsend += DELTA;
                                            proclist = new int[maxsend];
                                            //memory->grow(proclist, maxsend, "grid:proclist");
                                            sbuf = new Connect[maxsend];
                                              //memory->srealloc(sbuf, maxsend * sizeof(Connect), "grid:sbuf");
                                        }
                                        proclist[nsend] = cells[jcell].proc;
                                        sbuf[nsend].itype = itype;
                                        sbuf[nsend].marktype = marktype;
                                        sbuf[nsend].jcell = cells[jcell].ilocal;
                                        nsend++;
                                    }
                                }
                            }
                            else continue;
                        }
                    }

                    // swap set lists

                    nset = nsetnew;
                    if (set == set1)
                    {
                        set = set2;
                        setnew = set1;
                    }
                    else
                    {
                        set = set1;
                        setnew = set2;
                    }
                }

                // if no proc has info to communicate, then done iterating

                int anysend=0;
                sparta.mpi.MPI_Allreduce(ref nsend, ref anysend, 1, MPI.MPI_INT, MPI.MPI_MAX, sparta.world);
                if (anysend==0) break;

                // perform irregular comm of each proc's comm list
                // realloc rbuf as needed

                nrecv = irregular.create_data_uniform(nsend, proclist, sparta.comm.commsortflag);
                if (nrecv > maxrecv)
                {
                    //memory->sfree(rbuf);
                    maxrecv = nrecv;
                    rbuf = new Connect[maxrecv];
                    //rbuf = (Connect*)memory->smalloc(maxrecv * sizeof(Connect), "grid:rbuf");
                }

                irregular.exchange_uniform(sbuf.ToString(), 3*sizeof(int), rbuf.ToString());

                // this proc received info to mark neighbor cell it owns
                // perform same marking logic as above

                for (i = 0; i < nrecv; i++)
                {
                    itype = rbuf[i].itype;
                    marktype = rbuf[i].marktype;
                    jcell = rbuf[i].jcell;
                    jtype = cinfo[jcell].type;

                    if (jtype == (int)Enum4.UNKNOWN)
                    {
                        cinfo[jcell].type = marktype;
                        set[nset++] = jcell;
                    }
                    else if (jtype == (int)Enum4.OVERLAP)
                    {
                        // don't think this restriction needed (Mar 2017)
                        //if (itype == OVERLAP) continue;
                        jcorner = cinfo[jcell].corner;
                        if (jcorner[0] != (int)Enum4.UNKNOWN) continue;
                        for (icorner = 0; icorner < ncorner; icorner++)
                            jcorner[icorner] = marktype;
                        if (marktype == (int)Enum4.INSIDE) cinfo[jcell].volume = 0.0;
                        else if (marktype == (int)Enum4.OUTSIDE)
                        {
                            double[] lo = cells[jcell].lo;
                            double[] hi = cells[jcell].hi;
                            if (dimension == 3)
                                cinfo[jcell].volume =
                                  (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
                            else if (sparta.domain.axisymmetric!=0)
                                cinfo[jcell].volume =
                                  MyConst.MY_PI * (hi[1] * hi[1] - lo[1] * lo[1]) * (hi[0] - lo[0]);
                            else
                                cinfo[jcell].volume = (hi[0] - lo[0]) * (hi[1] - lo[1]);
                        }
                        // do not add to setnew, see comment above
                        //set[nset++] = jcell;
                    }
                    else
                    {
                        if (marktype != jtype)
                        {
                            Console.WriteLine("JCELL3 {0} id {1} marktype {2} jtype {3}\n",
                                   jcell, cells[jcell].id, marktype, jtype);
                            sparta.error.one("Cell type mis-match when marking on neigh proc");
                        }
                    }
                }
            }

            // NOTE: at this point could make a final attempt to mark
            //   any remaining UNKNOWN corner pts of an overlap cell
            //   to avoid warnings and errors in type_check()
            // when doing grid adaptation, a new cell could possibly still be unmarked,
            //   because its already marked neighbors which are non-OVERLAP cells
            //     will not try to mark it
            //   this is in contrast to first-time marking, when sweep entire grid
            // the final attempt logic could do this:
            //   instead of having marked cells mark their unmarked neighbors
            //   have unmarked cells look at their neighbors to acquire markings
            //   not necessary when doing initial full-grid sweep,
            //     but may be necessary when doing incremental adaptation

            // all done with marking
            // set type and cflags for all sub cells from split cell it belongs to

            int splitcell;

            for (icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit > 0) continue;
                splitcell = sinfo[cells[icell].isplit].icell;
                cinfo[icell].type = cinfo[splitcell].type;
                for (j = 0; j < ncorner; j++)
                    cinfo[icell].corner[j] = cinfo[splitcell].corner[j];
            }

            // set volume of cells that are now INSIDE to 0.0
            // this allows error check in Collide and FixGridCheck for particles
            //   in zero-volume cells

            for (icell = 0; icell < nlocal; icell++)
                if (cinfo[icell].type == (int)Enum4.INSIDE) cinfo[icell].volume = 0.0;

            /*
            printf("POST INOUT {}: {}\n",sparta.comm.me,grid->nlocal);
            for (int i = 0; i < grid->nlocal; i++) {
              Grid::ChildCell *g = &grid->cells[i];
              if (g->id == 52)
              printf("ICELL {}: {} id {} pid {} lo %g %g "
                     "hi %g %g type {} corners {} {} {} {} vol %g\n",
                     sparta.comm.me,i,g->id,pcells[g->iparent].id,
                     g->lo[0],
                     g->lo[1],
                     g->hi[0],
                     g->hi[1],
                     grid->cinfo[i].type,
                     grid->cinfo[i].corner[0],
                     grid->cinfo[i].corner[1],
                     grid->cinfo[i].corner[2],
                     grid->cinfo[i].corner[3],grid->cinfo[i].volume);
            }
            */

            // clean up

            //delete irregular;
            //memory->destroy(set1);
            //memory->destroy(set2);
            //memory->destroy(proclist);
            //memory->sfree(sbuf);
            //memory->sfree(rbuf);
        }
        public void check_uniform()
        {
            // maxlevel = max level of any child cell in grid

            maxlevel = 0;
            for (int i = 0; i < nparent; i++)
                maxlevel = Math.Max(maxlevel, pcells[i].level);
            maxlevel++;

            // grid is uniform only if parents of all child cells are at same level

            int plevel = -1;
            for (int i = 0; i < nlocal; i++)
                plevel = Math.Max(plevel, pcells[cells[i].iparent].level);

            int all=0;
            sparta.mpi.MPI_Allreduce(ref plevel, ref all, 1, MPI.MPI_INT, MPI.MPI_MAX, sparta.world);

            uniform = 1;
            for (int i = 0; i < nlocal; i++)
                if (pcells[cells[i].iparent].level != all) uniform = 0;

            sparta.mpi.MPI_Allreduce(ref uniform, ref all, 1,MPI.MPI_INT, MPI.MPI_MIN, sparta.world);
            if (all==0) uniform = 0;

            if (uniform!=0)
            {
                int[] lflag = new int[maxlevel];
                for (int i = 0; i < maxlevel; i++) lflag[i] = 0;

                int level;
                unx = uny = unz = 1;

                for (int i = 0; i < nparent; i++)
                {
                    level = pcells[i].level;
                    if (lflag[level]!=0) continue;
                    lflag[level] = 1;
                    unx *= pcells[i].nx;
                    uny *= pcells[i].ny;
                    unz *= pcells[i].nz;
                }
            }
        }
        public void type_check(int flag = 1)
        {
            int i;

            // check cell types

            int unknown = 0;
            for (int icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                if (cinfo[icell].type == (int)Enum4.UNKNOWN) unknown++;
            }
            int unknownall=0;
            sparta.mpi.MPI_Allreduce(ref unknown, ref unknownall, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);

            if (unknownall!=0)
            {
                string str=string.Format("Grid cells marked as unknown = {0}", unknownall);
                sparta.error.all(str);
            }

            // check corner flags of cells that are OVERLAP
            // warn if any interior corner flags are not set
            // error if any corner flags on global boundaries are unset

            double[] boxlo = sparta.domain.boxlo;
            double[] boxhi = sparta.domain.boxhi;
            int dimension = sparta.domain.dimension;

            int ncorner = 4;
            if (dimension == 3) ncorner = 8;

            double[] x=new double[3];
            int inside = 0;
            int outside = 0;

            for (int icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                if (cinfo[icell].type != (int)Enum4.OVERLAP) continue;
                for (i = 0; i < ncorner; i++)
                {
                    if (cinfo[icell].corner[i] != (int)Enum4.UNKNOWN) continue;
                    if (i % 2 == 0) x[0] = cells[icell].lo[0];
                    else x[0] = cells[icell].hi[0];
                    if ((i / 2) % 2 == 0) x[1] = cells[icell].lo[1];
                    else x[1] = cells[icell].hi[1];
                    if (dimension == 3)
                    {
                        if (i / 4 == 0) x[2] = cells[icell].lo[2];
                        else x[2] = cells[icell].hi[2];
                    }
                    else x[2] = 0.0;

                    if (Geometry.point_on_hex(x, boxlo, boxhi)!=0)
                    {
                        Console.WriteLine("BAD CORNER icell {0} id {1} type {3} icorner {4} x {5:G} {6:G} {7:G} cflags {8} {9} {10} {11} \n",
                               icell, cells[icell].id, cinfo[icell].type, i, x[0], x[1], x[2],
                               cinfo[icell].corner[0],
                               cinfo[icell].corner[1],
                               cinfo[icell].corner[2],
                               cinfo[icell].corner[3]);
                        outside++;
                    }
                    else inside++;
                }
            }

            int insideall=0;
            sparta.mpi.MPI_Allreduce(ref inside, ref insideall, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            if (insideall!=0)
            {
                string str = string.Format("Grid cell interior corner points marked as unknown (volume will be wrong if cell is effectively outside) = {0}",
                    insideall);
                if (sparta.comm.me == 0) sparta.error.all(str);
            }

            int outsideall=0;
            sparta.mpi.MPI_Allreduce(ref outside, ref outsideall, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            if (outsideall!=0)
            {
                string str = string.Format("Grid cell corner points on boundary marked as unknown = {0}",
                    outsideall);
                sparta.error.all(str);
            }

            int volzero = 0;
            for (int icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                if (cinfo[icell].type == (int)Enum4.OUTSIDE && cinfo[icell].volume == 0.0) volzero++;
            }
            int volzeroall=0;
            sparta.mpi.MPI_Allreduce(ref volzero, ref volzeroall, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            if (outsideall!=0)
            {
                string str = string.Format("Grid cells marked outside, but with zero volume = {0}",
                    volzeroall);
                sparta.error.all(str);
            }

            if (flag!=0) flow_stats();
        }






        //      public void refine_cell(int, int, int, int, int, int*,
        //                 class Cut2d *, class Cut3d *);
        //      public void coarsen_cell(int, int, int*, int*, int*, class AdaptGrid *,
        //                  class Cut2d *, class Cut3d *);

        //      public void group(int, char**);
        //      public int add_group(const char*);
        //      public int find_group(const char*);

        public virtual void grow_pcells(int n)
        {
            if (nparent + n >= maxparent)
            {
                int oldmax = maxparent;
                while (maxparent < nparent + n) maxparent += DELTA;
                pcells =new ParentCell[maxparent];
                
                //memset(ref pcells[oldmax], 0, (maxparent - oldmax) * sizeof(ParentCell));
            }
        }

        //      public void write_restart(FILE*);
        //      public void read_restart(FILE*);
        //      public int size_restart();
        //      public int pack_restart(char*);
        //      public int unpack_restart(char*);

        //      public bigint memory_usage();

        //      public void debug();

        //      // grid_comm.cpp



        //       // grid_surf.cpp



        //      // grid_id.cpp



        //      // extract/return neighbor flag for iface from per-cell nmask
        //      // inlined for efficiency

        int neigh_decode(int nmask, int iface)
        {
            return (nmask & neighmask[iface]) >> neighshift[iface];
        }

        //      // overwrite neighbor flag for iface in per-cell nmask
        //      // first line zeroes the iface bits via one's complement of mask
        //      // inlined for efficiency
        //      // return updated nmask

        int neigh_encode(int flag, int nmask, int iface)
        {
            nmask &= ~neighmask[iface];
            nmask |= flag << neighshift[iface];
            return nmask;
        }


        // protected:
        protected int me;
        protected int maxcell;             // size of cells
        protected int maxsplit;            // size of sinfo
        protected int maxparent;           // size of pcells
        protected int maxbits;             // max bits allowed in a cell ID
         
        protected int[] neighmask=new int[6];        // bit-masks for each face in nmask
        protected int[] neighshift=new int[6];       // bit-shifts for each face in nmask

        // connection between one of my cells and a neighbor cell on another proc

        protected struct Connect
        {
            public int itype;           // type of sending cell
            public int marktype;        // new type value (IN/OUT) for neighbor cell
            public int jcell;           // index of neighbor cell on receiving proc (owner)
        };

        // bounding box for a clump of grid cells

        protected struct Box
        {
            public double[] lo, hi;    // opposite corners of extended bbox
            public int proc;              // proc that owns it
        };

        // Particle class values used for packing/unpacking particles in grid comm

        protected int ncustom;
        protected int nbytes_particle, nbytes_custom, nbytes_total;

        // private methods

        void acquire_ghosts_all()
        {
            exist_ghost = 1;
            nempty = 0;

            // compute total # of ghosts so can pre-allocate cells array

            int nghost_new=0;
            sparta.mpi.MPI_Allreduce(ref nlocal,ref nghost_new, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            nghost_new -= nlocal;
            grow_cells(nghost_new, 0);

            // create buf for holding all of my cells, not including sub cells

            int sendsize = 0;
            for (int icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                StringBuilder buf = new StringBuilder();
                sendsize += pack_one(icell, 0, 0, 0,ref buf);
            }

            StringBuilder sbuf = new StringBuilder();
            

            // pack each unsplit or split cell
            // subcells will be packed by split cell

            sendsize = 0;
            for (int icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                
                sendsize += pack_one(icell, 0, 0, 1, ref sbuf);
            }

            // circulate buf of my grid cells around to all procs
            // unpack augments my ghost cells with info from other procs

            //gptr = this;
            sparta.comm.ring(sendsize, sizeof(char),sbuf, 1, unpack_ghosts, null, 0);

            //memory->destroy(sbuf);
        }
        void acquire_ghosts_near()
        {
            exist_ghost = 1;

            // bb lo/hi = bounding box for my owned cells

            int i;
            double[] bblo=new double[3], bbhi=new double[3];
            double[] lo,hi;

            for (i = 0; i < 3; i++)
            {
                bblo[i] = BIG;
                bbhi[i] = -BIG;
            }

            for (int icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                
                lo = new double[3];
                hi = new double[3];
                lo = cells[icell].lo;
                hi = cells[icell].hi;
                for (i = 0; i < 3; i++)
                {
                    bblo[i] = Math.Min(bblo[i], lo[i]);
                    bbhi[i] = Math.Max(bbhi[i], hi[i]);
                }
            }
            // ebb lo/hi = bbox + grid cutoff
            // trim to simulation box in non-periodic dims
            // if bblo/hi is at periodic boundary and cutoff is 0.0,
            //   add cell_epsilon to insure ghosts across periodic boundary acquired,
            //   else may be UNKNOWN to owned cell

            double[] boxlo = sparta.domain.boxlo;
            double[] boxhi = sparta.domain.boxhi;
            int[] bflag = sparta.domain.bflag;

            double[] ebblo=new double[3], ebbhi=new double[3];
            for (i = 0; i < 3; i++)
            {
                ebblo[i] = bblo[i] - cutoff;
                ebbhi[i] = bbhi[i] + cutoff;
                if (bflag[2 * i] != (int)Enum2.PERIODIC) ebblo[i] = Math.Max(ebblo[i], boxlo[i]);
                if (bflag[2 * i] != (int)Enum2.PERIODIC) ebbhi[i] = Math.Min(ebbhi[i], boxhi[i]);
                if (bflag[2 * i] == (int)Enum2.PERIODIC && bblo[i] == boxlo[i] && cutoff == 0.0)
                    ebblo[i] -= cell_epsilon;
                if (bflag[2 * i] == (int)Enum2.PERIODIC && bbhi[i] == boxhi[i] && cutoff == 0.0)
                    ebbhi[i] += cell_epsilon;
            }
            // box = ebbox split across periodic BC
            // 27 is max number of periodic images in 3d

            Box[] box=new Box[27];
            int nbox = box_periodic(ebblo, ebbhi, box);

            // boxall = collection of boxes from all procs

            int me = sparta.comm.me;
            int nprocs = sparta.comm.nprocs;

            int nboxall=0;
            sparta.mpi.MPI_Allreduce(ref nbox, ref nboxall, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);

            int[] recvcounts,displs;
            recvcounts = new int[nprocs];
            displs = new int[nprocs];
            //memory->create(recvcounts, nprocs, "grid:recvcounts");
            //memory->create(displs, nprocs, "grid:displs");

            int nsend = nbox * Marshal.SizeOf(typeof(Box));
            //sparta.mpi.MPI_Allgather(ref nsend, 1, MPI.MPI_INT, recvcounts, 1, MPI.MPI_INT, sparta.world);
            displs[0] = 0;
            for (i = 1; i < nprocs; i++) displs[i] = displs[i - 1] + recvcounts[i - 1];

            Box[] boxall = new Box[nboxall];
            //MPI_Allgatherv(box, nsend, MPI_CHAR, boxall, recvcounts, displs, MPI_CHAR, world);

            //memory->destroy(recvcounts);
            //memory->destroy(displs);

            // nlist = # of boxes that overlap with my bbox, skipping self boxes
            // list = indices into boxall of overlaps
            // overlap = true overlap or just touching

            int nlist = 0;
            int[] list;
            list = new int[nboxall];
            //memory->create(list, nboxall, "grid:list");

            for (i = 0; i < nboxall; i++)
            {
                if (boxall[i].proc == me) continue;
                if (box_overlap(bblo, bbhi, boxall[i].lo, boxall[i].hi)!=0) list[nlist++] = i;
            }

            // loop over my owned cells, not including sub cells
            // each may overlap with multiple boxes in list
            // on 1st pass, just tally memory to send copies of my cells
            // use lastproc to insure a cell only overlaps once per other proc
            // if oflag = 2 = my cell just touches box,
            // so flag grid cell as EMPTY ghost by setting nsurf = -1

            int j, oflag, lastproc, nsurf_hold=0;

            nsend = 0;
            int sendsize = 0;
            for (int icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                lo = cells[icell].lo;
                hi = cells[icell].hi;
                lastproc = -1;
                for (i = 0; i < nlist; i++)
                {
                    j = list[i];
                    oflag = box_overlap(lo, hi, boxall[j].lo, boxall[j].hi);
                    if (oflag!=0) continue;
                    if (boxall[j].proc == lastproc) continue;
                    lastproc = boxall[j].proc;

                    if (oflag == 2)
                    {
                        nsurf_hold = cells[icell].nsurf;
                        cells[icell].nsurf = -1;
                    }
                    StringBuilder bufempty = new StringBuilder();
                    sendsize += pack_one(icell, 0, 0, 0,ref bufempty);
                    if (oflag == 2) cells[icell].nsurf = nsurf_hold;
                    nsend++;
                }
            }

            StringBuilder sbuf = new StringBuilder();

            int[] proclist,sizelist;
            proclist = new int[nsend];
            sizelist = new int[nsend];
            //memory->create(proclist, nsend, "grid:proclist");
            //memory->create(sizelist, nsend, "grid:sizelist");

            // on 2nd pass over local cells, fill the send buf
            // use lastproc to insure a cell only overlaps once per other proc
            // if oflag = 2 = my cell just touches box,
            // so flag grid cell as EMPTY ghost by setting nsurf = -1

            nsend = 0;
            sendsize = 0;

            for (int icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                lo = cells[icell].lo;
                hi = cells[icell].hi;
                lastproc = -1;
                for (i = 0; i < nlist; i++)
                {
                    j = list[i];
                    oflag = box_overlap(lo, hi, boxall[j].lo, boxall[j].hi);
                    if (oflag!=0) continue;
                    if (boxall[j].proc == lastproc) continue;
                    lastproc = boxall[j].proc;

                    if (oflag == 2)
                    {
                        nsurf_hold = cells[icell].nsurf;
                        cells[icell].nsurf = -1;
                    }
                    sizelist[nsend] = pack_one(icell, 0, 0, 1, ref sbuf);
                    if (oflag == 2) cells[icell].nsurf = nsurf_hold;
                    proclist[nsend] = lastproc;
                    sendsize += sizelist[nsend];
                    nsend++;
                }
            }



            Irregular irregular = new Irregular(sparta);
            int recvsize;
            int nrecv = irregular.create_data_variable(nsend, proclist, sizelist,
                                                        out recvsize, sparta.comm.commsortflag);

            StringBuilder rbuf = new StringBuilder();
            //memory->create(rbuf, recvsize, "grid:rbuf");
            //memset(rbuf, 0, recvsize);
            byte[] sbytebuf = Encoding.UTF8.GetBytes(sbuf.ToString());
            byte[] rbytebuf = Encoding.UTF8.GetBytes(rbuf.ToString());
            irregular.exchange_variable(sbytebuf, sizelist, rbytebuf);
            //delete irregular;

            // unpack received grid cells as ghost cells

            int offset = 0;
            for (i = 0; i < nrecv; i++)
                offset += unpack_one(rbuf, 0, 0);

            // set nempty = # of EMPTY ghost cells I store

            nempty = 0;
            for (int icell = nlocal; icell < nlocal + nghost; icell++)
                if (cells[icell].nsurf < 0) nempty++;
        }

        //maybe ref?
        void box_intersect(double[] alo, double[] ahi, double[] blo, double[] bhi,double[] lo, double[] hi)
        {
            lo[0] = Math.Max(alo[0], blo[0]);
            hi[0] = Math.Min(ahi[0], bhi[0]);
            lo[1] = Math.Max(alo[1], blo[1]);
            hi[1] = Math.Min(ahi[1], bhi[1]);
            lo[2] = Math.Max(alo[2], blo[2]);
            hi[2] = Math.Min(ahi[2], bhi[2]);
        }
        int box_overlap(double[] alo, double[] ahi, double[] blo, double[] bhi)
        {
            double[] lo=new double[3], hi=new double[3];
            box_intersect(alo, ahi, blo, bhi, lo, hi);

            if (lo[0] > hi[0]) return 0;
            if (lo[1] > hi[1]) return 0;
            if (lo[2] > hi[2]) return 0;

            if (lo[0] == hi[0]) return 2;
            if (lo[1] == hi[1]) return 2;
            if (lo[2] == hi[2]) return 2;

            return 1;
        }
        int box_periodic(double[] lo, double[] hi, Box[] box)
        {
            int ilo, ihi, jlo, jhi, klo, khi;
            ilo = ihi = jlo = jhi = klo = khi = 0;

            int[] bflag = sparta.domain.bflag;
            if (bflag[0] ==(int)Enum2.PERIODIC)
            {
                ilo = -1; ihi = 1;
            }
            if (bflag[2] == (int)Enum2.PERIODIC)
            {
                jlo = -1; jhi = 1;
            }
            if (bflag[4] == (int)Enum2.PERIODIC && sparta.domain.dimension == 3)
            {
                klo = -1; khi = 1;
            }

            double[] boxlo = new double[3], boxhi = new double[3], prd = new double[3];
            boxlo = sparta.domain.boxlo;
            boxhi = sparta.domain.boxhi;
            prd =   sparta.domain.prd;

            double[] plo=new double[3], phi=new double[3], olo=new double[3], ohi=new double[3];

            int n = 0;

            int i, j, k;
            for (k = klo; k <= khi; k++)
            {
                for (j = jlo; j <= jhi; j++)
                {
                    for (i = ilo; i <= ihi; i++)
                    {
                        plo[0] = lo[0] + i * prd[0];
                        phi[0] = hi[0] + i * prd[0];
                        plo[1] = lo[1] + j * prd[1];
                        phi[1] = hi[1] + j * prd[1];
                        plo[2] = lo[2] + k * prd[2];
                        phi[2] = hi[2] + k * prd[2];
                        box_intersect(plo, phi, boxlo, boxhi, olo, ohi);
                        if (olo[0] >= ohi[0] || olo[1] >= ohi[1] || olo[2] >= ohi[2]) continue;

                        box[n].proc = me;
                        box[n].lo = olo;
                        box[n].hi = ohi;
                        //box[n].lo[0] = olo[0];
                        //box[n].lo[1] = olo[1];
                        //box[n].lo[2] = olo[2];
                        //box[n].hi[0] = ohi[0];
                        //box[n].hi[1] = ohi[1];
                        //box[n].hi[2] = ohi[2];
                        n++;
                    }
                }
            }

            return n;
        }

        protected virtual void grow_cells(int n, int m)
        {
            if (nlocal + nghost + n >= maxcell)
            {
                int oldmax = maxcell;
                while (maxcell < nlocal + nghost + n) maxcell += DELTA;
                cells = new ChildCell[maxcell];
                //memset(ref cells[oldmax], 0, (maxcell - oldmax) * sizeof(ChildCell));
            }

            if (nlocal + m >= maxlocal)
            {
                int oldmax = maxlocal;
                while (maxlocal < nlocal + m) maxlocal += DELTA;
                cinfo = new ChildInfo[maxlocal];
                  //memory->srealloc(cinfo, maxlocal * sizeof(ChildInfo), "grid:cinfo");
                //memset(&cinfo[oldmax], 0, (maxlocal - oldmax) * sizeof(ChildInfo));
            }
        }
        protected virtual void grow_sinfo(int n)
        {
            if (nsplitlocal + nsplitghost + n >= maxsplit)
            {
                int oldmax = maxsplit;
                while (maxsplit < nsplitlocal + nsplitghost + n) maxsplit += DELTA;
                sinfo = new SplitInfo[maxsplit];
                  //memory->srealloc(sinfo, maxsplit * sizeof(SplitInfo), "grid:sinfo");
                //memset(&sinfo[oldmax], 0, (maxsplit - oldmax) * sizeof(SplitInfo));
            }
        }


        
        

        //// grid_comm.cpp
        //// static variable for ring communication callback to access class data
        //// callback function for ring communication

        //static Grid* gptr;

        public void weight(int narg, string[] arg)
        {
            int i;

            if (exist != 0) sparta.error.all("Cannot weight cells before grid is defined");
            if (narg > 0 &&  narg != 1) sparta.error.all("Illegal weight command");

            // if called from read_restart with narg = -1, cellweightflag is already set

            if (narg == 1)
            {
                if (string.Equals(arg[0], "none")) cellweightflag = (int)Enum6.NOWEIGHT;
                else if (string.Equals(arg[0], "volume")) cellweightflag = (int)Enum6.VOLWEIGHT;
                else if (string.Equals(arg[0], "radius")) cellweightflag = (int)Enum6.RADWEIGHT;
                else sparta.error.all("Illegal weight command");
            }

            if (cellweightflag == (int)Enum6.RADWEIGHT &&  sparta.domain.axisymmetric != 0)
                sparta.error.all("Cannot use weight cell radius unless axisymmetric");

            // set per-cell weights

            for (i = 0; i < nlocal; i++) weight_one(i);
        }
        public void weight_one(int icell)
        {
            double[] lo, hi;

            int dimension = sparta.domain.dimension;
            int axisymmetric = sparta.domain.axisymmetric;

            if (cellweightflag == (int)Enum6.NOWEIGHT)
            {
                cinfo[icell].weight = 1.0;
            }
            else if (cellweightflag == (int)Enum6.VOLWEIGHT)
            {
                lo = cells[icell].lo;
                hi = cells[icell].hi;
                if (dimension == 3)
                    cinfo[icell].weight = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
                else if (axisymmetric!=0)
                    cinfo[icell].weight = MyConst.MY_PI * (hi[1] * hi[1] - lo[1] * lo[1]) * (hi[0] - lo[0]);
                else
                    cinfo[icell].weight = (hi[0] - lo[0]) * (hi[1] - lo[1]);
            }
            else if (cellweightflag == (int)Enum6.RADWEIGHT)
            {
                lo = cells[icell].lo;
                hi = cells[icell].hi;
                cinfo[icell].weight = 0.5 * (hi[1] + lo[1]) * (hi[0] - lo[0]);
            }
        }
    }

}