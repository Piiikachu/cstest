using System;
using System.Collections.Generic;
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
            int icell;                // index of split cell in cells this belongs to
            int xsub;                 // which sub cell (0 to Nsplit-1) xsplit is in
            double[] xsplit;         // coords of point in split cell
            int[] csplits;             // sub cell (0 to Nsplit-1) each Nsurf belongs to
            int[] csubs;               // indices in cells of Nsplit sub cells
        };

        // parent cell
        // global list of parent cells is stored by all procs

        public struct ParentCell
        {
            cellint id;               // cell ID in bitwise format, 0 = root
            int mask;                 // grid group mask
            int level;                // level in hierarchical grid, 0 = root
            int nbits;                // # of bits to encode my ID, also my siblings
            int newbits;              // # of additional bits to encode my children
            int iparent;              // index of parent, -1 if id=root
            int grandparent;          // 1 if this cell is a grandparent, 0 if not
            int nx, ny, nz;             // sub grid within cell
            double[] lo, hi;       // opposite corner pts of cell
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

            //gnames = (char**)memory->smalloc(MAXGROUP * sizeof(char*), "grid:gnames");
            //bitmask = (int*)memory->smalloc(MAXGROUP * sizeof(int), "grid:bitmask");
            //inversemask = (int*)memory->smalloc(MAXGROUP * sizeof(int),
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
        //      public void init();
        //      public void add_child_cell(cellint, int, double*, double*);
        //      public void add_parent_cell(cellint, int, int, int, int, double*, double*);
        //      public void add_split_cell(int);
        //      public void add_sub_cell(int, int);
        public void remove_ghosts()
        {
            exist_ghost = 0;
            nghost = nunsplitghost = nsplitghost = nsubghost = 0;
        }
        //      public void setup_owned();
        //      public void acquire_ghosts();
        //      public void rehash();
        //      public void find_neighbors();
        //      public void unset_neighbors();
        //      public void reset_neighbors();
        //      public void set_inout();
        //      public void check_uniform();
        //      public void type_check(int flag = 1);
       


    
    

    //      public void refine_cell(int, int, int, int, int, int*,
    //                 class Cut2d *, class Cut3d *);
    //      public void coarsen_cell(int, int, int*, int*, int*, class AdaptGrid *,
    //                  class Cut2d *, class Cut3d *);

    //      public void group(int, char**);
    //      public int add_group(const char*);
    //      public int find_group(const char*);

    //      public virtual void grow_pcells(int);

    //      public void write_restart(FILE*);
    //      public void read_restart(FILE*);
    //      public int size_restart();
    //      public int pack_restart(char*);
    //      public int unpack_restart(char*);

    //      public bigint memory_usage();

    //      public void debug();

    //      // grid_comm.cpp

    //      public int pack_one(int, char*, int, int, int);
    //      public int unpack_one(char*, int, int);
    //      public int pack_one_adapt(char*, char*, int);
    //      public int pack_particles(int, char*, int);
    //      public int unpack_particles(char*, int);
    //      public void unpack_particles_adapt(int, char*);
    //      public void compress();

    //       // grid_surf.cpp



    //      // grid_id.cpp

    //      int id_find_child(int, double*);
    //      int id_find_parent(cellint, cellint &);
    //      cellint id_str2num(char*);
    //      void id_num2str(cellint, char*);
    //      void id_pc_split(char*, char*, char*);
    //      void id_child_lohi(int, cellint, double*, double*);
    //      int id_bits(int, int, int);
    //      cellint id_find_face(double*, int, int, double*, double*);
    //      int id_child_from_parent_corner(int, int);

    //      // extract/return neighbor flag for iface from per-cell nmask
    //      // inlined for efficiency

    //      inline int neigh_decode(int nmask, int iface)
    //      {
    //          return (nmask & neighmask[iface]) >> neighshift[iface];
    //      }

    //      // overwrite neighbor flag for iface in per-cell nmask
    //      // first line zeroes the iface bits via one's complement of mask
    //      // inlined for efficiency
    //      // return updated nmask

    //      inline int neigh_encode(int flag, int nmask, int iface)
    //      {
    //          nmask &= ~neighmask[iface];
    //          nmask |= flag << neighshift[iface];
    //          return nmask;
    //      }


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
            int itype;           // type of sending cell
            int marktype;        // new type value (IN/OUT) for neighbor cell
            int jcell;           // index of neighbor cell on receiving proc (owner)
        };

        // bounding box for a clump of grid cells

        protected struct Box
        {
            double[] lo, hi;    // opposite corners of extended bbox
            int proc;              // proc that owns it
        };

        // Particle class values used for packing/unpacking particles in grid comm

        protected int ncustom;
        protected int nbytes_particle, nbytes_custom, nbytes_total;

        // private methods

        //void acquire_ghosts_all();
        //void acquire_ghosts_near();

        //void box_intersect(double*, double*, double*, double*,
        //                   double*, double*);
        //int box_overlap(double*, double*, double*, double*);
        //int box_periodic(double*, double*, Box*);

        //virtual void grow_cells(int, int);
        //virtual void grow_sinfo(int);

        //void surf2grid_stats();
        //void flow_stats();
        //double flow_volume();

        //// grid_comm.cpp
        //// static variable for ring communication callback to access class data
        //// callback function for ring communication

        //static Grid* gptr;
        //static void unpack_ghosts(int, char*);
        public void weight(int narg, string[] arg)
        {
            int i;

            if (exist != 0) sparta.error.all("Cannot weight cells before grid is defined");
            if (narg > 0 && narg != 1) sparta.error.all("Illegal weight command");

            // if called from read_restart with narg = -1, cellweightflag is already set

            if (narg == 1)
            {
                if (string.Equals(arg[0], "none")) cellweightflag = (int)Enum6.NOWEIGHT;
                else if (string.Equals(arg[0], "volume")) cellweightflag = (int)Enum6.VOLWEIGHT;
                else if (string.Equals(arg[0], "radius")) cellweightflag = (int)Enum6.RADWEIGHT;
                else sparta.error.all("Illegal weight command");
            }

            if (cellweightflag == (int)Enum6.RADWEIGHT && sparta.domain.axisymmetric != 0)
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