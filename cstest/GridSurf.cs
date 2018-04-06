using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public partial class Grid
    {
        //      public void surf2grid(int, int outflag = 1);
        //      public void surf2grid_one(int, int, int, int, class Cut3d *, class Cut2d *);
        //      public void clear_surf();
        //      public void clear_surf_restart();
        //      public void combine_split_cell_particles(int, int);
        //      public void assign_split_cell_particles(int);
        public void allocate_surf_arrays()
        {
            csurfs = new MyPage<int>(maxsurfpercell, Math.Max(100 * maxsurfpercell, 1024));
            csplits = new MyPage<int>(maxsurfpercell, Math.Max(100 * maxsurfpercell, 1024));
            csubs = new MyPage<int>(MAXSPLITPERCELL, 128);
        }
        //      int* csubs_request(int);
    }
}
