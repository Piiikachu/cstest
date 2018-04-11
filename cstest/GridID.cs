using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using bigint = System.Int64;
using cellint = System.Int32;

namespace cstest
{
    public partial class Grid
    {
        //      int id_find_child(int, double*);
        //      int id_find_parent(cellint, cellint &);
        //      cellint id_str2num(char*);
        //      void id_num2str(cellint, char*);
        //      void id_pc_split(char*, char*, char*);

        public void id_child_lohi(int iparent, cellint ichild, double[] lo, double[] hi)
        {
            ParentCell p = pcells[iparent];
            ichild--;

            int nx = p.nx;
            int ny = p.ny;
            int nz = p.nz;

            int ix = ichild % nx;
            int iy = (ichild / nx) % ny;
            int iz = ichild / (nx * ny);

            double[] plo = p.lo;
            double[] phi = p.hi;

            lo[0] = plo[0] + ix * (phi[0] - plo[0]) / nx;
            lo[1] = plo[1] + iy * (phi[1] - plo[1]) / ny;
            lo[2] = plo[2] + iz * (phi[2] - plo[2]) / nz;

            hi[0] = plo[0] + (ix + 1) * (phi[0] - plo[0]) / nx;
            hi[1] = plo[1] + (iy + 1) * (phi[1] - plo[1]) / ny;
            hi[2] = plo[2] + (iz + 1) * (phi[2] - plo[2]) / nz;

            if (ix == nx - 1) hi[0] = phi[0];
            if (iy == ny - 1) hi[1] = phi[1];
            if (iz == nz - 1) hi[2] = phi[2];
        }

        int id_bits(int nx, int ny, int nz)
        {
            bigint n = ((bigint)nx) * ny * nz;
            bigint nstore = 1;
            int nbits = 1;
            while (nstore < n)
            {
                nstore = 2 * nstore + 1;
                nbits++;
            }
            return nbits;
        }
        //      cellint id_find_face(double*, int, int, double*, double*);
        //      int id_child_from_parent_corner(int, int);


    }
}
