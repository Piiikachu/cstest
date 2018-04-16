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
        int id_find_child(int iparent, double[] x)
        {
            ParentCell p = pcells[iparent];
            double[] lo = p.lo;
            double[] hi = p.hi;
            int nx = p.nx;
            int ny = p.ny;
            int nz = p.nz;
            int ix =(int)((x[0] - lo[0]) * nx / (hi[0] - lo[0]));
            int iy =(int)((x[1] - lo[1]) * ny / (hi[1] - lo[1]));
            int iz =(int)((x[2] - lo[2]) * nz / (hi[2] - lo[2]));
            if (ix == nx) ix--;
            if (iy == ny) iy--;
            if (iz == nz) iz--;

            cellint ichild = iz * nx * ny + iy * nx + ix + 1;
            cellint idchild = p.id | (ichild << p.nbits);

            if (hash[idchild] == hash[hash.Keys.Last()]) return -1;
            int index = hash[idchild];
            if (index > 0) return index - 1;
            return id_find_child(-index - 1, x);
        }
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
        cellint id_find_face(double[] x, int icell, int dim, double[] olo, double[] ohi)
        {
            ParentCell p = pcells[icell];
            double[] lo = p.lo;
            double[] hi = p.hi;

            // go down a level to find (ix,iy,iz) of new cell that contains pt x

            int nx = p.nx;
            int ny = p.ny;
            int nz = p.nz;
            int ix =(int)((x[0] - lo[0]) * nx / (hi[0] - lo[0]));
            int iy =(int)((x[1] - lo[1]) * ny / (hi[1] - lo[1]));
            int iz =(int)((x[2] - lo[2]) * nz / (hi[2] - lo[2]));
            if (ix == nx) ix--;
            if (iy == ny) iy--;
            if (iz == nz) iz--;

            // calculate lo/hi of new cell
            // exact same math as in id_child_lohi()

            double[] newlo=new double[3], newhi=new double[3];

            newlo[0] = lo[0] + ix * (hi[0] - lo[0]) / nx;
            newlo[1] = lo[1] + iy * (hi[1] - lo[1]) / ny;
            newlo[2] = lo[2] + iz * (hi[2] - lo[2]) / nz;

            newhi[0] = lo[0] + (ix + 1) * (hi[0] - lo[0]) / nx;
            newhi[1] = lo[1] + (iy + 1) * (hi[1] - lo[1]) / ny;
            newhi[2] = lo[2] + (iz + 1) * (hi[2] - lo[2]) / nz;

            if (ix == nx - 1) newhi[0] = hi[0];
            if (iy == ny - 1) newhi[1] = hi[1];
            if (iz == nz - 1) newhi[2] = hi[2];

            // if new cell does not fully overlap olo/ohi face, return parent ID

            if (dim != 0 && (newlo[0] > olo[0] || newhi[0] < ohi[0])) return p.id;
            if (dim != 1 && (newlo[1] > olo[1] || newhi[1] < ohi[1])) return p.id;
            if (dim != 2 && (newlo[2] > olo[2] || newhi[2] < ohi[2])) return p.id;

            // id = ID of new cell
            // if I don't store new ID, it's a child ID, return it
            // if I do store new ID, determine if parent or child cell
            // if child, return it
            // if parent, recurse

            cellint ichild = ((cellint)iz) * nx * ny + ((cellint)iy) * nx + ix + 1;
            cellint id = p.id | (ichild << p.nbits);
            if (hash[id] == hash[hash.Keys.Last()]) return id;
            icell = hash[id];
            if (icell > 0) return id;
            return id_find_face(x, -icell - 1, dim, olo, ohi);
        }
        //      int id_child_from_parent_corner(int, int);


    }
}
