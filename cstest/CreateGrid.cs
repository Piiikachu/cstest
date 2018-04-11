using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using bigint = System.Int64;
using cellint = System.Int32;

namespace cstest
{
    public class CreateGrid
    {
        enum Enum1{ NONE, LEVEL, STRIDE, CLUMP, BLOCK, RANDOM };
        enum Enum2 { XYZ, XZY, YXZ, YZX, ZXY, ZYX };
        enum Enum3 { ANY, ALL };
        private SPARTA sparta;
        public CreateGrid(SPARTA sparta)
        {
            this.sparta = sparta;
        }
        public void command(int narg, string[] args)
        {
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);
            if (sparta.domain.box_exist == 0)
                sparta.error.all("Cannot create grid before simulation box is defined");
            if (sparta.grid.exist != 0)
                sparta.error.all("Cannot create grid when grid is already defined");

            sparta.grid.exist = 1;

            if (narg < 3) sparta.error.all("Illegal create_grid command");

            int nx = int.Parse(arg[0]);
            int ny = int.Parse(arg[1]);
            int nz = int.Parse(arg[2]);

            if (nx < 1 || ny < 1 || nz < 1)
                sparta.error.all("Illegal create_grid command");
            if (sparta.domain.dimension == 2 && nz != 1)
                sparta.error.all("Create_grid nz value must be 1 for a 2d simulation");

            // optional args

            dimension = sparta.domain.dimension;

            int nlevels = 1;
            int bstyle = (int)Enum1.NONE;
            int px = 0;
            int py = 0;
            int pz = 0;
            int order=0;
            int inside = (int)Enum3.ANY;

            int iarg = 3;
            while (iarg < narg)
            {
                if (string.Equals(arg[iarg], "level"))
                {
                    if (iarg + 8 > narg) sparta.error.all("Illegal create_grid command");
                    if (bstyle != (int)Enum1.NONE && bstyle != (int)Enum1.LEVEL)
                        sparta.error.all("Illegal create_grid command");
                    bstyle = (int)Enum1.LEVEL;
                    if (int.Parse(arg[iarg + 1]) != nlevels + 1)
                        sparta.error.all("Illegal create_grid command");
                    nlevels++;
                    iarg += 8;

                }
                else if (string.Equals(arg[iarg], "region"))
                {
                    if (iarg + 6 > narg) sparta.error.all("Illegal create_grid command");
                    if (bstyle != (int)Enum1.NONE && bstyle != (int)Enum1.LEVEL)
                        sparta.error.all("Illegal create_grid command");
                    bstyle = (int)Enum1.LEVEL;
                    if (int.Parse(arg[iarg + 1]) != nlevels + 1)
                        sparta.error.all("Illegal create_grid command");
                    if (sparta.domain.find_region(arg[iarg + 2]) < 0)
                        sparta.error.all("Create_grid region ID does not exist");
                    nlevels++;
                    iarg += 6;

                }
                else if (string.Equals(arg[iarg], "stride"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal create_grid command");
                    if (bstyle != (int)Enum1.NONE) sparta.error.all("Illegal create_grid command");
                    bstyle = (int)Enum1.STRIDE;
                    if (string.Equals(arg[iarg + 1], "xyz")) order = (int)Enum2.XYZ;
                    else if (string.Equals(arg[iarg + 1], "xzy")) order = (int)Enum2.XZY;
                    else if (string.Equals(arg[iarg + 1], "yxz")) order = (int)Enum2.YXZ;
                    else if (string.Equals(arg[iarg + 1], "yzx")) order = (int)Enum2.YZX;
                    else if (string.Equals(arg[iarg + 1], "zxy")) order = (int)Enum2.ZXY;
                    else if (string.Equals(arg[iarg + 1], "zyx")) order = (int)Enum2.ZYX;
                    else sparta.error.all("Illegal create_grid command");
                    iarg += 2;

                }
                else if (string.Equals(arg[iarg], "clump"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal create_grid command");
                    if (bstyle != (int)Enum1.NONE) sparta.error.all("Illegal create_grid command");
                    bstyle = (int)Enum1.CLUMP;
                    if (string.Equals(arg[iarg + 1], "xyz")) order = (int)Enum2.XYZ;
                    else if (string.Equals(arg[iarg + 1], "xzy")) order = (int)Enum2.XZY;
                    else if (string.Equals(arg[iarg + 1], "yxz")) order = (int)Enum2.YXZ;
                    else if (string.Equals(arg[iarg + 1], "yzx")) order = (int)Enum2.YZX;
                    else if (string.Equals(arg[iarg + 1], "zxy")) order = (int)Enum2.ZXY;
                    else if (string.Equals(arg[iarg + 1], "zyx")) order = (int)Enum2.ZYX;
                    else sparta.error.all("Illegal create_grid command");
                    iarg += 2;

                }
                else if (string.Equals(arg[iarg], "block"))
                {
                    if (iarg + 4 > narg) sparta.error.all("Illegal create_grid command");
                    if (bstyle != (int)Enum1.NONE) sparta.error.all("Illegal create_grid command");
                    bstyle = (int)Enum1.BLOCK;
                    if (string.Equals(arg[iarg + 1], "*")) px = 0;
                    else px = int.Parse(arg[iarg + 1]);
                    if (string.Equals(arg[iarg + 2], "*")) py = 0;
                    else py = int.Parse(arg[iarg + 2]);
                    if (string.Equals(arg[iarg + 3], "*")) pz = 0;
                    else pz = int.Parse(arg[iarg + 3]);
                    iarg += 4;

                }
                else if (string.Equals(arg[iarg], "random"))
                {
                    if (iarg + 1 > narg) sparta.error.all("Illegal create_grid command");
                    if (bstyle != (int)Enum1.NONE) sparta.error.all("Illegal create_grid command");
                    bstyle = (int)Enum1.RANDOM;
                    iarg += 1;

                }
                else if (string.Equals(arg[iarg], "inside"))
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal create_grid command");
                    if (string.Equals(arg[iarg + 1], "any")) inside = (int)Enum3.ANY;
                    else if (string.Equals(arg[iarg + 1], "all")) inside = (int)Enum3.ALL;
                    else sparta.error.all("Illegal create_grid command");
                    iarg += 2;

                }
                else sparta.error.all("Illegal create_grid command");

            }
            if (bstyle == (int)Enum1.NONE) bstyle = (int)Enum1.LEVEL;

            // partition domain across procs, only for BLOCK style

            if (bstyle == (int)Enum1.BLOCK)
            {
                procs2grid(nx, ny, nz, ref px, ref py, ref pz);
                if (px * py * pz != sparta.comm.nprocs)
                    sparta.error.all("Bad grid of processors for create_grid");
            }

            // create root parent cell
            // treat first 3 args as if specified as "level 1 * * * Nx Ny Nz"

            sparta.mpi.MPI_Barrier(sparta.world);
            double time1 = sparta.mpi.MPI_Wtime();

            int level = 1;
            int xlo, xhi, ylo, yhi, zlo, zhi;
            xlo = xhi = ylo = yhi = zlo = zhi = 1;
            iarg = 3;
            Region region = null;


            // loop over levels
            // new level determines assignment of previous-level cells
            //   to be parent cells vs child cells
            // parent cells are those which are further partitioned by new level
            // child cells are those that are not
            // add all cells in previous level as either parent or child cells
            // if this is last level, also add all current level cells as child cells

            int me = sparta.comm.me;
            int nprocs = sparta.comm.nprocs;
            bigint count = 0;

            int pnx, pny, pnz, ix, iy, iz, nbits, pflag, proc;
            cellint m, nth=0, idgrandparent, idparent, idchild;
            double[] lo = new double[3], hi = new double[3];
            Grid.ParentCell p;

            while (true)
            {
                // add previous level cells as parent or child cells
                // loop over all parent cells to find ones two levels up
                // use their info to generate parent or child cells at previous level
                // decision on parent vs child in previous level depends on 
                //   pxyz lo/hi bounds in this level
                //   or on whether parent is in/out of region

                if (level == 1)
                {
                    sparta.grid.add_parent_cell(0, -1, nx, ny, nz, sparta.domain.boxlo, sparta.domain.boxhi);

                }
                else
                {
                    int nparent1 = sparta.grid.nparent;
                    int prevlevel = level - 2;

                    for (int igrandparent = 0; igrandparent < nparent1; igrandparent++)
                    {
                        if (sparta.grid.pcells[igrandparent].level != prevlevel) continue;
                        p = sparta.grid.pcells[igrandparent];

                        idgrandparent = p.id;
                        nbits = p.nbits;
                        pnx = p.nx;
                        pny = p.ny;
                        pnz = p.nz;

                        m = 0;
                        for (iz = 0; iz < pnz; iz++)
                            for (iy = 0; iy < pny; iy++)
                                for (ix = 0; ix < pnx; ix++)
                                {
                                    m++;
                                    idparent = idgrandparent | (m << nbits);
                                    sparta.grid.id_child_lohi(igrandparent, m, lo, hi);
                                    if (region != null) pflag = cell_in_region(lo, hi, region, inside);
                                    else
                                    {
                                        pflag = 1;
                                        if (ix + 1 < xlo || ix + 1 > xhi) pflag = 0;
                                        if (iy + 1 < ylo || iy + 1 > yhi) pflag = 0;
                                        if (iz + 1 < zlo || iz + 1 > zhi) pflag = 0;
                                    }
                                    if (pflag != 0)
                                        sparta.grid.add_parent_cell(idparent, igrandparent, nx, ny, nz, lo, hi);
                                    else
                                    {
                                        if (count % nprocs == me)
                                            sparta.grid.add_child_cell(idparent, igrandparent, lo, hi);
                                        count++;
                                    }
                                }
                    }


                }

                // final level, add current level cells as child cells
                // loop over all parent cells to find ones at previous level
                // use their info to generate my child cells at this level
                // if BSTYLE is set, there is only 1 level, create proc's cells directly

                if (level == nlevels)
                {
                    Grid.ParentCell[] pcells = sparta.grid.pcells;
                    int nparent2 = sparta.grid.nparent;
                    int prevlevel = level - 1;

                    for (int iparent = 0; iparent < nparent2; iparent++)
                    {
                        if (pcells[iparent].level != prevlevel) continue;
                        p = pcells[iparent];
                        idparent = p.id;
                        nbits = p.nbits;
                        nx = p.nx;
                        ny = p.ny;
                        nz = p.nz;

                        if (bstyle == (int)Enum1.LEVEL)
                        {
                            cellint ntotal = (cellint)nx * ny * nz;
                            int firstproc = (int)count % nprocs;
                            cellint ifirst = me - firstproc + 1;
                            if (ifirst <= 0) ifirst += nprocs;
                            for (m = ifirst; m <= ntotal; m += nprocs)
                            {
                                idchild = idparent | (m << nbits);
                                sparta.grid.id_child_lohi(iparent, m, lo, hi);
                                sparta.grid.add_child_cell(idchild, iparent, lo, hi);
                            }
                            count += ntotal;

                            // loop over all child cells
                            // convert M to Nth based on order
                            // assign each cell to proc based on Nth and STRIDE or CLUMP

                        }
                        else if (bstyle == (int)Enum1.STRIDE || bstyle == (int)Enum1.CLUMP)
                        {
                            cellint ntotal = (cellint)nx * ny * nz;
                            for (m = 0; m < ntotal; m++)
                            {
                                ix = m % nx;
                                iy = (m / nx) % ny;
                                iz = m / (nx * ny);
                                if (order == (int)Enum2.XYZ) nth = (cellint)iz * nx * ny + iy * nx + ix;
                                else if (order == (int)Enum2.XZY) nth = (cellint)iy * nx * nz + iz * nx + ix;
                                else if (order == (int)Enum2.YXZ) nth = (cellint)iz * ny * nx + ix * ny + iy;
                                else if (order == (int)Enum2.YZX) nth = (cellint)ix * ny * nz + iz * ny + iy;
                                else if (order == (int)Enum2.ZXY) nth = (cellint)iy * nz * nx + ix * nz + iz;
                                else if (order == (int)Enum2.ZYX) nth = (cellint)ix * nz * ny + iy * nz + iz;
                                nth++;
                                if (bstyle == (int)Enum1.STRIDE) proc = nth % nprocs;
                                else proc = Convert.ToInt32(1.0 * nth / ntotal * nprocs);
                                if (proc != me) continue;
                                idchild = idparent | (nth << nbits);
                                sparta.grid.id_child_lohi(iparent, nth, lo, hi);
                                sparta.grid.add_child_cell(idchild, iparent, lo, hi);
                            }
                            count += ntotal;

                            // loop over subset of cells in my BLOCK

                        }
                        else if (bstyle == (int)Enum1.BLOCK)
                        {
                            int ipx = me % px;
                            int ipy = (me / px) % py;
                            int ipz = me / (px * py);

                            int ixstart = Convert.ToInt32(1.0 * ipx / px * nx);
                            int ixstop = Convert.ToInt32(1.0 * (ipx + 1) / px * nx);
                            int iystart = Convert.ToInt32(1.0 * ipy / py * ny);
                            int iystop = Convert.ToInt32(1.0 * (ipy + 1) / py * ny);
                            int izstart = Convert.ToInt32(1.0 * ipz / pz * nz);
                            int izstop = Convert.ToInt32(1.0 * (ipz + 1) / pz * nz);

                            for (iz = izstart; iz < izstop; iz++)
                            {
                                for (iy = iystart; iy < iystop; iy++)
                                {
                                    for (ix = ixstart; ix < ixstop; ix++)
                                    {
                                        m = (cellint)iz * nx * ny + iy * nx + ix;
                                        m++;
                                        idchild = idparent | (m << nbits);
                                        sparta.grid.id_child_lohi(iparent, m, lo, hi);
                                        sparta.grid.add_child_cell(idchild, iparent, lo, hi);
                                    }
                                }
                            }



                        }
                        else if (bstyle == (int) Enum1.RANDOM)
                        {
                            RanPark random = new RanPark(sparta.update.ranmaster.uniform());
                            cellint ntotal = (cellint)nx * ny * nz;
                            for (m = 0; m < ntotal; m++)
                            {
                                proc =Convert.ToInt32(nprocs * random.uniform());
                                if (proc != me) continue;
                                idchild = idparent | ((m + 1) << nbits);
                                sparta.grid.id_child_lohi(iparent, m + 1, lo, hi);
                                sparta.grid.add_child_cell(idchild, iparent, lo, hi);
                            }
                            count += ntotal;
                            
                        }

                    }
                    break;

                }
                if (level == nlevels) break;
                // args for next level
                // levels must be in ascending order starting at beginning of arg list

                level++;

                if (string.Equals(arg[iarg], "level"))
                {
                    bounds(arg[iarg + 2], nx,ref xlo,ref xhi);
                    bounds(arg[iarg + 3], ny,ref ylo,ref yhi);
                    bounds(arg[iarg + 4], nz, ref zlo,ref zhi);
                    nx = int.Parse(arg[iarg + 5]);
                    ny = int.Parse(arg[iarg + 6]);
                    nz = int.Parse(arg[iarg + 7]);
                    iarg += 8;
                }
                else if (string.Equals(arg[iarg], "region"))
                {
                    int iregion = sparta.domain.find_region(arg[iarg + 2]);
                    region = sparta.domain.regions[iregion];
                    nx = int.Parse(arg[iarg + 3]);
                    ny = int.Parse(arg[iarg + 4]);
                    nz = int.Parse(arg[iarg + 5]);
                    iarg += 6;
                }
                else sparta.error.all("Illegal create_grid command");

                if (nx < 1 || ny < 1 || nz < 1)
                    sparta.error.all("Illegal create_grid command");
                if (dimension == 2)
                {
                    if (zlo != 1 || zhi != 1 || nz != 1)
                        sparta.error.all("Illegal create_grid command");
                }

            }

            // set grandparent flag for all parent cells

            Grid.ParentCell[] pcells1 = sparta.grid.pcells;
            int nparent = sparta.grid.nparent;

            for (int i = 1; i < nparent; i++)
                pcells1[pcells1[i].iparent].grandparent = 1;

            // invoke grid methods to complete grid setup

            if (nprocs == 1 || bstyle == (int)Enum1.CLUMP || bstyle == (int)Enum1.BLOCK) sparta.grid.clumped = 1;
            else sparta.grid.clumped = 0;

            sparta.mpi.MPI_Barrier(sparta.world);
            double time2 = sparta.mpi.MPI_Wtime();

            sparta.grid.setup_owned();
            sparta.grid.acquire_ghosts();
            sparta.grid.find_neighbors();
            sparta.grid.check_uniform();
            sparta.comm.reset_neighbors();

            sparta.mpi.MPI_Barrier(sparta.world);
            double time3 = sparta.mpi.MPI_Wtime();

            // stats

            double time_total = time3 - time1;

            if (sparta.comm.me == 0)
            {

                string str1 = string.Format("Created {0:G} child grid cells\n", sparta.grid.ncell);
                string str2=string.Format("  parent cells = {0}\n", sparta.grid.nparent);
                string str3 = string.Format("  CPU time = {0:G} secs\n", time_total);
                string str4 = string.Format("  create/ghost percent = {0:G} {1:G}\n", 100.0 * (time2 - time1) / time_total, 100.0 * (time3 - time2) / time_total);
                if (sparta.screen!=null)
                {
                    StreamWriter sw=  new StreamWriter(sparta.screen);
                    sw.Write(str1);
                    sw.Write(str2);
                    sw.Write(str3);
                    sw.Write(str4);

                }

                if (sparta.logfile!=null)
                {
                    StreamWriter sw1 = new StreamWriter(sparta.logfile);
                    sw1.Write(str1);
                    sw1.Write(str2);
                    sw1.Write(str3);
                    sw1.Write(str4);
                }
            }

        }


        private int dimension;

        private void bounds(string str, int nmax, ref int nlo, ref int nhi)
        {
            Console.WriteLine("CreateGrid.bounds {0}",str);
            //char* ptr = strchr(str, '*');

            //if (ptr == NULL)
            //{
            //    nlo = nhi = atoi(str);
            //}
            //else if (strlen(str) == 1)
            //{
            //    nlo = 1;
            //    nhi = nmax;
            //}
            //else if (ptr == str)
            //{
            //    nlo = 1;
            //    nhi = atoi(ptr + 1);
            //}
            //else if (strlen(ptr + 1) == 0)
            //{
            //    nlo = atoi(str);
            //    nhi = nmax;
            //}
            //else
            //{
            //    nlo = atoi(str);
            //    nhi = atoi(ptr + 1);
            //}

            //printf("BOUNDS %s nmax %d nlo/hi %d %d\n",str,nmax,nlo,nhi);

            if (nlo < 1 || nhi > nmax || nlo > nhi)
                sparta.error.all("Numeric index is out of bounds");
        }
        private int cell_in_region(double[] lo, double[] hi,Region region, int inside)
        {
            double[] x=new double[3];

            int n = 0;

            x[0] = lo[0]; x[1] = lo[1]; x[2] = lo[2];
            n += region.match(x);
            x[0] = hi[0]; x[1] = lo[1]; x[2] = lo[2];
            n += region.match(x);
            x[0] = lo[0]; x[1] = hi[1]; x[2] = lo[2];
            n += region.match(x);
            x[0] = hi[0]; x[1] = hi[1]; x[2] = lo[2];
            n += region.match(x);

            if (dimension == 3)
            {
                x[0] = lo[0]; x[1] = lo[1]; x[2] = hi[2];
                n += region.match(x);
                x[0] = hi[0]; x[1] = lo[1]; x[2] = hi[2];
                n += region.match(x);
                x[0] = lo[0]; x[1] = hi[1]; x[2] = hi[2];
                n += region.match(x);
                x[0] = hi[0]; x[1] = hi[1]; x[2] = hi[2];
                n += region.match(x);
            }

            if (inside == (int)Enum3.ANY)
            {
                if (n!=0) return 1;
            }
            else if (inside == (int)Enum3.ALL)
            {
                if (dimension == 2 && n == 4) return 1;
                if (dimension == 3 && n == 8) return 1;
            }

            return 0;
        }
        void procs2grid(int nx, int ny, int nz,ref int px,ref int py, ref int pz)
        {
            int upx = px;
            int upy = py;
            int upz = pz;

            int nprocs = sparta.comm.nprocs;

            // all 3 proc counts are specified

            if (px!=0 && py != 0 && pz != 0) return;

            // 2 out of 3 proc counts are specified

            if (py > 0 && pz > 0)
            {
                px = nprocs / (py * pz);
                return;
            }
            else if (px > 0 && pz > 0)
            {
                py = nprocs / (px * pz);
                return;
            }
            else if (px > 0 && py > 0)
            {
                pz = nprocs / (px * py);
                return;
            }

            // determine cross-sectional areas
            // area[0] = xy, area[1] = xz, area[2] = yz

            double[] area=new double[3];
            area[0] = nx * ny;
            area[1] = nx * nz;
            area[2] = ny * nz;

            double bestsurf = 2.0 * (area[0] + area[1] + area[2]);

            // loop thru all possible factorizations of nprocs
            // only consider valid cases that match procgrid settings
            // surf = surface area of a proc sub-domain

            int ipx, ipy, ipz, valid;
            double surf;

            ipx = 1;

            while (ipx <= nprocs)
            {
                valid = 1;
                if (upx!=0 && ipx != upx) valid = 0;
                if (nprocs % ipx!=0) valid = 0;
                if (valid==0)
                {
                    ipx++;
                    continue;
                }

                ipy = 1;
                while (ipy <= nprocs / ipx)
                {
                    valid = 1;
                    if (upy!=0 && ipy != upy) valid = 0;
                    if ((nprocs / ipx) % ipy!=0) valid = 0;
                    if (valid==0)
                    {
                        ipy++;
                        continue;
                    }

                    ipz = nprocs / ipx / ipy;
                    valid = 1;
                    if (upz != 0 && ipz != upz) valid = 0;
                    if (dimension == 2 && ipz != 1) valid = 0;
                    if (valid==0)
                    {
                        ipy++;
                        continue;
                    }

                    surf = area[0] / ipx / ipy + area[1] / ipx / ipz + area[2] / ipy / ipz;
                    if (surf < bestsurf)
                    {
                        bestsurf = surf;
                        px = ipx;
                        py = ipy;
                        pz = ipz;
                    }
                    ipy++;
                }

                ipx++;
            }




        }
    }
}
