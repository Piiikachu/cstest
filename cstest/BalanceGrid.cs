using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using cellint = System.Int32;
using bigint = System.Int64;
using System.IO;

namespace cstest
{
    class BalanceGrid
    {
        //#define RCB_DEBUG 1     // un-comment to include RCB proc boxes in image

        enum Enum1{ NONE, STRIDE, CLUMP, BLOCK, RANDOM, PROC, BISECTION };
        enum Enum2 { XYZ, XZY, YXZ, YZX, ZXY, ZYX };
        enum Enum3 { CELL, PARTICLE };

        public const double ZEROPARTICLE = 0.1;
        private SPARTA sparta;
        public BalanceGrid(SPARTA sparta)
        {
            this.sparta = sparta;
        }
        public void command(int narg, string[] args, int outflag = 1)
        {
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);
            if (sparta.grid.exist ==0)
                sparta.error.all( "Cannot balance grid before grid is defined");

            if (narg < 1) sparta.error.all( "Illegal balance_grid command");

            int order=0, bstyle=0;
            int px = 0, py = 0, pz = 0;
            int rcbwt=0, rcbflip=0;

            if (string.Equals(arg[0], "none")  )
            {
                if (narg != 1) sparta.error.all( "Illegal balance_grid command");
                bstyle = (int)Enum1.NONE;

            }
            else if (string.Equals(arg[0], "stride")  )
            {
                if (narg != 2) sparta.error.all( "Illegal balance_grid command");
                bstyle = (int)Enum1.STRIDE;
                if (string.Equals(arg[1], "xyz")  ) order = (int)Enum2.XYZ;
                else if (string.Equals(arg[1], "xzy")  ) order = (int)Enum2.XZY;
                else if (string.Equals(arg[1], "yxz")  ) order = (int)Enum2.YXZ;
                else if (string.Equals(arg[1], "yzx")  ) order = (int)Enum2.YZX;
                else if (string.Equals(arg[1], "zxy")  ) order = (int)Enum2.ZXY;
                else if (string.Equals(arg[1], "zyx")  ) order = (int)Enum2.ZYX;
                else sparta.error.all( "Illegal balance_grid command");

            }
            else if (string.Equals(arg[0], "clump")  )
            {
                if (narg != 2) sparta.error.all( "Illegal balance_grid command");
                bstyle = (int)Enum1.CLUMP;
                if (string.Equals(arg[1], "xyz")  ) order = (int)Enum2.XYZ;
                else if (string.Equals(arg[1], "xzy")  ) order = (int)Enum2.XZY;
                else if (string.Equals(arg[1], "yxz")  ) order = (int)Enum2.YXZ;
                else if (string.Equals(arg[1], "yzx")  ) order = (int)Enum2.YZX;
                else if (string.Equals(arg[1], "zxy")  ) order = (int)Enum2.ZXY;
                else if (string.Equals(arg[1], "zyx")  ) order = (int)Enum2.ZYX;
                else sparta.error.all( "Illegal balance_grid command");

            }
            else if (string.Equals(arg[0], "block")  )
            {
                if (narg != 4) sparta.error.all( "Illegal balance_grid command");
                bstyle = (int)Enum1.BLOCK;
                if (string.Equals(arg[1], "*")  ) px = 0;
                else px = int.Parse(arg[1]);
                if (string.Equals(arg[2], "*")  ) py = 0;
                else py = int.Parse(arg[2]);
                if (string.Equals(arg[3], "*")  ) pz = 0;
                else pz = int.Parse(arg[3]);

            }
            else if (string.Equals(arg[0], "random")  )
            {
                if (narg != 1) sparta.error.all( "Illegal balance_grid command");
                bstyle = (int)Enum1.RANDOM;

            }
            else if (string.Equals(arg[0], "proc")  )
            {
                if (narg != 1) sparta.error.all( "Illegal balance_grid command");
                bstyle = (int)Enum1.PROC;

            }
            else if (string.Equals(arg[0], "rcb")  )
            {
                if (narg != 2 && narg != 3)
                    sparta.error.all( "Illegal balance_grid command");
                bstyle = (int)Enum1.BISECTION;
                if (string.Equals(arg[1], "cell")  ) rcbwt = (int)Enum3.CELL;
                else if (string.Equals(arg[1], "part")  ) rcbwt = (int)Enum3.PARTICLE;
                else sparta.error.all( "Illegal balance_grid command");
                // undocumented optional 3rd arg
                // rcbflip = 3rd arg = 1 forces rcb.compute() to flip sign
                //           of all grid cell "dots" to force totally different
                //           assignment of grid cells to procs and induce
                //           complete rebalance data migration
                rcbflip = 0;
                if (narg == 3) rcbflip = int.Parse(arg[2]);

            }
            else sparta.error.all( "Illegal balance_grid command");


            // error check on methods only allowed for a uniform grid

            if (bstyle == (int)Enum1.STRIDE || bstyle == (int)Enum1.CLUMP || bstyle == (int)Enum1.BLOCK)
                if (sparta.grid.uniform==0)
                    sparta.error.all( "Invalid balance_grid style for non-uniform grid");

            // re-assign each of my local child cells to a proc
            // only assign unsplit and split cells
            // do not assign sub-cells since they migrate with their split cell
            // set nmigrate = # of cells that will migrate to a new proc
            // reset proc field in cells for migrating cells
            // style NONE performs no re-assignment

            sparta.mpi.MPI_Barrier(sparta.world);
            double time1 = sparta.mpi.MPI_Wtime();

            Grid.ChildCell[] cells = sparta.grid.cells;
            Grid.ChildInfo[] cinfo = sparta.grid.cinfo;
            int nglocal = sparta.grid.nlocal;

            int nprocs = sparta.comm.nprocs;
            int newproc;
            int nmigrate = 0;

            if (bstyle ==(int)Enum1.STRIDE)
            {
                cellint idm1, ix, iy, iz, nth=0;

                cellint nx = sparta.grid.unx;
                cellint ny = sparta.grid.uny;
                cellint nz = sparta.grid.unz;

                for (int icell = 0; icell < nglocal; icell++)
                {
                    if (cells[icell].nsplit <= 0) continue;
                    idm1 = cells[icell].id - 1;
                    ix = idm1 % nx;
                    iy = (idm1 / nx) % ny;
                    iz = idm1 / (nx * ny);

                    if (order == (int)Enum2.XYZ) nth = iz * nx * ny + iy * nx + ix;
                    else if (order == (int)Enum2.XZY) nth = iy * nx * nz + iz * nx + ix;
                    else if (order == (int)Enum2.YXZ) nth = iz * ny * nx + ix * ny + iy;
                    else if (order == (int)Enum2.YZX) nth = ix * ny * nz + iz * ny + iy;
                    else if (order == (int)Enum2.ZXY) nth = iy * nz * nx + ix * nz + iz;
                    else if (order == (int)Enum2.ZYX) nth = ix * nz * ny + iy * nz + iz;

                    newproc = nth % nprocs;

                    if (newproc != cells[icell].proc) nmigrate++;
                    cells[icell].proc = newproc;
                }

            }
            else if (bstyle == (int)Enum1.BLOCK)
            {
                cellint idm1, ix, iy, iz;
                int ipx, ipy, ipz;

                int nx = sparta.grid.unx;
                int ny = sparta.grid.uny;
                int nz = sparta.grid.unz;

                procs2grid(nx, ny, nz,ref px,ref py,ref pz);
                if (px * py * pz != nprocs)
                    sparta.error.all("Bad grid of processors for balance_grid block");

                for (int icell = 0; icell < nglocal; icell++)
                {
                    if (cells[icell].nsplit <= 0) continue;
                    idm1 = cells[icell].id - 1;
                    ix = idm1 % nx;
                    iy = (idm1 / nx) % ny;
                    iz = idm1 / (nx * ny);
                    ipx = ix * px / nx;
                    ipy = iy * py / ny;
                    ipz = iz * pz / nz;

                    newproc = ipz * px * py + ipy * px + ipx;

                    if (newproc != cells[icell].proc) nmigrate++;
                    cells[icell].proc = newproc;
                }

            }
            else if (bstyle == (int)Enum1.RANDOM)
            {
                
                RanPark random = new RanPark(sparta.update.ranmaster.uniform());
                double seed = sparta.update.ranmaster.uniform();
                random.reset(seed, sparta.comm.me, 100);

                for (int icell = 0; icell < nglocal; icell++)
                {
                    if (cells[icell].nsplit <= 0) continue;
                    newproc =Convert.ToInt32( nprocs * random.uniform());
                    if (newproc != cells[icell].proc) nmigrate++;
                    cells[icell].proc = newproc;
                }

                //delete random;

            }
            else if (bstyle==(int)Enum1.BISECTION)
            {
                RCB rcb = new RCB(sparta);

                double[,] x=new double[nglocal,3];

                double[] lo,hi;

                int nbalance = 0;

                for (int icell = 0; icell < nglocal; icell++)
                {
                    if (cells[icell].nsplit <= 0) continue;
                    lo = cells[icell].lo;
                    hi = cells[icell].hi;
                    x[nbalance,0] = 0.5 * (lo[0] + hi[0]);
                    x[nbalance,1] = 0.5 * (lo[1] + hi[1]);
                    x[nbalance,2] = 0.5 * (lo[2] + hi[2]);
                    nbalance++;
                }

                double[] wt;

                if (rcbwt==(int)Enum3.PARTICLE)
                {
                    sparta.particle.sort();
                    int zero = 0;
                    int n;
                    wt = new double[nglocal];
                    nbalance = 0;
                    for (int icell = 0; icell < nglocal; icell++)
                    {
                        if (cells[icell].nsplit <= 0) continue;
                        n = cinfo[icell].count;
                        if (n!=0) wt[nbalance++] = n;
                        else
                        {
                            wt[nbalance++] = ZEROPARTICLE;
                            zero++;
                        }
                    }
                }
                else
                {
                    wt = new double[nglocal];
                }

                rcb.compute(nbalance, x, wt, rcbflip);
                rcb.invert();

                nbalance = 0;
                int[] sendproc = rcb.sendproc;
                for (int icell = 0; icell < nglocal; icell++)
                {
                    if (cells[icell].nsplit <= 0) continue;
                    cells[icell].proc = sendproc[nbalance++];
                }
                nmigrate = nbalance - rcb.nkeep;

                
            }

            // set clumped of not, depending on style
            // NONE style does not change clumping

            if (nprocs == 1 || bstyle == (int)Enum1.CLUMP || bstyle == (int)Enum1.BLOCK || bstyle == (int)Enum1.BISECTION)
                sparta.grid.clumped = 1;
            else if (bstyle != (int)Enum1.NONE) sparta.grid.clumped = 0;

            sparta.mpi.MPI_Barrier(sparta.world);
            double time2 = sparta.mpi.MPI_Wtime();

            // sort particles
            // NOTE: not needed again if rcbwt = PARTICLE for bstyle = BISECTION ??

            sparta.particle.sort();
            sparta.mpi.MPI_Barrier(sparta.world);
            double time3 = sparta.mpi.MPI_Wtime();


            // invoke init() so all grid cell info, including collide & fixes,
            //   is ready to migrate
            // for init, do not require surfs be assigned collision models
            //   this allows balance call early in script, e.g. from ReadRestart
            // migrate grid cells and their particles to new owners
            // invoke grid methods to complete grid setup

            int ghost_previous = sparta.grid.exist_ghost;

            sparta.domain.boundary_collision_check = 0;
            sparta.surf.surf_collision_check = 0;
            sparta.init();
            sparta.domain.boundary_collision_check = 1;
            sparta.surf.surf_collision_check = 1;

            sparta.grid.unset_neighbors();
            sparta.grid.remove_ghosts();
            sparta.comm.migrate_cells(nmigrate); 

            sparta.mpi.MPI_Barrier(sparta.world);
            double time4 = sparta.mpi.MPI_Wtime();

            sparta.grid.setup_owned();
            sparta.grid.acquire_ghosts();
            if (ghost_previous!=0)
            {
               sparta.grid.reset_neighbors(); 
            }
            else sparta.grid.find_neighbors();
            sparta.comm.reset_neighbors();


            // reallocate per grid arrays in per grid dumps

            for (int i = 0; i < sparta.output.ndump; i++)
                sparta.output.dump[i].reset_grid();

            sparta.mpi.MPI_Barrier(sparta.world);
            double time5 = sparta.mpi.MPI_Wtime();

            // stats on balance operation
            // only print if outflag = 1
            // some callers suppress output, e.g. ReadRestart

            bigint count = nmigrate;
            bigint nmigrate_all=0;
            sparta.mpi.MPI_Allreduce(ref count, ref nmigrate_all, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            double time_total = time5 - time1;

            if (sparta.comm.me == 0 && outflag != 0)
            {
                string str1 = string.Format("Balance grid migrated {0} cells\n", nmigrate_all);
                string str2 = string.Format("  CPU time = {0:G} secs\n", time_total);
                string str3 = string.Format("  reassign/sort/migrate/ghost percent = {0:G6} {1:G6} {2:G6} {3:G6}\n",
                            100.0 * (time2 - time1) / time_total, 100.0 * (time3 - time2) / time_total,
                            100.0 * (time4 - time3) / time_total, 100.0 * (time5 - time4) / time_total);

                Console.WriteLine(str1 + str2 + str3);
                if (sparta.screen != null)
                {
                    new StreamWriter(sparta.screen).WriteLine(str1 + str2 + str3);
                    //fprintf(screen, "Balance grid migrated " BIGINT_FORMAT " cells\n",
                    //        nmigrate_all);
                    //fprintf(screen, "  CPU time = %g secs\n", time_total);
                    //fprintf(screen, "  reassign/sort/migrate/ghost percent = %g %g %g %g\n",
                    //        100.0 * (time2 - time1) / time_total, 100.0 * (time3 - time2) / time_total,
                    //        100.0 * (time4 - time3) / time_total, 100.0 * (time5 - time4) / time_total);
                }
                if (sparta.logfile!=null)
                {
                    new StreamWriter(sparta.logfile).WriteLine(str1 + str2 + str3);
                }

            }

        }


        private void procs2grid(int nx, int ny, int nz,ref int px,ref int py,ref int pz)
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
                if (nprocs % ipx != 0) valid = 0;
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
                    if (sparta.domain.dimension == 2 && ipz != 1) valid = 0;
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
