using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public partial class Grid
    {
        public void surf2grid(int subflag, int outflag)
        {
            int i, isub, nsurf, nsplitone, xsub;
            int[] surfmap,ptr;
            double[] lo,hi,vols;
            double[] xsplit=new double[3];
            ChildCell c;
            SplitInfo s;
            Cut2d cut2d;
            Cut3d cut3d;

            int dim = sparta.domain.dimension;

            double[] slo = sparta.surf.bblo;
            double[] shi = sparta.surf.bbhi;

            cut3d = new Cut3d(sparta);
            cut2d = new Cut2d(sparta, sparta.domain.axisymmetric);

            // compute overlap of surfs with each cell I own
            // info stored in nsurf,csurfs

            //double t1 = MPI_Wtime();

            for (int icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;

                lo = cells[icell].lo;
                hi = cells[icell].hi;
                if (box_overlap(lo, hi, slo, shi)==0) continue;

                ptr = csurfs.vget();

                if (dim == 3)
                    nsurf = cut3d.surf2grid(cells[icell].id, cells[icell].lo, cells[icell].hi,
                                             ptr, maxsurfpercell);
                else
                    nsurf = cut2d.surf2grid(cells[icell].id, cells[icell].lo, cells[icell].hi,
                                             ptr, maxsurfpercell);

                if (nsurf < 0) sparta.error.one("Too many surfs in one cell");
                if (nsurf!=0)
                {
                    cinfo[icell].type = (int)Enum4.OVERLAP;
                    cells[icell].nsurf = nsurf;
                    cells[icell].csurfs = ptr;
                    csurfs.vgot(nsurf);
                }
            }

            //double t2 = MPI_Wtime();
            //printf("TIME %g\n",t2-t1);

            if (outflag!=0) surf2grid_stats();

            // compute cut volume and possible split of each grid cell by surfs
            // decrement nunsplitlocal if convert an unsplit cell to split cell
            // if nsplitone > 1, create new split cell sinfo and sub-cells

            int ncurrent = nlocal;
            for (int icell = 0; icell < ncurrent; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                if (cinfo[icell].type != OVERLAP) continue;

                surfmap = csplits.vget();
                c = cells[icell];

                if (dim == 3)
                    nsplitone = cut3d.split(c.id, c.lo, c.hi, c.nsurf, c.csurfs,
                                             vols, surfmap, cinfo[icell].corner, xsub, xsplit);
                else
                    nsplitone = cut2d.split(c.id, c.lo, c.hi, c.nsurf, c.csurfs,
                                             vols, surfmap, cinfo[icell].corner, xsub, xsplit);
                cells[icell] = c;

                if (nsplitone == 1)
                {
                    cinfo[icell].volume = vols[0];

                }
                else if (subflag!=0)
                {
                    cells[icell].nsplit = nsplitone;
                    nunsplitlocal--;

                    cells[icell].isplit = nsplitlocal;
                    add_split_cell(1);
                    s = sinfo[nsplitlocal - 1];
                    s.icell = icell;
                    s.csplits = surfmap;
                    s.xsub = xsub;
                    s.xsplit[0] = xsplit[0];
                    s.xsplit[1] = xsplit[1];
                    if (dim == 3) s.xsplit[2] = xsplit[2];
                    else s.xsplit[2] = 0.0;

                    sinfo[nsplitlocal - 1] = s;
                    ptr = s.csubs = csubs.vget();

                    for (i = 0; i < nsplitone; i++)
                    {
                        isub = nlocal;
                        add_sub_cell(icell, 1);
                        cells[isub].nsplit = -i;
                        cinfo[isub].volume = vols[i];
                        ptr[i] = isub;
                    }

                    csplits.vgot(cells[icell].nsurf);
                    csubs.vgot(nsplitone);

                }
                else
                {
                    if (cells[icell].nsplit != nsplitone)
                    {
                        Console.WriteLine("BAD {0} {1}: {2} {3}\n", icell, cells[icell].id,
                               nsplitone, cells[icell].nsplit);
                        sparta.error.one("Inconsistent surface to grid mapping in read_restart");
                    }

                    s = sinfo[cells[icell].isplit];
                    s.csplits = surfmap;
                    s.xsub = xsub;
                    s.xsplit[0] = xsplit[0];
                    s.xsplit[1] = xsplit[1];
                    if (dim == 3) s.xsplit[2] = xsplit[2];
                    else s.xsplit[2] = 0.0;
                    sinfo[cells[icell].isplit] = s;

                    ptr = s.csubs;
                    for (i = 0; i < nsplitone; i++)
                    {
                        isub = ptr[i];
                        cells[isub].nsurf = cells[icell].nsurf;
                        cells[isub].csurfs = cells[icell].csurfs;
                        cinfo[isub].volume = vols[i];
                    }

                    csplits.vgot(cells[icell].nsurf);
                }
            }

            //double t3 = MPI_Wtime();
            //printf("TIME %g\n",t3-t2);

            // stats on pushed cells and unmarked corner points in OVERLAP cells

            if (outflag!=0)
            {
                int npushmax;
                int[] npushcell;
                if (dim == 3)
                {
                    npushmax = cut3d.npushmax;
                    npushcell = cut3d.npushcell;
                }
                else
                {
                    npushmax = cut2d.npushmax;
                    npushcell = cut2d.npushcell;
                }
                int[] npushall = new int[npushmax + 1];
                //sparta.mpi.MPI_Allreduce(ref npushcell, ref npushall, npushmax + 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
                if (sparta.comm.me == 0)
                {
                    if (sparta.screen!=null)
                    {
                        StreamWriter sw = new StreamWriter(sparta.screen);
                        
                        sw.WriteLine( "  ");
                        Console.WriteLine("  ");
                        for (int ii = 1; ii <= npushmax; ii++)
                        {
                            sw.WriteLine( "{0} ", npushall[ii]);
                            Console.WriteLine("{0} ", npushall[ii]);
                        }

                        sw.WriteLine( "= number of pushed cells\n");
                        Console.WriteLine("= number of pushed cells\n");
                    }
                    if (sparta.logfile!=null)
                    {
                        StreamWriter sw = new StreamWriter(sparta.logfile);
                        sw.WriteLine( "  ");
                        for (int ii = 1; ii <= npushmax; ii++)
                            sw.WriteLine( "{0} ", npushall[ii]);
                        sw.WriteLine( "= number of pushed cells\n");
                    }
                }
                //delete[] npushall;

                int noverlap = 0;
                int ncorner = 0;
                for (int icell = 0; icell < nlocal; icell++)
                {
                    if (cells[icell].nsplit <= 0) continue;
                    if (cinfo[icell].type == (int)Enum4.OVERLAP)
                    {
                        noverlap++;
                        if (cinfo[icell].corner[0] == (int)Enum4.UNKNOWN) ncorner++;
                    }
                }
                int ncornerall=0, noverlapall=0;
                sparta.mpi.MPI_Allreduce(ref ncorner, ref ncornerall, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
                sparta.mpi.MPI_Allreduce(ref noverlap, ref noverlapall, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
                if (sparta.comm.me == 0)
                {
                    string str = string.Format("  {0} {1} = cells overlapping surfs,overlap cells with unmarked corner pts\n",
                                        noverlapall, ncornerall);
                    Console.WriteLine(str);
                    if (sparta.screen!=null)
                    {
                        new StreamWriter(sparta.screen).WriteLine(str);
                    }

                    if (sparta.logfile!=null)
                    {
                        new StreamWriter(sparta.logfile).WriteLine(str);
                    }
                }
            }

            // clean up

        }
        //      public void surf2grid_one(int, int, int, int, class Cut3d *, class Cut2d *);
        public void clear_surf()
        {
            int dimension = sparta.domain.dimension;
            int ncorner = 8;
            if (dimension == 2) ncorner = 4;
            double[] lo,hi;

            // if surfs no longer exist, set cell type to OUTSIDE, else UNKNOWN
            // set corner points of every cell to UNKNOWN

            int celltype = (int)Enum4.UNKNOWN;
            if (sparta.surf.exist==0) celltype = (int)Enum4.OUTSIDE;

            int nlocal_prev = nlocal;

            int icell = 0;
            while (icell < nlocal)
            {
                if (cells[icell].nsplit <= 0)
                {
                    cells[icell] = cells[nlocal - 1];
                    cinfo[icell] = cinfo[nlocal - 1];
                    //memcpy(&cells[icell], &cells[nlocal - 1], sizeof(ChildCell));
                    //memcpy(&cinfo[icell], &cinfo[nlocal - 1], sizeof(ChildInfo));
                    nlocal--;
                }
                else
                {
                    cells[icell].ilocal = icell;
                    cells[icell].nsurf = 0;
                    cells[icell].csurfs = null;
                    cells[icell].nsplit = 1;
                    cells[icell].isplit = -1;
                    cinfo[icell].type = celltype;
                    for (int m = 0; m < ncorner; m++) cinfo[icell].corner[m] = (int)Enum4.UNKNOWN;
                    lo = cells[icell].lo;
                    hi = cells[icell].hi;
                    if (dimension == 3)
                        cinfo[icell].volume = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
                    else if (sparta.domain.axisymmetric!=0)
                        cinfo[icell].volume = MyConst.MY_PI * (hi[1] * hi[1] - lo[1] * lo[1]) * (hi[0] - lo[0]);
                    else
                        cinfo[icell].volume = (hi[0] - lo[0]) * (hi[1] - lo[1]);
                    icell++;
                }
            }

            // if particles exist and local cell count changed
            // repoint particles to new icell indices
            // assumes particles are sorted,
            //   so count/first values in compressed cinfo data struct are still valid
            // when done, particles are still sorted

            if (sparta.particle.exist!=0 && nlocal < nlocal_prev)
            {
                Particle.OnePart[] particles = sparta.particle.particles;
                int[] next = sparta.particle.next;

                int ip;
                for (int aa = 0; aa < nlocal; aa++)
                {
                    ip = cinfo[aa].first;
                    while (ip >= 0)
                    {
                        particles[ip].icell = aa;
                        ip = next[ip];
                    }
                }
            }

            // reset csurfs and csplits and csubs so can refill

            csurfs .reset();
            csplits.reset();
            csubs  .reset();

            // reset all cell counters

            nunsplitlocal = nlocal;
            nsplitlocal = nsublocal = 0;
        }
        //      public void clear_surf_restart();
        public void combine_split_cell_particles(int icell, int relabel)
        {
            int ip, iplast=0, jcell;

            int nsplit = cells[icell].nsplit;
            int[] mycsubs = sinfo[cells[icell].isplit].csubs;
            int count = 0;
            int first = -1;

            int[] next = sparta.particle.next;

            for (int i = 0; i < nsplit; i++)
            {
                jcell = mycsubs[i];
                count += cinfo[jcell].count;
                if (cinfo[jcell].first < 0) continue;

                if (first < 0) first = cinfo[jcell].first;
                else next[iplast] = cinfo[jcell].first;

                ip = cinfo[jcell].first;
                while (ip >= 0)
                {
                    iplast = ip;
                    ip = next[ip];
                }
            }

            cinfo[icell].count = count;
            cinfo[icell].first = first;

            // repoint each particle now in parent split cell to the split cell

            if (relabel!=0)
            {
                Particle.OnePart[] particles = sparta.particle.particles;
                ip = first;
                while (ip >= 0)
                {
                    particles[ip].icell = icell;
                    ip = next[ip];
                }
            }
        }
        //      public void assign_split_cell_particles(int);
        public void allocate_surf_arrays()
        {
            csurfs = new MyPage<int>(maxsurfpercell, Math.Max(100 * maxsurfpercell, 1024));
            csplits = new MyPage<int>(maxsurfpercell, Math.Max(100 * maxsurfpercell, 1024));
            csubs = new MyPage<int>(MAXSPLITPERCELL, 128);
        }
        //      int* csubs_request(int);

        void surf2grid_stats(int subflag, int outflag)
        {
            int i, isub, nsurf, nsplitone, xsub;
            int[] surfmap,ptr;
            double[] lo,hi,vols;
            double[] xsplit=new double[3];
            ChildCell c;
            SplitInfo s;
            Cut2d cut2d;
            Cut3d cut3d;

            int dim = sparta.domain.dimension;

            double[] slo = sparta.surf.bblo;
            double[] shi = sparta.surf.bbhi;

            cut3d = new Cut3d(sparta);
            cut2d = new Cut2d(sparta, sparta.domain.axisymmetric);

            // compute overlap of surfs with each cell I own
            // info stored in nsurf,csurfs

            //double t1 = MPI_Wtime();

            for (int icell = 0; icell < nlocal; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;

                lo = cells[icell].lo;
                hi = cells[icell].hi;
                if (box_overlap(lo, hi, slo, shi)==0) continue;

                ptr = csurfs.vget();

                if (dim == 3)
                    nsurf = cut3d.surf2grid(cells[icell].id, cells[icell].lo, cells[icell].hi,
                                             ptr, maxsurfpercell);
                else
                    nsurf = cut2d.surf2grid(cells[icell].id, cells[icell].lo, cells[icell].hi,
                                             ptr, maxsurfpercell);

                if (nsurf < 0) sparta.error.one( "Too many surfs in one cell");
                if (nsurf!=0)
                {
                    cinfo[icell].type = (int)Enum4.OVERLAP;
                    cells[icell].nsurf = nsurf;
                    cells[icell].csurfs = ptr;
                    csurfs.vgot(nsurf);
                }
            }

            //double t2 = MPI_Wtime();
            //printf("TIME %g\n",t2-t1);

            if (outflag!=0) surf2grid_stats();

            // compute cut volume and possible split of each grid cell by surfs
            // decrement nunsplitlocal if convert an unsplit cell to split cell
            // if nsplitone > 1, create new split cell sinfo and sub-cells

            int ncurrent = nlocal;
            for (int icell = 0; icell < ncurrent; icell++)
            {
                if (cells[icell].nsplit <= 0) continue;
                if (cinfo[icell].type != (int)Enum4.OVERLAP) continue;

                surfmap = csplits.vget();
                c = cells[icell];

                if (dim == 3)
                    nsplitone = cut3d.split(c.id, c.lo, c.hi, c.nsurf, c.csurfs,
                                             vols, surfmap, cinfo[icell].corner, xsub, xsplit);
                else
                    nsplitone = cut2d.split(c.id, c.lo, c.hi, c.nsurf, c.csurfs,
                                             vols, surfmap, cinfo[icell].corner, xsub, xsplit);

                if (nsplitone == 1)
                {
                    cinfo[icell].volume = vols[0];

                }
                else if (subflag)
                {
                    cells[icell].nsplit = nsplitone;
                    nunsplitlocal--;

                    cells[icell].isplit = nsplitlocal;
                    add_split_cell(1);
                    s = &sinfo[nsplitlocal - 1];
                    s.icell = icell;
                    s.csplits = surfmap;
                    s.xsub = xsub;
                    s.xsplit[0] = xsplit[0];
                    s.xsplit[1] = xsplit[1];
                    if (dim == 3) s.xsplit[2] = xsplit[2];
                    else s.xsplit[2] = 0.0;

                    ptr = s.csubs = csubs.vget();

                    for (i = 0; i < nsplitone; i++)
                    {
                        isub = nlocal;
                        add_sub_cell(icell, 1);
                        cells[isub].nsplit = -i;
                        cinfo[isub].volume = vols[i];
                        ptr[i] = isub;
                    }

                    csplits.vgot(cells[icell].nsurf);
                    csubs.vgot(nsplitone);

                }
                else
                {
                    if (cells[icell].nsplit != nsplitone)
                    {
                        printf("BAD %d %d: %d %d\n", icell, cells[icell].id,
                               nsplitone, cells[icell].nsplit);
                        sparta.error.one(
                                   "Inconsistent surface to grid mapping in read_restart");
                    }

                    s = &sinfo[cells[icell].isplit];
                    s.csplits = surfmap;
                    s.xsub = xsub;
                    s.xsplit[0] = xsplit[0];
                    s.xsplit[1] = xsplit[1];
                    if (dim == 3) s.xsplit[2] = xsplit[2];
                    else s.xsplit[2] = 0.0;

                    ptr = s.csubs;
                    for (i = 0; i < nsplitone; i++)
                    {
                        isub = ptr[i];
                        cells[isub].nsurf = cells[icell].nsurf;
                        cells[isub].csurfs = cells[icell].csurfs;
                        cinfo[isub].volume = vols[i];
                    }

                    csplits.vgot(cells[icell].nsurf);
                }
            }

            //double t3 = MPI_Wtime();
            //printf("TIME %g\n",t3-t2);

            // stats on pushed cells and unmarked corner points in OVERLAP cells

            if (outflag!=0)
            {
                int npushmax;
                int[] npushcell;
                if (dim == 3)
                {
                    npushmax = cut3d.npushmax;
                    npushcell = cut3d.npushcell;
                }
                else
                {
                    npushmax = cut2d.npushmax;
                    npushcell = cut2d.npushcell;
                }
                int[] npushall = new int[npushmax + 1];
                MPI_Allreduce(npushcell, npushall, npushmax + 1, MPI_INT, MPI_SUM, world);
                if (sparta.comm.me == 0)
                {
                    if (screen)
                    {
                        fprintf(screen, "  ");
                        for (int i = 1; i <= npushmax; i++)
                            fprintf(screen, "%d ", npushall[i]);
                        fprintf(screen, "= number of pushed cells\n");
                    }
                    if (logfile)
                    {
                        fprintf(logfile, "  ");
                        for (int i = 1; i <= npushmax; i++)
                            fprintf(logfile, "%d ", npushall[i]);
                        fprintf(logfile, "= number of pushed cells\n");
                    }
                }
                delete[] npushall;

                int noverlap = 0;
                int ncorner = 0;
                for (int icell = 0; icell < nlocal; icell++)
                {
                    if (cells[icell].nsplit <= 0) continue;
                    if (cinfo[icell].type == OVERLAP)
                    {
                        noverlap++;
                        if (cinfo[icell].corner[0] == UNKNOWN) ncorner++;
                    }
                }
                int ncornerall, noverlapall;
                MPI_Allreduce(&ncorner, &ncornerall, 1, MPI_INT, MPI_SUM, world);
                MPI_Allreduce(&noverlap, &noverlapall, 1, MPI_INT, MPI_SUM, world);
                if (comm.me == 0)
                {
                    if (screen) fprintf(screen, "  %d %d = cells overlapping surfs, "
              
                                        "overlap cells with unmarked corner pts\n",
                                        noverlapall, ncornerall);
                    if (logfile) fprintf(logfile, "  %d %d = cells overlapping surfs, "
              
                                         "overlap cells with unmarked corner pts\n",
                                         noverlapall, ncornerall);
                }
            }

            // clean up

            //if (dim == 3) delete cut3d;
            //else delete cut2d;
        }

    }
}
