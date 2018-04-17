using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public partial class Grid
    {
        //      public int pack_one(int, char*, int, int, int);
        public int unpack_one(StringBuilder buf, int ownflag, int molflag)
        {
            //string str = buf.ToString();
            //int icell;
            //int ptr = 0;
            //if (ownflag!=0) icell = nlocal;
            //else icell = nlocal + nghost;
            //grow_cells(1, ownflag);
            //if (ownflag!=0) nlocal++;
            //else nghost++;

            Console.WriteLine("gridcomm.unpack_one" + buf);
            return 16;
        }
        //      public int pack_one_adapt(char*, char*, int);
        public int pack_particles(int icell, ref StringBuilder buf, int memflag)
        {
            int n = 0;
            int np = cinfo[icell].count;
            if (memflag != 0)
            {
                //*((int*)ptr) = np;
            }

            n += sizeof(int);

            if (np == 0) return n;

            if (memflag != 0)
            {
                Particle.OnePart[] particles = sparta.particle.particles;
                int[] next = sparta.particle.next;
                int ip = cinfo[icell].first;

                while (ip >= 0)
                {
                    //memcpy(ptr, &particles[ip], nbytes_particle);
                    buf.Append(particles[ip]);
                    n += nbytes_particle;
                    if (ncustom != 0)
                    {
                        sparta.particle.pack_custom(ip, ref buf);
                        n += nbytes_custom;
                    }
                    particles[ip].icell = -1;
                    ip = next[ip];
                }
            }
            else n += np * nbytes_total;

            return n;
        }
        //      public int unpack_particles(char*, int);
        //      public void unpack_particles_adapt(int, char*);
        public void compress()
        {
            // must compress per-cell arrays in collide and fixes before
            // cells data structure changes

            if (sparta.collide != null) sparta.collide.compress_grid();
            if (sparta.modify.n_pergrid != 0) sparta.modify.compress_grid(0);

            // copy of integer lists
            // create new lists

            MyPage<int> csurfs_old = csurfs;
            MyPage<int> csplits_old = csplits;
            MyPage<int> csubs_old = csubs;

            csurfs = null; csplits = null; csubs = null;
            allocate_surf_arrays();

            // compress cells and cinfo and sinfo arrays
            // 4 cases:
            // discard unsplit or sub cell:
            //   continue
            // keep unsplit or sub cell:
            //   memcpy of cells and cinfo, setup csurfs
            //   reset cells.ilocal and sinfo.csubs
            //   increment nlocal, nunsplitlocal, nsublocal
            // discard split cell:
            //   flag all sub cells with other proc
            //   continue
            // keep split cell:
            //   memcpy of cells and cinfo and sinfo, setup csurfs/csplits/csubs
            //   reset sinfo.icell
            //   reset cells.isplit, ditto for sub cells
            //   increment nlocal, nsplitlocal
            // relies on sub cells appearing in cells array after their split cell
            //   no need for sub cells of one split cell to be contiguous before or after

            int ncurrent = nlocal;
            nlocal = nunsplitlocal = nsplitlocal = nsublocal = 0;

            for (int icell = 0; icell < ncurrent; icell++)
            {

                // unsplit and sub cells

                if (cells[icell].nsplit <= 1)
                {
                    if (cells[icell].proc != me) continue;

                    if (icell != nlocal)
                    {
                        cells[nlocal] = cells[icell];
                        cinfo[nlocal] = cinfo[icell];
                        //memcpy(&cells[nlocal], &cells[icell], sizeof(ChildCell));
                        //memcpy(&cinfo[nlocal], &cinfo[icell], sizeof(ChildInfo));
                    }

                    cells[nlocal].ilocal = nlocal;

                    // for unsplit cell, new copy of all csurfs indices
                    // for sub cell, just copy csurfs ptr from its split cell

                    if (cells[nlocal].nsurf != 0)
                    {
                        if (cells[nlocal].nsplit == 1)
                        {
                            int[] oldcsurfs = cells[icell].csurfs;
                            cells[nlocal].csurfs = csurfs.vget();
                            Array.Copy(oldcsurfs, cells[nlocal].csurfs, cells[nlocal].nsurf);
                            //memcpy(cells[nlocal].csurfs, oldcsurfs,
                            //       cells[nlocal].nsurf * sizeof(int));
                            csurfs.vgot(cells[nlocal].nsurf);
                        }
                        else
                            cells[nlocal].csurfs =
                              cells[sinfo[cells[nlocal].isplit].icell].csurfs;
                    }

                    if (cells[nlocal].nsplit <= 0)
                    {
                        sinfo[cells[nlocal].isplit].csubs[-cells[nlocal].nsplit] = nlocal;
                        nsublocal++;
                    }
                    else nunsplitlocal++;

                    nlocal++;

                    // split cells

                }
                else
                {
                    if (cells[icell].proc != me)
                    {
                        int isplits = cells[icell].isplit;
                        int nsplits = cells[icell].nsplit;
                        for (int i = 0; i < nsplits; i++)
                        {
                            int m = sinfo[isplits].csubs[i];
                            cells[m].proc = cells[icell].proc;
                        }
                        continue;
                    }

                    if (icell != nlocal)
                    {
                        cells[nlocal] = cells[icell];
                        cinfo[nlocal] = cinfo[icell];
                        // memcpy(&cells[nlocal], &cells[icell], sizeof(ChildCell));
                        // memcpy(&cinfo[nlocal], &cinfo[icell], sizeof(ChildInfo));
                    }

                    cells[nlocal].ilocal = nlocal;

                    // new copy of all csurfs indices

                    int[] oldcsurfs = cells[icell].csurfs;
                    cells[nlocal].csurfs = csurfs.vget();
                    Array.Copy(oldcsurfs, cells[nlocal].csurfs, cells[nlocal].nsurf);
                    //memcpy(cells[nlocal].csurfs, oldcsurfs,
                    //cells[nlocal].nsurf * sizeof(int));
                    csurfs.vgot(cells[nlocal].nsurf);

                    // compress sinfo

                    int isplit = cells[nlocal].isplit;
                    if (isplit != nsplitlocal)
                    {
                        sinfo[nsplitlocal] = sinfo[isplit];
                        //memcpy(&sinfo[nsplitlocal], &sinfo[isplit], sizeof(SplitInfo));
                    }

                    cells[nlocal].isplit = nsplitlocal;
                    sinfo[nsplitlocal].icell = nlocal;

                    // new copy of all csplits indices

                    int[] oldcsplits = sinfo[isplit].csplits;
                    sinfo[nsplitlocal].csplits = csplits.vget();
                    Array.Copy(oldcsplits, sinfo[nsplitlocal].csplits, cells[nlocal].nsurf);
                    //memcpy(sinfo[nsplitlocal].csplits, oldcsplits,
                    //               cells[nlocal].nsurf * sizeof(int));
                    csplits.vgot(cells[nlocal].nsurf);

                    // new csubs list of length nsplit
                    // values unset for now, set one-by-one when sub cells are compressed

                    int[] oldcsubs = sinfo[isplit].csubs;
                    sinfo[nsplitlocal].csubs = csubs.vget();
                    csubs.vgot(cells[nlocal].nsplit);

                    // point each sub cell in old sinfo csubs at new compressed sinfo index

                    int nsplit = cells[nlocal].nsplit;
                    for (int i = 0; i < nsplit; i++)
                        cells[oldcsubs[i]].isplit = nsplitlocal;

                    nsplitlocal++;
                    nlocal++;
                }
            }
            hashfilled = 0;

            // delete old integer lists

            //delete csurfs_old;
            //delete csplits_old;
            //delete csubs_old;

            // repoint particles in all remaining grid cells to new icell indices
            // assumes particles are sorted and have not yet been compressed,
            //   so count/first values in compressed cinfo data struct are still valid
            // when done, particles are still sorted

            Particle.OnePart[] particles = sparta.particle.particles;
            int[] next = sparta.particle.next;

            int ip;
            for (int icell = 0; icell < nlocal; icell++)
            {
                ip = cinfo[icell].first;
                while (ip >= 0)
                {
                    particles[ip].icell = icell;
                    ip = next[ip];
                }
            }

            // some fixes have post-compress operations to perform

            if (sparta.modify.n_pergrid!=0) sparta.modify.compress_grid(1);
        }

        void unpack_ghosts(int nsize, StringBuilder buf)
        {
            int n = 0;
            while (n < nsize)
            {
                n += unpack_one(buf, 0, 0);
            }
        }
              
        public int pack_one(int icell,int ownflag, int molflag, int memflag,ref StringBuilder buf)
        {
            int n = 0;
            

            if (memflag != 0)
            {
                buf.Append(cells[icell]);
            }
            n +=Marshal.SizeOf(cells[icell]);

            if (cells[icell].nsurf < 0)
            {
                return n;
            }

            if (cells[icell].nsurf == 0)
            {
                int nsurf = cells[icell].nsurf;
                if (memflag != 0)
                {
                    buf.Append(cells[icell].csurfs);
                }
                n += nsurf * sizeof(int);
            }

            if (ownflag!=0)
            {
                if (memflag != 0)
                {
                    buf.Append(cinfo[icell]);
                }
                n += Marshal.SizeOf(cinfo[icell]);
            }

            if (cells[icell].nsplit>1)
            {
                int isplit = cells[icell].isplit;
                if (memflag!=0)
                {
                    buf.Append(sinfo[isplit]);
                }
                n += Marshal.SizeOf(sinfo[isplit]);

                int nsurf = cells[icell].nsurf;
                if (memflag != 0)
                {
                    buf.Append(sinfo[isplit].csplits);
                }
                n += nsurf*sizeof(int);

                if (ownflag!=0)
                {
                    int[] csubs = sinfo[isplit].csubs;
                                                             
                    double[] dptr = new double[nsplit];
                    if (memflag != 0)
                    {
                        for (int i = 0; i < nsplit; i++)
                        {
                            dptr[i] = cinfo[csubs[i]].volume;
                        }
                        n += nsplit * sizeof(double);
                    }
                }

            }

            if (ownflag!=0)
            {
                if (sparta.collide!=null)
                {
                    n += sparta.collide.pack_grid_one(icell,ref buf, memflag);
                }
                if (sparta.modify.n_pergrid!=0)
                {
                    n += sparta.modify.pack_grid_one(icell, ref buf, memflag);
                }
            }

            if (molflag==0)
            {
                return n;
            }

            n += pack_particles(icell, ref buf, memflag);

            if (cells[icell].nsplit>1)
            {
                int isplit = cells[icell].isplit;
                int nsplit = cells[icell].nsplit;
                for (int i = 0; i < nsplit; i++)
                {
                    int m = sinfo[isplit].csubs[i];
                    n += pack_particles(m,ref buf, memflag);
                }



            }

            return n;

        }


    }
}
