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

            Console.WriteLine("gridcomm.unpack_one"+buf);
            return 16;
        }
        //      public int pack_one_adapt(char*, char*, int);
        public int pack_particles(int icell,ref StringBuilder buf, int memflag)
        {
            int n = 0;
            int np = cinfo[icell].count;
            if (memflag!=0)
            {
                //*((int*)ptr) = np;
            }

            n += sizeof(int);

            if (np==0) return n;

            if (memflag!=0)
            {
                Particle.OnePart[] particles = sparta.particle.particles;
                int[] next = sparta.particle.next;
                int ip = cinfo[icell].first;

                while (ip >= 0)
                {
                    //memcpy(ptr, &particles[ip], nbytes_particle);
                    buf.Append(particles[ip]);
                    n += nbytes_particle;
                    if (ncustom!=0)
                    {
                        sparta.particle.pack_custom(ip,ref buf);
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
        //      public void compress();

        void unpack_ghosts(int nsize, StringBuilder buf)
        {
            int n = 0;
            while (n < nsize)
            {
                n += unpack_one(buf, 0, 0);
            }
        }

        int pack_one(int icell,int ownflag, int molflag, int memflag,ref StringBuilder buf)
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
