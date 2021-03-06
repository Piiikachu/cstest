﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;
using bigint = System.Int64;
using cellint = System.Int32;
namespace cstest
{
    public class Particle
    {
        enum Enum1 { PKEEP, PINSERT, PDONE, PDISCARD, PENTRY, PEXIT, PSURF };  // several files
        enum Enum2 { NONE, DISCRETE, SMOOTH };            // several files
        enum Enum3 { INT, DOUBLE };                      // several files

        public const int DELTA = 16384;
        public const int DELTASPECIES = 16;
        public const int DELTAMIXTURE = 8;
        public const int MAXLINE = 1024;
        public const string AIR = "N O NO";


        public int exist;                // 1 if particles exist
        public int sorted;               // 1 if particles are sorted by grid cell

        public struct Species
        {          // info on each particle species
            public string id;
            public double molwt;
            public double mass;
            public double rotrel;
            public double vibrel;
            public double vibtemp;
            public double specwt;
            public double charge;
            public int rotdof, vibdof;
            public int internaldof;
        };

        public List<Species> species;         // list of particle species info
        public int nspecies;             // # of defined species

        public List<Mixture> mixture;
        public int nmixture;
        public int maxmixture;

        [StructLayout(LayoutKind.Sequential,Pack =1)]
        public struct OnePart
        {
            public int id;                 // particle ID
            public int ispecies;           // particle species index
            public int icell;              // which local Grid::cells the particle is in
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3, ArraySubType = UnmanagedType.SysInt)]
            public double[] x;            // particle position
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3, ArraySubType = UnmanagedType.SysInt)]
            public double[] v;            // particle velocity
            public double erot;            // rotational energy
            public double evib;            // vibrational energy
            public int flag;               // used for migration status
            public double dtremain;        // portion of move timestep remaining
            public double weight;          // particle or cell weight, if weighting enabled

            
        }
        [StructLayout(LayoutKind.Sequential, Pack = 1)]
        public struct OnePartRestart
        {
            public int id;                 // particle ID
            public int ispecies;           // particle species index
            public cellint icell;          // cell ID the particle is in
            public int nsplit;             // 1 for unsplit cell
                                           // else neg of sub cell index (0 to Nsplit-1)
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3, ArraySubType = UnmanagedType.SysInt)]
            public double[] x;            // particle position
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3, ArraySubType = UnmanagedType.SysInt)]
            public double[] v;            // particle velocity
            public double erot;            // rotational energy
            public double evib;            // vibrational energy
        }

        public bigint nglobal;           // global # of particles
        public int nlocal;               // # of particles I own
        public int maxlocal;             // max # particles list can hold
        public OnePart[] particles;       // list of particles I own

        // currently stored in grid.h for every cell, whether I own it or not
        // not sure why storing it here is slower

        //int *cellcount;           // count of particles in each grid cell I own
        //int *first;               // index of first particle in each grid cell

        public int[] next;                // index of next particle in each grid cell

        // extra custom vectors/arrays for per-particle data
        // ncustom > 0 if there are any extra arrays
        // custom attributes are created by various commands
        // these variables are public, others below are private

        public int ncustom;              // # of custom attributes, some may be deleted
        public int[] etype;               // type = INT/DOUBLE of each attribute
        public int[] esize;               // size = 0 for vector, N for array columns
        public int[] ewhich;              // index into eivec,eiarray,edvec,edarray for data
         
        public int[][] eivec;              // pointer to each integer vector
        public int[][][] eiarray;           // pointer to each integer array
        public double[][] edvec;           // pointer to each double vector
        public double[][][] edarray;        // pointer to each double array
         
         // restart buffers, filled by read_restart
         
        public int nlocal_restart;
        public string particle_restart;
        
        public int copy, copymode;        // 1 if copy of class (prevents deallocation of
                                          //  base class when child copy is destroyed)

        // methods

        public void init()
        {
            for (int i = 0; i < nmixture; i++) mixture[i].init();

            // RNG for particle weighting

            if (wrandom==null)
            {
                wrandom = new RanPark(sparta.update.ranmaster.uniform());
                double seed = sparta.update.ranmaster.uniform();
                wrandom.reset(seed, me, 100);
            }

            // if first run after reading a restart file,
            // delete any custom particle attributes that have not been re-defined
            // use nactive since remove_custom() may alter ncustom

            if (custom_restart_flag!= null)
            {
                int nactive = ncustom;
                for (int i = 0; i < nactive; i++)
                    if (custom_restart_flag[i] == 0) remove_custom(i);
                //delete[] custom_restart_flag;
                custom_restart_flag = null;
            }

            // reallocate cellcount and first lists as needed
            // NOTE: when grid becomes dynamic, will need to do this in sort()

            //if (maxgrid < sparta.grid.nlocal) {
            //  maxgrid = sparta.grid.nlocal;
            //    memory.destroy(cellcount);
            //memory.destroy(first);
            //memory.create(first,maxgrid,"particle:first");
            //memory.create(cellcount,maxgrid,"particle:cellcount");
            // }
        }
        public virtual void compress_migrate(int nmigrate, int[] mlist)
        {
            int i, j, k;
            int nbytes = Marshal.SizeOf(typeof(OnePart));

            if (ncustom==0)
            {
                for (i = 0; i < nmigrate; i++)
                {
                    j = mlist[i];
                    k = nlocal - 1;
                    while (k == mlist[nmigrate - 1] && k > j)
                    {
                        nmigrate--;
                        nlocal--;
                        k--;
                    }
                    nlocal--;
                    if (j == k) continue;
                    //memcpy(&particles[j], &particles[k], nbytes);
                    Array.Copy(particles, k, particles, j, 1);
                }

            }
            else
            {
                for (i = 0; i < nmigrate; i++)
                {
                    j = mlist[i];
                    k = nlocal - 1;
                    while (k == mlist[nmigrate - 1] && k > j)
                    {
                        nmigrate--;
                        nlocal--;
                        k--;
                    }
                    nlocal--;
                    if (j == k) continue;
                    //memcpy(&particles[j], &particles[k], nbytes);
                    Array.Copy(particles, k, particles, j, 1);
                    copy_custom(j, k);
                }
            }

            sorted = 0;
        }
        public void compress_rebalance()
        {
            int nbytes = 96;

            if (ncustom==0)
            {
                int i = 0;
                while (i < nlocal)
                {
                    if (particles[i].icell < 0)
                    {
                        particles[i] = particles[nlocal - 1];
                        //memcpy(&particles[i], &particles[nlocal - 1], nbytes);
                        nlocal--;
                    }
                    else i++;
                }

            }
            else
            {
                int i = 0;
                while (i < nlocal)
                {
                    if (particles[i].icell < 0)
                    {
                        particles[i] = particles[nlocal - 1];
                        //memcpy(&particles[i], &particles[nlocal - 1], nbytes);
                        copy_custom(i, nlocal - 1);
                        nlocal--;
                    }
                    else i++;
                }
            }

            sorted = 0;
        }
        public void compress_reactions(int ndelete, int[] dellist)
        {
            int i;

            int nbytes = Marshal.SizeOf(typeof(OnePart));

            // reallocate next list as needed

            if (maxsort < maxlocal)
            {
                maxsort = maxlocal;
                //memory.destroy(next);
                //memory.create(next, maxsort, "particle:next");
                next = new int[maxsort];
            }

            // use next as a scratch vector
            // is always defined when performing collisions and thus gas reactions
            // next is only used for upper locs from nlocal-ndelete to nlocal
            // next[i] = current index of atom originally at index i, when i >= nlocal
            // next[i] = original index of atom currently at index i, when i <= nlocal

            int upper = nlocal - ndelete;
            for (i = upper; i < nlocal; i++) next[i] = i;

            // i = current index of atom to remove, even if it previously moved

            if (ncustom==0)
            {
                for (int m = 0; m < ndelete; m++)
                {
                    i = dellist[m];
                    if (i >= nlocal) i = next[i];
                    nlocal--;
                    if (i == nlocal) continue;
                    //memcpy(&particles[i], &particles[nlocal], nbytes);
                    particles[i] = particles[nlocal];
                    if (i >= upper) next[i] = next[nlocal];
                    next[next[nlocal]] = i;
                }

            }
            else
            {
                for (int m = 0; m < ndelete; m++)
                {
                    i = dellist[m];
                    if (i >= nlocal) i = next[i];
                    nlocal--;
                    if (i == nlocal) continue;
                    //memcpy(&particles[i], &particles[nlocal], nbytes);
                    particles[i] = particles[nlocal];
                    copy_custom(i, nlocal);
                    if (i >= upper) next[i] = next[nlocal];
                    next[next[nlocal]] = i;
                }
            }
        }
        public void sort()
        {
            sorted = 1;

            // reallocate next list as needed
            // NOTE: why not just compare maxsort to nlocal?
            //       then could realloc less often?
            //       ditto for all reallocs of next in related methods

            if (maxsort < maxlocal)
            {
                maxsort = maxlocal;
                next = new cellint[maxsort];

            }
            // initialize linked list of particles in cells I own

            Grid.ChildInfo[] cinfo = sparta.grid.cinfo;
            int nglocal = sparta.grid.nlocal;

            for (int icells = 0; icells < nglocal; icells++)
            {
                cinfo[icells].first = -1;
                cinfo[icells].count = 0;
                //cellcount[i] = 0;
                //first[i] = -1;
            }
            // reverse loop over partlcles to store linked lists in forward order
            // icell = global cell the particle is in

            int icell;
            for (int i = nlocal - 1; i >= 0; i--)
            {
                icell = particles[i].icell;
                next[i] = cinfo[icell].first;
                cinfo[icell].first = i;
                cinfo[icell].count++;

                // NOTE: this method seems much slower for some reason
                // uses separate, smaller vectors for first & cellcount
                //icell = cells[particles[i].icell].local;
                //next[i] = first[icell];
                //first[icell] = i;
                //cellcount[icell]++;
            }
        }
        public void sort_allocate()
        {
            if (maxsort < maxlocal)
            {
                maxsort = maxlocal;
                //memory.destroy(next);
                //memory.create(next, maxsort, "particle:next");
                next = new int[maxsort];
            }
        }
        public void remove_all_from_cell(int ip)
        {
            while (ip >= 0)
            {
                particles[ip].icell = -1;
                ip = next[ip];
            }
        }
        public virtual void grow(int nextra)
        {
            bigint target = (bigint)nlocal + nextra;
            if (target <= maxlocal) return;

            int oldmax = maxlocal;
            bigint newmax = maxlocal;
            while (newmax < target) newmax += DELTA;

            if (newmax > Run.MAXSMALLINT)
                sparta.error.one("Per-processor particle count is too big");

            maxlocal = (int)newmax;
            //particles = (OnePart*)
            //  memory.srealloc(particles, maxlocal * sizeof(OnePart),
            //           "particle:particles");
            //memset(&particles[oldmax], 0, (maxlocal - oldmax) * sizeof(OnePart));
            if (particles==null)
            {
                particles = new OnePart[maxlocal];
            }
            else
            {
                Array.Resize<OnePart>(ref particles, maxlocal);
            }
            if (ncustom == 0) return;

            for (int i = 0; i < ncustom; i++)
            {
                if (ename[i] == null) continue;
                grow_custom(i, oldmax, maxlocal);
            }
        }
        public virtual void grow_species()
        {
            species = new List<Species>(maxspecies);
            //Console.WriteLine("particle.add_species().grow_species");
        }
        //public void grow_next();
        public virtual void pre_weight()
        {
            int icell;
            Grid.ChildInfo[] cinfo = sparta.grid.cinfo;

            for (int i = 0; i < nlocal; i++)
            {
                icell = particles[i].icell;
                particles[i].weight = cinfo[icell].weight;
            }
        }
        public virtual void post_weight()
        {
            int m, icell, nclone;
            double ratio, fraction;

            int nbytes = Marshal.SizeOf(typeof(OnePart));
            Grid.ChildInfo[] cinfo = sparta.grid.cinfo;

            // nlocal_original-1 = index of last original particle

            int nlocal_original = nlocal;
            int i = 0;

            while (i < nlocal_original)
            {
                icell = particles[i].icell;

                // next particle will be an original particle
                // skip it if no weight change

                if (particles[i].weight == cinfo[icell].weight)
                {
                    i++;
                    continue;
                }

                // ratio < 1.0 is candidate for deletion
                // if deleted and particle that takes its place is cloned (Nloc > Norig)
                //   then skip it via i++, else will examine it on next iteration

                ratio = particles[i].weight / cinfo[icell].weight;

                if (ratio < 1.0)
                {
                    if (wrandom.uniform() > ratio)
                    {
                        //memcpy(&particles[i], &particles[nlocal - 1], nbytes);
                        particles[i] =particles[nlocal - 1];
                        if (ncustom!=0) copy_custom(i, nlocal - 1);
                        if (nlocal > nlocal_original) i++;
                        else nlocal_original--;
                        nlocal--;
                    }
                    else i++;
                    continue;
                }

                // ratio > 1.0 is candidate for cloning
                // create Nclone new particles each with unique ID 

                nclone = Convert.ToInt32(ratio);
                fraction = ratio - nclone;
                nclone--;
                if (wrandom.uniform() < fraction) nclone++;

                for (m = 0; m < nclone; m++)
                {
                    clone_particle(i);
                    particles[nlocal - 1].id = (int)(Run.MAXSMALLINT * wrandom.uniform());
                }
                i++;
            }
        }

        public virtual int add_particle(int id, int ispecies, int icell,
                           double[] x, double[] v, double erot, double evib)
        {
            int reallocflag = 0;
            if (nlocal == maxlocal)
            {
                grow(1);
                reallocflag = 1;
            }

            OnePart p = particles[nlocal];

            p.id = id;
            p.ispecies = ispecies;
            p.icell = icell;
            p.x = new double[3];
            p.x[0] = x[0];
            p.x[1] = x[1];
            p.x[2] = x[2];
            p.v = new double[3];
            p.v[0] = v[0];
            p.v[1] = v[1];
            p.v[2] = v[2];
            p.erot = erot;
            p.evib = evib;
            p.flag = (int)Enum1.PKEEP;

            //p.dtremain = 0.0;    not needed due to memset in grow() ??
            //p.weight = 1.0;      not needed due to memset in grow() ??
            particles[nlocal] = p;
            nlocal++;
            return reallocflag;
        }
        public int clone_particle(int index)
        {
            int reallocflag = 0;
            if (nlocal == maxlocal)
            {
                grow(1);
                reallocflag = 1;
            }

            //memcpy(&particles[nlocal], &particles[index], sizeof(OnePart));
            particles[nlocal] = particles[index];
            if (ncustom!=0) copy_custom(nlocal, index);

            nlocal++;
            return reallocflag;
        }
        public void add_species(int narg, string[] args)
        {
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);
            int i, j, n; ;

            if (narg < 2) sparta.error.all("Illegal species command");

            if (me == 0)
            {
                try
                {
                    fp = new FileStream(arg[0], FileMode.Open, FileAccess.Read);
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message.ToString());
                } 
                if (fp == null)
                {
                    string str = string.Format("Cannot open species file {0}", arg[0]);
                    sparta.error.one(str);
                }
            }
            // nfilespecies = # of species defined in file
            // filespecies = list of species defined in file

            nfilespecies = maxfilespecies = 0;
            filespecies = null;

            if (me == 0) read_species_file();
            sparta.mpi.MPI_Bcast(ref nfilespecies, 1, MPI.MPI_INT, 0, sparta.world);
            if (sparta.comm.me != 0)
            {
                filespecies = new List<Species>(nfilespecies);
            }
            sparta.mpi.MPI_Bcast(ref filespecies, nfilespecies * Marshal.SizeOf(typeof(Species)), MPI.MPI_BYTE, 0, sparta.world);

            // newspecies = # of new user-requested species
            // names = list of new species IDs
            // customize abbreviations by adding new keyword in 2 places

            //char line[MAXLINE];
            string line;

            int newspecies = 0;
            for (int iarg = 1; iarg < narg; iarg++)
            {
                if (string.Equals(arg[iarg], "air"))
                {
                    line = string.Copy(AIR);
                    newspecies += wordcount(line);
                }
                else newspecies++;
            }
            string[] names = new string[newspecies];
            newspecies = 0;

            for (int iarg = 1; iarg < narg; iarg++)
            {
                if (string.Equals(arg[iarg], "air"))
                {
                    line = string.Copy(AIR);
                    newspecies += wordcount(line, ref names);
                }
                else names[newspecies++] = arg[iarg];
            }
            // species ID must be all alphanumeric chars, underscore, plus/minus

            for (i = 0; i < newspecies; i++)
            {
                n = names[i].Length;
                for (j = 0; j < n - 1; j++)
                    if (!char.IsLetterOrDigit(names[i][j]) && names[i][j] != '_' &&
                        names[i][j] != '+' && names[i][j] != '-')
                        sparta.error.all("Invalid character in species ID");
            }
            // extend species list if necessary

            if (nspecies + newspecies > maxspecies)
            {
                while (nspecies + newspecies > maxspecies) maxspecies += DELTASPECIES;
                grow_species();

                //Console.ReadKey();
            }
            // extract info on user-requested species from file species list
            // add new species to default mixtures "all" and "species"

            int imix_all = find_mixture("all");
            int imix_species = find_mixture("species");

            for (i = 0; i < newspecies; i++)
            {
                for (j = 0; j < nspecies; j++)
                    if (string.Equals(names[i], species[j].id)) break;
                if (j < nspecies) sparta.error.all("Species ID is already defined");
                for (j = 0; j < nfilespecies; j++)
                    if (string.Equals(names[i], filespecies[j].id)) break;
                if (j == nfilespecies)
                    sparta.error.all("Species ID does not appear in species file");
                //memcpy(&species[nspecies], &filespecies[j], sizeof(Species));
                species.Add( filespecies[j]);
                nspecies++;


                mixture[imix_all].add_species_default(species[nspecies - 1].id);
                mixture[imix_species].add_species_default(species[nspecies - 1].id);
            }

        }
        public void add_mixture(int narg, string[] arg)
        {
            if (narg < 1) sparta.error.all("Illegal mixture command");

            // imix = index if mixture ID already exists
            // else instantiate a new mixture

            int imix = find_mixture(arg[0]);

            if (imix < 0)
            {
                if (nmixture == maxmixture)
                {
                    maxmixture += DELTAMIXTURE;
                    mixture = new List<Mixture>(maxmixture);

                }
                imix = nmixture;
                nmixture++;
                mixture.Add(new Mixture(sparta, arg[0]));
            }

            mixture[imix].command(narg, arg);
        }
        public int find_species(string id)
        {
            for (int i = 0; i < nspecies; i++)
                if (string.Equals(id, species[i].id)) return i;
            return -1;
        }
        public int find_mixture(string id)
        {
            for (int i = 0; i < nmixture; i++)
                if (string.Equals(id, mixture[i].id)) return i;
            return -1;
        }
        public double erot(int isp, double temp_thermal, RanPark erandom)
        {
            double eng, a, erm, b;

            if (sparta.collide==null || sparta.collide.rotstyle == (int)Enum2.NONE) return 0.0;
            if (species[isp].rotdof < 2) return 0.0;

            if (species[isp].rotdof == 2)
                eng = -Math.Log(erandom.uniform()) * sparta.update.boltz * temp_thermal;
            else
            {
                a = 0.5 * sparta.particle.species[isp].rotdof - 1.0;
                while (true)
                {
                    // energy cut-off at 10 kT
                    erm = 10.0 * erandom.uniform();
                    b = Math.Pow(erm / a, a) * Math.Exp(a - erm);
                    if (b > erandom.uniform()) break;
                }
                eng = erm * sparta.update.boltz * temp_thermal;
            }

            return eng;
        }
        public double evib(int isp, double temp_thermal, RanPark erandom)
        {
            double eng, a, erm, b;

            int vibstyle = (int)Enum2.NONE;
            if (sparta.collide!=null) vibstyle = sparta.collide.vibstyle;
            if (vibstyle == (int)Enum2.NONE || species[isp].vibdof < 2) return 0.0;

            eng = 0.0;
            if (vibstyle == (int)Enum2.DISCRETE && species[isp].vibdof == 2)
            {
                int ivib = Convert.ToInt32(-Math.Log(erandom.uniform()) * temp_thermal /
                                             sparta.particle.species[isp].vibtemp);
                eng = ivib * sparta.update.boltz * sparta.particle.species[isp].vibtemp;
            }
            else if (vibstyle == (int)Enum2.SMOOTH || species[isp].vibdof >= 2)
            {
                if (species[isp].vibdof == 2)
                    eng = -Math.Log(erandom.uniform()) * sparta.update.boltz * temp_thermal;
                else if (species[isp].vibdof > 2)
                {
                    a = 0.5 * sparta.particle.species[isp].vibdof - 1;
                    while (true)
                    {
                        // energy cut-off at 10 kT
                        erm = 10.0 * erandom.uniform();
                        b = Math.Pow(erm / a, a) * Math.Exp(a - erm);
                        if (b > erandom.uniform()) break;
                    }
                    eng = erm * sparta.update.boltz * temp_thermal;
                }
            }

            return eng;
        }

        //public  void write_restart_species(FILE* fp);
        //public void read_restart_species(FILE* fp);
        //public void write_restart_mixture(FILE* fp);
        //public void read_restart_mixture(FILE* fp);

        //public int size_restart();
        //public int pack_restart(char*);
        //public int unpack_restart(char*);

        public int find_custom(string name)
        {
            for (int i = 0; i < ncustom; i++)
                if (ename[i]!=null && string.Equals(ename[i], name)) return i;
            return -1;
        }
        //public int add_custom(char*, int, int);
        public void grow_custom(int index, int nold, int nnew)
        {
            if (etype[index] == (int)Enum3.INT)
            {
                if (esize[index] == 0)
                {
                    int[] ivector = eivec[ewhich[index]];
                    //memory.grow(ivector, nnew, "particle:eivec");
                    int[] tempvector = new int[nnew];
                    if (ivector!=null)
                    {
                        Array.Copy(ivector, tempvector, ivector.Length);
                        //memset(&ivector[nold], 0, (nnew - nold) * sizeof(int));
                        for (int i = nold+1; i < nnew; i++)
                        {
                            tempvector[i] = 0;
                        }
                        ivector = new int[nnew];
                        Array.Copy(tempvector, ivector, nnew);
                    }
                    else
                    {
                        ivector = new int[nnew];
                    }

                    eivec[ewhich[index]] = ivector;
                }
                else
                {
                    int[][] iarray = eiarray[ewhich[index]];
                    //memory.grow(iarray, nnew, esize[index], "particle:eiarray");
                    
                    if (iarray != null)
                    {
                        int[][] temparray = new int[nnew][];
                        for (int i = 0; i < nold; i++)
                        {
                            //temparray[i] = iarray[i];
                            Array.Copy(iarray[i], temparray[i], iarray[i].Length);
                        }
                        for (int i = nold; i < nnew; i++)
                        {
                            temparray[i] = new int[esize[index]];
                        }
                        iarray = temparray;
                    }
                    else
                    {
                        iarray = new int[nnew][];
                        for (int i = 0; i < nnew; i++)
                        {
                            iarray[i] = new int[esize[index]];
                        }
                    }
                        
                    
                        //memset(&iarray[nold][0], 0, (nnew - nold) * esize[index] * sizeof(int));
                    eiarray[ewhich[index]] = iarray;
                }

            }
            else
            {
                if (esize[index] == 0)
                {
                    double[] dvector = edvec[ewhich[index]];
                    double[] tempvector = new double[nnew];
                    Array.Copy(dvector, tempvector, dvector.Length);
                    for (int i = nold+1; i < nnew; i++)
                    {
                        tempvector[i] = 0;
                    }
                    //memory.grow(dvector, nnew, "particle:edvec");
                    dvector = new double[nnew];
                    Array.Copy(tempvector, dvector, nnew);
                    //if (dvector!=null) memset(&dvector[nold], 0, (nnew - nold) * sizeof(double));
                    edvec[ewhich[index]] = dvector;
                }
                else
                {
                    double[][] darray = edarray[ewhich[index]];
                    //memory.grow(darray, nnew, esize[index], "particle:edarray");
                    if (darray != null)
                    {
                        double[][] tempdarray = new double[nnew][];
                        for (int i = 0; i < nold; i++)
                        {
                            //tempdarray[i] = darray[i];
                            Array.Copy(darray[i], tempdarray[i], darray[i].Length);
                        }
                        for (int i = nold+1; i < nnew; i++)
                        {
                            tempdarray[i] = new double[esize[index]];
                        }
                    }
                        //memset(&darray[nold][0], 0, (nnew - nold) * esize[index] * sizeof(double));
                    edarray[ewhich[index]] = darray;
                }
            }
        }
        public void remove_custom(int index)
        {
            //delete[] ename[index];
            ename[index] = null;

            if (etype[index] == (int)Enum3.INT)
            {
                if (esize[index] == 0)
                {
                    //memory.destroy(eivec[ewhich[index]]);
                    eivec[ewhich[index]] = null;
                    ncustom_ivec--;
                    for (int i = ewhich[index]; i < ncustom_ivec; i++)
                    {
                        icustom_ivec[i] = icustom_ivec[i + 1];
                        ewhich[icustom_ivec[i]] = i;
                        eivec[i] = eivec[i + 1];
                    }
                }
                else
                {
                    //memory.destroy(eiarray[ewhich[index]]);
                    eiarray[ewhich[index]] = null;
                    ncustom_iarray--;
                    for (int i = ewhich[index]; i < ncustom_iarray; i++)
                    {
                        icustom_iarray[i] = icustom_iarray[i + 1];
                        ewhich[icustom_iarray[i]] = i;
                        eiarray[i] = eiarray[i + 1];
                        eicol[i] = eicol[i + 1];
                    }
                }
            }
            else if (etype[index] == (int)Enum3.DOUBLE)
            {
                if (esize[index] == 0)
                {
                    //memory.destroy(edvec[ewhich[index]]);
                    edvec[ewhich[index]] = null;
                    ncustom_dvec--;
                    for (int i = ewhich[index]; i < ncustom_dvec; i++)
                    {
                        icustom_dvec[i] = icustom_dvec[i + 1];
                        ewhich[icustom_dvec[i]] = i;
                        edvec[i] = edvec[i + 1];
                    }
                }
                else
                {
                    //memory.destroy(edarray[ewhich[index]]);
                    edarray[ewhich[index]] = null;
                    ncustom_darray--;
                    for (int i = ewhich[index]; i < ncustom_darray; i++)
                    {
                        icustom_darray[i] = icustom_darray[i + 1];
                        ewhich[icustom_darray[i]] = i;
                        edarray[i] = edarray[i + 1];
                        edcol[i] = edcol[i + 1];
                    }
                }
            }

            // set ncustom = 0 if custom list is now entirely empty

            int empty = 1;
            for (int i = 0; i < ncustom; i++)
                if (ename[i]!=null) empty = 0;
            if (empty != 0) ncustom = 0;
        }
        public void copy_custom(int i, int j)
        {
            int m;

            // caller does not always check this
            // shouldn't be a problem, but valgrind can complain if memcpy to self
            // oddly memcpy(&particles[i],&particles[j],sizeof(OnePart)) seems OK

            if (i == j) return;

            // 4 flavors of vectors/arrays

            if (ncustom_ivec!=0)
            {
                for (m = 0; m < ncustom_ivec; m++) eivec[m][i] = eivec[m][j];
            }
            if (ncustom_iarray != 0)
            {
                for (m = 0; m < ncustom_iarray; m++)
                {
                    Array.Copy(eiarray[m][j], eiarray[m][i], eicol[m]);
                    //memcpy(eiarray[m][i], eiarray[m][j], eicol[m] * sizeof(int));
                }
            }
            if (ncustom_dvec != 0)
            {
                for (m = 0; m < ncustom_dvec; m++) edvec[m][i] = edvec[m][j];
            }
            if (ncustom_darray != 0)
            {
                for (m = 0; m < ncustom_darray; m++)
                {
                    Array.Copy(edarray[m][j], edarray[m][i], edcol[m]);
                    //memcpy(edarray[m][i], edarray[m][j], edcol[m] * sizeof(double));
                }
            }
        }
        public int sizeof_custom()
        {
            int n = 0;

            n += ncustom_ivec * sizeof(int);
            if (ncustom_iarray!=0)
                for (int i = 0; i < ncustom_iarray; i++)
                    n += eicol[i] * sizeof(int);

            n = IROUNDUP(n);

            n += ncustom_dvec * sizeof(double);
            if (ncustom_darray != 0)
                for (int i = 0; i < ncustom_darray; i++)
                    n += edcol[i] * sizeof(double);

            return n;
        }
        //public void write_restart_custom(FILE* fp);
        //public void read_restart_custom(FILE* fp);
        public void pack_custom(int n,ref StringBuilder buf)
        {
            int i, ptr = 0;

            if (ncustom_ivec!=0)
            {
                for (i = 0; i < ncustom_ivec; i++)
                {
                    //memcpy(ptr, &eivec[i][n], sizeof(int));
                    buf.Append(eivec[i][n]);
                    ptr += sizeof(int);
                }
            }
            if (ncustom_iarray!=0)
            {
                for (i = 0; i < ncustom_iarray; i++)
                {
                    //memcpy(ptr, eiarray[i][n], eicol[i] * sizeof(int));
                    buf.Append(eiarray[i][n]);
                    ptr += eicol[i] * sizeof(int);
                }
            }


            if (ncustom_dvec!=0)
            {
                for (i = 0; i < ncustom_dvec; i++)
                {
                    //memcpy(ptr, &edvec[i][n], sizeof(double));
                    buf.Append(edvec[i][n]);
                    ptr += sizeof(double);
                }
            }
            if (ncustom_darray!=0)
            {
                for (i = 0; i < ncustom_darray; i++)
                {
                    //memcpy(ptr, edarray[i][n], edcol[i] * sizeof(double));
                    buf.Append(edarray[i][n]);
                    ptr += edcol[i] * sizeof(double);
                }
            }

        }
        public void unpack_custom(byte[] buf, int n)
        {
            Console.WriteLine("Particle.unpack_custom()");
            //int i;
            //char* ptr = buf;

            //if (ncustom_ivec!=0)
            //{
            //    for (i = 0; i < ncustom_ivec; i++)
            //    {
            //        memcpy(&eivec[i][n], ptr, sizeof(int));
            //        ptr += sizeof(int);
            //    }
            //}
            //if (ncustom_iarray)
            //{
            //    for (i = 0; i < ncustom_iarray; i++)
            //    {
            //        memcpy(eiarray[i][n], ptr, eicol[i] * sizeof(int));
            //        ptr += eicol[i] * sizeof(int);
            //    }
            //}

            //ptr = ROUNDUP(ptr);

            //if (ncustom_dvec)
            //{
            //    for (i = 0; i < ncustom_dvec; i++)
            //    {
            //        memcpy(&edvec[i][n], ptr, sizeof(double));
            //        ptr += sizeof(double);
            //    }
            //}
            //if (ncustom_darray)
            //{
            //    for (i = 0; i < ncustom_darray; i++)
            //    {
            //        memcpy(edarray[i][n], ptr, edcol[i] * sizeof(double));
            //        ptr += edcol[i] * sizeof(double);
            //    }
            //}
        }

        public bigint memory_usage()
        {
            bigint bytes = (bigint)maxlocal * Marshal.SizeOf(typeof(OnePart));
            bytes += (bigint)maxlocal * sizeof(int);
            for (int i = 0; i < ncustom_ivec; i++)
                bytes += (bigint)maxlocal * sizeof(int);
            for (int i = 0; i < ncustom_iarray; i++)
                bytes += (bigint)maxlocal * eicol[i] * sizeof(int);
            for (int i = 0; i < ncustom_dvec; i++)
                bytes += (bigint)maxlocal * sizeof(double);
            for (int i = 0; i < ncustom_darray; i++)
                bytes += (bigint)maxlocal * edcol[i] * sizeof(double);
            return bytes;
        }
        private SPARTA sparta;
        public Particle(SPARTA sparta)
        {
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);

            exist = sorted = 0;
            nglobal = 0;
            nlocal = maxlocal = 0;
            particles = null;

            nspecies = maxspecies = 0;
            species = null;

            //maxgrid = 0;
            //cellcount = null;
            //first = null;
            maxsort = 0;
            next = null;

            // create two default mixtures

            nmixture = maxmixture = 0;
            mixture = null;

            string[] newarg = new string[1];
            newarg[0] = "all";
            add_mixture(1, newarg);
            newarg[0] = "species";
            add_mixture(1, newarg);
            //delete[] newarg;

            // custom per-particle vectors/arrays

            ncustom = 0;
            ename = null;
            etype = esize = ewhich = null;

            ncustom_ivec = ncustom_iarray = 0;
            icustom_ivec = icustom_iarray = null;
            eivec = null;
            eiarray = null;
            eicol = null;

            ncustom_dvec = ncustom_darray = 0;
            icustom_dvec = icustom_darray = null;
            edvec = null;
            edarray = null;
            edcol = null;

            custom_restart_flag = null;

            // RNG for particle weighting

            wrandom = null;

            copy = copymode = 0;
        }

        protected int me;
        protected int maxgrid;              // max # of indices first can hold
        protected int maxsort;              // max # of particles next can hold
        protected int maxspecies;           // max size of species list

        protected List<Species> filespecies;     // list of species read from file
        protected int nfilespecies;         // # of species read from file
        protected int maxfilespecies;       // max size of filespecies list
        protected FileStream fp;                 // species file pointer

        protected RanPark wrandom;   // RNG for particle weighting

        // extra custom vectors/arrays for per-particle data
        // ncustom > 0 if there are any extra arrays
        // these varaiables are private, others above are public

        protected string[] ename;             // name of each attribute

        protected int ncustom_ivec;         // # of integer vector attributes
        protected int ncustom_iarray;       // # of integer array attributes
        protected int[] icustom_ivec;        // index into ncustom for each integer vector
        protected int[] icustom_iarray;      // index into ncustom for each integer array
        protected int[] eicol;               // # of columns in each integer array (esize)

        protected int ncustom_dvec;         // # of double vector attributes
        protected int ncustom_darray;       // # of double array attributes
        protected int[] icustom_dvec;        // index into ncustom for each double vector
        protected int[] icustom_darray;      // index into ncustom for each double array
        protected int[] edcol;               // # of columns in each double array (esize)

        protected int[] custom_restart_flag; // flag on each custom vec/array read from restart
                                             // used to delete them if not redefined in 
                                             // restart script

        // private methods

        private void read_species_file()
        {
            using (StreamReader sr = new StreamReader(fp))
            {
                int NWORDS = 10;
                string line;
                while ((line=sr.ReadLine())!=null)
                {
                    ParseSpecies(line);
                    
                }

            }

        }

        private void ParseSpecies(string line)
        {
            if (string.IsNullOrWhiteSpace(line))
            {
                return;
            }
            if (line.StartsWith("#"))
            {
                return;
            }

            string[] words = line.Split();
            List<string> wordlist = new List<string>();
            foreach (string word in words)
            {
                if (!string.IsNullOrWhiteSpace(word))
                {
                    wordlist.Add(word);
                }
            }

            if (wordlist.Count!=10)
            {
                sparta.error.one("Incorrect line format in species file");
            }

            if (nfilespecies==maxfilespecies)
            {
                maxfilespecies += DELTASPECIES;
                filespecies = new List<Species>(maxfilespecies);
            }
            Species fsp = new Species();

            if (wordlist[0].Length+1>16)
            {
                sparta.error.one("Invalid species ID in species file");
            }
            fsp.id = string.Copy(wordlist[0]);
            fsp.molwt = double.Parse(wordlist[1]);
            fsp.mass= double.Parse(wordlist[2]);
            fsp.rotdof = int.Parse(wordlist[3]);
            fsp.rotrel = double.Parse(wordlist[4]);
            fsp.vibdof = int.Parse(wordlist[5]);
            fsp.vibrel = double.Parse(wordlist[6]);
            fsp.vibtemp = double.Parse(wordlist[7]);
            fsp.specwt = double.Parse(wordlist[8]);
            fsp.charge= double.Parse(wordlist[9]);
            if (fsp.rotdof > 0 || fsp.vibdof > 0)
            {
                fsp.internaldof = 1;
            }
            else
            {
                fsp.internaldof = 0;
            }
            filespecies.Add(fsp);
            nfilespecies++;
        }

        private int wordcount(string line, ref string[] words)
        {
            //int nwords = 0;

            //string word = strtok(line, " \t");

            //while (word)
            //{
            //    if (words) words[nwords] = word;
            //    nwords++;
            //    word = strtok(NULL, " \t");
            //}
            words = line.Split('\t');

            return words.Length;
        }
        private int wordcount(string line)
        {
            string[] words = line.Split('\t');

            return words.Length;
        }
        public static int IROUNDUP(int A)
        {
            return ((((A) + 7) / 8) * 8);
        }
    }
}