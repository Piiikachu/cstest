using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;
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

        public struct OnePart
        {
            public int id;                 // particle ID
            public int ispecies;           // particle species index
            public int icell;              // which local Grid::cells the particle is in
            public double[] x;            // particle position
            public double[] v;            // particle velocity
            public double erot;            // rotational energy
            public double evib;            // vibrational energy
            public int flag;               // used for migration status
            public double dtremain;        // portion of move timestep remaining
            public double weight;          // particle or cell weight, if weighting enabled
        }

        public struct OnePartRestart
        {
            public int id;                 // particle ID
            public int ispecies;           // particle species index
            public cellint icell;          // cell ID the particle is in
            public int nsplit;             // 1 for unsplit cell
                                           // else neg of sub cell index (0 to Nsplit-1)
            public double[] x;            // particle position
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
         
        public int[,] eivec;              // pointer to each integer vector
        public int[,,] eiarray;           // pointer to each integer array
        public double[,] edvec;           // pointer to each double vector
        public double[,,] edarray;        // pointer to each double array
         
         // restart buffers, filled by read_restart
         
        public int nlocal_restart;
        public string particle_restart;
        
        public int copy, copymode;        // 1 if copy of class (prevents deallocation of
                                   //  base class when child copy is destroyed)

        // methods

        //public void init();
        //public virtual void compress_migrate(int, int*);
        //public void compress_rebalance();
        //public void compress_reactions(int, int*);
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
        //public void sort_allocate();
        //public void remove_all_from_cell(int);
        //public virtual void grow(int);
        public virtual void grow_species()
        {
            species = new List<Species>(maxspecies);
            //Console.WriteLine("particle.add_species()->grow_species");
        }
        //public void grow_next();
        //public virtual void pre_weight();
        //public virtual void post_weight();

        //public virtual int add_particle(int, int, int, double*, double*, double, double);
        //public int clone_particle(int);
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
        //public double erot(int, double, class RanPark *);
        //public double evib(int, double, class RanPark *);

        //public  void write_restart_species(FILE* fp);
        //public void read_restart_species(FILE* fp);
        //public void write_restart_mixture(FILE* fp);
        //public void read_restart_mixture(FILE* fp);

        //public int size_restart();
        //public int pack_restart(char*);
        //public int unpack_restart(char*);

        //public int find_custom(char*);
        //public int add_custom(char*, int, int);
        //public void grow_custom(int, int, int);
        //public void remove_custom(int);
        //public void copy_custom(int, int);
        //public int sizeof_custom();
        //public void write_restart_custom(FILE* fp);
        //public void read_restart_custom(FILE* fp);
        //public void pack_custom(int, char*);
        //public void unpack_custom(char*, int);

        //public bigint memory_usage();
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
    }
}