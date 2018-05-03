using System;
using System.Collections.Generic;
using System.Text;
using bigint = System.Int64;
namespace cstest
{
    public class Collide
    {
        public const int DELTAPART = 128;

        enum Enum1{ NONE, DISCRETE, SMOOTH };       // several files
        enum Enum2{ PKEEP, PINSERT, PDONE, PDISCARD, PENTRY, PEXIT, PSURF };   // several files

        public const int DELTAGRID = 1000;        // must be bigger than split cells per cell
        public const int DELTADELETE = 1024;
        public const int DELTAELECTRON = 128;

        public const double BIG = 1.0e20;

        public string style;
        public int rotstyle;       // none/smooth rotational modes
        public int vibstyle;       // none/discrete/smooth vibrational modes
        public int nearcp;         // 1 for near neighbor collisions
        public int nearlimit;      // limit on neighbor serach for near neigh collisions
         
        public int ncollide_one, nattempt_one, nreact_one;
        public bigint ncollide_running, nattempt_running, nreact_running;
        private SPARTA sparta;

        public Collide(SPARTA sparta, int narg, string[] arg)
        {
            this.sparta = sparta;
            int n = arg[0].Length + 1;
            style = string.Copy(arg[0]);

            n = arg[1].Length + 1;
            mixID = string.Copy(arg[1]);

            random = new RanPark(sparta.update.ranmaster.uniform());
            double seed = sparta.update.ranmaster.uniform();
            random.reset(seed, sparta.comm.me, 100);

            ngroups = 0;

            npmax = 0;
            plist = null;

            nglocal = nglocalmax = 0;

            ngroup = null;
            maxgroup = null;
            glist = null;
            gpair = null;

            maxdelete = 0;
            dellist = null;

            vre_first = 1;
            vre_start = 1;
            vre_every = 0;
            remainflag = 1;
            vremax = null;
            vremax_initial = null;
            remain = null;
            rotstyle =(int)Enum1.SMOOTH;
            vibstyle = (int)Enum1.NONE;
            nearcp = 0;
            nearlimit = 10;

            recomb_ijflag = null;

            ambiflag = 0;
            maxelectron = 0;
            elist = null;
        }
        public virtual void init()
        {
            // error check

            if (ambiflag!=0 && nearcp != 0)
                sparta.error.all("Ambipolar collision model does not yet support near-neighbor collisions");

            // require mixture to contain all species

            int imix = sparta.particle.find_mixture(mixID);
            if (imix < 0) sparta.error.all("Collision mixture does not exist");
            mixture = sparta.particle.mixture[imix];

            if (mixture.nspecies != sparta.particle.nspecies)
                sparta.error.all("Collision mixture does not contain all species");

            if (sparta.kokkos != null && kokkosable == 0)
                sparta.error.all("Must use Kokkos-supported collision style if Kokkos is enabled");

            // reallocate one-cell data structs for one or many groups

            oldgroups = ngroups;
            ngroups = mixture.ngroup;

            if (ngroups != oldgroups)
            {
                if (oldgroups == 1)
                {
                    //memory.destroy(plist);
                    npmax = 0;
                    plist = null;
                }
                if (oldgroups > 1)
                {
                    //delete[] ngroup;
                    //delete[] maxgroup;
                    //for (int i = 0; i < oldgroups; i++) memory.destroy(glist[i]);
                    //delete[] glist;
                    //memory.destroy(gpair);
                    ngroup = null;
                    maxgroup = null;
                    glist = null;
                    gpair = null;
                }

                if (ngroups == 1)
                {
                    npmax = DELTAPART;
                    plist = new int[npmax];
                    //memory.create(plist, npmax, "collide:plist");
                }
                if (ngroups > 1)
                {
                    ngroup = new int[ngroups];
                    maxgroup = new int[ngroups];
                    glist = new int[ngroups][];
                    for (int i = 0; i < ngroups; i++)
                    {
                        maxgroup[i] = DELTAPART;
                        glist[i] = new int[DELTAPART];
                        //memory.create(glist[i], DELTAPART, "collide:glist");
                    }
                    gpair = new int[ngroups * ngroups, 3];
                    //memory.create(gpair, ngroups * ngroups, 3, "collide:gpair");
                }
            }

            // allocate vremax,remain if group count changed
            // will always be allocated on first run since oldgroups = 0
            // set vremax_intitial via values calculated by collide style

            if (ngroups != oldgroups)
            {
                //memory.destroy(vremax);
                //memory.destroy(vremax_initial);
                //memory.destroy(remain);
                nglocal = sparta.grid.nlocal;
                nglocalmax = nglocal;
                vremax = new double[nglocalmax][,];
                for (int i = 0; i < nglocalmax; i++)
                {
                    vremax[i] = new double[ngroups, ngroups];
                }
                
                vremax_initial = new double[ngroups, ngroups];
                //memory.create(vremax, nglocalmax, ngroups, ngroups, "collide:vremax");
                //memory.create(vremax_initial, ngroups, ngroups, "collide:vremax_initial");
                if (remainflag != 0)
                {
                    remain = new double[nglocalmax][,]; 
                    for (int i = 0; i < nglocalmax; i++)
                    {
                        remain[i] = new double[ngroups, ngroups];
                    }
                }
                //memory.create(remain, nglocalmax, ngroups, ngroups, "collide:remain");

                for (int igroup = 0; igroup < ngroups; igroup++)
                    for (int jgroup = 0; jgroup < ngroups; jgroup++)
                        vremax_initial[igroup,jgroup] = vremax_init(igroup, jgroup);
            }
            // if recombination reactions exist, set flags per species pair

            recombflag = 0;
            if (sparta.react!=null)
            {
                System.Console.WriteLine("Collide init react");
                //recombflag = sparta.react.recombflag;
                //recomb_boost_inverse = sparta.react.recomb_boost_inverse;
            }

            if (recombflag != 0)
            {
                int nspecies = sparta.particle.nspecies;
                //memory.destroy(recomb_ijflag);
                recomb_ijflag = new int[nspecies, nspecies];
                //memory.create(recomb_ijflag, nspecies, nspecies, "collide:recomb_ijflag");
                for (int i = 0; i < nspecies; i++)
                    for (int j = 0; j < nspecies; j++)
                        //recomb_ijflag[i,j] = sparta.react.recomb_exist(i, j);
                        System.Console.WriteLine("collide init react");
            }

            // find ambipolar fix
            // set ambipolar vector/array indices
            // if reactions defined, check that they are valid ambipolar reactions

            if (ambiflag!=0)
            {
                index_ionambi = sparta.particle.find_custom("ionambi");
                index_velambi = sparta.particle.find_custom("velambi");
                if (index_ionambi < 0 || index_velambi < 0)
                    sparta.error.all("Collision ambipolar without fix ambipolar");
                if (sparta.react != null)
                {
                    System.Console.WriteLine("collide init react");
                    //sparta.react.ambi_check();
                }

                int ifix;
                for (ifix = 0; ifix < sparta.modify.nfix; ifix++)
                    if (string.Equals(sparta.modify.fix[ifix].style, "ambipolar")) break;
                FixAmbipolar afix = (FixAmbipolar)sparta.modify.fix[ifix];
                ambispecies = afix.especies;
                ions = afix.ions;
            }

            // vre_next = next timestep to zero vremax & remain, based on vre_every

            if (vre_every!=0) vre_next = (sparta.update.ntimestep / vre_every) * vre_every + vre_every;
            else vre_next = sparta.update.laststep + 1;

            // if requested reset vremax & remain
            // must be after per-species vremax_initial is setup

            if (vre_first != 0 || vre_start != 0)
            {
                reset_vremax();
                vre_first = 0;
            }

            // initialize running stats before each run

            ncollide_running = nattempt_running = nreact_running = 0;
            
        }
        //public void modify_params(int, char**);
        public void reset_vremax()
        {
            for (int icell = 0; icell < nglocal; icell++)
                for (int igroup = 0; igroup < ngroups; igroup++)
                    for (int jgroup = 0; jgroup < ngroups; jgroup++)
                    {
                        vremax[icell][igroup,jgroup] = vremax_initial[igroup,jgroup];
                        if (remainflag!=0) remain[icell][igroup,jgroup] = 0.0;
                    }
        }
        public virtual void collisions()
        {
            // if requested, reset vrwmax & remain

            if (sparta.update.ntimestep == vre_next)
            {
                reset_vremax();
                vre_next += vre_every;
            }

            // counters

            ncollide_one = nattempt_one = nreact_one = 0;
            ndelete = 0;

            // perform collisions without or with ambipolar approximation
            // one variant is optimized for a single group

            if (ambiflag==0)
            {
                //Console.WriteLine("Collide.collisions().ambiflag==0");
                if (nearcp == 0)
                {
                    if (ngroups == 1) collisions_one(0);
                    else collisions_group(0);
                }
                else
                {
                    if (ngroups == 1) collisions_one (1);
                    else collisions_group(1);
                }
            }
            else
            {
                //if (ngroups == 1) collisions_one_ambipolar();
                //else collisions_group_ambipolar();
                Console.WriteLine("Collide.collisions().ambiflag!=0");
            }

            // remove any particles deleted in chemistry reactions
            // if reactions occurred, particles are no longer sorted
            // e.g. compress_reactions may have reallocated particle.next vector

            if (ndelete!=0) sparta.particle.compress_reactions(ndelete, dellist);
            if (sparta.react !=null) sparta.particle.sorted = 0;

            // accumulate running totals

            nattempt_running += nattempt_one;
            ncollide_running += ncollide_one;
            nreact_running += nreact_one;
            //Console.WriteLine("Collide virtual collisions");
        }

        public virtual double vremax_init(int igroup, int jgroup)
        {
            return 0;
        }
        public virtual double attempt_collision(int cell, int np, double volume)
        {
            Console.WriteLine("Collide virtual attempt_collision1");
            return 0;
        }
        public virtual double attempt_collision(int icell, int igroup, int jgroup, double volume)
        {
            Console.WriteLine("Collide virtual attempt_collision2");
            return 0;
        }
        public virtual int test_collision(int icell, int igroup, int jgroup,
                        Particle.OnePart ip, Particle.OnePart jp)
        {
            Console.WriteLine("Collide virtual test_collision");

            return 0;
        }
        public virtual void setup_collision(Particle.OnePart ip, Particle.OnePart jp)
        {
            Console.WriteLine("Collide virtual setup_collision");
        }
        public virtual int perform_collision(ref Particle.OnePart? ip, ref Particle.OnePart? jp,
                                       ref Particle.OnePart? kp)
        {
            Console.WriteLine("Collide virtual perform_collision");
            return 0;
        }

        //public virtual double extract(int, const char*) {return 0.0;


        public virtual int pack_grid_one(int icell,ref StringBuilder buf, int memflag)
        {
            int nbytes = ngroups * ngroups * sizeof(double);
            
            Grid.ChildCell[] cells = sparta.grid.cells;

            int n;
            if (remainflag!=0)
            {
                if (memflag!=0)
                {
                    //memcpy(buf, &vremax[icell][0][0], nbytes);
                    buf.Append(vremax[icell][0,0]);
                    //memcpy(&buf[nbytes], &remain[icell][0][0], nbytes);
                    buf.Append(remain[icell][0,0]);
                }
                n = 2 * nbytes;
            }
            else
            {
                if (memflag!=0)
                {
                    //memcpy(buf, &vremax[icell][0][0], nbytes);
                    buf.Append(vremax[icell][0,0]);
                }

                n = nbytes;
            }

            if (cells[icell].nsplit > 1)
            {
                int isplit = cells[icell].isplit;
                int nsplit = cells[icell].nsplit;
                for (int i = 0; i < nsplit; i++)
                {
                    int m = sparta.grid.sinfo[isplit].csubs[i];
                    if (remainflag!=0)
                    {
                        if (memflag!=0)
                        {
                            //memcpy(&buf[n], &vremax[m][0][0], nbytes);
                            buf.Append(vremax[m][0,0]);
                            n += nbytes;
                            //memcpy(&buf[n], &remain[m][0][0], nbytes);
                            buf.Append(remain[m][0,0]);
                            n += nbytes;
                        }
                        else n += 2 * nbytes;
                    }
                    else
                    {
                        if (memflag!=0)
                        {
                            //memcpy(&buf[n], &vremax[m][0][0], nbytes);
                            buf.Append(vremax[m][ 0, 0]);
                        }

                        n += nbytes;
                    }
                }
            }

            return n;
        }
        //public virtual int unpack_grid_one(int, char*);
        public virtual void compress_grid()
        {

            int nbytes = ngroups * ngroups * sizeof(double);

            int me = sparta.comm.me;
            Grid.ChildCell[] cells = sparta.grid.cells;

            // keep an unsplit or split cell if staying on this proc
            // keep a sub cell if its split cell is staying on this proc

            int ncurrent = nglocal;
            nglocal = 0;
            for (int icell = 0; icell < ncurrent; icell++)
            {
                if (cells[icell].nsplit >= 1)
                {
                    if (cells[icell].proc != me) continue;
                }
                else
                {
                    int isplit = cells[icell].isplit;
                    if (cells[sparta.grid.sinfo[isplit].icell].proc != me) continue;
                }

                if (nglocal != icell)
                {
                    vremax[nglocal] = vremax[icell];
                    //memcpy(&vremax[nglocal][0][0], &vremax[icell][0][0], nbytes);
                    if (remainflag!=0 )
                    {
                        remain[nglocal] = remain[icell];
                       //memcpy(&remain[nglocal][0][0], &remain[icell][0][0], nbytes);
                    }
                }
                nglocal++;
            }
        }
        //public virtual void adapt_grid();

        protected int npmax;          // max # of particles in plist
        protected int[] plist;         // list of particles in a single group
         
        protected int nglocal;        // current size of per-cell arrays
        protected int nglocalmax;     // max allocated size of per-cell arrays
         
        protected int ngroups;        // # of groups
        protected int[] ngroup;        // # of particles in one cell of each group
        protected int[] maxgroup;      // max # of glist indices allocated per group
        protected int[][] glist;        // indices of particles in one cell of each group
         
        protected int npair;          // # of group pairs to do collisions for
        protected int[,] gpair;        // Nx3 list of species pairs to do collisions for
                             // 0 = igroup, 1 = jgroup, 2 = # of attempt collisions
         
        protected int max_nn;             // allocated size of nn_last_partner
        protected int[] nn_last_partner;   // index+1 of last collision partner for each particle
                                 // 0 = no collision yet (on this step)
        protected int[] nn_last_partner_igroup;   // ditto for igroup and jgroup particles
        protected int[] nn_last_partner_jgroup;
         
        protected int ndelete, maxdelete;      // # of particles removed by chemsitry
        protected int[] dellist;               // list of particle indices to delete

        protected string mixID;               // ID of mixture to use for groups
        protected Mixture  mixture;    // ptr to mixture
        protected RanPark  random;     // RNG for collision generation

        protected int vre_first;      // 1 for first run after collision style is defined
        protected int vre_start;      // 1 if reset vre params at start of each run
        protected int vre_every;      // reset vre params every this many steps
        protected bigint vre_next;    // next timestep to reset vre params on
        protected int remainflag;     // 1 if remain defined, else use random fraction
         
        protected double[][,] vremax;   // max relative velocity, per cell, per group pair
        protected double[][,] remain;   // collision number remainder, per cell, per group pair
        protected double[,] vremax_initial;   // initial vremax value, per group pair

        // recombination reactions

        protected int recombflag;               // 1 if recomb reactions enabled, 0 if not
        protected double recomb_boost_inverse;  // recombination rate boost factor from React
        protected int[,] recomb_ijflag;          // 1 if species I,J have recomb reaction(s)

        protected int oldgroups;         // pass from parent to child class
        protected int copymode;          // 1 if copy of class (prevents deallocation of
                                //  base class when child copy is destroyed)
        protected int kokkosable;        // 1 if collide method supports Kokkos

        // ambipolar approximation data structs

        protected int ambiflag;       // 1 if ambipolar option is enabled
        protected int ambispecies;    // species for ambipolar electrons
        protected int index_ionambi;  // 2 custom ambipolar vectors
        protected int index_velambi;
        protected int[] ions;          // ptr to fix ambipolar list of ions

        protected int nelectron;                // # of ambipolar electrons in elist
        protected int maxelectron;              // max # elist can hold
        protected List<Particle.OnePart> elist;     // list of ambipolar electrons
                                                    // for one grid cell or pair of groups in cell

        //inline void addgroup(int igroup, int n)
        //{
        //    if (ngroup[igroup] == maxgroup[igroup])
        //    {
        //        maxgroup[igroup] += DELTAPART;
        //        memory.grow(glist[igroup], maxgroup[igroup], "collide:grouplist");
        //    }
        //    glist[igroup][ngroup[igroup]++] = n;
        //}

        public void collisions_one(int NEARCP)
        {
            int i, j, k, m, n, ip, np;
            int nattempt, reactflag;
            double attempt, volume;
            Particle.OnePart? ipart=null,jpart=null,kpart=null;

            // loop over cells I own

            Grid.ChildInfo[] cinfo = sparta.grid.cinfo;

            Particle.OnePart[] particles = sparta.particle.particles;
            int[] next = sparta.particle.next;

            for (int icell = 0; icell < nglocal; icell++)
            {
                np = cinfo[icell].count;
                if (np <= 1) continue;

                if (NEARCP!=0)
                {
                    if (np > max_nn)
                    {
                        //realloc_nn(np, nn_last_partner);
                        nn_last_partner = new int[np];
                    }

                    //memset(nn_last_partner, 0, np * sizeof(int));
                }

                ip = cinfo[icell].first;
                volume = cinfo[icell].volume / cinfo[icell].weight;
                if (volume == 0.0) sparta.error.one( "Collision cell volume is zero");

                // setup particle list for this cell

                if (np > npmax)
                {
                    npmax = np + DELTAPART;
                    //memory.destroy(plist);
                    //memory.create(plist, npmax, "collide:plist");
                    plist = new int[npmax];
                }

                n = 0;
                while (ip >= 0)
                {
                    plist[n++] = ip;
                    ip = next[ip];
                }

                // attempt = exact collision attempt count for a pair of groups
                // nattempt = rounded attempt with RN

                attempt = attempt_collision(icell, np, volume);
                nattempt = Convert.ToInt32(attempt);

                if (nattempt==0) continue;
                nattempt_one += nattempt;

                // perform collisions
                // select random pair of particles, cannot be same
                // test if collision actually occurs

                for (m = 0; m < nattempt; m++)
                {
                    i = (int)(np * random.uniform());
                    if (NEARCP!=0) j = find_nn(i, np);
                    else
                    {
                        j = (int)(np * random.uniform());
                        while (i == j) j = (int)(np * random.uniform());
                    }

                    ipart = particles[plist[i]];
                    jpart = particles[plist[j]];

                    // test if collision actually occurs
                    // continue to next collision if no reaction

                    if (test_collision(icell, 0, 0, ipart.Value, jpart.Value)==0) continue;

                    if (NEARCP!=0)
                    {
                        nn_last_partner[i] = j + 1;
                        nn_last_partner[j] = i + 1;
                    }

                    // if recombination reaction is possible for this IJ pair
                    // pick a 3rd particle to participate and set cell number density
                    // unless boost factor turns it off, or there is no 3rd particle

                    if (recombflag!=0 && recomb_ijflag[ipart.Value.ispecies,jpart.Value.ispecies]!=0)
                    {
                        //if (random.uniform() > sparta.react.recomb_boost_inverse)
                        //    sparta.react.recomb_species = -1;
                        //else if (np <= 2)
                        //    sparta.react.recomb_species = -1;
                        //else
                        //{
                        //    k = (int)(np * random.uniform());
                        //    while (k == i || k == j) k = (int)(np * random.uniform());
                        //    sparta.react.recomb_part3 = particles[plist[k]];
                        //    sparta.react.recomb_species = sparta.react.recomb_part3.ispecies;
                        //    sparta.react.recomb_density = np * sparta.update.fnum / volume;
                        //}
                        Console.WriteLine("Collide collisions_one()-> react");
                    }

                    // perform collision and possible reaction
                    
                    setup_collision(ipart.Value, jpart.Value);
                    reactflag = perform_collision(ref ipart,ref jpart,ref kpart);
                    ncollide_one++;
                    if (reactflag!=0) nreact_one++;
                    else continue;

                    // if jpart destroyed, delete from plist
                    // also add particle to deletion list
                    // exit attempt loop if only single particle left

                    if (jpart==null)
                    {
                        if (ndelete == maxdelete)
                        {
                            maxdelete += DELTADELETE;
                            //memory.grow(dellist, maxdelete, "collide:dellist");
                            dellist = new int[maxdelete];
                        }
                        dellist[ndelete++] = plist[j];
                        np--;
                        plist[j] = plist[np];
                        if (NEARCP!=0) nn_last_partner[j] = nn_last_partner[np];
                        if (np < 2) break;
                    }

                    // if kpart created, add to plist
                    // kpart was just added to particle list, so index = nlocal-1
                    // particle data structs may have been realloced by kpart

                    if (kpart!=null)
                    {
                        if (np == npmax)
                        {
                            npmax = np + DELTAPART;
                            //memory.grow(plist, npmax, "collide:plist");
                            plist = new int[npmax];
                        }
                        if (NEARCP!=0) set_nn(np);
                        plist[np++] = sparta.particle.nlocal - 1;
                        particles = sparta.particle.particles;
                    }
                }
            }
        }
        void collisions_group(int NEARCP)
        {
            Console.WriteLine("Collide  collisions_group");
        }
        //void collisions_one_ambipolar();
        //void collisions_group_ambipolar();
        //void ambi_reset(int, int, int, int, Particle.OnePart, Particle.OnePart,
        //                Particle.OnePart, int[]);
        //void ambi_check();
        //void grow_percell(int);

        int find_nn(int i, int np)
        {
            int jneigh;
            double dx, dy, dz, rsq;
            double[] xj;

            // if np = 2, just return J = non-I particle
            // np is never < 2

            if (np == 2) return (i + 1) % 2;

            Particle.OnePart ipart,jpart;
            Particle.OnePart[] particles = sparta.particle.particles;
            double dt = sparta.update.dt;

            // thresh = distance particle I moves in this timestep

            ipart = particles[plist[i]];
            double[] vi = ipart.v;
            double[] xi = ipart.x;
            double threshsq = dt * dt * (vi[0] * vi[0] + vi[1] * vi[1] + vi[2] * vi[2]);
            double minrsq = BIG;

            // nlimit = max # of J candidates to consider

            int nlimit = Math.Min(nearlimit, np - 1);
            int count = 0;

            // pick a random starting J
            // jneigh = collision partner when exit loop
            //   set to initial J as default in case no Nlimit J meets criteria

            int j = (int)(np * random.uniform());
            while (i == j) j = (int)(np * random.uniform());
            jneigh = j;

            while (count < nlimit)
            {
                count++;

                // skip this J if I,J last collided with each other

                if (nn_last_partner[i] == j + 1 && nn_last_partner[j] == i + 1)
                {
                    j++;
                    if (j == np) j = 0;
                    continue;
                }

                // rsq = squared distance between particles I and J
                // if rsq = 0.0, skip this J
                //   could be I = J, or a cloned J at same position as I
                // if rsq <= threshsq, this J is collision partner
                // if rsq = smallest yet seen, this J is tentative collision partner

                jpart = particles[plist[j]];
                xj = jpart.x;
                dx = xi[0] - xj[0];
                dy = xi[1] - xj[1];
                dz = xi[2] - xj[2];
                rsq = dx * dx + dy * dy + dz * dz;

                if (rsq > 0.0)
                {
                    if (rsq <= threshsq)
                    {
                        jneigh = j;
                        break;
                    }
                    if (rsq < minrsq)
                    {
                        minrsq = rsq;
                        jneigh = j;
                    }
                }

                j++;
                if (j == np) j = 0;
            }

            return jneigh;
        }
        //int find_nn_group(int, int[], int, int[], int[], int[]);
        //void realloc_nn(int, int[]&);
        void set_nn(int n)
        {
            if (n == max_nn)
            {
                max_nn *= 2;
                //memory->grow(nn_last_partner, max_nn, "collide:nn_last_partner");
                nn_last_partner = new int[max_nn];
            }
            nn_last_partner[n] = 0;
        }
        //void set_nn_group(int);

    }
}
