using System;

namespace cstest
{
    public class Mixture
    {
        public const int DELTA = 8;

        public string id;                   // ID of mixture
        public int nspecies;               // # of species in mixture
        public int[] species;               // species[i] = particle species index of 
                                     //              mixture species I
         
        public int ngroup;                 // # of defined groups
        public string[] groups;              // group IDs
        public int[] mix2group;             // m2g[i] = group that mixture species I is in
         
         // global attributes
        public double nrho;                // number density
        public int nrho_flag;              // 1 if user set nrho
        public double nrho_user;           // user value
        public double[] vstream=new double[3];          // stream velocity
        public int vstream_flag;           // 1 if user set vstream
        public double[] vstream_user=new double[3];     // user value
        public double temp_thermal;        // thermal temperature
        public double temp_rot;            // rotational temperature
        public double temp_vib;            // vibrational temperature
        public int temp_thermal_flag;      // 1 if user set thermal temp
        public int temp_rot_flag;          // 1 if user set rotational temp
        public int temp_vib_flag;          // 1 if user set vibrational temp
        public double temp_thermal_user;   // user value
        public double temp_rot_user;       // user value
        public double temp_vib_user;       // user value
         
         // per-species attributes
        public double[] fraction;           // relative fraction of each species
        public int[] fraction_flag;         // 1 if user set fraction for a species
        public double[] fraction_user;      // user fractional value
         
         // set by init()
         
        public double[] cummulative;        // cummulative fraction for each species
        public int[] groupsize;             // # of species in each group
        public int[][] groupspecies;         // list of particle species indices in each group
        public int[] species2group;         // s2g[i] = group that particle species I is in
                                     // -1 if species I not in mixture
        public int[] species2species;       // s2s[i] = mixture species that 
                                     //   particle species I is
                                     // -1 if species I not in mixture
        public double[] vscale;             // pre-computed velocity scale factor

        public void copy(Mixture old)
        {
            nrho = old.nrho;
            nrho_flag = old.nrho_flag;
            nrho_user = old.nrho_user;
            vstream[0] = old.vstream[0];
            vstream[1] = old.vstream[1];
            vstream[2] = old.vstream[2];
            vstream_flag = old.vstream_flag;
            vstream_user[0] = old.vstream_user[0];
            vstream_user[1] = old.vstream_user[1];
            vstream_user[2] = old.vstream_user[2];
            temp_thermal = old.temp_thermal;
            temp_thermal_flag = old.temp_thermal_flag;
            temp_thermal_user = old.temp_thermal_user;

            nspecies = maxspecies = old.nspecies;
            allocate();

            for (int i = 0; i < nspecies; i++)
            {
                species[i] = old.species[i];
                fraction[i] = old.fraction[i];
                fraction_flag[i] = old.fraction_flag[i];
                fraction_user[i] = old.fraction_user[i];
                mix2group[i] = old.mix2group[i];
            }

            for (int i = 0; i < old.ngroup; i++)
                add_group(old.groups[i]);
        }
        public void command(int narg, string[] arg)
        {
            if (narg < 1) sparta.error.all( "Illegal mixture command");

            // nsp = # of species listed before optional keywords
            // iarg = start of optional keywords

            int iarg;
            for (iarg = 1; iarg < narg; iarg++)
            {
                if (string.Equals(arg[iarg], "nrho") ) break;
                if (string.Equals(arg[iarg], "vstream") ) break;
                if (string.Equals(arg[iarg], "temp") ) break;
                if (string.Equals(arg[iarg], "trot") ) break;
                if (string.Equals(arg[iarg], "tvib") ) break;
                if (string.Equals(arg[iarg], "frac") ) break;
                if (string.Equals(arg[iarg], "group")) break;
                if (string.Equals(arg[iarg], "copy") ) break;
                if (string.Equals(arg[iarg], "delete") ) break;
            }
            int nsp = iarg - 1;

            // add_species() processes list of species
            // params() processes remaining optional keywords
            string[] arg1;
            string[] arg2;
            string[] arg3;
            if (arg.Length>1)
            {
                arg1 = new string[arg.Length - 1];
                arg2 = new string[arg.Length - iarg];
                arg3 = new string[arg.Length - iarg - copyarg];
                Array.Copy(arg, 1, arg1, 0, arg1.Length);
                Array.Copy(arg, iarg, arg2, 0, arg2.Length );
                Array.Copy(arg, iarg + copyarg, arg3, 0, arg3.Length);
                add_species(nsp, arg1);
                myparams(narg - iarg, arg2);
            }
            else
            {
                arg3 = new string[arg.Length];
            }




            // if copy keyword was used, create a new mixture via add_mixture()
            // then invoke its copy() method, passing it this mixture

            if (copyflag!=0)
            {
                sparta.particle.add_mixture(1, arg3);
                sparta.particle.mixture[sparta.particle.nmixture - 1].copy(this);
            }
        }
        //public void init();
        //public int init_fraction(int[], double[], double[], double[]);
        public void add_species_default(string name)
        {
            int index = sparta.particle.find_species(name);
            if (nspecies == maxspecies) allocate();
            species[nspecies] = index;

            if (all_default!=0 && ngroup == 0) add_group("all");
            if (species_default!=0) add_group(name);
            mix2group[nspecies] = ngroup - 1;

            nspecies++;
        }
        public int find_group(string idgroup)
        {
            for (int i = 0; i < ngroup; i++)
                if (string.Equals(groups[i], idgroup)) return i;
            return -1;
        }
        //public void write_restart(FILE* fp);
        //public void read_restart(FILE* fp);
        private SPARTA sparta;

        public Mixture(SPARTA sparta,string userid)
        {
            this.sparta = sparta;
            // mixture ID must be all alphanumeric chars or underscores
            int n = userid.Length + 1;
            id =string.Copy(userid);

            for (int i = 0; i < n - 1; i++)
                if (!char.IsLetterOrDigit(id[i]) && id[i] != '_')
                    sparta.error.all(
                       "Mixture ID must be alphanumeric or underscore characters");

            // special default mixtures

            all_default = species_default = 0;
            if (string.Equals(id, "all") ) all_default = 1;
            if (string.Equals(id, "species") ) species_default = 1;

            // initialize mixture values

            nspecies = maxspecies = 0;
            species = null;

            nrho_flag = 0;
            vstream_flag = 0;
            temp_thermal_flag = 0;
            temp_rot_flag = 0;
            temp_vib_flag = 0;

            fraction = null;
            fraction_user = null;
            fraction_flag = null;
            cummulative = null;

            ngroup = maxgroup = 0;
            groups = null;
            groupsize = null;
            groupspecies = null;

            mix2group = null;
            species2group = null;
            species2species = null;

            vscale = null;
            active = null;

            allocate();
        }


        private int maxspecies, maxgroup;
        private int copyflag, copyarg;
         
        private int activeflag;             // 1 if species are listed in mixture command
        private int[] active;                // flags for species listed in mixture command
        private int all_default;            // 1 if this is default mixture "all"
        private int species_default;        // 1 if this is default mixture "species"

        private void add_species(int narg, string[] arg)
        {
            int i, j, index;

            // activeflag = 1 if species are listed
            // active[i] = 0 if current mixture species I is not in the list
            // active[i] = 1 if species I is in the list but already existed in mixture
            // active[i] = 2 if species I was just added b/c it was in the list
            // active[i] = 3 if species is being removed from mixture via delete keyword
            //             this flag is set in params()

            if (narg!=0) activeflag = 1;
            else activeflag = 0;
            for (i = 0; i < nspecies; i++) active[i] = 0;

            for (i = 0; i < narg; i++)
            {
                index = sparta.particle.find_species(arg[i]);
                if (index < 0) sparta.error.all( "Mixture species is not defined");
                for (j = 0; j < nspecies; j++)
                    if (species[j] == index) break;
                if (j < nspecies) active[j] = 1;
                else
                {
                    if (all_default!=0 || species_default!=0)
                        sparta.error.all( "Cannot add new species to mixture all or species");
                    if (nspecies == maxspecies) allocate();
                    active[nspecies] = 2;
                    species[nspecies++] = index;
                }
            }
        }
        private void myparams(int narg, string[] arg)
        {
            // for global attributes, set immediately
            // for per-species attributes, store flags

            copyflag = 0;
            int deleteflag = 0;
            int fracflag = 0;
            int groupflag = 0;
            double fracvalue=0;
            int grouparg=0;

            int iarg = 0;

            while (iarg<narg)
            {
                switch (arg[iarg])
                {
                    case "nrho":
                        if (iarg + 2 > narg) sparta.error.all("Illegal mixture command");
                        nrho_flag = 1;
                        nrho_user = double.Parse(arg[iarg + 1]);
                        if (nrho_user <= 0.0) sparta.error.all("Illegal mixture command");
                        iarg += 2;
                        break;
                    case "vstream":
                        if (iarg + 4 > narg) sparta.error.all("Illegal mixture command");
                        vstream_flag = 1;
                        vstream_user[0] = double.Parse(arg[iarg + 1]);
                        vstream_user[1] = double.Parse(arg[iarg + 2]);
                        vstream_user[2] = double.Parse(arg[iarg + 3]);
                        iarg += 4;
                        break;
                    case "temp":
                        if (iarg + 2 > narg) sparta.error.all("Illegal mixture command");
                        temp_thermal_flag = 1;
                        temp_thermal_user = double.Parse(arg[iarg + 1]);
                        if (temp_thermal_user <= 0.0)
                            sparta.error.all("Illegal mixture command");
                        iarg += 2;
                        break;
                    case "trot":
                        if (iarg + 2 > narg) sparta.error.all("Illegal mixture command");
                        temp_rot_flag = 1;
                        temp_rot_user = double.Parse(arg[iarg + 1]);
                        if (temp_rot_user <= 0.0)
                            sparta.error.all("Illegal mixture command");
                        iarg += 2;
                        break;
                    case "tvib":
                        if (iarg + 2 > narg) sparta.error.all("Illegal mixture command");
                        temp_vib_flag = 1;
                        temp_vib_user = double.Parse(arg[iarg + 1]);
                        if (temp_vib_user <= 0.0)
                            sparta.error.all("Illegal mixture command");
                        iarg += 2;
                        break;
                    case "frac":
                        if (iarg + 2 > narg) sparta.error.all("Illegal mixture command");
                        fracflag = 1;
                        fracvalue = double.Parse(arg[iarg + 1]);
                        if (fracvalue < 0.0 || fracvalue > 1.0)
                            sparta.error.all("Illegal mixture command");
                        iarg += 2;
                        break;
                    case "group":
                        if (iarg + 2 > narg) sparta.error.all("Illegal mixture command");
                        groupflag = 1;
                        grouparg = iarg + 1;
                        int n = arg[grouparg].Length;
                        for (int i = 0; i < n; i++)
                            if (!char.IsLetterOrDigit(arg[grouparg][i]) && arg[grouparg][i] != '_')
                                sparta.error.all("Mixture group ID must be alphanumeric or underscore characters");
                        if (all_default!=0 || species_default!=0)
                            sparta.error.all("Cannot use group keyword with mixture all or species");
                        iarg += 2;
                        break;
                    case "copy":
                        if (iarg + 2 > narg) sparta.error.all("Illegal mixture command");
                        if (sparta.particle.find_mixture(arg[iarg + 1]) >= 0)
                            sparta.error.all("New mixture copy mixture already exists");
                        copyflag = 1;
                        copyarg = iarg + 1;
                        iarg += 2;
                        break;
                    case "delete":
                        if (iarg + 1 > narg) sparta.error.all("Illegal mixture command");
                        deleteflag = 1;
                        if (activeflag!=0)
                            sparta.error.all("Mixture delete cannot list species before keyword");
                        if (iarg != 0)
                            sparta.error.all("Mixture delete must be only keyword");
                        iarg = narg;
                        break;

                    default:
                        sparta.error.all( "Illegal mixture command");
                        break;
                }
            }
            if (deleteflag!=0)
            {
                int m, ispecies;
                for (iarg = 1; iarg < narg; iarg++)
                {
                    ispecies = sparta.particle.find_species(arg[iarg]);
                    if (ispecies < 0)
                        sparta.error.all( "Mixture delete species is not recognized");
                    for (m = 0; m < nspecies; m++)
                        if (species[m] == ispecies) break;
                    if (m == nspecies)
                        sparta.error.all( "Mixture delete species is not in mixture");
                    active[m] = 1;
                }

                int nspecies_original = nspecies;
                m = 0;
                for (int i = 0; i < nspecies_original; i++)
                {
                    species[m] = species[i];
                    mix2group[m] = mix2group[i];
                    if (active[i]!=0) nspecies--;
                    else m++;
                }
            }
            // assign per-species attributes

            if (fracflag!=0)
            {
                for (int i = 0; i < nspecies; i++)
                    if (active[i]!=0)
                    {
                        fraction_flag[i] = 1;
                        fraction_user[i] = fracvalue;
                    }
            }
            // assign species to groups
            // end up with:
            //   every species assigned to exactly one group via mix2group
            //   no empty groups via shrink_groups()
            // if group-ID = SELF:
            //   no listed species: delete groups, assign each species to own group
            //   listed species: assign each listed species to own group
            // else if group-ID = user name:
            //   no listed species: delete groups, assign all species to group-ID
            //   listed species: assign listed species to group-ID
            // else if group keyword not specified:
            //   assign any listed species that are new each to group "default"


            if (groupflag!=0)
            {
                if (string.Equals(arg[grouparg], "SELF"))
                {
                    if (activeflag==0)
                    {
                        delete_groups();
                        for (int i = 0; i < nspecies; i++)
                        {
                            add_group(sparta.particle.species[species[i]].id);
                            mix2group[i] = ngroup - 1;
                        }
                    }
                    else
                    {
                        for (int i = 0; i < nspecies; i++)
                        {
                            if (active[i]==0) continue;
                            int igroup = find_group(sparta.particle.species[species[i]].id);
                            if (igroup < 0)
                            {
                                add_group(sparta.particle.species[species[i]].id);
                                igroup = ngroup - 1;
                            }
                            mix2group[i] = igroup;
                        }
                    }

                }
                else
                {
                    if (activeflag==0)
                    {
                        delete_groups();
                        add_group(arg[grouparg]);
                        for (int i = 0; i < nspecies; i++) mix2group[i] = ngroup - 1;
                    }
                    else
                    {
                        int igroup = find_group(arg[grouparg]);
                        if (igroup < 0)
                        {
                            add_group(arg[grouparg]);
                            igroup = ngroup - 1;
                        }
                        for (int i = 0; i < nspecies; i++)
                            if (active[i]!=0) mix2group[i] = igroup;
                    }
                }

            }
            else if (activeflag!=0)
            {
                for (int i = 0; i < nspecies; i++)
                {
                    if (active[i] != 2) continue;
                    int igroup = find_group("default");
                    if (igroup < 0)
                    {
                        add_group("default");
                        igroup = ngroup - 1;
                    }
                    mix2group[i] = igroup;
                }
            }

            // remove empty groups due to deleteflag or groupflag operations

            shrink_groups();



        }
        private void allocate()
        {
            int old = maxspecies;
            maxspecies += DELTA;
            species = new int[maxspecies];
            fraction = new double[maxspecies];
            fraction_flag = new int[maxspecies];
            fraction_user = new double[maxspecies];
            cummulative = new double[maxspecies];
            mix2group = new int[maxspecies];
            vscale = new double[maxspecies];
            active = new int[maxspecies];


            for (int i = old; i < maxspecies; i++)
            {
                fraction_flag[i] = 0;
                fraction_user[i] = 0.0;
            }
        }
        private void delete_groups()
        {
            //for (int i = 0; i < ngroup; i++) delete[] groups[i];
            ngroup = 0;
        }
        private void shrink_groups()
        {
            int i, nsp;

            int igroup = 0;
            while (igroup < ngroup)
            {
                nsp = 0;
                for (i = 0; i < nspecies; i++)
                    if (mix2group[i] == igroup) nsp++;
                if (nsp == 0)
                {
                    //delete[] groups[igroup];
                    for (i = igroup; i < ngroup - 1; i++)
                        groups[i] = groups[i + 1];
                    for (i = 0; i < nspecies; i++)
                        if (mix2group[i] > igroup) mix2group[i]--;
                    ngroup--;
                }
                else igroup++;
            }
        }
        private void add_group(string idgroup)
        {
            if (ngroup == maxgroup)
            {
                maxgroup += DELTA;
                groups = new string[maxgroup];

            }

            int n = idgroup.Length + 1;
            groups[ngroup] = string.Copy(idgroup);
            ngroup++;
        }

    }
}