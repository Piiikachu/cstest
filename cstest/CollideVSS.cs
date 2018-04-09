using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    class CollideVSS : Collide
    {
        enum Enum1 { NONE, DISCRETE, SMOOTH };            // several files
        enum Enum2 { CONSTANT, VARIABLE };

        public const int MAXLINE = 1024;
        //public virtual void init();

        //public double vremax_init(int, int);
        //public virtual double attempt_collision(int, int, double);
        //public double attempt_collision(int, int, int, double);
        //public virtual int test_collision(int, int, int, Particle::OnePart*, Particle::OnePart*);
        //public virtual void setup_collision(Particle::OnePart*, Particle::OnePart*);
        //public virtual int perform_collision(Particle::OnePart*&, Particle::OnePart*&,
        //                       Particle::OnePart*&);
        //public double extract(int, const char*);
        private SPARTA sparta;
        public CollideVSS(SPARTA sparta, int narg, string[] arg) : base(sparta, narg, arg)
        {
            this.sparta = sparta;
            if (narg < 3) sparta.error.all("Illegal collide command");

            // optional args

            relaxflag = (int)Enum2.CONSTANT;

            int iarg = 3;
            while (iarg < narg)
            {
                if (string.Compare(arg[iarg], "relax") == 0)
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal collide command");
                    if (string.Compare(arg[iarg + 1], "constant") == 0) relaxflag = (int)Enum2.CONSTANT;
                    else if (string.Compare(arg[iarg + 1], "variable") == 0) relaxflag = (int)Enum2.VARIABLE;
                    else sparta.error.all("Illegal collide command");
                    iarg += 2;
                }
                else sparta.error.all("Illegal collide command");
            }

            // proc 0 reads file to extract params for current species
            // broadcasts params to all procs
            nparams = sparta.particle.nspecies;
            if (nparams == 0)
                sparta.error.all("Cannot use collide command with no species defined");

            mparams = new List<Params>(nparams);
            if (sparta.comm.me == 0) read_param_file(arg[2]);
            //sparta.mpi.MPI_Bcast(mparams, nparams * sizeof(Params), MPI_BYTE, 0, world);

            for (int i = 0; i < nparams; i++)
            {
                if (mparams[i].diam < 0.0)
                {
                    string str = string.Format("Species {0} did not appear in VSS parameter file",
                        sparta.particle.species[i].id);
                    sparta.error.one(str);

                }
            }
            prefactor = new double[nparams, nparams];
        }
        public struct State
        {      // two-particle state
            public double vr2;
            public double vr;
            public double imass, jmass;
            public double mr;
            public double ave_rotdof;
            public double ave_vibdof;
            public double ave_dof;
            public double etrans;
            public double erot;
            public double evib;
            public double eexchange;
            public double eint;
            public double etotal;
            public double ucmf;
            public double vcmf;
            public double wcmf;
        }

        public struct Params
        {             // VSS model parameters
            public double diam;
            public double omega;
            public double tref;
            public double alpha;
            public double rotc1;
            public double rotc2;
            public double rotc3;
            public double vibc1;
            public double vibc2;
        }
        protected int relaxflag, eng_exchange;
        protected double vr_indice;
        protected double[,] prefactor; // static portion of collision attempt frequency

        protected State precoln;       // state before collision
        protected State postcoln;      // state after collision

        protected List<Params> mparams;             // VSS params for each species
        protected int nparams;                // # of per-species params read in

        //protected void SCATTER_TwoBodyScattering(Particle::OnePart*,
        //                Particle::OnePart*);
        //protected void EEXCHANGE_NonReactingEDisposal(Particle::OnePart*,
        //                     Particle::OnePart*);
        //protected void SCATTER_ThreeBodyScattering(Particle::OnePart*,
        //                                  Particle::OnePart*,
        //                                  Particle::OnePart*);
        //protected void EEXCHANGE_ReactingEDisposal(Particle::OnePart*,
        //                                 Particle::OnePart*,
        //                                 Particle::OnePart*);

        //protected double sample_bl(RanPark*, double, double);
        //protected double rotrel(int, double);
        //protected double vibrel(int, double);

        protected void read_param_file(string fname)
        {
            FileStream fp = new FileStream(fname, FileMode.Open, FileAccess.Read);
            if (fp == null)
            {
                string str = string.Format("Cannot open VSS parameter file {0}", fname);
                sparta.error.one(str);
            }
            // set all diameters to -1, so can detect if not read

            for (int i = 0; i < nparams; i++)
            {
                Params mparam = new Params();
                if (mparams.Count!=nparams)
                {
                    mparam.diam = -1.0;
                    mparams.Add(mparam);
                }
                else
                {
                    mparam.diam = -1.0;

                    mparams[i] = mparam;
                }
                
                
                
            }

            // read file line by line
            // skip blank lines or comment lines starting with '#'
            // all other lines must have at least NWORDS 


            //todo: readvssfile()
            Console.WriteLine("read vss file");
        }
        //protected int wordcount(char*);
        //protected void wordparse(int, char*, char**);
    }
}
