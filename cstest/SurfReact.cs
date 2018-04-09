using System;
using System.Collections.Generic;
using System.IO;

namespace cstest
{
    public class SurfReact
    {
        enum Enum1{ DISSOCIATION, EXCHANGE, RECOMBINATION };        // other surf react files
        enum Enum2 { SIMPLE };                                     // other surf react files

        public const int MAXREACTANT = 1;
        public const int MAXPRODUCT = 2;
        public const int MAXCOEFF = 2;

        public const int MAXLINE = 1024;
        public const int DELTALIST = 16;


        public string id;
        public string style;
         
        public int vector_flag;          // 0/1 if compute_vector() function exists
        public int size_vector;          // length of global vector
        private SPARTA sparta;
         
        public SurfReact(SPARTA sparta, int narg, string[] args)
        {
            this.sparta = sparta;
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);

            // ID and style
            // ID must be all alphanumeric chars or underscores

            int n = arg[0].Length + 1;
            id = string.Copy( arg[0]);

            for (int i = 0; i < n - 1; i++)
                if (!char.IsLetterOrDigit(id[i]) && id[i] != '_')
                    sparta.error.all("Surf_react ID must be alphanumeric or underscore characters");

            n = arg[0].Length + 1;
            style = string.Copy(arg[0]);

            vector_flag = 1;
            size_vector = 2;

            nsingle = ntotal = 0;

            // surface reaction data structs

            nlist = maxlist = 0;
            rlist = null;

            reactions = null;
            indices = null;
        }
        //public virtual void init();
        //public virtual int react(Particle::OnePart*&, double*, Particle::OnePart*&) = 0;

        //public void tally_update();
        //public double compute_vector(int i);

        protected FileStream fp;
        protected int nsingle, ntotal;
        protected double[] one=new double[2], all=new double[2];
         
        protected struct OneReaction
        {
            int active;                    // 1 if reaction is active
            int type;                      // reaction type = DISSOCIATION, etc
            int style;                     // reaction style = ARRHENIUS, etc
            int ncoeff;                    // # of numerical coeffs
            int nreactant, nproduct;        // # of reactants and products
            string[] id_reactants,id_products;  // species IDs of reactants/products
            int[] reactants,products;      // species indices of reactants/products
            double coeff;                 // numerical coeffs for reaction
        };

        protected List<OneReaction> rlist;              // list of all reactions read from file
        protected int nlist;                       // # of reactions read from file
        protected int maxlist;                     // max # of reactions in rlist

        // possible reactions a reactant species is part of

        protected struct ReactionI
        {
            int[] list;           // list of indices into rlist, ptr into indices
            int n;               // # of reactions in list
        };

        protected List<ReactionI> reactions;       // reactions for all species
        protected int[] indices;               // master list of indices

        //protected void init_reactions();
        //protected void readfile(char*);
        //protected int readone(char*, char*, int &, int &);

    }
}