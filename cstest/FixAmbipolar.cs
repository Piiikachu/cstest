namespace cstest
{
    public class FixAmbipolar : Fix
    {
        public int especies;               // index of electron species
        public int[] ions;                  // 1 if a particle species is an ionx

        //public int setmask();
        //public void init();
        //public void add_particle(int, double, double, double, double*);
        //public void surf_react(Particle::OnePart*, int &, int &);


        private int maxion;                 // length of ions vector
        private int ionindex, velindex;      // indices into particle custom data structs
        private RanPark random;
        public FixAmbipolar(SPARTA sparta, int narg, string[] arg) : base(sparta, narg, arg)
        {
        }
    }
}