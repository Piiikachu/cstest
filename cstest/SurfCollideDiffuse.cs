using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    class SurfCollideDiffuse : SurfCollide
    {
        private SPARTA sparta;
        public SurfCollideDiffuse(SPARTA sparta, int narg, string[] arg) : base(sparta, narg, arg)
        {
            this.sparta = sparta;
            if (narg < 4) sparta.error.all("Illegal surf_collide diffuse command");

            tstr = null;

            if (arg[2].Contains("v_"))
            {
                //int n = arg[2,2].Length + 1;
                //tstr = new char[n];
                //strcpy(tstr, &arg[2][2]);
                tstr = string.Copy(arg[2]);

            }
            else
            {
                twall = double.Parse(arg[2]);
                if (twall <= 0.0) sparta.error.all("Surf_collide diffuse temp <= 0.0");
            }

            acc = double.Parse (arg[3]);
            if (acc < 0.0 || acc > 1.0)
                sparta.error.all("Illegal surf_collide diffuse command");

            // optional args
            tflag = rflag = 0;

            int iarg = 4;

            while (iarg < narg)
            {
                switch (arg[iarg])
                {
                    case "translate":
                        if (iarg + 4 > narg)
                            sparta.error.all("Illegal surf_collide diffuse command");
                        tflag = 1;
                        vx = double.Parse(arg[iarg + 1]);
                        vy = double.Parse(arg[iarg + 2]);
                        vz = double.Parse(arg[iarg + 3]);
                        iarg += 4;
                        break;
                    case "rotate":
                        if (iarg + 7 > narg)
                            sparta.error.all("Illegal surf_collide diffuse command");
                        rflag = 1;
                        px = double.Parse(arg[iarg + 1]);
                        py = double.Parse(arg[iarg + 2]);
                        pz = double.Parse(arg[iarg + 3]);
                        wx = double.Parse(arg[iarg + 4]);
                        wy = double.Parse(arg[iarg + 5]);
                        wz = double.Parse(arg[iarg + 6]);
                        if (sparta.domain.dimension == 2 && pz != 0.0)
                            sparta.error.all("Surf_collide diffuse rotation invalid for 2d");
                        if (sparta.domain.dimension == 2 && (wx != 0.0 || wy != 0.0))
                            sparta.error.all("Surf_collide diffuse rotation invalid for 2d");
                        iarg += 7;
                        break;
                    default:
                        sparta.error.all("Illegal surf_collide diffuse command");
                        break;
                }




            }
            if (tflag!=0 && rflag!=0) sparta.error.all("Illegal surf_collide diffuse command");
            if (tflag!=0 || rflag!=0) trflag = 1;
            else trflag = 0;

            // initialize RNG

            random = new RanPark(sparta.update.ranmaster.uniform());
            double seed = sparta.update.ranmaster.uniform();
            random.reset(seed, sparta.comm.me, 100);

        }
        //public void init();
        public override Particle.OnePart? collide(ref Particle.OnePart? ip, double[] norm,double dtremain, int isr)
        {
            nsingle++;

            // if surface chemistry defined, attempt reaction
            // reaction = 1 if reaction took place

            Particle.OnePart iorig;
            Particle.OnePart? jp=null ;
            int reaction = 0;

            if (isr >= 0)
            {
                if (sparta.modify.n_surf_react!=0)
                {
                    //memcpy(&iorig, ip, sizeof(Particle.OnePart));
                    iorig = ip.Value;
                }

                reaction = sparta.surf.sr[isr].react(ref ip, norm, out jp);
                if (reaction!=0)
                {
                    sparta.surf.nreact_one++;
                }
            }
           

            // diffuse reflection for each particle
            if (ip!=null)
            {
                diffuse(ip.Value, norm);
            }
            if (jp!=null)
            {
                diffuse(jp.Value, norm);
            }

            // call any fixes with a surf_react() method
            // they may reset j to -1, e.g. fix ambipolar
            //   in which case newly created j is deleted

            if (reaction!=0 && sparta.modify.n_surf_react!=0)
            {
                //int i = -1;
                //if (ip) i = ip - sparta.particle.particles[0];
                //int j = -1;
                //if (jp) j = jp - particle.particles;
                //sparta.modify.surf_react(&iorig, i, j);
                //if (jp && j < 0)
                //{
                //jp = NULL;
                //particle.nlocal--;
                //}
                Console.WriteLine("SurfCollideDiffuse.collide-> react");
            }

            return jp;
        }

        //public void dynamic();
        protected double twall;              // surface temperature
        protected double acc;                // surface accomodation coeff
        protected double vx, vy, vz;           // translational velocity of surface
        protected double wx, wy, wz;           // angular velocity of surface
        protected double px, py, pz;           // point to rotate surface around
        protected int tflag, rflag;           // flags for translation and rotation
        protected int trflag;                // 1 if either tflag or rflag is set
         
        protected string tstr;                // temperature variable name (NULL if constant)
        protected int tvar;                  // index of equal-style variable
         
        protected  RanPark random;     // RNG for particle reflection

        void diffuse(Particle.OnePart p, double[] norm)
        {
            // specular reflection
            // reflect incident v around norm

            if (random.uniform() > acc)
            {
                MathExtra.reflect3(p.v, norm);

                // diffuse reflection
                // vrm = most probable speed of species, eqns (4.1) and (4.7)
                // vperp = velocity component perpendicular to surface along norm, eqn (12.3)
                // vtan12 = 2 velocity components tangential to surface
                // tangent1 = component of particle v tangential to surface,
                //   check if tangent1 = 0 (normal collision), set randomly
                // tangent2 = norm x tangent1 = orthogonal tangential direction
                // tangent12 are both unit vectors

            }
            else
            {
                double[] tangent1=new double[3], tangent2=new double[3];
                List<Particle.Species> species = sparta.particle.species;
                int ispecies = p.ispecies;

                double vrm = Math.Sqrt(2.0 * sparta.update.boltz * twall / species[ispecies].mass);
                double vperp = vrm * Math.Sqrt(-Math.Log(random.uniform()));

                double theta = MyConst.MY_2PI * random.uniform();
                double vtangent = vrm * Math.Sqrt(-Math.Log(random.uniform()));
                double vtan1 = vtangent * Math.Sin(theta);
                double vtan2 = vtangent * Math.Cos(theta);

                double[] v = p.v;
                double dot = MathExtra.dot3(v, norm);

                double beta_un, normalized_distbn_fn;

                tangent1[0] = v[0] - dot * norm[0];
                tangent1[1] = v[1] - dot * norm[1];
                tangent1[2] = v[2] - dot * norm[2];

                if (MathExtra.lensq3(tangent1) == 0.0)
                {
                    tangent2[0] = random.uniform();
                    tangent2[1] = random.uniform();
                    tangent2[2] = random.uniform();
                    MathExtra.cross3(norm, tangent2, tangent1);
                }

                MathExtra.norm3(tangent1);
                MathExtra.cross3(norm, tangent1, tangent2);

                // add in translation or rotation vector if specified
                // only keep portion of vector tangential to surface element

                if (trflag!=0)
                {
                    double vxdelta, vydelta, vzdelta;
                    if (tflag != 0)
                    {
                        vxdelta = vx; vydelta = vy; vzdelta = vz;
                        double adot = vxdelta * norm[0] + vydelta * norm[1] + vzdelta * norm[2];

                        if (Math.Abs(adot) > 0.001)
                        {
                            adot /= vrm;
                            do
                            {
                                do
                                {
                                    beta_un = (6.0 * random.uniform() - 3.0);
                                } while (beta_un + adot < 0.0);
                                normalized_distbn_fn = 2.0 * (beta_un + adot) /
                                  (adot + Math.Sqrt(adot * adot + 2.0)) *
                                  Math.Exp(0.5 + (0.5 * adot) * (adot - Math.Sqrt(adot * adot + 2.0)) -
                                      beta_un * beta_un);
                            } while (normalized_distbn_fn < random.uniform());
                            vperp = beta_un * vrm;
                        }

                    }
                    else
                    {
                        double[] x = p.x;
                        vxdelta = wy * (x[2] - pz) - wz * (x[1] - py);
                        vydelta = wz * (x[0] - px) - wx * (x[2] - pz);
                        vzdelta = wx * (x[1] - py) - wy * (x[0] - px);
                        double adot = vxdelta * norm[0] + vydelta * norm[1] + vzdelta * norm[2];
                        vxdelta -= adot * norm[0];
                        vydelta -= adot * norm[1];
                        vzdelta -= adot * norm[2];
                    }

                    v[0] = vperp * norm[0] + vtan1 * tangent1[0] + vtan2 * tangent2[0] + vxdelta;
                    v[1] = vperp * norm[1] + vtan1 * tangent1[1] + vtan2 * tangent2[1] + vydelta;
                    v[2] = vperp * norm[2] + vtan1 * tangent1[2] + vtan2 * tangent2[2] + vzdelta;

                    // no translation or rotation

                }
                else
                {
                    v[0] = vperp * norm[0] + vtan1 * tangent1[0] + vtan2 * tangent2[0];
                    v[1] = vperp * norm[1] + vtan1 * tangent1[1] + vtan2 * tangent2[1];
                    v[2] = vperp * norm[2] + vtan1 * tangent1[2] + vtan2 * tangent2[2];
                }

                p.erot = sparta.particle.erot(ispecies, twall, random);
                p.evib = sparta.particle.evib(ispecies, twall, random);
            }
        }
    }
}
