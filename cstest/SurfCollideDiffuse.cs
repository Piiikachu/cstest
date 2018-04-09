using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    class SurfCollideDiffuse : SurfCollide
    {
        public SurfCollideDiffuse(SPARTA sparta, int narg, string[] arg) : base(sparta, narg, arg)
        {
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
        //public Particle::OnePart* collide(Particle::OnePart*&, double*, double &, int);

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

  //void diffuse(Particle::OnePart*, double*);
    }
}
