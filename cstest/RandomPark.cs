using System;
namespace cstest
{
    public class RanPark
    {
        public const int IA = 16807;
        public const int IM = 2147483647;
        public const double AM = (1.0 / IM);
        public const int IQ = 127773;
        public const int IR = 2836;


        public RanPark(int iseed)
        {
            seed = iseed;
            save = 0;
        }
        public RanPark(double rseed)
        {
            seed = Convert.ToInt32(rseed * IM);
            if (seed == 0) seed = 1;
            save = 0;
        }
        public void reset(double rseed, int offset, int warmup)
        {
            seed = Convert.ToInt32(rseed * IM + offset% IM);
            if (seed < 0) seed = -seed;
            if (seed == 0) seed = 1;
            for (int i = 0; i < warmup; i++) uniform();
        }
        public double uniform()
        {
            int k = seed / IQ;
            seed = IA * (seed - k * IQ) - IR * k;
            if (seed < 0) seed += IM;
            double ans = AM * seed;
            return ans;
        }
        public double gaussian()
        {
            double first, v1, v2, rsq, fac;

            if (save==0)
            {
                while (true)
                {
                    v1 = 2.0 * uniform() - 1.0;
                    v2 = 2.0 * uniform() - 1.0;
                    rsq = v1 * v1 + v2 * v2;
                    if (rsq < 1.0 && rsq != 0.0) break;
                }
                fac = Math.Sqrt(-2.0 * Math.Log(rsq) / rsq);
                second = v1 * fac;
                first = v2 * fac;
                save = 1;
            }
            else
            {
                first = second;
                save = 0;
            }
            return first;
        }


        private int seed, save;
        private double second;
    }
}