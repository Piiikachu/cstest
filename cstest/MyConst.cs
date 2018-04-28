using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    class MyConst
    {
        public const double THIRD = 1.0 / 3.0;
        public const double MY_PI = 3.14159265358979323846; // pi
        public const double MY_2PI = 6.28318530717958647692; // 2pi
        public const double MY_3PI = 9.42477796076937971538; // 3pi
        public const double MY_PI2 = 1.57079632679489661923; // pi/2
        public const double MY_PI3 = 1.04719755119659774615; // pi/3
        public const double MY_PI4 = 0.78539816339744830962; // pi/4
        public const double MY_PIS = 1.77245385090551602729; // sqrt(pi)
        public static double Erf(double x)
        {
            // constants
            double a1 = 0.254829592;
            double a2 = -0.284496736;
            double a3 = 1.421413741;
            double a4 = -1.453152027;
            double a5 = 1.061405429;
            double p = 0.3275911;

            // Save the sign of x
            int sign = 1;
            if (x < 0)
                sign = -1;
            x = Math.Abs(x);

            // A&S formula 7.1.26
            double t = 1.0 / (1.0 + p * x);
            double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);

            return sign * y;
        }

        const double LOG_2PI = 1.83787706640934548;  /* $\log 2\pi$ */
        const int N = 8;

        const double B0 = 1;                 /* Bernoulli numbers */
        const double B1 = (-1.0 / 2.0);
        const double B2 = (1.0 / 6.0);
        const double B4 = (-1.0 / 30.0);
        const double B6 = (1.0 / 42.0);
        const double B8 = (-1.0 / 30.0);
        const double B10 = (5.0 / 66.0);
        const double B12 = (-691.0 / 2730.0);
        const double B14 = (7.0 / 6.0);
        const double B16 = (-3617.0 / 510.0);
                    
        static double loggamma(double x)  /* the natural logarithm of the Gamma function. */
        {
            double v, w;

            v = 1;
            while (x < N) { v *= x; x++; }
            w = 1 / (x * x);
            return ((((((((B16 / (16 * 15)) * w + (B14 / (14 * 13))) * w
                        + (B12 / (12 * 11))) * w + (B10 / (10 * 9))) * w
                        + (B8 / (8 * 7))) * w + (B6 / (6 * 5))) * w
                        + (B4 / (4 * 3))) * w + (B2 / (2 * 1))) / x
                        + 0.5 * LOG_2PI - Math.Log(v) - x + (x - 0.5) * Math.Log(x);
        }

        public static double Gamma(double x)  /* Gamma function */
        {
            if (x == 0.0)
            { /* Pole Error */
                //errno = ERANGE;
                
                return 1 / x < 0 ? double.NegativeInfinity : double.PositiveInfinity;
                throw new Exception();
            }
            if (x < 0)
            {
                int sign;
                double zero = 0.0;
                double i, f;
                i = Math.Floor(-x);
                f = -x - i;
                //f = modf(-x, &i);
                if (f == 0.0)
                { /* Domain Error */
                    //errno = EDOM;
                    throw new Exception() ;
                    return zero / zero;
                }
                sign = (i% 2.0 != 0.0) ? 1 : -1;
                return sign * MY_PI / (Math.Sin(MY_PI * f) * Math.Exp(loggamma(1 - x)));
            }
            return Math.Exp(loggamma(x));
        }

    }
}
