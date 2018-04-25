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
    }
}
