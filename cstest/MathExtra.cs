using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using bigint = System.Int64;

namespace cstest
{
    public  static class MathExtra
    {
        /* ----------------------------------------------------------------------
   normalize a vector in place
------------------------------------------------------------------------- */

        public static void norm3(double[] v)
        {
            double scale = 1.0 / Math.Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;
        }

        /* ----------------------------------------------------------------------
           normalize a vector, return in ans
        ------------------------------------------------------------------------- */

        public static void normalize3(double[] v, double[] ans)
        {
            double scale = 1.0 / Math.Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            ans[0] = v[0] * scale;
            ans[1] = v[1] * scale;
            ans[2] = v[2] * scale;
        }

        /* ----------------------------------------------------------------------
           scale a vector to length in place
        ------------------------------------------------------------------------- */

        public static void snorm3(double length, double[] v)
        {
            double scale = length / Math.Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;
        }

        /* ----------------------------------------------------------------------
           scale a vector to length
        ------------------------------------------------------------------------- */

        public static void snormalize3(double length, double[] v, double[] ans)
        {
            double scale = length / Math.Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            ans[0] = v[0] * scale;
            ans[1] = v[1] * scale;
            ans[2] = v[2] * scale;
        }

        /* ----------------------------------------------------------------------
           negate vector v in place
        ------------------------------------------------------------------------- */

        public static void negate3(double[] v)
        {
            v[0] = -v[0];
            v[1] = -v[1];
            v[2] = -v[2];
        }

        /* ----------------------------------------------------------------------
           scale vector v by s in place
        ------------------------------------------------------------------------- */

        public static void scale3(double s, double[] v)
        {
            v[0] *= s;
            v[1] *= s;
            v[2] *= s;
        }

        /* ----------------------------------------------------------------------
           scale vector v by s, return in ans
        ------------------------------------------------------------------------- */

        public static void scale3(double s, double[] v, double[] ans)
        {
            ans[0] = s * v[0];
            ans[1] = s * v[1];
            ans[2] = s * v[2];
        }

        /* ----------------------------------------------------------------------
           axpy: y = alpha*x + y
           y is replaced by result
        ------------------------------------------------------------------------- */

        public static void axpy3(double alpha, double[] x, double[] y)
        {
            y[0] += alpha * x[0];
            y[1] += alpha * x[1];
            y[2] += alpha * x[2];
        }

        /* ----------------------------------------------------------------------
           axpy: ynew = alpha*x + y
        ------------------------------------------------------------------------- */

        public static void axpy3(double alpha, double[] x, double[] y,
                              double[] ynew)
        {
            ynew[0] += alpha * x[0] + y[0];
            ynew[1] += alpha * x[1] + y[1];
            ynew[2] += alpha * x[2] + y[2];
        }

        /* ----------------------------------------------------------------------
           ans = v1 + v2
        ------------------------------------------------------------------------- */

        public static void add3(double[] v1, double[] v2, double[] ans)
        {
            ans[0] = v1[0] + v2[0];
            ans[1] = v1[1] + v2[1];
            ans[2] = v1[2] + v2[2];
        }

        /* ----------------------------------------------------------------------
           ans = v1 - v2
        ------------------------------------------------------------------------- */

        public static void sub3(double[] v1, double[] v2, double[] ans)
        {
            ans[0] = v1[0] - v2[0];
            ans[1] = v1[1] - v2[1];
            ans[2] = v1[2] - v2[2];
        }

        /* ----------------------------------------------------------------------
           length of vector v
        ------------------------------------------------------------------------- */

        public static double len3(double[] v)
        {
            return Math.Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }

        /* ----------------------------------------------------------------------
           squared length of vector v, or dot product of v with itself
        ------------------------------------------------------------------------- */

        public static double lensq3(double[] v)
        {
            return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        }

        /* ----------------------------------------------------------------------
           dot product of 2 vectors
        ------------------------------------------------------------------------- */

        public static double dot3(double[] v1, double[] v2)
        {
            return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        }

        /* ----------------------------------------------------------------------
           cross product of 2 vectors
        ------------------------------------------------------------------------- */

        public static void cross3(double[] v1, double[] v2, double[] ans)
        {
            ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
            ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
            ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
        }

        /* ----------------------------------------------------------------------
           reflect vector v around unit normal n
           return updated v of same length = v - 2(v dot n)n
        ------------------------------------------------------------------------- */

        public static void reflect3(double[] v, double[] n)
        {
            double dot = dot3(v, n);
            v[0] -= 2.0 * dot * n[0];
            v[1] -= 2.0 * dot * n[1];
            v[2] -= 2.0 * dot * n[2];
        }

        /* ----------------------------------------------------------------------
           determinant of a matrix
        ------------------------------------------------------------------------- */

        public static double det3(double[,] m)
        {
            double ans = m[0, 0] * m[1, 1] * m[2, 2] - m[0, 0] * m[1, 2] * m[2, 1] -
              m[1, 0] * m[0, 1] * m[2, 2] + m[1, 0] * m[0, 2] * m[2, 1] +
              m[2, 0] * m[0, 1] * m[1, 2] - m[2, 0] * m[0, 2] * m[1, 1];
            return ans;
        }

        /* ----------------------------------------------------------------------
           diagonal matrix times a full matrix
        ------------------------------------------------------------------------- */

        public static void diag_times3(double[] d, double[,] m, double[,] ans)
        {
            ans[0, 0] = d[0] * m[0, 0];
            ans[0, 1] = d[0] * m[0, 1];
            ans[0, 2] = d[0] * m[0, 2];
            ans[1, 0] = d[1] * m[1, 0];
            ans[1, 1] = d[1] * m[1, 1];
            ans[1, 2] = d[1] * m[1, 2];
            ans[2, 0] = d[2] * m[2, 0];
            ans[2, 1] = d[2] * m[2, 1];
            ans[2, 2] = d[2] * m[2, 2];
        }

        /* ----------------------------------------------------------------------
           add two matrices
        ------------------------------------------------------------------------- */

        public static void plus3(double[,] m, double[,] m2,
                      double[,] ans)
        {
            ans[0, 0] = m[0, 0] + m2[0, 0];
            ans[0, 1] = m[0, 1] + m2[0, 1];
            ans[0, 2] = m[0, 2] + m2[0, 2];
            ans[1, 0] = m[1, 0] + m2[1, 0];
            ans[1, 1] = m[1, 1] + m2[1, 1];
            ans[1, 2] = m[1, 2] + m2[1, 2];
            ans[2, 0] = m[2, 0] + m2[2, 0];
            ans[2, 1] = m[2, 1] + m2[2, 1];
            ans[2, 2] = m[2, 2] + m2[2, 2];
        }

        /* ----------------------------------------------------------------------
           multiply mat1 times mat2
        ------------------------------------------------------------------------- */

        public static void times3(double[,] m, double[,] m2,
                               double[,] ans)
        {
            ans[0, 0] = m[0, 0] * m2[0, 0] + m[0, 1] * m2[1, 0] + m[0, 2] * m2[2, 0];
            ans[0, 1] = m[0, 0] * m2[0, 1] + m[0, 1] * m2[1, 1] + m[0, 2] * m2[2, 1];
            ans[0, 2] = m[0, 0] * m2[0, 2] + m[0, 1] * m2[1, 2] + m[0, 2] * m2[2, 2];
            ans[1, 0] = m[1, 0] * m2[0, 0] + m[1, 1] * m2[1, 0] + m[1, 2] * m2[2, 0];
            ans[1, 1] = m[1, 0] * m2[0, 1] + m[1, 1] * m2[1, 1] + m[1, 2] * m2[2, 1];
            ans[1, 2] = m[1, 0] * m2[0, 2] + m[1, 1] * m2[1, 2] + m[1, 2] * m2[2, 2];
            ans[2, 0] = m[2, 0] * m2[0, 0] + m[2, 1] * m2[1, 0] + m[2, 2] * m2[2, 0];
            ans[2, 1] = m[2, 0] * m2[0, 1] + m[2, 1] * m2[1, 1] + m[2, 2] * m2[2, 1];
            ans[2, 2] = m[2, 0] * m2[0, 2] + m[2, 1] * m2[1, 2] + m[2, 2] * m2[2, 2];
        }

        /* ----------------------------------------------------------------------
           multiply the transpose of mat1 times mat2
        ------------------------------------------------------------------------- */

        public static void transpose_times3(double[,] m, double[,] m2,
                                         double[,] ans)
        {
            ans[0, 0] = m[0, 0] * m2[0, 0] + m[1, 0] * m2[1, 0] + m[2, 0] * m2[2, 0];
            ans[0, 1] = m[0, 0] * m2[0, 1] + m[1, 0] * m2[1, 1] + m[2, 0] * m2[2, 1];
            ans[0, 2] = m[0, 0] * m2[0, 2] + m[1, 0] * m2[1, 2] + m[2, 0] * m2[2, 2];
            ans[1, 0] = m[0, 1] * m2[0, 0] + m[1, 1] * m2[1, 0] + m[2, 1] * m2[2, 0];
            ans[1, 1] = m[0, 1] * m2[0, 1] + m[1, 1] * m2[1, 1] + m[2, 1] * m2[2, 1];
            ans[1, 2] = m[0, 1] * m2[0, 2] + m[1, 1] * m2[1, 2] + m[2, 1] * m2[2, 2];
            ans[2, 0] = m[0, 2] * m2[0, 0] + m[1, 2] * m2[1, 0] + m[2, 2] * m2[2, 0];
            ans[2, 1] = m[0, 2] * m2[0, 1] + m[1, 2] * m2[1, 1] + m[2, 2] * m2[2, 1];
            ans[2, 2] = m[0, 2] * m2[0, 2] + m[1, 2] * m2[1, 2] + m[2, 2] * m2[2, 2];
        }

        /* ----------------------------------------------------------------------
           multiply mat1 times transpose of mat2
        ------------------------------------------------------------------------- */

        public static void times3_transpose(double[,] m, double[,] m2,
                                         double[,] ans)
        {
            ans[0, 0] = m[0, 0] * m2[0, 0] + m[0, 1] * m2[0, 1] + m[0, 2] * m2[0, 2];
            ans[0, 1] = m[0, 0] * m2[1, 0] + m[0, 1] * m2[1, 1] + m[0, 2] * m2[1, 2];
            ans[0, 2] = m[0, 0] * m2[2, 0] + m[0, 1] * m2[2, 1] + m[0, 2] * m2[2, 2];
            ans[1, 0] = m[1, 0] * m2[0, 0] + m[1, 1] * m2[0, 1] + m[1, 2] * m2[0, 2];
            ans[1, 1] = m[1, 0] * m2[1, 0] + m[1, 1] * m2[1, 1] + m[1, 2] * m2[1, 2];
            ans[1, 2] = m[1, 0] * m2[2, 0] + m[1, 1] * m2[2, 1] + m[1, 2] * m2[2, 2];
            ans[2, 0] = m[2, 0] * m2[0, 0] + m[2, 1] * m2[0, 1] + m[2, 2] * m2[0, 2];
            ans[2, 1] = m[2, 0] * m2[1, 0] + m[2, 1] * m2[1, 1] + m[2, 2] * m2[1, 2];
            ans[2, 2] = m[2, 0] * m2[2, 0] + m[2, 1] * m2[2, 1] + m[2, 2] * m2[2, 2];
        }

        /* ----------------------------------------------------------------------
           invert a matrix
           does NOT checks for singular or badly scaled matrix
        ------------------------------------------------------------------------- */

        public static void invert3(double[,] m, double[,] ans)
        {
            double den = m[0, 0] * m[1, 1] * m[2, 2] - m[0, 0] * m[1, 2] * m[2, 1];
            den += -m[1, 0] * m[0, 1] * m[2, 2] + m[1, 0] * m[0, 2] * m[2, 1];
            den += m[2, 0] * m[0, 1] * m[1, 2] - m[2, 0] * m[0, 2] * m[1, 1];

            ans[0, 0] = (m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1]) / den;
            ans[0, 1] = -(m[0, 1] * m[2, 2] - m[0, 2] * m[2, 1]) / den;
            ans[0, 2] = (m[0, 1] * m[1, 2] - m[0, 2] * m[1, 1]) / den;
            ans[1, 0] = -(m[1, 0] * m[2, 2] - m[1, 2] * m[2, 0]) / den;
            ans[1, 1] = (m[0, 0] * m[2, 2] - m[0, 2] * m[2, 0]) / den;
            ans[1, 2] = -(m[0, 0] * m[1, 2] - m[0, 2] * m[1, 0]) / den;
            ans[2, 0] = (m[1, 0] * m[2, 1] - m[1, 1] * m[2, 0]) / den;
            ans[2, 1] = -(m[0, 0] * m[2, 1] - m[0, 1] * m[2, 0]) / den;
            ans[2, 2] = (m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]) / den;
        }

        /* ----------------------------------------------------------------------
           matrix times vector
        ------------------------------------------------------------------------- */

        public static void matvec(double[,] m, double[] v, double[] ans)
        {
            ans[0] = m[0, 0] * v[0] + m[0, 1] * v[1] + m[0, 2] * v[2];
            ans[1] = m[1, 0] * v[0] + m[1, 1] * v[1] + m[1, 2] * v[2];
            ans[2] = m[2, 0] * v[0] + m[2, 1] * v[1] + m[2, 2] * v[2];
        }

        /* ----------------------------------------------------------------------
           matrix times vector
        ------------------------------------------------------------------------- */

        public static void matvec(double[] ex, double[] ey, double[] ez,
                       double[] v, double[] ans)
        {
            ans[0] = ex[0] * v[0] + ey[0] * v[1] + ez[0] * v[2];
            ans[1] = ex[1] * v[0] + ey[1] * v[1] + ez[1] * v[2];
            ans[2] = ex[2] * v[0] + ey[2] * v[1] + ez[2] * v[2];
        }

        /* ----------------------------------------------------------------------
           transposed matrix times vector
        ------------------------------------------------------------------------- */

        public static void transpose_matvec(double[,] m, double[] v,
                         double[] ans)
        {
            ans[0] = m[0, 0] * v[0] + m[1, 0] * v[1] + m[2, 0] * v[2];
            ans[1] = m[0, 1] * v[0] + m[1, 1] * v[1] + m[2, 1] * v[2];
            ans[2] = m[0, 2] * v[0] + m[1, 2] * v[1] + m[2, 2] * v[2];
        }

        /* ----------------------------------------------------------------------
           transposed matrix times vector
        ------------------------------------------------------------------------- */

        public static void transpose_matvec(double[] ex, double[] ey,
                         double[] ez, double[] v,
                         double[] ans)
        {
            ans[0] = ex[0] * v[0] + ex[1] * v[1] + ex[2] * v[2];
            ans[1] = ey[0] * v[0] + ey[1] * v[1] + ey[2] * v[2];
            ans[2] = ez[0] * v[0] + ez[1] * v[1] + ez[2] * v[2];
        }

        /* ----------------------------------------------------------------------
           transposed matrix times diagonal matrix
        ------------------------------------------------------------------------- */

        public static void transpose_diag3(double[,] m, double[] d,
                        double[,] ans)
        {
            ans[0, 0] = m[0, 0] * d[0];
            ans[0, 1] = m[1, 0] * d[1];
            ans[0, 2] = m[2, 0] * d[2];
            ans[1, 0] = m[0, 1] * d[0];
            ans[1, 1] = m[1, 1] * d[1];
            ans[1, 2] = m[2, 1] * d[2];
            ans[2, 0] = m[0, 2] * d[0];
            ans[2, 1] = m[1, 2] * d[1];
            ans[2, 2] = m[2, 2] * d[2];
        }

        /* ----------------------------------------------------------------------
           row vector times matrix
        ------------------------------------------------------------------------- */

        public static void vecmat(double[] v, double[,] m, double[] ans)
        {
            ans[0] = v[0] * m[0, 0] + v[1] * m[1, 0] + v[2] * m[2, 0];
            ans[1] = v[0] * m[0, 1] + v[1] * m[1, 1] + v[2] * m[2, 1];
            ans[2] = v[0] * m[0, 2] + v[1] * m[1, 2] + v[2] * m[2, 2];
        }

        /* ----------------------------------------------------------------------
           matrix times scalar, in place
        ------------------------------------------------------------------------- */

        public static void scalar_times3(double f, double[,] m)
        {
            m[0, 0] *= f; m[0, 1] *= f; m[0, 2] *= f;
            m[1, 0] *= f; m[1, 1] *= f; m[1, 2] *= f;
            m[2, 0] *= f; m[2, 1] *= f; m[2, 2] *= f;
        }

        /* ----------------------------------------------------------------------
           compute quaternion from axis-angle rotation
           v MUST be a unit vector
        ------------------------------------------------------------------------- */

        public static void axisangle_to_quat(double[] v, double angle,
                          double[] quat)
        {
            double halfa = 0.5 * angle;
            double sina = Math.Sin(halfa);
            quat[0] = Math.Cos(halfa);
            quat[1] = v[0] * sina;
            quat[2] = v[1] * sina;
            quat[3] = v[2] * sina;
        }

        /* ----------------------------------------------------------------------
           compute rotation matrix from quaternion
           quat = [w i j k]
        ------------------------------------------------------------------------- */

        public static void quat_to_mat(double[] quat, double[,] mat)
        {
            double w2 = quat[0] * quat[0];
            double i2 = quat[1] * quat[1];
            double j2 = quat[2] * quat[2];
            double k2 = quat[3] * quat[3];
            double twoij = 2.0 * quat[1] * quat[2];
            double twoik = 2.0 * quat[1] * quat[3];
            double twojk = 2.0 * quat[2] * quat[3];
            double twoiw = 2.0 * quat[1] * quat[0];
            double twojw = 2.0 * quat[2] * quat[0];
            double twokw = 2.0 * quat[3] * quat[0];

            mat[0, 0] = w2 + i2 - j2 - k2;
            mat[0, 1] = twoij - twokw;
            mat[0, 2] = twojw + twoik;

            mat[1, 0] = twoij + twokw;
            mat[1, 1] = w2 - i2 + j2 - k2;
            mat[1, 2] = twojk - twoiw;

            mat[2, 0] = twoik - twojw;
            mat[2, 1] = twojk + twoiw;
            mat[2, 2] = w2 - i2 - j2 + k2;
        }
        /* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   1 = lower bound, Nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = 1 to Nmax,
     (3) i* = i to Nmax, (4) *j = 1 to j, (5) i*j = i to j
   return nlo,nhi
   return 0 if successful
   return 1 if numeric values are out of lower/upper bounds
------------------------------------------------------------------------- */

        //public static int bounds(char* str, int nmax, int &nlo, int &nhi)
        //{
        //    char* ptr = strchr(str, '*');

        //    if (ptr == NULL)
        //    {
        //        nlo = nhi = atoi(str);
        //    }
        //    else if (strlen(str) == 1)
        //    {
        //        nlo = 1;
        //        nhi = nmax;
        //    }
        //    else if (ptr == str)
        //    {
        //        nlo = 1;
        //        nhi = atoi(ptr + 1);
        //    }
        //    else if (strlen(ptr + 1) == 0)
        //    {
        //        nlo = atoi(str);
        //        nhi = nmax;
        //    }
        //    else
        //    {
        //        nlo = atoi(str);
        //        nhi = atoi(ptr + 1);
        //    }

        //    if (nlo < 1 || nhi > nmax) return 1;
        //    return 0;
        //}

        /* ----------------------------------------------------------------------
           convert a 64-bit integer into a string like "9.36B"
           K = thousand, M = million, B = billion, T = trillion, P = peta, E = exa
           for easier-to-understand output
        ------------------------------------------------------------------------- */

        public static string num2str(bigint n)
        {
            string outstr;
            if (n < 100000) outstr=string.Format( "({0:G3}K)", 1.0e-3 * n);
            else if (n < 1000000000) outstr=string.Format("({0:G3}M)", 1.0e-6 * n);
            else if (n < 1000000000000) outstr=string.Format("({0:G3}B)", 1.0e-9 * n);
            else if (n < 1000000000000000) outstr=string.Format("({0:G3}T)", 1.0e-12 * n);
            else if (n < 1000000000000000000) outstr=string.Format("({0:G3}P)", 1.0e-15 * n);
            else outstr=string.Format("({0:G3}E)", 1.0e-18 * n);
            return outstr;
        }

    }
}
