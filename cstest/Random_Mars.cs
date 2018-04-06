﻿using System;

namespace cstest
{
    public class RanMars
    {
        private int initflag, save;
        private int i97, j97;
        private double c, cd, cm;
        private double second;
        private double[] u;
        private SPARTA sparta;
        public RanMars(SPARTA sparta)
        {
            this.sparta = sparta;
            initflag = 0;
            u = null;
        }
        public void init(int seed)
        {
            int ij, kl, i, j, k, l, ii, jj, m;
            double s, t;

            initflag = 1;

            // assume input seed is positive value > 0
            // insure seed is from 1 to 900,000,000 inclusive

            while (seed > 900000000) seed -= 900000000;

            u = new double[97 + 1];

            ij = (seed - 1) / 30082;
            kl = (seed - 1) - 30082 * ij;
            i = (ij / 177) % 177 + 2;
            j = ij % 177 + 2;
            k = (kl / 169) % 178 + 1;
            l = kl % 169;
            for (ii = 1; ii <= 97; ii++)
            {
                s = 0.0;
                t = 0.5;
                for (jj = 1; jj <= 24; jj++)
                {
                    m = ((i * j) % 179) * k % 179;
                    i = j;
                    j = k;
                    k = m;
                    l = (53 * l + 1) % 169;
                    if ((l * m) % 64 >= 32) s = s + t;
                    t = 0.5 * t;
                }
                u[ii] = s;
            }
            c = 362436.0 / 16777216.0;
            cd = 7654321.0 / 16777216.0;
            cm = 16777213.0 / 16777216.0;
            i97 = 97;
            j97 = 33;
            uniform();
        }
        public double uniform()
        {
            if (initflag==0) sparta.error.all( "Seed command has not been used");

            double uni = u[i97] - u[j97];
            if (uni < 0.0) uni += 1.0;
            u[i97] = uni;
            i97--;
            if (i97 == 0) i97 = 97;
            j97--;
            if (j97 == 0) j97 = 97;
            c -= cd;
            if (c < 0.0) c += cm;
            uni -= c;
            if (uni < 0.0) uni += 1.0;
            return uni;
        }
        public double gaussian()
        {
            double first, v1, v2, rsq, fac;

            if (save==0)
            {
                int again = 1;
                do
                {
                    v1 = 2.0 * uniform() - 1.0;
                    v2 = 2.0 * uniform() - 1.0;
                    rsq = v1 * v1 + v2 * v2;
                    if (rsq < 1.0 && rsq != 0.0) again = 0;
                } while (again!=0);
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


    }
}