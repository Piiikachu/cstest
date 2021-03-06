﻿using System;

namespace cstest
{
    public class SurfCollide
    {
        public string id;
        public string style;
         
        public int vector_flag;          // 0/1 if compute_vector() function exists
        public int size_vector;          // length of global vector

        public SurfCollide(SPARTA sparta, int narg, string[] arg)
        {
            // ID and style
            // ID must be all alphanumeric chars or underscores

            int n = arg[0].Length + 1;
            id = string.Copy(arg[0]);

            for (int i = 0; i < n - 1; i++)
                if (!char.IsLetterOrDigit(id[i]) && id[i] != '_')
                    sparta.error.all("Surf_collide ID must be alphanumeric or underscore characters");

            n = arg[1].Length + 1;
            style = string.Copy(arg[1]);

            vector_flag = 1;
            size_vector = 2;

            nsingle = ntotal = 0;

            copy = 0;
        }
        //public SurfCollide(SPARTA sparta) 
        //{

        //}
        public virtual void init()
        {
            nsingle = ntotal = 0;
        }
        public virtual Particle.OnePart? collide(ref Particle.OnePart? ip, double[] norm,double dtremain, int isr)
        {
            Console.WriteLine("SurfCollide virtual collide");
            return null;
        }

        //virtual void dynamic() { }
        public void tally_update()
        {
            ntotal += nsingle;
            nsingle = 0;
        }
        //      double compute_vector(int i);

        public int copy;

        
        protected int nsingle, ntotal;
        protected double[] one=new double[2], all=new double[2];
    }
}