﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class Timer
    {
        public enum Enum1
        { TIME_LOOP, TIME_MOVE, TIME_COLLIDE,
            TIME_SORT, TIME_COMM, TIME_MODIFY, TIME_OUTPUT, TIME_N };

        public double[] array;

        public void init()
        {
            for (int i = 0; i < (int)Enum1.TIME_N; i++) array[i] = 0.0;
        }
        public void stamp()
        {
            // uncomment if want synchronized timing
            // MPI_Barrier(world);
            previous_time = sparta.mpi.MPI_Wtime();
        }
        public void stamp(int which)
        {
            // uncomment if want synchronized timing
            // MPI_Barrier(world);
            double current_time = sparta.mpi.MPI_Wtime();
            array[which] += current_time - previous_time;
            previous_time = current_time;
        }
        public void barrier_start(int which)
        {
            sparta.mpi.MPI_Barrier(sparta.world);
            array[which] = sparta.mpi.MPI_Wtime();
        }
        public void barrier_stop(int which)
        {
            sparta.mpi.MPI_Barrier(sparta.world);
            double current_time = sparta.mpi.MPI_Wtime();
            array[which] = current_time - array[which];
        }
        public double elapsed(int which)
        {
            double current_time = sparta.mpi.MPI_Wtime();
            return (current_time - array[which]);
        }


        private double previous_time;
        private SPARTA sparta;

        public Timer(SPARTA sparta)
        {
            this.sparta = sparta;
            array = new double[(int)Enum1.TIME_N];
        }
    }
}
