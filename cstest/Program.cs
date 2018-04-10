using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    //todo:update.cpp line 277
    class Program
    {

        static void Main(string[] args)
        {
            args =new string[]{ "-in","in.circle"};
            MPI mpi = new MPI();
            mpi.MPI_Init(args);

            SPARTA sparta = new SPARTA(args, mpi,MPI.MPI_COMM_WORLD);
            sparta.input.file();

            mpi.MPI_Finalize();
            Console.ReadKey();
        }



    }
}
