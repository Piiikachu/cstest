using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    class Program
    {

        static void Main(string[] args)
        {
            string str = Console.ReadLine();
            if (string.IsNullOrWhiteSpace(str))
            {
                args = new string[] { "-in", "in.circle" };
            }
            else if(new FileInfo(str).Exists)
            {
                args = new string[] { "-in", str };
            }
            else
            {
                Console.WriteLine("File {0} doesn't exits",str);
                Console.ReadKey();
                Environment.Exit(1);
            }
            
            MPI mpi = new MPI();
            mpi.MPI_Init(args);

            SPARTA sparta = new SPARTA(args, mpi,MPI.MPI_COMM_WORLD);
            sparta.input.file();
           
            mpi.MPI_Finalize();
            Console.ReadKey();
        }



    }
}
