using System.IO;

//todo: 需要优化  委托

namespace cstest
{
    public class Error
    {
        private SPARTA sparta;
        public Error(SPARTA sparta)
        {
            this.sparta = sparta;
        }

        public void universe_all(string str)
        {
            sparta.mpi.MPI_Barrier(sparta.universe.uworld);
            if (sparta.universe.me == 0)
            {
                string strstr = string.Format("ERROR: {0} \n", str);
                System.Console.WriteLine(strstr);
                if (sparta.universe.uscreen != null)
                {
                   
                    new StreamWriter(sparta.universe.uscreen).WriteLine(strstr);
                 
                }

                if (sparta.universe.ulogfile != null)
                {
                    new StreamWriter(sparta.universe.ulogfile).WriteLine(strstr);

                }
            }

            //if (output!=null) delete output;
            if (sparta.universe.nworlds > 1)
            {
                if (sparta.screen != null) sparta.screen.Close();
                if (sparta.logfile != null) sparta.logfile.Close();
            }
            if (sparta.universe.ulogfile != null) sparta.universe.ulogfile.Close();

            sparta.mpi.MPI_Finalize();
            System.Console.WriteLine(" exit(1);");
            System.Console.ReadKey();
        }

        public void universe_all(string file, int line, string str)
        {
            sparta.mpi.MPI_Barrier(sparta.universe.uworld);

            if (sparta.universe.me == 0)
            {
                string strstr = string.Format("ERROR: {0} ({1}:{2})\n", str, file, line);
                System.Console.WriteLine(strstr);
                if (sparta.universe.uscreen != null)
                {
                    System.Console.WriteLine(strstr);
                    //fprintf(sparta.universe.uscreen,
                    //           "ERROR: %s (%s:%d)\n", str, file, line);
                }

                if (sparta.universe.ulogfile != null)
                {
                    StreamWriter sw = new StreamWriter(sparta.universe.ulogfile);
                    sw.WriteLine(strstr);

                    //fprintf(sparta.universe.ulogfile,
                    //            "ERROR: %s (%s:%d)\n", str, file, line);
                }
            }

            //if (output!=null) delete output;
            if (sparta.universe.nworlds > 1)
            {
                if (sparta.screen != null) sparta.screen.Close();
                if (sparta.logfile != null) sparta.logfile.Close();
            }
            if (sparta.universe.ulogfile != null) sparta.universe.ulogfile.Close();

            sparta.mpi.MPI_Finalize();
            System.Console.WriteLine(" exit(1);");
            System.Console.ReadKey();
        }


        public void universe_one(string str)
        {
            if (sparta.universe.uscreen != null)
            {
                string tempstr = string.Format("ERROR on proc {0}: {1} \n", sparta.universe.me, str);
                System.Console.WriteLine(tempstr);
                new StreamWriter(sparta.universe.uscreen).Write(tempstr);
            }

            sparta.mpi.MPI_Abort(sparta.universe.uworld, 1);
        }

        public void universe_one(string file, int line, string str)
        {
            if (sparta.universe.uscreen != null)
            {
                string tempstr = string.Format("ERROR on proc {0}: {1} ({2}:{3})\n", sparta.universe.me, str, file, line);
                System.Console.WriteLine(tempstr);
                //fprintf(universe->uscreen, "ERROR on proc %d: %s (%s:%d)\n",                    universe->me, str, file, line);
                StreamWriter sw = new StreamWriter(sparta.universe.uscreen);
                sw.Write(tempstr);


            }

            sparta.mpi.MPI_Abort(sparta.universe.uworld, 1);
        }

        public void all(string str)
        {
            sparta.mpi.MPI_Barrier(sparta.world);

            int me = 0;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);

            if (me == 0)
            {
                string strstr = string.Format("ERROR: {0} \n", str);
                System.Console.WriteLine(strstr);
                if (sparta.screen != null)
                {
                    new StreamWriter(sparta.screen).Write(strstr);
                }

                if (sparta.logfile != null)
                {
                    new StreamWriter(sparta.logfile).Write(strstr);
                }
            }

            //if (output) delete output;
            if (sparta.screen != null) sparta.screen.Close();
            if (sparta.logfile != null) sparta.logfile.Close();

            sparta.mpi.MPI_Finalize();
            System.Console.WriteLine("exit(1)");
            System.Console.ReadKey();
        }

        public void all(string file, int line, string str)
        {
            sparta.mpi.MPI_Barrier(sparta.world);

            int me = 0;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);

            if (me == 0)
            {
                string strstr = string.Format("ERROR: {0} ({1}:{2})\n", str, file, line);
                System.Console.WriteLine(strstr);
                if (sparta.screen != null)
                {
                    new StreamWriter(sparta.screen).Write(strstr);
                }

                if (sparta.logfile != null)
                {
                    new StreamWriter(sparta.logfile).Write(strstr);
                }
            }

            //if (output) delete output;
            if (sparta.screen != null) sparta.screen.Close();
            if (sparta.logfile != null) sparta.logfile.Close();

            sparta.mpi.MPI_Finalize();
            System.Console.WriteLine("exit(1)");
            System.Console.ReadKey();
        }


        public void one(string str)
        {
            int me = 0;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);
            string strstr = string.Format("ERROR on proc {0}: {1} \n", me, str);
            System.Console.WriteLine(strstr);
            if (sparta.screen != null)
            {

                new StreamWriter(sparta.screen).Write(strstr);
            }

            if (sparta.universe.nworlds > 1)
            {
                new StreamWriter(sparta.universe.uscreen).Write(strstr);
            }

            sparta.mpi.MPI_Abort(sparta.world, 1);
            System.Console.ReadKey();
        }


        public void one(string file, int line, string str)
        {
            int me = 0;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);
            string strstr = string.Format("ERROR on proc {0}: {1} ({2}:{3})\n", me, str, file, line);
            System.Console.WriteLine(strstr);
            if (sparta.screen != null)
            {

                new StreamWriter(sparta.screen).Write(strstr);
            }

            if (sparta.universe.nworlds > 1)
            {
                new StreamWriter(sparta.universe.uscreen).Write(strstr);
            }

            sparta.mpi.MPI_Abort(sparta.world, 1);
        }
        public void warning(string file, int line, string str, int logflag)
        {
            string strstr = string.Format("WARNING: {0} ({1}:{2})\n", str, file, line);
            System.Console.WriteLine(strstr);
            if (sparta.screen!=null)
            {
                new StreamWriter(sparta.screen).Write(strstr);
            }

            if (logflag!=0 && sparta.logfile!=null)
            {
                new StreamWriter(sparta.logfile).Write(strstr);
            }
        }
        public void message(string file, int line, string str, int logflag)
        {
            string strstr = string.Format("{0} ({1}:{2})\n", str, file, line);
            System.Console.WriteLine(strstr);
            if (sparta.screen != null)
            {
                new StreamWriter(sparta.screen).Write(strstr);
            }

            if (logflag != 0 && sparta.logfile != null)
            {
                new StreamWriter(sparta.logfile).Write(strstr);
            }
            
        }
        public void done()
        {
            sparta.mpi.MPI_Barrier(sparta.world);

            //if (output) delete output;
            if (sparta.screen!=null) sparta.screen.Close();
            if (sparta.logfile!=null) sparta.logfile.Close();

            sparta.mpi.MPI_Finalize();
            System.Console.WriteLine("exit(1);");
            System.Console.ReadKey();
        }
    }
}