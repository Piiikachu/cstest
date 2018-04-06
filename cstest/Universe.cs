using System;
using System.IO;
using MPI_Comm = System.Int32;
namespace cstest
{
    public class Universe
    {
        public string version;          // SPARTA version string = date
         
        public MPI_Comm uworld;        // communicator for entire universe
        public int me, nprocs;          // my place in universe
        
        public FileStream uscreen;          // universe screen output
        public FileStream ulogfile;         // universe logfile
        
        public int existflag;          // 1 if universe exists due to -partition flag
        public int nworlds;            // # of worlds in universe
        public int iworld;             // which world I am in
        public int[] procs_per_world;   // # of procs in each world
        public int[] root_proc;         // root proc in each world
        
        public MPI_Comm uorig;         // original communicator passed to SPARTA instance
        public int[] uni2orig;          // proc I in universe uworld is 
                                // proc uni2orig[I] in original communicator

        public Universe(SPARTA sparta, MPI_Comm communicator)
        {
            version = Version.SPARTA_VERSION;

            uworld = uorig = communicator;
            sparta.mpi.MPI_Comm_rank(uworld,ref me);
            sparta.mpi.MPI_Comm_size(uworld, ref nprocs);

            //uscreen = stdout;
            ulogfile = null;

            existflag = 0;
            nworlds = 0;
            procs_per_world = null;
            root_proc = null;

            //memory->create(uni2orig, nprocs, "universe:uni2orig");
            int[] uni2orig = new int[nprocs];
            for (int i = 0; i < nprocs; i++) uni2orig[i] = i;
        }
        public void Add_world(string str)
        {
            int n, nper;
            string ptr;

            if (str == null)
            {
                n = 1;
                nper = nprocs;
            }
            //else if ((ptr = strchr(str, 'x')) != null)
            //{
            //    ptr = '\0';
            //    n = Int32.Parse(str);
            //    nper = Int32.Parse(ptr + 1);
            //}
            else if (str.Contains("x"))
            {
                string[] split= str.Split(new char[] { 'x'}, 2);
                n =Int32.Parse(split[0]);
                nper = Int32.Parse(split[1]);

            }
            else
            {
                n = 1;
                nper = Int32.Parse(str);
            }

            //memory->grow(procs_per_world, nworlds + n, "universe:procs_per_world");
            //memory->grow(root_proc, (nworlds + n), "universe:root_proc");
            procs_per_world = new int[nworlds + n];
            root_proc = new int[nworlds + n];

            for (int i = 0; i < n; i++)
            {
                procs_per_world[nworlds] = nper;
                if (nworlds == 0) root_proc[nworlds] = 0;
                else
                    root_proc[nworlds] = root_proc[nworlds - 1] + procs_per_world[nworlds - 1];
                if (me >= root_proc[nworlds]) iworld = nworlds;
                nworlds++;
            }
        }
        public int Consistent()
        {
            int n = 0;
            for (int i = 0; i < nworlds; i++) n += procs_per_world[i];
            if (n == nprocs) return 1;
            else return 0;
        }
    }
}