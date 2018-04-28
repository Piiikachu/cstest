using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class Finish
    {
        private SPARTA sparta;
        enum EnumTime
        {
            Time_LOOP, Time_MOVE,  Time_COLLIDE,
            Time_SORT, Time_COMM,  Time_MODIFY,  Time_OUTPUT, Time_N
        };
        public Finish(SPARTA sparta)
        {
            this.sparta = sparta;
        }
        public void end(int flag, double time_multiple_runs)
        {
            int i;
            int[] histo=new int[10];
            int loopflag, statsflag, timeflag, histoflag;
            double time, tmp=0, ave=0, max = 0, min = 0;
            double time_loop=0, time_other=0;

            int me=0, nprocs=0;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);
            sparta.mpi.MPI_Comm_size(sparta.world, ref nprocs);

            // choose flavors of statistical output
            // flag = 0 = just loop summary
            // flag = 1 = dynamics or minimization

            loopflag = 1;
            statsflag = timeflag = histoflag = 0;
            if (flag == 1) statsflag = timeflag = histoflag = 1;

            // loop stats
            // time_multiple_runs used for moves/CPU/proc statistic below

            if (loopflag!=0)
            {
                time_other = sparta.timer.array[(int)EnumTime.Time_LOOP] -
                  (sparta.timer.array[(int)EnumTime.Time_MOVE] + sparta.timer.array[(int)EnumTime.Time_COLLIDE] +
                   sparta.timer.array[(int)EnumTime.Time_SORT] + sparta.timer.array[(int)EnumTime.Time_COMM] +
                   sparta.timer.array[(int)EnumTime.Time_MODIFY] + sparta.timer.array[(int)EnumTime.Time_OUTPUT]);

                time_loop = sparta.timer.array[(int)EnumTime.Time_LOOP];
                sparta.mpi.MPI_Allreduce(ref time_loop, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time_loop = tmp / nprocs;

                if (time_multiple_runs == 0.0) time_multiple_runs = time_loop;
                else
                {
                    tmp = time_multiple_runs;
                    sparta.mpi.MPI_Allreduce(ref tmp, ref time_multiple_runs, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                    time_multiple_runs /= nprocs;
                }
            }

            // recalculate nglobal

            Int64 n = sparta.particle.nlocal;
            sparta.mpi.MPI_Allreduce(ref n, ref sparta.particle.nglobal, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);

            // overall loop time

            if (me == 0)
            {
                string str = string.Format("Loop time of {0} on {1} procs for {2} steps with {3} particles\n", time_loop, nprocs, sparta.update.nsteps, sparta.particle.nglobal);
                Console.WriteLine(str);
                if (sparta.screen!=null) new StreamWriter(sparta.screen).WriteLine(str);
                if (sparta.logfile!=null) new StreamWriter(sparta.logfile).WriteLine(str);
            }

            // cummulative stats over entire run

            if (statsflag!=0)
            {
                Int64 nmove_total=0, ntouch_total=0, ncomm_total=0;
                Int64 nboundary_total=0, nexit_total=0;
                Int64 nscheck_total=0, nscollide_total=0, nsreact_total=0;
                Int64 nattempt_total = 0;
                Int64 ncollide_total = 0;
                Int64 nreact_total = 0;
                int stuck_total=0;

                sparta.mpi.MPI_Allreduce(ref sparta.update.nmove_running, ref nmove_total, 1,
                              MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                sparta.mpi.MPI_Allreduce(ref sparta.update.ntouch_running, ref ntouch_total, 1,
                              MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                sparta.mpi.MPI_Allreduce(ref sparta.update.ncomm_running, ref ncomm_total, 1,
                              MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                sparta.mpi.MPI_Allreduce(ref sparta.update.nboundary_running, ref nboundary_total, 1,
                              MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                sparta.mpi.MPI_Allreduce(ref sparta.update.nexit_running, ref nexit_total, 1,
                              MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                sparta.mpi.MPI_Allreduce(ref sparta.update.nscheck_running, ref nscheck_total, 1,
                              MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                sparta.mpi.MPI_Allreduce(ref sparta.update.nscollide_running, ref nscollide_total, 1,
                              MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                sparta.mpi.MPI_Allreduce(ref sparta.surf.nreact_running, ref nsreact_total, 1,
                              MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                if (sparta.collide!=null)
                {
                    sparta.mpi.MPI_Allreduce(ref sparta.collide.nattempt_running, ref nattempt_total, 1,
                                  MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                    sparta.mpi.MPI_Allreduce(ref sparta.collide.ncollide_running, ref ncollide_total, 1,
                                  MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                    sparta.mpi.MPI_Allreduce(ref sparta.collide.nreact_running, ref nreact_total, 1,
                                  MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
                }
                sparta.mpi.MPI_Allreduce(ref sparta.update.nstuck, ref stuck_total, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);

                double pms, pmsp, ctps, cis, pfc, pfcwb, pfeb, schps, sclps, srps, caps, cps, rps;
                pms = pmsp = ctps = cis = pfc = pfcwb = pfeb =
                  schps = sclps = srps = caps = cps = rps = 0.0;

                Int64 elapsed = sparta.update.ntimestep - sparta.update.first_running_step;
                if (elapsed!=0) pms = 1.0 * nmove_total / elapsed;

                if (nmove_total!=0)
                {
                    pmsp = 1.0 * nmove_total / time_multiple_runs / nprocs;
                    ctps = 1.0 * ntouch_total / nmove_total;
                    cis = 1.0 * sparta.update.niterate_running / elapsed;
                    pfc = 1.0 * ncomm_total / nmove_total;
                    pfcwb = 1.0 * nboundary_total / nmove_total;
                    pfeb = 1.0 * nexit_total / nmove_total;
                    schps = 1.0 * nscheck_total / nmove_total;
                    sclps = 1.0 * nscollide_total / nmove_total;
                    srps = 1.0 * nsreact_total / nmove_total;
                    caps = 1.0 * nattempt_total / nmove_total;
                    cps = 1.0 * ncollide_total / nmove_total;
                    rps = 1.0 * nreact_total / nmove_total;
                }

                string str;

                if (me == 0)
                {
                    str = "\n";
                    str += string.Format("Particle moves    = {0} {1}\n",
                                nmove_total, MathExtra.num2str(nmove_total));
                    str += string.Format("Cells touched     = {0} {1}\n",
                                ntouch_total, MathExtra.num2str(ntouch_total));
                    str += string.Format("Particle comms    = {0} {1}\n",
                                ncomm_total, MathExtra.num2str(ncomm_total));
                    str += string.Format("Boundary collides = {0} {1}\n",
                            nboundary_total, MathExtra.num2str(nboundary_total));
                    str += string.Format("Boundary exits    = {0} {1}\n",
                            nboundary_total, MathExtra.num2str(nexit_total));
                    str += string.Format("SurfColl checks   = {0} {1}\n",
                            nscheck_total, MathExtra.num2str(nscheck_total));
                    str += string.Format("SurfColl occurs   = {0} {1}\n",
                            nscollide_total, MathExtra.num2str(nscollide_total));
                    str += string.Format("Surf reactions    = {0} {1}\n",
                            nsreact_total, MathExtra.num2str(nsreact_total));
                    str += string.Format("Collide attempts  = {0} {1}\n",
                            nattempt_total, MathExtra.num2str(nattempt_total));
                    str += string.Format("Collide occurs    = {0} {1}\n",
                            ncollide_total, MathExtra.num2str(ncollide_total));
                    str += string.Format("Gas reactions     = {0} {1}\n",
                            nreact_total, MathExtra.num2str(nreact_total));
                    str += string.Format("Particles stuck   = {0}\n", stuck_total);

                    str += string.Format("\n");
                    str += string.Format("Particle-moves/CPUsec/proc: {0}\n", pmsp);
                    str += string.Format("Particle-moves/step: {0}\n", pms);
                    str += string.Format("Cell-touches/particle/step: {0}\n", ctps);
                    str += string.Format("Particle comm iterations/step: {0}\n", cis);
                    str += string.Format("Particle fraction communicated: {0}\n", pfc);
                    str += string.Format("Particle fraction colliding with boundary: {0}\n", pfcwb);
                    str += string.Format("Particle fraction exiting boundary: {0}\n", pfeb);
                    str += string.Format("Surface-checks/particle/step: {0}\n", schps);
                    str += string.Format("Surface-collisions/particle/step: {0}\n", sclps);
                    str += string.Format("Surface-reactions/particle/step: {0}\n", srps);
                    str += string.Format("Collision-attempts/particle/step: {0}\n", caps);
                    str += string.Format("Collisions/particle/step: {0}\n", cps);
                    str += string.Format("Gas-reactions/particle/step: {0}\n", rps);
                    Console.WriteLine(str);
                    if (sparta.screen!=null)
                    {
                        new StreamWriter(sparta.screen).WriteLine(str);
                    }
                    if (sparta.logfile!=null)
                    {
                        new StreamWriter(sparta.logfile).WriteLine( str);
                        
                    }
                }
            }

            // timing breakdowns

            if (timeflag != 0)
            {
                if (me == 0)
                {
                    if (sparta.screen!=null) new StreamWriter(sparta.screen).WriteLine( "\n");
                    if (sparta.logfile!=null) new StreamWriter(sparta.logfile).WriteLine( "\n");
                }

                time = sparta.timer.array[(int)EnumTime.Time_MOVE];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    Console.WriteLine("Move  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.screen!=null)
                        new StreamWriter(sparta.screen).WriteLine( "Move  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.logfile!=null)
                        new StreamWriter(sparta.logfile).WriteLine( "Move  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                }

                time = sparta.timer.array[(int)EnumTime.Time_COLLIDE];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    Console.WriteLine("Coll  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.screen!=null)
                        new StreamWriter(sparta.screen).WriteLine( "Coll  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.logfile!=null)
                        new StreamWriter(sparta.logfile).WriteLine( "Coll  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                }

                time = sparta.timer.array[(int)EnumTime.Time_SORT];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    Console.WriteLine("Sort  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.screen!=null)
                        new StreamWriter(sparta.screen).WriteLine( "Sort  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.logfile!=null)
                        new StreamWriter(sparta.logfile).WriteLine( "Sort  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                }

                time = sparta.timer.array[(int)EnumTime.Time_COMM];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    Console.WriteLine("Comm  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.screen!=null)
                        new StreamWriter(sparta.screen).WriteLine( "Comm  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.logfile!=null)
                        new StreamWriter(sparta.logfile).WriteLine( "Comm  time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                }

                time = sparta.timer.array[(int)EnumTime.Time_MODIFY];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    Console.WriteLine("Modfy time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.screen!=null)
                        new StreamWriter(sparta.screen).WriteLine( "Modfy time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.logfile!=null)
                        new StreamWriter(sparta.logfile).WriteLine( "Modfy time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                }

                time = sparta.timer.array[(int)EnumTime.Time_OUTPUT];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    Console.WriteLine("Outpt time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.screen!=null)
                        new StreamWriter(sparta.screen).WriteLine( "Outpt time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.logfile!=null)
                        new StreamWriter(sparta.logfile).WriteLine( "Outpt time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                }

                time = time_other;
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    Console.WriteLine("Other time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.screen!=null)
                        new StreamWriter(sparta.screen).WriteLine( "Other time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                    if (sparta.logfile!=null)
                        new StreamWriter(sparta.logfile).WriteLine( "Other time (%%) = {0} ({1})\n",
                            time, time / time_loop * 100.0);
                }
            }

            // histograms

            if (histoflag != 0)
            {
                if (me == 0)
                {
                    Console.WriteLine("\n");
                    if (sparta.screen!=null) new StreamWriter(sparta.screen).WriteLine( "\n");
                    if (sparta.logfile!=null) new StreamWriter(sparta.logfile).WriteLine( "\n");
                }

                tmp = sparta.particle.nlocal;
                stats(1, ref tmp, ref ave, ref max, ref min, 10, ref histo);
                if (me == 0)
                {
                    string str=string.Format("Particles: {0} ave {1} max {2} min\n", ave, max, min);
                    str += "Histogram:";
                    for (int a = 0; a < 10; a++)
                    {
                        str+=string.Format(" {0}", histo[a]);
                    }
                    str += "\n";
                    Console.WriteLine(str);
                    if (sparta.screen!=null)
                    {
                        new StreamWriter(sparta.screen).WriteLine( str);
                    }
                    if (sparta.logfile!=null)
                    {
                        new StreamWriter(sparta.logfile).WriteLine( str);
                    }
                }

                tmp = sparta.grid.nlocal;
                stats(1, ref tmp, ref ave, ref max, ref min, 10, ref histo);
                if (me == 0)
                {
                    string str=string.Format("Cells:     {0} ave {1} max {2} min\n", ave, max, min);
                    str += "Histogram:";
                    for (int a = 0; a < 10; a++)
                    {
                        str += string.Format(" {0}", histo[a]);
                    }
                    str += "\n";
                    Console.WriteLine(str);
                    if (sparta.screen!=null)
                    {
                        new StreamWriter(sparta.screen).WriteLine( str);
                    }
                    if (sparta.logfile!=null)
                    {
                        new StreamWriter(sparta.logfile).WriteLine(str);
                    }
                }

                tmp = sparta.grid.nghost;
                stats(1, ref tmp, ref ave, ref max, ref min, 10, ref histo);
                if (me == 0)
                {
                    string str = string.Format("GhostCell: {0} ave {1} max {2} min\n", ave, max, min);
                    str += "Histogram:";
                    for (int a = 0; a < 10; a++)
                    {
                        str += string.Format(" {0}", histo[a]);
                    }
                    str += "\n";
                    Console.WriteLine(str);
                    if (sparta.screen!=null)
                    {
                        new StreamWriter(sparta.screen).WriteLine( str);
                    }
                    if (sparta.logfile!=null)
                    {
                        new StreamWriter(sparta.logfile).WriteLine( str);
                    }
                }

                tmp = sparta.grid.nempty;
                stats(1, ref tmp, ref ave, ref max, ref min, 10,ref histo);
                if (me == 0)
                {
                    string str = string.Format("EmptyCell: {0} ave {1} max {2} min\n", ave, max, min);
                    str += "Histogram:";
                    for (int a = 0; a < 10; a++)
                    {
                        str += string.Format(" {0}", histo[a]);
                    }
                    str += "\n";
                    Console.WriteLine(str);
                    if (sparta.screen!=null)
                    {
                        new StreamWriter(sparta.screen).WriteLine(str);
                    }
                    if (sparta.logfile!=null)
                    {
                        new StreamWriter(sparta.logfile).WriteLine( str);
                    }
                }
            }

            if (sparta.logfile!=null) sparta.logfile.Flush();
        }


        private void stats(int n,ref double data,
           ref double pave, ref double pmax,ref double pmin,
           int nhisto,ref int[] histo)
        {
            int i, m;
            int[] histotmp;

            double min = 1.0e20;
            double max = -1.0e20;
            double ave = 0.0;
            //for (i = 0; i < n; i++)
            //{
            //    ave += data[i];
            //    if (data[i] < min) min = data[i];
            //    if (data[i] > max) max = data[i];
            //}
            ave += data;
            int ntotal=0;
            sparta.mpi.MPI_Allreduce(ref n, ref ntotal, 1, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            double tmp=0;
            sparta.mpi.MPI_Allreduce(ref ave, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
            ave = tmp / ntotal;
            sparta.mpi.MPI_Allreduce(ref min, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_MIN, sparta.world);
            min = tmp;
            sparta.mpi.MPI_Allreduce(ref max, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_MAX, sparta.world);
            max = tmp;

            for (i = 0; i < nhisto; i++) histo[i] = 0;


            double del = max - min;
            //for (i = 0; i < n; i++)
            //{
            //    if (del == 0.0) m = 0;
            //    else m = Convert.ToInt32((data - min) / del * nhisto);
            //    if (m > nhisto - 1) m = nhisto - 1;
            //    histo++;
            //}
            if (del == 0.0) m = 0;
            else m = Convert.ToInt32((data - min) / del * nhisto);
            if (m > nhisto - 1) m = nhisto - 1;

            histotmp = new int[nhisto];
            //memory->create(histotmp, nhisto, "finish:histotmp");
            //sparta.mpi.MPI_Allreduce(ref histo,ref histotmp, nhisto, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            histotmp = histo;
            //for (i = 0; i < nhisto; i++) histo[i] = histotmp[i];
            //memory->destroy(histotmp);

            pave = ave;
            pmax = max;
            pmin = min;
        }
    }
}
