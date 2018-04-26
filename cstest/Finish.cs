using System;
using System.Collections.Generic;
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
            double time, tmp=0, ave, max, min;
            double time_loop, time_other;

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
                if (screen) fprintf(screen,
                        "Loop time of %g on %d procs for %d steps with "
            
                        BIGINT_FORMAT " particles\n",
                        time_loop, nprocs, sparta.update.nsteps, sparta.particle.nglobal);
                if (logfile) fprintf(logfile,
                         "Loop time of %g on %d procs for %d steps with "
            
                         BIGINT_FORMAT " particles\n",
                         time_loop, nprocs, sparta.update.nsteps, sparta.particle.nglobal);
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
                sparta.mpi.MPI_Allreduce(ref sparta.update.nstuck, ref stuck_total, 1, MPI_INT, MPI.MPI_SUM, sparta.world);

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

                char str[32];

                if (me == 0)
                {
                    if (screen)
                    {
                        fprintf(screen, "\n");
                        fprintf(screen, "Particle moves    = " BIGINT_FORMAT " %s\n",
                                nmove_total, MathExtra::num2str(nmove_total, str));
                        fprintf(screen, "Cells touched     = " BIGINT_FORMAT " %s\n",
                                ntouch_total, MathExtra::num2str(ntouch_total, str));
                        fprintf(screen, "Particle comms    = " BIGINT_FORMAT " %s\n",
                                ncomm_total, MathExtra::num2str(ncomm_total, str));
                        fprintf(screen, "Boundary collides = " BIGINT_FORMAT " %s\n",
                                nboundary_total, MathExtra::num2str(nboundary_total, str));
                        fprintf(screen, "Boundary exits    = " BIGINT_FORMAT " %s\n",
                                nexit_total, MathExtra::num2str(nexit_total, str));
                        fprintf(screen, "SurfColl checks   = " BIGINT_FORMAT " %s\n",
                                nscheck_total, MathExtra::num2str(nscheck_total, str));
                        fprintf(screen, "SurfColl occurs   = " BIGINT_FORMAT " %s\n",
                                nscollide_total, MathExtra::num2str(nscollide_total, str));
                        fprintf(screen, "Surf reactions    = " BIGINT_FORMAT " %s\n",
                                nsreact_total, MathExtra::num2str(nsreact_total, str));
                        fprintf(screen, "Collide attempts  = " BIGINT_FORMAT " %s\n",
                                nattempt_total, MathExtra::num2str(nattempt_total, str));
                        fprintf(screen, "Collide occurs    = " BIGINT_FORMAT " %s\n",
                                ncollide_total, MathExtra::num2str(ncollide_total, str));
                        fprintf(screen, "Gas reactions     = " BIGINT_FORMAT " %s\n",
                                nreact_total, MathExtra::num2str(nreact_total, str));
                        fprintf(screen, "Particles stuck   = %d\n", stuck_total);

                        fprintf(screen, "\n");
                        fprintf(screen, "Particle-moves/CPUsec/proc: %g\n", pmsp);
                        fprintf(screen, "Particle-moves/step: %g\n", pms);
                        fprintf(screen, "Cell-touches/particle/step: %g\n", ctps);
                        fprintf(screen, "Particle comm iterations/step: %g\n", cis);
                        fprintf(screen, "Particle fraction communicated: %g\n", pfc);
                        fprintf(screen, "Particle fraction colliding with boundary: %g\n", pfcwb);
                        fprintf(screen, "Particle fraction exiting boundary: %g\n", pfeb);
                        fprintf(screen, "Surface-checks/particle/step: %g\n", schps);
                        fprintf(screen, "Surface-collisions/particle/step: %g\n", sclps);
                        fprintf(screen, "Surface-reactions/particle/step: %g\n", srps);
                        fprintf(screen, "Collision-attempts/particle/step: %g\n", caps);
                        fprintf(screen, "Collisions/particle/step: %g\n", cps);
                        fprintf(screen, "Gas-reactions/particle/step: %g\n", rps);
                    }
                    if (logfile)
                    {
                        fprintf(logfile, "\n");
                        fprintf(logfile, "Particle moves    = " BIGINT_FORMAT " %s\n",
                                nmove_total, MathExtra::num2str(nmove_total, str));
                        fprintf(logfile, "Cells touched     = " BIGINT_FORMAT " %s\n",
                                ntouch_total, MathExtra::num2str(ntouch_total, str));
                        fprintf(logfile, "Particle comms    = " BIGINT_FORMAT " %s\n",
                                ncomm_total, MathExtra::num2str(ncomm_total, str));
                        fprintf(logfile, "Boundary collides = " BIGINT_FORMAT " %s\n",
                                nboundary_total, MathExtra::num2str(nboundary_total, str));
                        fprintf(logfile, "Boundary exits    = " BIGINT_FORMAT " %s\n",
                                nexit_total, MathExtra::num2str(nexit_total, str));
                        fprintf(logfile, "SurfColl checks   = " BIGINT_FORMAT " %s\n",
                                nscheck_total, MathExtra::num2str(nscheck_total, str));
                        fprintf(logfile, "SurfColl occurs   = " BIGINT_FORMAT " %s\n",
                                nscollide_total, MathExtra::num2str(nscollide_total, str));
                        fprintf(logfile, "Surf reactions    = " BIGINT_FORMAT " %s\n",
                                nsreact_total, MathExtra::num2str(nsreact_total, str));
                        fprintf(logfile, "Collide attempts  = " BIGINT_FORMAT " %s\n",
                                nattempt_total, MathExtra::num2str(nattempt_total, str));
                        fprintf(logfile, "Collide occurs    = " BIGINT_FORMAT " %s\n",
                                ncollide_total, MathExtra::num2str(ncollide_total, str));
                        fprintf(logfile, "Reactions         = " BIGINT_FORMAT " %s\n",
                                nreact_total, MathExtra::num2str(nreact_total, str));
                        fprintf(logfile, "Particles stuck   = %d\n", stuck_total);

                        fprintf(logfile, "\n");
                        fprintf(logfile, "Particle-moves/CPUsec/proc: %g\n", pmsp);
                        fprintf(logfile, "Particle-moves/step: %g\n", pms);
                        fprintf(logfile, "Cell-touches/particle/step: %g\n", ctps);
                        fprintf(logfile, "Particle comm iterations/step: %g\n", cis);
                        fprintf(logfile, "Particle fraction communicated: %g\n", pfc);
                        fprintf(logfile, "Particle fraction colliding with boundary: %g\n",
                                pfcwb);
                        fprintf(logfile, "Particle fraction exiting boundary: %g\n", pfeb);
                        fprintf(logfile, "Surface-checks/particle/step: %g\n", schps);
                        fprintf(logfile, "Surface-collisions/particle/step: %g\n", sclps);
                        fprintf(logfile, "Surf-reactions/particle/step: %g\n", srps);
                        fprintf(logfile, "Collision-attempts/particle/step: %g\n", caps);
                        fprintf(logfile, "Collisions/particle/step: %g\n", cps);
                        fprintf(logfile, "Reactions/particle/step: %g\n", rps);
                    }
                }
            }

            // timing breakdowns

            if (timeflag != 0)
            {
                if (me == 0)
                {
                    if (screen) fprintf(screen, "\n");
                    if (logfile) fprintf(logfile, "\n");
                }

                time = sparta.timer.array[(int)EnumTime.Time_MOVE];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    if (screen)
                        fprintf(screen, "Move  time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                    if (logfile)
                        fprintf(logfile, "Move  time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                }

                time = sparta.timer.array[(int)EnumTime.Time_COLLIDE];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    if (screen)
                        fprintf(screen, "Coll  time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                    if (logfile)
                        fprintf(logfile, "Coll  time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                }

                time = sparta.timer.array[(int)EnumTime.Time_SORT];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    if (screen)
                        fprintf(screen, "Sort  time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                    if (logfile)
                        fprintf(logfile, "Sort  time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                }

                time = sparta.timer.array[(int)EnumTime.Time_COMM];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    if (screen)
                        fprintf(screen, "Comm  time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                    if (logfile)
                        fprintf(logfile, "Comm  time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                }

                time = sparta.timer.array[(int)EnumTime.Time_MODIFY];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    if (screen)
                        fprintf(screen, "Modfy time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                    if (logfile)
                        fprintf(logfile, "Modfy time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                }

                time = sparta.timer.array[(int)EnumTime.Time_OUTPUT];
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    if (screen)
                        fprintf(screen, "Outpt time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                    if (logfile)
                        fprintf(logfile, "Outpt time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                }

                time = time_other;
                sparta.mpi.MPI_Allreduce(ref time, ref tmp, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                time = tmp / nprocs;
                if (me == 0)
                {
                    if (screen)
                        fprintf(screen, "Other time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                    if (logfile)
                        fprintf(logfile, "Other time (%%) = %g (%g)\n",
                            time, time / time_loop * 100.0);
                }
            }

            // histograms

            if (histoflag != 0)
            {
                if (me == 0)
                {
                    if (screen) fprintf(screen, "\n");
                    if (logfile) fprintf(logfile, "\n");
                }

                tmp = sparta.particle.nlocal;
                stats(1, &tmp, &ave, &max, &min, 10, histo);
                if (me == 0)
                {
                    if (screen)
                    {
                        fprintf(screen, "Particles: %g ave %g max %g min\n", ave, max, min);
                        fprintf(screen, "Histogram:");
                        for (i = 0; i < 10; i++) fprintf(screen, " %d", histo[i]);
                        fprintf(screen, "\n");
                    }
                    if (logfile)
                    {
                        fprintf(logfile, "Particles: %g ave %g max %g min\n", ave, max, min);
                        fprintf(logfile, "Histogram:");
                        for (i = 0; i < 10; i++) fprintf(logfile, " %d", histo[i]);
                        fprintf(logfile, "\n");
                    }
                }

                tmp = sparta.grid.nlocal;
                stats(1, &tmp, &ave, &max, &min, 10, histo);
                if (me == 0)
                {
                    if (screen)
                    {
                        fprintf(screen, "Cells:     %g ave %g max %g min\n", ave, max, min);
                        fprintf(screen, "Histogram:");
                        for (i = 0; i < 10; i++) fprintf(screen, " %d", histo[i]);
                        fprintf(screen, "\n");
                    }
                    if (logfile)
                    {
                        fprintf(logfile, "Cells:      %g ave %g max %g min\n", ave, max, min);
                        fprintf(logfile, "Histogram:");
                        for (i = 0; i < 10; i++) fprintf(logfile, " %d", histo[i]);
                        fprintf(logfile, "\n");
                    }
                }

                tmp = sparta.grid.nghost;
                stats(1, &tmp, &ave, &max, &min, 10, histo);
                if (me == 0)
                {
                    if (screen)
                    {
                        fprintf(screen, "GhostCell: %g ave %g max %g min\n", ave, max, min);
                        fprintf(screen, "Histogram:");
                        for (i = 0; i < 10; i++) fprintf(screen, " %d", histo[i]);
                        fprintf(screen, "\n");
                    }
                    if (logfile)
                    {
                        fprintf(logfile, "GhostCell: %g ave %g max %g min\n", ave, max, min);
                        fprintf(logfile, "Histogram:");
                        for (i = 0; i < 10; i++) fprintf(logfile, " %d", histo[i]);
                        fprintf(logfile, "\n");
                    }
                }

                tmp = sparta.grid.nempty;
                stats(1, &tmp, &ave, &max, &min, 10, histo);
                if (me == 0)
                {
                    if (screen)
                    {
                        fprintf(screen, "EmptyCell: %g ave %g max %g min\n", ave, max, min);
                        fprintf(screen, "Histogram:");
                        for (i = 0; i < 10; i++) fprintf(screen, " %d", histo[i]);
                        fprintf(screen, "\n");
                    }
                    if (logfile)
                    {
                        fprintf(logfile, "EmptyCell: %g ave %g max %g min\n", ave, max, min);
                        fprintf(logfile, "Histogram:");
                        for (i = 0; i < 10; i++) fprintf(logfile, " %d", histo[i]);
                        fprintf(logfile, "\n");
                    }
                }
            }

            if (logfile) fflush(logfile);
        }


        private void stats(int n, double[] data,
           double[] pave, double[] pmax, double[] pmin,
           int nhisto, int[] histo)
        {
            int i, m;
            int[] histotmp;

            double min = 1.0e20;
            double max = -1.0e20;
            double ave = 0.0;
            for (i = 0; i < n; i++)
            {
                ave += data[i];
                if (data[i] < min) min = data[i];
                if (data[i] > max) max = data[i];
            }

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
            for (i = 0; i < n; i++)
            {
                if (del == 0.0) m = 0;
                else m = Convert.ToInt32((data[i] - min) / del * nhisto);
                if (m > nhisto - 1) m = nhisto - 1;
                histo[m]++;
            }
            histotmp = new int[nhisto];
            //memory->create(histotmp, nhisto, "finish:histotmp");
            //sparta.mpi.MPI_Allreduce(ref histo,ref histotmp, nhisto, MPI.MPI_INT, MPI.MPI_SUM, sparta.world);
            histotmp = histo;
            for (i = 0; i < nhisto; i++) histo[i] = histotmp[i];
            //memory->destroy(histotmp);

            pave = ave;
            pmax = max;
            pmin = min;
        }
    }
}
