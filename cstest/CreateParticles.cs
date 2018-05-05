using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class CreateParticles
    {
        private SPARTA sparta;
        enum Enum1{ UNKNOWN, OUTSIDE, INSIDE, OVERLAP };   // same as Grid

        public const double EPSZERO = 1.0e-14;
        public CreateParticles(SPARTA sparta)
        {
            this.sparta = sparta;
        }
        public void command(int narg, string[] arg)
        {
            if (sparta.grid.exist==0)
                sparta.error.all("Cannot create particles before grid is defined");

            sparta.particle.exist = 1;

            if (narg < 1) sparta.error.all("Illegal create_particles command");

            imix = sparta.particle.find_mixture(arg[0]);
            if (imix < 0) sparta.error.all("Create_particles mixture ID does not exist");
            sparta.particle.mixture[imix].init();

            // style arg

            Int64 np = 0;
            single = 0;

            int iarg = 1;
            if (string.Compare(arg[iarg], "n") == 0)
            {
                if (iarg + 2 > narg) sparta.error.all("Illegal create_particles command");
                np = Int64.Parse(arg[iarg + 1]);
                if (np < 0) sparta.error.all("Illegal create_particles command");
                iarg += 2;
            }
            else if (string.Compare(arg[iarg], "single") == 0)
            {
                if (iarg + 8 > narg) sparta.error.all("Illegal create_particles command");
                single = 1;
                mspecies = sparta.particle.find_species(arg[iarg + 1]);
                if (mspecies < 0)
                    sparta.error.all("Create_particles species ID does not exist");
                xp = double.Parse(arg[iarg + 2]);
                yp = double.Parse(arg[iarg + 3]);
                zp = double.Parse(arg[iarg + 4]);
                vx = double.Parse(arg[iarg + 5]);
                vy = double.Parse(arg[iarg + 6]);
                vz = double.Parse(arg[iarg + 7]);
                iarg += 8;
            }
            else sparta.error.all("Illegal create_particles command");

            // optional args

            int globalflag = 0;
            twopass = 0;
            region = null;
            speciesflag = densflag = velflag = tempflag = 0;
            sstr = sxstr = systr = szstr = null;
            dstr = dxstr = dystr = dzstr = null;
            tstr = txstr = tystr = tzstr = null;
            vxstr = vystr = vzstr = vstrx = vstry = vstrz = null;

            while (iarg < narg)
            {
                if (string.Compare(arg[iarg], "global") == 0)
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal create_particles command");
                    if (string.Compare(arg[iarg + 1], "no") == 0) globalflag = 0;
                    else if (string.Compare(arg[iarg + 1], "yes") == 0) globalflag = 1;
                    else sparta.error.all("Illegal create_particles command");
                    iarg += 2;
                }
                else if (string.Compare(arg[iarg], "region") == 0)
                {
                    if (iarg + 2 > narg) sparta.error.all("Illegal create_particles command");
                    int iregion = sparta.domain.find_region(arg[iarg + 1]);
                    if (iregion < 0)
                        sparta.error.all("Create_particles region does not exist");
                    region = sparta.domain.regions[iregion];
                    iarg += 2;
                }
                else if (string.Compare(arg[iarg], "species") == 0)
                {
                    if (iarg + 5 > narg) sparta.error.all("Illegal create_particles command");
                    speciesflag = 1;
                    sstr = arg[iarg + 1];
                    if (string.Compare(arg[iarg + 2], "null") == 0) sxstr = null;
                    else sxstr = arg[iarg + 2];
                    if (string.Compare(arg[iarg + 3], "null") == 0) systr = null;
                    else systr = arg[iarg + 3];
                    if (string.Compare(arg[iarg + 4], "null") == 0) szstr = null;
                    else szstr = arg[iarg + 4];
                    iarg += 5;
                }
                else if (string.Compare(arg[iarg], "density") == 0)
                {
                    if (iarg + 5 > narg) sparta.error.all("Illegal create_particles command");
                    densflag = 1;
                    dstr = arg[iarg + 1];
                    if (string.Compare(arg[iarg + 2], "null") == 0) dxstr = null;
                    else dxstr = arg[iarg + 2];
                    if (string.Compare(arg[iarg + 3], "null") == 0) dystr = null;
                    else dystr = arg[iarg + 3];
                    if (string.Compare(arg[iarg + 4], "null") == 0) dzstr = null;
                    else dzstr = arg[iarg + 4];
                    iarg += 5;
                }
                else if (string.Compare(arg[iarg], "temperature") == 0)
                {
                    if (iarg + 5 > narg) sparta.error.all("Illegal create_particles command");
                    tempflag = 1;
                    tstr = arg[iarg + 1];
                    if (string.Compare(arg[iarg + 2], "null") == 0) txstr = null;
                    else txstr = arg[iarg + 2];
                    if (string.Compare(arg[iarg + 3], "null") == 0) tystr = null;
                    else tystr = arg[iarg + 3];
                    if (string.Compare(arg[iarg + 4], "null") == 0) tzstr = null;
                    else tzstr = arg[iarg + 4];
                    iarg += 5;
                }
                else if (string.Compare(arg[iarg], "velocity") == 0)
                {
                    if (iarg + 7 > narg) sparta.error.all("Illegal create_particles command");
                    velflag = 1;
                    if (string.Compare(arg[iarg + 1], "null") == 0) vxstr = null;
                    else vxstr = arg[iarg + 1];
                    if (string.Compare(arg[iarg + 2], "null") == 0) vystr = null;
                    else vystr = arg[iarg + 2];
                    if (string.Compare(arg[iarg + 3], "null") == 0) vzstr = null;
                    else vzstr = arg[iarg + 3];
                    if (string.Compare(arg[iarg + 4], "null") == 0) vstrx = null;
                    else vstrx = arg[iarg + 4];
                    if (string.Compare(arg[iarg + 5], "null") == 0) vstry = null;
                    else vstry = arg[iarg + 5];
                    if (string.Compare(arg[iarg + 6], "null") == 0) vstrz = null;
                    else vstrz = arg[iarg + 6];
                    iarg += 7;
                }
                else if (string.Compare(arg[iarg], "twopass") == 0)
                {
                    if (iarg + 1 > narg) sparta.error.all("Illegal create_particles command");
                    twopass = 1;
                    iarg += 1;
                }
                else sparta.error.all("Illegal create_particles command");
            }

            if (globalflag!=0)
                sparta.error.all("Create_particles global option not yet implemented");

            // error checks and further setup for variables

            if (speciesflag != 0)
            {
                svar = sparta.input.variable.find(sstr);
                if (svar < 0)
                    sparta.error.all("Variable name for create_particles does not exist");
                if (sparta.input.variable.equal_style(svar)==0)
                    sparta.error.all("Variable for create_particles is invalid style");
                if (sxstr != null)
                {
                    sxvar = sparta.input.variable.find(sxstr);
                    if (sxvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(sxvar)==0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (systr != null)
                {
                    syvar = sparta.input.variable.find(systr);
                    if (syvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(syvar)==0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (szstr != null)
                {
                    szvar = sparta.input.variable.find(szstr);
                    if (szvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(szvar)==0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
            }

            if (densflag != 0)
            {
                dvar = sparta.input.variable.find(dstr);
                if (dvar < 0)
                    sparta.error.all("Variable name for create_particles does not exist");
                if (sparta.input.variable.equal_style(dvar)==0)
                    sparta.error.all("Variable for create_particles is invalid style");
                if (dxstr != null)
                {
                    dxvar = sparta.input.variable.find(dxstr);
                    if (dxvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(dxvar)==0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (dystr != null)
                {
                    dyvar = sparta.input.variable.find(dystr);
                    if (dyvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(dyvar)==0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (dzstr != null)
                {
                    dzvar = sparta.input.variable.find(dzstr);
                    if (dzvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(dzvar) == 0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
            }

            if (tempflag != 0)
            {
                tvar = sparta.input.variable.find(tstr);
                if (tvar < 0)
                    sparta.error.all("Variable name for create_particles does not exist");
                if (sparta.input.variable.equal_style(tvar)==0)
                    sparta.error.all("Variable for create_particles is invalid style");
                if (txstr != null)
                {
                    txvar = sparta.input.variable.find(txstr);
                    if (txvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(txvar) == 0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (tystr != null)
                {
                    tyvar = sparta.input.variable.find(tystr);
                    if (tyvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(tyvar) == 0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (tzstr != null)
                {
                    tzvar = sparta.input.variable.find(tzstr);
                    if (tzvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(tzvar) == 0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
            }

            if (velflag != 0)
            {
                if (vxstr != null)
                {
                    vxvar = sparta.input.variable.find(vxstr);
                    if (vxvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.equal_style(vxvar)==0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (vystr != null)
                {
                    vyvar = sparta.input.variable.find(vystr);
                    if (vyvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.equal_style(vyvar)==0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (vzstr != null)
                {
                    vzvar = sparta.input.variable.find(vzstr);
                    if (vzvar < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.equal_style(vzvar)==0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (vstrx != null)
                {
                    vvarx = sparta.input.variable.find(vstrx);
                    if (vvarx < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(vvarx) == 0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (vstry != null)
                {
                    vvary = sparta.input.variable.find(vstry);
                    if (vvary < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(vvary) == 0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
                if (vstrz != null)
                {
                    vvarz = sparta.input.variable.find(vstrz);
                    if (vvarz < 0)
                        sparta.error.all("Variable name for create_particles does not exist");
                    if (sparta.input.variable.internal_style(vvarz) == 0)
                        sparta.error.all("Variable for create_particles is invalid style");
                }
            }

            // calculate Np if not set explicitly

            if (single!=0) np = 1;
            else if (np == 0)
            {
                Grid.ChildCell[] cells = sparta.grid.cells;
                Grid.ChildInfo[] cinfo = sparta.grid.cinfo;
                int nglocal = sparta.grid.nlocal;

                double flowvolme = 0.0;
                for (int icell = 0; icell < nglocal; icell++)
                {
                    if (cells[icell].nsplit > 1) continue;
                    if (cinfo[icell].type != (int)Enum1.INSIDE)
                        flowvolme += cinfo[icell].volume / cinfo[icell].weight;
                }
                double flowvol=0;
                sparta.mpi.MPI_Allreduce(ref flowvolme, ref flowvol, 1, MPI.MPI_DOUBLE, MPI.MPI_SUM, sparta.world);
                np = (long)(sparta.particle.mixture[imix].nrho * flowvol / sparta.update.fnum);
            }

            // generate particles
            // NOTE: invoke local or global option here

            if (sparta.comm.me == 0)
                if (sparta.screen!=null)
                {
                    fprintf(screen, "Creating particles ...\n");
                }

            sparta.mpi.MPI_Barrier(sparta.world);
            double time1 = sparta.mpi.MPI_Wtime();

            Int64 nprevious = sparta.particle.nglobal;
            if (single!=0) create_single();
            else if (globalflag==0)
            {
                if (twopass!=0) create_local_twopass(np);
                else create_local(np);
            } //else create_global(np);

            sparta.mpi.MPI_Barrier(sparta.world);
            double time2 = sparta.mpi.MPI_Wtime();

            // error check
            // only if no region and no variable species/density specified

            Int64 nglobal=0;
            Int64 nme = sparta.particle.nlocal;
            sparta.mpi.MPI_Allreduce(ref nme, ref nglobal, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);
            if (region==null && speciesflag == 0 && densflag == 0 && tempflag == 0 &&
                nglobal - nprevious != np)
            {
                string str=string.Format( "Created incorrect # of particles: {0} versus {1}" ,
                    nglobal - nprevious, np);
                sparta.error.all(str);
            }
            Int64 ncreated = nglobal - nprevious;
            sparta.particle.nglobal = nglobal;

            // print stats

            if (sparta.comm.me == 0)
            {
                if (screen)
                {
                    fprintf(screen, "Created " BIGINT_FORMAT " particles\n", ncreated);
                    fprintf(screen, "  CPU time = %g secs\n", time2 - time1);
                }
                if (logfile)
                {
                    fprintf(logfile, "Created " BIGINT_FORMAT " particles\n", ncreated);
                    fprintf(logfile, "  CPU time = %g secs\n", time2 - time1);
                }
            }
        }
        //public int evib(int);
        //public double erot(int);

        protected int imix, single, mspecies, twopass;
        protected double xp, yp, zp, vx, vy, vz;
        protected Region region;

        protected int speciesflag, densflag, velflag, tempflag, normflag;
        protected string sstr,sxstr,systr,szstr;
        protected string dstr,dxstr,dystr,dzstr;
        protected string tstr,txstr,tystr,tzstr;
        protected string vxstr,vystr,vzstr,vstrx,vstry,vstrz;
        protected int svar, sxvar, syvar, szvar;
        protected int dvar, dxvar, dyvar, dzvar;
        protected int tvar, txvar, tyvar, tzvar;
        protected int vxvar, vyvar, vzvar, vvarx, vvary, vvarz;
        protected string sxstr_copy,systr_copy,szstr_copy;
        protected string dxstr_copy,dystr_copy,dzstr_copy;
        protected string txstr_copy,tystr_copy,tzstr_copy;
        protected string vstrx_copy,vstry_copy,vstrz_copy;

        //protected virtual void create_single();
        //protected virtual void create_local(Int64);
        //protected void create_local_twopass(Int64);
        //protected int species_variable(double);
        //protected double density_variable(double, double);
        //protected double temperature_variable(double);
        //protected void velocity_variable(double, double, double);
        //protected int outside_region(int, double, double);
    }
}
