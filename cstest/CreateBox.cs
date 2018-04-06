using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class CreateBox
    {
        private SPARTA sparta;
        public CreateBox(SPARTA sparta)
        {
            this.sparta = sparta;
        }

        public void command(int narg,string[] args)
        {
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);
            if (sparta.domain.box_exist!=0)
            {
                sparta.error.all("Cannot create_box after simulation box is defined");
            }
            //if (sparta.domain.dimension == 2 && sparta.domain.zperiodic == 0)
            //  sparta.error.all("Cannot run 2d simulation with nonperiodic Z dimension");

            sparta.domain.box_exist = 1;

            if (narg != 6) sparta.error.all( "Illegal create_box command");

            sparta.domain.boxlo[0] = double.Parse(arg[0]);
            sparta.domain.boxhi[0] = double.Parse(arg[1]);
            sparta.domain.boxlo[1] = double.Parse(arg[2]);
            sparta.domain.boxhi[1] = double.Parse(arg[3]);
            sparta.domain.boxlo[2] = double.Parse(arg[4]);
            sparta.domain.boxhi[2] = double.Parse(arg[5]);

            if (sparta.domain.dimension == 2)
            {
                if (sparta.domain.boxlo[2] >= 0.0 || sparta.domain.boxhi[2] <= 0.0)
                    sparta.error.all(
                       "Create_box z box bounds must straddle 0.0 for 2d simulations");
            }
            if (sparta.domain.axisymmetric!=0 && sparta.domain.boxlo[1] != 0.0)
                sparta.error.all( "Box ylo must be 0.0 for axi-symmetric model");

            // problem setup using info from header

            sparta.update.ntimestep = 0;

            sparta.domain.print_box("Created ");
            sparta.domain.set_initial_box();
            sparta.domain.set_global_box();
        }

    }
}
