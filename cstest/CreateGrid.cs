using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class CreateGrid
    {
        enum Enum1{ NONE, LEVEL, STRIDE, CLUMP, BLOCK, RANDOM };
        enum Enum2 { XYZ, XZY, YXZ, YZX, ZXY, ZYX };
        enum Enum3 { ANY, ALL };
        private SPARTA sparta;
        public CreateGrid(SPARTA sparta)
        {
            this.sparta = sparta;
        }
        public void command(int narg, string[] args)
        {
            string[] arg = new string[narg];
            Array.Copy(args, 1, arg, 0, narg);
            if (sparta.domain.box_exist ==0 )
                sparta.error.all("Cannot create grid before simulation box is defined");
            if (sparta.grid.exist != 0)
                sparta.error.all("Cannot create grid when grid is already defined");

            sparta.grid.exist = 1;

            if (narg < 3) sparta.error.all("Illegal create_grid command");

            int nx = int.Parse(arg[0]);
            int ny = int.Parse(arg[1]);
            int nz = int.Parse(arg[2]);

            if (nx < 1 || ny < 1 || nz < 1)
                sparta.error.all("Illegal create_grid command");
            if (sparta.domain.dimension == 2 && nz != 1)
                sparta.error.all("Create_grid nz value must be 1 for a 2d simulation");

            // optional args

            dimension = sparta.domain.dimension;

            int nlevels = 1;
            int bstyle = (int)Enum1.NONE;
            int px = 0;
            int py = 0;
            int pz = 0;
            int order;
            int inside = (int)Enum3.ANY;

            int iarg = 3;
            while (iarg < narg)
            {
                if (string.Equals(arg[iarg], "level"))
                {
                    if (iarg + 8 > narg) sparta.error.all( "Illegal create_grid command");
                    if (bstyle != (int)Enum1.NONE && bstyle != (int)Enum1.LEVEL)
                        sparta.error.all( "Illegal create_grid command");
                    bstyle = (int)Enum1.LEVEL;
                    if (int.Parse(arg[iarg + 1]) != nlevels + 1)
                        sparta.error.all( "Illegal create_grid command");
                    nlevels++;
                    iarg += 8;

                }
                else if (string.Equals(arg[iarg], "region"))
                {
                    if (iarg + 6 > narg) sparta.error.all( "Illegal create_grid command");
                    if (bstyle !=(int)Enum1.NONE && bstyle !=(int)Enum1.LEVEL)
                        sparta.error.all( "Illegal create_grid command");
                    bstyle =(int)Enum1.LEVEL;
                    if (int.Parse(arg[iarg + 1]) != nlevels + 1)
                        sparta.error.all( "Illegal create_grid command");
                    if (sparta.domain.find_region(arg[iarg + 2]) < 0)
                        sparta.error.all( "Create_grid region ID does not exist");
                    nlevels++;
                    iarg += 6;

                }
                else if (string.Equals(arg[iarg], "stride"))
                {
                    if (iarg + 2 > narg) sparta.error.all( "Illegal create_grid command");
                    if (bstyle !=(int)Enum1.NONE) sparta.error.all( "Illegal create_grid command");
                    bstyle = (int)Enum1.STRIDE;
                    if (string.Equals(arg[iarg + 1], "xyz")  ) order = (int)Enum2.XYZ;
                    else if (string.Equals(arg[iarg + 1], "xzy")  ) order = (int)Enum2.XZY;
                    else if (string.Equals(arg[iarg + 1], "yxz")  ) order = (int)Enum2.YXZ;
                    else if (string.Equals(arg[iarg + 1], "yzx")  ) order = (int)Enum2.YZX;
                    else if (string.Equals(arg[iarg + 1], "zxy")  ) order = (int)Enum2.ZXY;
                    else if (string.Equals(arg[iarg + 1], "zyx")  ) order = (int)Enum2.ZYX;
                    else sparta.error.all( "Illegal create_grid command");
                    iarg += 2;

                }
                else if (string.Equals(arg[iarg], "clump"))
                {
                    if (iarg + 2 > narg) sparta.error.all( "Illegal create_grid command");
                    if (bstyle !=(int)Enum1.NONE) sparta.error.all( "Illegal create_grid command");
                    bstyle = (int)Enum1.CLUMP;
                    if (string.Equals(arg[iarg + 1], "xyz")  ) order = (int)Enum2.XYZ;
                    else if (string.Equals(arg[iarg + 1], "xzy")  ) order = (int)Enum2.XZY;
                    else if (string.Equals(arg[iarg + 1], "yxz")  ) order = (int)Enum2.YXZ;
                    else if (string.Equals(arg[iarg + 1], "yzx")  ) order = (int)Enum2.YZX;
                    else if (string.Equals(arg[iarg + 1], "zxy")  ) order = (int)Enum2.ZXY;
                    else if (string.Equals(arg[iarg + 1], "zyx")  ) order = (int)Enum2.ZYX;
                    else sparta.error.all( "Illegal create_grid command");
                    iarg += 2;

                }
                else if (string.Equals(arg[iarg], "block")  )
                {
                    if (iarg + 4 > narg) sparta.error.all( "Illegal create_grid command");
                    if (bstyle !=(int)Enum1.NONE) sparta.error.all( "Illegal create_grid command");
                    bstyle = (int)Enum1.BLOCK;
                    if (string.Equals(arg[iarg + 1], "*")  ) px = 0;
                    else px = int.Parse(arg[iarg + 1]);
                    if (string.Equals(arg[iarg + 2], "*")  ) py = 0;
                    else py = int.Parse(arg[iarg + 2]);
                    if (string.Equals(arg[iarg + 3], "*")  ) pz = 0;
                    else pz = int.Parse(arg[iarg + 3]);
                    iarg += 4;

                }
                else if (string.Equals(arg[iarg], "random")  )
                {
                    if (iarg + 1 > narg) sparta.error.all( "Illegal create_grid command");
                    if (bstyle !=(int)Enum1.NONE) sparta.error.all( "Illegal create_grid command");
                    bstyle = (int)Enum1.RANDOM;
                    iarg += 1;

                }
                else if (string.Equals(arg[iarg], "inside")  )
                {
                    if (iarg + 2 > narg) sparta.error.all( "Illegal create_grid command");
                    if (string.Equals(arg[iarg + 1], "any")  ) inside = (int)Enum3.ANY;
                    else if (string.Equals(arg[iarg + 1], "all")  ) inside = (int)Enum3.ALL;
                    else sparta.error.all( "Illegal create_grid command");
                    iarg += 2;

                }
                else sparta.error.all( "Illegal create_grid command");




            }
        }


        private int dimension;

        //private void bounds(string str, int nmax, ref int nlo, ref int nhi)
        //{

        //}
        //      private int cell_in_region(double*, double*, class Region *, int);
        //void procs2grid(int, int, int, int &, int &, int &);
    }
}
