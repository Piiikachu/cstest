using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public partial class Grid
    {
        //      public int pack_one(int, char*, int, int, int);
        public int unpack_one(string buf, int ownflag, int molflag)
        {

        }
        //      public int pack_one_adapt(char*, char*, int);
        //      public int pack_particles(int, char*, int);
        //      public int unpack_particles(char*, int);
        //      public void unpack_particles_adapt(int, char*);
        //      public void compress();

        static void unpack_ghosts(int nsize, string buf)
        {
            int n = 0;
            while (n < nsize)
            {
                n += unpack_one(, 0, 0);
            }
        }
    }
}
