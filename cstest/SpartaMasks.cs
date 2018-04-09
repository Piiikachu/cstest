using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class SpartaMasks
    {
        // data masks

        public const uint EMPTY_MASK = 0x00000000;
        public const uint ALL_MASK = 0xffffffff;

        // particles

        public const uint PARTICLE_MASK = 0x00000001;
        public const uint SPECIES_MASK = 0x00000002;

        // grid

        public const uint CELL_MASK = 0x00000004;
        public const uint CINFO_MASK = 0x00000008;
        public const uint PCELL_MASK = 0x00000016;
        public const uint SINFO_MASK = 0x00000032;

        // collide

        public const uint VREMAX_MASK = 0x00000064;
        public const uint REMAIN_MASK = 0x00000128;

        // surf

        public const uint PT_MASK = 0x00000256;
        public const uint LINE_MASK = 0x00000512;
        public const uint TRI_MASK = 0x00001024;

    }
}
