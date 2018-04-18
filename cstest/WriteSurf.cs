using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class WriteSurf
    {
        public const int MAXLINE = 256;

        private SPARTA sparta;

        public WriteSurf(SPARTA sparta)
        {
            this.sparta = sparta;
        }
        //public void command(int, char**);
        public void write_file(FileStream fp)
        {
            int dim = sparta.domain.dimension;

            List<Surf.Point> pts = sparta.surf.pts;
            List<Surf.Line> lines = sparta.surf.lines;
            List<Surf.Tri> tris = sparta.surf.tris;

            int npoint = sparta.surf.npoint;
            int nline = sparta.surf.nline;
            int ntri = sparta.surf.ntri;

            // header section
            StreamWriter sw = new StreamWriter(fp);
            
            sw.WriteLine("# Surface element file written by SPARTA\n\n");
            sw.WriteLine("{0} points\n", npoint);
            if (dim == 2) sw.WriteLine("{0} lines\n", nline);
            else sw.WriteLine("{0} triangles\n", ntri);
            sw.WriteLine("\n");

            // points

            sw.WriteLine("Points\n\n");
            if (dim == 2)
            {
                for (int i = 0; i < npoint; i++)
                    sw.WriteLine("{0} %20.15g %20.15g\n", i + 1, pts[i].x[0], pts[i].x[1]);
            }
            else
            {
                for (int i = 0; i < npoint; i++)
                    sw.WriteLine("{0} %20.15g %20.15g %20.15g\n", i + 1,
                        pts[i].x[0], pts[i].x[1], pts[i].x[2]);
            }

            // lines

            if (dim == 2)
            {
                sw.WriteLine("\nLines\n\n");
                for (int i = 0; i < nline; i++)
                    sw.WriteLine("{0} {1} {2} {3}\n", i + 1, lines[i].type,
                        lines[i].p1 + 1, lines[i].p2 + 1);
            }

            // triangles

            if (dim == 3)
            {
                sw.WriteLine("\nTriangles\n\n");
                for (int i = 0; i < ntri; i++)
                    sw.WriteLine("{0} {1} {2} {3} {4}\n", i + 1, tris[i].type,
                        tris[i].p1 + 1, tris[i].p2 + 1, tris[i].p3 + 1);
            }
        }
    }
}
