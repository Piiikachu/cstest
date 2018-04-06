namespace cstest
{
    public class Region
    {

        public string id, style;
        public int interior;                     // 1 for interior, 0 for exterior
        public double extent_xlo, extent_xhi;     // bounding box on region
        public double extent_ylo, extent_yhi;
        public double extent_zlo, extent_zhi;
        public int bboxflag;                     // 1 if bounding box is computable

        public Region(SPARTA sparta, string[] arg)
        {
            id = string.Copy(arg[0]);
            style = string.Copy(arg[1]);
        }

        // called by other classes to check point versus region

        //       int match(double*);

        //       // implemented by each region, not called by other classes

        //       virtual int inside(double*) = 0;

        //protected:
        // void options(int, char**);
    }
}