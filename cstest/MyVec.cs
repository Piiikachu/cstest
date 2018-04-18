using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace cstest
{
    public class MyVec<T>
    {

        public int n;        // current length of vector

        // construct an empty vector
        // no need for destructor, since STL vec is a class member

        public MyVec()
        {
            n = nmax = 0;
        }

        // access an element of vector
        // caller must insure index is in range 0 to Nmax-1
        // inlined so performance should be same as STL vector access

        //inline T& operator[] (int index) {return vec[index];}

        // resize vector to nnew, only if necessary

        void grow(int nnew)
        {
            if (nnew <= nmax) return;
            vec=new List<T>(nnew);
            nmax = nnew;
        }


        private List<T> vec;    // underlying STL vector
        private int nmax;              // length of STL vector allocation
    }
}
