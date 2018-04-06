using System;
using System.Collections.Generic;

namespace cstest
{
    public class MyPage<T> where T :struct,IComparable
    {
        public int ndatum;      // total # of stored datums
        public int nchunk;      // total # of stored chunks
        public int errorflag;   // flag > 1 if error has occurred
                                // 1 = invalid inputs
                                // 2 = memory allocation error
                                // 3 = chunk size exceeded maxchunk
        private int maxchunk;   // max # of datums in one requested chunk
        private int pagesize;   // # of datums in one page, default = 1024
        private int pagedelta;  // # of pages to allocate at once, default = 1
         
        private List<T[]> pages;      // list of allocated pages
        private T[] page;        // ptr to current page
        private int npage;      // # of allocated pages
        private int ipage;      // index of current page
        private int index;      // current index on current page
         
        private void allocate()
        {
            npage += pagedelta;
            //pages = (T**)realloc(pages, npage * sizeof(T*));
            pages = new List<T[]>(npage);
            if (pages!=null)
            {
                errorflag = 2;
                return;
            }
            for (int i = npage - pagedelta; i < npage; i++)
            {
                //pages[i] = (T*)malloc(pagesize * sizeof(T));
                pages[i] = new T[pagesize];
                if (pages[i]!=null) errorflag = 2;
            }
        }

        public MyPage(int user_maxchunk = 1, int user_pagesize = 1024,
               int user_pagedelta = 1)
        {
            maxchunk = user_maxchunk;
            pagesize = user_pagesize;
            pagedelta = user_pagedelta;

            errorflag = 0;
            if (maxchunk <= 0 || pagesize <= 0 || pagedelta <= 0) errorflag = 1;
            if (maxchunk > pagesize) errorflag = 1;
            if (errorflag!=0) return;

            // initial page allocation

            ndatum = nchunk = 0;
            pages = null;
            npage = 0;
            allocate();
            if (errorflag!=0) return;
            ipage = index = 0;
            page = pages[ipage];
        }

        public T get()
        {
            ndatum++;
            nchunk++;
            if (index < pagesize) return page[index++];
            ipage++;
            if (ipage == npage)
            {
                allocate();
                if (errorflag!=0) return default(T);
            }
            page = pages[ipage];
            index = 0;
            return page[index++];
        }

        // get ptr to location that can store N datums
        // error if N > maxchunk
        // return null if run out of memory

        public T get(int n)
        {
            if (n > maxchunk)
            {
                errorflag = 3;
                return default(T);
            }
            ndatum += n;
            nchunk++;
            if (index + n <= pagesize)
            {
                int start = index;
                index += n;
                return page[start];
            }
            ipage++;
            if (ipage == npage)
            {
                allocate();
                if (errorflag!=0) return default(T);
            }
            page = pages[ipage];
            index = n;
            return page[0];
        }

        // get ptr to location that can store maxchunk datums
        // will return same ptr as previous call if vgot() not called
        // return null if run out of memory

        public T vget()
        {
            if (index + maxchunk <= pagesize) return page[index];
            ipage++;
            if (ipage == npage)
            {
                allocate();
                if (errorflag!=0) return default(T);
            }
            page = pages[ipage];
            index = 0;
            return page[index];
        }

        // increment by N = # of values stored in loc returned by vget()
        // OK to not call if vget() ptr was not used
        // error if N > maxchunk

        public void vgot(int n)
        {
            if (n > maxchunk) errorflag = 3;
            ndatum += n;
            nchunk++;
            index += n;
        }

        // reset index to beginning of first page
        // effectively clears all pages, without freeing any memory

        public void reset()
        {
            ndatum = nchunk = 0;
            index = ipage = 0;
            page = pages[ipage];
        }

        // return total size of all allocated pages

        public int size()
        {
            Console.WriteLine("mypage.size()->fake size");
            return npage * pagesize * sizeof(int); //sizeof(T);
        }
    }
}