
#define MPI_STUBS


using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//using MPI_STATUS_IGNORE =null;

using MPI_Comm = System.Int32;
using MPI_Request = System.Int32;
using MPI_Datatype = System.Int32;
using MPI_Fint = System.Int32;
using MPI_Group = System.Int32;
using MPI_Offset = System.Int32;
using MPI_Op = System.Int32;

namespace cstest
{
    public class MPI
    {
        public const int MPI_COMM_WORLD = 0;

        public const int MPI_SUCCESS = 0;
        public const int MPI_ERR_ARG = -1;

        public const int MPI_INT = 1;
        public const int MPI_FLOAT = 2;
        public const int MPI_DOUBLE = 3;
        public const int MPI_CHAR = 4;
        public const int MPI_BYTE = 5;
        public const int MPI_LONG = 6;
        public const int MPI_LONG_LONG = 7;
        public const int MPI_DOUBLE_INT = 8;

        public const int MPI_SUM = 1;
        public const int MPI_MAX = 2;
        public const int MPI_MIN = 3;
        public const int MPI_MAXLOC = 4;
        public const int MPI_MINLOC = 5;
        public const int MPI_LOR = 6;

        public const int MPI_UNDEFINED = -1;
        public const int MPI_COMM_NULL = -1;
        public const int MPI_GROUP_EMPTY = -1;

        public const int MPI_ANY_SOURCE = -1;


        public int MPI_IN_PLACE;

        public const int MPI_MAX_PROCESSOR_NAME = 128;

        const int MAXEXTRA_DATATYPE = 16;

        public int nextra_datatype = new Int32();
        public MPI_Datatype[] ptr_datatype = new MPI_Datatype[MAXEXTRA_DATATYPE];
        public int[] index_datatype = new int[MAXEXTRA_DATATYPE];
        public int[] size_datatype = new int[MAXEXTRA_DATATYPE];

        static int _mpi_is_initialized = 0;


        public struct _MPI_Status
        {
            int MPI_SOURCE;
        }

        struct _mpi_double_int
        {
            Double value;
            int proc;
        }

        /* Function prototypes for MPI stubs */

        public int MPI_Init(string[] args)
        {
            if (_mpi_is_initialized > 0)
            {
                Console.WriteLine();
                Console.WriteLine("MPI Stub WARNING: MPI already initialized\n");
                return 1;
            }
            if (_mpi_is_initialized < 0)
            {
                Console.WriteLine("MPI Stub WARNING: MPI already finalized\n");
                return 1;
            }
            _mpi_is_initialized = 1;
            return 0;
        }
        public int MPI_Initialized(ref int flag)
        {
            flag = (_mpi_is_initialized > 0) ? 1 : 0;
            return 0;
        }
        public int MPI_Finalized(ref int flag)
        {
            flag = (_mpi_is_initialized < 0) ? 1 : 0;
            return 0;
        }
        public int MPI_Get_processor_name(string name, ref int resultlen)
        {
            const string host = "localhost";
            int len;

            if (name != null || resultlen != 0) return MPI_ERR_ARG;

            len = host.Length;
            name = String.Copy(host);
            resultlen = len;
            return MPI_SUCCESS;
        }
        public int MPI_Get_version(ref int major, ref int minor)
        {
            if (major != 0 || minor != 0) return MPI_ERR_ARG;

            major = 1;
            minor = 2;
            return MPI_SUCCESS;
        }
        public int MPI_Comm_rank(MPI_Comm comm, ref int me)
        {
            me = 0;
            return 0;
        }
        public int MPI_Comm_size(MPI_Comm comm, ref int nprocs)
        {
            nprocs = 1;
            return 0;
        }
        public int MPI_Abort(MPI_Comm comm, int errorcode)
        {
            Console.WriteLine("exit in MPI_Abort");
            //exit(1);
            return 0;
        }
        public int MPI_Finalize()
        {
            if (_mpi_is_initialized == 0)
            {
                Console.WriteLine("MPI Stub WARNING: MPI not yet initialized\n");
                return 1;
            }
            if (_mpi_is_initialized < 0)
            {
                Console.WriteLine("MPI Stub WARNING: MPI already finalized\n");
                return 1;
            }
            _mpi_is_initialized = -1;
            return 0;
        }
        public double MPI_Wtime()
        {
            return DateTime.Now.Millisecond;
        }

        //static 
        int stubtypesize(MPI_Datatype datatype)
        {
            if (datatype == MPI_INT) return sizeof(int);
            else if (datatype == MPI_FLOAT) return sizeof(float);
            else if (datatype == MPI_DOUBLE) return sizeof(double);
            else if (datatype == MPI_CHAR) return sizeof(char);
            else if (datatype == MPI_BYTE) return sizeof(char);
            else if (datatype == MPI_LONG) return sizeof(long);
            else if (datatype == MPI_LONG_LONG) return sizeof(System.UInt64);
#if unsafe
            else if (datatype == MPI_DOUBLE_INT) return sizeof(_mpi_double_int);
#endif      

            else
            {
                int i;
                for (i = 0; i < nextra_datatype; i++)
                {
                    if (datatype == index_datatype[i]) return size_datatype[i];
                }
            }
            return 0;
        }

        public int MPI_Type_size(MPI_Datatype datatype, ref int size)
        {
            if (size == 0) return MPI_ERR_ARG;

            size = stubtypesize(datatype);
            return 0;
        }
        //public int MPI_Send(const void* buf, int count, MPI_Datatype datatype,int dest, int tag, MPI_Comm comm)
        //{

        //}
        //public int MPI_Isend(const void* buf, int count, MPI_Datatype datatype,int source, int tag, MPI_Comm comm, MPI_Request* request)
        //{

        //}
        //public int MPI_Rsend(const void* buf, int count, MPI_Datatype datatype,int dest, int tag, MPI_Comm comm)
        //{

        //}
        //int MPI_Recv(void* buf, int count, MPI_Datatype datatype,
        //     int source, int tag, MPI_Comm comm, MPI_Status* status){}
        //int MPI_Irecv(void* buf, int count, MPI_Datatype datatype,
        //              int source, int tag, MPI_Comm comm, MPI_Request* request){}
        //int MPI_Wait(MPI_Request* request, MPI_Status* status){}
        //int MPI_Waitall(int n, MPI_Request* request, MPI_Status* status){}
        //int MPI_Waitany(int count, MPI_Request* request, int* index,
        //                MPI_Status* status){}
        //int MPI_Sendrecv(const void* sbuf, int scount, MPI_Datatype sdatatype,
        //          int dest, int stag, void* rbuf, int rcount,
        //          MPI_Datatype rdatatype, int source, int rtag,
        //          MPI_Comm comm, MPI_Status* status){}
        //int MPI_Get_count(MPI_Status* status, MPI_Datatype datatype, int* count){}

        public int MPI_Comm_split(MPI_Comm comm, int color, int key,ref MPI_Comm comm_out)
        {
            comm_out = comm;
            return 0;
        }
        //int MPI_Comm_dup(MPI_Comm comm, MPI_Comm* comm_out){}
        //int MPI_Comm_free(MPI_Comm* comm){}
        //MPI_Fint MPI_Comm_c2f(MPI_Comm comm){}
        //MPI_Comm MPI_Comm_f2c(MPI_Fint comm){}
        //int MPI_Comm_group(MPI_Comm comm, MPI_Group* group){}
        //int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm* newcomm){}
        //int MPI_Group_incl(MPI_Group group, int n, int* ranks, MPI_Group* newgroup){}

        //int MPI_Cart_create(MPI_Comm comm_old, int ndims, int* dims, int* periods,
        //                    int reorder, MPI_Comm* comm_cart){}
        //int MPI_Cart_get(MPI_Comm comm, int maxdims, int* dims, int* periods,
        //                 int* coords){}
        //int MPI_Cart_shift(MPI_Comm comm, int direction, int displ,
        //                   int* source, int* dest){}
        //int MPI_Cart_rank(MPI_Comm comm, int* coords, int* rank){}

        public int MPI_Type_contiguous(int count, MPI_Datatype oldtype,ref MPI_Datatype newtype)
        {
            if (nextra_datatype == MAXEXTRA_DATATYPE) return -1;
            ptr_datatype[nextra_datatype] = newtype;
            index_datatype[nextra_datatype] = -(nextra_datatype + 1);
            size_datatype[nextra_datatype] = count * stubtypesize(oldtype);
            nextra_datatype++;
            return 0;
        }
        public int MPI_Type_commit(ref MPI_Datatype datatype)
        {
            int i;
            for (i = 0; i < nextra_datatype; i++)
                if (datatype == ptr_datatype[i]) datatype = index_datatype[i];
            return 0;
        }
        //int MPI_Type_free(MPI_Datatype* datatype){}

        public delegate void MPI_User_function(IntPtr inn,IntPtr inout,ref int a,ref MPI_Datatype b);


        public int MPI_Op_create(MPI_User_function function, int commute,ref MPI_Op op)
        {
            return 0;
        }
        //int MPI_Op_free(MPI_Op* op){}

        public int MPI_Barrier(MPI_Comm comm)
        {
            return 0;
        }
        public int MPI_Bcast(System.IntPtr buf, int count, MPI_Datatype datatype,int root, MPI_Comm comm)
        {
            return 0;
        }
        //int MPI_Allreduce(void* sendbuf, void* recvbuf, int count,
        //                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){}
        //int MPI_Reduce(void* sendbuf, void* recvbuf, int count,
        //                   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){}
        //int MPI_Scan(void* sendbuf, void* recvbuf, int count,
        //             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){}
        //int MPI_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        //                  void* recvbuf, int recvcount, MPI_Datatype recvtype,
        //                  MPI_Comm comm){}
        //int MPI_Allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        //                   void* recvbuf, int* recvcounts, int* displs,
        //                   MPI_Datatype recvtype, MPI_Comm comm){}
        //int MPI_Reduce_scatter(void* sendbuf, void* recvbuf, int* recvcounts,
        //                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){}
        //int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        //               void* recvbuf, int recvcount, MPI_Datatype recvtype,
        //               int root, MPI_Comm comm){}
        //int MPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        //                void* recvbuf, int* recvcounts, int* displs,
        //                MPI_Datatype recvtype, int root, MPI_Comm comm){}
        //int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        //                void* recvbuf, int recvcount, MPI_Datatype recvtype,
        //                int root, MPI_Comm comm){}
        //int MPI_Scatterv(void* sendbuf, int* sendcounts, int* displs,
        //                 MPI_Datatype sendtype, void* recvbuf, int recvcount,
        //                 MPI_Datatype recvtype, int root, MPI_Comm comm){}
        //int MPI_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        //                 void* recvbuf, int recvcount, MPI_Datatype recvtype,
        //                 MPI_Comm comm){}
        //int MPI_Alltoallv(void* sendbuf, int* sendcounts, int* sdispls,
        //                  MPI_Datatype sendtype,
        //                  void* recvbuf, int* recvcounts, int* rdispls,
        //                  MPI_Datatype recvtype, MPI_Comm comm){}
        /* ---------------------------------------------------------------------- */
    }

}

