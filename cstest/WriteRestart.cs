﻿using System;
using System.IO;

namespace cstest
{
    public class WriteRestart
    {
        public const string MAGIC_STRING = "SpartA RestartT";
        public const int ENDIAN = 0x0001;
        public const int ENDIANSWAP = 0x1000;
        public const int VERSION_NUMERIC = 0;

        enum Enum1{
            VERSION, SMALLINT, CELLINT, BIGINT,
            UNITS, NTIMESTEP, NPROCS,
            FNUM, NRHO, VSTREAM, TEMP_THERMAL, GRAVITY, SURFMAX, GRIDCUT, GRID_WEIGHT,
            COMM_SORT, COMM_STYLE,
            DIMENSION, AXISYMMETRIC, BOXLO, BOXHI, BFLAG,
            NPARTICLE, NUNSPLIT, NSPLIT, NSUB, NPOINT, NSURF,
            SPECIES, MIXTURE, PARTICLE_CUSTOM, GRID, SURF,
            MULTIPROC, PROCSPERFILE, PERPROC
        };    // new fields added after PERPROC
        private SPARTA sparta;
        public WriteRestart(SPARTA sparta)
        {
            this.sparta = sparta;
            sparta.mpi.MPI_Comm_rank(sparta.world, ref me);
            sparta.mpi.MPI_Comm_size(sparta.world, ref nprocs);
            multiproc = 0;
        }
        //public void command(int, char**);
        //public void multiproc_options(int, int, char**);
        public void write(string file)
        {
            // open single restart file or base file for multiproc case

            if (me == 0)
            {
                string hfile=null;
                
                if (multiproc!=0)
                {
                    if (file.Contains("%"))
                    {
                        string[] tmpstr = file.Split('%');
                        hfile = string.Format(" {0} {1} {2}", tmpstr[0], "base", tmpstr[1]);
                    }
                    //hfile = new char[strlen(file) + 16];
                    //char* ptr = strchr(file, '%');
                    //*ptr = '\0';
                    //sprintf(hfile, "%s%s%s", file, "base", ptr + 1);
                    //*ptr = '%';
                }
                else hfile = file;
                fp = new FileStream(hfile, FileMode.OpenOrCreate, FileAccess.Write);
                //fp = fopen(hfile, "wb");
                if (fp == null)
                {
                    string str=string.Format( "Cannot open restart file {0}", hfile);
                    sparta.error.one(str);
                }
                //if (multiproc!=0) delete[] hfile;
            }

            // proc 0 writes magic string, endian flag, numeric version

            if (me == 0)
            {
                magic_string();
                endian();
                version_numeric();
            }

            // proc 0 writes header info
            // also simulation box, particle species, parent grid cells, surf info

            Int64 btmp = sparta.particle.nlocal;
            sparta.mpi.MPI_Allreduce(ref btmp, ref sparta.particle.nglobal, 1, MPI.MPI_LONG_LONG, MPI.MPI_SUM, sparta.world);

            if (me == 0)
            {
                header();
                box_params();
                particle_params();
                grid_params();
                surf_params();
            }

            // communication buffer for my per-proc info = child grid cells and particles
            // max_size = largest buffer needed by any proc

            int send_size = grid->size_restart();
            send_size += sparta.particle.size_restart();

            int max_size;
            MPI_Allreduce(&send_size, &max_size, 1, MPI_INT, MPI_MAX, world);

            char* buf;
            memory->create(buf, max_size, "write_restart:buf");
            memset(buf, 0, max_size);

            // all procs write file layout info which may include per-proc sizes

            file_layout(send_size);

            // header info is complete
            // if multiproc output:
            //   close header file, open multiname file on each writing proc,
            //   write PROCSPERFILE into new file

            if (multiproc)
            {
                if (me == 0) fclose(fp);

                char* multiname = new char[strlen(file) + 16];
                char* ptr = strchr(file, '%');
                *ptr = '\0';
                sprintf(multiname, "%s%d%s", file, icluster, ptr + 1);
                *ptr = '%';

                if (filewriter)
                {
                    fp = fopen(multiname, "wb");
                    if (fp == NULL)
                    {
                        char str[128];
                        sprintf(str, "Cannot open restart file %s", multiname);
                        sparta.error.one(str);
                    }
                    write_int(PROCSPERFILE, nclusterprocs);
                }

                delete[] multiname;
            }

            // pack my child grid and particle data into buf

            int n = grid->pack_restart(buf);
            n += sparta.particle.pack_restart(&buf[n]);

            // output of one or more native files
            // filewriter = 1 = this proc writes to file
            // ping each proc in my cluster, receive its data, write data to file
            // else wait for ping from fileproc, send my data to fileproc

            int tmp, recv_size;
            MPI_Status status;
            MPI_Request request;

            if (filewriter)
            {
                for (int iproc = 0; iproc < nclusterprocs; iproc++)
                {
                    if (iproc)
                    {
                        MPI_Irecv(buf, max_size, MPI_CHAR, me + iproc, 0, world, &request);
                        MPI_Send(&tmp, 0, MPI_INT, me + iproc, 0, world);
                        MPI_Wait(&request, &status);
                        MPI_Get_count(&status, MPI_CHAR, &recv_size);
                    }
                    else recv_size = send_size;

                    write_char_vec(PERPROC, recv_size, buf);
                }
                fclose(fp);

            }
            else
            {
                MPI_Recv(&tmp, 0, MPI_INT, fileproc, 0, world, &status);
                MPI_Rsend(buf, send_size, MPI_CHAR, fileproc, 0, world);
            }

            // clean up

            memory->destroy(buf);
        }


        int me, nprocs;
        FileStream fp;

        int multiproc;             // 0 = proc 0 writes for all
                                   // else # of procs writing files
        int nclusterprocs;         // # of procs in my cluster that write to one file
        int filewriter;            // 1 if this proc writes a file, else 0
        int fileproc;              // ID of proc in my cluster who writes to file
        int icluster;              // which cluster I am in

        void header()
        {
            write_string(VERSION, universe->version);
            write_int(SMALLINT, sizeof(smallint));
            write_int(CELLINT, sizeof(cellint));
            write_int(BIGINT, sizeof(bigint));
            write_string(UNITS, update->unit_style);
            write_bigint(NTIMESTEP, update->ntimestep);
            write_int(NPROCS, nprocs);

            write_double(FNUM, update->fnum);
            write_double(NRHO, update->nrho);
            write_double_vec(VSTREAM, 3, update->vstream);
            write_double(TEMP_THERMAL, update->temp_thermal);
            write_double_vec(GRAVITY, 3, update->gravity);
            write_int(SURFMAX, grid->maxsurfpercell);
            write_double(GRIDCUT, grid->cutoff);
            write_int(COMM_SORT, comm->commsortflag);
            write_int(COMM_STYLE, comm->commpartstyle);
            write_int(GRID_WEIGHT, grid->cellweightflag);

            write_bigint(NPARTICLE, particle->nglobal);
            write_bigint(NUNSPLIT, grid->nunsplit);
            write_int(NSPLIT, grid->nsplit);
            write_int(NSUB, grid->nsub);
            write_int(NPOINT, surf->npoint);
            if (domain->dimension == 2) write_int(NSURF, surf->nline);
            else write_int(NSURF, surf->ntri);

            // -1 flag signals end of header

            int flag = -1;
            fwrite(&flag, sizeof(int), 1, fp);
        }
        //void box_params();
        //void particle_params();
        //void grid_params();
        //void surf_params();
        //void file_layout(int);

        void magic_string()
        {
            using (StreamWriter sw=new StreamWriter(fp))
            {
                sw.Write(MAGIC_STRING);
                
                sw.Close();
            }
            //int n = strlen(MAGIC_STRING) + 1;
            //char* str = new char[n];
            //strcpy(str, MAGIC_STRING);
            //fwrite(str, sizeof(char), n, fp);
            //delete[] str;
        }
        void endian()
        {
            int endian = ENDIAN;
            //fwrite(&endian, sizeof(int), 1, fp);
            using (StreamWriter sw = new StreamWriter(fp))
            {
                sw.Write(endian);

                sw.Close();
            }
        }
        void version_numeric()
        {
            int vn = VERSION_NUMERIC;
            //fwrite(&vn, sizeof(int), 1, fp);
            using (StreamWriter sw = new StreamWriter(fp))
            {
                sw.Write(vn);

                sw.Close();
            }
        }

        //void write_int(int, int);
        //void write_bigint(int, Int64);
        //void write_double(int, double);
        //void write_string(int, char*);
        //void write_int_vec(int, int, int*);
        //void write_double_vec(int, int, double*);
        //void write_char_vec(int, int, char*);
    }
}