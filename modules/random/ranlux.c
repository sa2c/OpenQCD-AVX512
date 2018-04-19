
/*******************************************************************************
*
* File ranlux.c
*
* Copyright (C) 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Initialization of the ranlux generators
*
* The externally accessible functions are
*
*   void start_ranlux(int level,int seed)
*     Initializes the random number generators ranlxs and ranlxd on all
*     processes in different ways. The luxury level should be 0 (recommended)
*     or 1 (exceptional) and the seed can be any positive integer less than
*     or equal to INT_MAX/NPROC. An error occurs if the seed is not in this
*     range.
*
*   void export_ranlux(int tag,char *out)
*     Writes the tag, the lattice sizes, the process grid and the state of
*     the random number generators ranlxs and ranlxd to the file "out" from
*     process 0. The state of the generators is retrieved from all processes
*     and written to the file in the order specified in the notes.
*
*   int import_ranlux(char *in)
*     Reads the state of the random number generators ranlxs and ranlxd from
*     the file "in". The file is read from process 0 only and the data on
*     the file are expected in the form written by export_ranlux(). An error
*     occurs if the lattice sizes and process grid read from the file do not
*     coincide with the actual values of these parameters. The program then
*     resets the generators on all processes to the state read from the file.
*     The value returned is the tag read from the file.
*
* Notes:
*
* The program start_ranlux() guarantees that all generators are initialized
* with different seed values. Moreover, the initialization is guaranteed to
* be pairwise different from the previous one when start_ranlux() is called
* a second time with another value of "seed".
*
* The functions in this module assign a lexicographic index
*
*  id=n3+NPROC3*n2+NPROC2*NPROC3*n1+NPROC1*NPROC2*NPROC3*n0
*
* to the MPI process with Cartesian grid coordinates (n0,n1,n2,n3). On a
* given process, the initialization of the generators depends only on the
* index and the parameters "level" and "seed".
*
* The function export_ranlux() writes the state of the generators to the
* specified file in the order of increasing index. Independently of the
* machine, the data are written in little-endian byte order, using a 4 byte
* integer type. The import function assumes the data on the input file to
* be of this kind and converts them to big-endian byte order if the machine
* is big endian.
*
* The programs in this module act globally and must be called simultaneously
* on all processes.
*
*******************************************************************************/

#define RANLUX_C

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "mpi.h"
#include "utils.h"
#include "lattice.h"
#include "random.h"
#include "global.h"

static int *rlxs_state=NULL,*rlxd_state;
static stdint_t *state;


static int check_machine(void)
{
   int np,ie;

   MPI_Comm_size(MPI_COMM_WORLD,&np);

   error_root(np!=NPROC,1,"check_machine [ranlux.c]",
              "Actual number of processes does not match NPROC");

   error_root(sizeof(stdint_t)!=4,1,"check_machine [ranlux.c]",
              "Size of a stdint_t integer is not 4");

   ie=endianness();
   error_root(ie==UNKNOWN_ENDIAN,1,"check_machine [ranlux.c]",
              "Unkown endianness");

   return ie;
}


static int alloc_state(void)
{
   int nlxs,nlxd,n;

   nlxs=rlxs_size();
   nlxd=rlxd_size();
   n=nlxs+nlxd;

   if (rlxs_state==NULL)
   {
      rlxs_state=malloc(n*sizeof(int));
      rlxd_state=rlxs_state+nlxs;
      state=malloc(n*sizeof(stdint_t));
      error((rlxs_state==NULL)||(state==NULL),1,"alloc_state [ranlux.c]",
            "Unable to allocate state arrays");
   }

   return n;
}


static int get_ip(int n)
{
   int np[4];

   np[3]=n%NPROC3;
   n/=NPROC3;
   np[2]=n%NPROC2;
   n/=NPROC2;
   np[1]=n%NPROC1;
   n/=NPROC1;
   np[0]=n;

   return ipr_global(np);
}


static void get_state(void)
{
   rlxs_get(rlxs_state);
   rlxd_get(rlxd_state);
}


static void reset_state(void)
{
   rlxs_reset(rlxs_state);
   rlxd_reset(rlxd_state);
}


void start_ranlux(int level,int seed)
{
   int my_rank,max_seed,loc_seed;
   int n,iprms[2];

   if (NPROC>1)
   {
      iprms[0]=level;
      iprms[1]=seed;

      MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

      error((iprms[0]!=level)||(iprms[1]!=seed),1,
            "start_ranlux [ranlux.c]","Input parameters are not global");
   }

   max_seed=INT_MAX/NPROC;

   error_root((level<0)||(level>1)||(seed<1)||(seed>max_seed),1,
              "start_ranlux [ranlux.c]","Parameters are out of range");

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   loc_seed=0;

   for (n=0;n<NPROC;n++)
   {
      if (get_ip(n)==my_rank)
         loc_seed=seed+n*max_seed;
   }

   rlxs_init(level,loc_seed);
   rlxd_init(level+1,loc_seed);
}


void export_ranlux(int tag,char *out)
{
   int my_rank,ie,iw;
   int dmy,tag0,tag1;
   int n,ip,ns,k;
   stdint_t lsize[9];
   MPI_Status stat;
   FILE *fout=NULL;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   dmy=1;
   tag0=mpi_tag();
   tag1=mpi_tag();
   iw=0;
   ie=check_machine();
   ns=alloc_state();

   if (my_rank==0)
   {
      fout=fopen(out,"wb");
      error_root(fout==NULL,1,"export_ranlux [ranlux.c]",
                 "Unable to open output file");

      lsize[0]=(stdint_t)(tag);
      lsize[1]=(stdint_t)(NPROC0);
      lsize[2]=(stdint_t)(NPROC1);
      lsize[3]=(stdint_t)(NPROC2);
      lsize[4]=(stdint_t)(NPROC3);

      lsize[5]=(stdint_t)(L0);
      lsize[6]=(stdint_t)(L1);
      lsize[7]=(stdint_t)(L2);
      lsize[8]=(stdint_t)(L3);

      if (ie==BIG_ENDIAN)
         bswap_int(9,lsize);

      iw+=fwrite(lsize,sizeof(stdint_t),9,fout);
   }

   for (n=0;n<NPROC;n++)
   {
      ip=get_ip(n);

      if (ip>0)
      {
         if (my_rank==0)
         {
            MPI_Send(&dmy,1,MPI_INT,ip,tag0,MPI_COMM_WORLD);
            MPI_Recv(rlxs_state,ns,MPI_INT,ip,tag1,MPI_COMM_WORLD,&stat);
         }
         else if (my_rank==ip)
         {
            get_state();
            MPI_Recv(&dmy,1,MPI_INT,0,tag0,MPI_COMM_WORLD,&stat);
            MPI_Send(rlxs_state,ns,MPI_INT,0,tag1,MPI_COMM_WORLD);
         }
      }
      else if (my_rank==0)
         get_state();

      if (my_rank==0)
      {
         for (k=0;k<ns;k++)
            state[k]=(stdint_t)(rlxs_state[k]);

         if (ie==BIG_ENDIAN)
            bswap_int(ns,state);

         iw+=fwrite(state,sizeof(stdint_t),ns,fout);
      }
   }

   if (my_rank==0)
   {
      error_root(iw!=(9+NPROC*ns),1,"export_ranlux [ranlux.c]",
                 "Incorrect write count");
      fclose(fout);
   }
}


int import_ranlux(char *in)
{
   int my_rank,ie,ir;
   int dmy,tag0,tag1;
   int n,ip,ns,k;
   stdint_t lsize[9];
   MPI_Status stat;
   FILE *fin=NULL;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   dmy=1;
   tag0=mpi_tag();
   tag1=mpi_tag();
   ir=0;
   ie=check_machine();
   ns=alloc_state();

   if (my_rank==0)
   {
      fin=fopen(in,"rb");
      error_root(fin==NULL,1,"import_ranlux [ranlux.c]",
                 "Unable to open input file");

      ir+=fread(lsize,sizeof(stdint_t),9,fin);

      if (ie==BIG_ENDIAN)
         bswap_int(9,lsize);

      n=0;
      n|=(lsize[1]!=(stdint_t)(NPROC0));
      n|=(lsize[2]!=(stdint_t)(NPROC1));
      n|=(lsize[3]!=(stdint_t)(NPROC2));
      n|=(lsize[4]!=(stdint_t)(NPROC3));

      n|=(lsize[5]!=(stdint_t)(L0));
      n|=(lsize[6]!=(stdint_t)(L1));
      n|=(lsize[7]!=(stdint_t)(L2));
      n|=(lsize[8]!=(stdint_t)(L3));

      error_root(n!=0,1,"import_ranlux [ranlux.c]",
                 "Lattice size or process grid mismatch");
   }

   for (n=0;n<NPROC;n++)
   {
      ip=get_ip(n);

      if (my_rank==0)
      {
         ir+=fread(state,sizeof(stdint_t),ns,fin);

         if (ie==BIG_ENDIAN)
            bswap_int(ns,state);

         for (k=0;k<ns;k++)
            rlxs_state[k]=(int)(state[k]);
      }

      if (ip>0)
      {
         if (my_rank==0)
         {
            MPI_Send(rlxs_state,ns,MPI_INT,ip,tag1,MPI_COMM_WORLD);
            MPI_Recv(&dmy,1,MPI_INT,ip,tag0,MPI_COMM_WORLD,&stat);
         }
         else if (my_rank==ip)
         {
            MPI_Recv(rlxs_state,ns,MPI_INT,0,tag1,MPI_COMM_WORLD,&stat);
            MPI_Send(&dmy,1,MPI_INT,0,tag0,MPI_COMM_WORLD);
            reset_state();
         }
      }
      else if (my_rank==0)
         reset_state();
   }

   if (my_rank==0)
   {
      n=(int)(lsize[0]);
      error_root(ir!=(9+NPROC*ns),1,"import_ranlux [ranlux.c]",
                 "Incorrect read count");
      fclose(fin);
   }

   MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);

   return n;
}
