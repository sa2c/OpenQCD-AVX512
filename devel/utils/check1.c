
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005, 2008, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Copying of files. After running this program, one can verify that all
* bytes have been copied correctly using the diff utility.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "utils.h"
#include "archive.h"
#include "global.h"

#define NRAN 10000

static float r[NRAN];


int main(int argc,char *argv[])
{
   int my_rank,n,iw;
   FILE *flog=NULL,*fdat=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Copying of .log and .dat files from process 0\n");
      printf("---------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,1234);
   ranlxs(r,NRAN);

   if (my_rank==0)
   {
      printf("Write 10 random numbers to test.log (in asci format)\n");
      printf("and %d numbers to test.dat (in binary format).\n\n",NRAN);

      fdat=fopen("test.log","w");
      error_root(fdat==NULL,1,"main [check1.c]","Unable to open file");
      iw=0;
      for (n=0;n<10;n++)
         iw|=(fprintf(fdat,"r[%d] = %.6e\n",n,r[n])<0);
      error_root(iw!=0,1,"main [check1.c]","Incorrect write count");
      fclose(fdat);

      fdat=fopen("test.dat","wb");
      error_root(fdat==NULL,1,"main [check1.c]","Unable to open file");
      iw=fwrite(&r[0],sizeof(float),NRAN,fdat);
      error_root(iw!=NRAN,1,"main [check1.c]","Incorrect write count");
      fclose(fdat);

      printf("Copy the files to test.log~ and test.dat~, respectively.\n");
      printf("The copying may then be verified using the diff utility.\n\n");

      copy_file("test.log","test.log~");
      copy_file("test.dat","test.dat~");
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
