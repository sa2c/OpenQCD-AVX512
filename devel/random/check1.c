
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of import/export functions for the ranlux generators
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
#include "global.h"

#define NRAN 10000

static float r[2*NRAN];
static double rd[2*NRAN];


int main(int argc,char *argv[])
{
   int my_rank,tag,k,ie,ied;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Check of import/export functions for the ranlux generators\n");
      printf("----------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,1234);
   ranlxs(r,NRAN);
   ranlxd(rd,NRAN);
   tag=98029;

   export_ranlux(tag,"check1.dat");
   ranlxs(r,NRAN);
   ranlxd(rd,NRAN);

   k=import_ranlux("check1.dat");
   error (k!=tag,1,"main [check1.c]",
          "Import_ranlux() returns incorrect tag");

   ranlxs(r+NRAN,NRAN);
   ranlxd(rd+NRAN,NRAN);

   ie=0;
   ied=0;

   for (k=0;k<NRAN;k++)
   {
      ie|=(r[k]!=r[NRAN+k]);
      ied|=(rd[k]!=rd[NRAN+k]);
   }

   error((ie!=0)||(ied!=0),1,"main [check1.c]",
         "Export/import of the generator states failed");

   if (my_rank==0)
   {
      remove("check1.dat");
      printf("No errors detected --- all is fine\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
