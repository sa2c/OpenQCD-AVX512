
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2007, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Writing and reading spinor fields.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "sflds.h"
#include "linalg.h"
#include "archive.h"
#include "global.h"


int main(int argc,char *argv[])
{
   int my_rank,nsize,k;
   double d,dmax;
   spinor_dble **psd;
   char loc_dir[NAME_SIZE],name[NAME_SIZE];
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);
      fin=freopen("check4.in","r",stdin);

      printf("\n");
      printf("Writing and reading spinor fields\n");
      printf("---------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("loc_dir","%s\n",loc_dir);
      fclose(fin);
   }

   MPI_Bcast(loc_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   start_ranlux(0,123456);
   geometry();
   alloc_wsd(6);
   psd=reserve_wsd(6);

   check_dir(loc_dir);
   nsize=name_size("%s/testsfld_%d.%d",loc_dir,NPROC,6);
   error_root(nsize>=NAME_SIZE,1,"main [check4.c]","loc_dir name is too long");

   for (k=0;k<3;k++)
   {
      random_sd(VOLUME,psd[k],1.0);
      sprintf(name,"%s/testsfld_%d.%d",loc_dir,my_rank,k);
      write_sfld(name,psd[k]);
   }

   for (k=0;k<3;k++)
   {
      sprintf(name,"%s/testsfld_%d.%d",loc_dir,my_rank,k);
      read_sfld(name,psd[k+3]);
      remove(name);
   }

   dmax=0.0;

   for (k=0;k<3;k++)
   {
      mulr_spinor_add_dble(VOLUME,psd[k],psd[k+3],-1.0);
      d=norm_square_dble(VOLUME,1,psd[k]);

      if (d>dmax)
         dmax=d;
   }

   if (my_rank==0)
   {
      printf("Wrote 3 spinor fields to the files\n"
             "%s/testsfld_*\n"
             "on the local disks. ",loc_dir);
      printf("Then read the fields from there and removed\n"
             "the files.\n\n");
      printf("Maximal deviation = %.1e ",sqrt(dmax));
      printf("(should be exactly equal to 0.0)\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
