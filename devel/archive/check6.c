
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2007, 2008, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Importing a previously exported spinor field.
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
#include "lattice.h"
#include "sflds.h"
#include "linalg.h"
#include "archive.h"
#include "global.h"

static spinor_dble **psd;
static const spinor_dble sd0={{{0.0}}};


static void ptfld(int k)
{
   int x0,x1,x2,x3,y0,y1,y2,y3,ix;
   spinor_dble *s;

   y0=L0*cpr[0];
   y1=L1*cpr[1];
   y2=L2*cpr[2];
   y3=L3*cpr[3];

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
               s=psd[k]+ix;

               (*s)=sd0;
               (*s).c1.c1.re=(double)(y0+x0);
               (*s).c2.c1.re=(double)(y1+x1);
               (*s).c3.c1.re=(double)(y2+x2);
               (*s).c4.c1.re=(double)(y3+x3);
            }
         }
      }
   }
}


int main(int argc,char *argv[])
{
   int my_rank,nsize;
   double d;
   char sfld_dir[NAME_SIZE],name[NAME_SIZE];
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check6.log","w",stdout);
      fin=freopen("check6.in","r",stdin);

      printf("\n");
      printf("Importing a previously exported spinor field\n");
      printf("--------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("sfld_dir","%s\n",sfld_dir);
      fclose(fin);
   }

   start_ranlux(0,123456);
   geometry();
   alloc_wsd(2);
   psd=reserve_wsd(2);
   ptfld(0);

   check_dir_root(sfld_dir);
   nsize=name_size("%s/testsfld",sfld_dir);
   error_root(nsize>=NAME_SIZE,1,"main [check6.c]","sfld_dir name is too long");
   sprintf(name,"%s/testsfld",sfld_dir);

   import_sfld(name,psd[1]);

   mulr_spinor_add_dble(VOLUME,psd[0],psd[1],-1.0);
   d=norm_square_dble(VOLUME,1,psd[0]);

   if (my_rank==0)
   {
      printf("Imported field from file\n"
             "%s\n\n",name);
      printf("Deviation = %.1e ",sqrt(d));
      printf("(should be exactly equal to 0.0)\n\n");
      remove(name);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
