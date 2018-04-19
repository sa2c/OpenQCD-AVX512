
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2007, 2008, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Exporting and importing spinor fields.
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

static const spinor_dble sd0={{{0.0}}};
static spinor_dble **psd;


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
   int my_rank,nsize,k;
   double d,dmax;
   char sfld_dir[NAME_SIZE],name[NAME_SIZE];
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check5.log","w",stdout);
      fin=freopen("check5.in","r",stdin);

      printf("\n");
      printf("Exporting and importing spinor fields\n");
      printf("-------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("sfld_dir","%s\n",sfld_dir);
      fclose(fin);
   }

   MPI_Bcast(sfld_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   start_ranlux(0,123456);
   geometry();
   alloc_wsd(6);
   psd=reserve_wsd(6);

   check_dir_root(sfld_dir);
   nsize=name_size("%s/testsfld%d",sfld_dir,6);
   error_root(nsize>=NAME_SIZE,1,"main [check5.c]","sfld_dir name is too long");

   for (k=0;k<3;k++)
   {
      random_sd(VOLUME,psd[k],1.0);
      sprintf(name,"%s/testsfld%d",sfld_dir,k);
      export_sfld(name,psd[k]);
   }

   for (k=0;k<3;k++)
   {
      sprintf(name,"%s/testsfld%d",sfld_dir,k);
      import_sfld(name,psd[k+3]);
      remove(name);
   }

   dmax=0.0;

   for (k=0;k<3;k++)
   {
      mulr_spinor_add_dble(VOLUME,psd[k],psd[k+3],-1.0);
      d=norm_square_dble(VOLUME,0,psd[k]);

      if (d>dmax)
         dmax=d;
   }

   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Exported 3 spinor fields to the directory\n"
             "%s\n",sfld_dir);
      printf("Then reimported and deleted them\n\n");
      printf("Maximal deviation = %.1e ",sqrt(dmax));
      printf("(should be exactly equal to 0.0)\n\n");
   }

   ptfld(4);
   sprintf(name,"%s/testsfld",sfld_dir);
   export_sfld(name,psd[4]);

   if (my_rank==0)
   {
      printf("Point source field exported to file\n"
             "%s\n\n",name);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
