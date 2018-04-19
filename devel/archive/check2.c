
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005, 2007, 2010-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Exporting and importing gauge configurations.
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
#include "uflds.h"
#include "su3fcts.h"
#include "linalg.h"
#include "archive.h"
#include "global.h"


static int cmp_ud(su3_dble *u,su3_dble *v)
{
   int it;

   it =((*u).c11.re!=(*v).c11.re);
   it|=((*u).c11.im!=(*v).c11.im);
   it|=((*u).c12.re!=(*v).c12.re);
   it|=((*u).c12.im!=(*v).c12.im);
   it|=((*u).c13.re!=(*v).c13.re);
   it|=((*u).c13.im!=(*v).c13.im);

   it|=((*u).c21.re!=(*v).c21.re);
   it|=((*u).c21.im!=(*v).c21.im);
   it|=((*u).c22.re!=(*v).c22.re);
   it|=((*u).c22.im!=(*v).c22.im);
   it|=((*u).c23.re!=(*v).c23.re);
   it|=((*u).c23.im!=(*v).c23.im);

   it|=((*u).c31.re!=(*v).c31.re);
   it|=((*u).c31.im!=(*v).c31.im);
   it|=((*u).c32.re!=(*v).c32.re);
   it|=((*u).c32.im!=(*v).c32.im);
   it|=((*u).c33.re!=(*v).c33.re);
   it|=((*u).c33.im!=(*v).c33.im);

   return it;
}


static int check_ud(su3_dble *usv)
{
   int it;
   su3_dble *u,*um;

   u=udfld();
   um=u+4*VOLUME;
   it=0;

   for (;u<um;u++)
   {
      it|=cmp_ud(u,usv);
      usv+=1;
   }

   return it;
}


int main(int argc,char *argv[])
{
   int my_rank,bc,nsize,ie;
   double phi[2],phi_prime[2],theta[3];
   su3_dble *udb,**usv;
   char cnfg_dir[NAME_SIZE],cnfg[NAME_SIZE];
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      fin=freopen("check2.in","r",stdin);

      printf("\n");
      printf("Exporting and importing gauge configurations\n");
      printf("--------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("cnfg_dir","%s\n",cnfg_dir);
      fclose(fin);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }

   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);

   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.5;
   theta[1]=1.0;
   theta[2]=-0.5;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(1);

   start_ranlux(0,123456);
   geometry();
   alloc_wud(1);

   check_dir_root(cnfg_dir);
   nsize=name_size("%s/testcnfg",cnfg_dir);
   error_root(nsize>=NAME_SIZE,1,"main [check2.c]","cnfg_dir name is too long");

   if (my_rank==0)
   {
      printf("Export random field configurations to the file\n"
             "%s/testcnfg.\n",cnfg_dir);
      printf("Then read the fields from there and compare with the saved "
             "fields.\n\n");
   }

   udb=udfld();
   usv=reserve_wud(1);
   random_ud();
   cm3x3_assign(4*VOLUME,udb,usv[0]);

   sprintf(cnfg,"%s/testcnfg",cnfg_dir);
   export_cnfg(cnfg);

   random_ud();
   import_cnfg(cnfg);

   ie=(check_bc(0.0)^0x1);
   ie|=check_ud(usv[0]);
   error(ie!=0,1,"main [check2.c]","The gauge field is not properly restored");
   print_flags();

   if (my_rank==0)
   {
      printf("No errors detected --- the fields are correctly exported\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
