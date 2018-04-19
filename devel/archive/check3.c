
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005, 2007, 2008, 2010-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Importing a configuration previously exported by check2.
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


static double avg_plaq(void)
{
   double plaq;

   plaq=plaq_sum_dble(1);

   return plaq/((double)(6*NPROC)*(double)(VOLUME));
}


int main(int argc,char *argv[])
{
   int my_rank,bc,nsize,ir,ie;
   stdint_t l[4];
   double phi[2],phi_prime[2],theta[3];
   double plaq0,plaq1,plaq2;
   char cnfg_dir[NAME_SIZE],cnfg[NAME_SIZE];
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      fin=freopen("check3.in","r",stdin);

      printf("\n");
      printf("Importing gauge fields exported by check2\n");
      printf("-----------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("cnfg_dir","%s\n",cnfg_dir);
      fclose(fin);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>]");
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
   print_bc_parms(0);

   start_ranlux(0,9876);
   geometry();
   random_ud();
   plaq0=avg_plaq();

   check_dir_root(cnfg_dir);
   nsize=name_size("%s/testcnfg",cnfg_dir);
   error_root(nsize>=NAME_SIZE,1,"main [check3.c]","cnfg_dir name is too long");
   sprintf(cnfg,"%s/testcnfg",cnfg_dir);

   if (my_rank==0)
   {
      fin=fopen(cnfg,"rb");
      error_root(fin==NULL,1,"main [check3.c]","Unable to open input file");

      ir=fread(l,sizeof(stdint_t),4,fin);
      ir+=fread(&plaq1,sizeof(double),1,fin);
      error_root(ir!=5,1,"main [check3.c]","Incorrect read count");
      fclose(fin);

      if (endianness()==BIG_ENDIAN)
      {
         bswap_int(4,l);
         bswap_double(1,&plaq1);
      }

      printf("Random gauge field, average plaquette = %.15e\n\n",plaq0);
      printf("Now read gauge field from file\n"
             "%s:\n",cnfg);
      printf("%dx%dx%dx%d lattice\n",
             (int)(l[0]),(int)(l[1]),(int)(l[2]),(int)(l[3]));
      printf("Average plaquette = %.15e\n",plaq1);
   }

   import_cnfg(cnfg);
   ie=check_bc(0.0);
   error(ie!=1,1,"main [check3.c]","Boundary conditions are not preserved");
   plaq2=avg_plaq();

   if (my_rank==0)
   {
      printf("Should be         = %.15e\n\n",plaq2);
      remove(cnfg);
   }

   print_flags();

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
