
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2009-2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the program tcharge_slices().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "forces.h"
#include "wflow.h"
#include "tcharge.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int bc,n,dn;
static double eps,Q1,Q2,Q[N0],Q0[N0];


int main(int argc,char *argv[])
{
   int my_rank,i,imax,t;
   double phi[2],phi_prime[2],theta[3];
   double nplaq,act,dev;
   FILE *fin=NULL,*flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      fin=freopen("check3.in","r",stdin);

      printf("\n");
      printf("Check of the program tcharge_slices()\n");
      printf("-------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("n","%d",&n);
      read_line("dn","%d",&dn);
      read_line("eps","%lf",&eps);
      fclose(fin);

      printf("n = %d\n",n);
      printf("dn = %d\n",dn);
      printf("eps = %.2e\n\n",eps);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dn,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&eps,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   set_lat_parms(6.0,1.0,0,NULL,1.0);

   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(0);

   start_ranlux(0,123456);
   geometry();
   alloc_wfd(2);

   if (bc==0)
      nplaq=(double)(6*N0-6)*(double)(N1*N2*N3);
   else
      nplaq=(double)(6*N0)*(double)(N1*N2*N3);

   random_ud();
   imax=n/dn;

   for (i=0;i<imax;i++)
   {
      fwd_euler(dn,eps);

      act=action0(1)/nplaq;
      Q1=tcharge();
      Q2=tcharge_slices(Q);
      dev=fabs(Q1-Q2);

      for (t=0;t<N0;t++)
      {
         Q2-=Q[t];
         Q0[t]=Q[t];
      }

      dev+=fabs(Q2);

      if (my_rank==0)
      {
         printf("n=%3d, act=%.4e, Q=% .2e, dev=%.1e, Q[0...%d]=% .2e",
                (i+1)*dn,act,Q1,dev,N0-1,Q[0]);

         for (t=1;t<N0;t++)
            printf(", % .2e",Q[t]);

         printf("\n");
      }

      MPI_Bcast(Q0,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);

      for (t=0;t<N0;t++)
      {
         if ((Q[t]-Q0[t])!=0.0)
            break;
      }

      error(t!=N0,1,"main [check3.c]",
            "Charge slices are not globally the same");
   }

   if (my_rank==0)
   {
      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
