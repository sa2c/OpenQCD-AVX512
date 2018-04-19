
/*******************************************************************************
*
* File time1.c
*
* Copyright (C) 2005, 2008, 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of Dw() and Dwhat().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "global.h"


int main(int argc,char *argv[])
{
   int my_rank,bc,count,nt;
   int i,nflds;
   float mu;
   double phi[2],phi_prime[2],theta[3];
   double wt1,wt2,wdt;
   spinor **ps;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("time1.log","w",stdout);

      printf("\n");
      printf("Timing of Dw() and Dwhat()\n");
      printf("--------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      if (NPROC>1)
         printf("There are %d MPI processes\n",NPROC);
      else
         printf("There is 1 MPI process\n");

      if ((VOLUME*sizeof(float))<(64*1024))
      {
         printf("The local size of the gauge field is %d KB\n",
                (int)((72*VOLUME*sizeof(float))/(1024)));
         printf("The local size of a quark field is %d KB\n",
                (int)((24*VOLUME*sizeof(float))/(1024)));
      }
      else
      {
         printf("The local size of the gauge field is %d MB\n",
                (int)((72*VOLUME*sizeof(float))/(1024*1024)));
         printf("The local size of a quark field is %d MB\n",
                (int)((24*VOLUME*sizeof(float))/(1024*1024)));
      }

#if (defined x64)
#if (defined AVX)
      printf("Using AVX instructions\n");
#else
      printf("Using SSE3 instructions and 16 xmm registers\n");
#endif
#if (defined P3)
      printf("Assuming SSE prefetch instructions fetch 32 bytes\n");
#elif (defined PM)
      printf("Assuming SSE prefetch instructions fetch 64 bytes\n");
#elif (defined P4)
      printf("Assuming SSE prefetch instructions fetch 128 bytes\n");
#else
      printf("SSE prefetch instructions are not used\n");
#endif
#endif
      printf("\n");

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [time1.c]",
                    "Syntax: time1 [-bc <type>]");
   }

   set_lat_parms(5.5,1.0,0,NULL,1.978);
   print_lat_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.35;
   theta[1]=-1.25;
   theta[2]=0.78;
   set_bc_parms(bc,0.55,0.78,0.9012,1.2034,phi,phi_prime,theta);
   print_bc_parms(2);

   start_ranlux(0,12345);
   geometry();

   set_sw_parms(-0.0123);
   mu=0.0785f;

   random_ud();
   set_ud_phase();
   sw_term(NO_PTS);
   assign_ud2u();
   assign_swd2sw();

   nflds=(int)((4*1024*1024)/(VOLUME*sizeof(float)))+1;
   if ((nflds%2)==1)
      nflds+=1;
   alloc_ws(nflds);
   ps=reserve_ws(nflds);

   for (i=0;i<nflds;i++)
      random_s(VOLUME,ps[i],1.0f);

   nt=(int)(1.0e6f/(double)(nflds*VOLUME));
   if (nt<2)
      nt=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<nt;count++)
      {
         for (i=0;i<nflds;i+=2)
            Dw(mu,ps[i],ps[i+1]);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      nt*=2;
   }

   wdt=4.0e6*wdt/((double)(nt)*(double)(nflds*VOLUME));

   if (my_rank==0)
   {
      printf("Time per lattice point for Dw():\n");
      printf("%4.3f micro sec (%d Mflops)\n\n",wdt,(int)(1920.0/wdt));
   }

   nt=(int)(1.0e6/(double)(nflds*VOLUME));
   if (nt<2)
      nt=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<nt;count++)
      {
         for (i=0;i<nflds;i+=2)
            Dwhat(mu,ps[i],ps[i+1]);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      nt*=2;
   }

   wdt=4.0e6*wdt/((double)(nt)*(double)(nflds*VOLUME));

   if (my_rank==0)
   {
      printf("Time per lattice point for Dwhat():\n");
      printf("%4.3f micro sec (%d Mflops)\n\n",wdt,(int)(1908.0/wdt));
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
