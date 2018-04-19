
/*******************************************************************************
*
* File time2.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of the SAP preconditioner.
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
#include "sap.h"
#include "global.h"


int main(int argc,char *argv[])
{
   int my_rank,bc,count,nt;
   int ncy,nmr,bs[4];
   int n,ie;
   float mu;
   double phi[2],phi_prime[2],theta[3];
   double rbb,wt1,wt2,wdt;
   spinor **ps;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("time2.log","w",stdout);
      fin=freopen("time2.in","r",stdin);

      printf("\n");
      printf("Timing of the SAP preconditioner\n");
      printf("--------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

#if (defined x64)
#if (defined AVX)
#if (defined FMA3)
   printf("Using AVX and FMA3 instructions\n");
#else
   printf("Using AVX instructions\n");
#endif
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

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      read_line("ncy","%d",&ncy);
      read_line("nmr","%d",&nmr);
      fclose(fin);

      printf("bs = %d %d %d %d\n",bs[0],bs[1],bs[2],bs[3]);
      printf("ncy = %d\n",ncy);
      printf("nmr = %d\n\n",nmr);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [time2.c]",
                    "Syntax: time2 [-bc <type>]");
   }

   set_lat_parms(5.5,1.0,0,NULL,1.978);
   print_lat_parms();

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncy,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.34;
   theta[1]=-1.25;
   theta[2]=0.58;
   set_bc_parms(bc,1.0,1.0,0.9012,1.2034,phi,phi_prime,theta);
   print_bc_parms(2);

   start_ranlux(0,12345);
   geometry();
   alloc_ws(3);
   set_sap_parms(bs,0,1,1);
   alloc_bgr(SAP_BLOCKS);

   set_sw_parms(0.0123);
   mu=0.0785f;
   rbb=2.0*(1.0/(double)(bs[0])+1.0/(double)(bs[1])+
            1.0/(double)(bs[2])+1.0/(double)(bs[3]));

   random_ud();
   set_ud_phase();
   sw_term(NO_PTS);
   assign_ud2ubgr(SAP_BLOCKS);
   assign_swd2swbgr(SAP_BLOCKS,NO_PTS);

   ps=reserve_ws(3);
   random_s(VOLUME,ps[2],1.0f);
   bnd_s2zero(ALL_PTS,ps[2]);
   normalize(VOLUME,1,ps[2]);

   nt=(int)(2.0e6/(double)(ncy*nmr*VOLUME));
   if (nt<2)
      nt=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<nt;count++)
      {
         set_s2zero(VOLUME,ps[0]);
         assign_s2s(VOLUME,ps[2],ps[1]);

         for (n=0;n<ncy;n++)
            sap(mu,0,nmr,ps[0],ps[1]);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      nt*=2;
   }

   wdt=2.0e6*wdt/((double)(nt)*(double)(ncy*VOLUME));

   if (my_rank==0)
   {
      printf("Using the MinRes block solver:\n");
      printf("Time per lattice point:   %.3f micro sec\n",(double)(ncy)*wdt);
      printf("Time per point and cycle: %.3f micro sec",wdt);
      printf(" (about %d Mflops)\n\n",
             (int)(((double)(nmr*2256+24)+112.0*rbb)/wdt));
   }

   ie=assign_swd2swbgr(SAP_BLOCKS,ODD_PTS);
   error_root(ie,1,"main [time2.c]",
              "The inversion of the SW term was not safe");

   nt=(int)(2.0e6/(double)(ncy*nmr*VOLUME));
   if (nt<2)
      nt=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<nt;count++)
      {
         set_s2zero(VOLUME,ps[0]);
         assign_s2s(VOLUME,ps[2],ps[1]);

         for (n=0;n<ncy;n++)
            sap(mu,1,nmr,ps[0],ps[1]);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      nt*=2;
   }

   wdt=2.0e6*wdt/((double)(nt)*(double)(ncy*VOLUME));

   if (my_rank==0)
   {
      printf("Using the even-odd preconditioned MinRes block solver:\n");
      printf("Time per lattice point:   %.3f micro sec\n",(double)(ncy)*wdt);
      printf("Time per point and cycle: %.3f micro sec",wdt);
      printf(" (about %d Mflops)\n\n",
             (int)(((double)((nmr+1)*2076+48)+112.0*rbb)/wdt));
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
