
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005, 2008, 2012, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Direct test of the Schwarz alternating procedure.
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
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "block.h"
#include "dirac.h"
#include "sap.h"
#include "global.h"


int main(int argc,char *argv[])
{
   int my_rank,bc;
   int n,ie,itm;
   int bs[4],nmr;
   float mu,res,del[3];
   double phi[2],phi_prime[2],theta[3];
   spinor **ps;
   tm_parms_t tm;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      fin=freopen("check2.in","r",stdin);

      printf("\n");
      printf("Direct test of the Schwarz alternating procedure\n");
      printf("------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      read_line("mu","%f",&mu);
      read_line("nmr","%d",&nmr);
      fclose(fin);

      printf("bs = %d %d %d %d\n",bs[0],bs[1],bs[2],bs[3]);
      printf("mu = %.3e\n",mu);
      printf("nmr = %d\n\n",nmr);
      fflush(flog);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }

   set_lat_parms(5.5,1.0,0,NULL,1.978);
   print_lat_parms();

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&mu,1,MPI_FLOAT,0,MPI_COMM_WORLD);
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
   set_sap_parms(bs,0,1,1);
   alloc_bgr(SAP_BLOCKS);
   alloc_ws(4);
   ps=reserve_ws(4);

   set_sw_parms(0.05);

   for (itm=0;itm<2;itm++)
   {
      if (itm==0)
         set_tm_parms(1);
      else
         set_tm_parms(0);

      random_ud();
      set_ud_phase();
      sw_term(NO_PTS);
      assign_ud2u();
      assign_swd2sw();
      assign_ud2ubgr(SAP_BLOCKS);
      assign_swd2swbgr(SAP_BLOCKS,NO_PTS);

      set_s2zero(VOLUME,ps[0]);
      random_s(VOLUME,ps[1],1.0f);
      bnd_s2zero(ALL_PTS,ps[1]);
      normalize(VOLUME,1,ps[1]);
      assign_s2s(VOLUME,ps[1],ps[2]);

      if (my_rank==0)
      {
         tm=tm_parms();
         printf("Twisted-mass flag = %d\n",tm.eoflg);
         printf("MinRes block solver:\n");
      }

      for (n=0;n<8;n++)
      {
         sap(mu,0,nmr,ps[0],ps[1]);
         res=norm_square(VOLUME,1,ps[1]);
         res=(float)(sqrt((double)(res)));

         if (my_rank==0)
            printf("n = %d: \t residue = %.2e\t ",n+1,res);

         Dw(mu,ps[0],ps[3]);
         mulr_spinor_add(VOLUME,ps[3],ps[2],-1.0f);
         mulr_spinor_add(VOLUME,ps[3],ps[1],1.0f);
         del[0]=norm_square(VOLUME,1,ps[3]);
         del[0]=(float)(sqrt((double)(del[0])));

         assign_s2s(VOLUME,ps[0],ps[3]);
         bnd_s2zero(ALL_PTS,ps[3]);
         mulr_spinor_add(VOLUME,ps[3],ps[0],-1.0f);
         del[1]=norm_square(VOLUME,1,ps[3]);
         del[1]=(float)(sqrt((double)(del[1])));

         assign_s2s(VOLUME,ps[1],ps[3]);
         bnd_s2zero(ALL_PTS,ps[3]);
         mulr_spinor_add(VOLUME,ps[3],ps[1],-1.0f);
         del[2]=norm_square(VOLUME,1,ps[3]);
         del[2]=(float)(sqrt((double)(del[1])));

         if (my_rank==0)
            printf("check = %.2e, bnd checks = %.1e,%.1e\n",
                   del[0],del[1],del[2]);
      }

      ie=assign_swd2swbgr(SAP_BLOCKS,ODD_PTS);
      error_root(ie,1,"main [check2.c]",
                 "The inversion of the SW term was not safe");

      set_s2zero(VOLUME,ps[0]);
      random_s(VOLUME,ps[1],1.0f);
      bnd_s2zero(ALL_PTS,ps[1]);
      normalize(VOLUME,1,ps[1]);
      assign_s2s(VOLUME,ps[1],ps[2]);

      if (my_rank==0)
      {
         printf("\n");
         printf("Even-odd preconditioned MinRes block solver:\n");
      }

      for (n=0;n<8;n++)
      {
         sap(mu,1,nmr,ps[0],ps[1]);
         res=norm_square(VOLUME,1,ps[1]);
         res=(float)(sqrt((double)(res)));

         if (my_rank==0)
            printf("n = %d: \t residue = %.2e\t ",n+1,res);

         Dw(mu,ps[0],ps[3]);
         mulr_spinor_add(VOLUME,ps[3],ps[2],-1.0f);
         mulr_spinor_add(VOLUME,ps[3],ps[1],1.0f);
         del[0]=norm_square(VOLUME,1,ps[3]);
         del[0]=(float)(sqrt((double)(del[0])));

         assign_s2s(VOLUME,ps[0],ps[3]);
         bnd_s2zero(ALL_PTS,ps[3]);
         mulr_spinor_add(VOLUME,ps[3],ps[0],-1.0f);
         del[1]=norm_square(VOLUME,1,ps[3]);
         del[1]=(float)(sqrt((double)(del[1])));

         assign_s2s(VOLUME,ps[1],ps[3]);
         bnd_s2zero(ALL_PTS,ps[3]);
         mulr_spinor_add(VOLUME,ps[3],ps[1],-1.0f);
         del[2]=norm_square(VOLUME,1,ps[3]);
         del[2]=(float)(sqrt((double)(del[1])));

         if (my_rank==0)
            printf("check = %.2e, bnd checks = %.1e,%.1e\n",
                   del[0],del[1],del[2]);
      }

      if (my_rank==0)
         printf("\n");
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
