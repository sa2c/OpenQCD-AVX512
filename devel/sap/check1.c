
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the block solver programs.
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
   int nb,isw,ie,itm;
   int bs[4],n,k,vol,volh;
   float mu,res0,res[8],res_max[8];
   double phi[2],phi_prime[2],theta[3];
   spinor **ps;
   block_t *b;
   tm_parms_t tm;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Check of the block solver programs\n");
      printf("----------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n\n",bs[0],bs[1],bs[2],bs[3]);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>]");
   }

   set_lat_parms(5.5,1.0,0,NULL,1.978);
   print_lat_parms();

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
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

   start_ranlux(0,1234);
   geometry();
   set_sap_parms(bs,0,1,1);
   alloc_bgr(SAP_BLOCKS);
   alloc_ws(4);

   set_sw_parms(0.05);
   mu=0.123f;
   ps=reserve_ws(4);

   for (itm=0;itm<2;itm++)
   {
      if (itm==1)
         set_tm_parms(1);

      random_ud();
      set_ud_phase();
      sw_term(NO_PTS);
      assign_ud2ubgr(SAP_BLOCKS);
      assign_swd2swbgr(SAP_BLOCKS,NO_PTS);

      b=blk_list(SAP_BLOCKS,&nb,&isw);
      vol=(*b).vol;
      volh=vol/2;

      for (k=0;k<8;k++)
         res_max[k]=0.0f;

      random_s(VOLUME,ps[0],1.0f);
      bnd_s2zero(ALL_PTS,ps[0]);
      set_s2zero(VOLUME,ps[1]);

      for (n=0;n<nb;n++)
      {
         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],1);
         res0=norm_square(vol,0,b[n].s[1]);

         for (k=0;k<8;k++)
         {
            blk_mres(n,mu,4);

            assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[1],2);
            mulr_spinor_add(vol,b[n].s[2],b[n].s[0],1.0f);
            assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,2,ps[1]);

            Dw_blk(SAP_BLOCKS,n,mu,2,0);
            assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],2);
            mulr_spinor_add(vol,b[n].s[2],b[n].s[0],-1.0f);
            res[k]=norm_square(vol,0,b[n].s[2])/res0;

            if (res[k]>res_max[k])
               res_max[k]=res[k];
         }
      }

      if (NPROC>1)
      {
         MPI_Reduce(res_max,res,8,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);

         for (k=0;k<8;k++)
            res_max[k]=res[k];
      }

      if (my_rank==0)
      {
         tm=tm_parms();
         printf("Twisted-mass flag = %d\n",tm.eoflg);
         printf("Check of blk_mres():\n");

         for (k=0;k<8;k++)
            printf("nmr = %2d, res_max = %.1e\n",
                   4*(k+1),sqrt((double)(res_max[k])));
      }

      for (k=0;k<8;k++)
         res_max[k]=0.0f;

      ie=assign_swd2swbgr(SAP_BLOCKS,ODD_PTS);
      error_root(ie,1,"main [check1.c]",
                 "The inversion of the SW term was not safe");

      random_s(VOLUME,ps[0],1.0f);
      bnd_s2zero(ALL_PTS,ps[0]);
      set_s2zero(VOLUME,ps[1]);

      for (n=0;n<nb;n++)
      {
         assign_s2sblk(SAP_BLOCKS,n,EVEN_PTS,ps[0],1);
         res0=norm_square(volh,0,b[n].s[1]);

         for (k=0;k<8;k++)
         {
            blk_eo_mres(n,mu,3);

            assign_s2sblk(SAP_BLOCKS,n,EVEN_PTS,ps[1],2);
            mulr_spinor_add(volh,b[n].s[2],b[n].s[0],1.0f);
            assign_sblk2s(SAP_BLOCKS,n,EVEN_PTS,2,ps[1]);

            Dwhat_blk(SAP_BLOCKS,n,mu,2,0);
            assign_s2sblk(SAP_BLOCKS,n,EVEN_PTS,ps[0],2);
            mulr_spinor_add(volh,b[n].s[2],b[n].s[0],-1.0f);
            res[k]=norm_square(volh,0,b[n].s[2])/res0;

            if (res[k]>res_max[k])
               res_max[k]=res[k];
         }
      }

      if (NPROC>1)
      {
         MPI_Reduce(res_max,res,8,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);

         for (k=0;k<8;k++)
            res_max[k]=res[k];
      }

      if (my_rank==0)
      {
         printf("Check of blk_eo_mres():\n");

         for (k=0;k<8;k++)
            printf("nmr = %2d, res_max = %.1e\n",
                   3*(k+1),sqrt((double)(res_max[k])));

         printf("\n");
      }
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
