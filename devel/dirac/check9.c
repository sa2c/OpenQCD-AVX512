
/*******************************************************************************
*
* File check9.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Comparison of Dw_bnd() with Dw().
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
#include "sap.h"
#include "dirac.h"
#include "global.h"

typedef union
{
   weyl w;
   float r[12];
} spin_t;


static void blk_s2zero(int ic,spinor *s)
{
   int nb,isw;
   int nbh,n,nm,vol;
   block_t *b;

   b=blk_list(SAP_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;

   if (ic^isw)
      n=nbh;
   else
      n=0;

   nm=n+nbh;

   for (;n<nm;n++)
   {
      set_s2zero(vol,b[n].s[0]);
      assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,0,s);
   }
}


static void random_weyl(int vol,weyl *w)
{
   spin_t *ws,*wm;

   ws=(spin_t*)(w);
   wm=ws+vol;

   for (;ws<wm;ws++)
      gauss((*ws).r,12);
}


static void random_bnd(int ic)
{
   int bc,nb,isw;
   int nbh,n,nm;
   block_t *b;
   bndry_t *bb;

   bc=bc_type();
   b=blk_list(SAP_BLOCKS,&nb,&isw);
   nbh=nb/2;

   if (ic^isw)
      n=nbh;
   else
      n=0;

   nm=n+nbh;

   for (;n<nm;n++)
   {
      bb=b[n].bb;

      if ((cpr[0]==0)&&(b[n].bo[0]==0)&&(bc!=3))
         random_weyl(bb[0].vol,bb[0].w[0]);

      if ((cpr[0]==(NPROC0-1))&&((b[n].bo[0]+b[n].bs[0])==L0)&&(bc==0))
         random_weyl(bb[1].vol,bb[1].w[0]);
   }
}


int main(int argc,char *argv[])
{
   int my_rank,bc;
   int nb,isw,nbh,ic;
   int bs[4],n,nm,vol,ie;
   float mu,d,dmax;
   double phi[2],phi_prime[2],theta[3];
   spinor **ps;
   block_t *b;
   sw_parms_t swp;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check9.log","w",stdout);
      fin=freopen("check7.in","r",stdin);

      printf("\n");
      printf("Comparison of Dw_bnd() with Dw()\n");
      printf("--------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n",bs[0],bs[1],bs[2],bs[3]);


      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check9.c]",
                    "Syntax: check9 [-bc <type>]");
   }

   set_lat_parms(5.5,1.0,0,NULL,1.978);
   print_lat_parms();

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
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

   start_ranlux(0,1234);
   geometry();
   set_sap_parms(bs,0,1,1);
   alloc_bgr(SAP_BLOCKS);
   alloc_ws(4);

   swp=set_sw_parms(0.05);
   mu=0.123f;

   if (my_rank==0)
      printf("m0 = %.4e, mu = %.4e, csw = %.4e, cF = %.4e, cF' = %.4e\n\n",
             swp.m0,mu,swp.csw,swp.cF[0],swp.cF[1]);

   random_ud();
   set_ud_phase();
   sw_term(NO_PTS);

   assign_ud2u();
   assign_swd2sw();
   assign_ud2ubgr(SAP_BLOCKS);

   ps=reserve_ws(4);
   b=blk_list(SAP_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;

   ie=0;
   dmax=0.0f;

   for (ic=0;ic<2;ic++)
   {
      random_s(VOLUME,ps[0],1.0f);
      assign_s2s(VOLUME,ps[0],ps[3]);

      if (ic^isw)
         n=nbh;
      else
         n=0;
      nm=n+nbh;

      for (;n<nm;n++)
      {
         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],0);
         Dw_bnd(SAP_BLOCKS,n,0,0);
      }

      random_bnd(ic);
      random_s(VOLUME,ps[1],1.0f);
      assign_s2s(VOLUME,ps[1],ps[2]);
      sap_com(ic,ps[1]);
      mulr_spinor_add(VOLUME,ps[2],ps[1],-1.0f);

      blk_s2zero(ic^0x1,ps[0]);
      Dw(mu,ps[0],ps[1]);
      blk_s2zero(ic,ps[1]);

      mulr_spinor_add(VOLUME,ps[1],ps[2],-1.0f);
      d=norm_square(VOLUME,1,ps[1])/
         norm_square(VOLUME,1,ps[2]);

      if (d>dmax)
         dmax=d;

      if (ic^isw)
         n=nbh;
      else
         n=0;
      nm=n+nbh;

      for (;n<nm;n++)
      {
         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[3],0);
         Dw_bnd(SAP_BLOCKS,n,0,0);
         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],1);
         mulr_spinor_add(vol,b[n].s[0],b[n].s[1],-1.0f);
         d=norm_square(vol,0,b[n].s[0]);

         if (d!=0.0f)
            ie=1;
      }
   }

   error(ie,1,"main [check9.c]",
     "Dw_bnd() changes the input field where it should not");

   dmax=(float)(sqrt((double)(dmax)));

   if (my_rank==0)
   {
      printf("The maximal relative deviation is %.1e\n\n",dmax);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
