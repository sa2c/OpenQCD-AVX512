
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2007, 2008, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs in the module dfl_subspace.c.
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
#include "flags.h"
#include "lattice.h"
#include "block.h"
#include "linalg.h"
#include "sflds.h"
#include "vflds.h"
#include "dfl.h"
#include "global.h"


static void check_basis(int Ns,double *dev0,double *dev1)
{
   int nb,isw,i,j;
   double dev,x[2],y[2];
   complex_dble z;
   block_t *b,*bm;

   b=blk_list(DFL_BLOCKS,&nb,&isw);
   bm=b+nb;

   x[0]=0.0;
   x[1]=0.0;

   for (;b<bm;b++)
   {
      for (i=1;i<=Ns;i++)
      {
         for (j=1;j<=i;j++)
         {
            z=spinor_prod_dble((*b).vol,0,(*b).sd[i],(*b).sd[j]);
            dev=sqrt(z.re*z.re+z.im*z.im);

            if (i==j)
               dev=fabs(1.0-dev);

            if (dev>x[0])
               x[0]=dev;
         }

         assign_s2sd((*b).vol,(*b).s[i],(*b).sd[0]);
         mulr_spinor_add_dble((*b).vol,(*b).sd[0],(*b).sd[i],-1.0);
         dev=norm_square_dble((*b).vol,0,(*b).sd[0]);
         dev=sqrt(dev);

         if (dev>x[1])
            x[1]=dev;
      }
   }

   if (NPROC>1)
   {
      MPI_Reduce(x,y,2,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(y,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

      (*dev0)=y[0];
      (*dev1)=y[1];
   }
   else
   {
      (*dev0)=x[0];
      (*dev1)=x[1];
   }
}


int main(int argc,char *argv[])
{
   int my_rank,bc,i;
   int bs[4],Ns,nv;
   double phi[2],phi_prime[2],theta[3];
   double dev,dev0,dev1;
   complex **vm,**wv,z;
   complex_dble **wvd;
   spinor **ws;
   spinor_dble **wsd;
   FILE *fin=NULL,*flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      fin=freopen("check2.in","r",stdin);

      printf("\n");
      printf("Check of the programs in the module dfl_subspace.c\n");
      printf("--------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      read_line("Ns","%d",&Ns);
      fclose(fin);

      printf("bs = %d %d %d %d\n",bs[0],bs[1],bs[2],bs[3]);
      printf("Ns = %d\n\n",Ns);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(0);

   start_ranlux(0,123456);
   geometry();
   set_dfl_parms(bs,Ns);

   alloc_ws(Ns+1);
   alloc_wsd(1);
   alloc_wv(2);
   alloc_wvd(2);

   ws=reserve_ws(Ns+1);
   wsd=reserve_wsd(1);
   vm=vflds()+Ns;
   wv=reserve_wv(2);
   wvd=reserve_wvd(2);
   nv=Ns*VOLUME/(bs[0]*bs[1]*bs[2]*bs[3]);

   for (i=0;i<Ns;i++)
   {
      random_s(VOLUME,ws[i],1.0f);
      bnd_s2zero(ALL_PTS,ws[i]);
   }

   dfl_subspace(ws);
   check_basis(Ns,&dev0,&dev1);

   message("Orthonormality of the basis vectors: %.1e\n",dev0);
   message("Single-precision basis vectors:      %.1e\n\n",dev1);

   dev0=0.0;
   dev1=0.0;

   for (i=0;i<Ns;i++)
   {
      dfl_v2s(vm[i],ws[Ns]);
      mulr_spinor_add(VOLUME,ws[Ns],ws[i],-1.0f);
      dev=(double)(norm_square(VOLUME,1,ws[Ns])/
                   norm_square(VOLUME,1,ws[i]));
      if (dev>dev0)
         dev0=dev;

      assign_s2s(VOLUME,ws[i],ws[Ns]);
      dfl_sub_v2s(vm[i],ws[Ns]);
      dev=(double)(norm_square(VOLUME,1,ws[Ns])/
                   norm_square(VOLUME,1,ws[i]));
      if (dev>dev1)
         dev1=dev;
   }

   if (my_rank==0)
   {
      printf("Check of the single-precision vector modes:\n");
      printf("Using dfl_v2s:     %.1e\n",sqrt(dev0));
      printf("Using dfl_sub_v2s: %.1e\n\n",sqrt(dev1));
   }

   dev0=0.0;
   dev1=0.0;

   random_v(nv,wv[0],1.0f);
   random_vd(nv,wvd[0],1.0);

   dfl_v2s(wv[0],ws[Ns]);
   dfl_s2v(ws[Ns],wv[1]);
   z.re=-1.0f;
   z.im=0.0f;
   mulc_vadd(nv,wv[0],wv[1],z);
   dev0=(double)(vnorm_square(nv,1,wv[0])/
                 vnorm_square(nv,1,wv[1]));

   dfl_vd2sd(wvd[0],wsd[0]);
   dfl_sd2vd(wsd[0],wvd[1]);
   diff_vd2v(nv,wvd[0],wvd[1],wv[0]);
   assign_vd2v(nv,wvd[1],wv[1]);
   dev1=(double)(vnorm_square(nv,1,wv[0])/
                 vnorm_square(nv,1,wv[1]));

   dfl_sub_vd2sd(wvd[1],wsd[0]);
   dev=norm_square_dble(VOLUME,1,wsd[0])/vnorm_square_dble(nv,1,wvd[1]);
   if (dev>dev1)
      dev1=dev;

   if (my_rank==0)
   {
      printf("Check of\n");
      printf("dfl_s2v,..:   %.1e\n",sqrt(dev0));
      printf("dfl_sd2vd,..: %.1e\n\n",sqrt(dev1));
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
