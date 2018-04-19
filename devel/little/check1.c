
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2007, 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs in the module Aw_gen.c.
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
#include "utils.h"
#include "flags.h"
#include "sflds.h"
#include "linalg.h"
#include "little.h"
#include "global.h"

#define NPTS 2048

static int imb[NPTS];
static su3_dble ud[NPTS] ALIGNED16;
static su3_dble vd[NPTS] ALIGNED16;
static spinor_dble sd[3][NPTS] ALIGNED16;


static void random_imb(int vol)
{
   int i,j,a,b;
   float r[2],rvol;

   for (i=0;i<vol;i++)
      imb[i]=i;

   rvol=(float)(vol);

   for (i=0;i<(16*vol);i++)
   {
      ranlxs(r,2);

      a=(int)(rvol*r[0]);
      if (a>=vol)
         a=vol-1;
      b=(int)(rvol*r[1]);
      if (b>=vol)
         b=vol-1;

      j=imb[a];
      imb[a]=imb[b];
      imb[b]=j;
   }
}


static void random_ufld(int vol)
{
   int i;

   for (i=0;i<vol;i++)
   {
      random_su3_dble(ud+i);
      random_su3_dble(vd+i);
   }
}


static void random_sfld(int vol)
{
   int i;

   for (i=0;i<3;i++)
   {
      random_sd(vol,sd[i],1.0);
      normalize_dble(vol,0,sd[i]);
   }
}


static void permute_sd(int vol,spinor_dble *s,spinor_dble *r)
{
   int i;

   for (i=0;i<vol;i++)
      r[i]=s[imb[i]];
}


static void apply_ud(int vol,su3_dble *u,spinor_dble *s,spinor_dble *r)
{
   int i;

   for (i=0;i<vol;i++)
   {
      _su3_multiply(r[i].c1,u[i],s[i].c1);
      _su3_multiply(r[i].c2,u[i],s[i].c2);
      _su3_multiply(r[i].c3,u[i],s[i].c3);
      _su3_multiply(r[i].c4,u[i],s[i].c4);
   }
}


static su3_vector_dble mul_cplx(complex_dble z,su3_vector_dble s)
{
   su3_vector_dble r;

   r.c1.re=z.re*s.c1.re-z.im*s.c1.im;
   r.c1.im=z.im*s.c1.re+z.re*s.c1.im;
   r.c2.re=z.re*s.c2.re-z.im*s.c2.im;
   r.c2.im=z.im*s.c2.re+z.re*s.c2.im;
   r.c3.re=z.re*s.c3.re-z.im*s.c3.im;
   r.c3.im=z.im*s.c3.re+z.re*s.c3.im;

   return r;
}


static spinor_dble mul_gamma(int mu,spinor_dble s)
{
   spinor_dble r;
   complex_dble i,m_i,m_1;

   i.re=0.0;
   i.im=1.0;

   m_i.re=0.0;
   m_i.im=-1.0;

   m_1.re=-1.0;
   m_1.im=0.0;

   if (mu==0)
   {
      r.c1=mul_cplx(m_1,s.c3);
      r.c2=mul_cplx(m_1,s.c4);
      r.c3=mul_cplx(m_1,s.c1);
      r.c4=mul_cplx(m_1,s.c2);
   }
   else if (mu==1)
   {
      r.c1=mul_cplx(m_i,s.c4);
      r.c2=mul_cplx(m_i,s.c3);
      r.c3=mul_cplx(i,s.c2);
      r.c4=mul_cplx(i,s.c1);
   }
   else if (mu==2)
   {
      r.c1=mul_cplx(m_1,s.c4);
      r.c2=s.c3;
      r.c3=s.c2;
      r.c4=mul_cplx(m_1,s.c1);
   }
   else if (mu==3)
   {
      r.c1=mul_cplx(m_i,s.c3);
      r.c2=mul_cplx(i,s.c4);
      r.c3=mul_cplx(i,s.c1);
      r.c4=mul_cplx(m_i,s.c2);
   }
   else
   {
      r.c1=s.c1;
      r.c2=s.c2;
      r.c3=mul_cplx(m_1,s.c3);
      r.c4=mul_cplx(m_1,s.c4);
   }

   return r;
}


int main(int argc,char *argv[])
{
   int my_rank,ix,mu,vol,iv;
   double rv,dev0[2],dev1[2];
   complex_dble sp0[2],sp1[2];
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Check of the programs in the module Aw_gen.c\n");
      printf("--------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n",L0,L1,L2,L3);
      fflush(flog);
   }

   start_ranlux(0,123456);

   for (iv=0;iv<3;iv++)
   {
      if (iv==0)
         vol=256;
      else if (iv==1)
         vol=111;
      else
         vol=2023;

      if (my_rank==0)
      {
         printf("\n");
         printf("vol = %3d, absolute deviations:\n",vol);
      }

      rv=sqrt((double)(vol));
      random_imb(vol);
      random_ufld(vol);

      random_sfld(vol);
      gather_sd(vol,imb,sd[0],sd[1]);
      permute_sd(vol,sd[0],sd[2]);
      mulr_spinor_add_dble(vol,sd[1],sd[2],-1.0);
      dev0[0]=norm_square_dble(vol,1,sd[1]);

      random_sfld(vol);
      permute_sd(vol,sd[0],sd[1]);
      gather_ud(vol,imb,ud,vd);
      apply_ud(vol,vd,sd[1],sd[2]);
      apply_ud(vol,ud,sd[0],sd[1]);
      permute_sd(vol,sd[1],sd[0]);
      mulr_spinor_add_dble(vol,sd[0],sd[2],-1.0);
      dev0[1]=norm_square_dble(vol,1,sd[0]);

      if (my_rank==0)
      {
         printf("gather_sd():     %.1e\n",sqrt(dev0[0]));
         printf("gather_ud():     %.1e\n",sqrt(dev0[1]));
      }

      random_sfld(vol);
      apply_u2sd(vol,imb,ud,sd[0],sd[1]);
      permute_sd(vol,sd[0],sd[2]);
      apply_ud(vol,ud,sd[2],sd[0]);
      mulr_spinor_add_dble(vol,sd[1],sd[0],-1.0);
      dev0[0]=norm_square_dble(vol,1,sd[1]);

      random_sfld(vol);
      apply_udag2sd(vol,imb,ud,sd[0],sd[1]);
      permute_sd(vol,sd[0],sd[2]);
      apply_ud(vol,ud,sd[1],sd[0]);
      mulr_spinor_add_dble(vol,sd[0],sd[2],-1.0);
      dev0[1]=norm_square_dble(vol,1,sd[0]);

      if (my_rank==0)
      {
         printf("apply_u2sd():    %.1e\n",sqrt(dev0[0]));
         printf("apply_udag2sd(): %.1e\n",sqrt(dev0[1]));
      }

      for (mu=0;mu<4;mu++)
      {
         random_sd(vol,sd[0],1.0);
         random_sd(vol,sd[1],1.0);
         normalize_dble(vol,0,sd[0]);
         normalize_dble(vol,0,sd[1]);

         spinor_prod_gamma[mu](vol,sd[0],sd[1],sp0);
         sp1[0]=spinor_prod_dble(vol,0,sd[0],sd[1]);

         for (ix=0;ix<vol;ix++)
            sd[2][ix]=mul_gamma(mu,sd[1][ix]);

         sp1[1]=spinor_prod_dble(vol,0,sd[0],sd[2]);

         dev0[0]=fabs(sp0[0].re-sp1[0].re)+fabs(sp0[0].im-sp1[0].im);
         dev0[1]=fabs(sp0[1].re-sp1[1].re)+fabs(sp0[1].im-sp1[1].im);

         MPI_Reduce(dev0,dev1,2,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

         if (my_rank==0)
         {
            printf("spinor_prod_gamma[%d]: sp[0] = %.1e, sp[1] = %.1e\n",
                   mu,rv*dev1[0],rv*dev1[1]);
         }
      }
   }

   if (my_rank==0)
   {
      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
