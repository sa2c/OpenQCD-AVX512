
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2009, 2011, 2016 Filippo Palombi, Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Verifies that exp(X) is in SU(3) and that exp(X)*exp(-X)=1
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "random.h"
#include "su3fcts.h"

#define NTEST 100000
#define SEED 1234

static double mu[3];
static su3_alg_dble *X;
static su3_dble *r,*u,*v,*w;
static ch_drv0_t *sp;


static void alloc_Xu(void)
{
   X=amalloc(1*sizeof(*X),4);
   r=amalloc(4*sizeof(*r),4);
   sp=amalloc(1*sizeof(*sp),4);

   error((X==NULL)||(r==NULL)||(sp==NULL),1,
         "alloc_Xu [check6.c]","Unable to allocate matrices");

   u=r+1;
   v=r+2;
   w=r+3;
}


static void random_Xu(void)
{
   for (;;)
   {
      ranlxd(mu,2);
      mu[0]=2.0*mu[0]-1.0;
      mu[1]=2.0*mu[1]-1.0;
      mu[2]=-mu[0]-mu[1];

      if (fabs(mu[2])<=1.0)
         break;
   }

   cm3x3_zero(1,u);
   (*u).c11.im=mu[0];
   (*u).c22.im=mu[1];
   (*u).c33.im=mu[2];

   random_su3_dble(r);
   su3xsu3(r,u,w);
   su3xsu3dag(w,r,u);

   (*X).c1=((*u).c11.im-(*u).c22.im)/3.0;
   (*X).c2=((*u).c11.im-(*u).c33.im)/3.0;
   (*X).c3=(*u).c12.re;
   (*X).c4=(*u).c12.im;
   (*X).c5=(*u).c13.re;
   (*X).c6=(*u).c13.im;
   (*X).c7=(*u).c23.re;
   (*X).c8=(*u).c23.im;
}


static void flip_signX(void)
{
   (*X).c1=-(*X).c1;
   (*X).c2=-(*X).c2;
   (*X).c3=-(*X).c3;
   (*X).c4=-(*X).c4;
   (*X).c5=-(*X).c5;
   (*X).c6=-(*X).c6;
   (*X).c7=-(*X).c7;
   (*X).c8=-(*X).c8;
}


static double dev_uv(void)
{
   int i;
   double r[18],dev,dmax;

   r[ 0]=(*u).c11.re-(*v).c11.re;
   r[ 1]=(*u).c11.im-(*v).c11.im;
   r[ 2]=(*u).c12.re-(*v).c12.re;
   r[ 3]=(*u).c12.im-(*v).c12.im;
   r[ 4]=(*u).c13.re-(*v).c13.re;
   r[ 5]=(*u).c13.im-(*v).c13.im;

   r[ 6]=(*u).c21.re-(*v).c21.re;
   r[ 7]=(*u).c21.im-(*v).c21.im;
   r[ 8]=(*u).c22.re-(*v).c22.re;
   r[ 9]=(*u).c22.im-(*v).c22.im;
   r[10]=(*u).c23.re-(*v).c23.re;
   r[11]=(*u).c23.im-(*v).c23.im;

   r[12]=(*u).c31.re-(*v).c31.re;
   r[13]=(*u).c31.im-(*v).c31.im;
   r[14]=(*u).c32.re-(*v).c32.re;
   r[15]=(*u).c32.im-(*v).c32.im;
   r[16]=(*u).c33.re-(*v).c33.re;
   r[17]=(*u).c33.im-(*v).c33.im;

   dmax=0.0;

   for (i=0;i<18;i++)
   {
      dev=fabs(r[i]);
      if (dev>dmax)
         dmax=dev;
   }

   return dmax;
}


int main(void)
{
   int i;
   double dev,dmax1,dmax2;

   printf("\n");
   printf("Simple checks of exp(X) as calculated by chexp_drv0()\n");
   printf("-----------------------------------------------------\n\n");

   printf("Test performed on %d random matrices X\n",NTEST);

   rlxd_init(1,SEED);
   alloc_Xu();

   dmax1=0.0;
   dmax2=0.0;

   for (i=0;i<NTEST;i++)
   {
      random_Xu();
      chexp_drv0(X,sp);
      ch2mat((*sp).p,X,u);
      (*v)=(*u);
      project_to_su3_dble(v);

      dev=dev_uv();
      if (dev>dmax1)
         dmax1=dev;

      flip_signX();
      chexp_drv0(X,sp);
      ch2mat((*sp).p,X,w);
      su3xsu3(u,w,v);
      cm3x3_unity(1,u);

      dev=dev_uv();
      if (dev>dmax2)
         dmax2=dev;
   }

   printf ("Maximal deviation of exp(X) from SU(3)     = %.1e\n",dmax1);
   printf ("Maximal deviation of exp(X)*exp(-X) from 1 = %.1e\n\n",dmax2);

   exit(0);
}
