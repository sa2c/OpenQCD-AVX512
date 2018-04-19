
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2009, 2011, 2016 Filippo Palombi, Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Comparison of chexp_drv2() with chexp_drv0() and chexp_drv1() and
* invariance of the calculated coefficients under rotations of X
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "random.h"
#include "su3fcts.h"

#define NTEST 100000
#define SEED 773

static double mu[3];
static su3_alg_dble *X;
static su3_dble *r,*u,*w;
static ch_drv0_t *sp;
static ch_drv1_t *sg;
static ch_drv2_t *sf;


static void alloc_Xu(void)
{
   X=amalloc(1*sizeof(*X),4);
   r=amalloc(3*sizeof(*r),4);
   sp=amalloc(1*sizeof(*sp),4);
   sg=amalloc(1*sizeof(*sg),4);
   sf=amalloc(2*sizeof(*sf),4);

   error((X==NULL)||(r==NULL)||(sp==NULL)||(sg==NULL)||(sf==NULL),1,
         "alloc_Xu [check4.c]","Unable to allocate matrices");

   u=r+1;
   w=r+2;
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


static double dev_sp(void)
{
   int i;
   double r[8],dev,dmax;

   r[0]=(*sp).t-(*sf).t;
   r[1]=(*sp).d-(*sf).d;

   for (i=0;i<3;i++)
   {
      r[2*i+2]=(*sp).p[i].re-(*sf).p[i].re;
      r[2*i+3]=(*sp).p[i].im-(*sf).p[i].im;
   }

   dmax=0.0;

   for (i=0;i<8;i++)
   {
      dev=fabs(r[i]);
      if (dev>dmax)
         dmax=dev;
   }

   return dmax;
}


static double dev_sg(void)
{
   int i;
   double r[20],dev,dmax;

   r[0]=(*sp).t-(*sf).t;
   r[1]=(*sp).d-(*sf).d;

   for (i=0;i<3;i++)
   {
      r[6*i+2]=(*sg).p[i].re-(*sf).p[i].re;
      r[6*i+3]=(*sg).p[i].im-(*sf).p[i].im;

      r[6*i+4]=(*sg).pt[i].re-(*sf).pt[i].re;
      r[6*i+5]=(*sg).pt[i].im-(*sf).pt[i].im;

      r[6*i+6]=(*sg).pd[i].re-(*sf).pd[i].re;
      r[6*i+7]=(*sg).pd[i].im-(*sf).pd[i].im;
   }

   dmax=0.0;

   for (i=0;i<20;i++)
   {
      dev=fabs(r[i]);
      if (dev>dmax)
         dmax=dev;
   }

   return dmax;
}


static double dev_sf(void)
{
   int i;
   double r[38],dev,dmax;
   ch_drv2_t *sf1,*sf2;

   sf1=sf;
   sf2=sf+1;

   r[0]=(*sf1).t-(*sf2).t;
   r[1]=(*sf1).d-(*sf2).d;

   for (i=0;i<3;i++)
   {
      r[12*i+2]=(*sf1).p[i].re-(*sf2).p[i].re;
      r[12*i+3]=(*sf1).p[i].im-(*sf2).p[i].im;

      r[12*i+4]=(*sf1).pt[i].re-(*sf2).pt[i].re;
      r[12*i+5]=(*sf1).pt[i].im-(*sf2).pt[i].im;

      r[12*i+6]=(*sf1).pd[i].re-(*sf2).pd[i].re;
      r[12*i+7]=(*sf1).pd[i].im-(*sf2).pd[i].im;

      r[12*i+8]=(*sf1).ptt[i].re-(*sf2).ptt[i].re;
      r[12*i+9]=(*sf1).ptt[i].im-(*sf2).ptt[i].im;

      r[12*i+10]=(*sf1).ptd[i].re-(*sf2).ptd[i].re;
      r[12*i+11]=(*sf1).ptd[i].im-(*sf2).ptd[i].im;

      r[12*i+12]=(*sf1).pdd[i].re-(*sf2).pdd[i].re;
      r[12*i+13]=(*sf1).pdd[i].im-(*sf2).pdd[i].im;
   }

   dmax=0.0;

   for (i=0;i<38;i++)
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
   double dev,dmax1,dmax2,dmax3;

   printf("\n");
   printf("Invariance of chexp_drv2() under rotations of X\n");
   printf("-----------------------------------------------\n\n");

   printf("Test performed on %d random matrices X\n\n",NTEST);

   rlxd_init(1,SEED);
   alloc_Xu();

   dmax1=0.0;
   dmax2=0.0;
   dmax3=0.0;

   for (i=0;i<NTEST;i++)
   {
      random_Xu();
      chexp_drv0(X,sp);
      chexp_drv1(X,sg);
      chexp_drv2(X,sf);

      dev=dev_sp();
      if (dev>dmax1)
         dmax1=dev;

      dev=dev_sg();
      if (dev>dmax2)
         dmax2=dev;

      random_su3_dble(r);
      rotate_su3alg(r,X);
      chexp_drv2(X,sf+1);

      dev=dev_sf();
      if (dev>dmax3)
         dmax3=dev;
   }

   printf ("Comparision of chexp_drv0 and chexp_drv2 = %.1e\n",dmax1);
   printf ("Comparision of chexp_drv1 and chexp_drv2 = %.1e\n",dmax2);
   printf ("Rotation invariance of chexp_drv2        = %.1e\n\n",dmax3);

   exit(0);
}
