
/*******************************************************************************
*
* File check7.c
*
* Copyright (C) 2009, 2011, 2016 Martin Luescher, Filippo Palombi
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of expXsu3() using the spectral representation of X
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "random.h"
#include "su3fcts.h"

#define NTEST 50000
#define SEED 38579

static double mu[3],t,d,eps;
static su3_alg_dble *X;
static su3_dble *r,*u,*w,*y,*z;


static void alloc_Xu(void)
{
   X=amalloc(1*sizeof(*X),4);
   r=amalloc(5*sizeof(*r),4);

   error((X==NULL)||(r==NULL),1,
         "alloc_Xu [check7.c]","Unable to allocate matrices");

   u=r+1;
   w=r+2;
   y=r+3;
   z=r+4;
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

   t=0.5*(mu[0]*mu[0]+mu[1]*mu[1]+mu[2]*mu[2]);
   d=mu[0]*mu[1]*mu[2];

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

   ranlxd(&eps,1);
   eps*=20.0;

   cm3x3_zero(1,u);
   (*u).c11.re=cos(eps*mu[0]);
   (*u).c22.re=cos(eps*mu[1]);
   (*u).c33.re=cos(eps*mu[2]);
   (*u).c11.im=sin(eps*mu[0]);
   (*u).c22.im=sin(eps*mu[1]);
   (*u).c33.im=sin(eps*mu[2]);

   su3xsu3(r,u,w);
   su3xsu3dag(w,r,u);

   random_su3_dble(z);
   su3xsu3(u,z,w);
}


static double dev_yw(void)
{
   int i;
   double r[18],dev,dmax;

   r[ 0]=(*y).c11.re-(*w).c11.re;
   r[ 1]=(*y).c11.im-(*w).c11.im;
   r[ 2]=(*y).c12.re-(*w).c12.re;
   r[ 3]=(*y).c12.im-(*w).c12.im;
   r[ 4]=(*y).c13.re-(*w).c13.re;
   r[ 5]=(*y).c13.im-(*w).c13.im;

   r[ 6]=(*y).c21.re-(*w).c21.re;
   r[ 7]=(*y).c21.im-(*w).c21.im;
   r[ 8]=(*y).c22.re-(*w).c22.re;
   r[ 9]=(*y).c22.im-(*w).c22.im;
   r[10]=(*y).c23.re-(*w).c23.re;
   r[11]=(*y).c23.im-(*w).c23.im;

   r[12]=(*y).c31.re-(*w).c31.re;
   r[13]=(*y).c31.im-(*w).c31.im;
   r[14]=(*y).c32.re-(*w).c32.re;
   r[15]=(*y).c32.im-(*w).c32.im;
   r[16]=(*y).c33.re-(*w).c33.re;
   r[17]=(*y).c33.im-(*w).c33.im;

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
   double dev,dmax;

   printf("\n");
   printf("Check of expXsu3()\n");
   printf("------------------\n\n");

   printf("Test performed on %d random matrices X and u using the\n",NTEST);
   printf("spectral representation of X\n\n");

   rlxd_init(1,SEED);
   alloc_Xu();

   dmax=0.0;

   for (i=0;i<NTEST;i++)
   {
      random_Xu();
      (*y)=(*z);
      expXsu3(eps,X,y);

      dev=dev_yw();
      if (dev>dmax)
	dmax=dev;
   }

   printf ("Maximal deviation of exp(X)*u = %.1e\n\n",dmax);

   exit(0);
}
