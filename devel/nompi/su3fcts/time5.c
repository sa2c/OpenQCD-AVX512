
/*******************************************************************************
*
* File time5.c
*
* Copyright (C) 2009, 2011, 2013, 2016 Filippo Palombi, Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of chexp_drv*(), ch2mat() and expXsu3()
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "utils.h"
#include "random.h"
#include "su3fcts.h"

static double mu[3];
static su3_alg_dble *X;
static su3_dble *r,*u,*v,*w,*uu;
static ch_drv0_t *sp;
static ch_drv1_t *sg;
static ch_drv2_t *sf;
static double eps;


static void alloc_Xu(void)
{
   X=amalloc(2*sizeof(*X),4);
   r=amalloc(4*sizeof(*r),4);
   sp=amalloc(2*sizeof(*sp),4);
   sg=amalloc(2*sizeof(*sg),4);
   sf=amalloc(2*sizeof(*sf),4);
   uu=amalloc(2*sizeof(*uu),4);

   error((X==NULL)||(r==NULL)||(sp==NULL)||(sf==NULL)||(sg==NULL)||(uu==NULL),1,
         "alloc_Xu [time5.c]","Unable to allocate matrices");

   u=r+1;
   v=r+2;
   w=r+3;
}


static void random_X(void)
{
   int i;

   ranlxd(&eps,1);
   eps*=0.5;

   for (i=0;i<2;i++)
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

      X[i].c1=((*u).c11.im-(*u).c22.im)/3.0;
      X[i].c2=((*u).c11.im-(*u).c33.im)/3.0;
      X[i].c3=(*u).c12.re;
      X[i].c4=(*u).c12.im;
      X[i].c5=(*u).c13.re;
      X[i].c6=(*u).c13.im;
      X[i].c7=(*u).c23.re;
      X[i].c8=(*u).c23.im;
   }

   random_su3_dble(uu);
   random_su3_dble(uu+1);
}


static int eval_nsplt(double eps,su3_alg_dble *X)
{
   double nfrb;
   int n;

   nfrb=4.0*(3.0*((*X).c1*(*X).c1+(*X).c2*(*X).c2-(*X).c1*(*X).c2)+
             (*X).c3*(*X).c3+(*X).c4*(*X).c4+(*X).c5*(*X).c5+
             (*X).c6*(*X).c6+(*X).c7*(*X).c7+(*X).c8*(*X).c8);

   nfrb*=eps*eps;
   n=0;
   while(nfrb>3.0)
   {
      nfrb*=0.25;
      n++;
   }

   return n;
}


static int find_N(void)
{
   int i;
   double r;

   r=1.0;

   for (i=1;r>DBL_EPSILON;i++)
      r/=(double)(i);

   i+=7;

   return i+(i%2);
}


int main(void)
{
   int k,n,count,ns,nsplt,nop;
   double t1,t2,dt;

   printf("\n");
   printf("Timing of chexp_drv*(), ch2mat() and expXsu3()\n");
   printf("----------------------------------------------\n\n");

#if (defined AVX)
#if (defined FMA3)
   printf("Using AVX and FMA3 instructions\n");
#else
   printf("Using AVX instructions\n");
#endif
#elif (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n");
#endif

   printf("Measurement made with all data in cache\n\n");

   alloc_Xu();
   rlxd_init(1,12345);
   random_X();
   ns=find_N();

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
      {
         chexp_drv0(X,sp);
         chexp_drv0(X+1,sp+1);
      }
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=(1.0e6/(double)(n));

   printf("Time per call of chexp_drv0():\n");
   printf("%4.3f nsec [%d Mflops]\n\n",
          1.0e3*dt,(int)((double)(76+(ns-6)*7)/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
      {
         chexp_drv1(X,sg);
         chexp_drv1(X+1,sg+1);
      }
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=(1.0e6/(double)(n));

   printf("Time per call of chexp_drv1():\n");
   printf("%4.3f nsec [%d Mflops])\n\n",
          1.0e3*dt,(int)((double)(82+(ns-3)*15)/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
      {
         chexp_drv2(X,sf);
         chexp_drv2(X+1,sf+1);
      }
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=(1.0e6/(double)(n));

   printf("Time per call of chexp_drv2():\n");
   printf("%4.3f nsec [%d Mflops])\n\n",
          1.0e3*dt,(int)((double)(106+ns*25)/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
      {
         ch2mat(sp[0].p,X,u);
         ch2mat(sp[1].p,X+1,v);
      }
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=(1.0e6/(double)(n));

   printf("Time per call of ch2mat():\n");
   printf("%4.3f nsec [%d Mflops])\n\n",1.0e3*dt,(int)(212.0/dt));
   printf("Time per call of expXsu3():\n");

   for (k=0;k<8;k++)
   {
      eps*=2.0;
      nsplt=eval_nsplt(eps,X);

      n=(int)(1.0e6);
      dt=0.0;

      while (dt<1.5)
      {
         t1=(double)clock();
         for (count=0;count<n;count++)
         {
            expXsu3(eps,X,uu);
            expXsu3(eps,X,uu+1);
         }
         t2=(double)clock();
         dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
         n*=2;
      }

      dt*=(1.0e6/(double)(n));
      nop=8+76+(ns-6)*7+212+192*(nsplt+1);

      printf("no. splits: %2d,  %4.3f nsec [%d Mflops]\n",
             nsplt,1.0e3*dt,(int)(nop/dt));
   }

   printf("\n");
   exit(0);
}
