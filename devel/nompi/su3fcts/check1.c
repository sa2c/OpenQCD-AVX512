
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005, 2009, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of su3xsu3, su3dagxsu3, ...
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "random.h"
#include "utils.h"
#include "su3fcts.h"


static double max_dev(su3_dble *u,su3_dble *v)
{
   int i;
   double r[18],s[18];
   double nrm,d,dmax;

   r[ 0]=(*u).c11.re;
   r[ 1]=(*u).c11.im;
   r[ 2]=(*u).c12.re;
   r[ 3]=(*u).c12.im;
   r[ 4]=(*u).c13.re;
   r[ 5]=(*u).c13.im;

   r[ 6]=(*u).c21.re;
   r[ 7]=(*u).c21.im;
   r[ 8]=(*u).c22.re;
   r[ 9]=(*u).c22.im;
   r[10]=(*u).c23.re;
   r[11]=(*u).c23.im;

   r[12]=(*u).c31.re;
   r[13]=(*u).c31.im;
   r[14]=(*u).c32.re;
   r[15]=(*u).c32.im;
   r[16]=(*u).c33.re;
   r[17]=(*u).c33.im;

   s[ 0]=(*v).c11.re;
   s[ 1]=(*v).c11.im;
   s[ 2]=(*v).c12.re;
   s[ 3]=(*v).c12.im;
   s[ 4]=(*v).c13.re;
   s[ 5]=(*v).c13.im;

   s[ 6]=(*v).c21.re;
   s[ 7]=(*v).c21.im;
   s[ 8]=(*v).c22.re;
   s[ 9]=(*v).c22.im;
   s[10]=(*v).c23.re;
   s[11]=(*v).c23.im;

   s[12]=(*v).c31.re;
   s[13]=(*v).c31.im;
   s[14]=(*v).c32.re;
   s[15]=(*v).c32.im;
   s[16]=(*v).c33.re;
   s[17]=(*v).c33.im;

   nrm=0.0;
   dmax=0.0;

   for (i=0;i<18;i++)
   {
      nrm+=r[i]*r[i];
      d=(r[i]-s[i])*(r[i]-s[i]);
      i+=1;
      nrm+=r[i]*r[i];
      d+=(r[i]-s[i])*(r[i]-s[i]);

      if (d>dmax)
         dmax=d;
   }

   return sqrt(dmax/nrm);
}


static void random_u3alg(u3_alg_dble *X)
{
   double r[9];

   ranlxd(r,9);

   (*X).c1=r[0]-0.5;
   (*X).c2=r[1]-0.5;
   (*X).c3=r[2]-0.5;
   (*X).c4=r[3]-0.5;
   (*X).c5=r[4]-0.5;
   (*X).c6=r[5]-0.5;
   (*X).c7=r[6]-0.5;
   (*X).c8=r[7]-0.5;
   (*X).c9=r[8]-0.5;
}


static void X2u(u3_alg_dble *X,su3_dble *u)
{
   (*u).c11.re=0.0;
   (*u).c11.im= (*X).c1;
   (*u).c22.re=0.0;
   (*u).c22.im= (*X).c2;
   (*u).c33.re=0.0;
   (*u).c33.im= (*X).c3;

   (*u).c12.re= (*X).c4;
   (*u).c12.im= (*X).c5;
   (*u).c21.re=-(*X).c4;
   (*u).c21.im= (*X).c5;

   (*u).c13.re= (*X).c6;
   (*u).c13.im= (*X).c7;
   (*u).c31.re=-(*X).c6;
   (*u).c31.im= (*X).c7;

   (*u).c23.re= (*X).c8;
   (*u).c23.im= (*X).c9;
   (*u).c32.re=-(*X).c8;
   (*u).c32.im= (*X).c9;
}


int main(void)
{
   double d1,d2,d3,d4;
   su3_dble *u,*v,*w1,*w2;
   u3_alg_dble *X;

   printf("\n");
   printf("Check of su3xsu3, su3dagxsu3, ...\n");
   printf("---------------------------------\n\n");

#if (defined AVX)
#if (defined FMA3)
   printf("Using AVX and FMA3 instructions\n\n");
#else
   printf("Using AVX instructions\n\n");
#endif
#elif (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   u=amalloc(4*sizeof(su3_dble),4);
   X=amalloc(sizeof(u3_alg_dble),3);
   error((u==NULL)||(X==NULL),1,"main [check1.c]",
         "Unable to allocate auxiliary arrays");

   v=u+1;
   w1=u+2;
   w2=u+3;

   rlxd_init(1,23456);

   random_su3_dble(u);
   random_su3_dble(v);
   su3xsu3(u,v,w1);
   _su3_times_su3(*w2,*u,*v);
   d1=max_dev(w1,w2);

   random_su3_dble(u);
   random_su3_dble(v);
   su3dagxsu3(u,v,w1);
   _su3_dagger(*w2,*u);
   *u=*w2;
   _su3_times_su3(*w2,*u,*v);
   d2=max_dev(w1,w2);

   random_su3_dble(u);
   random_su3_dble(v);
   su3xsu3dag(u,v,w1);
   _su3_dagger(*w2,*v);
   *v=*w2;
   _su3_times_su3(*w2,*u,*v);
   d3=max_dev(w1,w2);

   random_su3_dble(u);
   random_su3_dble(v);
   su3dagxsu3dag(u,v,w1);
   _su3_dagger(*w2,*u);
   *u=*w2;
   _su3_dagger(*w2,*v);
   *v=*w2;
   _su3_times_su3(*w2,*u,*v);
   d4=max_dev(w1,w2);

   printf("su3xsu3:       %.2e\n",d1);
   printf("su3dagxsu3:    %.2e\n",d2);
   printf("su3xsu3dag:    %.2e\n",d3);
   printf("su3dagxsu3dag: %.2e\n",d4);

   random_su3_dble(u);
   random_u3alg(X);
   su3xu3alg(u,X,w1);
   X2u(X,v);
   _su3_times_su3(*w2,*u,*v);
   d1=max_dev(w1,w2);

   random_su3_dble(u);
   random_u3alg(X);
   su3dagxu3alg(u,X,w1);
   _su3_dagger(*w2,*u);
   *u=*w2;
   X2u(X,v);
   _su3_times_su3(*w2,*u,*v);
   d2=max_dev(w1,w2);

   random_su3_dble(v);
   random_u3alg(X);
   u3algxsu3(X,v,w1);
   X2u(X,u);
   _su3_times_su3(*w2,*u,*v);
   d3=max_dev(w1,w2);

   random_su3_dble(v);
   random_u3alg(X);
   u3algxsu3dag(X,v,w1);
   X2u(X,u);
   _su3_dagger(*w2,*v);
   *v=*w2;
   _su3_times_su3(*w2,*u,*v);
   d4=max_dev(w1,w2);

   printf("su3xu3alg:     %.2e\n",d1);
   printf("su3dagxu3alg:  %.2e\n",d2);
   printf("u3algxsu3:     %.2e\n",d3);
   printf("u3algxsu3dag:  %.2e\n\n",d4);

   exit(0);
}
