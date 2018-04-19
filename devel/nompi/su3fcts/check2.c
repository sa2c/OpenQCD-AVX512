
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005, 2009, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of prod2su3alg, prod2u3alg and rotate_su3alg
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "utils.h"
#include "su3fcts.h"

static su3_dble Q ALIGNED16;
static su3_dble u ALIGNED16;
static su3_dble v ALIGNED16;
static su3_alg_dble X ALIGNED16;
static u3_alg_dble Y;


static void random_su3alg(su3_alg_dble *X)
{
   double r[8];

   ranlxd(r,8);

   (*X).c1=r[0]-0.5;
   (*X).c2=r[1]-0.5;
   (*X).c3=r[2]-0.5;
   (*X).c4=r[3]-0.5;
   (*X).c5=r[4]-0.5;
   (*X).c6=r[5]-0.5;
   (*X).c7=r[6]-0.5;
   (*X).c8=r[7]-0.5;
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


static void X2u(su3_alg_dble *X,su3_dble *u)
{
   (*u).c11.re=0.0;
   (*u).c11.im= (*X).c1+(*X).c2;
   (*u).c22.re=0.0;
   (*u).c22.im= (*X).c2-2.0*(*X).c1;
   (*u).c33.re=0.0;
   (*u).c33.im= (*X).c1-2.0*(*X).c2;

   (*u).c12.re= (*X).c3;
   (*u).c12.im= (*X).c4;
   (*u).c21.re=-(*X).c3;
   (*u).c21.im= (*X).c4;

   (*u).c13.re= (*X).c5;
   (*u).c13.im= (*X).c6;
   (*u).c31.re=-(*X).c5;
   (*u).c31.im= (*X).c6;

   (*u).c23.re= (*X).c7;
   (*u).c23.im= (*X).c8;
   (*u).c32.re=-(*X).c7;
   (*u).c32.im= (*X).c8;
}


int main(void)
{
   double tr,d;

   printf("\n");
   printf("Check of prod2su3alg and rotate_su3alg\n");
   printf("--------------------------------------\n\n");

#if (defined AVX)
#if (defined FMA3)
   printf("Using AVX and FMA3 instructions\n\n");
#else
   printf("Using AVX instructions\n\n");
#endif
#elif (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   rlxd_init(1,23456);

   printf("prod2su3alg:\n");
   random_su3_dble(&u);
   random_su3_dble(&v);
   random_su3alg(&X);

   tr=prod2su3alg(&u,&v,&X);
   _su3_times_su3(Q,u,v);
   tr-=(Q.c11.re+Q.c22.re+Q.c33.re);

   Q.c11.re=0.5*(Q.c11.re-Q.c11.re);
   Q.c11.im=0.5*(Q.c11.im+Q.c11.im);
   Q.c12.re=0.5*(Q.c12.re-Q.c21.re);
   Q.c12.im=0.5*(Q.c12.im+Q.c21.im);
   Q.c13.re=0.5*(Q.c13.re-Q.c31.re);
   Q.c13.im=0.5*(Q.c13.im+Q.c31.im);

   Q.c22.re=0.5*(Q.c22.re-Q.c22.re);
   Q.c22.im=0.5*(Q.c22.im+Q.c22.im);
   Q.c23.re=0.5*(Q.c23.re-Q.c32.re);
   Q.c23.im=0.5*(Q.c23.im+Q.c32.im);

   Q.c33.re=0.5*(Q.c33.re-Q.c33.re);
   Q.c33.im=0.5*(Q.c33.im+Q.c33.im);

   d=(Q.c11.im+Q.c22.im+Q.c33.im)/3.0;
   Q.c11.im-=d;
   Q.c22.im-=d;
   Q.c33.im-=d;

   d=fabs(Q.c11.im-X.c1-X.c2);
   printf("X.c11.im: %.2e\n",d);
   d=fabs(Q.c22.im+2.0*X.c1-X.c2);
   printf("X.c22.im: %.2e\n",d);
   d=fabs(Q.c33.im-X.c1+2.0*X.c2);
   printf("X.c33.im: %.2e\n",d);

   d=fabs(Q.c12.re-X.c3);
   printf("X.c12.re: %.2e\n",d);
   d=fabs(Q.c12.im-X.c4);
   printf("X.c12.im: %.2e\n",d);

   d=fabs(Q.c13.re-X.c5);
   printf("X.c13.re: %.2e\n",d);
   d=fabs(Q.c13.im-X.c6);
   printf("X.c13.im: %.2e\n",d);

   d=fabs(Q.c23.re-X.c7);
   printf("X.c23.re: %.2e\n",d);
   d=fabs(Q.c23.im-X.c8);
   printf("X.c23.im: %.2e\n",d);
   d=fabs(tr);
   printf("Return value: %.2e\n\n",d);

   printf("prod2u3alg:\n");
   random_su3_dble(&u);
   random_su3_dble(&v);
   random_u3alg(&Y);

   prod2u3alg(&u,&v,&Y);
   _su3_times_su3(Q,u,v);

   Q.c11.re=Q.c11.re-Q.c11.re;
   Q.c11.im=Q.c11.im+Q.c11.im;
   Q.c12.re=Q.c12.re-Q.c21.re;
   Q.c12.im=Q.c12.im+Q.c21.im;
   Q.c13.re=Q.c13.re-Q.c31.re;
   Q.c13.im=Q.c13.im+Q.c31.im;

   Q.c22.re=Q.c22.re-Q.c22.re;
   Q.c22.im=Q.c22.im+Q.c22.im;
   Q.c23.re=Q.c23.re-Q.c32.re;
   Q.c23.im=Q.c23.im+Q.c32.im;

   Q.c33.re=Q.c33.re-Q.c33.re;
   Q.c33.im=Q.c33.im+Q.c33.im;

   d=fabs(Q.c11.im-Y.c1);
   printf("X.c11.im: %.2e\n",d);
   d=fabs(Q.c22.im-Y.c2);
   printf("X.c22.im: %.2e\n",d);
   d=fabs(Q.c33.im-Y.c3);
   printf("X.c33.im: %.2e\n",d);

   d=fabs(Q.c12.re-Y.c4);
   printf("X.c12.re: %.2e\n",d);
   d=fabs(Q.c12.im-Y.c5);
   printf("X.c12.im: %.2e\n",d);

   d=fabs(Q.c13.re-Y.c6);
   printf("X.c13.re: %.2e\n",d);
   d=fabs(Q.c13.im-Y.c7);
   printf("X.c13.im: %.2e\n",d);

   d=fabs(Q.c23.re-Y.c8);
   printf("X.c23.re: %.2e\n",d);
   d=fabs(Q.c23.im-Y.c9);
   printf("X.c23.im: %.2e\n\n",d);

   printf("rotate_su3alg:\n");
   random_su3_dble(&u);
   random_su3alg(&X);
   X2u(&X,&v);

   rotate_su3alg(&u,&X);

   _su3_times_su3(Q,u,v);
   _su3_dagger(v,u);
   _su3_times_su3(u,Q,v);

   d=fabs(u.c11.im-X.c1-X.c2);
   printf("X.c11.im: %.2e\n",d);
   d=fabs(u.c22.im+2.0*X.c1-X.c2);
   printf("X.c22.im: %.2e\n",d);
   d=fabs(u.c33.im-X.c1+2.0*X.c2);
   printf("X.c33.im: %.2e\n",d);

   d=fabs(u.c12.re-X.c3);
   printf("X.c12.re: %.2e\n",d);
   d=fabs(u.c12.im-X.c4);
   printf("X.c12.im: %.2e\n",d);

   d=fabs(u.c13.re-X.c5);
   printf("X.c13.re: %.2e\n",d);
   d=fabs(u.c13.im-X.c6);
   printf("X.c13.im: %.2e\n",d);

   d=fabs(u.c23.re-X.c7);
   printf("X.c23.re: %.2e\n",d);
   d=fabs(u.c23.im-X.c8);
   printf("X.c23.im: %.2e\n\n",d);

   exit(0);
}
