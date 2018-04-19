
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2009, 2011, 2016 Filippo Palombi, Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of chexp_drv2() in the case of diagonal X
*
* This program verifies that eqs. (4.1)-(4.6) of the notes "SU(3) matrix
* functions" are satisfied by the coefficients obtained by chexp_drv2().
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "random.h"
#include "su3fcts.h"

#define NTEST 100000
#define SEED 8923

static double mu[3],xt[2],yt[2];
static complex_dble t[2][3],x[3],y[3];
static complex_dble xsq[3],stx[2][3],stt[2][2][3];
static complex_dble ex[3],dex[2][3],ddex[2][2][3];
static complex_dble df[2][3],ddf[2][2][3];
static su3_alg_dble *X;
static ch_drv2_t *sf;


static void mul_vec(complex_dble *a,complex_dble *b,complex_dble *c)
{
   int i;

   for (i=0;i<3;i++)
   {
      c[i].re=a[i].re*b[i].re-a[i].im*b[i].im;
      c[i].im=a[i].re*b[i].im+a[i].im*b[i].re;
   }
}


static void add_vec(complex_dble z,complex_dble *a,complex_dble *b)
{
   int i;

   for (i=0;i<3;i++)
   {
      b[i].re+=(z.re*a[i].re-z.im*a[i].im);
      b[i].im+=(z.re*a[i].im+z.im*a[i].re);
   }
}


static void alloc_X(void)
{
   X=amalloc(1*sizeof(*X),4);
   sf=amalloc(1*sizeof(*sf),4);

   error((X==NULL)||(sf==NULL),1,
         "alloc_X [check5.c]","Unable to allocate matrices");
}


static void set_tk(void)
{
   double r;

   t[0][0].re=0.0;
   t[0][0].im=0.5;
   t[0][1].re=0.0;
   t[0][1].im=-0.5;
   t[0][2].re=0.0;
   t[0][2].im=0.0;

   r=1.0/(2.0*sqrt(3.0));

   t[1][0].re=0.0;
   t[1][0].im=r;
   t[1][1].re=0.0;
   t[1][1].im=r;
   t[1][2].re=0.0;
   t[1][2].im=-2.0*r;
}


static void random_X(void)
{
   double s;

   for (;;)
   {
      ranlxd(mu,2);
      mu[0]=2.0*mu[0]-1.0;
      mu[1]=2.0*mu[1]-1.0;
      mu[2]=-mu[0]-mu[1];

      if (fabs(mu[2])<=1.0)
         break;
   }

   (*X).c1=(mu[0]-mu[1])/3.0;
   (*X).c2=(mu[0]-mu[2])/3.0;
   (*X).c3=0.0;
   (*X).c4=0.0;
   (*X).c5=0.0;
   (*X).c6=0.0;
   (*X).c7=0.0;
   (*X).c8=0.0;

   x[0].re=0.0;
   x[0].im=mu[0];
   x[1].re=0.0;
   x[1].im=mu[1];
   x[2].re=0.0;
   x[2].im=mu[2];

   xt[0]=x[0].im-x[1].im;
   xt[1]=sqrt(3.0)*(x[0].im+x[1].im);

   s=(mu[0]*mu[0]+mu[1]*mu[1]+mu[2]*mu[2])/3.0;

   y[0].re=0.0;
   y[0].im=mu[0]*mu[0]-s;
   y[1].re=0.0;
   y[1].im=mu[1]*mu[1]-s;
   y[2].re=0.0;
   y[2].im=mu[2]*mu[2]-s;

   yt[0]=y[0].im-y[1].im;
   yt[1]=sqrt(3.0)*(y[0].im+y[1].im);

   ex[0].re=cos(mu[0]);
   ex[0].im=sin(mu[0]);
   ex[1].re=cos(mu[1]);
   ex[1].im=sin(mu[1]);
   ex[2].re=cos(mu[2]);
   ex[2].im=sin(mu[2]);
}


static void diff_exp(void)
{
   int i,j;

   for (i=0;i<2;i++)
      mul_vec(t[i],ex,dex[i]);

   for (i=0;i<2;i++)
   {
      for (j=0;j<2;j++)
         mul_vec(t[i],dex[j],ddex[i][j]);
   }
}


static void diff_fk(void)
{
   int i,j,k;
   double d;

   for (i=0;i<2;i++)
   {
      for (k=0;k<3;k++)
      {
         df[i][k].re=0.5*(xt[i]*(*sf).pt[k].re+yt[i]*(*sf).pd[k].re);
         df[i][k].im=0.5*(xt[i]*(*sf).pt[k].im+yt[i]*(*sf).pd[k].im);
      }
   }

   d=1.0/sqrt(3.0);

   for (k=0;k<3;k++)
   {
      ddf[0][0][k].re=0.5*((*sf).pt[k].re+d*xt[1]*(*sf).pd[k].re);
      ddf[0][0][k].im=0.5*((*sf).pt[k].im+d*xt[1]*(*sf).pd[k].im);

      ddf[1][1][k].re=0.5*((*sf).pt[k].re-d*xt[1]*(*sf).pd[k].re);
      ddf[1][1][k].im=0.5*((*sf).pt[k].im-d*xt[1]*(*sf).pd[k].im);

      ddf[0][1][k].re=0.5*d*xt[0]*(*sf).pd[k].re;
      ddf[0][1][k].im=0.5*d*xt[0]*(*sf).pd[k].im;

      ddf[1][0][k].re=0.5*d*xt[0]*(*sf).pd[k].re;
      ddf[1][0][k].im=0.5*d*xt[0]*(*sf).pd[k].im;
   }

   for (i=0;i<2;i++)
   {
      for (j=0;j<2;j++)
      {
         for (k=0;k<3;k++)
         {
            ddf[i][j][k].re+=0.25*(xt[i]*xt[j]*(*sf).ptt[k].re+
                                   xt[i]*yt[j]*(*sf).ptd[k].re+
                                   yt[i]*xt[j]*(*sf).ptd[k].re+
                                   yt[i]*yt[j]*(*sf).pdd[k].re);

            ddf[i][j][k].im+=0.25*(xt[i]*xt[j]*(*sf).ptt[k].im+
                                   xt[i]*yt[j]*(*sf).ptd[k].im+
                                   yt[i]*xt[j]*(*sf).ptd[k].im+
                                   yt[i]*yt[j]*(*sf).pdd[k].im);
         }
      }
   }
}


static void set_prods(void)
{
   int i,j;

   mul_vec(x,x,xsq);

   for (i=0;i<2;i++)
   {
      mul_vec(t[i],x,stx[i]);

      for (j=0;j<2;j++)
         mul_vec(t[i],t[j],stt[i][j]);
   }
}


static void subtract_chexp(void)
{
   int i,j,k;
   complex_dble z;

   for (i=0;i<2;i++)
   {
      for (k=0;k<3;k++)
      {
         dex[i][k].re-=df[i][0].re;
         dex[i][k].im-=df[i][0].im;
      }

      z.re=-df[i][1].re;
      z.im=-df[i][1].im;
      add_vec(z,x,dex[i]);

      z.re=-df[i][2].re;
      z.im=-df[i][2].im;
      add_vec(z,xsq,dex[i]);

      z.re=-(*sf).p[1].re;
      z.im=-(*sf).p[1].im;
      add_vec(z,t[i],dex[i]);

      z.re=-2.0*(*sf).p[2].re;
      z.im=-2.0*(*sf).p[2].im;
      add_vec(z,stx[i],dex[i]);

      for (j=0;j<2;j++)
      {
         for (k=0;k<3;k++)
         {
            ddex[i][j][k].re-=ddf[i][j][0].re;
            ddex[i][j][k].im-=ddf[i][j][0].im;
         }

         z.re=-ddf[i][j][1].re;
         z.im=-ddf[i][j][1].im;
         add_vec(z,x,ddex[i][j]);

         z.re=-ddf[i][j][2].re;
         z.im=-ddf[i][j][2].im;
         add_vec(z,xsq,ddex[i][j]);

         z.re=-df[i][1].re;
         z.im=-df[i][1].im;
         add_vec(z,t[j],ddex[i][j]);

         z.re=-df[j][1].re;
         z.im=-df[j][1].im;
         add_vec(z,t[i],ddex[i][j]);

         z.re=-2.0*df[i][2].re;
         z.im=-2.0*df[i][2].im;
         add_vec(z,stx[j],ddex[i][j]);

         z.re=-2.0*df[j][2].re;
         z.im=-2.0*df[j][2].im;
         add_vec(z,stx[i],ddex[i][j]);

         z.re=-2.0*(*sf).p[2].re;
         z.im=-2.0*(*sf).p[2].im;
         add_vec(z,stt[i][j],ddex[i][j]);
      }
   }
}


static double dev_dex(void)
{
   int i,k;
   double dev,dmax;

   dmax=0.0;

   for (i=0;i<2;i++)
   {
      for (k=0;k<3;k++)
      {
         dev=dex[i][k].re*dex[i][k].re+dex[i][k].im*dex[i][k].im;
         if (dev>dmax)
            dmax=dev;
      }
   }

   return sqrt(dmax);
}


static double dev_ddex(void)
{
   int i,j,k;
   double dev,dmax;

   dmax=0.0;

   for (i=0;i<2;i++)
   {
      for (j=0;j<2;j++)
      {
         for (k=0;k<3;k++)
         {
            dev=ddex[i][j][k].re*ddex[i][j][k].re+
                ddex[i][j][k].im*ddex[i][j][k].im;
            if (dev>dmax)
               dmax=dev;
         }
      }
   }

   return sqrt(dmax);
}


int main(void)
{
   int i;
   double dev,dmax1,dmax2;

   printf("\n");
   printf("Check of chexp_drv2() for diagonal X\n");
   printf("------------------------------------\n\n");

   printf("Test performed on %d random matrices X\n\n",NTEST);

   rlxd_init(1,SEED);
   alloc_X();
   set_tk();

   dmax1=0.0;
   dmax2=0.0;

   for (i=0;i<NTEST;i++)
   {
      random_X();
      chexp_drv2(X,sf);
      diff_exp();
      diff_fk();
      set_prods();
      subtract_chexp();

      dev=dev_dex();
      if (dev>dmax1)
         dmax1=dev;

      dev=dev_ddex();
      if (dev>dmax2)
         dmax2=dev;
   }

   printf ("Maximal deviation of 1st derivatives = %.1e\n",dmax1);
   printf ("Maximal deviation of 2nd derivatives = %.1e\n\n",dmax2);

   exit(0);
}
