
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2007, 2009, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of cmat_vec_dble, cmat_add_dble, ...
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "random.h"
#include "utils.h"
#include "linalg.h"

#define NMAX 32

#define cadd(u,v,w) \
   (u).re=(v).re+(w).re;\
   (u).im=(v).im+(w).im

#define csub(u,v,w) \
   (u).re=(v).re-(w).re;\
   (u).im=(v).im-(w).im

#define cmul(u,v,w) \
   (u).re=(v).re*(w).re-(v).im*(w).im;\
   (u).im=(v).re*(w).im+(v).im*(w).re

#define cmul_assign(u,v,w) \
   (u).re+=((v).re*(w).re-(v).im*(w).im);\
   (u).im+=((v).re*(w).im+(v).im*(w).re)


static void mvec(int n,complex_dble *a,complex_dble *v,complex_dble *w)
{
   int i,j;
   complex_dble z;

   for (i=0;i<n;i++)
   {
      z.re=0.0;
      z.im=0.0;

       for (j=0;j<n;j++)
       {
          cmul_assign(z,a[i*n+j],v[j]);
       }

       w[i]=z;
   }
}


static void mvec_assign(int n,complex_dble *a,complex_dble *v,complex_dble *w)
{
   int i,j;
   complex_dble z;

   for (i=0;i<n;i++)
   {
      z.re=0.0;
      z.im=0.0;

       for (j=0;j<n;j++)
       {
          cmul_assign(z,a[i*n+j],v[j]);
       }

       w[i].re+=z.re;
       w[i].im+=z.im;
   }
}


static void madd(int n,complex_dble *a,complex_dble *b,complex_dble *c)
{
   int i,j;

   for (i=0;i<n;i++)
   {
       for (j=0;j<n;j++)
       {
          cadd(c[i*n+j],a[i*n+j],b[i*n+j]);
       }
   }
}


static void msub(int n,complex_dble *a,complex_dble *b,complex_dble *c)
{
   int i,j;

   for (i=0;i<n;i++)
   {
       for (j=0;j<n;j++)
       {
          csub(c[i*n+j],a[i*n+j],b[i*n+j]);
       }
   }
}


static void mdag(int n,complex_dble *a,complex_dble *b)
{
   int i,j;

   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
      {
         b[i*n+j].re=a[j*n+i].re;
         b[i*n+j].im=-a[j*n+i].im;
      }
   }
}


static void mmul(int n,complex_dble *a,complex_dble *b,complex_dble *c)
{
   int i,j,l;
   complex_dble z;

   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
      {
         z.re=0.0;
         z.im=0.0;

         for (l=0;l<n;l++)
         {
            cmul_assign(z,a[i*n+l],b[l*n+j]);
         }

         c[i*n+j]=z;

      }
   }
}


static double fnorm(int n,complex_dble *a)
{
   int i,j;
   double sm;

   sm=0.0;

   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
         sm+=(a[i*n+j].re*a[i*n+j].re+a[i*n+j].im*a[i*n+j].im);
   }

   return sqrt(sm);
}


static double vdev(int n,complex_dble *v,complex_dble *w)
{
   int i;
   double d,dmax;

   dmax=0.0;

   for (i=0;i<n;i++)
   {
      d=fabs(v[i].re-w[i].re)+fabs(v[i].im-w[i].im);
      if (d>dmax)
         dmax=d;
   }

   return dmax;
}


static double mdev(int n,complex_dble *a,complex_dble *b)
{
   int i,j;
   double d,dmax;

   dmax=0.0;

   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
      {
         d=fabs(a[i*n+j].re-b[i*n+j].re)+fabs(a[i*n+j].im-b[i*n+j].im);
         if (d>dmax)
            dmax=d;
      }
   }

   return dmax;
}


static void rvec(int n,complex_dble *v)
{
   gauss_dble((double*)(v),2*n);
}


static void rmat(int n,complex_dble *a)
{
   int i,j;
   double r;

   r=0.1/(double)(n*n);

   gauss_dble((double*)(a),2*n*n);

   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
      {
         a[i*n+j].re*=r;
         a[i*n+j].im*=r;
      }

      a[i*n+i].re+=1.0;
   }
}


static void vec2vec(int n,complex_dble *v,complex_dble *w)
{
   int i;

   for (i=0;i<n;i++)
      w[i]=v[i];
}


static void mat2mat(int n,complex_dble *a,complex_dble *b)
{
   int i,j;

   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
         b[i*n+j]=a[i*n+j];
   }
}


int main(void)
{
   int n,ie;
   double d,d1,d2,d3,d4,d5,d6,d7;
   double k1,k2,kmax;
   complex_dble *a1,*b1,*c1,*v1,*w1;
   complex_dble *a2,*b2,*c2,*v2,*w2;

   printf("\n");
   printf("Check of cmat_vec_dble, cmat_add_dble, ...\n");
   printf("------------------------------------------\n\n");

#if (defined AVX)
#if (defined FMA3)
   printf("Using AVX and FMA3 instructions\n\n");
#else
   printf("Using AVX instructions\n\n");
#endif
#elif (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   a1=amalloc(2*NMAX*(3*NMAX+2)*sizeof(*a1),6);
   error(a1==NULL,1,"main [check2.c]","Unable to allocate auxiliary arrays");

   b1=a1+NMAX*NMAX;
   c1=b1+NMAX*NMAX;
   v1=c1+NMAX*NMAX;
   w1=v1+NMAX;
   a2=w1+NMAX;
   b2=a2+NMAX*NMAX;
   c2=b2+NMAX*NMAX;
   v2=c2+NMAX*NMAX;
   w2=v2+NMAX;

   d1=0.0;
   d2=0.0;
   d3=0.0;
   d4=0.0;
   d5=0.0;
   d6=0.0;
   d7=0.0;
   kmax=0.0;

   for (n=1;n<=NMAX;n++)
   {
      rvec(n,v1);
      rmat(n,a1);
      vec2vec(n,v1,v2);
      mat2mat(n,a1,a2);

      cmat_vec_dble(n,a1,v1,w1);
      mvec(n,a2,v2,w2);

      d=vdev(n,w1,w2);
      if (d>d1)
         d1=d;

      error((mdev(n,a1,a2)!=0.0)||(vdev(n,v1,v2)!=0.0),1,"main [check2.c]",
            "cmat_vec_dble: input values have changed");

      rvec(n,v1);
      rvec(n,w1);
      rmat(n,a1);
      vec2vec(n,v1,v2);
      vec2vec(n,w1,w2);
      mat2mat(n,a1,a2);

      cmat_vec_assign_dble(n,a1,v1,w1);
      mvec_assign(n,a2,v2,w2);

      d=vdev(n,w1,w2);
      if (d>d1)
         d1=d;

      error((mdev(n,a1,a2)!=0.0)||(vdev(n,v1,v2)!=0.0),1,"main [check2.c]",
            "cmat_vec_assign_dble: input values have changed");

      rmat(n,a1);
      rmat(n,b1);
      rmat(n,c1);
      mat2mat(n,a1,a2);
      mat2mat(n,b1,b2);
      rmat(n,c2);

      cmat_add_dble(n,a1,b1,c1);
      madd(n,a2,b2,c2);

      d=mdev(n,c1,c2);
      if (d>d2)
         d2=d;

      error((mdev(n,a1,a2)!=0.0)||(mdev(n,b1,b2)!=0.0),1,"main [check2.c]",
            "cmat_add_dble: input values have changed");

      rmat(n,a1);
      rmat(n,b1);
      rmat(n,c1);
      mat2mat(n,a1,a2);
      mat2mat(n,b1,b2);
      rmat(n,c2);

      cmat_sub_dble(n,a1,b1,c1);
      msub(n,a2,b2,c2);

      d=mdev(n,c1,c2);
      if (d>d3)
         d3=d;

      error((mdev(n,a1,a2)!=0.0)||(mdev(n,b1,b2)!=0.0),1,"main [check2.c]",
            "cmat_sub_dble: input values have changed");

      rmat(n,a1);
      rmat(n,b1);
      rmat(n,c1);
      mat2mat(n,a1,a2);
      mat2mat(n,b1,b2);
      rmat(n,c2);

      cmat_mul_dble(n,a1,b1,c1);
      mmul(n,a2,b2,c2);

      d=mdev(n,c1,c2);
      if (d>d4)
         d4=d;

      error((mdev(n,a1,a2)!=0.0)||(mdev(n,b1,b2)!=0.0),1,"main [check2.c]",
            "cmat_mul_dble: input values have changed");

      rmat(n,a1);
      rmat(n,b1);
      mat2mat(n,a1,a2);
      rmat(n,b2);

      cmat_dag_dble(n,a1,b1);
      mdag(n,a2,b2);

      d=mdev(n,b1,b2);
      if (d>d5)
         d5=d;

      error(mdev(n,a1,a2)!=0.0,1,"main [check2.c]",
            "cmat_dag_dble: input values have changed");

      rmat(n,a1);
      rmat(n,b1);
      mat2mat(n,a1,a2);
      rmat(n,b2);
      rmat(n,c2);

      ie=cmat_inv_dble(n,a1,b1,&k1);
      mmul(n,a2,b1,b2);
      mmul(n,a2,b2,c2);

      d=mdev(n,a2,c2);
      if (d>d6)
         d6=d;

      if (k1>kmax)
         kmax=k1;

      k2=fnorm(n,a1)*fnorm(n,b1);
      d=fabs(k2/k1-1.0);
      if (d>d7)
         d7=d;

      error(ie!=0,1,"main [check2.c]",
            "cmat_inv_dble: singular matrix encountered");

      error(mdev(n,a1,a2)!=0.0,1,"main [check2.c]",
            "cmat_inv_dble: input values have changed");
   }

   printf("Consider matrices of size up to %dx%d\n\n",NMAX,NMAX);

   printf("The maximal observed deviations are:\n\n");
   printf("cmat_vec_dble:     %.1e\n",d1);
   printf("cmat_add_dble:     %.1e\n",d2);
   printf("cmat_sub_dble:     %.1e\n",d3);
   printf("cmat_mul_dble:     %.1e\n",d4);
   printf("cmat_dag_dble:     %.1e\n",d5);
   printf("cmat_inv_dble:     %.1e,  condition number:  max=%.1e, dev=%.1e\n\n",
          d6,kmax,d7);

   exit(0);
}
