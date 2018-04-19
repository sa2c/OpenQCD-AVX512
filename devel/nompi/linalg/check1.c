
/*******************************************************************************
*
* File check1.c
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


static void mvec(int n,complex *a,complex *v,complex *w)
{
   int i,j;
   complex z;

   for (i=0;i<n;i++)
   {
      z.re=0.0f;
      z.im=0.0f;

       for (j=0;j<n;j++)
       {
          cmul_assign(z,a[i*n+j],v[j]);
       }

       w[i]=z;
   }
}


static void mvec_assign(int n,complex *a,complex *v,complex *w)
{
   int i,j;
   complex z;

   for (i=0;i<n;i++)
   {
      z.re=0.0f;
      z.im=0.0f;

       for (j=0;j<n;j++)
       {
          cmul_assign(z,a[i*n+j],v[j]);
       }

       w[i].re+=z.re;
       w[i].im+=z.im;
   }
}


static void madd(int n,complex *a,complex *b,complex *c)
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


static void msub(int n,complex *a,complex *b,complex *c)
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


static void mdag(int n,complex *a,complex *b)
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


static void mmul(int n,complex *a,complex *b,complex *c)
{
   int i,j,l;
   complex z;

   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
      {
         z.re=0.0f;
         z.im=0.0f;

         for (l=0;l<n;l++)
         {
            cmul_assign(z,a[i*n+l],b[l*n+j]);
         }

         c[i*n+j]=z;

      }
   }
}


static float vdev(int n,complex *v,complex *w)
{
   int i;
   float d,dmax;

   dmax=0.0f;

   for (i=0;i<n;i++)
   {
      d= (float)(fabs((double)(v[i].re-w[i].re)))+
         (float)(fabs((double)(v[i].im-w[i].im)));
      if (d>dmax)
         dmax=d;
   }

   return dmax;
}


static float mdev(int n,complex *a,complex *b)
{
   int i,j;
   float d,dmax;

   dmax=0.0f;

   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
      {
         d= (float)(fabs((double)(a[i*n+j].re-b[i*n+j].re)))+
            (float)(fabs((double)(a[i*n+j].im-b[i*n+j].im)));
         if (d>dmax)
            dmax=d;
      }
   }

   return dmax;
}


static void rvec(int n,complex *v)
{
   gauss((float*)(v),2*n);
}


static void rmat(int n,complex *a)
{
   int i,j;
   float r;

   r=0.1f/(float)(n*n);

   gauss((float*)(a),2*n*n);

   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
      {
         a[i*n+j].re*=r;
         a[i*n+j].im*=r;
      }

      a[i*n+i].re+=1.0f;
   }
}


static void vec2vec(int n,complex *v,complex *w)
{
   int i;

   for (i=0;i<n;i++)
      w[i]=v[i];
}


static void mat2mat(int n,complex *a,complex *b)
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
   int n;
   float d,d1,d2,d3,d4,d5;
   complex *a1,*b1,*c1,*v1,*w1;
   complex *a2,*b2,*c2,*v2,*w2;

   printf("\n");
   printf("Check of cmat_vec, cmat_add, ...\n");
   printf("--------------------------------\n\n");

#if (defined AVX)
#if (defined FMA3)
   printf("Using AVX and FMA3 instructions\n\n");
#else
   printf("Using AVX instructions\n\n");
#endif
#elif (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   a1=amalloc(2*NMAX*(3*NMAX+2)*sizeof(*a1),4);
   error(a1==NULL,1,"main [check1.c]","Unable to allocate auxiliary arrays");

   b1=a1+NMAX*NMAX;
   c1=b1+NMAX*NMAX;
   v1=c1+NMAX*NMAX;
   w1=v1+NMAX;
   a2=w1+NMAX;
   b2=a2+NMAX*NMAX;
   c2=b2+NMAX*NMAX;
   v2=c2+NMAX*NMAX;
   w2=v2+NMAX;

   d1=0.0f;
   d2=0.0f;
   d3=0.0f;
   d4=0.0f;
   d5=0.0f;

   for (n=1;n<=NMAX;n++)
   {
      rvec(n,v1);
      rmat(n,a1);
      vec2vec(n,v1,v2);
      mat2mat(n,a1,a2);

      cmat_vec(n,a1,v1,w1);
      mvec(n,a2,v2,w2);

      d=vdev(n,w1,w2);
      if (d>d1)
         d1=d;

      error((mdev(n,a1,a2)!=0.0f)||(vdev(n,v1,v2)!=0.0f),1,"main [check1.c]",
            "cmat_vec: input values have changed");

      rvec(n,v1);
      rvec(n,w1);
      rmat(n,a1);
      vec2vec(n,v1,v2);
      vec2vec(n,w1,w2);
      mat2mat(n,a1,a2);

      cmat_vec_assign(n,a1,v1,w1);
      mvec_assign(n,a2,v2,w2);

      d=vdev(n,w1,w2);
      if (d>d1)
         d1=d;

      error((mdev(n,a1,a2)!=0.0f)||(vdev(n,v1,v2)!=0.0f),1,"main [check1.c]",
            "cmat_vec_assign: input values have changed");

      rmat(n,a1);
      rmat(n,b1);
      rmat(n,c1);
      mat2mat(n,a1,a2);
      mat2mat(n,b1,b2);
      rmat(n,c2);

      cmat_add(n,a1,b1,c1);
      madd(n,a2,b2,c2);

      d=mdev(n,c1,c2);
      if (d>d2)
         d2=d;

      error((mdev(n,a1,a2)!=0.0f)||(mdev(n,b1,b2)!=0.0f),1,"main [check1.c]",
            "cmat_add: input values have changed");

      rmat(n,a1);
      rmat(n,b1);
      rmat(n,c1);
      mat2mat(n,a1,a2);
      mat2mat(n,b1,b2);
      rmat(n,c2);

      cmat_sub(n,a1,b1,c1);
      msub(n,a2,b2,c2);

      d=mdev(n,c1,c2);
      if (d>d3)
         d3=d;

      error((mdev(n,a1,a2)!=0.0f)||(mdev(n,b1,b2)!=0.0f),1,"main [check1.c]",
            "cmat_sub: input values have changed");

      rmat(n,a1);
      rmat(n,b1);
      rmat(n,c1);
      mat2mat(n,a1,a2);
      mat2mat(n,b1,b2);
      rmat(n,c2);

      cmat_mul(n,a1,b1,c1);
      mmul(n,a2,b2,c2);

      d=mdev(n,c1,c2);
      if (d>d4)
         d4=d;

      error((mdev(n,a1,a2)!=0.0f)||(mdev(n,b1,b2)!=0.0f),1,"main [check1.c]",
            "cmat_mul: input values have changed");

      rmat(n,a1);
      rmat(n,b1);
      mat2mat(n,a1,a2);
      rmat(n,b2);

      cmat_dag(n,a1,b1);
      mdag(n,a2,b2);

      d=mdev(n,b1,b2);
      if (d>d5)
         d5=d;

      error(mdev(n,a1,a2)!=0.0f,1,"main [check1.c]",
            "cmat_dag: input values have changed");
   }

   printf("Consider matrices of size up to %dx%d\n\n",NMAX,NMAX);

   printf("The maximal observed deviations are:\n\n");
   printf("cmat_vec:     %.1e\n",d1);
   printf("cmat_add:     %.1e\n",d2);
   printf("cmat_sub:     %.1e\n",d3);
   printf("cmat_mul:     %.1e\n",d4);
   printf("cmat_dag:     %.1e\n\n",d5);

   exit(0);
}
