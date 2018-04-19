
/*******************************************************************************
*
* File check7.c
*
* Copyright (C) 2005, 2009, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of det_pauli_dble()
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "utils.h"
#include "random.h"
#include "linalg.h"
#include "sw_term.h"

#define NM 10000

static double dd[6] ALIGNED16;
static complex_dble aa[36] ALIGNED16;
static complex_dble bb[36] ALIGNED16;
static complex_dble vv[36] ALIGNED16;
static complex_dble ww[36] ALIGNED16;


static complex_dble random_dd(double mu)
{
   int i;
   complex_dble det,z;

   ranlxd(dd,6);
   det.re=1.0;
   det.im=0.0;

   for (i=0;i<6;i++)
   {
      if (dd[i]<0.5)
         dd[i]-=0.6;
      else
         dd[i]-=0.4;

      z.re=det.re*dd[i]-det.im*mu;
      z.im=det.re*mu+det.im*dd[i];

      det.re=z.re;
      det.im=z.im;
   }

   return det;
}


static double norm(complex_dble *v)
{
   int i;
   double r;

   r=0.0;

   for (i=0;i<6;i++)
      r+=(v[i].re*v[i].re+v[i].im*v[i].im);

   return sqrt(r);
}


static complex_dble prod(complex_dble *v,complex_dble *w)
{
   int i;
   complex_dble z;

   z.re=0.0;
   z.im=0.0;

   for (i=0;i<6;i++)
   {
      z.re+=(v[i].re*w[i].re+v[i].im*w[i].im);
      z.im+=(v[i].re*w[i].im-v[i].im*w[i].re);
   }

   return z;
}


static void proj(complex_dble *v,complex_dble *w)
{
   int i;
   complex_dble z;

   z=prod(v,w);

   for (i=0;i<6;i++)
   {
      w[i].re-=(z.re*v[i].re-z.im*v[i].im);
      w[i].im-=(z.re*v[i].im+z.im*v[i].re);
   }
}


static void random_vv(void)
{
   int i,j;
   double r,ri[12];
   complex_dble *vi;

   for (i=0;i<6;i++)
   {
      vi=vv+6*i;
      r=0.0;

      while (r<1.0)
      {
         gauss_dble(ri,12);

         for (j=0;j<6;j++)
         {
            vi[j].re=ri[2*j];
            vi[j].im=ri[2*j+1];
         }

         for (j=0;j<i;j++)
            proj(vv+6*j,vi);

         for (j=0;j<i;j++)
            proj(vv+6*j,vi);

         r=norm(vi);
      }

      for (j=0;j<6;j++)
      {
         vi[j].re/=r;
         vi[j].im/=r;
      }
   }
}


static complex_dble random_aa(double mu)
{
   int i,j;
   complex_dble det;

   det=random_dd(mu);
   random_vv();

   for (i=0;i<6;i++)
   {
      for (j=0;j<6;j++)
      {
         if (i==j)
         {
            aa[6*i+j].re=dd[i];
            aa[6*i+j].im=mu;
         }
         else
         {
            aa[6*i+j].re=0.0;
            aa[6*i+j].im=0.0;
         }
      }
   }

   cmat_mul_dble(6,aa,vv,bb);
   cmat_dag_dble(6,vv,ww);
   cmat_mul_dble(6,ww,bb,aa);

   return det;
}


static void aa2pauli(pauli_dble *m)
{
   int i,j,k;
   double *u;

   u=(*m).u;
   k=6;

   for (i=0;i<6;i++)
   {
      u[i]=aa[6*i+i].re;

      for (j=i+1;j<6;j++)
      {
         u[k]=aa[6*i+j].re;
         k+=1;
         u[k]=aa[6*i+j].im;
         k+=1;
      }
   }
}


int main(void)
{
   int n;
   double mu,d,dmax;
   complex_dble det1,det2;
   pauli_dble *md;

   printf("\n");
   printf("Check of det_pauli_dble()\n");
   printf("-------------------------\n\n");

#if (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   rlxd_init(1,1234);
   md=amalloc(NM*sizeof(pauli_dble),4);
   error(md==NULL,1,"main [check7.c]",
         "Unable to allocate auxiliary arrays");
   mu=0.1234;
   dmax=0.0;

   for (n=0;n<NM;n++)
   {
      det1=random_aa(mu);
      aa2pauli(md);
      det2=det_pauli_dble(mu,md);

      det2.re-=det1.re;
      det2.im-=det1.im;

      d=sqrt((det2.re*det2.re+det2.im*det2.im)/
             (det1.re*det1.re+det1.im*det1.im));

      if (d>dmax)
         dmax=d;

      md+=1;
   }

   printf("%d Gaussian random matrices M, mu=%.4f\n",NM,mu);
   printf("Maximal relative deviation of det(M+i*mu) = %.1e\n\n",dmax);
   exit(0);
}
