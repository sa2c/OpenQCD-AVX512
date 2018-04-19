
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005, 2009, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of mul_pauli_dble()
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

typedef union
{
   weyl_dble w;
   complex_dble c[6];
   double r[12];
} spin_t;

typedef union
{
   complex_dble c[36];
   double r[72];
} mat_t;

#if (defined AVX)
static pauli_dble mp ALIGNED32;
static spin_t s1 ALIGNED32;
static spin_t s2 ALIGNED32;
static spin_t r1 ALIGNED32;
static spin_t r2 ALIGNED32;
static mat_t mv ALIGNED32;
#else
static pauli_dble mp ALIGNED16;
static spin_t s1 ALIGNED16;
static spin_t s2 ALIGNED16;
static spin_t r1 ALIGNED16;
static spin_t r2 ALIGNED16;
static mat_t mv ALIGNED16;
#endif


static void cpvec(int n,complex_dble *s,complex_dble *r)
{
   int i;

   for (i=0;i<n;i++)
      r[i]=s[i];
}


static int diffvec(int n,complex_dble *s,complex_dble *r)
{
   int i;

   for (i=0;i<n;i++)
   {
      if ((s[i].re!=r[i].re)||(s[i].im!=r[i].im))
         return 1;
   }

   return 0;
}


int main(void)
{
   int i,j,k;
   double mu,d,dmax;

   printf("\n");
   printf("Check of mul_pauli_dble()\n");
   printf("-------------------------\n\n");

#if (defined AVX)
   printf("Using AVX instructions\n\n");
#elif (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   rlxd_init(1,3898);
   gauss_dble(s1.r,12);
   gauss_dble(mv.r,72);
   mu=0.1234;

   for (i=0;i<6;i++)
   {
      for (j=i;j<6;j++)
      {
         if (j>i)
         {
            mv.c[6*i+j].re= mv.c[6*j+i].re;
            mv.c[6*i+j].im=-mv.c[6*j+i].im;
         }
         else
            mv.c[6*i+j].im=0.0;
      }
   }

   k=6;

   for (i=0;i<6;i++)
   {
      mp.u[i]=mv.c[6*i+i].re;

      for (j=i+1;j<6;j++)
      {
         mp.u[k]=mv.c[6*i+j].re;
         k+=1;
         mp.u[k]=mv.c[6*i+j].im;
         k+=1;
      }
   }

   cpvec(6,s1.c,s2.c);
   mul_pauli_dble(mu,&mp,&(s1.w),&(r1.w));

   error(diffvec(6,s1.c,s2.c),1,"main [check2.c]",
         "mul_pauli_dble() modifies the source spinor");

   cmat_vec_dble(6,mv.c,s2.c,r2.c);

   for (i=0;i<6;i++)
   {
      r2.c[i].re-=mu*s2.c[i].im;
      r2.c[i].im+=mu*s2.c[i].re;
   }

   printf("r1: result, r2: expected result\n\n");

   for (i=0;i<2;i++)
   {
      for (j=0;j<3;j++)
      {
         k=3*i+j;
         printf("r1.c%d.c%d=(% .7e,% .7e)\n",i+1,j+1,r1.c[k].re,r1.c[k].im);
         printf("r2.c%d.c%d=(% .7e,% .7e)\n",i+1,j+1,r2.c[k].re,r2.c[k].im);
         printf("\n");
      }
   }

   dmax=0.0;

   for (i=0;i<12;i++)
   {
      d=fabs(r1.r[i]-r2.r[i]);
      if (d>dmax)
         dmax=d;
   }

   printf("Maximal absolute deviation = %.1e\n",dmax);

   mul_pauli_dble(mu,&mp,&(s1.w),&(s1.w));
   error(diffvec(6,s1.c,r1.c),1,"main [check2.c]",
         "mul_pauli_dble() is incorrect when r=s");
   printf("Works correctly if input and output spinors coincide\n\n");

   exit(0);
}
