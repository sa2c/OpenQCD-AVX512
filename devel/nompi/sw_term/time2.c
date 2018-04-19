
/*******************************************************************************
*
* File time2.c
*
* Copyright (C) 2005, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of mul_pauli_dble()
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "utils.h"
#include "random.h"
#include "linalg.h"
#include "sw_term.h"

typedef union
{
   weyl_dble w;
   double r[12];
} spin_t;

#if (defined AVX)
static pauli_dble mp1 ALIGNED32;
static pauli_dble mp2 ALIGNED32;
static spin_t s1 ALIGNED32;
static spin_t s2 ALIGNED32;
static spin_t r1 ALIGNED32;
static spin_t r2 ALIGNED32;
#else
static pauli_dble mp1 ALIGNED16;
static pauli_dble mp2 ALIGNED16;
static spin_t s1 ALIGNED16;
static spin_t s2 ALIGNED16;
static spin_t r1 ALIGNED16;
static spin_t r2 ALIGNED16;
#endif


int main(void)
{
   int n,count;
   double mu1,mu2;
   double t1,t2,dt;

   printf("\n");
   printf("Timing of mul_pauli_dble()\n");
   printf("--------------------------\n\n");

#if (defined AVX)
   printf("Using AVX instructions\n\n");
#elif (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   printf("Measurement made with all data in cache\n\n");

   rlxd_init(1,23456);
   ranlxd(mp1.u,36);
   ranlxd(mp2.u,36);
   ranlxd(s1.r,12);
   ranlxd(s2.r,12);
   mu1=0.1234;
   mu2=0.5678;

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
      {
         mul_pauli_dble(mu1,&mp1,&(s1.w),&(r1.w));
         mul_pauli_dble(mu2,&mp2,&(s2.w),&(r2.w));
      }
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=(1.0e6/(double)(n));

   printf("Time per call of mul_pauli_dble():\n");
   printf("%.4f usec (%d Mflops [%d bit arithmetic])\n\n",
          dt,(int)(276.0/dt),(int)(sizeof(spinor_dble)/3));

   exit(0);
}
