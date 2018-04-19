
/*******************************************************************************
*
* File time2.c
*
* Copyright (C) 2005, 2008, 2009, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of the SU(3) x SU(3)-vector multiplication (double-precision programs)
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "random.h"
#include "su3fcts.h"

static su3_dble u[4] ALIGNED16;
static su3_vector_dble s[8] ALIGNED16;
static su3_vector_dble r[8] ALIGNED16;
static su3_vector_dble t[8] ALIGNED16;

#if (defined x64)
#if (defined AVX)
#include "avx.h"

#define _su3_fast_multiply(r1,r2,u,s1,s2) \
   _avx_pair_load_dble(s1,s2); \
   _avx_su3_multiply_pair_dble(u); \
   _avx_pair_store_up_dble(r1,r2)

#define _su3_fast_inverse_multiply(r1,r2,u,s1,s2) \
   _avx_pair_load_dble(s1,s2); \
   _avx_su3_inverse_multiply_pair_dble(u); \
   _avx_pair_store_up_dble(r1,r2)

static void fast_multiply(su3_dble *ua,su3_vector_dble *sa,
                          su3_vector_dble *ra)
{
   _su3_fast_multiply((*(ra  )),(*(ra+1)),(*(ua  )),(*(sa  )),(*(sa+1)));
   _su3_fast_multiply((*(ra+2)),(*(ra+3)),(*(ua+1)),(*(sa+2)),(*(sa+3)));
   _su3_fast_multiply((*(ra+4)),(*(ra+5)),(*(ua+2)),(*(sa+4)),(*(sa+5)));
   _su3_fast_multiply((*(ra+6)),(*(ra+7)),(*(ua+3)),(*(sa+6)),(*(sa+7)));
}


static void fast_inverse_multiply(su3_dble *ua,su3_vector_dble *sa,
                                  su3_vector_dble *ra)
{
   _su3_fast_inverse_multiply((*(ra  )),(*(ra+1)),(*(ua  )),(*(sa  )),(*(sa+1)));
   _su3_fast_inverse_multiply((*(ra+2)),(*(ra+3)),(*(ua+1)),(*(sa+2)),(*(sa+3)));
   _su3_fast_inverse_multiply((*(ra+4)),(*(ra+5)),(*(ua+2)),(*(sa+4)),(*(sa+5)));
   _su3_fast_inverse_multiply((*(ra+6)),(*(ra+7)),(*(ua+3)),(*(sa+6)),(*(sa+7)));
}

#else
#include "sse2.h"

#define _su3_fast_multiply(r,u,s) \
   _sse_load_dble(s); \
   _sse_su3_multiply_dble(u); \
   _sse_store_up_dble(r)

#define _su3_fast_inverse_multiply(r,u,s) \
   _sse_load_dble(s); \
   _sse_su3_inverse_multiply_dble(u); \
   _sse_store_up_dble(r)


static void fast_multiply(su3_dble *ua,su3_vector_dble *sa,
                          su3_vector_dble *ra)
{
   _su3_fast_multiply((*(ra  )),(*(ua  )),(*(sa  )));
   _su3_fast_multiply((*(ra+1)),(*(ua  )),(*(sa+1)));
   _su3_fast_multiply((*(ra+2)),(*(ua+1)),(*(sa+2)));
   _su3_fast_multiply((*(ra+3)),(*(ua+1)),(*(sa+3)));
   _su3_fast_multiply((*(ra+4)),(*(ua+2)),(*(sa+4)));
   _su3_fast_multiply((*(ra+5)),(*(ua+2)),(*(sa+5)));
   _su3_fast_multiply((*(ra+6)),(*(ua+3)),(*(sa+6)));
   _su3_fast_multiply((*(ra+7)),(*(ua+3)),(*(sa+7)));
}


static void fast_inverse_multiply(su3_dble *ua,su3_vector_dble *sa,
                                  su3_vector_dble *ra)
{
   _su3_fast_inverse_multiply((*(ra  )),(*(ua  )),(*(sa  )));
   _su3_fast_inverse_multiply((*(ra+1)),(*(ua  )),(*(sa+1)));
   _su3_fast_inverse_multiply((*(ra+2)),(*(ua+1)),(*(sa+2)));
   _su3_fast_inverse_multiply((*(ra+3)),(*(ua+1)),(*(sa+3)));
   _su3_fast_inverse_multiply((*(ra+4)),(*(ua+2)),(*(sa+4)));
   _su3_fast_inverse_multiply((*(ra+5)),(*(ua+2)),(*(sa+5)));
   _su3_fast_inverse_multiply((*(ra+6)),(*(ua+3)),(*(sa+6)));
   _su3_fast_inverse_multiply((*(ra+7)),(*(ua+3)),(*(sa+7)));
}

#endif

static void slow_inverse_multiply(su3_dble *ua,su3_vector_dble *sa,
                                  su3_vector_dble *ra)
{
   _su3_inverse_multiply((*(ra  )),(*(ua  )),(*(sa  )));
   _su3_inverse_multiply((*(ra+1)),(*(ua  )),(*(sa+1)));
   _su3_inverse_multiply((*(ra+2)),(*(ua+1)),(*(sa+2)));
   _su3_inverse_multiply((*(ra+3)),(*(ua+1)),(*(sa+3)));
   _su3_inverse_multiply((*(ra+4)),(*(ua+2)),(*(sa+4)));
   _su3_inverse_multiply((*(ra+5)),(*(ua+2)),(*(sa+5)));
   _su3_inverse_multiply((*(ra+6)),(*(ua+3)),(*(sa+6)));
   _su3_inverse_multiply((*(ra+7)),(*(ua+3)),(*(sa+7)));
}

#endif

static void slow_multiply(su3_dble *ua,su3_vector_dble *sa,
                          su3_vector_dble *ra)
{
   _su3_multiply((*(ra  )),(*(ua  )),(*(sa  )));
   _su3_multiply((*(ra+1)),(*(ua  )),(*(sa+1)));
   _su3_multiply((*(ra+2)),(*(ua+1)),(*(sa+2)));
   _su3_multiply((*(ra+3)),(*(ua+1)),(*(sa+3)));
   _su3_multiply((*(ra+4)),(*(ua+2)),(*(sa+4)));
   _su3_multiply((*(ra+5)),(*(ua+2)),(*(sa+5)));
   _su3_multiply((*(ra+6)),(*(ua+3)),(*(sa+6)));
   _su3_multiply((*(ra+7)),(*(ua+3)),(*(sa+7)));
}


int main(void)
{
   int k,n,count;
   double t1,t2,dt;
#if (defined x64)
   double delta,diff,norm;
#endif

   printf("\n");
   printf("Time per double-precision SU(3) x SU(3)-vector multiplication\n");
   printf("-------------------------------------------------------------\n\n");

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

   rlxd_init(1,123456);

   for (k=0;k<4;k++)
      random_su3_dble(u+k);

   gauss_dble((double*)(s),48);
   gauss_dble((double*)(r),48);
   gauss_dble((double*)(t),48);

#if (defined x64)

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         fast_multiply(u,s,r);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=1.0e6/(double)(4*n);

   printf("The time per v=U*w is     %4.3f nsec",1.0e3*dt);
   printf(" [%d Mflops]\n",(int)(66.0/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         fast_inverse_multiply(u,s,r);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=1.0e6/(double)(4*n);

   printf("The time per v=U^dag*w is %4.3f nsec",1.0e3*dt);
   printf(" [%d Mflops]\n",(int)(66.0/dt));

#endif

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         slow_multiply(u,s,t);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=1.0e6/(double)(4*n);

   printf("Using x87 FPU instructions:\n");
   printf("The time per v=U*w is     %4.3f nsec",1.0e3*dt);
   printf(" [%d Mflops]\n",(int)(66.0/dt));

#if (defined x64)

   fast_multiply(u,s,r);
   slow_multiply(u,s,t);
   delta=0.0;

   for (k=0;k<8;k++)
   {
      _vector_sub_assign(r[k],t[k]);
      diff=_vector_prod_re(r[k],r[k]);
      norm=_vector_prod_re(s[k],s[k]);
      diff=sqrt(diff/norm);
      if (diff>delta)
         delta=diff;
   }

#if (defined AVX)
   printf("||U*w_AVX-U*w_FPU||<= %.1e*||w||\n",delta);
#else
   printf("||U*w_SSE-U*w_FPU||<= %.1e*||w||\n",delta);
#endif

   fast_inverse_multiply(u,s,r);
   slow_inverse_multiply(u,s,t);
   delta=0.0;

   for (k=0;k<8;k++)
   {
      _vector_sub_assign(r[k],t[k]);
      diff=_vector_prod_re(r[k],r[k]);
      norm=_vector_prod_re(s[k],s[k]);
      diff=sqrt(diff/norm);
      if (diff>delta)
         delta=diff;
   }

#if (defined AVX)
   printf("||U^dag*w_AVX-U^dag*w_FPU||<= %.1e*||w||\n",delta);
#else
   printf("||U^dag*w_SSE-U^dag*w_FPU||<= %.1e*||w||\n",delta);
#endif
#endif

   exit(0);
}
