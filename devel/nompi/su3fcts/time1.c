
/*******************************************************************************
*
* File time1.c
*
* Copyright (C) 2005, 2008, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of the SU(3) x SU(3)-vector multiplication (single-precision programs)
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "random.h"
#include "su3fcts.h"

static su3 u[4] ALIGNED16;
static su3_vector s[8] ALIGNED16;
static su3_vector r[8] ALIGNED16;
static su3_vector t[8] ALIGNED16;

#if (defined x64)
#if (defined AVX)
#include "avx.h"

#define _avx_vector_quartet_load(s)              \
__asm__ __volatile__ ("vmovaps %0, %%xmm0 \n\t" \
                      "vmovaps %2, %%xmm1 \n\t" \
                      "vmovaps %4, %%xmm2" \
                      : \
                      : \
                      "m" ((s[0]).c1), \
                      "m" ((s[0]).c2), \
                      "m" ((s[0]).c3), \
                      "m" ((s[1]).c1), \
                      "m" ((s[1]).c2), \
                      "m" ((s[1]).c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm0, %%ymm0 \n\t" \
                      "vinsertf128 $0x1, %2, %%ymm1, %%ymm1 \n\t" \
                      "vinsertf128 $0x1, %4, %%ymm2, %%ymm2" \
                      : \
                      : \
                      "m" ((s[2]).c1), \
                      "m" ((s[2]).c2), \
                      "m" ((s[2]).c3), \
                      "m" ((s[3]).c1), \
                      "m" ((s[3]).c2), \
                      "m" ((s[3]).c3) \
                      : \
                      "xmm6", "xmm7", "xmm8")

#define _avx_vector_quartet_store_up(r) \
__asm__ __volatile__ ("vmovaps %%xmm3, %0 \n\t" \
                      "vmovaps %%xmm4, %2 \n\t" \
                      "vmovaps %%xmm5, %4" \
                      : \
                      : \
                      "m" ((r[0]).c1), \
                      "m" ((r[0]).c2), \
                      "m" ((r[0]).c3), \
                      "m" ((r[1]).c1), \
                      "m" ((r[1]).c2), \
                      "m" ((r[1]).c3)); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm3, %0 \n\t" \
                      "vextractf128 $0x1, %%ymm4, %2 \n\t" \
                      "vextractf128 $0x1, %%ymm5, %4" \
                      : \
                      : \
                      "m" ((r[2]).c1), \
                      "m" ((r[2]).c2), \
                      "m" ((r[2]).c3), \
                      "m" ((r[3]).c1), \
                      "m" ((r[3]).c2), \
                      "m" ((r[3]).c3))


#define _avx_vector_quartet_save(s) \
__asm__ __volatile__ ("vshufps $0x44, %%ymm1, %%ymm0, %%ymm9 \n\t" \
                      "vshufps $0xe4, %%ymm0, %%ymm2, %%ymm10 \n\t" \
                      "vshufps $0xee, %%ymm2, %%ymm1, %%ymm11" \
                      : \
                      : \
                      : \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vmovaps %%xmm9, %0 \n\t" \
                      "vmovaps %%xmm10, %2 \n\t" \
                      "vmovaps %%xmm11, %4" \
                      : \
                      : \
                      "m" ((s[0]).c1), \
                      "m" ((s[0]).c2), \
                      "m" ((s[0]).c3), \
                      "m" ((s[1]).c1), \
                      "m" ((s[1]).c2), \
                      "m" ((s[1]).c3)); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm9, %0 \n\t" \
                      "vextractf128 $0x1, %%ymm10, %2 \n\t" \
                      "vextractf128 $0x1, %%ymm11, %4" \
                      : \
                      : \
                      "m" ((s[2]).c1), \
                      "m" ((s[2]).c2), \
                      "m" ((s[2]).c3), \
                      "m" ((s[3]).c1), \
                      "m" ((s[3]).c2), \
                      "m" ((s[3]).c3))


static void swizzle_spinors(su3_vector *sa,su3_vector *ra)
{
   _avx_vector_quartet_load(sa);
   _avx_vector_quartet_save(sa);
   _avx_vector_quartet_load(ra);
   _avx_vector_quartet_save(ra);
   sa+=4;
   ra+=4;
   _avx_vector_quartet_load(sa);
   _avx_vector_quartet_save(sa);
   _avx_vector_quartet_load(ra);
   _avx_vector_quartet_save(ra);
}

static void fast_multiply(su3 *ua,su3_vector *sa,su3_vector *ra)
{
   _avx_vector_quartet_load(sa);
   _avx_su3_pair_multiply(ua[0],ua[1]);
   _avx_vector_quartet_store_up(ra);
   ua+=2;
   sa+=4;
   ra+=4;
   _avx_vector_quartet_load(sa);
   _avx_su3_pair_multiply(ua[0],ua[1]);
   _avx_vector_quartet_store_up(ra);
}


static void fast_inverse_multiply(su3 *ua,su3_vector *sa,su3_vector *ra)
{
   _avx_vector_quartet_load(sa);
   _avx_su3_pair_inverse_multiply(ua[0],ua[1]);
   _avx_vector_quartet_store_up(ra);
   ua+=2;
   sa+=4;
   ra+=4;
   _avx_vector_quartet_load(sa);
   _avx_su3_pair_inverse_multiply(ua[0],ua[1]);
   _avx_vector_quartet_store_up(ra);
}


static void fast_mixed_multiply(su3 *ua,su3_vector *sa,su3_vector *ra)
{
   _avx_vector_quartet_load(sa);
   _avx_su3_pair_mixed_multiply(ua[0],ua[1]);
   _avx_vector_quartet_store_up(ra);
   ua+=2;
   sa+=4;
   ra+=4;
   _avx_vector_quartet_load(sa);
   _avx_su3_pair_mixed_multiply(ua[0],ua[1]);
   _avx_vector_quartet_store_up(ra);
}


static void slow_mixed_multiply(su3 *ua,su3_vector *sa,su3_vector *ra)
{
   _su3_multiply((*(ra  )),(*(ua  )),(*(sa  )));
   _su3_multiply((*(ra+1)),(*(ua  )),(*(sa+1)));
   _su3_inverse_multiply((*(ra+2)),(*(ua+1)),(*(sa+2)));
   _su3_inverse_multiply((*(ra+3)),(*(ua+1)),(*(sa+3)));
   _su3_multiply((*(ra+4)),(*(ua+2)),(*(sa+4)));
   _su3_multiply((*(ra+5)),(*(ua+2)),(*(sa+5)));
   _su3_inverse_multiply((*(ra+6)),(*(ua+3)),(*(sa+6)));
   _su3_inverse_multiply((*(ra+7)),(*(ua+3)),(*(sa+7)));
}

#else
#include "sse2.h"

#define _su3_fast_multiply(r1,r2,u,s1,s2) \
   _sse_pair_load(s1,s2); \
   _sse_su3_multiply(u); \
   _sse_pair_store_up(r1,r2)

#define _su3_fast_inverse_multiply(r1,r2,u,s1,s2) \
   _sse_pair_load(s1,s2); \
   _sse_su3_inverse_multiply(u); \
   _sse_pair_store_up(r1,r2)


static void fast_multiply(su3 *ua,su3_vector *sa,su3_vector *ra)
{
   _su3_fast_multiply((*(ra  )),(*(ra+1)),(*(ua  )),(*(sa  )),(*(sa+1)));
   _su3_fast_multiply((*(ra+2)),(*(ra+3)),(*(ua+1)),(*(sa+2)),(*(sa+3)));
   _su3_fast_multiply((*(ra+4)),(*(ra+5)),(*(ua+2)),(*(sa+4)),(*(sa+5)));
   _su3_fast_multiply((*(ra+6)),(*(ra+7)),(*(ua+3)),(*(sa+6)),(*(sa+7)));
}


static void fast_inverse_multiply(su3 *ua,su3_vector *sa,su3_vector *ra)
{
   _su3_fast_inverse_multiply((*(ra  )),(*(ra+1)),(*(ua  )),(*(sa  )),(*(sa+1)));
   _su3_fast_inverse_multiply((*(ra+2)),(*(ra+3)),(*(ua+1)),(*(sa+2)),(*(sa+3)));
   _su3_fast_inverse_multiply((*(ra+4)),(*(ra+5)),(*(ua+2)),(*(sa+4)),(*(sa+5)));
   _su3_fast_inverse_multiply((*(ra+6)),(*(ra+7)),(*(ua+3)),(*(sa+6)),(*(sa+7)));
}

#endif

static void slow_inverse_multiply(su3 *ua,su3_vector *sa,su3_vector *ra)
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

static void slow_multiply(su3 *ua,su3_vector *sa,su3_vector *ra)
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
   printf("Time per single-precision SU(3) x SU(3)-vector multiplication\n");
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

   rlxs_init(0,123456);

   for (k=0;k<4;k++)
      random_su3(u+k);

   gauss((float*)(s),48);
   gauss((float*)(r),48);
   gauss((float*)(t),48);

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

#if (defined AVX)

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         fast_mixed_multiply(u,s,r);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=1.0e6/(double)(4*n);

   printf("The time per v=U/U^dag*w is %4.3f nsec",1.0e3*dt);
   printf(" [%d Mflops]\n",(int)(66.0/dt));

#endif
#endif

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<(n/2);count++)
         slow_multiply(u,s,t);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=1.0e6/(double)(2*n);

   printf("Using x87 FPU instructions:\n");
   printf("The time per v=U*w is     %4.3f nsec",1.0e3*dt);
   printf(" [%d Mflops]\n",(int)(66.0/dt));

#if (defined x64)

   fast_multiply(u,s,r);
#if (defined AVX)
   swizzle_spinors(s,r);
#endif
   slow_multiply(u,s,t);
   delta=0.0;

   for (k=0;k<8;k++)
   {
      _vector_sub_assign(r[k],t[k]);
      diff=(double)(_vector_prod_re(r[k],r[k]));
      norm=(double)(_vector_prod_re(s[k],s[k]));
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
#if (defined AVX)
   swizzle_spinors(s,r);
#endif
   slow_inverse_multiply(u,s,t);
   delta=0.0;

   for (k=0;k<8;k++)
   {
      _vector_sub_assign(r[k],t[k]);
      diff=(double)(_vector_prod_re(r[k],r[k]));
      norm=(double)(_vector_prod_re(s[k],s[k]));
      diff=sqrt(diff/norm);
      if (diff>delta)
         delta=diff;
   }

#if (defined AVX)
   printf("||U^dag*w_AVX-U^dag*w_FPU||<= %.1e*||w||\n",delta);
#else
   printf("||U^dag*w_SSE-U^dag*w_FPU||<= %.1e*||w||\n",delta);
#endif

#if (defined AVX)

   fast_mixed_multiply(u,s,r);
   swizzle_spinors(s,r);
   slow_mixed_multiply(u,s,t);
   delta=0.0;

   for (k=0;k<8;k++)
   {
      _vector_sub_assign(r[k],t[k]);
      diff=(double)(_vector_prod_re(r[k],r[k]));
      norm=(double)(_vector_prod_re(s[k],s[k]));
      diff=sqrt(diff/norm);
      if (diff>delta)
         delta=diff;
   }

   printf("||U/U^dag*w_AVX-U/U^dag*w_FPU||<= %.1e*||w||\n",delta);

#endif
#endif

   exit(0);
}
