
/*******************************************************************************
*
* File sse2.h
*
* Copyright (C) 2005, 2008, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Macros for Dirac spinors, SU(3) vectors and SU(3) matrices using inline
* assembly SSE3 instructions. The machine is assumed to comply with the
* x86-64 instruction set.
*
*******************************************************************************/

#ifndef SSE2_H
#define SSE2_H

#ifndef SSE_H
#include "sse.h"
#endif

typedef struct __attribute__ ((aligned (16)))
{
   double c1,c2;
} sse_double;

static sse_double _sse_sgn1_dble __attribute__ ((unused)) ={-1.0,1.0};
static sse_double _sse_sgn2_dble __attribute__ ((unused)) ={1.0,-1.0};
static sse_double _sse_sgn_dble __attribute__ ((unused)) ={-1.0,-1.0};

/*******************************************************************************
*
* Macros for double-precision su3 vectors
*
* Most of these macros operate on su3 vectors that are stored
* in  xmm0,xmm1,xmm2 or xmm3,xmm4,xmm5. For example,
*
* xmm0 -> s.c1.re,s.c1.im
* xmm1 -> s.c2.re,s.c2.im
* xmm2 -> s.c3.re,s.c3.im
*
* where s is of type su3_vector_dble
*
*******************************************************************************/

/*
* Loads an su3 vector s to xmm0,xmm1,xmm2
*/

#define _sse_load_dble(s) \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Loads an su3 vector s to xmm3,xmm4,xmm5
*/

#define _sse_load_up_dble(s) \
__asm__ __volatile__ ("movapd %0, %%xmm3 \n\t" \
                      "movapd %1, %%xmm4 \n\t" \
                      "movapd %2, %%xmm5" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Stores xmm0,xmm1,xmm2 to the components r.c1,r.c2,r.c3 of an su3 vector
*/

#define _sse_store_dble(r) \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((r).c1), \
                      "=m" ((r).c2), \
                      "=m" ((r).c3))

/*
* Stores xmm3,xmm4,xmm5 to the components r.c1,r.c2,r.c3 of an su3 vector
*/

#define _sse_store_up_dble(r) \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((r).c1), \
                      "=m" ((r).c2), \
                      "=m" ((r).c3))

/*
* Multiplies xmm0,xmm1,xmm2 with a constant sse_double c
*/

#define _sse_vector_mul_dble(c) \
__asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t" \
                      "mulpd %0, %%xmm1 \n\t" \
                      "mulpd %0, %%xmm2" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Multiplies xmm3,xmm4,xmm5 with a constant sse_double c
*/

#define _sse_vector_mul_up_dble(c) \
__asm__ __volatile__ ("mulpd %0, %%xmm3 \n\t" \
                      "mulpd %0, %%xmm4 \n\t" \
                      "mulpd %0, %%xmm5" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Adds xmm3,xmm4,xmm5 to xmm0,xmm1,xmm2
*/

#define _sse_vector_add_dble() \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Subtracts xmm3,xmm4,xmm5 from xmm0,xmm1,xmm2
*/

#define _sse_vector_sub_dble() \
__asm__ __volatile__ ("subpd %%xmm3, %%xmm0 \n\t" \
                      "subpd %%xmm4, %%xmm1 \n\t" \
                      "subpd %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Multiplies xmm3,xmm4,xmm5 with i
*/

#define _sse_vector_i_mul_dble() \
__asm__ __volatile__ ("shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
                      "shufpd $0x1, %%xmm4, %%xmm4 \n\t" \
                      "shufpd $0x1, %%xmm5, %%xmm5 \n\t" \
                      "mulpd %0, %%xmm3 \n\t" \
                      "mulpd %0, %%xmm4 \n\t" \
                      "mulpd %0, %%xmm5" \
                      : \
                      : \
                      "m" (_sse_sgn1_dble) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies xmm3,xmm4,xmm5 with i and adds them to xmm0,xmm1,xmm2
*/

#define _sse_vector_i_add_dble() \
__asm__ __volatile__ ("shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
                      "shufpd $0x1, %%xmm4, %%xmm4 \n\t" \
                      "shufpd $0x1, %%xmm5, %%xmm5 \n\t" \
                      "addsubpd %%xmm3, %%xmm0 \n\t" \
                      "addsubpd %%xmm4, %%xmm1 \n\t" \
                      "addsubpd %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
*  Loads (z.re,z.re) to xmm6 and (-z.im,z.im) to xmm7
*/

#define _sse_load_cmplx_dble(z) \
__asm__ __volatile__ ("movddup %0, %%xmm6 \n\t" \
                      "movddup %1, %%xmm7 \n\t" \
                      "mulpd %2, %%xmm7" \
                      : \
                      : \
                      "m" ((z).re), \
                      "m" ((z).im), \
                      "m" (_sse_sgn1_dble) \
                      : \
                      "xmm6", "xmm7")

/*
*  Multiplies the complex numbers in xmm0,xmm1,xmm2 by z, assuming z has
*  been loaded to xmm6,xmm7 by _sse_load_cmplx_dble(z). The result appears
*  in xmm0,xmm1,xmm2 and xmm3,xmm4,xmm5,xmm6,xmm7 are unchanged
*/

#define _sse_mulc_vector_dble() \
__asm__ __volatile__ ("movapd %%xmm0, %%xmm8 \n\t" \
                      "movapd %%xmm1, %%xmm9 \n\t" \
                      "movapd %%xmm2, %%xmm10 \n\t" \
                      "mulpd %%xmm6, %%xmm0 \n\t" \
                      "mulpd %%xmm6, %%xmm1 \n\t" \
                      "mulpd %%xmm6, %%xmm2 \n\t" \
                      "shufpd $0x1, %%xmm8, %%xmm8 \n\t" \
                      "shufpd $0x1, %%xmm9, %%xmm9 \n\t" \
                      "shufpd $0x1, %%xmm10, %%xmm10 \n\t" \
                      "mulpd %%xmm7, %%xmm8 \n\t" \
                      "mulpd %%xmm7, %%xmm9 \n\t" \
                      "mulpd %%xmm7, %%xmm10 \n\t" \
                      "addpd %%xmm8, %%xmm0 \n\t" \
                      "addpd %%xmm9, %%xmm1 \n\t" \
                      "addpd %%xmm10, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm8", "xmm9", "xmm10")

/*
*  Multiplies the complex numbers in xmm3,xmm4,xmm5 by z, assuming z has
*  been loaded to xmm6,xmm7 by _sse_load_cmplx_dble(z). The result appears
*  in xmm3,xmm4,xmm5 and xmm0,xmm1,xmm2,xmm6,xmm7 are unchanged
*/

#define _sse_mulc_vector_up_dble() \
__asm__ __volatile__ ("movapd %%xmm3, %%xmm8 \n\t" \
                      "movapd %%xmm4, %%xmm9 \n\t" \
                      "movapd %%xmm5, %%xmm10 \n\t" \
                      "mulpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm6, %%xmm4 \n\t" \
                      "mulpd %%xmm6, %%xmm5 \n\t" \
                      "shufpd $0x1, %%xmm8, %%xmm8 \n\t" \
                      "shufpd $0x1, %%xmm9, %%xmm9 \n\t" \
                      "shufpd $0x1, %%xmm10, %%xmm10 \n\t" \
                      "mulpd %%xmm7, %%xmm8 \n\t" \
                      "mulpd %%xmm7, %%xmm9 \n\t" \
                      "mulpd %%xmm7, %%xmm10 \n\t" \
                      "addpd %%xmm8, %%xmm3 \n\t" \
                      "addpd %%xmm9, %%xmm4 \n\t" \
                      "addpd %%xmm10, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm8", "xmm9", "xmm10")

/*
*  Computes s+z*r assuming s is stored in xmm0,xmm1,xmm2 and that z
*  has been loaded to xmm6,xmm7 by _sse_load_cmplx_dble(z). The result
*  appears in xmm0,xmm1,xmm2 and xmm3,xmm4,xmm5,xmm6,xmm7 are unchanged
*/

#define _sse_mulc_vector_add_dble(r) \
__asm__ __volatile__ ("movapd %0, %%xmm8 \n\t" \
                      "movapd %1, %%xmm9 \n\t" \
                      "movapd %2, %%xmm10 \n\t" \
                      "movapd %%xmm8, %%xmm11 \n\t" \
                      "movapd %%xmm9, %%xmm12 \n\t" \
                      "movapd %%xmm10, %%xmm13" \
                      : \
                      : \
                      "m" ((r).c1), \
                      "m" ((r).c2), \
                      "m" ((r).c3) \
                      : \
                      "xmm8", "xmm9", "xmm10", \
                      "xmm11", "xmm12", "xmm13"); \
__asm__ __volatile__ ("mulpd %%xmm6, %%xmm8 \n\t" \
                      "mulpd %%xmm6, %%xmm9 \n\t" \
                      "mulpd %%xmm6, %%xmm10 \n\t" \
                      "shufpd $0x1, %%xmm11, %%xmm11 \n\t" \
                      "shufpd $0x1, %%xmm12, %%xmm12 \n\t" \
                      "shufpd $0x1, %%xmm13, %%xmm13 \n\t" \
                      "addpd %%xmm8, %%xmm0 \n\t" \
                      "addpd %%xmm9, %%xmm1 \n\t" \
                      "addpd %%xmm10, %%xmm2 \n\t" \
                      "mulpd %%xmm7, %%xmm11 \n\t" \
                      "mulpd %%xmm7, %%xmm12 \n\t" \
                      "mulpd %%xmm7, %%xmm13 \n\t" \
                      "addpd %%xmm11, %%xmm0 \n\t" \
                      "addpd %%xmm12, %%xmm1 \n\t" \
                      "addpd %%xmm13, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm8", "xmm9", "xmm10", \
                      "xmm11", "xmm12", "xmm13")

/*
*  Computes s+z*r assuming s is stored in xmm3,xmm4,xmm5 and that z
*  has been loaded to xmm6,xmm7 by _sse_load_cmplx_dble(z). The result
*  appears in xmm4,xmm5,xmm6 and xmm0,xmm1,xmm2,xmm6,xmm7 are unchanged
*/

#define _sse_mulc_vector_add_up_dble(r) \
__asm__ __volatile__ ("movapd %0, %%xmm8 \n\t" \
                      "movapd %1, %%xmm9 \n\t" \
                      "movapd %2, %%xmm10 \n\t" \
                      "movapd %%xmm8, %%xmm11 \n\t" \
                      "movapd %%xmm9, %%xmm12 \n\t" \
                      "movapd %%xmm10, %%xmm13" \
                      : \
                      : \
                      "m" ((r).c1), \
                      "m" ((r).c2), \
                      "m" ((r).c3) \
                      : \
                      "xmm8", "xmm9", "xmm10", \
                      "xmm11", "xmm12", "xmm13"); \
__asm__ __volatile__ ("mulpd %%xmm6, %%xmm8 \n\t" \
                      "mulpd %%xmm6, %%xmm9 \n\t" \
                      "mulpd %%xmm6, %%xmm10 \n\t" \
                      "shufpd $0x1, %%xmm11, %%xmm11 \n\t" \
                      "shufpd $0x1, %%xmm12, %%xmm12 \n\t" \
                      "shufpd $0x1, %%xmm13, %%xmm13 \n\t" \
                      "addpd %%xmm8, %%xmm3 \n\t" \
                      "addpd %%xmm9, %%xmm4 \n\t" \
                      "addpd %%xmm10, %%xmm5 \n\t" \
                      "mulpd %%xmm7, %%xmm11 \n\t" \
                      "mulpd %%xmm7, %%xmm12 \n\t" \
                      "mulpd %%xmm7, %%xmm13 \n\t" \
                      "addpd %%xmm11, %%xmm3 \n\t" \
                      "addpd %%xmm12, %%xmm4 \n\t" \
                      "addpd %%xmm13, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm8", "xmm9", "xmm10", \
                      "xmm11", "xmm12", "xmm13")

/*
*  Loads (c,c) to xmm6 and xmm7
*/

#define _sse_load_real_dble(c) \
__asm__ __volatile__ ("movddup %0, %%xmm6 \n\t" \
                      "movddup %0, %%xmm7" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm6", "xmm7")

/*
*  Multiplies the complex numbers in xmm0,xmm1,xmm2 by c, assuming c has
*  been loaded to xmm6,xmm7 by _sse_load_real_dble(c). The result appears
*  in xmm0,xmm1,xmm2 all other xmm registers are unchanged
*/

#define _sse_mulr_vector_dble() \
__asm__ __volatile__ ("mulpd %%xmm6, %%xmm0 \n\t" \
                      "mulpd %%xmm7, %%xmm1 \n\t" \
                      "mulpd %%xmm6, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
*  Multiplies the complex numbers in xmm3,xmm4,xmm5 by c, assuming c has
*  been loaded to xmm6,xmm7 by _sse_load_real_dble(z). The result appears
*  in xmm3,xmm4,xmm5 all other xmm registers are unchanged
*/

#define _sse_mulr_vector_up_dble() \
__asm__ __volatile__ ("mulpd %%xmm7, %%xmm3 \n\t" \
                      "mulpd %%xmm6, %%xmm4 \n\t" \
                      "mulpd %%xmm7, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
*  Computes s+c*r assuming r is stored in xmm0,xmm1,xmm2 and that c
*  has been loaded to xmm6,xmm7 by _sse_load_real_dble(z). The result
*  appears in xmm0,xmm1,xmm2 all other xmm registers are unchanged
*/

#define _sse_mulr_vector_add_dble(s) \
__asm__ __volatile__ ("mulpd %%xmm6, %%xmm0 \n\t" \
                      "mulpd %%xmm7, %%xmm1 \n\t" \
                      "mulpd %%xmm6, %%xmm2 \n\t" \
                      "addpd %0, %%xmm0 \n\t" \
                      "addpd %1, %%xmm1 \n\t" \
                      "addpd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
*  Computes s+c*r assuming r is stored in xmm3,xmm4,xmm5 and that c
*  has been loaded to xmm6,xmm7 by _sse_load_real_dble(c). The result
*  appears in xmm4,xmm5,xmm6 and all other xmm registers are unchanged
*/

#define _sse_mulr_vector_add_up_dble(s) \
__asm__ __volatile__ ("mulpd %%xmm7, %%xmm3 \n\t" \
                      "mulpd %%xmm6, %%xmm4 \n\t" \
                      "mulpd %%xmm7, %%xmm5 \n\t" \
                      "addpd %0, %%xmm3 \n\t" \
                      "addpd %1, %%xmm4 \n\t" \
                      "addpd %2, %%xmm5" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/******************************************************************************
*
*  Action of su3 matrices on su3 vectors
*
******************************************************************************/

/*
* Multiplies an su3 vector s with an su3 matrix u, assuming s is
* stored in  xmm0,xmm1,xmm2.
*
* On output the result is in xmm3,xmm4,xmm5 and all registers except
* for xmm15 are changed.
*/

#define _sse_su3_multiply_dble(u) \
__asm__ __volatile__ ("movddup %0, %%xmm3 \n\t" \
                      "movddup %1, %%xmm6 \n\t" \
                      "movddup %2, %%xmm4 \n\t" \
                      "movddup %3, %%xmm7 \n\t" \
                      "movddup %4, %%xmm5 \n\t" \
                      "movddup %5, %%xmm8 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm1, %%xmm8 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "addpd %%xmm8, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c22.re), \
                      "m" ((u).c31.re), \
                      "m" ((u).c32.re) \
                      : \
                      "xmm0", "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("movddup %0, %%xmm9 \n\t" \
                      "movddup %1, %%xmm10 \n\t" \
                      "movddup %2, %%xmm11 \n\t" \
                      "movddup %3, %%xmm12 \n\t" \
                      "movddup %4, %%xmm13 \n\t" \
                      "movddup %5, %%xmm14 \n\t" \
                      "mulpd %%xmm2, %%xmm9 \n\t" \
                      "mulpd %%xmm0, %%xmm10 \n\t" \
                      "mulpd %%xmm2, %%xmm11 \n\t" \
                      "mulpd %%xmm0, %%xmm12 \n\t" \
                      "addpd %%xmm9, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm13 \n\t" \
                      "addsubpd %%xmm10, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm14 \n\t" \
                      "addpd %%xmm11, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c13.re), \
                      "m" ((u).c21.im), \
                      "m" ((u).c33.re), \
                      "m" ((u).c11.im), \
                      "m" ((u).c23.re), \
                      "m" ((u).c31.im) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm9", \
                      "xmm10", "xmm11", "xmm12", "xmm13", \
                      "xmm14"); \
__asm__ __volatile__ ("shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "addsubpd %%xmm12, %%xmm3 \n\t" \
                      "addpd %%xmm13, %%xmm4 \n\t" \
                      "addsubpd %%xmm14, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm1", "xmm2", "xmm3", "xmm4", \
                      "xmm5"); \
__asm__ __volatile__ ("movddup %0, %%xmm6 \n\t" \
                      "movddup %1, %%xmm7 \n\t" \
                      "movddup %2, %%xmm8 \n\t" \
                      "movddup %3, %%xmm9 \n\t" \
                      "movddup %4, %%xmm10 \n\t" \
                      "movddup %5, %%xmm11 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm8 \n\t" \
                      "mulpd %%xmm2, %%xmm9 \n\t" \
                      "addsubpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm1, %%xmm10 \n\t" \
                      "addsubpd %%xmm7, %%xmm4 \n\t" \
                      "mulpd %%xmm2, %%xmm11 \n\t" \
                      "addsubpd %%xmm8, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c12.im), \
                      "m" ((u).c23.im), \
                      "m" ((u).c32.im), \
                      "m" ((u).c13.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c33.im) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm6", \
                      "xmm7", "xmm8", "xmm9", "xmm10", \
                      "xmm11"); \
__asm__ __volatile__ ("addsubpd %%xmm9, %%xmm3 \n\t" \
                      "addsubpd %%xmm10, %%xmm4 \n\t" \
                      "addsubpd %%xmm11, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies an su3 vector s with an su3 matrix u^dagger, assuming s is
* stored in  xmm0,xmm1,xmm2.
*
* On output the result is in xmm3,xmm4,xmm5 and all registers except
* for xmm15 are changed.
*/

#define _sse_su3_inverse_multiply_dble(u) \
__asm__ __volatile__ ("movddup %0, %%xmm3 \n\t" \
                      "movddup %1, %%xmm6 \n\t" \
                      "movddup %2, %%xmm4 \n\t" \
                      "movddup %3, %%xmm7 \n\t" \
                      "movddup %4, %%xmm5 \n\t" \
                      "movddup %5, %%xmm8 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %6, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm8 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "addpd %%xmm8, %%xmm5 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c22.re), \
                      "m" ((u).c13.re), \
                      "m" ((u).c23.re), \
                      "m" (_sse_sgn1_dble) \
                      : \
                      "xmm0", "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("movddup %0, %%xmm9 \n\t" \
                      "movddup %1, %%xmm10 \n\t" \
                      "movddup %2, %%xmm11 \n\t" \
                      "movddup %3, %%xmm12 \n\t" \
                      "movddup %4, %%xmm13 \n\t" \
                      "movddup %5, %%xmm14 \n\t" \
                      "mulpd %%xmm2, %%xmm9 \n\t" \
                      "mulpd %6, %%xmm1 \n\t" \
                      "mulpd %%xmm0, %%xmm10 \n\t" \
                      "mulpd %%xmm2, %%xmm11 \n\t" \
                      "mulpd %%xmm0, %%xmm12 \n\t" \
                      "addpd %%xmm9, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm13 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "addpd %%xmm10, %%xmm4 \n\t" \
                      "mulpd %6, %%xmm2 \n\t" \
                      "addpd %%xmm11, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm14 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2" \
                      : \
                      : \
                      "m" ((u).c31.re), \
                      "m" ((u).c12.im), \
                      "m" ((u).c33.re), \
                      "m" ((u).c11.im), \
                      "m" ((u).c32.re), \
                      "m" ((u).c13.im), \
                      "m" (_sse_sgn1_dble) \
                      : \
                      "xmm1", "xmm2", "xmm3", "xmm4", \
                      "xmm5", "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("movddup %0, %%xmm6 \n\t" \
                      "movddup %1, %%xmm7 \n\t" \
                      "movddup %2, %%xmm8 \n\t" \
                      "movddup %3, %%xmm9 \n\t" \
                      "movddup %4, %%xmm10 \n\t" \
                      "movddup %5, %%xmm11 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "addpd %%xmm12, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm13, %%xmm4 \n\t" \
                      "mulpd %%xmm1, %%xmm8 \n\t" \
                      "addpd %%xmm14, %%xmm5 \n\t" \
                      "mulpd %%xmm2, %%xmm9 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm1, %%xmm10 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "mulpd %%xmm2, %%xmm11 \n\t" \
                      "addpd %%xmm8, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c21.im), \
                      "m" ((u).c32.im), \
                      "m" ((u).c23.im), \
                      "m" ((u).c31.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c33.im) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm6", \
                      "xmm7", "xmm8", "xmm9", "xmm10", \
                      "xmm11"); \
__asm__ __volatile__ ("addpd %%xmm9, %%xmm3 \n\t" \
                      "addpd %%xmm10, %%xmm4 \n\t" \
                      "addpd %%xmm11, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

#endif
