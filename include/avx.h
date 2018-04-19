
/*******************************************************************************
*
* File avx.h
*
* Copyright (C) 2013, 2016 Martin Luescher, Isabel Campos
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Macros for Dirac spinors, SU(3) vectors and SU(3) matrices using inline
* assembly AVX instructions. The machine is assumed to comply with the
* x86-64 instruction set.
*
*******************************************************************************/

#ifndef AVX_H
#define AVX_H

#ifndef SSE2_H
#include "sse2.h"
#endif

typedef struct __attribute__ ((aligned (32)))
{
   float c1,c2,c3,c4;
   float c5,c6,c7,c8;
} avx_float;

typedef struct __attribute__ ((aligned (32)))
{
   double c1,c2,c3,c4;
} avx_double;

static avx_double _avx_sgn12_dble __attribute__ ((unused)) ={-1.0,-1.0,1.0,1.0};
static avx_double _avx_sgn13_dble __attribute__ ((unused)) ={-1.0,1.0,-1.0,1.0};
static avx_double _avx_sgn14_dble __attribute__ ((unused)) ={-1.0,1.0,1.0,-1.0};
static avx_double _avx_sgn23_dble __attribute__ ((unused)) ={1.0,-1.0,-1.0,1.0};
static avx_double _avx_sgn24_dble __attribute__ ((unused)) ={1.0,-1.0,1.0,-1.0};
static avx_double _avx_sgn34_dble __attribute__ ((unused)) ={1.0,1.0,-1.0,-1.0};
static avx_double _avx_sgn_dble __attribute__ ((unused)) ={-1.0,-1.0,-1.0,-1.0};

static avx_float _avx_sgn_add __attribute__ ((unused))
={1.0f,1.0f,1.0f,1.0f,-1.0f,-1.0f,-1.0f,-1.0f};
static avx_float _avx_sgn_i_add __attribute__ ((unused))
={-1.0f,1.0f,-1.0f,1.0f,1.0f,-1.0f,1.0f,-1.0f};
static avx_float _avx_sgn_addsub __attribute__ ((unused))
={1.0f,1.0f,-1.0f,-1.0f,-1.0f,-1.0f,1.0f,1.0f};
static avx_float _avx_sgn_i_addsub __attribute__ ((unused))
={-1.0f,1.0f,1.0f,-1.0f,1.0f,-1.0f,-1.0f,1.0f};

#define _avx_zeroall() \
__asm__ __volatile__ ("vzeroall")

#define _avx_zeroupper() \
__asm__ __volatile__ ("vzeroupper")

/*******************************************************************************
*
* Macros operating on single precision data
*
*******************************************************************************/

/*******************************************************************************
*
* Macros for spinors in su3_vector order
*
*******************************************************************************/

/*
* Loads two spinors sl and sh to the low and high lanes of ymm0,..,ymm5. The
* ordering of the spinor components in the low lane is
*
* xmm0 <- sl.c1.c1,sl.c2.c1
* xmm1 <- sl.c1.c2,sl.c2.c2
* xmm2 <- sl.c1.c3,sl.c2.c3
* xmm3 <- sl.c3.c1,sl.c4.c1
* xmm4 <- sl.c3.c2,sl.c4.c2
* xmm5 <- sl.c3.c3,sl.c4.c3
*
* and those in the high lane are arranged in the same way. The registers
* ymm6,..,ymm11 are changed on exit.
*/

#define _avx_spinor_pair_load34(sl,sh) \
__asm__ __volatile__ ("vmovaps %0, %%xmm6 \n\t" \
                      "vmovaps %2, %%xmm7 \n\t" \
                      "vmovaps %4, %%xmm8" \
                      : \
                      : \
                      "m" ((sl).c1.c1), \
                      "m" ((sl).c1.c2), \
                      "m" ((sl).c1.c3), \
                      "m" ((sl).c2.c1), \
                      "m" ((sl).c2.c2), \
                      "m" ((sl).c2.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vmovaps %0, %%xmm9 \n\t" \
                      "vmovaps %2, %%xmm10 \n\t" \
                      "vmovaps %4, %%xmm11" \
                      : \
                      : \
                      "m" ((sl).c3.c1), \
                      "m" ((sl).c3.c2), \
                      "m" ((sl).c3.c3), \
                      "m" ((sl).c4.c1), \
                      "m" ((sl).c4.c2), \
                      "m" ((sl).c4.c3) \
                      : \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %2, %%ymm7, %%ymm7 \n\t" \
                      "vinsertf128 $0x1, %4, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((sh).c1.c1), \
                      "m" ((sh).c1.c2), \
                      "m" ((sh).c1.c3), \
                      "m" ((sh).c2.c1), \
                      "m" ((sh).c2.c2), \
                      "m" ((sh).c2.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm9, %%ymm9 \n\t" \
                      "vinsertf128 $0x1, %2, %%ymm10, %%ymm10 \n\t" \
                      "vinsertf128 $0x1, %4, %%ymm11, %%ymm11" \
                      : \
                      : \
                      "m" ((sh).c3.c1), \
                      "m" ((sh).c3.c2), \
                      "m" ((sh).c3.c3), \
                      "m" ((sh).c4.c1), \
                      "m" ((sh).c4.c2), \
                      "m" ((sh).c4.c3) \
                      : \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vshufps $0xe4, %%ymm7, %%ymm6, %%ymm0 \n\t" \
                      "vshufps $0xe4, %%ymm10, %%ymm9, %%ymm3 \n\t" \
                      "vshufps $0x4e, %%ymm8, %%ymm6, %%ymm1 \n\t" \
                      "vshufps $0x4e, %%ymm11, %%ymm9, %%ymm4 \n\t" \
                      "vshufps $0xe4, %%ymm8, %%ymm7, %%ymm2 \n\t" \
                      "vshufps $0xe4, %%ymm11, %%ymm10, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Loads two spinors sl and sh to the low and high lanes of ymm0,..,ymm5. The
* ordering of the spinor components in the low lane is
*
* xmm0 <- sl.c1.c1,sl.c2.c1
* xmm1 <- sl.c1.c2,sl.c2.c2
* xmm2 <- sl.c1.c3,sl.c2.c3
* xmm3 <- sl.c4.c1,sl.c3.c1       (note: unusual order)
* xmm4 <- sl.c4.c2,sl.c3.c2
* xmm5 <- sl.c4.c3,sl.c3.c3
*
* and those in the high lane are arranged in the same way. The registers
* ymm6,..,ymm11 are changed on exit.
*/

#define _avx_spinor_pair_load43(sl,sh) \
__asm__ __volatile__ ("vmovaps %0, %%xmm6 \n\t" \
                      "vmovaps %2, %%xmm7 \n\t" \
                      "vmovaps %4, %%xmm8" \
                      : \
                      : \
                      "m" ((sl).c1.c1), \
                      "m" ((sl).c1.c2), \
                      "m" ((sl).c1.c3), \
                      "m" ((sl).c2.c1), \
                      "m" ((sl).c2.c2), \
                      "m" ((sl).c2.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vmovaps %0, %%xmm9 \n\t" \
                      "vmovaps %2, %%xmm10 \n\t" \
                      "vmovaps %4, %%xmm11" \
                      : \
                      : \
                      "m" ((sl).c3.c1), \
                      "m" ((sl).c3.c2), \
                      "m" ((sl).c3.c3), \
                      "m" ((sl).c4.c1), \
                      "m" ((sl).c4.c2), \
                      "m" ((sl).c4.c3) \
                      : \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %2, %%ymm7, %%ymm7 \n\t" \
                      "vinsertf128 $0x1, %4, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((sh).c1.c1), \
                      "m" ((sh).c1.c2), \
                      "m" ((sh).c1.c3), \
                      "m" ((sh).c2.c1), \
                      "m" ((sh).c2.c2), \
                      "m" ((sh).c2.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm9, %%ymm9 \n\t" \
                      "vinsertf128 $0x1, %2, %%ymm10, %%ymm10 \n\t" \
                      "vinsertf128 $0x1, %4, %%ymm11, %%ymm11" \
                      : \
                      : \
                      "m" ((sh).c3.c1), \
                      "m" ((sh).c3.c2), \
                      "m" ((sh).c3.c3), \
                      "m" ((sh).c4.c1), \
                      "m" ((sh).c4.c2), \
                      "m" ((sh).c4.c3) \
                      : \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vshufps $0xe4, %%ymm7, %%ymm6, %%ymm0 \n\t" \
                      "vshufps $0x4e, %%ymm9, %%ymm10, %%ymm3 \n\t" \
                      "vshufps $0x4e, %%ymm8, %%ymm6, %%ymm1 \n\t" \
                      "vshufps $0xe4, %%ymm9, %%ymm11, %%ymm4 \n\t" \
                      "vshufps $0xe4, %%ymm8, %%ymm7, %%ymm2 \n\t" \
                      "vshufps $0x4e, %%ymm10, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Loads the spinor s to xmm0,..,xmm5 in the order
*
* xmm0 <- s.c1.c1,s.c2.c1
* xmm1 <- s.c1.c2,s.c2.c2
* xmm2 <- s.c1.c3,s.c2.c3
* xmm3 <- s.c3.c1,s.c4.c1
* xmm4 <- s.c3.c2,s.c4.c2
* xmm5 <- s.c3.c3,s.c4.c3
*
* and duplicates these values to the upper lanes of ymm0,..ymm5. The registers
* ymm6,..,ymm11 are changed on exit.
*/

#define _avx_spinor_load_dup(s) \
__asm__ __volatile__ ("vbroadcastf128 %0, %%ymm6 \n\t" \
                      "vbroadcastf128 %2, %%ymm7 \n\t" \
                      "vbroadcastf128 %4, %%ymm8" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vbroadcastf128 %0, %%ymm9 \n\t" \
                      "vbroadcastf128 %2, %%ymm10 \n\t" \
                      "vbroadcastf128 %4, %%ymm11" \
                      : \
                      : \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2), \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vshufps $0xe4, %%ymm7, %%ymm6, %%ymm0 \n\t" \
                      "vshufps $0xe4, %%ymm10, %%ymm9, %%ymm3 \n\t" \
                      "vshufps $0x4e, %%ymm8, %%ymm6, %%ymm1 \n\t" \
                      "vshufps $0x4e, %%ymm11, %%ymm9, %%ymm4 \n\t" \
                      "vshufps $0xe4, %%ymm8, %%ymm7, %%ymm2 \n\t" \
                      "vshufps $0xe4, %%ymm11, %%ymm10, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Stores the low and high lanes of ymm0,..,ymm5 to the spinors rl and rh,
* assuming the spinor components are ordered as if they were loaded with
* _avx_spinor_pair_load34(rl,rh). The registers ymm6,..,ymm11 are changed
* on exit.
*/

#define _avx_spinor_pair_store34(rl,rh) \
__asm__ __volatile__ ("vshufps $0x44, %%ymm1, %%ymm0, %%ymm6 \n\t" \
                      "vshufps $0x44, %%ymm4, %%ymm3, %%ymm9 \n\t" \
                      "vshufps $0xe4, %%ymm0, %%ymm2, %%ymm7 \n\t" \
                      "vshufps $0xe4, %%ymm3, %%ymm5, %%ymm10 \n\t" \
                      "vshufps $0xee, %%ymm2, %%ymm1, %%ymm8 \n\t" \
                      "vshufps $0xee, %%ymm5, %%ymm4, %%ymm11" \
                      : \
                      : \
                      : \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vmovaps %%xmm6, %0 \n\t" \
                      "vmovaps %%xmm7, %2 \n\t" \
                      "vmovaps %%xmm8, %4" \
                      : \
                      "=m" ((rl).c1.c1), \
                      "=m" ((rl).c1.c2), \
                      "=m" ((rl).c1.c3), \
                      "=m" ((rl).c2.c1), \
                      "=m" ((rl).c2.c2), \
                      "=m" ((rl).c2.c3)); \
__asm__ __volatile__ ("vmovaps %%xmm9, %0 \n\t" \
                      "vmovaps %%xmm10, %2 \n\t" \
                      "vmovaps %%xmm11, %4" \
                      : \
                      "=m" ((rl).c3.c1), \
                      "=m" ((rl).c3.c2), \
                      "=m" ((rl).c3.c3), \
                      "=m" ((rl).c4.c1), \
                      "=m" ((rl).c4.c2), \
                      "=m" ((rl).c4.c3)); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm6, %0 \n\t" \
                      "vextractf128 $0x1, %%ymm7, %2 \n\t" \
                      "vextractf128 $0x1, %%ymm8, %4" \
                      : \
                      "=m" ((rh).c1.c1), \
                      "=m" ((rh).c1.c2), \
                      "=m" ((rh).c1.c3), \
                      "=m" ((rh).c2.c1), \
                      "=m" ((rh).c2.c2), \
                      "=m" ((rh).c2.c3)); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm9, %0 \n\t" \
                      "vextractf128 $0x1, %%ymm10, %2 \n\t" \
                      "vextractf128 $0x1, %%ymm11, %4" \
                      : \
                      "=m" ((rh).c3.c1), \
                      "=m" ((rh).c3.c2), \
                      "=m" ((rh).c3.c3), \
                      "=m" ((rh).c4.c1), \
                      "=m" ((rh).c4.c2), \
                      "=m" ((rh).c4.c3))

/*
* Stores the low and high lanes of ymm0,..,ymm5 to the spinors rl and rh,
* assuming the spinor components are ordered as if they were loaded with
* _avx_spinor_pair_load43(rl,rh). The registers ymm6,..,ymm11 are changed
* on exit.
*/

#define _avx_spinor_pair_store43(rl,rh) \
__asm__ __volatile__ ("vshufps $0x44, %%ymm1, %%ymm0, %%ymm6 \n\t" \
                      "vshufps $0xee, %%ymm4, %%ymm3, %%ymm9 \n\t" \
                      "vshufps $0xe4, %%ymm0, %%ymm2, %%ymm7 \n\t" \
                      "vshufps $0x4e, %%ymm3, %%ymm5, %%ymm10 \n\t" \
                      "vshufps $0xee, %%ymm2, %%ymm1, %%ymm8 \n\t" \
                      "vshufps $0x44, %%ymm5, %%ymm4, %%ymm11" \
                      : \
                      : \
                      : \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vmovaps %%xmm6, %0 \n\t" \
                      "vmovaps %%xmm7, %2 \n\t" \
                      "vmovaps %%xmm8, %4" \
                      : \
                      "=m" ((rl).c1.c1), \
                      "=m" ((rl).c1.c2), \
                      "=m" ((rl).c1.c3), \
                      "=m" ((rl).c2.c1), \
                      "=m" ((rl).c2.c2), \
                      "=m" ((rl).c2.c3)); \
__asm__ __volatile__ ("vmovaps %%xmm9, %0 \n\t" \
                      "vmovaps %%xmm10, %2 \n\t" \
                      "vmovaps %%xmm11, %4" \
                      : \
                      "=m" ((rl).c3.c1), \
                      "=m" ((rl).c3.c2), \
                      "=m" ((rl).c3.c3), \
                      "=m" ((rl).c4.c1), \
                      "=m" ((rl).c4.c2), \
                      "=m" ((rl).c4.c3)); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm6, %0 \n\t" \
                      "vextractf128 $0x1, %%ymm7, %2 \n\t" \
                      "vextractf128 $0x1, %%ymm8, %4" \
                      : \
                      "=m" ((rh).c1.c1), \
                      "=m" ((rh).c1.c2), \
                      "=m" ((rh).c1.c3), \
                      "=m" ((rh).c2.c1), \
                      "=m" ((rh).c2.c2), \
                      "=m" ((rh).c2.c3)); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm9, %0 \n\t" \
                      "vextractf128 $0x1, %%ymm10, %2 \n\t" \
                      "vextractf128 $0x1, %%ymm11, %4" \
                      : \
                      "=m" ((rh).c3.c1), \
                      "=m" ((rh).c3.c2), \
                      "=m" ((rh).c3.c3), \
                      "=m" ((rh).c4.c1), \
                      "=m" ((rh).c4.c2), \
                      "=m" ((rh).c4.c3))

/*
* Loads the lower Weyl spinors of the Dirac spinors sl and sh to the low and
* high lanes of ymm0,..,ymm3. The ordering of the spinor components in the
* low lane is
*
* xmm0 <- sl.c1.c1,sl.c2.c1
* xmm1 <- sl.c1.c2,sl.c2.c2
* xmm2 <- sl.c1.c3,sl.c2.c3
*
* and those in the high lane are arranged in the same way. The registers
* ymm6,..,ymm8 are changed on exit. Also applies if sl and sh are Weyl
* spinors.
*/

#define _avx_weyl_pair_load12(sl,sh) \
__asm__ __volatile__ ("vmovaps %0, %%xmm6 \n\t" \
                      "vmovaps %2, %%xmm7 \n\t" \
                      "vmovaps %4, %%xmm8" \
                      : \
                      : \
                      "m" ((sl).c1.c1), \
                      "m" ((sl).c1.c2), \
                      "m" ((sl).c1.c3), \
                      "m" ((sl).c2.c1), \
                      "m" ((sl).c2.c2), \
                      "m" ((sl).c2.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %2, %%ymm7, %%ymm7 \n\t" \
                      "vinsertf128 $0x1, %4, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((sh).c1.c1), \
                      "m" ((sh).c1.c2), \
                      "m" ((sh).c1.c3), \
                      "m" ((sh).c2.c1), \
                      "m" ((sh).c2.c2), \
                      "m" ((sh).c2.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vshufps $0xe4, %%ymm7, %%ymm6, %%ymm0 \n\t" \
                      "vshufps $0x4e, %%ymm8, %%ymm6, %%ymm1 \n\t" \
                      "vshufps $0xe4, %%ymm8, %%ymm7, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Loads the upper Weyl spinors of the Dirac spinors sl and sh to the low and
* high lanes of ymm0,..,ymm3. The ordering of the spinor components in the
* low lane is
*
* xmm0 <- sl.c3.c1,sl.c4.c1
* xmm1 <- sl.c3.c2,sl.c4.c2
* xmm2 <- sl.c3.c3,sl.c4.c3
*
* and those in the high lane are arranged in the same way. The registers
* ymm6,..,ymm8 are changed on exit.
*/

#define _avx_weyl_pair_load34(sl,sh) \
__asm__ __volatile__ ("vmovaps %0, %%xmm6 \n\t" \
                      "vmovaps %2, %%xmm7 \n\t" \
                      "vmovaps %4, %%xmm8" \
                      : \
                      : \
                      "m" ((sl).c3.c1), \
                      "m" ((sl).c3.c2), \
                      "m" ((sl).c3.c3), \
                      "m" ((sl).c4.c1), \
                      "m" ((sl).c4.c2), \
                      "m" ((sl).c4.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %2, %%ymm7, %%ymm7 \n\t" \
                      "vinsertf128 $0x1, %4, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((sh).c3.c1), \
                      "m" ((sh).c3.c2), \
                      "m" ((sh).c3.c3), \
                      "m" ((sh).c4.c1), \
                      "m" ((sh).c4.c2), \
                      "m" ((sh).c4.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vshufps $0xe4, %%ymm7, %%ymm6, %%ymm0 \n\t" \
                      "vshufps $0x4e, %%ymm8, %%ymm6, %%ymm1 \n\t" \
                      "vshufps $0xe4, %%ymm8, %%ymm7, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Stores the low and high lanes of ymm0,..,ymm3 to the lower Weyl spinors
* of the Dirac spinors rl and rh, assuming the spinor components are ordered
* as if they were loaded with _avx_weyl_pair_load12(rl,rh). The registers
* ymm6,..,ymm8 are changed on exit. Also applies if rl and rh are Weyl
* spinors.
*/

#define _avx_weyl_pair_store12(rl,rh) \
__asm__ __volatile__ ("vshufps $0x44, %%ymm1, %%ymm0, %%ymm6 \n\t" \
                      "vshufps $0xe4, %%ymm0, %%ymm2, %%ymm7 \n\t" \
                      "vshufps $0xee, %%ymm2, %%ymm1, %%ymm8" \
                      : \
                      : \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vmovaps %%xmm6, %0 \n\t" \
                      "vmovaps %%xmm7, %2 \n\t" \
                      "vmovaps %%xmm8, %4" \
                      : \
                      "=m" ((rl).c1.c1), \
                      "=m" ((rl).c1.c2), \
                      "=m" ((rl).c1.c3), \
                      "=m" ((rl).c2.c1), \
                      "=m" ((rl).c2.c2), \
                      "=m" ((rl).c2.c3)); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm6, %0 \n\t" \
                      "vextractf128 $0x1, %%ymm7, %2 \n\t" \
                      "vextractf128 $0x1, %%ymm8, %4" \
                      : \
                      "=m" ((rh).c1.c1), \
                      "=m" ((rh).c1.c2), \
                      "=m" ((rh).c1.c3), \
                      "=m" ((rh).c2.c1), \
                      "=m" ((rh).c2.c2), \
                      "=m" ((rh).c2.c3))

/*
* Stores the low and high lanes of ymm0,..,ymm3 to the upper Weyl spinors
* of the Dirac spinors rl and rh, assuming the spinor components are ordered
* as if they were loaded with _avx_weyl_pair_load34(rl,rh). The registers
* ymm6,..,ymm8 are changed on exit.
*/

#define _avx_weyl_pair_store34(rl,rh) \
__asm__ __volatile__ ("vshufps $0x44, %%ymm1, %%ymm0, %%ymm6 \n\t" \
                      "vshufps $0xe4, %%ymm0, %%ymm2, %%ymm7 \n\t" \
                      "vshufps $0xee, %%ymm2, %%ymm1, %%ymm8" \
                      : \
                      : \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vmovaps %%xmm6, %0 \n\t" \
                      "vmovaps %%xmm7, %2 \n\t" \
                      "vmovaps %%xmm8, %4" \
                      : \
                      "=m" ((rl).c3.c1), \
                      "=m" ((rl).c3.c2), \
                      "=m" ((rl).c3.c3), \
                      "=m" ((rl).c4.c1), \
                      "=m" ((rl).c4.c2), \
                      "=m" ((rl).c4.c3)); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm6, %0 \n\t" \
                      "vextractf128 $0x1, %%ymm7, %2 \n\t" \
                      "vextractf128 $0x1, %%ymm8, %4" \
                      : \
                      "=m" ((rh).c3.c1), \
                      "=m" ((rh).c3.c2), \
                      "=m" ((rh).c3.c3), \
                      "=m" ((rh).c4.c1), \
                      "=m" ((rh).c4.c2), \
                      "=m" ((rh).c4.c3))

/*
* Splits the registers ymm3,..,ymm5 according to
*
*  xmm3 <- ymm3_lo + ymm3_hi
*  xmm4 <- ymm4_lo + ymm4_hi
*  xmm5 <- ymm5_lo + ymm5_hi
*
*  xmm6 <- ymm3_lo - ymm3_hi
*  xmm7 <- ymm4_lo - ymm4_hi
*  xmm8 <- ymm5_lo - ymm5_hi
*
* where *_lo and *_hi are the low and high lanes of the registers. The
* registers ymm9,..,ymm11 are used as workspace.
*/

#define _avx_spinor_split() \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm3, %%xmm9 \n\t" \
                      "vextractf128 $0x1, %%ymm4, %%xmm10 \n\t" \
                      "vextractf128 $0x1, %%ymm5, %%xmm11 \n\t" \
                      "vsubps %%xmm9, %%xmm3, %%xmm6 \n\t" \
                      "vsubps %%xmm10, %%xmm4, %%xmm7 \n\t" \
                      "vsubps %%xmm11, %%xmm5, %%xmm8 \n\t" \
                      "vaddps %%xmm9, %%xmm3, %%xmm3 \n\t" \
                      "vaddps %%xmm10, %%xmm4, %%xmm4 \n\t" \
                      "vaddps %%xmm11, %%xmm5, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11")

/*
* Moves the lower lanes of ymm6,..,ymm8 to the upper lanes of ymm3,..,ymm5.
*/

#define _avx_spinor_unsplit() \
__asm__ __volatile__ ("vinsertf128 $0x1, %%xmm6, %%ymm3, %%ymm3 \n\t" \
                      "vinsertf128 $0x1, %%xmm7, %%ymm4, %%ymm4 \n\t" \
                      "vinsertf128 $0x1, %%xmm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies ymm3,..,ymm5 by the avx_float c. The register ymm15 is used as
* workspace.
*/

#define _avx_spinor_mul_up(c) \
__asm__ __volatile__ ("vmovaps %0, %%ymm15 \n\t" \
                      "vmulps %%ymm15, %%ymm3, %%ymm3 \n\t" \
                      "vmulps %%ymm15, %%ymm4, %%ymm4 \n\t" \
                      "vmulps %%ymm15, %%ymm5, %%ymm5" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm15")

/*
* Exchanges real and imaginary parts of the double words in ymm3,..,ymm5
* and multiplies these registers by the avx_float c. The register ymm15 is
* used as workspace.
*/

#define _avx_spinor_imul_up(c) \
__asm__ __volatile__ ("vmovaps %0, %%ymm15 \n\t" \
                      "vpermilps $0xb1, %%ymm3, %%ymm3 \n\t" \
                      "vpermilps $0xb1, %%ymm4, %%ymm4 \n\t" \
                      "vpermilps $0xb1, %%ymm5, %%ymm5 \n\t" \
                      "vmulps %%ymm15, %%ymm3, %%ymm3 \n\t" \
                      "vmulps %%ymm15, %%ymm4, %%ymm4 \n\t" \
                      "vmulps %%ymm15, %%ymm5, %%ymm5" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm15")

/*
* Exchanges the high and low words in the two lanes of ymm3,..,ymm5.
*/

#define _avx_spinor_xch_up() \
__asm__ __volatile__ ("vpermilps $0x4e, %%ymm3, %%ymm3 \n\t" \
                      "vpermilps $0x4e, %%ymm4, %%ymm4 \n\t" \
                      "vpermilps $0x4e, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Exchanges the high and low words in the two lanes of ymm3,..,ymm5, then the
* real and imaginary parts of the words and finally multiplies the registers
* by the avx_float c. The register ymm15 is used as workspace.
*/

#define _avx_spinor_xch_imul_up(c) \
__asm__ __volatile__ ("vmovaps %0, %%ymm15 \n\t" \
                      "vpermilps $0x1b, %%ymm3, %%ymm3 \n\t" \
                      "vpermilps $0x1b, %%ymm4, %%ymm4 \n\t" \
                      "vpermilps $0x1b, %%ymm5, %%ymm5 \n\t" \
                      "vmulps %%ymm15, %%ymm3, %%ymm3 \n\t" \
                      "vmulps %%ymm15, %%ymm4, %%ymm4 \n\t" \
                      "vmulps %%ymm15, %%ymm5, %%ymm5" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm15")

/*
* Multiplies xmm6,..,xmm8 by the sse_float c. The register ymm15 is used as
* workspace.
*/

#define _avx_weyl_mul(c) \
__asm__ __volatile__ ("vmovaps %0, %%xmm15 \n\t" \
                      "vmulps %%xmm15, %%xmm6, %%xmm6 \n\t" \
                      "vmulps %%xmm15, %%xmm7, %%xmm7 \n\t" \
                      "vmulps %%xmm15, %%xmm8, %%xmm8" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm6", "xmm7", "xmm8", "xmm15")

/*
* Exchanges real and imaginary parts of the double words in xmm6,..,xmm8
* and multiplies these registers by the sse_float c. The register ymm15 is
* used as workspace.
*/

#define _avx_weyl_imul(c) \
__asm__ __volatile__ ("vmovaps %0, %%xmm15 \n\t" \
                      "vpermilps $0xb1, %%xmm6, %%xmm6 \n\t" \
                      "vpermilps $0xb1, %%xmm7, %%xmm7 \n\t" \
                      "vpermilps $0xb1, %%xmm8, %%xmm8 \n\t" \
                      "vmulps %%xmm15, %%xmm6, %%xmm6 \n\t" \
                      "vmulps %%xmm15, %%xmm7, %%xmm7 \n\t" \
                      "vmulps %%xmm15, %%xmm8, %%xmm8" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm6", "xmm7", "xmm8", "xmm15")

/*
* Exchanges the high and low words of xmm6,..,xmm8.
*/

#define _avx_weyl_xch() \
__asm__ __volatile__ ("vpermilps $0x4e, %%xmm6, %%xmm6 \n\t" \
                      "vpermilps $0x4e, %%xmm7, %%xmm7 \n\t" \
                      "vpermilps $0x4e, %%xmm8, %%xmm8" \
                      : \
                      : \
                      : \
                      "xmm6", "xmm7", "xmm8")

/*
* Exchanges the high and low words of xmm6,..,xmm8, then the real and
* imaginary parts of the words and finally multiplies the registers by
* the sse_float c. The register ymm15 is used as workspace.
*/

#define _avx_weyl_xch_imul(c) \
__asm__ __volatile__ ("vmovaps %0, %%xmm15 \n\t" \
                      "vpermilps $0x1b, %%xmm6, %%xmm6 \n\t" \
                      "vpermilps $0x1b, %%xmm7, %%xmm7 \n\t" \
                      "vpermilps $0x1b, %%xmm8, %%xmm8 \n\t" \
                      "vmulps %%xmm15, %%xmm6, %%xmm6 \n\t" \
                      "vmulps %%xmm15, %%xmm7, %%xmm7 \n\t" \
                      "vmulps %%xmm15, %%xmm8, %%xmm8" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm6", "xmm7", "xmm8", "xmm15")

/*
* Adds ymm3,..,ymm5 to ymm0,..,ymm2
*/

#define _avx_spinor_add() \
__asm__ __volatile__ ("vaddps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddps %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Subtracts ymm3,..,ymm5 from ymm0,..,ymm2
*/

#define _avx_spinor_sub() \
__asm__ __volatile__ ("vsubps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vsubps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vsubps %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Adds (subtracts) the low (high) words in the two lanes of ymm3,..,ymm5
* to (from) ymm0,..,ymm2. The registers ymm6,ymm7,ymm8 are changed on exit.
*/

#define _avx_spinor_addsub() \
__asm__ __volatile__ ("vaddps %%ymm3, %%ymm0, %%ymm6 \n\t" \
                      "vaddps %%ymm4, %%ymm1, %%ymm7 \n\t" \
                      "vaddps %%ymm5, %%ymm2, %%ymm8 \n\t" \
                      "vsubps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vsubps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vsubps %%ymm5, %%ymm2, %%ymm2 \n\t" \
                      "vblendps $0x33, %%ymm6, %%ymm0, %%ymm0 \n\t" \
                      "vblendps $0x33, %%ymm7, %%ymm1, %%ymm1 \n\t" \
                      "vblendps $0x33, %%ymm8, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm6", "xmm7", "xmm8")

/*
* Adds (subtracts) the high (low) words in the two lanes of ymm3,..,ymm5
* to (from) ymm0,..,ymm2. The registers ymm6,..,ymm8 are changed on exit.
*/

#define _avx_spinor_subadd() \
__asm__ __volatile__ ("vaddps %%ymm3, %%ymm0, %%ymm6 \n\t" \
                      "vaddps %%ymm4, %%ymm1, %%ymm7 \n\t" \
                      "vaddps %%ymm5, %%ymm2, %%ymm8 \n\t" \
                      "vsubps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vsubps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vsubps %%ymm5, %%ymm2, %%ymm2 \n\t" \
                      "vblendps $0xcc, %%ymm6, %%ymm0, %%ymm0 \n\t" \
                      "vblendps $0xcc, %%ymm7, %%ymm1, %%ymm1 \n\t" \
                      "vblendps $0xcc, %%ymm8, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm6", "xmm7", "xmm8")

/*
* Multiplies ymm3,..,ymm5 with i and adds them to ymm0,..,ymm2. The
* registers ymm3,..,ymm5 are changed on exit.
*/

#define _avx_spinor_i_add() \
__asm__ __volatile__ ("vpermilps $0xb1, %%ymm3, %%ymm3 \n\t" \
                      "vpermilps $0xb1, %%ymm4, %%ymm4 \n\t" \
                      "vpermilps $0xb1, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddsubps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddsubps %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies ymm3,..,ymm5 with i and subtracts them from ymm0,..,ymm2.
*/

#define _avx_spinor_i_sub() \
__asm__ __volatile__ ("vpermilps $0xb1, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0xb1, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xb1, %%ymm2, %%ymm2 \n\t" \
                      "vaddsubps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddsubps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddsubps %%ymm5, %%ymm2, %%ymm2 \n\t" \
                      "vpermilps $0xb1, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0xb1, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xb1, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Exchanges the high and low words of ymm3,..,ymm5, multiplies them with i
* and adds the result to ymm0,..,ymm2. The registers ymm3,..,ymm5 are
* changed on exit.
*/

#define _avx_spinor_xch_i_add() \
__asm__ __volatile__ ("vpermilps $0x1b, %%ymm3, %%ymm3 \n\t" \
                      "vpermilps $0x1b, %%ymm4, %%ymm4 \n\t" \
                      "vpermilps $0x1b, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddsubps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddsubps %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Exchanges the high and low words of ymm3,..,ymm5, multiplies them with i
* and subtracts the result from ymm0,..,ymm2.
*/

#define _avx_spinor_xch_i_sub() \
__asm__ __volatile__ ("vpermilps $0x1b, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0x1b, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0x1b, %%ymm2, %%ymm2 \n\t" \
                      "vaddsubps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddsubps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddsubps %%ymm5, %%ymm2, %%ymm2 \n\t" \
                      "vpermilps $0x1b, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0x1b, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0x1b, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Multiplies the low and high words in the two lanes of ymm3,..,ymm5 with
* i and -i respectively and adds these registers to ymm0,..,ymm2. The
* registers ymm3,..,ymm5 are changed on exit.
*/

#define _avx_spinor_i_addsub() \
__asm__ __volatile__ ("vpermilps $0xb4, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0xb4, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xb4, %%ymm2, %%ymm2 \n\t" \
                      "vpermilps $0xe1, %%ymm3, %%ymm3 \n\t" \
                      "vpermilps $0xe1, %%ymm4, %%ymm4 \n\t" \
                      "vpermilps $0xe1, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddsubps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddsubps %%ymm5, %%ymm2, %%ymm2 \n\t" \
                      "vpermilps $0xb4, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0xb4, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xb4, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies the low and high words in the two lanes of ymm3,..,ymm5 with
* -i and i respectively and adds these registers to ymm0,..,ymm2. The
* registers ymm3,..,ymm5 are changed on exit.
*/

#define _avx_spinor_i_subadd() \
__asm__ __volatile__ ("vpermilps $0xe1, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0xe1, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xe1, %%ymm2, %%ymm2 \n\t" \
                      "vpermilps $0xb4, %%ymm3, %%ymm3 \n\t" \
                      "vpermilps $0xb4, %%ymm4, %%ymm4 \n\t" \
                      "vpermilps $0xb4, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddsubps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddsubps %%ymm5, %%ymm2, %%ymm2 \n\t" \
                      "vpermilps $0xe1, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0xe1, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xe1, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Exchanges the high and low words in each lane of ymm3,..,ymm5.
*/

#define _avx_spinor_xch() \
__asm__ __volatile__ ("vpermilps $0x4e, %%ymm3, %%ymm3 \n\t" \
                      "vpermilps $0x4e, %%ymm4, %%ymm4 \n\t" \
                      "vpermilps $0x4e, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

/******************************************************************************
*
*  Action of su3 matrices on su3 vectors
*
******************************************************************************/

/*
* Multiplies pairs of su3 vectors, stored in the low and high lanes of
* ymm0,..,ymm2, with su3 matrices ul and uh, respectively. The vectors
* are assumed to be in vertical order and the products are returned in the
* same order in the registers ymm3,..,ymm5. All registers except for
* ymm15 are changed on exit.
*/

#if (defined FMA3)

#define _avx_su3_pair_multiply(ul,uh) \
__asm__ __volatile__ ("vpermilps $0xf5, %%ymm0, %%ymm6 \n\t" \
                      "vpermilps $0xf5, %%ymm1, %%ymm7 \n\t" \
                      "vpermilps $0xf5, %%ymm2, %%ymm8 \n\t" \
                      "vpermilps $0xa0, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0xa0, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xa0, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm12 \n\t" \
                      "vbroadcastsd %4, %%ymm13 \n\t" \
                      "vbroadcastsd %5, %%ymm14 \n\t" \
                      "vblendps $0xf0, %%ymm12, %%ymm9, %%ymm9 \n\t" \
                      "vblendps $0xf0, %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vblendps $0xf0, %%ymm14, %%ymm11, %%ymm11 \n\t" \
                      "vpermilps $0xb1, %%ymm9, %%ymm12 \n\t" \
                      "vpermilps $0xb1, %%ymm10, %%ymm13 \n\t" \
                      "vpermilps $0xb1, %%ymm11, %%ymm14" \
                      : \
                      : \
                      "m" ((ul).c11), \
                      "m" ((ul).c22), \
                      "m" ((ul).c33), \
                      "m" ((uh).c11), \
                      "m" ((uh).c22), \
                      "m" ((uh).c33) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm6, %%ymm12, %%ymm3 \n\t" \
                      "vmulps %%ymm7, %%ymm13, %%ymm4 \n\t" \
                      "vmulps %%ymm8, %%ymm14, %%ymm5 \n\t" \
                      "vfmaddsub231ps %%ymm0, %%ymm9, %%ymm3 \n\t" \
                      "vfmaddsub231ps %%ymm1, %%ymm10, %%ymm4 \n\t" \
                      "vfmaddsub231ps %%ymm2, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm12 \n\t" \
                      "vbroadcastsd %1, %%ymm13 \n\t" \
                      "vbroadcastsd %2, %%ymm14 \n\t" \
                      "vbroadcastsd %3, %%ymm9 \n\t" \
                      "vbroadcastsd %4, %%ymm10 \n\t" \
                      "vbroadcastsd %5, %%ymm11 \n\t" \
                      "vblendps $0xf0, %%ymm9, %%ymm12, %%ymm12 \n\t" \
                      "vblendps $0xf0, %%ymm10, %%ymm13, %%ymm13 \n\t" \
                      "vblendps $0xf0, %%ymm11, %%ymm14, %%ymm14 \n\t" \
                      "vpermilps $0xb1, %%ymm12, %%ymm9 \n\t" \
                      "vpermilps $0xb1, %%ymm13, %%ymm10 \n\t" \
                      "vpermilps $0xb1, %%ymm14, %%ymm11" \
                      : \
                      : \
                      "m" ((ul).c12), \
                      "m" ((ul).c23), \
                      "m" ((ul).c31), \
                      "m" ((uh).c12), \
                      "m" ((uh).c23), \
                      "m" ((uh).c31) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vfmaddsub231ps %%ymm7, %%ymm9, %%ymm3 \n\t" \
                      "vfmaddsub231ps %%ymm8, %%ymm10, %%ymm4 \n\t" \
                      "vfmaddsub231ps %%ymm6, %%ymm11, %%ymm5 \n\t" \
                      "vfmaddsub231ps %%ymm1, %%ymm12, %%ymm3 \n\t" \
                      "vfmaddsub231ps %%ymm2, %%ymm13, %%ymm4 \n\t" \
                      "vfmaddsub231ps %%ymm0, %%ymm14, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm12 \n\t" \
                      "vbroadcastsd %4, %%ymm13 \n\t" \
                      "vbroadcastsd %5, %%ymm14 \n\t" \
                      "vblendps $0xf0, %%ymm12, %%ymm9, %%ymm9 \n\t" \
                      "vblendps $0xf0, %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vblendps $0xf0, %%ymm14, %%ymm11, %%ymm11 \n\t" \
                      "vpermilps $0xb1, %%ymm9, %%ymm12 \n\t" \
                      "vpermilps $0xb1, %%ymm10, %%ymm13 \n\t" \
                      "vpermilps $0xb1, %%ymm11, %%ymm14" \
                      : \
                      : \
                      "m" ((ul).c13), \
                      "m" ((ul).c21), \
                      "m" ((ul).c32), \
                      "m" ((uh).c13), \
                      "m" ((uh).c21), \
                      "m" ((uh).c32) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vfmaddsub231ps %%ymm8, %%ymm12, %%ymm3 \n\t" \
                      "vfmaddsub231ps %%ymm6, %%ymm13, %%ymm4 \n\t" \
                      "vfmaddsub231ps %%ymm7, %%ymm14, %%ymm5 \n\t" \
                      "vfmaddsub231ps %%ymm2, %%ymm9, %%ymm3 \n\t" \
                      "vfmaddsub231ps %%ymm0, %%ymm10, %%ymm4 \n\t" \
                      "vfmaddsub231ps %%ymm1, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_su3_pair_multiply(ul,uh) \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm3 \n\t" \
                      "vbroadcastss %1, %%xmm6 \n\t" \
                      "vbroadcastss %2, %%xmm4 \n\t" \
                      "vbroadcastss %3, %%xmm9 \n\t" \
                      "vbroadcastss %4, %%xmm10 \n\t" \
                      "vbroadcastss %5, %%xmm11 \n\t" \
                      "vinsertf128 $0x1, %%xmm9, %%ymm3, %%ymm3 \n\t" \
                      "vinsertf128 $0x1, %%xmm10, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %%xmm11, %%ymm4, %%ymm4" \
                      : \
                      : \
                      "m" ((ul).c11.re), \
                      "m" ((ul).c12.re), \
                      "m" ((ul).c21.re), \
                      "m" ((uh).c11.re), \
                      "m" ((uh).c12.re), \
                      "m" ((uh).c21.re) \
                      : \
                      "xmm3", "xmm4", "xmm6", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm7 \n\t" \
                      "vbroadcastss %1, %%xmm5 \n\t" \
                      "vbroadcastss %2, %%xmm8 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm7, %%ymm7 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm5, %%ymm5 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((ul).c22.re), \
                      "m" ((ul).c31.re), \
                      "m" ((ul).c32.re), \
                      "m" ((uh).c22.re), \
                      "m" ((uh).c31.re), \
                      "m" ((uh).c32.re) \
                      : \
                      "xmm5", "xmm7", "xmm8", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm0, %%ymm3, %%ymm3 \n\t" \
                      "vmulps %%ymm1, %%ymm6, %%ymm6 \n\t" \
                      "vmulps %%ymm0, %%ymm4, %%ymm4 \n\t" \
                      "vmulps %%ymm1, %%ymm7, %%ymm7 \n\t" \
                      "vmulps %%ymm0, %%ymm5, %%ymm5 \n\t" \
                      "vmulps %%ymm1, %%ymm8, %%ymm8 \n\t" \
                      "vaddps %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddps %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddps %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm9 \n\t"  \
                      "vbroadcastss %1, %%xmm10 \n\t" \
                      "vbroadcastss %2, %%xmm11 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm9, %%ymm9 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm10, %%ymm10 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm11, %%ymm11 \n\t" \
                      "vpermilps $0xb1, %%ymm0, %%ymm0" \
                      : \
                      : \
                      "m" ((ul).c13.re), \
                      "m" ((ul).c21.im), \
                      "m" ((ul).c33.re), \
                      "m" ((uh).c13.re), \
                      "m" ((uh).c21.im), \
                      "m" ((uh).c33.re) \
                      : \
                      "xmm0", "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm6 \n\t"  \
                      "vbroadcastss %1, %%xmm7 \n\t" \
                      "vbroadcastss %2, %%xmm8 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm7, %%ymm7 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((ul).c11.im), \
                      "m" ((ul).c23.re), \
                      "m" ((ul).c31.im), \
                      "m" ((uh).c11.im), \
                      "m" ((uh).c23.re), \
                      "m" ((uh).c31.im) \
                      : \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm2, %%ymm9, %%ymm9 \n\t" \
                      "vmulps %%ymm0, %%ymm10, %%ymm10 \n\t" \
                      "vmulps %%ymm2, %%ymm11, %%ymm11 \n\t" \
                      "vmulps %%ymm0, %%ymm6, %%ymm6 \n\t" \
                      "vmulps %%ymm2, %%ymm7, %%ymm7 \n\t" \
                      "vmulps %%ymm0, %%ymm8, %%ymm8 \n\t" \
                      "vaddps %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubps %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddps %%ymm11, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubps %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddps %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubps %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8",  \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vpermilps $0xb1, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xb1, %%ymm2, %%ymm2 \n\t" \
                      "vbroadcastss %0, %%xmm12 \n\t" \
                      "vbroadcastss %1, %%xmm13 \n\t" \
                      "vbroadcastss %2, %%xmm14 \n\t" \
                      "vbroadcastss %3, %%xmm9 \n\t" \
                      "vbroadcastss %4, %%xmm10 \n\t" \
                      "vbroadcastss %5, %%xmm11 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm9, %%ymm9 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm10, %%ymm10 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm11, %%ymm11" \
                      : \
                      : \
                      "m" ((uh).c12.im), \
                      "m" ((uh).c23.im), \
                      "m" ((uh).c32.im), \
                      "m" ((ul).c12.im), \
                      "m" ((ul).c23.im), \
                      "m" ((ul).c32.im) \
                      : \
                      "xmm1", "xmm2", "xmm9", "xmm10", \
                      "xmm11", "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm6 \n\t" \
                      "vbroadcastss %1, %%xmm7 \n\t" \
                      "vbroadcastss %2, %%xmm8 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm7, %%ymm7 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((ul).c13.im), \
                      "m" ((ul).c22.im), \
                      "m" ((ul).c33.im), \
                      "m" ((uh).c13.im), \
                      "m" ((uh).c22.im), \
                      "m" ((uh).c33.im) \
                      : \
                      "xmm6", "xmm7", "xmm8", "xmm12", \
                      "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm1, %%ymm9, %%ymm9 \n\t"  \
                      "vmulps %%ymm2, %%ymm10, %%ymm10 \n\t" \
                      "vmulps %%ymm1, %%ymm11, %%ymm11 \n\t" \
                      "vmulps %%ymm2, %%ymm6, %%ymm6 \n\t" \
                      "vmulps %%ymm1, %%ymm7, %%ymm7 \n\t" \
                      "vmulps %%ymm2, %%ymm8, %%ymm8 \n\t" \
                      "vaddsubps %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubps %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubps %%ymm11, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubps %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubps %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubps %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11")

#endif

/*
* Multiplies pairs of su3 vectors, stored in the low and high lanes of
* ymm0,..,ymm2, by the su3 matrices ul^dagger and uh^dagger, respectively.
* The vectors are assumed to be in vertical order and the products are returned
* in the same order in the registers ymm3,..,ymm5. All registers except for
* ymm15 are changed on exit.
*/

#if (defined FMA3)

#define _avx_su3_pair_inverse_multiply(ul,uh) \
__asm__ __volatile__ ("vpermilps $0xf5, %%ymm0, %%ymm6 \n\t" \
                      "vpermilps $0xf5, %%ymm1, %%ymm7 \n\t" \
                      "vpermilps $0xf5, %%ymm2, %%ymm8 \n\t" \
                      "vpermilps $0xa0, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0xa0, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xa0, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm12 \n\t" \
                      "vbroadcastsd %4, %%ymm13 \n\t" \
                      "vbroadcastsd %5, %%ymm14 \n\t" \
                      "vblendps $0xf0, %%ymm12, %%ymm9, %%ymm9 \n\t" \
                      "vblendps $0xf0, %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vblendps $0xf0, %%ymm14, %%ymm11, %%ymm11 \n\t" \
                      "vpermilps $0xb1, %%ymm9, %%ymm12 \n\t" \
                      "vpermilps $0xb1, %%ymm10, %%ymm13 \n\t" \
                      "vpermilps $0xb1, %%ymm11, %%ymm14" \
                      : \
                      : \
                      "m" ((ul).c11), \
                      "m" ((ul).c22), \
                      "m" ((ul).c33), \
                      "m" ((uh).c11), \
                      "m" ((uh).c22), \
                      "m" ((uh).c33) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm0, %%ymm9, %%ymm3 \n\t" \
                      "vmulps %%ymm1, %%ymm10, %%ymm4 \n\t" \
                      "vmulps %%ymm2, %%ymm11, %%ymm5 \n\t" \
                      "vfmsubadd231ps %%ymm6, %%ymm12, %%ymm3 \n\t" \
                      "vfmsubadd231ps %%ymm7, %%ymm13, %%ymm4 \n\t" \
                      "vfmsubadd231ps %%ymm8, %%ymm14, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm12 \n\t" \
                      "vbroadcastsd %4, %%ymm13 \n\t" \
                      "vbroadcastsd %5, %%ymm14 \n\t" \
                      "vblendps $0xf0, %%ymm12, %%ymm9, %%ymm9 \n\t" \
                      "vblendps $0xf0, %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vblendps $0xf0, %%ymm14, %%ymm11, %%ymm11 \n\t" \
                      "vpermilps $0xb1, %%ymm9, %%ymm12 \n\t" \
                      "vpermilps $0xb1, %%ymm10, %%ymm13 \n\t" \
                      "vpermilps $0xb1, %%ymm11, %%ymm14" \
                      : \
                      : \
                      "m" ((ul).c21), \
                      "m" ((ul).c32), \
                      "m" ((ul).c13), \
                      "m" ((uh).c21), \
                      "m" ((uh).c32), \
                      "m" ((uh).c13) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vfmsubadd231ps %%ymm1, %%ymm9, %%ymm3 \n\t" \
                      "vfmsubadd231ps %%ymm2, %%ymm10, %%ymm4 \n\t" \
                      "vfmsubadd231ps %%ymm0, %%ymm11, %%ymm5 \n\t" \
                      "vfmsubadd231ps %%ymm7, %%ymm12, %%ymm3 \n\t" \
                      "vfmsubadd231ps %%ymm8, %%ymm13, %%ymm4 \n\t" \
                      "vfmsubadd231ps %%ymm6, %%ymm14, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm12 \n\t" \
                      "vbroadcastsd %4, %%ymm13 \n\t" \
                      "vbroadcastsd %5, %%ymm14 \n\t" \
                      "vblendps $0xf0, %%ymm12, %%ymm9, %%ymm9 \n\t" \
                      "vblendps $0xf0, %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vblendps $0xf0, %%ymm14, %%ymm11, %%ymm11 \n\t" \
                      "vpermilps $0xb1, %%ymm9, %%ymm12 \n\t" \
                      "vpermilps $0xb1, %%ymm10, %%ymm13 \n\t" \
                      "vpermilps $0xb1, %%ymm11, %%ymm14" \
                      : \
                      : \
                      "m" ((ul).c31), \
                      "m" ((ul).c12), \
                      "m" ((ul).c23), \
                      "m" ((uh).c31), \
                      "m" ((uh).c12), \
                      "m" ((uh).c23) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vfmsubadd231ps %%ymm2, %%ymm9, %%ymm3 \n\t" \
                      "vfmsubadd231ps %%ymm0, %%ymm10, %%ymm4 \n\t" \
                      "vfmsubadd231ps %%ymm1, %%ymm11, %%ymm5 \n\t" \
                      "vfmsubadd231ps %%ymm8, %%ymm12, %%ymm3 \n\t" \
                      "vfmsubadd231ps %%ymm6, %%ymm13, %%ymm4 \n\t" \
                      "vfmsubadd231ps %%ymm7, %%ymm14, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_su3_pair_inverse_multiply(ul,uh)          \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm6 \n\t" \
                      "vbroadcastss %1, %%xmm9 \n\t" \
                      "vbroadcastss %2, %%xmm7 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm9, %%ymm9 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm7, %%ymm7" \
                      : \
                      : \
                      "m" ((ul).c11.im), \
                      "m" ((ul).c21.im), \
                      "m" ((ul).c12.im), \
                      "m" ((uh).c11.im), \
                      "m" ((uh).c21.im), \
                      "m" ((uh).c12.im) \
                      : \
                      "xmm6", "xmm7", "xmm9", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm10 \n\t" \
                      "vbroadcastss %1, %%xmm8 \n\t" \
                      "vbroadcastss %2, %%xmm11 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm10, %%ymm10 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm8, %%ymm8 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm11, %%ymm11" \
                      : \
                      : \
                      "m" ((ul).c22.im), \
                      "m" ((ul).c13.im), \
                      "m" ((ul).c23.im), \
                      "m" ((uh).c22.im), \
                      "m" ((uh).c13.im), \
                      "m" ((uh).c23.im) \
                      : \
                      "xmm8", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm0, %%ymm6, %%ymm6 \n\t" \
                      "vmulps %%ymm1, %%ymm9, %%ymm9 \n\t" \
                      "vmulps %%ymm0, %%ymm7, %%ymm7 \n\t" \
                      "vmulps %%ymm1, %%ymm10, %%ymm10 \n\t" \
                      "vmulps %%ymm0, %%ymm8, %%ymm8 \n\t" \
                      "vmulps %%ymm1, %%ymm11, %%ymm11 \n\t" \
                      "vaddps %%ymm6, %%ymm9, %%ymm9 \n\t" \
                      "vaddps %%ymm7, %%ymm10, %%ymm10 \n\t" \
                      "vaddps %%ymm8, %%ymm11, %%ymm11" \
                      : \
                      : \
                      : \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm3 \n\t"  \
                      "vbroadcastss %1, %%xmm4 \n\t" \
                      "vbroadcastss %2, %%xmm5 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm3, %%ymm3 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm4, %%ymm4 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm5, %%ymm5 \n\t" \
                      "vpermilps $0xb1, %%ymm0, %%ymm0" \
                      : \
                      : \
                      "m" ((ul).c11.re), \
                      "m" ((ul).c12.re), \
                      "m" ((ul).c13.re), \
                      "m" ((uh).c11.re), \
                      "m" ((uh).c12.re), \
                      "m" ((uh).c13.re) \
                      : \
                      "xmm0", "xmm3", "xmm4", "xmm5", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm6 \n\t" \
                      "vbroadcastss %1, %%xmm7 \n\t" \
                      "vbroadcastss %2, %%xmm8 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm7, %%ymm7 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((ul).c31.im), \
                      "m" ((ul).c32.im), \
                      "m" ((ul).c33.im), \
                      "m" ((uh).c31.im), \
                      "m" ((uh).c32.im), \
                      "m" ((uh).c33.im) \
                      : \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm0, %%ymm3, %%ymm3 \n\t" \
                      "vmulps %%ymm0, %%ymm4, %%ymm4 \n\t" \
                      "vmulps %%ymm0, %%ymm5, %%ymm5 \n\t" \
                      "vmulps %%ymm2, %%ymm6, %%ymm6 \n\t" \
                      "vmulps %%ymm2, %%ymm7, %%ymm7 \n\t" \
                      "vmulps %%ymm2, %%ymm8, %%ymm8 \n\t" \
                      "vaddsubps %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubps %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubps %%ymm11, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubps %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubps %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubps %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vpermilps $0xb1, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xb1, %%ymm2, %%ymm2 \n\t" \
                      "vbroadcastss %0, %%xmm12 \n\t" \
                      "vbroadcastss %1, %%xmm13 \n\t" \
                      "vbroadcastss %2, %%xmm14 \n\t" \
                      "vbroadcastss %3, %%xmm9 \n\t" \
                      "vbroadcastss %4, %%xmm10 \n\t" \
                      "vbroadcastss %5, %%xmm11 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm9, %%ymm9 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm10, %%ymm10 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm11, %%ymm11" \
                      : \
                      : \
                      "m" ((uh).c21.re), \
                      "m" ((uh).c32.re), \
                      "m" ((uh).c23.re), \
                      "m" ((ul).c21.re), \
                      "m" ((ul).c32.re), \
                      "m" ((ul).c23.re) \
                      : \
                      "xmm1", "xmm2", "xmm9", "xmm10", \
                      "xmm11", "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm6 \n\t" \
                      "vbroadcastss %1, %%xmm7 \n\t" \
                      "vbroadcastss %2, %%xmm8 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm7, %%ymm7 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((ul).c31.re), \
                      "m" ((ul).c22.re), \
                      "m" ((ul).c33.re), \
                      "m" ((uh).c31.re), \
                      "m" ((uh).c22.re), \
                      "m" ((uh).c33.re) \
                      : \
                      "xmm6", "xmm7", "xmm8", "xmm12", \
                      "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm1, %%ymm9, %%ymm9 \n\t" \
                      "vmulps %%ymm2, %%ymm10, %%ymm10 \n\t" \
                      "vmulps %%ymm1, %%ymm11, %%ymm11 \n\t" \
                      "vmulps %%ymm2, %%ymm6, %%ymm6 \n\t" \
                      "vmulps %%ymm1, %%ymm7, %%ymm7 \n\t" \
                      "vmulps %%ymm2, %%ymm8, %%ymm8 \n\t" \
                      "vaddps %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddps %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddps %%ymm11, %%ymm5, %%ymm5 \n\t" \
                      "vaddps %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddps %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddps %%ymm8, %%ymm5, %%ymm5 \n\t" \
                      "vpermilps $0xb1, %%ymm3, %%ymm3 \n\t" \
                      "vpermilps $0xb1, %%ymm4, %%ymm4 \n\t" \
                      "vpermilps $0xb1, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11")

#endif

/*
* Multiplies pairs of su3 vectors, stored in the low and high lanes of
* ymm0,..,ymm2, by the su3 matrices ul and uh^dagger, respectively. The
* vectors are assumed to be in vertical order and the products are returned
* in the same order in the registers ymm3,..,ymm5. All registers except
* for ymm15 are changed on exit.
*/

#if (defined FMA3)

#define _avx_su3_pair_mixed_multiply(ul,uh) \
__asm__ __volatile__ ("vpermilps $0xf5, %%ymm0, %%ymm6 \n\t" \
                      "vpermilps $0xf5, %%ymm1, %%ymm7 \n\t" \
                      "vpermilps $0xf5, %%ymm2, %%ymm8 \n\t" \
                      "vpermilps $0xa0, %%ymm0, %%ymm0 \n\t" \
                      "vpermilps $0xa0, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xa0, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm12 \n\t" \
                      "vbroadcastsd %4, %%ymm13 \n\t" \
                      "vbroadcastsd %5, %%ymm14 \n\t" \
                      "vmulps %6, %%xmm9, %%xmm9 \n\t" \
                      "vmulps %6, %%xmm10, %%xmm10 \n\t" \
                      "vmulps %6, %%xmm11, %%xmm11 \n\t" \
                      "vblendps $0xf0, %%ymm12, %%ymm9, %%ymm9 \n\t" \
                      "vblendps $0xf0, %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vblendps $0xf0, %%ymm14, %%ymm11, %%ymm11 \n\t" \
                      "vpermilps $0xb1, %%ymm9, %%ymm12 \n\t" \
                      "vpermilps $0xb1, %%ymm10, %%ymm13 \n\t" \
                      "vpermilps $0xb1, %%ymm11, %%ymm14" \
                      : \
                      : \
                      "m" ((ul).c11), \
                      "m" ((ul).c22), \
                      "m" ((ul).c33), \
                      "m" ((uh).c11), \
                      "m" ((uh).c22), \
                      "m" ((uh).c33), \
                      "m" (_sse_sgn24) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm0, %%ymm9, %%ymm3 \n\t" \
                      "vmulps %%ymm1, %%ymm10, %%ymm4 \n\t" \
                      "vmulps %%ymm2, %%ymm11, %%ymm5 \n\t" \
                      "vfmsubadd231ps %%ymm6, %%ymm12, %%ymm3 \n\t" \
                      "vfmsubadd231ps %%ymm7, %%ymm13, %%ymm4 \n\t" \
                      "vfmsubadd231ps %%ymm8, %%ymm14, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm12 \n\t" \
                      "vbroadcastsd %4, %%ymm13 \n\t" \
                      "vbroadcastsd %5, %%ymm14 \n\t" \
                      "vmulps %6, %%xmm9, %%xmm9 \n\t" \
                      "vmulps %6, %%xmm10, %%xmm10 \n\t" \
                      "vmulps %6, %%xmm11, %%xmm11 \n\t" \
                      "vblendps $0xf0, %%ymm12, %%ymm9, %%ymm9 \n\t" \
                      "vblendps $0xf0, %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vblendps $0xf0, %%ymm14, %%ymm11, %%ymm11 \n\t" \
                      "vpermilps $0xb1, %%ymm9, %%ymm12 \n\t" \
                      "vpermilps $0xb1, %%ymm10, %%ymm13 \n\t" \
                      "vpermilps $0xb1, %%ymm11, %%ymm14" \
                      : \
                      : \
                      "m" ((ul).c12), \
                      "m" ((ul).c23), \
                      "m" ((ul).c31), \
                      "m" ((uh).c21), \
                      "m" ((uh).c32), \
                      "m" ((uh).c13), \
                      "m" (_sse_sgn24) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vfmsubadd231ps %%ymm1, %%ymm9, %%ymm3 \n\t" \
                      "vfmsubadd231ps %%ymm2, %%ymm10, %%ymm4 \n\t" \
                      "vfmsubadd231ps %%ymm0, %%ymm11, %%ymm5 \n\t" \
                      "vfmsubadd231ps %%ymm7, %%ymm12, %%ymm3 \n\t" \
                      "vfmsubadd231ps %%ymm8, %%ymm13, %%ymm4 \n\t" \
                      "vfmsubadd231ps %%ymm6, %%ymm14, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm12 \n\t" \
                      "vbroadcastsd %4, %%ymm13 \n\t" \
                      "vbroadcastsd %5, %%ymm14 \n\t" \
                      "vmulps %6, %%xmm9, %%xmm9 \n\t" \
                      "vmulps %6, %%xmm10, %%xmm10 \n\t" \
                      "vmulps %6, %%xmm11, %%xmm11 \n\t" \
                      "vblendps $0xf0, %%ymm12, %%ymm9, %%ymm9 \n\t" \
                      "vblendps $0xf0, %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vblendps $0xf0, %%ymm14, %%ymm11, %%ymm11 \n\t" \
                      "vpermilps $0xb1, %%ymm9, %%ymm12 \n\t" \
                      "vpermilps $0xb1, %%ymm10, %%ymm13 \n\t" \
                      "vpermilps $0xb1, %%ymm11, %%ymm14" \
                      : \
                      : \
                      "m" ((ul).c13), \
                      "m" ((ul).c21), \
                      "m" ((ul).c32), \
                      "m" ((uh).c31), \
                      "m" ((uh).c12), \
                      "m" ((uh).c23), \
                      "m" (_sse_sgn24) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vfmsubadd231ps %%ymm2, %%ymm9, %%ymm3 \n\t" \
                      "vfmsubadd231ps %%ymm0, %%ymm10, %%ymm4 \n\t" \
                      "vfmsubadd231ps %%ymm1, %%ymm11, %%ymm5 \n\t" \
                      "vfmsubadd231ps %%ymm8, %%ymm12, %%ymm3 \n\t" \
                      "vfmsubadd231ps %%ymm6, %%ymm13, %%ymm4 \n\t" \
                      "vfmsubadd231ps %%ymm7, %%ymm14, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_su3_pair_mixed_multiply(ul,uh) \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm3 \n\t" \
                      "vbroadcastss %1, %%xmm6 \n\t" \
                      "vbroadcastss %2, %%xmm4 \n\t" \
                      "vbroadcastss %3, %%xmm9 \n\t" \
                      "vbroadcastss %4, %%xmm10 \n\t" \
                      "vbroadcastss %5, %%xmm11 \n\t" \
                      "vinsertf128 $0x1, %%xmm9, %%ymm3, %%ymm3 \n\t" \
                      "vinsertf128 $0x1, %%xmm10, %%ymm6, %%ymm6 \n\t" \
                      "vinsertf128 $0x1, %%xmm11, %%ymm4, %%ymm4" \
                      : \
                      : \
                      "m" ((ul).c11.re), \
                      "m" ((ul).c12.re), \
                      "m" ((ul).c21.re), \
                      "m" ((uh).c11.re), \
                      "m" ((uh).c21.re), \
                      "m" ((uh).c12.re) \
                      : \
                      "xmm3", "xmm4", "xmm6", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm7 \n\t" \
                      "vbroadcastss %1, %%xmm5 \n\t" \
                      "vbroadcastss %2, %%xmm8 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm7, %%ymm7 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm5, %%ymm5 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((ul).c22.re), \
                      "m" ((ul).c31.re), \
                      "m" ((ul).c32.re), \
                      "m" ((uh).c22.re), \
                      "m" ((uh).c13.re), \
                      "m" ((uh).c23.re) \
                      : \
                      "xmm5", "xmm7", "xmm8", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm0, %%ymm3, %%ymm3 \n\t" \
                      "vmulps %%ymm1, %%ymm6, %%ymm6 \n\t" \
                      "vmulps %%ymm0, %%ymm4, %%ymm4 \n\t" \
                      "vmulps %%ymm1, %%ymm7, %%ymm7 \n\t" \
                      "vmulps %%ymm0, %%ymm5, %%ymm5 \n\t" \
                      "vmulps %%ymm1, %%ymm8, %%ymm8 \n\t" \
                      "vaddps %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddps %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddps %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm9 \n\t"  \
                      "vbroadcastss %1, %%xmm10 \n\t" \
                      "vbroadcastss %2, %%xmm11 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vinsertf128 $0x1, %%xmm12, %%ymm9, %%ymm9 \n\t" \
                      "vperm2f128 $0x1, %%ymm13, %%ymm13, %%ymm13 \n\t" \
                      "vinsertf128 $0x1, %%xmm14, %%ymm11, %%ymm11 \n\t" \
                      "vsubps %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vpermilps $0xb1, %%ymm0, %%ymm0" \
                      : \
                      : \
                      "m" ((ul).c13.re), \
                      "m" ((ul).c21.im), \
                      "m" ((ul).c33.re), \
                      "m" ((uh).c31.re), \
                      "m" ((uh).c12.im), \
                      "m" ((uh).c33.re) \
                      : \
                      "xmm0", "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm6 \n\t"  \
                      "vbroadcastss %1, %%xmm7 \n\t" \
                      "vbroadcastss %2, %%xmm8 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vperm2f128 $0x1, %%ymm12, %%ymm12, %%ymm12 \n\t" \
                      "vinsertf128 $0x1, %%xmm13, %%ymm7, %%ymm7 \n\t" \
                      "vperm2f128 $0x1, %%ymm14, %%ymm14, %%ymm14 \n\t" \
                      "vsubps %%ymm12, %%ymm6, %%ymm6 \n\t" \
                      "vsubps %%ymm14, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((ul).c11.im), \
                      "m" ((ul).c23.re), \
                      "m" ((ul).c31.im), \
                      "m" ((uh).c11.im), \
                      "m" ((uh).c32.re), \
                      "m" ((uh).c13.im) \
                      : \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm2, %%ymm9, %%ymm9 \n\t" \
                      "vmulps %%ymm0, %%ymm10, %%ymm10 \n\t" \
                      "vmulps %%ymm2, %%ymm11, %%ymm11 \n\t" \
                      "vmulps %%ymm0, %%ymm6, %%ymm6 \n\t" \
                      "vmulps %%ymm2, %%ymm7, %%ymm7 \n\t" \
                      "vmulps %%ymm0, %%ymm8, %%ymm8 \n\t" \
                      "vaddps %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubps %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddps %%ymm11, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubps %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddps %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubps %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8",  \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vpermilps $0xb1, %%ymm1, %%ymm1 \n\t" \
                      "vpermilps $0xb1, %%ymm2, %%ymm2 \n\t" \
                      "vbroadcastss %0, %%xmm12 \n\t" \
                      "vbroadcastss %1, %%xmm13 \n\t" \
                      "vbroadcastss %2, %%xmm14 \n\t" \
                      "vbroadcastss %3, %%xmm9 \n\t" \
                      "vbroadcastss %4, %%xmm10 \n\t" \
                      "vbroadcastss %5, %%xmm11 \n\t" \
                      "vperm2f128 $0x1, %%ymm12, %%ymm12, %%ymm12 \n\t" \
                      "vperm2f128 $0x1, %%ymm13, %%ymm13, %%ymm13 \n\t" \
                      "vperm2f128 $0x1, %%ymm14, %%ymm14, %%ymm14 \n\t" \
                      "vsubps %%ymm12, %%ymm9, %%ymm9 \n\t" \
                      "vsubps %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vsubps %%ymm14, %%ymm11, %%ymm11" \
                      : \
                      : \
                      "m" ((uh).c21.im), \
                      "m" ((uh).c32.im), \
                      "m" ((uh).c23.im), \
                      "m" ((ul).c12.im), \
                      "m" ((ul).c23.im), \
                      "m" ((ul).c32.im) \
                      : \
                      "xmm1", "xmm2", "xmm9", "xmm10", \
                      "xmm11", "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vbroadcastss %0, %%xmm6 \n\t" \
                      "vbroadcastss %1, %%xmm7 \n\t" \
                      "vbroadcastss %2, %%xmm8 \n\t" \
                      "vbroadcastss %3, %%xmm12 \n\t" \
                      "vbroadcastss %4, %%xmm13 \n\t" \
                      "vbroadcastss %5, %%xmm14 \n\t" \
                      "vperm2f128 $0x1, %%ymm12, %%ymm12, %%ymm12 \n\t" \
                      "vperm2f128 $0x1, %%ymm13, %%ymm13, %%ymm13 \n\t" \
                      "vperm2f128 $0x1, %%ymm14, %%ymm14, %%ymm14 \n\t" \
                      "vsubps %%ymm12, %%ymm6, %%ymm6 \n\t" \
                      "vsubps %%ymm13, %%ymm7, %%ymm7 \n\t" \
                      "vsubps %%ymm14, %%ymm8, %%ymm8" \
                      : \
                      : \
                      "m" ((ul).c13.im), \
                      "m" ((ul).c22.im), \
                      "m" ((ul).c33.im), \
                      "m" ((uh).c31.im), \
                      "m" ((uh).c22.im), \
                      "m" ((uh).c33.im) \
                      : \
                      "xmm6", "xmm7", "xmm8", "xmm12", \
                      "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulps %%ymm1, %%ymm9, %%ymm9 \n\t"  \
                      "vmulps %%ymm2, %%ymm10, %%ymm10 \n\t" \
                      "vmulps %%ymm1, %%ymm11, %%ymm11 \n\t" \
                      "vmulps %%ymm2, %%ymm6, %%ymm6 \n\t" \
                      "vmulps %%ymm1, %%ymm7, %%ymm7 \n\t" \
                      "vmulps %%ymm2, %%ymm8, %%ymm8 \n\t" \
                      "vaddsubps %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubps %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubps %%ymm11, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubps %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubps %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubps %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11")

#endif

/******************************************************************************
*
*  Macros for single precision Dirac spinors in linear order
*
******************************************************************************/

/*
*  Loads the spinor s to the registers ymm0,..,ymm2 in linear order.
*/

#define _avx_spinor_load(s) \
__asm__ __volatile__ ("vmovaps %0, %%ymm0" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1) \
                      : \
                      "xmm0"); \
__asm__ __volatile__ ("vmovaps %0, %%ymm1" \
                      : \
                      : \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3), \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2) \
                      : \
                      "xmm1"); \
__asm__ __volatile__ ("vmovaps %0, %%ymm2" \
                      : \
                      : \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm2")

/*
*  Loads the spinor s to the registers ymm3,..,ymm5 in linear order.
*/

#define _avx_spinor_load_up(s) \
__asm__ __volatile__ ("vmovaps %0, %%ymm3" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1) \
                      : \
                      "xmm3"); \
__asm__ __volatile__ ("vmovaps %0, %%ymm4" \
                      : \
                      : \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3), \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2) \
                      : \
                      "xmm4"); \
__asm__ __volatile__ ("vmovaps %0, %%ymm5" \
                      : \
                      : \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm5")

/*
*  Stores the registers ymm0,..,ymm2 to the spinor s in linear order.
*/

#define _avx_spinor_store(s) \
__asm__ __volatile__ ("vmovaps %%ymm0, %0 \n\t" \
                      : \
                      "=m" ((s).c1.c1), \
                      "=m" ((s).c1.c2), \
                      "=m" ((s).c1.c3), \
                      "=m" ((s).c2.c1)); \
__asm__ __volatile__ ("vmovaps %%ymm1, %0 \n\t" \
                      : \
                      "=m" ((s).c2.c2), \
                      "=m" ((s).c2.c3), \
                      "=m" ((s).c3.c1), \
                      "=m" ((s).c3.c2)); \
__asm__ __volatile__ ("vmovaps %%ymm2, %0 \n\t" \
                      : \
                      "=m" ((s).c3.c3), \
                      "=m" ((s).c4.c1), \
                      "=m" ((s).c4.c2), \
                      "=m" ((s).c4.c3))

/*
*  Stores the registers ymm3,..,ymm5 to the spinor s in linear order.
*/

#define _avx_spinor_store_up(s) \
__asm__ __volatile__ ("vmovaps %%ymm3, %0 \n\t" \
                      : \
                      "=m" ((s).c1.c1), \
                      "=m" ((s).c1.c2), \
                      "=m" ((s).c1.c3), \
                      "=m" ((s).c2.c1)); \
__asm__ __volatile__ ("vmovaps %%ymm4, %0 \n\t" \
                      : \
                      "=m" ((s).c2.c2), \
                      "=m" ((s).c2.c3), \
                      "=m" ((s).c3.c1), \
                      "=m" ((s).c3.c2)); \
__asm__ __volatile__ ("vmovaps %%ymm5, %0 \n\t" \
                      : \
                      "=m" ((s).c3.c3), \
                      "=m" ((s).c4.c1), \
                      "=m" ((s).c4.c2), \
                      "=m" ((s).c4.c3))

/*
*  Loads (z.re,z.re,..,z.re) to ymm12 and (-z.im,z.im,..,z.im) to ymm13.
*/

#define _avx_load_cmplx(z) \
__asm__ __volatile__ ("vxorps %%ymm13, %%ymm13, %%ymm13 \n\t" \
                      "vbroadcastss %0, %%ymm12 \n\t" \
                      "vaddsubps %%ymm12, %%ymm13, %%ymm13 \n\t" \
                      "vbroadcastss %1, %%ymm12" \
                      : \
                      : \
                      "m" ((z).im), \
                      "m" ((z).re) \
                      : \
                      "xmm12", "xmm13")

/*
*  Loads (z.re,z.re,..,z.re) to ymm14 and (-z.im,z.im,..,z.im) to ymm15
*/

#define _avx_load_cmplx_up(z) \
__asm__ __volatile__ ("vxorps %%ymm15, %%ymm15, %%ymm15 \n\t" \
                      "vbroadcastss %0, %%ymm14 \n\t" \
                      "vaddsubps %%ymm14, %%ymm15, %%ymm15 \n\t" \
                      "vbroadcastss %1, %%ymm14" \
                      : \
                      : \
                      "m" ((z).im), \
                      "m" ((z).re) \
                      : \
                      "xmm14", "xmm15")

/*
*  Multiplies the spinor s by the complex number z and assigns the result to
*  ymm0,..,ymm2, assuming z was loaded using _avx_load_cmplx(z). The registers
*  ymm3,..,ymm5 are used as workspace.
*/

#if (defined FMA3)

#define _avx_mulc_spinor(s) \
_avx_spinor_load(s); \
__asm__ __volatile__ ("vpermilps $0xb1, %%ymm0, %%ymm3 \n\t" \
                      "vpermilps $0xb1, %%ymm1, %%ymm4 \n\t" \
                      "vpermilps $0xb1, %%ymm2, %%ymm5 \n\t" \
                      "vmulps %%ymm12, %%ymm0, %%ymm0 \n\t" \
                      "vmulps %%ymm12, %%ymm1, %%ymm1 \n\t" \
                      "vmulps %%ymm12, %%ymm2, %%ymm2 \n\t" \
                      "vfmadd231ps %%ymm13, %%ymm3, %%ymm0 \n\t" \
                      "vfmadd231ps %%ymm13, %%ymm4, %%ymm1 \n\t" \
                      "vfmadd231ps %%ymm13, %%ymm5, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_mulc_spinor(s) \
_avx_spinor_load(s); \
__asm__ __volatile__ ("vpermilps $0xb1, %%ymm0, %%ymm3 \n\t" \
                      "vpermilps $0xb1, %%ymm1, %%ymm4 \n\t" \
                      "vpermilps $0xb1, %%ymm2, %%ymm5 \n\t" \
                      "vmulps %%ymm12, %%ymm0, %%ymm0 \n\t" \
                      "vmulps %%ymm13, %%ymm3, %%ymm3 \n\t" \
                      "vmulps %%ymm12, %%ymm1, %%ymm1 \n\t" \
                      "vmulps %%ymm13, %%ymm4, %%ymm4 \n\t" \
                      "vmulps %%ymm12, %%ymm2, %%ymm2 \n\t" \
                      "vmulps %%ymm13, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vaddps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddps %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

#endif

/*
*  Multiplies the spinor s by the complex number z and adds the result to
*  ymm0,..,ymm2, assuming z was loaded using _avx_load_cmplx_up(z). The
*  registers ymm3,..,ymm8 are used as workspace.
*/

#if (defined FMA3)

#define _avx_mulc_spinor_add(s) \
_avx_spinor_load_up(s); \
__asm__ __volatile__ ("vpermilps $0xb1, %%ymm3, %%ymm6 \n\t" \
                      "vpermilps $0xb1, %%ymm4, %%ymm7 \n\t" \
                      "vpermilps $0xb1, %%ymm5, %%ymm8 \n\t" \
                      "vfmadd231ps %%ymm14, %%ymm3, %%ymm0 \n\t" \
                      "vfmadd231ps %%ymm14, %%ymm4, %%ymm1 \n\t" \
                      "vfmadd231ps %%ymm14, %%ymm5, %%ymm2 \n\t" \
                      "vfmadd231ps %%ymm15, %%ymm6, %%ymm0 \n\t" \
                      "vfmadd231ps %%ymm15, %%ymm7, %%ymm1 \n\t" \
                      "vfmadd231ps %%ymm15, %%ymm8, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm6", "xmm7", "xmm8")

#else

#define _avx_mulc_spinor_add(s) \
_avx_spinor_load_up(s); \
__asm__ __volatile__ ("vpermilps $0xb1, %%ymm3, %%ymm6 \n\t" \
                      "vpermilps $0xb1, %%ymm4, %%ymm7 \n\t" \
                      "vpermilps $0xb1, %%ymm5, %%ymm8 \n\t" \
                      "vmulps %%ymm14, %%ymm3, %%ymm3 \n\t" \
                      "vmulps %%ymm15, %%ymm6, %%ymm6 \n\t" \
                      "vmulps %%ymm14, %%ymm4, %%ymm4 \n\t" \
                      "vmulps %%ymm15, %%ymm7, %%ymm7 \n\t" \
                      "vmulps %%ymm14, %%ymm5, %%ymm5 \n\t" \
                      "vmulps %%ymm15, %%ymm8, %%ymm8 \n\t" \
                      "vaddps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddps %%ymm5, %%ymm2, %%ymm2 \n\t" \
                      "vaddps %%ymm6, %%ymm0, %%ymm0 \n\t" \
                      "vaddps %%ymm7, %%ymm1, %%ymm1 \n\t" \
                      "vaddps %%ymm8, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8")

#endif

/*
*  Broadcasts the real number c to ymm12 and ymm13.
*/

#define _avx_load_real(c) \
__asm__ __volatile__ ("vbroadcastss %0, %%ymm12 \n\t" \
                      "vbroadcastss %0, %%ymm13" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm12", "xmm13")

/*
*  Broadcasts the real number c to ymm14 and ymm15.
*/

#define _avx_load_real_up(c) \
__asm__ __volatile__ ("vbroadcastss %0, %%ymm14 \n\t" \
                      "vbroadcastss %0, %%ymm15" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm14", "xmm15")

/*
*  Multiplies the spinor s by the real number c and assigns the result to
*  ymm0,..,ymm2, assuming c was loaded using _avx_load_real(c).
*/

#define _avx_mulr_spinor(s) \
_avx_spinor_load(s); \
__asm__ __volatile__ ("vmulps %%ymm12, %%ymm0, %%ymm0 \n\t" \
                      "vmulps %%ymm13, %%ymm1, %%ymm1 \n\t" \
                      "vmulps %%ymm12, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
*  Multiplies the spinor s by the real number c and adds the result to
*  ymm0,..,ymm2, assuming c was loaded using _avx_load_real_up(c). The
*  registers ymm3,..,ymm5 are used as workspace.
*/

#if (defined FMA3)

#define _avx_mulr_spinor_add(s) \
_avx_spinor_load_up(s); \
__asm__ __volatile__ ("vfmadd231ps %%ymm14, %%ymm3, %%ymm0 \n\t" \
                      "vfmadd231ps %%ymm15, %%ymm4, %%ymm1 \n\t" \
                      "vfmadd231ps %%ymm14, %%ymm5, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

#else

#define _avx_mulr_spinor_add(s) \
_avx_spinor_load_up(s); \
__asm__ __volatile__ ("vmulps %%ymm14, %%ymm3, %%ymm3 \n\t" \
                      "vmulps %%ymm15, %%ymm4, %%ymm4 \n\t" \
                      "vmulps %%ymm14, %%ymm5, %%ymm5 \n\t" \
                      "vaddps %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddps %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddps %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#endif

/*******************************************************************************
*
* Macros operating on double precision data
*
*******************************************************************************/

/*******************************************************************************
*
* Macros for su3_vector data
*
* Most of these macros operate on pairs of su3 vectors that are stored
* in the low and high lanes of ymm0,..,ymm2 or ymm3,..,ymm5. For example,
*
* ymm0 <- sl.c1.re,sl.c1.im,sh.c1.re,sh.c1.im
* ymm1 <- sl.c2.re,sl.c2.im,sh.c2.re,sh.c2.im
* ymm2 <- sl.c3.re,sl.c3.im,sh.c3.re,sh.c3.im
*
* (where sl and sh are of type su3_vector).
*
*******************************************************************************/

/*
* Loads two su3 vectors sl and sh to the low and high lanes of ymm0,..,ymm2.
*/

#define _avx_pair_load_dble(sl,sh) \
__asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t" \
                      "vmovapd %1, %%xmm1 \n\t" \
                      "vmovapd %2, %%xmm2 \n\t" \
                      "vinsertf128 $0x1, %3, %%ymm0, %%ymm0 \n\t" \
                      "vinsertf128 $0x1, %4, %%ymm1, %%ymm1 \n\t" \
                      "vinsertf128 $0x1, %5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      "m" ((sl).c1), \
                      "m" ((sl).c2), \
                      "m" ((sl).c3), \
                      "m" ((sh).c1), \
                      "m" ((sh).c2), \
                      "m" ((sh).c3) \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Loads two su3 vectors sl and sh to the low and high lanes of ymm3,..,ymm5.
*/

#define _avx_pair_load_up_dble(sl,sh) \
__asm__ __volatile__ ("vmovapd %0, %%xmm3 \n\t" \
                      "vmovapd %1, %%xmm4 \n\t" \
                      "vmovapd %2, %%xmm5 \n\t" \
                      "vinsertf128 $0x1, %3, %%ymm3, %%ymm3 \n\t" \
                      "vinsertf128 $0x1, %4, %%ymm4, %%ymm4 \n\t" \
                      "vinsertf128 $0x1, %5, %%ymm5, %%ymm5" \
                      : \
                      : \
                      "m" ((sl).c1), \
                      "m" ((sl).c2), \
                      "m" ((sl).c3), \
                      "m" ((sh).c1), \
                      "m" ((sh).c2), \
                      "m" ((sh).c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Stores the low and high lanes of ymm0,..,ymm2 to the su3 vectors rl and rh.
*/

#define _avx_pair_store_dble(rl,rh) \
__asm__ __volatile__ ("vmovapd %%xmm0, %0 \n\t" \
                      "vmovapd %%xmm1, %1 \n\t" \
                      "vmovapd %%xmm2, %2 \n\t" \
                      "vextractf128 $0x1, %%ymm0, %3 \n\t" \
                      "vextractf128 $0x1, %%ymm1, %4 \n\t" \
                      "vextractf128 $0x1, %%ymm2, %5" \
                      : \
                      "=m" ((rl).c1), \
                      "=m" ((rl).c2), \
                      "=m" ((rl).c3), \
                      "=m" ((rh).c1), \
                      "=m" ((rh).c2), \
                      "=m" ((rh).c3))

/*
* Stores the low and high lanes of ymm3,..,ymm5 to the su3 vectors rl and rh.
*/

#define _avx_pair_store_up_dble(rl,rh) \
__asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t" \
                      "vmovapd %%xmm4, %1 \n\t" \
                      "vmovapd %%xmm5, %2 \n\t" \
                      "vextractf128 $0x1, %%ymm3, %3 \n\t" \
                      "vextractf128 $0x1, %%ymm4, %4 \n\t" \
                      "vextractf128 $0x1, %%ymm5, %5" \
                      : \
                      "=m" ((rl).c1), \
                      "=m" ((rl).c2), \
                      "=m" ((rl).c3), \
                      "=m" ((rh).c1), \
                      "=m" ((rh).c2), \
                      "=m" ((rh).c3))

/*
* Loads the components of a Weyl spinor s to ymm0,..,ymm2 in linear order.
*/

#define _avx_weyl_load_dble(s) \
__asm__ __volatile__ ("vmovapd %0, %%ymm0 \n\t" \
                      "vmovapd %2, %%ymm1 \n\t" \
                      "vmovapd %4, %%ymm2" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Loads the components of a Weyl spinor s to ymm3,..,ymm5 in linear order.
*/

#define _avx_weyl_load_up_dble(s) \
__asm__ __volatile__ ("vmovapd %0, %%ymm3 \n\t" \
                      "vmovapd %2, %%ymm4 \n\t" \
                      "vmovapd %4, %%ymm5" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Stores ymm0,..,ymm2 to the components of a Weyl spinor s in linear order.
*/

#define _avx_weyl_store_dble(s) \
__asm__ __volatile__ ("vmovapd %%ymm0, %0 \n\t" \
                      "vmovapd %%ymm1, %2 \n\t" \
                      "vmovapd %%ymm2, %4" \
                      : \
                      "=m" ((s).c1.c1), \
                      "=m" ((s).c1.c2), \
                      "=m" ((s).c1.c3), \
                      "=m" ((s).c2.c1), \
                      "=m" ((s).c2.c2), \
                      "=m" ((s).c2.c3))

/*
* Stores ymm3,..,ymm5 to the components of a Weyl spinor s in linear order.
*/

#define _avx_weyl_store_up_dble(s) \
__asm__ __volatile__ ("vmovapd %%ymm3, %0 \n\t" \
                      "vmovapd %%ymm4, %2 \n\t" \
                      "vmovapd %%ymm5, %4" \
                      : \
                      "=m" ((s).c1.c1), \
                      "=m" ((s).c1.c2), \
                      "=m" ((s).c1.c3), \
                      "=m" ((s).c2.c1), \
                      "=m" ((s).c2.c2), \
                      "=m" ((s).c2.c3))

/*
* Adds ymm3,..,ymm5 to ymm0,..,ymm2.
*/

#define _avx_vector_add_dble() \
__asm__ __volatile__ ("vaddpd %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Subtracts ymm3,..,ymm5 from ymm0,..,ymm2.
*/

#define _avx_vector_sub_dble() \
__asm__ __volatile__ ("vsubpd %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vsubpd %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vsubpd %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Multiplies the high lanes of ymm3,..,ymm5 by -1 and adds these registers
* to ymm0,..,ymm2.
*/

#if (defined FMA3)

#define _avx_vector_addsub_dble() \
__asm__ __volatile__ ("vfmadd231pd %0, %%ymm3, %%ymm0 \n\t" \
                      "vfmadd231pd %0, %%ymm4, %%ymm1 \n\t" \
                      "vfmadd231pd %0, %%ymm5, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn34_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2")

#else

#define _avx_vector_addsub_dble() \
__asm__ __volatile__ ("vmulpd %0, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %0, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %0, %%ymm5, %%ymm5 \n\t" \
                      "vaddpd %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn34_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#endif

/*
* Multiplies the low lanes of ymm3,..,ymm5 by -1 and adds these registers
* to ymm0,..,ymm2.
*/

#if (defined FMA3)

#define _avx_vector_subadd_dble() \
__asm__ __volatile__ ("vfmadd231pd %0, %%ymm3, %%ymm0 \n\t" \
                      "vfmadd231pd %0, %%ymm4, %%ymm1 \n\t" \
                      "vfmadd231pd %0, %%ymm5, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn12_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2")

#else

#define _avx_vector_subadd_dble() \
__asm__ __volatile__ ("vmulpd %0, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %0, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %0, %%ymm5, %%ymm5 \n\t" \
                      "vaddpd %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn12_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#endif

/*
* Multiplies the registers ymm3,..,ymm5 by i and adds them to ymm0,..,ymm2.
*/

#define _avx_vector_i_add_dble() \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm3, %%ymm3 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm4 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubpd %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddsubpd %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddsubpd %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies the registers ymm3,..,ymm5 by i and subtracts them from
* ymm0,..,ymm2.
*/

#if (defined FMA3)

#define _avx_vector_i_sub_dble() \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm3, %%ymm3 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm4 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm5 \n\t" \
                      "vfmadd231pd %0, %%ymm3, %%ymm0 \n\t" \
                      "vfmadd231pd %0, %%ymm4, %%ymm1 \n\t" \
                      "vfmadd231pd %0, %%ymm5, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn24_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_vector_i_sub_dble() \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm3, %%ymm3 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm4 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm5 \n\t" \
                      "vmulpd %0, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %0, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %0, %%ymm5, %%ymm5 \n\t" \
                      "vaddpd %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn24_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#endif

/*
* Exchanges the high and low lanes of ymm3,..,ymm5, multiplies them by i
* and adds the result to ymm0,..,ymm2.
*/

#define _avx_vector_xch_i_add_dble() \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm3, %%ymm3 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm4 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm5 \n\t" \
                      "vperm2f128 $0x1, %%ymm3, %%ymm3, %%ymm3 \n\t" \
                      "vperm2f128 $0x1, %%ymm4, %%ymm4, %%ymm4 \n\t" \
                      "vperm2f128 $0x1, %%ymm5, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubpd %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddsubpd %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddsubpd %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Exchanges the high and low lanes of ymm3,..,ymm5, multiplies them by i
* and subtracts the result from ymm0,..,ymm2.
*/

#if (defined FMA3)

#define _avx_vector_xch_i_sub_dble() \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm3, %%ymm3 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm4 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm5 \n\t" \
                      "vperm2f128 $0x1, %%ymm3, %%ymm3, %%ymm3 \n\t" \
                      "vperm2f128 $0x1, %%ymm4, %%ymm4, %%ymm4 \n\t" \
                      "vperm2f128 $0x1, %%ymm5, %%ymm5, %%ymm5 \n\t" \
                      "vfmadd231pd %0, %%ymm3, %%ymm0 \n\t" \
                      "vfmadd231pd %0, %%ymm4, %%ymm1 \n\t" \
                      "vfmadd231pd %0, %%ymm5, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn24_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_vector_xch_i_sub_dble() \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm3, %%ymm3 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm4 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm5 \n\t" \
                      "vperm2f128 $0x1, %%ymm3, %%ymm3, %%ymm3 \n\t" \
                      "vperm2f128 $0x1, %%ymm4, %%ymm4, %%ymm4 \n\t" \
                      "vperm2f128 $0x1, %%ymm5, %%ymm5, %%ymm5 \n\t" \
                      "vmulpd %0, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %0, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %0, %%ymm5, %%ymm5 \n\t" \
                      "vaddpd %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn24_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#endif

/*
* Multiplies the low and high lanes of ymm3,..,ymm5 by i and -i
* respectively and adds these registers to ymm0,..,ymm2.
*/

#if (defined FMA3)

#define _avx_vector_i_addsub_dble() \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm3, %%ymm3 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm4 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm5 \n\t" \
                      "vfmadd231pd %0, %%ymm3, %%ymm0 \n\t" \
                      "vfmadd231pd %0, %%ymm4, %%ymm1 \n\t" \
                      "vfmadd231pd %0, %%ymm5, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn14_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_vector_i_addsub_dble() \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm3, %%ymm3 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm4 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm5 \n\t" \
                      "vmulpd %0, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %0, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %0, %%ymm5, %%ymm5 \n\t" \
                      "vaddpd %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn14_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#endif

/*
* Multiplies the low and high words of ymm3,..,ymm5 by -i and i
* respectively and adds these registers to ymm0,..,ymm2.
*/

#if (defined FMA3)

#define _avx_vector_i_subadd_dble() \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm3, %%ymm3 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm4 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm5 \n\t" \
                      "vfmadd231pd %0, %%ymm3, %%ymm0 \n\t" \
                      "vfmadd231pd %0, %%ymm4, %%ymm1 \n\t" \
                      "vfmadd231pd %0, %%ymm5, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn23_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_vector_i_subadd_dble() \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm3, %%ymm3 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm4 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm5 \n\t" \
                      "vmulpd %0, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %0, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %0, %%ymm5, %%ymm5 \n\t" \
                      "vaddpd %%ymm3, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm4, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      "m" (_avx_sgn23_dble) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#endif

/*
* Exchanges the high and low lanes of ymm3,..,ymm5.
*/

#define _avx_vector_xch_dble() \
__asm__ __volatile__ ("vperm2f128 $0x1, %%ymm3, %%ymm3, %%ymm3 \n\t" \
                      "vperm2f128 $0x1, %%ymm4, %%ymm4, %%ymm4 \n\t" \
                      "vperm2f128 $0x1, %%ymm5, %%ymm5, %%ymm5" \
                      : \
                      : \
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

#if (defined FMA3)

#define _avx_su3_multiply_dble(u) \
__asm__ __volatile__ ("vpermpd $0x50, %%ymm0, %%ymm0 \n\t" \
                      "vpermpd $0x50, %%ymm1, %%ymm1 \n\t" \
                      "vpermpd $0x50, %%ymm2, %%ymm2 \n\t" \
                      "vbroadcastf128 %0, %%ymm6 \n\t" \
                      "vbroadcastf128 %1, %%ymm7 \n\t" \
                      "vbroadcastf128 %2, %%ymm8 \n\t" \
                      "vmulpd %%ymm0, %%ymm6, %%ymm3 \n\t" \
                      "vmulpd %%ymm1, %%ymm7, %%ymm4 \n\t" \
                      "vmulpd %%ymm2, %%ymm8, %%ymm5" \
                      : \
                      : \
                      "m" ((u).c11), \
                      "m" ((u).c22), \
                      "m" ((u).c33) \
                      : \
                      "xmm0", "xmm1", "xmm2", "xmm3", \
                      "xmm4", "xmm5", "xmm6", "xmm7", \
                      "xmm8"); \
__asm__ __volatile__ ("vbroadcastf128 %0, %%ymm9 \n\t"  \
                      "vbroadcastf128 %1, %%ymm10 \n\t" \
                      "vbroadcastf128 %2, %%ymm11 \n\t" \
                      "vbroadcastf128 %3, %%ymm12 \n\t" \
                      "vbroadcastf128 %4, %%ymm13 \n\t" \
                      "vbroadcastf128 %5, %%ymm14 \n\t" \
                      "vfmadd231pd %%ymm1, %%ymm9, %%ymm3 \n\t" \
                      "vfmadd231pd %%ymm2, %%ymm10, %%ymm4 \n\t" \
                      "vfmadd231pd %%ymm0, %%ymm11, %%ymm5 \n\t" \
                      "vfmadd231pd %%ymm2, %%ymm12, %%ymm3 \n\t" \
                      "vfmadd231pd %%ymm0, %%ymm13, %%ymm4 \n\t" \
                      "vfmadd231pd %%ymm1, %%ymm14, %%ymm5" \
                      : \
                      : \
                      "m" ((u).c12), \
                      "m" ((u).c23), \
                      "m" ((u).c31), \
                      "m" ((u).c13), \
                      "m" ((u).c21), \
                      "m" ((u).c32) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm9", \
                      "xmm10", "xmm11", "xmm12", "xmm13", \
                      "xmm14"); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm3, %%xmm6 \n\t" \
                      "vextractf128 $0x1, %%ymm4, %%xmm7 \n\t" \
                      "vextractf128 $0x1, %%ymm5, %%xmm8 \n\t" \
                      "vpermilpd $0x1, %%xmm6, %%xmm6 \n\t" \
                      "vpermilpd $0x1, %%xmm7, %%xmm7 \n\t" \
                      "vpermilpd $0x1, %%xmm8, %%xmm8 \n\t" \
                      "vaddsubpd %%xmm6, %%xmm3, %%xmm3 \n\t" \
                      "vaddsubpd %%xmm7, %%xmm4, %%xmm4 \n\t" \
                      "vaddsubpd %%xmm8, %%xmm5, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm6", \
                      "xmm7", "xmm8")

#else

#define _avx_su3_multiply_dble(u) \
__asm__ __volatile__ ("vinsertf128 $0x1, %%xmm0, %%ymm0, %%ymm0 \n\t" \
                      "vinsertf128 $0x1, %%xmm1, %%ymm1, %%ymm1 \n\t" \
                      "vinsertf128 $0x1, %%xmm2, %%ymm2, %%ymm2 \n\t" \
                      "vpermilpd $0xc, %%ymm0, %%ymm0 \n\t" \
                      "vpermilpd $0xc, %%ymm1, %%ymm1 \n\t" \
                      "vpermilpd $0xc, %%ymm2, %%ymm2 \n\t" \
                      "vbroadcastf128 %0, %%ymm6 \n\t" \
                      "vbroadcastf128 %1, %%ymm7 \n\t" \
                      "vbroadcastf128 %2, %%ymm8 \n\t" \
                      "vmulpd %%ymm0, %%ymm6, %%ymm3 \n\t" \
                      "vmulpd %%ymm1, %%ymm7, %%ymm4 \n\t" \
                      "vmulpd %%ymm2, %%ymm8, %%ymm5" \
                      : \
                      : \
                      "m" ((u).c11), \
                      "m" ((u).c22), \
                      "m" ((u).c33) \
                      : \
                      "xmm0", "xmm1", "xmm2", "xmm3", \
                      "xmm4", "xmm5", "xmm6", "xmm7", \
                      "xmm8"); \
__asm__ __volatile__ ("vbroadcastf128 %0, %%ymm9 \n\t"  \
                      "vbroadcastf128 %1, %%ymm10 \n\t" \
                      "vbroadcastf128 %2, %%ymm11 \n\t" \
                      "vmulpd %%ymm1, %%ymm9, %%ymm12 \n\t" \
                      "vmulpd %%ymm2, %%ymm10, %%ymm13 \n\t" \
                      "vmulpd %%ymm0, %%ymm11, %%ymm14 \n\t" \
                      "vaddpd %%ymm12, %%ymm3, %%ymm3 \n\t" \
                      "vaddpd %%ymm13, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm14, %%ymm5, %%ymm5" \
                      : \
                      : \
                      "m" ((u).c12), \
                      "m" ((u).c23), \
                      "m" ((u).c31) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm9", \
                      "xmm10", "xmm11", "xmm12", "xmm13", \
                      "xmm14"); \
__asm__ __volatile__ ("vbroadcastf128 %0, %%ymm6 \n\t"  \
                      "vbroadcastf128 %1, %%ymm7 \n\t" \
                      "vbroadcastf128 %2, %%ymm8 \n\t" \
                      "vmulpd %%ymm2, %%ymm6, %%ymm9 \n\t" \
                      "vmulpd %%ymm0, %%ymm7, %%ymm10 \n\t" \
                      "vmulpd %%ymm1, %%ymm8, %%ymm11 \n\t" \
                      "vaddpd %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddpd %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm11, %%ymm5, %%ymm5" \
                      : \
                      : \
                      "m" ((u).c13), \
                      "m" ((u).c21), \
                      "m" ((u).c32) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm6", \
                      "xmm7", "xmm8", "xmm9", "xmm10", \
                      "xmm11"); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm3, %%xmm6 \n\t" \
                      "vextractf128 $0x1, %%ymm4, %%xmm7 \n\t" \
                      "vextractf128 $0x1, %%ymm5, %%xmm8 \n\t" \
                      "vpermilpd $0x1, %%xmm6, %%xmm6 \n\t" \
                      "vpermilpd $0x1, %%xmm7, %%xmm7 \n\t" \
                      "vpermilpd $0x1, %%xmm8, %%xmm8 \n\t" \
                      "vaddsubpd %%xmm6, %%xmm3, %%xmm3 \n\t" \
                      "vaddsubpd %%xmm7, %%xmm4, %%xmm4 \n\t" \
                      "vaddsubpd %%xmm8, %%xmm5, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm6", \
                      "xmm7", "xmm8")

#endif

/*
* Multiplies an su3 vector s with an su3 matrix u^dagger, assuming s is
* stored in  xmm0,xmm1,xmm2.
*
* On output the result is in xmm3,xmm4,xmm5 and all registers except
* for xmm15 are changed.
*/

#if (defined FMA3)

#define _avx_su3_inverse_multiply_dble(u) \
__asm__ __volatile__ ("vpermpd $0x50, %%ymm0, %%ymm0 \n\t" \
                      "vpermpd $0x50, %%ymm1, %%ymm1 \n\t" \
                      "vpermpd $0x50, %%ymm2, %%ymm2 \n\t" \
                      "vbroadcastf128 %0, %%ymm6 \n\t" \
                      "vbroadcastf128 %1, %%ymm7 \n\t" \
                      "vbroadcastf128 %2, %%ymm8 \n\t" \
                      "vmulpd %%ymm0, %%ymm6, %%ymm3 \n\t" \
                      "vmulpd %%ymm1, %%ymm7, %%ymm4 \n\t" \
                      "vmulpd %%ymm2, %%ymm8, %%ymm5" \
                      : \
                      : \
                      "m" ((u).c11), \
                      "m" ((u).c22), \
                      "m" ((u).c33) \
                      : \
                      "xmm0", "xmm1", "xmm2", "xmm3", \
                      "xmm4", "xmm5", "xmm6", "xmm7", \
                      "xmm8"); \
__asm__ __volatile__ ("vbroadcastf128 %0, %%ymm9 \n\t"  \
                      "vbroadcastf128 %1, %%ymm10 \n\t" \
                      "vbroadcastf128 %2, %%ymm11 \n\t" \
                      "vbroadcastf128 %3, %%ymm12 \n\t" \
                      "vbroadcastf128 %4, %%ymm13 \n\t" \
                      "vbroadcastf128 %5, %%ymm14 \n\t" \
                      "vfmadd231pd %%ymm1, %%ymm9, %%ymm3 \n\t" \
                      "vfmadd231pd %%ymm2, %%ymm10, %%ymm4 \n\t" \
                      "vfmadd231pd %%ymm0, %%ymm11, %%ymm5 \n\t" \
                      "vfmadd231pd %%ymm2, %%ymm12, %%ymm3 \n\t" \
                      "vfmadd231pd %%ymm0, %%ymm13, %%ymm4 \n\t" \
                      "vfmadd231pd %%ymm1, %%ymm14, %%ymm5" \
                      : \
                      : \
                      "m" ((u).c21), \
                      "m" ((u).c32), \
                      "m" ((u).c13), \
                      "m" ((u).c31), \
                      "m" ((u).c12), \
                      "m" ((u).c23) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm9", \
                      "xmm10", "xmm11", "xmm12", "xmm13", \
                      "xmm14"); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm3, %%xmm6 \n\t" \
                      "vextractf128 $0x1, %%ymm4, %%xmm7 \n\t" \
                      "vextractf128 $0x1, %%ymm5, %%xmm8 \n\t" \
                      "vpermilpd $0x1, %%xmm3, %%xmm3 \n\t" \
                      "vpermilpd $0x1, %%xmm4, %%xmm4 \n\t" \
                      "vpermilpd $0x1, %%xmm5, %%xmm5 \n\t" \
                      "vaddsubpd %%xmm3, %%xmm6, %%xmm6 \n\t" \
                      "vaddsubpd %%xmm4, %%xmm7, %%xmm7 \n\t" \
                      "vaddsubpd %%xmm5, %%xmm8, %%xmm8 \n\t" \
                      "vpermilpd $0x1, %%xmm6, %%xmm3 \n\t" \
                      "vpermilpd $0x1, %%xmm7, %%xmm4 \n\t" \
                      "vpermilpd $0x1, %%xmm8, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm6", \
                      "xmm7", "xmm8")

#else

#define _avx_su3_inverse_multiply_dble(u) \
__asm__ __volatile__ ("vinsertf128 $0x1, %%xmm0, %%ymm0, %%ymm0 \n\t" \
                      "vinsertf128 $0x1, %%xmm1, %%ymm1, %%ymm1 \n\t" \
                      "vinsertf128 $0x1, %%xmm2, %%ymm2, %%ymm2 \n\t" \
                      "vpermilpd $0xc, %%ymm0, %%ymm0 \n\t" \
                      "vpermilpd $0xc, %%ymm1, %%ymm1 \n\t" \
                      "vpermilpd $0xc, %%ymm2, %%ymm2 \n\t" \
                      "vbroadcastf128 %0, %%ymm6 \n\t" \
                      "vbroadcastf128 %1, %%ymm7 \n\t" \
                      "vbroadcastf128 %2, %%ymm8 \n\t" \
                      "vmulpd %%ymm0, %%ymm6, %%ymm3 \n\t" \
                      "vmulpd %%ymm1, %%ymm7, %%ymm4 \n\t" \
                      "vmulpd %%ymm2, %%ymm8, %%ymm5" \
                      : \
                      : \
                      "m" ((u).c11), \
                      "m" ((u).c22), \
                      "m" ((u).c33) \
                      : \
                      "xmm0", "xmm1", "xmm2", "xmm3", \
                      "xmm4", "xmm5", "xmm6", "xmm7", \
                      "xmm8"); \
__asm__ __volatile__ ("vbroadcastf128 %0, %%ymm9 \n\t"  \
                      "vbroadcastf128 %1, %%ymm10 \n\t" \
                      "vbroadcastf128 %2, %%ymm11 \n\t" \
                      "vmulpd %%ymm1, %%ymm9, %%ymm12 \n\t" \
                      "vmulpd %%ymm2, %%ymm10, %%ymm13 \n\t" \
                      "vmulpd %%ymm0, %%ymm11, %%ymm14 \n\t" \
                      "vaddpd %%ymm12, %%ymm3, %%ymm3 \n\t" \
                      "vaddpd %%ymm13, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm14, %%ymm5, %%ymm5" \
                      : \
                      : \
                      "m" ((u).c21), \
                      "m" ((u).c32), \
                      "m" ((u).c13) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm9", \
                      "xmm10", "xmm11", "xmm12", "xmm13", \
                      "xmm14"); \
__asm__ __volatile__ ("vbroadcastf128 %0, %%ymm6 \n\t"  \
                      "vbroadcastf128 %1, %%ymm7 \n\t" \
                      "vbroadcastf128 %2, %%ymm8 \n\t" \
                      "vmulpd %%ymm2, %%ymm6, %%ymm9 \n\t" \
                      "vmulpd %%ymm0, %%ymm7, %%ymm10 \n\t" \
                      "vmulpd %%ymm1, %%ymm8, %%ymm11 \n\t" \
                      "vaddpd %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddpd %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm11, %%ymm5, %%ymm5" \
                      : \
                      : \
                      "m" ((u).c31), \
                      "m" ((u).c12), \
                      "m" ((u).c23) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm6", \
                      "xmm7", "xmm8", "xmm9", "xmm10", \
                      "xmm11"); \
__asm__ __volatile__ ("vextractf128 $0x1, %%ymm3, %%xmm6 \n\t" \
                      "vextractf128 $0x1, %%ymm4, %%xmm7 \n\t" \
                      "vextractf128 $0x1, %%ymm5, %%xmm8 \n\t" \
                      "vpermilpd $0x1, %%xmm3, %%xmm3 \n\t" \
                      "vpermilpd $0x1, %%xmm4, %%xmm4 \n\t" \
                      "vpermilpd $0x1, %%xmm5, %%xmm5 \n\t" \
                      "vaddsubpd %%xmm3, %%xmm6, %%xmm6 \n\t" \
                      "vaddsubpd %%xmm4, %%xmm7, %%xmm7 \n\t" \
                      "vaddsubpd %%xmm5, %%xmm8, %%xmm8 \n\t" \
                      "vpermilpd $0x1, %%xmm6, %%xmm3 \n\t" \
                      "vpermilpd $0x1, %%xmm7, %%xmm4 \n\t" \
                      "vpermilpd $0x1, %%xmm8, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm6", \
                      "xmm7", "xmm8")

#endif

/*
* Multiplies a pair sl,sh of su3 vectors by an su3 matrix u, assuming sl and
* sh are in the low and high lanes of ymm0,..,ymm2. On output the result is
* in ymm3,..,ymm5 and all registers except for ymm15 are changed.
*/

#if (defined FMA3)

#define _avx_su3_multiply_pair_dble(u) \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm0, %%ymm6 \n\t" \
                      "vpermilpd $0x5, %%ymm1, %%ymm7 \n\t" \
                      "vpermilpd $0x5, %%ymm2, %%ymm8 \n\t" \
                      "vbroadcastsd %0, %%ymm12 \n\t" \
                      "vbroadcastsd %1, %%ymm13 \n\t" \
                      "vbroadcastsd %2, %%ymm14 \n\t" \
                      "vbroadcastsd %3, %%ymm9 \n\t" \
                      "vbroadcastsd %4, %%ymm10 \n\t" \
                      "vbroadcastsd %5, %%ymm11" \
                      : \
                      : \
                      "m" ((u).c11.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c33.im), \
                      "m" ((u).c11.re), \
                      "m" ((u).c22.re), \
                      "m" ((u).c33.re) \
                      : \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulpd %%ymm6, %%ymm12, %%ymm3 \n\t" \
                      "vmulpd %%ymm7, %%ymm13, %%ymm4 \n\t" \
                      "vmulpd %%ymm8, %%ymm14, %%ymm5 \n\t" \
                      "vfmaddsub231pd %%ymm0, %%ymm9, %%ymm3 \n\t" \
                      "vfmaddsub231pd %%ymm1, %%ymm10, %%ymm4 \n\t" \
                      "vfmaddsub231pd %%ymm2, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm12 \n\t" \
                      "vbroadcastsd %1, %%ymm13 \n\t" \
                      "vbroadcastsd %2, %%ymm14 \n\t" \
                      "vbroadcastsd %3, %%ymm9 \n\t" \
                      "vbroadcastsd %4, %%ymm10 \n\t" \
                      "vbroadcastsd %5, %%ymm11" \
                      : \
                      : \
                      "m" ((u).c12.im), \
                      "m" ((u).c23.im), \
                      "m" ((u).c31.im), \
                      "m" ((u).c12.re), \
                      "m" ((u).c23.re), \
                      "m" ((u).c31.re) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vfmaddsub231pd %%ymm7, %%ymm12, %%ymm3 \n\t" \
                      "vfmaddsub231pd %%ymm8, %%ymm13, %%ymm4 \n\t" \
                      "vfmaddsub231pd %%ymm6, %%ymm14, %%ymm5 \n\t" \
                      "vfmaddsub231pd %%ymm1, %%ymm9, %%ymm3 \n\t" \
                      "vfmaddsub231pd %%ymm2, %%ymm10, %%ymm4 \n\t" \
                      "vfmaddsub231pd %%ymm0, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm12 \n\t" \
                      "vbroadcastsd %1, %%ymm13 \n\t" \
                      "vbroadcastsd %2, %%ymm14 \n\t" \
                      "vbroadcastsd %3, %%ymm9 \n\t" \
                      "vbroadcastsd %4, %%ymm10 \n\t" \
                      "vbroadcastsd %5, %%ymm11" \
                      : \
                      : \
                      "m" ((u).c13.im), \
                      "m" ((u).c21.im), \
                      "m" ((u).c32.im), \
                      "m" ((u).c13.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c32.re) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vfmaddsub231pd %%ymm8, %%ymm12, %%ymm3 \n\t" \
                      "vfmaddsub231pd %%ymm6, %%ymm13, %%ymm4 \n\t" \
                      "vfmaddsub231pd %%ymm7, %%ymm14, %%ymm5 \n\t" \
                      "vfmaddsub231pd %%ymm2, %%ymm9, %%ymm3 \n\t" \
                      "vfmaddsub231pd %%ymm0, %%ymm10, %%ymm4 \n\t" \
                      "vfmaddsub231pd %%ymm1, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_su3_multiply_pair_dble(u) \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm3 \n\t" \
                      "vbroadcastsd %1, %%ymm6 \n\t" \
                      "vbroadcastsd %2, %%ymm4 \n\t" \
                      "vbroadcastsd %3, %%ymm7 \n\t" \
                      "vbroadcastsd %4, %%ymm5 \n\t" \
                      "vbroadcastsd %5, %%ymm8" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c22.re), \
                      "m" ((u).c31.re), \
                      "m" ((u).c32.re) \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vmulpd %%ymm0, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %%ymm1, %%ymm6, %%ymm6 \n\t" \
                      "vmulpd %%ymm0, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %%ymm1, %%ymm7, %%ymm7 \n\t" \
                      "vmulpd %%ymm0, %%ymm5, %%ymm5 \n\t" \
                      "vmulpd %%ymm1, %%ymm8, %%ymm8 \n\t" \
                      "vaddpd %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddpd %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm6 \n\t" \
                      "vbroadcastsd %4, %%ymm7 \n\t" \
                      "vbroadcastsd %5, %%ymm8 \n\t" \
                      "vpermilpd $0x5, %%ymm0, %%ymm0 \n\t" \
                      : \
                      : \
                      "m" ((u).c13.re), \
                      "m" ((u).c21.im), \
                      "m" ((u).c33.re), \
                      "m" ((u).c11.im), \
                      "m" ((u).c23.re), \
                      "m" ((u).c31.im) \
                      : \
                      "xmm0", "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vmulpd %%ymm2, %%ymm9, %%ymm9 \n\t" \
                      "vmulpd %%ymm0, %%ymm10, %%ymm10 \n\t" \
                      "vmulpd %%ymm2, %%ymm11, %%ymm11 \n\t" \
                      "vmulpd %%ymm0, %%ymm6, %%ymm6 \n\t" \
                      "vmulpd %%ymm2, %%ymm7, %%ymm7 \n\t" \
                      "vmulpd %%ymm0, %%ymm8, %%ymm8 \n\t" \
                      "vaddpd %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubpd %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm11, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubpd %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddpd %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubpd %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8",  \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm6 \n\t" \
                      "vbroadcastsd %4, %%ymm7 \n\t" \
                      "vbroadcastsd %5, %%ymm8 \n\t" \
                      "vpermilpd $0x5, %%ymm1, %%ymm1 \n\t" \
                      "vpermilpd $0x5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      "m" ((u).c12.im), \
                      "m" ((u).c23.im), \
                      "m" ((u).c32.im), \
                      "m" ((u).c13.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c33.im) \
                      : \
                      "xmm1", "xmm2", "xmm6", "xmm7", \
                      "xmm8", "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vmulpd %%ymm1, %%ymm9, %%ymm9 \n\t" \
                      "vmulpd %%ymm2, %%ymm10, %%ymm10 \n\t" \
                      "vmulpd %%ymm1, %%ymm11, %%ymm11 \n\t" \
                      "vmulpd %%ymm2, %%ymm6, %%ymm6 \n\t" \
                      "vmulpd %%ymm1, %%ymm7, %%ymm7 \n\t" \
                      "vmulpd %%ymm2, %%ymm8, %%ymm8 \n\t" \
                      "vaddsubpd %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubpd %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubpd %%ymm11, %%ymm5, %%ymm5 \n\t" \
                      "vaddsubpd %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubpd %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubpd %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8",  \
                      "xmm9", "xmm10", "xmm11")

#endif

/*
* Multiplies a pair sl,sh of su3 vectors by an su3 matrix u^dagger, assuming
* sl and sh are in the low and high lanes of ymm0,..,ymm2. On output the
* result is in ymm3,..,ymm5 and all registers are changed.
*/

#if (defined FMA3)

#define _avx_su3_inverse_multiply_pair_dble(u) \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm0, %%ymm6 \n\t" \
                      "vpermilpd $0x5, %%ymm1, %%ymm7 \n\t" \
                      "vpermilpd $0x5, %%ymm2, %%ymm8 \n\t" \
                      "vbroadcastsd %0, %%ymm12 \n\t" \
                      "vbroadcastsd %1, %%ymm13 \n\t" \
                      "vbroadcastsd %2, %%ymm14 \n\t" \
                      "vbroadcastsd %3, %%ymm9 \n\t" \
                      "vbroadcastsd %4, %%ymm10 \n\t" \
                      "vbroadcastsd %5, %%ymm11" \
                      : \
                      : \
                      "m" ((u).c11.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c33.im), \
                      "m" ((u).c11.re), \
                      "m" ((u).c22.re), \
                      "m" ((u).c33.re) \
                      : \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vmulpd %%ymm6, %%ymm12, %%ymm3 \n\t" \
                      "vmulpd %%ymm7, %%ymm13, %%ymm4 \n\t" \
                      "vmulpd %%ymm8, %%ymm14, %%ymm5 \n\t" \
                      "vfmsubadd231pd %%ymm0, %%ymm9, %%ymm3 \n\t" \
                      "vfmsubadd231pd %%ymm1, %%ymm10, %%ymm4 \n\t" \
                      "vfmsubadd231pd %%ymm2, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm12 \n\t" \
                      "vbroadcastsd %1, %%ymm13 \n\t" \
                      "vbroadcastsd %2, %%ymm14 \n\t" \
                      "vbroadcastsd %3, %%ymm9 \n\t" \
                      "vbroadcastsd %4, %%ymm10 \n\t" \
                      "vbroadcastsd %5, %%ymm11" \
                      : \
                      : \
                      "m" ((u).c21.im), \
                      "m" ((u).c32.im), \
                      "m" ((u).c13.im), \
                      "m" ((u).c21.re), \
                      "m" ((u).c32.re), \
                      "m" ((u).c13.re) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vfmsubadd231pd %%ymm7, %%ymm12, %%ymm3 \n\t" \
                      "vfmsubadd231pd %%ymm8, %%ymm13, %%ymm4 \n\t" \
                      "vfmsubadd231pd %%ymm6, %%ymm14, %%ymm5 \n\t" \
                      "vfmsubadd231pd %%ymm1, %%ymm9, %%ymm3 \n\t" \
                      "vfmsubadd231pd %%ymm2, %%ymm10, %%ymm4 \n\t" \
                      "vfmsubadd231pd %%ymm0, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm12 \n\t" \
                      "vbroadcastsd %1, %%ymm13 \n\t" \
                      "vbroadcastsd %2, %%ymm14 \n\t" \
                      "vbroadcastsd %3, %%ymm9 \n\t" \
                      "vbroadcastsd %4, %%ymm10 \n\t" \
                      "vbroadcastsd %5, %%ymm11" \
                      : \
                      : \
                      "m" ((u).c31.im), \
                      "m" ((u).c12.im), \
                      "m" ((u).c23.im), \
                      "m" ((u).c31.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c23.re) \
                      : \
                      "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vfmsubadd231pd %%ymm8, %%ymm12, %%ymm3 \n\t" \
                      "vfmsubadd231pd %%ymm6, %%ymm13, %%ymm4 \n\t" \
                      "vfmsubadd231pd %%ymm7, %%ymm14, %%ymm5 \n\t" \
                      "vfmsubadd231pd %%ymm2, %%ymm9, %%ymm3 \n\t" \
                      "vfmsubadd231pd %%ymm0, %%ymm10, %%ymm4 \n\t" \
                      "vfmsubadd231pd %%ymm1, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_su3_inverse_multiply_pair_dble(u) \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm3 \n\t" \
                      "vbroadcastsd %1, %%ymm6 \n\t" \
                      "vbroadcastsd %2, %%ymm4 \n\t" \
                      "vbroadcastsd %3, %%ymm7 \n\t" \
                      "vbroadcastsd %4, %%ymm5 \n\t" \
                      "vbroadcastsd %5, %%ymm8" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c22.re), \
                      "m" ((u).c13.re), \
                      "m" ((u).c23.re) \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm6", \
                      "xmm7", "xmm8"); \
__asm__ __volatile__ ("vxorpd %%ymm15, %%ymm15, %%ymm15 \n\t" \
                      "vmulpd %%ymm0, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %%ymm1, %%ymm6, %%ymm6 \n\t" \
                      "vmulpd %%ymm0, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %%ymm1, %%ymm7, %%ymm7 \n\t" \
                      "vmulpd %%ymm0, %%ymm5, %%ymm5 \n\t" \
                      "vmulpd %%ymm1, %%ymm8, %%ymm8 \n\t" \
                      "vaddpd %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vsubpd %%ymm0, %%ymm15, %%ymm0 \n\t" \
                      "vaddpd %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm8, %%ymm5, %%ymm5 \n\t" \
                      "vpermilpd $0x5, %%ymm0, %%ymm0" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8", "xmm15"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm9 \n\t" \
                      "vbroadcastsd %1, %%ymm10 \n\t" \
                      "vbroadcastsd %2, %%ymm11 \n\t" \
                      "vbroadcastsd %3, %%ymm12 \n\t" \
                      "vbroadcastsd %4, %%ymm13 \n\t" \
                      "vbroadcastsd %5, %%ymm14" \
                      : \
                      : \
                      "m" ((u).c31.re), \
                      "m" ((u).c12.im), \
                      "m" ((u).c33.re), \
                      "m" ((u).c11.im), \
                      "m" ((u).c32.re), \
                      "m" ((u).c13.im) \
                      : \
                      "xmm9", "xmm10", "xmm11", "xmm12", \
                      "xmm13", "xmm14"); \
__asm__ __volatile__ ("vsubpd %%ymm1, %%ymm15, %%ymm1 \n\t" \
                      "vmulpd %%ymm2, %%ymm9, %%ymm9 \n\t" \
                      "vmulpd %%ymm0, %%ymm10, %%ymm10 \n\t" \
                      "vmulpd %%ymm2, %%ymm11, %%ymm11 \n\t" \
                      "vmulpd %%ymm0, %%ymm12, %%ymm12 \n\t" \
                      "vaddpd %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %%ymm2, %%ymm13, %%ymm13 \n\t" \
                      "vpermilpd $0x5, %%ymm1, %%ymm1 \n\t" \
                      "vaddsubpd %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vsubpd %%ymm2, %%ymm15, %%ymm2 \n\t" \
                      "vaddpd %%ymm11, %%ymm5, %%ymm5 \n\t" \
                      "vmulpd %%ymm0, %%ymm14, %%ymm14 \n\t" \
                      "vpermilpd $0x5, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm1", "xmm2", "xmm3", "xmm4", \
                      "xmm5", "xmm9", "xmm10", "xmm11", \
                      "xmm12", "xmm13", "xmm14"); \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm6 \n\t" \
                      "vbroadcastsd %1, %%ymm7 \n\t" \
                      "vbroadcastsd %2, %%ymm8 \n\t" \
                      "vbroadcastsd %3, %%ymm9 \n\t" \
                      "vbroadcastsd %4, %%ymm10 \n\t" \
                      "vbroadcastsd %5, %%ymm11" \
                      : \
                      : \
                      "m" ((u).c21.im), \
                      "m" ((u).c32.im), \
                      "m" ((u).c23.im), \
                      "m" ((u).c31.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c33.im) \
                      : \
                      "xmm6", "xmm7", "xmm8", "xmm9", \
                      "xmm10", "xmm11"); \
__asm__ __volatile__ ("vmulpd %%ymm1, %%ymm6, %%ymm6 \n\t" \
                      "vaddsubpd %%ymm12, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %%ymm2, %%ymm7, %%ymm7 \n\t" \
                      "vaddpd %%ymm13, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %%ymm1, %%ymm8, %%ymm8 \n\t" \
                      "vaddsubpd %%ymm14, %%ymm5, %%ymm5 \n\t" \
                      "vmulpd %%ymm2, %%ymm9, %%ymm9 \n\t" \
                      "vaddsubpd %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %%ymm1, %%ymm10, %%ymm10 \n\t" \
                      "vaddsubpd %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %%ymm2, %%ymm11, %%ymm11 \n\t" \
                      "vaddsubpd %%ymm8, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", "xmm6", \
                      "xmm7", "xmm8", "xmm9", "xmm10", \
                      "xmm11"); \
__asm__ __volatile__ ("vaddsubpd %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddsubpd %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddsubpd %%ymm11, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

#endif

/******************************************************************************
*
*  Macros for double precision Dirac spinors in linear order.
*
******************************************************************************/

/*
* Loads the spinor s to the registers ymm0,..,ymm5 in linear order.
*/

#define _avx_spinor_load_dble(s) \
__asm__ __volatile__ ("vmovapd %0, %%ymm0 \n\t" \
                      "vmovapd %2, %%ymm1 \n\t" \
                      "vmovapd %4, %%ymm2" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm0", "xmm1", "xmm2"); \
__asm__ __volatile__ ("vmovapd %0, %%ymm3 \n\t" \
                      "vmovapd %2, %%ymm4 \n\t" \
                      "vmovapd %4, %%ymm5" \
                      : \
                      : \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2), \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Loads the spinor s to the registers ymm6,..,ymm11 in linear order.
*/

#define _avx_spinor_load_up_dble(s) \
__asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t" \
                      "vmovapd %2, %%ymm7 \n\t" \
                      "vmovapd %4, %%ymm8" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vmovapd %0, %%ymm9 \n\t" \
                      "vmovapd %2, %%ymm10 \n\t" \
                      "vmovapd %4, %%ymm11" \
                      : \
                      : \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2), \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm9", "xmm10", "xmm11")

/*
* Stores the registers ymm0,..,ymm5 to the spinor s in linear order.
*/

#define _avx_spinor_store_dble(s) \
__asm__ __volatile__ ("vmovapd %%ymm0, %0 \n\t" \
                      "vmovapd %%ymm1, %2 \n\t" \
                      "vmovapd %%ymm2, %4" \
                      : \
                      "=m" ((s).c1.c1), \
                      "=m" ((s).c1.c2), \
                      "=m" ((s).c1.c3), \
                      "=m" ((s).c2.c1), \
                      "=m" ((s).c2.c2), \
                      "=m" ((s).c2.c3)); \
__asm__ __volatile__ ("vmovapd %%ymm3, %0 \n\t" \
                      "vmovapd %%ymm4, %2 \n\t" \
                      "vmovapd %%ymm5, %4" \
                      : \
                      "=m" ((s).c3.c1), \
                      "=m" ((s).c3.c2), \
                      "=m" ((s).c3.c3), \
                      "=m" ((s).c4.c1), \
                      "=m" ((s).c4.c2), \
                      "=m" ((s).c4.c3))

/*
* Stores the registers ymm6,..,ymm11 to the spinor s in linear order.
*/

#define _avx_spinor_store_up_dble(s) \
__asm__ __volatile__ ("vmovapd %%ymm6, %0 \n\t" \
                      "vmovapd %%ymm7, %2 \n\t" \
                      "vmovapd %%ymm8, %4" \
                      : \
                      "=m" ((s).c1.c1), \
                      "=m" ((s).c1.c2), \
                      "=m" ((s).c1.c3), \
                      "=m" ((s).c2.c1), \
                      "=m" ((s).c2.c2), \
                      "=m" ((s).c2.c3)); \
__asm__ __volatile__ ("vmovapd %%ymm9, %0 \n\t" \
                      "vmovapd %%ymm10, %2 \n\t" \
                      "vmovapd %%ymm11, %4" \
                      : \
                      "=m" ((s).c3.c1), \
                      "=m" ((s).c3.c2), \
                      "=m" ((s).c3.c3), \
                      "=m" ((s).c4.c1), \
                      "=m" ((s).c4.c2), \
                      "=m" ((s).c4.c3))

/*
* Loads (z.re,z.re,z.re,z.re) to ymm12 and (-z.im,z.im,-z.im,z.im) to ymm13.
*/

#define _avx_load_cmplx_dble(z) \
__asm__ __volatile__ ("vxorpd %%ymm13, %%ymm13, %%ymm13 \n\t" \
                      "vbroadcastsd %0, %%ymm12 \n\t" \
                      "vaddsubpd %%ymm12, %%ymm13, %%ymm13 \n\t" \
                      "vbroadcastsd %1, %%ymm12" \
                      : \
                      : \
                      "m" ((z).im), \
                      "m" ((z).re) \
                      : \
                      "xmm12", "xmm13")

/*
* Loads (z.re,z.re,z.re,z.re) to ymm14 and (-z.im,z.im,-z.im,z.im) to ymm15.
*/

#define _avx_load_cmplx_up_dble(z) \
__asm__ __volatile__ ("vxorpd %%ymm15, %%ymm15, %%ymm15 \n\t" \
                      "vbroadcastsd %0, %%ymm14 \n\t" \
                      "vaddsubpd %%ymm14, %%ymm15, %%ymm15 \n\t" \
                      "vbroadcastsd %1, %%ymm14" \
                      : \
                      : \
                      "m" ((z).im), \
                      "m" ((z).re) \
                      : \
                      "xmm14", "xmm15")

/*
* Multiplies the spinor s by the complex number z and assigns the result to
* ymm0,..,ymm5, assuming z was loaded using _avx_load_cmplx_dble(z). The
* registers ymm6,..,ymm11 are used as workspace.
*/

#if (defined FMA3)

#define _avx_mulc_spinor_dble(s) \
_avx_spinor_load_dble(s); \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm0, %%ymm6 \n\t" \
                      "vpermilpd $0x5, %%ymm1, %%ymm7 \n\t" \
                      "vpermilpd $0x5, %%ymm2, %%ymm8 \n\t" \
                      "vpermilpd $0x5, %%ymm3, %%ymm9 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm10 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm11 \n\t" \
                      "vmulpd %%ymm12, %%ymm0, %%ymm0 \n\t" \
                      "vmulpd %%ymm12, %%ymm1, %%ymm1 \n\t" \
                      "vmulpd %%ymm12, %%ymm2, %%ymm2 \n\t" \
                      "vmulpd %%ymm12, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %%ymm12, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %%ymm12, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vfmadd231pd %%ymm13, %%ymm6, %%ymm0 \n\t" \
                      "vfmadd231pd %%ymm13, %%ymm7, %%ymm1 \n\t" \
                      "vfmadd231pd %%ymm13, %%ymm8, %%ymm2 \n\t" \
                      "vfmadd231pd %%ymm13, %%ymm9, %%ymm3 \n\t" \
                      "vfmadd231pd %%ymm13, %%ymm10, %%ymm4 \n\t" \
                      "vfmadd231pd %%ymm13, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_mulc_spinor_dble(s) \
_avx_spinor_load_dble(s); \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm0, %%ymm6 \n\t" \
                      "vpermilpd $0x5, %%ymm1, %%ymm7 \n\t" \
                      "vpermilpd $0x5, %%ymm2, %%ymm8 \n\t" \
                      "vpermilpd $0x5, %%ymm3, %%ymm9 \n\t" \
                      "vpermilpd $0x5, %%ymm4, %%ymm10 \n\t" \
                      "vpermilpd $0x5, %%ymm5, %%ymm11 \n\t" \
                      "vmulpd %%ymm12, %%ymm0, %%ymm0 \n\t" \
                      "vmulpd %%ymm13, %%ymm6, %%ymm6 \n\t" \
                      "vmulpd %%ymm12, %%ymm1, %%ymm1 \n\t" \
                      "vmulpd %%ymm13, %%ymm7, %%ymm7 \n\t" \
                      "vmulpd %%ymm12, %%ymm2, %%ymm2 \n\t" \
                      "vmulpd %%ymm13, %%ymm8, %%ymm8" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vmulpd %%ymm12, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %%ymm13, %%ymm9, %%ymm9 \n\t" \
                      "vmulpd %%ymm12, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %%ymm13, %%ymm10, %%ymm10 \n\t" \
                      "vmulpd %%ymm12, %%ymm5, %%ymm5 \n\t" \
                      "vmulpd %%ymm13, %%ymm11, %%ymm11 \n\t" \
                      "vaddpd %%ymm6, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm7, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm8, %%ymm2, %%ymm2 \n\t" \
                      "vaddpd %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddpd %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm11, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm9", "xmm10", "xmm11")

#endif

/*
* Multiplies the spinor s by the complex number z and adds the result to
* ymm0,..,ymm5, assuming z was loaded using _avx_load_cmplx_up_dble(z). The
* registers ymm6,..,ymm11 are used as workspace.
*/

#if (defined FMA3)

#define _avx_mulc_spinor_add_dble(s) \
_avx_spinor_load_up_dble(s); \
__asm__ __volatile__ ("vfmadd231pd %%ymm14, %%ymm6, %%ymm0 \n\t" \
                      "vfmadd231pd %%ymm14, %%ymm7, %%ymm1 \n\t" \
                      "vfmadd231pd %%ymm14, %%ymm8, %%ymm2 \n\t" \
                      "vfmadd231pd %%ymm14, %%ymm9, %%ymm3 \n\t" \
                      "vfmadd231pd %%ymm14, %%ymm10, %%ymm4 \n\t" \
                      "vfmadd231pd %%ymm14, %%ymm11, %%ymm5 \n\t" \
                      "vpermilpd $0x5, %%ymm6, %%ymm6 \n\t"  \
                      "vpermilpd $0x5, %%ymm7, %%ymm7 \n\t" \
                      "vpermilpd $0x5, %%ymm8, %%ymm8 \n\t" \
                      "vpermilpd $0x5, %%ymm9, %%ymm9 \n\t" \
                      "vpermilpd $0x5, %%ymm10, %%ymm10 \n\t" \
                      "vpermilpd $0x5, %%ymm11, %%ymm11 \n\t" \
                      "vfmadd231pd %%ymm15, %%ymm6, %%ymm0 \n\t" \
                      "vfmadd231pd %%ymm15, %%ymm7, %%ymm1 \n\t" \
                      "vfmadd231pd %%ymm15, %%ymm8, %%ymm2 \n\t" \
                      "vfmadd231pd %%ymm15, %%ymm9, %%ymm3 \n\t" \
                      "vfmadd231pd %%ymm15, %%ymm10, %%ymm4 \n\t" \
                      "vfmadd231pd %%ymm15, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11")

#else

#define _avx_mulc_spinor_add_dble(s) \
__asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t" \
                      "vmovapd %2, %%ymm7 \n\t" \
                      "vmovapd %4, %%ymm8" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm6, %%ymm9 \n\t" \
                      "vpermilpd $0x5, %%ymm7, %%ymm10 \n\t" \
                      "vpermilpd $0x5, %%ymm8, %%ymm11 \n\t" \
                      "vmulpd %%ymm14, %%ymm6, %%ymm6 \n\t" \
                      "vmulpd %%ymm15, %%ymm9, %%ymm9 \n\t" \
                      "vmulpd %%ymm14, %%ymm7, %%ymm7 \n\t" \
                      "vmulpd %%ymm15, %%ymm10, %%ymm10 \n\t" \
                      "vmulpd %%ymm14, %%ymm8, %%ymm8 \n\t" \
                      "vmulpd %%ymm15, %%ymm11, %%ymm11 \n\t" \
                      "vaddpd %%ymm6, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm7, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm8, %%ymm2, %%ymm2 \n\t" \
                      "vaddpd %%ymm9, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm10, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm11, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t" \
                      "vmovapd %2, %%ymm7 \n\t" \
                      "vmovapd %4, %%ymm8" \
                      : \
                      : \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2), \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("vpermilpd $0x5, %%ymm6, %%ymm9 \n\t" \
                      "vpermilpd $0x5, %%ymm7, %%ymm10 \n\t" \
                      "vpermilpd $0x5, %%ymm8, %%ymm11 \n\t" \
                      "vmulpd %%ymm14, %%ymm6, %%ymm6 \n\t" \
                      "vmulpd %%ymm15, %%ymm9, %%ymm9 \n\t" \
                      "vmulpd %%ymm14, %%ymm7, %%ymm7 \n\t" \
                      "vmulpd %%ymm15, %%ymm10, %%ymm10 \n\t" \
                      "vmulpd %%ymm14, %%ymm8, %%ymm8 \n\t" \
                      "vmulpd %%ymm15, %%ymm11, %%ymm11 \n\t" \
                      "vaddpd %%ymm6, %%ymm3, %%ymm3 \n\t" \
                      "vaddpd %%ymm7, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm8, %%ymm5, %%ymm5 \n\t" \
                      "vaddpd %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddpd %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm11, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11")

#endif

/*
* Broadcasts the real number c to ymm12 and ymm13.
*/

#define _avx_load_real_dble(c) \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm12 \n\t" \
                      "vbroadcastsd %0, %%ymm13" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm12", "xmm13")

/*
* Broadcasts the real number c to ymm14 and ymm15.
*/

#define _avx_load_real_up_dble(c) \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm14 \n\t" \
                      "vbroadcastsd %0, %%ymm15" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm14", "xmm15")

/*
* Multiplies the spinor s by the real number c and assigns the result to
* ymm0,..,ymm5, assuming c was loaded using _avx_load_real_dble(c).
*/

#define _avx_mulr_spinor_dble(s) \
_avx_spinor_load_dble(s); \
__asm__ __volatile__ ("vmulpd %%ymm12, %%ymm0, %%ymm0 \n\t" \
                      "vmulpd %%ymm13, %%ymm1, %%ymm1 \n\t" \
                      "vmulpd %%ymm12, %%ymm2, %%ymm2 \n\t" \
                      "vmulpd %%ymm13, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %%ymm12, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %%ymm13, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")


/*
* Multiplies the spinor s by the real number c and adds the result to
* ymm0,..,ymm5, assuming c was loaded using _avx_load_real_up_dble(c).
* The registers ymm6,..,ymm11 are used as workspace.
*/

#if (defined FMA3)

#define _avx_mulr_spinor_add_dble(s) \
_avx_spinor_load_up_dble(s); \
__asm__ __volatile__ ("vfmadd231pd %%ymm14, %%ymm6, %%ymm0 \n\t" \
                      "vfmadd231pd %%ymm15, %%ymm7, %%ymm1 \n\t" \
                      "vfmadd231pd %%ymm14, %%ymm8, %%ymm2 \n\t" \
                      "vfmadd231pd %%ymm15, %%ymm9, %%ymm3 \n\t" \
                      "vfmadd231pd %%ymm14, %%ymm10, %%ymm4 \n\t" \
                      "vfmadd231pd %%ymm15, %%ymm11, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#else

#define _avx_mulr_spinor_add_dble(s) \
_avx_spinor_load_up_dble(s); \
__asm__ __volatile__ ("vmulpd %%ymm14, %%ymm6, %%ymm6 \n\t" \
                      "vmulpd %%ymm15, %%ymm7, %%ymm7 \n\t" \
                      "vmulpd %%ymm14, %%ymm8, %%ymm8 \n\t" \
                      "vmulpd %%ymm15, %%ymm9, %%ymm9 \n\t" \
                      "vmulpd %%ymm14, %%ymm10, %%ymm10 \n\t" \
                      "vmulpd %%ymm15, %%ymm11, %%ymm11 \n\t" \
                      "vaddpd %%ymm6, %%ymm0, %%ymm0 \n\t" \
                      "vaddpd %%ymm7, %%ymm1, %%ymm1 \n\t" \
                      "vaddpd %%ymm8, %%ymm2, %%ymm2 \n\t" \
                      "vaddpd %%ymm9, %%ymm3, %%ymm3 \n\t" \
                      "vaddpd %%ymm10, %%ymm4, %%ymm4 \n\t" \
                      "vaddpd %%ymm11, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11")

#endif

#endif
