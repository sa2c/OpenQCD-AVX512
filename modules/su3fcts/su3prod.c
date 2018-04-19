
/*******************************************************************************
*
* File su3prod.c
*
* Copyright (C) 2005, 2009-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Products of double-precision 3x3 matrices.
*
* The externally accessible functions are
*
*   void su3xsu3(su3_dble *u,su3_dble *v,su3_dble *w)
*     Computes w=u*v assuming that w is different from u.
*
*   void su3dagxsu3(su3_dble *u,su3_dble *v,su3_dble *w)
*     Computes w=u^dag*v assuming that w is different from u.
*
*   void su3xsu3dag(su3_dble *u,su3_dble *v,su3_dble *w)
*     Computes w=u*v^dag assuming that w is different from u and v.
*
*   void su3dagxsu3dag(su3_dble *u,su3_dble *v,su3_dble *w)
*     Computes w=u^dag*v^dag assuming that w is different from u and v.
*
*   void su3xu3alg(su3_dble *u,u3_alg_dble *X,su3_dble *v)
*     Computes v=u*X assuming that v is different from u.
*
*   void su3dagxu3alg(su3_dble *u,u3_alg_dble *X,su3_dble *v)
*     Computes v=u^dag*X assuming that v is different from u.
*
*   void u3algxsu3(u3_alg_dble *X,su3_dble *u,su3_dble *v)
*     Computes v=X*u assuming that v is different from u.
*
*   void u3algxsu3dag(u3_alg_dble *X,su3_dble *u,su3_dble *v)
*     Computes v=X*u^dag assuming that v is different from u.
*
*   double prod2su3alg(su3_dble *u,su3_dble *v,su3_alg_dble *X)
*     Computes the product w=u*v and assigns its traceless antihermitian
*     part (1/2)*[w-w^dag-(1/3)*tr{w-w^dag}] to X. The program returns
*     the real part of tr{w}.
*
*   void prod2u3alg(su3_dble *u,su3_dble *v,u3_alg_dble *X)
*     Computes the product w=u*v and assigns w-w^dag to X.
*
*   void rotate_su3alg(su3_dble *u,su3_alg_dble *X)
*     Replaces X by u*X*u^dag. The matrix u must be unitary but its
*     determinant may be different from 1.
*
* Notes:
*
* Unless stated otherwise, the matrices of type su3_dble are not assumed to
* be unitary or unimodular. They are just treated as general 3x3 complex
* matrices and the operations are applied to them as described.
*
* The elements X of the Lie algebra of U(3) are antihermitian 3x3 matrices
* that are represented by structures X with real entries X.c1,...,X.c9
* through
*
*  X_11=i*X.c1, X_22=i*X.c2, X_33=i*X.c3,
*
*  X_12=X.c4+i*X.c5, X_13=X.c6+i*X.c7, X_23=X.c8+i*X.c9
*
* The type su3_alg_dble [which represents elements of the Lie algebra of SU(3)]
* is described in the file linalg/liealg.c.
*
* If SSE2 or AVX instructions are used, all su3_dble and su3_alg_dble matrices
* are assumed to be aligned to 16 byte boundaries.
*
*******************************************************************************/

#define SU3PROD_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "su3fcts.h"

#if (defined AVX)
#include "avx.h"

static const sse_double c0={0.5,0.5},c1={-1.0/3.0,-1.0/3.0};
static su3_dble uX ALIGNED16;
static double tr ALIGNED8;


static void su3xsu3vec(su3_dble *u)
{
   _avx_su3_multiply_dble(*u);
}


static void su3xsu3vec_pair(su3_dble *u)
{
   _avx_su3_multiply_pair_dble(*u);
}


static void su3dagxsu3vec(su3_dble *u)
{
   _avx_su3_inverse_multiply_dble(*u);
}


static void su3dagxsu3vec_pair(su3_dble *u)
{
   _avx_su3_inverse_multiply_pair_dble(*u);
}


void su3xsu3(su3_dble *u,su3_dble *v,su3_dble *w)
{
   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2 \n\t"
                         "vinsertf128 $0x1, %3, %%ymm0, %%ymm0 \n\t"
                         "vinsertf128 $0x1, %4, %%ymm1, %%ymm1 \n\t"
                         "vinsertf128 $0x1, %5, %%ymm2, %%ymm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c21),
                         "m" ((*v).c31),
                         "m" ((*v).c12),
                         "m" ((*v).c22),
                         "m" ((*v).c32)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec_pair(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2 \n\t"
                         "vextractf128 $0x1, %%ymm3, %3 \n\t"
                         "vextractf128 $0x1, %%ymm4, %4 \n\t"
                         "vextractf128 $0x1, %%ymm5, %5"
                         :
                         "=m" ((*w).c11),
                         "=m" ((*w).c21),
                         "=m" ((*w).c31),
                         "=m" ((*w).c12),
                         "=m" ((*w).c22),
                         "=m" ((*w).c32));

   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c13),
                         "m" ((*v).c23),
                         "m" ((*v).c33)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2"
                         :
                         "=m" ((*w).c13),
                         "=m" ((*w).c23),
                         "=m" ((*w).c33));

   _avx_zeroupper();
}


void su3dagxsu3(su3_dble *u,su3_dble *v,su3_dble *w)
{
   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2 \n\t"
                         "vinsertf128 $0x1, %3, %%ymm0, %%ymm0 \n\t"
                         "vinsertf128 $0x1, %4, %%ymm1, %%ymm1 \n\t"
                         "vinsertf128 $0x1, %5, %%ymm2, %%ymm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c21),
                         "m" ((*v).c31),
                         "m" ((*v).c12),
                         "m" ((*v).c22),
                         "m" ((*v).c32)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3dagxsu3vec_pair(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2 \n\t"
                         "vextractf128 $0x1, %%ymm3, %3 \n\t"
                         "vextractf128 $0x1, %%ymm4, %4 \n\t"
                         "vextractf128 $0x1, %%ymm5, %5"
                         :
                         "=m" ((*w).c11),
                         "=m" ((*w).c21),
                         "=m" ((*w).c31),
                         "=m" ((*w).c12),
                         "=m" ((*w).c22),
                         "=m" ((*w).c32));

   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c13),
                         "m" ((*v).c23),
                         "m" ((*v).c33)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3dagxsu3vec(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2"
                         :
                         "=m" ((*w).c13),
                         "=m" ((*w).c23),
                         "=m" ((*w).c33));

   _avx_zeroupper();
}


void su3xsu3dag(su3_dble *u,su3_dble *v,su3_dble *w)
{
   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2 \n\t"
                         "vinsertf128 $0x1, %3, %%ymm0, %%ymm0 \n\t"
                         "vinsertf128 $0x1, %4, %%ymm1, %%ymm1 \n\t"
                         "vinsertf128 $0x1, %5, %%ymm2, %%ymm2 \n\t"
                         "vmulpd %6, %%ymm0, %%ymm0 \n\t"
                         "vmulpd %6, %%ymm1, %%ymm1 \n\t"
                         "vmulpd %6, %%ymm2, %%ymm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c12),
                         "m" ((*v).c13),
                         "m" ((*v).c21),
                         "m" ((*v).c22),
                         "m" ((*v).c23),
                         "m" (_avx_sgn24_dble)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec_pair(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2 \n\t"
                         "vextractf128 $0x1, %%ymm3, %3 \n\t"
                         "vextractf128 $0x1, %%ymm4, %4 \n\t"
                         "vextractf128 $0x1, %%ymm5, %5"
                         :
                         "=m" ((*w).c11),
                         "=m" ((*w).c21),
                         "=m" ((*w).c31),
                         "=m" ((*w).c12),
                         "=m" ((*w).c22),
                         "=m" ((*w).c32));

   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2 \n\t"
                         "vmulpd %3, %%xmm0, %%xmm0 \n\t"
                         "vmulpd %3, %%xmm1, %%xmm1 \n\t"
                         "vmulpd %3, %%xmm2, %%xmm2"
                         :
                         :
                         "m" ((*v).c31),
                         "m" ((*v).c32),
                         "m" ((*v).c33),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2"
                         :
                         "=m" ((*w).c13),
                         "=m" ((*w).c23),
                         "=m" ((*w).c33));

   _avx_zeroupper();
}


void su3dagxsu3dag(su3_dble *u,su3_dble *v,su3_dble *w)
{
   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2 \n\t"
                         "vinsertf128 $0x1, %3, %%ymm0, %%ymm0 \n\t"
                         "vinsertf128 $0x1, %4, %%ymm1, %%ymm1 \n\t"
                         "vinsertf128 $0x1, %5, %%ymm2, %%ymm2 \n\t"
                         "vmulpd %6, %%ymm0, %%ymm0 \n\t"
                         "vmulpd %6, %%ymm1, %%ymm1 \n\t"
                         "vmulpd %6, %%ymm2, %%ymm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c12),
                         "m" ((*v).c13),
                         "m" ((*v).c21),
                         "m" ((*v).c22),
                         "m" ((*v).c23),
                         "m" (_avx_sgn24_dble)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3dagxsu3vec_pair(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2 \n\t"
                         "vextractf128 $0x1, %%ymm3, %3 \n\t"
                         "vextractf128 $0x1, %%ymm4, %4 \n\t"
                         "vextractf128 $0x1, %%ymm5, %5"
                         :
                         "=m" ((*w).c11),
                         "=m" ((*w).c21),
                         "=m" ((*w).c31),
                         "=m" ((*w).c12),
                         "=m" ((*w).c22),
                         "=m" ((*w).c32));

   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2 \n\t"
                         "vmulpd %3, %%xmm0, %%xmm0 \n\t"
                         "vmulpd %3, %%xmm1, %%xmm1 \n\t"
                         "vmulpd %3, %%xmm2, %%xmm2"
                         :
                         :
                         "m" ((*v).c31),
                         "m" ((*v).c32),
                         "m" ((*v).c33),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3dagxsu3vec(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2"
                         :
                         "=m" ((*w).c13),
                         "=m" ((*w).c23),
                         "=m" ((*w).c33));

   _avx_zeroupper();
}


void su3xu3alg(su3_dble *u,u3_alg_dble *X,su3_dble *v)
{
   __asm__ __volatile__ ("vxorpd %%ymm5, %%ymm5, %%ymm5 \n\t"
                         "vbroadcastf128 %8, %%ymm6 \n\t"
                         "vmovupd %0, %%xmm1 \n\t"
                         "vmovupd %2, %%xmm2 \n\t"
                         "vmovhpd %4, %%xmm5, %%xmm0 \n\t"
                         "vmovhpd %5, %%xmm5, %%xmm3 \n\t"
                         "vinsertf128 $0x1, %6, %%ymm2, %%ymm2 \n\t"
                         "vinsertf128 $0x1, %%xmm1, %%ymm0, %%ymm0 \n\t"
                         "vmulpd %%xmm6, %%xmm1, %%xmm1 \n\t"
                         "vmulpd %%ymm6, %%ymm2, %%ymm2 \n\t"
                         "vinsertf128 $0x1, %%xmm3, %%ymm1, %%ymm1"
                         :
                         :
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c1),
                         "m" ((*X).c2),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2",
                         "xmm3", "xmm5", "xmm6");

   su3xsu3vec_pair(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2 \n\t"
                         "vextractf128 $0x1, %%ymm3, %3 \n\t"
                         "vextractf128 $0x1, %%ymm4, %4 \n\t"
                         "vextractf128 $0x1, %%ymm5, %5"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c21),
                         "=m" ((*v).c31),
                         "=m" ((*v).c12),
                         "=m" ((*v).c22),
                         "=m" ((*v).c32));

   __asm__ __volatile__ ("vxorpd %%xmm2, %%xmm2, %%xmm2\n\t"
                         "vmovupd %0, %%xmm0 \n\t"
                         "vmovupd %2, %%xmm1 \n\t"
                         "vmovhpd %4, %%xmm2, %%xmm2"
                         :
                         :
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" ((*X).c3)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2"
                         :
                         "=m" ((*v).c13),
                         "=m" ((*v).c23),
                         "=m" ((*v).c33));

   _avx_zeroupper();
}


void su3dagxu3alg(su3_dble *u,u3_alg_dble *X,su3_dble *v)
{
   __asm__ __volatile__ ("vxorpd %%ymm5, %%ymm5, %%ymm5 \n\t"
                         "vbroadcastf128 %8, %%ymm6 \n\t"
                         "vmovupd %0, %%xmm1 \n\t"
                         "vmovupd %2, %%xmm2 \n\t"
                         "vmovhpd %4, %%xmm5, %%xmm0 \n\t"
                         "vmovhpd %5, %%xmm5, %%xmm3 \n\t"
                         "vinsertf128 $0x1, %6, %%ymm2, %%ymm2 \n\t"
                         "vinsertf128 $0x1, %%xmm1, %%ymm0, %%ymm0 \n\t"
                         "vmulpd %%xmm6, %%xmm1, %%xmm1 \n\t"
                         "vmulpd %%ymm6, %%ymm2, %%ymm2 \n\t"
                         "vinsertf128 $0x1, %%xmm3, %%ymm1, %%ymm1"
                         :
                         :
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c1),
                         "m" ((*X).c2),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2",
                         "xmm3", "xmm5", "xmm6");

   su3dagxsu3vec_pair(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2 \n\t"
                         "vextractf128 $0x1, %%ymm3, %3 \n\t"
                         "vextractf128 $0x1, %%ymm4, %4 \n\t"
                         "vextractf128 $0x1, %%ymm5, %5"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c21),
                         "=m" ((*v).c31),
                         "=m" ((*v).c12),
                         "=m" ((*v).c22),
                         "=m" ((*v).c32));

   __asm__ __volatile__ ("vxorpd %%xmm2, %%xmm2, %%xmm2\n\t"
                         "vmovupd %0, %%xmm0 \n\t"
                         "vmovupd %2, %%xmm1 \n\t"
                         "vmovhpd %4, %%xmm2, %%xmm2"
                         :
                         :
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" ((*X).c3)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3dagxsu3vec(u);

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2"
                         :
                         "=m" ((*v).c13),
                         "=m" ((*v).c23),
                         "=m" ((*v).c33));

   _avx_zeroupper();
}


void u3algxsu3(u3_alg_dble *X,su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("vxorpd %%ymm5, %%ymm5, %%ymm5 \n\t"
                         "vbroadcastf128 %8, %%ymm6 \n\t"
                         "vmovupd %0, %%xmm1 \n\t"
                         "vmovupd %2, %%xmm2 \n\t"
                         "vmovhpd %4, %%xmm5, %%xmm0 \n\t"
                         "vmovhpd %5, %%xmm5, %%xmm3 \n\t"
                         "vinsertf128 $0x1, %6, %%ymm2, %%ymm2 \n\t"
                         "vinsertf128 $0x1, %%xmm1, %%ymm0, %%ymm0 \n\t"
                         "vmulpd %%xmm6, %%xmm1, %%xmm1 \n\t"
                         "vmulpd %%ymm6, %%ymm2, %%ymm2 \n\t"
                         "vinsertf128 $0x1, %%xmm3, %%ymm1, %%ymm1"
                         :
                         :
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c1),
                         "m" ((*X).c2),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2",
                         "xmm3", "xmm5", "xmm6");

   su3dagxsu3vec_pair(u);

   __asm__ __volatile__ ("vxorpd %%ymm0, %%ymm0, %%ymm0 \n\t"
                         "vaddsubpd %%ymm3, %%ymm0, %%ymm3 \n\t"
                         "vaddsubpd %%ymm4, %%ymm0, %%ymm4 \n\t"
                         "vaddsubpd %%ymm5, %%ymm0, %%ymm5 \n\t"
                         "vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2 \n\t"
                         "vextractf128 $0x1, %%ymm3, %3 \n\t"
                         "vextractf128 $0x1, %%ymm4, %4 \n\t"
                         "vextractf128 $0x1, %%ymm5, %5"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c12),
                         "=m" ((*v).c13),
                         "=m" ((*v).c21),
                         "=m" ((*v).c22),
                         "=m" ((*v).c23)
                         :
                         :
                         "xmm0", "xmm3", "xmm4", "xmm5");

   __asm__ __volatile__ ("vxorpd %%xmm2, %%xmm2, %%xmm2\n\t"
                         "vmovupd %0, %%xmm0 \n\t"
                         "vmovupd %2, %%xmm1 \n\t"
                         "vmovhpd %4, %%xmm2, %%xmm2"
                         :
                         :
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" ((*X).c3)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3dagxsu3vec(u);

   __asm__ __volatile__ ("vmulpd %3, %%xmm3, %%xmm3\n\t"
                         "vmulpd %3, %%xmm4, %%xmm4\n\t"
                         "vmulpd %3, %%xmm5, %%xmm5\n\t"
                         "vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2"
                         :
                         "=m" ((*v).c31),
                         "=m" ((*v).c32),
                         "=m" ((*v).c33)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm4", "xmm5");

   _avx_zeroupper();
}


void u3algxsu3dag(u3_alg_dble *X,su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("vxorpd %%ymm5, %%ymm5, %%ymm5 \n\t"
                         "vbroadcastf128 %8, %%ymm6 \n\t"
                         "vmovupd %0, %%xmm1 \n\t"
                         "vmovupd %2, %%xmm2 \n\t"
                         "vmovhpd %4, %%xmm5, %%xmm0 \n\t"
                         "vmovhpd %5, %%xmm5, %%xmm3 \n\t"
                         "vinsertf128 $0x1, %6, %%ymm2, %%ymm2 \n\t"
                         "vinsertf128 $0x1, %%xmm1, %%ymm0, %%ymm0 \n\t"
                         "vmulpd %%xmm6, %%xmm1, %%xmm1 \n\t"
                         "vmulpd %%ymm6, %%ymm2, %%ymm2 \n\t"
                         "vinsertf128 $0x1, %%xmm3, %%ymm1, %%ymm1"
                         :
                         :
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c1),
                         "m" ((*X).c2),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2",
                         "xmm3", "xmm5", "xmm6");

   su3xsu3vec_pair(u);

   __asm__ __volatile__ ("vxorpd %%ymm0, %%ymm0, %%ymm0 \n\t"
                         "vaddsubpd %%ymm3, %%ymm0, %%ymm3 \n\t"
                         "vaddsubpd %%ymm4, %%ymm0, %%ymm4 \n\t"
                         "vaddsubpd %%ymm5, %%ymm0, %%ymm5 \n\t"
                         "vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2 \n\t"
                         "vextractf128 $0x1, %%ymm3, %3 \n\t"
                         "vextractf128 $0x1, %%ymm4, %4 \n\t"
                         "vextractf128 $0x1, %%ymm5, %5"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c12),
                         "=m" ((*v).c13),
                         "=m" ((*v).c21),
                         "=m" ((*v).c22),
                         "=m" ((*v).c23)
                         :
                         :
                         "xmm0", "xmm3", "xmm4", "xmm5");

   __asm__ __volatile__ ("vxorpd %%xmm2, %%xmm2, %%xmm2\n\t"
                         "vmovupd %0, %%xmm0 \n\t"
                         "vmovupd %2, %%xmm1 \n\t"
                         "vmovhpd %4, %%xmm2, %%xmm2"
                         :
                         :
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" ((*X).c3)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec(u);

   __asm__ __volatile__ ("vmulpd %3, %%xmm3, %%xmm3\n\t"
                         "vmulpd %3, %%xmm4, %%xmm4\n\t"
                         "vmulpd %3, %%xmm5, %%xmm5\n\t"
                         "vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %1 \n\t"
                         "vmovapd %%xmm5, %2"
                         :
                         "=m" ((*v).c31),
                         "=m" ((*v).c32),
                         "=m" ((*v).c33)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm4", "xmm5");

   _avx_zeroupper();
}


double prod2su3alg(su3_dble *u,su3_dble *v,su3_alg_dble *X)
{
   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2 \n\t"
                         "vinsertf128 $0x1, %3, %%ymm0, %%ymm0 \n\t"
                         "vinsertf128 $0x1, %4, %%ymm1, %%ymm1 \n\t"
                         "vinsertf128 $0x1, %5, %%ymm2, %%ymm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c21),
                         "m" ((*v).c31),
                         "m" ((*v).c12),
                         "m" ((*v).c22),
                         "m" ((*v).c32)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec_pair(u);

   __asm__ __volatile__ ("vmovapd %%xmm5, %0 \n\t"
                         "vextractf128 $0x1, %%ymm3, %%xmm6 \n\t"
                         "vextractf128 $0x1, %%ymm4, %%xmm7 \n\t"
                         "vextractf128 $0x1, %%ymm5, %2 \n\t"
                         "vaddsubpd %%xmm4, %%xmm6, %%xmm6"
                         :
                         "=m" ((*X).c5),
                         "=m" ((*X).c6),
                         "=m" ((*X).c7),
                         "=m" ((*X).c8)
                         :
                         :
                         "xmm6", "xmm7");

   __asm__ __volatile__ ("vmovhpd %%xmm3, %0 \n\t"
                         "vmovhpd %%xmm7, %1 \n\t"
                         "vaddsd %%xmm7, %%xmm3, %%xmm15 \n\t"
                         "vmulpd %4, %%xmm6, %%xmm6 \n\t"
                         "vmovapd %%xmm6, %2"
                         :
                         "=m" ((*X).c1),
                         "=m" ((*X).c2),
                         "=m" ((*X).c3),
                         "=m" ((*X).c4)
                         :
                         "m" (c0)
                         :
                         "xmm6", "xmm15");

   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c13),
                         "m" ((*v).c23),
                         "m" ((*v).c33)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec(u);

   __asm__ __volatile__ ("vaddsd %%xmm15, %%xmm5, %%xmm5\n\t"
                         "vmovapd %%xmm5, %%xmm15 \n\t"
                         "vaddsubpd %0, %%xmm3, %%xmm3 \n\t"
                         "vaddsubpd %2, %%xmm4, %%xmm4 \n\t"
                         "vmovlpd %4, %%xmm5, %%xmm5 \n\t"
                         "vmovddup %5, %%xmm6 \n\t"
                         "vmulpd %6, %%xmm3, %%xmm3 \n\t"
                         "vsubpd %%xmm6, %%xmm5, %%xmm5 \n\t"
                         "vmulpd %6, %%xmm4, %%xmm4 \n\t"
                         "vmulpd %7, %%xmm5, %%xmm5"
                         :
                         :
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" ((*X).c2),
                         "m" ((*X).c1),
                         "m" (c0),
                         "m" (c1)
                         :
                         "xmm3", "xmm4", "xmm5", "xmm6",
                         "xmm15");

   __asm__ __volatile__ ("vmovapd %%xmm3, %0 \n\t"
                         "vmovapd %%xmm4, %2 \n\t"
                         "vmovapd %%xmm5, %4 \n\t"
                         "vmovlpd %%xmm15, %6"
                         :
                         "=m" ((*X).c5),
                         "=m" ((*X).c6),
                         "=m" ((*X).c7),
                         "=m" ((*X).c8),
                         "=m" ((*X).c1),
                         "=m" ((*X).c2),
                         "=m" (tr));

   _avx_zeroupper();

   return tr;
}


void prod2u3alg(su3_dble *u,su3_dble *v,u3_alg_dble *X)
{
   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2 \n\t"
                         "vinsertf128 $0x1, %3, %%ymm0, %%ymm0 \n\t"
                         "vinsertf128 $0x1, %4, %%ymm1, %%ymm1 \n\t"
                         "vinsertf128 $0x1, %5, %%ymm2, %%ymm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c21),
                         "m" ((*v).c31),
                         "m" ((*v).c12),
                         "m" ((*v).c22),
                         "m" ((*v).c32)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec_pair(u);

   __asm__ __volatile__ ("vmovupd %%xmm5, %0 \n\t"
                         "vextractf128 $0x1, %%ymm3, %%xmm6 \n\t"
                         "vextractf128 $0x1, %%ymm4, %%xmm7 \n\t"
                         "vextractf128 $0x1, %%ymm5, %2 \n\t"
                         "vaddsubpd %%xmm4, %%xmm6, %%xmm6 \n\t"
                         "vaddsubpd %%xmm3, %%xmm3, %%xmm3 \n\t"
                         "vaddsubpd %%xmm7, %%xmm7, %%xmm7"
                         :
                         "=m" ((*X).c6),
                         "=m" ((*X).c7),
                         "=m" ((*X).c8),
                         "=m" ((*X).c9)
                         :
                         :
                         "xmm3", "xmm6", "xmm7");

   __asm__ __volatile__ ("vmovupd %%xmm6, %0 \n\t"
                         "vmovhpd %%xmm3, %2 \n\t"
                         "vmovhpd %%xmm7, %3"
                         :
                         "=m" ((*X).c4),
                         "=m" ((*X).c5),
                         "=m" ((*X).c1),
                         "=m" ((*X).c2));

   __asm__ __volatile__ ("vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c13),
                         "m" ((*v).c23),
                         "m" ((*v).c33)
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec(u);

   __asm__ __volatile__ ("vmovupd %0, %%xmm0 \n\t"
                         "vmovupd %2, %%xmm1 \n\t"
                         "vaddsubpd %%xmm5, %%xmm5, %%xmm5 \n\t"
                         "vaddsubpd %%xmm0, %%xmm3, %%xmm3 \n\t"
                         "vaddsubpd %%xmm1, %%xmm4, %%xmm4"
                         :
                         :
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" ((*X).c9)
                         :
                         "xmm0", "xmm1", "xmm3",
                         "xmm4", "xmm5");

   __asm__ __volatile__ ("vmovhpd %%xmm5, %0 \n\t"
                         "vmovupd %%xmm3, %1 \n\t"
                         "vmovupd %%xmm4, %3"
                         :
                         "=m" ((*X).c3),
                         "=m" ((*X).c6),
                         "=m" ((*X).c7),
                         "=m" ((*X).c8),
                         "=m" ((*X).c9));

   _avx_zeroupper();
}


void rotate_su3alg(su3_dble *u,su3_alg_dble *X)
{
   __asm__ __volatile__ ("vxorpd %%xmm5, %%xmm5, %%xmm5 \n\t"
                         "vmovsd %0, %%xmm0 \n\t"
                         "vmovsd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2 \n\t"
                         "vmovapd %4, %%xmm3 \n\t"
                         "vmovapd %6, %%xmm4"
                         :
                         :
                         "m" ((*X).c1),
                         "m" ((*X).c2),
                         "m" ((*X).c3),
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8)
                         :
                         "xmm0", "xmm1", "xmm2",
                         "xmm3", "xmm4", "xmm5");

   __asm__ __volatile__ ("vaddsubpd %%xmm2, %%xmm5, %%xmm6 \n\t"
                         "vaddsubpd %%xmm3, %%xmm5, %%xmm7 \n\t"
                         "vaddsubpd %%xmm4, %%xmm5, %%xmm8 \n\t"
                         "vmovapd %%xmm2, %0 \n\t"
                         "vmovapd %%xmm3, %1 \n\t"
                         "vmovapd %%xmm4, %2 \n\t"
                         "vmovapd %%xmm6, %3 \n\t"
                         "vmovapd %%xmm7, %4 \n\t"
                         "vmovapd %%xmm8, %5"
                         :
                         "=m" (uX.c12),
                         "=m" (uX.c13),
                         "=m" (uX.c23),
                         "=m" (uX.c21),
                         "=m" (uX.c31),
                         "=m" (uX.c32)
                         :
                         :
                         "xmm6", "xmm7", "xmm8");

   __asm__ __volatile__ ("vaddsd %%xmm0, %%xmm1, %%xmm2 \n\t"
                         "vsubsd %%xmm0, %%xmm1, %%xmm3 \n\t"
                         "vsubsd %%xmm1, %%xmm0, %%xmm4 \n\t"
                         "vsubsd %%xmm0, %%xmm3, %%xmm3 \n\t"
                         "vsubsd %%xmm1, %%xmm4, %%xmm4 \n\t"
                         "vpermilpd $0x1, %%xmm2, %%xmm2 \n\t"
                         "vpermilpd $0x1, %%xmm3, %%xmm3 \n\t"
                         "vpermilpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "vmovapd %%xmm2, %0 \n\t"
                         "vmovapd %%xmm3, %1 \n\t"
                         "vmovapd %%xmm4, %2"
                         :
                         "=m" (uX.c11),
                         "=m" (uX.c22),
                         "=m" (uX.c33)
                         :
                         :
                         "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("vbroadcastf128 %6, %%ymm3 \n\t"
                         "vmovapd %0, %%xmm0 \n\t"
                         "vmovapd %1, %%xmm1 \n\t"
                         "vmovapd %2, %%xmm2 \n\t"
                         "vinsertf128 $0x1, %3, %%ymm0, %%ymm0 \n\t"
                         "vinsertf128 $0x1, %4, %%ymm1, %%ymm1 \n\t"
                         "vinsertf128 $0x1, %5, %%ymm2, %%ymm2 \n\t"
                         "vmulpd %%ymm3, %%ymm0, %%ymm0 \n\t"
                         "vmulpd %%ymm3, %%ymm1, %%ymm1 \n\t"
                         "vmulpd %%ymm3, %%ymm2, %%ymm2"
                         :
                         :
                         "m" ((*u).c11),
                         "m" ((*u).c12),
                         "m" ((*u).c13),
                         "m" ((*u).c21),
                         "m" ((*u).c22),
                         "m" ((*u).c23),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3");

   su3xsu3vec_pair(&uX);

   __asm__ __volatile__ ("vmovapd %%ymm3, %%ymm0 \n\t"
                         "vmovapd %%ymm4, %%ymm1 \n\t"
                         "vmovapd %%ymm5, %%ymm2"
                         :
                         :
                         :
                         "xmm0", "xmm1", "xmm2");

   su3xsu3vec_pair(u);

   __asm__ __volatile__ ("vxorpd %%ymm0, %%ymm0, %%ymm0 \n\t"
                         "vextractf128 $0x1, %%ymm4, %%xmm2 \n\t"
                         "vextractf128 $0x1, %%ymm3, %0 \n\t"
                         "vaddsubpd %%ymm5, %%ymm0, %%ymm0 \n\t"
                         "vpermilpd $0x1, %%xmm2, %%xmm2 \n\t"
                         "vpermilpd $0x1, %%xmm3, %%xmm3 \n\t"
                         "vmovapd %%xmm0, %2 \n\t"
                         "vextractf128 $0x1, %%ymm0, %4 \n\t"
                         "vsubsd %%xmm3, %%xmm2, %%xmm1 \n\t"
                         "vmulsd %8, %%xmm1, %%xmm1 \n\t"
                         "vsubsd %%xmm1, %%xmm3, %%xmm2 \n\t"
                         "vmovsd %%xmm1, %6 \n\t"
                         "vmovsd %%xmm2, %7"
                         :
                         "=m" ((*X).c3),
                         "=m" ((*X).c4),
                         "=m" ((*X).c5),
                         "=m" ((*X).c6),
                         "=m" ((*X).c7),
                         "=m" ((*X).c8),
                         "=m" ((*X).c1),
                         "=m" ((*X).c2)
                         :
                         "m" (c1.c1)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3");

   _avx_zeroupper();
}

#elif (defined x64)
#include "sse2.h"

static const sse_double c0={0.5,0.5},c1={-1.0/3.0,-1.0/3.0};
static su3_dble uX ALIGNED16;
static double tr ALIGNED8;


static void su3xsu3vec(su3_dble *u)
{
   _sse_su3_multiply_dble(*u);
}


static void su3dagxsu3vec(su3_dble *u)
{
   _sse_su3_inverse_multiply_dble(*u);
}


void su3xsu3(su3_dble *u,su3_dble *v,su3_dble *w)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c21),
                         "m" ((*v).c31)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c11),
                         "=m" ((*w).c21),
                         "=m" ((*w).c31));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c12),
                         "m" ((*v).c22),
                         "m" ((*v).c32)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c12),
                         "=m" ((*w).c22),
                         "=m" ((*w).c32));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c13),
                         "m" ((*v).c23),
                         "m" ((*v).c33)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c13),
                         "=m" ((*w).c23),
                         "=m" ((*w).c33));
}


void su3dagxsu3(su3_dble *u,su3_dble *v,su3_dble *w)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c21),
                         "m" ((*v).c31)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c11),
                         "=m" ((*w).c21),
                         "=m" ((*w).c31));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c12),
                         "m" ((*v).c22),
                         "m" ((*v).c32)
                        :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c12),
                         "=m" ((*w).c22),
                         "=m" ((*w).c32));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c13),
                         "m" ((*v).c23),
                         "m" ((*v).c33)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c13),
                         "=m" ((*w).c23),
                         "=m" ((*w).c33));
}


void su3xsu3dag(su3_dble *u,su3_dble *v,su3_dble *w)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm0 \n\t"
                         "mulpd %3, %%xmm1 \n\t"
                         "mulpd %3, %%xmm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c12),
                         "m" ((*v).c13),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c11),
                         "=m" ((*w).c21),
                         "=m" ((*w).c31));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm0 \n\t"
                         "mulpd %3, %%xmm1 \n\t"
                         "mulpd %3, %%xmm2"
                         :
                         :
                         "m" ((*v).c21),
                         "m" ((*v).c22),
                         "m" ((*v).c23),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c12),
                         "=m" ((*w).c22),
                         "=m" ((*w).c32));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm0 \n\t"
                         "mulpd %3, %%xmm1 \n\t"
                         "mulpd %3, %%xmm2"
                         :
                         :
                         "m" ((*v).c31),
                         "m" ((*v).c32),
                         "m" ((*v).c33),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c13),
                         "=m" ((*w).c23),
                         "=m" ((*w).c33));
}


void su3dagxsu3dag(su3_dble *u,su3_dble *v,su3_dble *w)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm0 \n\t"
                         "mulpd %3, %%xmm1 \n\t"
                         "mulpd %3, %%xmm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c12),
                         "m" ((*v).c13),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c11),
                         "=m" ((*w).c21),
                         "=m" ((*w).c31));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm0 \n\t"
                         "mulpd %3, %%xmm1 \n\t"
                         "mulpd %3, %%xmm2"
                         :
                         :
                         "m" ((*v).c21),
                         "m" ((*v).c22),
                         "m" ((*v).c23),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c12),
                         "=m" ((*w).c22),
                         "=m" ((*w).c32));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm0 \n\t"
                         "mulpd %3, %%xmm1 \n\t"
                         "mulpd %3, %%xmm2"
                         :
                         :
                         "m" ((*v).c31),
                         "m" ((*v).c32),
                         "m" ((*v).c33),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c13),
                         "=m" ((*w).c23),
                         "=m" ((*w).c33));
}


void su3xu3alg(su3_dble *u,u3_alg_dble *X,su3_dble *v)
{
   __asm__ __volatile__ ("xorpd %%xmm0, %%xmm0\n\t"
                         "movhpd %0, %%xmm0 \n\t"
                         "movupd %1, %%xmm1 \n\t"
                         "movupd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm1 \n\t"
                         "mulpd %3, %%xmm2"
                         :
                         :
                         "m" ((*X).c1),
                         "m" ((*X).c4),
                         "m" ((*X).c6),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c21),
                         "=m" ((*v).c31));

   __asm__ __volatile__ ("movupd %0, %%xmm0 \n\t"
                         "xorpd %%xmm1, %%xmm1\n\t"
                         "movhpd %1, %%xmm1 \n\t"
                         "movupd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm2"
                         :
                         :
                         "m" ((*X).c4),
                         "m" ((*X).c2),
                         "m" ((*X).c8),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c12),
                         "=m" ((*v).c22),
                         "=m" ((*v).c32));


   __asm__ __volatile__ ("movupd %0, %%xmm0 \n\t"
                         "movupd %1, %%xmm1 \n\t"
                         "xorpd %%xmm2, %%xmm2\n\t"
                         "movhpd %2, %%xmm2"
                         :
                         :
                         "m" ((*X).c6),
                         "m" ((*X).c8),
                         "m" ((*X).c3)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c13),
                         "=m" ((*v).c23),
                         "=m" ((*v).c33));
}


void su3dagxu3alg(su3_dble *u,u3_alg_dble *X,su3_dble *v)
{
   __asm__ __volatile__ ("xorpd %%xmm0, %%xmm0\n\t"
                         "movhpd %0, %%xmm0 \n\t"
                         "movupd %1, %%xmm1 \n\t"
                         "movupd %3, %%xmm2 \n\t"
                         "mulpd %5, %%xmm1 \n\t"
                         "mulpd %5, %%xmm2"
                         :
                         :
                         "m" ((*X).c1),
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c21),
                         "=m" ((*v).c31));

   __asm__ __volatile__ ("movupd %0, %%xmm0 \n\t"
                         "xorpd %%xmm1, %%xmm1\n\t"
                         "movhpd %2, %%xmm1 \n\t"
                         "movupd %3, %%xmm2 \n\t"
                         "mulpd %5, %%xmm2"
                         :
                         :
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c2),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c12),
                         "=m" ((*v).c22),
                         "=m" ((*v).c32));


   __asm__ __volatile__ ("movupd %0, %%xmm0 \n\t"
                         "movupd %2, %%xmm1 \n\t"
                         "xorpd %%xmm2, %%xmm2\n\t"
                         "movhpd %4, %%xmm2"
                         :
                         :
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" ((*X).c3)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c13),
                         "=m" ((*v).c23),
                         "=m" ((*v).c33));
}


void u3algxsu3(u3_alg_dble *X,su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("xorpd %%xmm0, %%xmm0\n\t"
                         "movhpd %0, %%xmm0 \n\t"
                         "movupd %1, %%xmm1 \n\t"
                         "movupd %3, %%xmm2 \n\t"
                         "mulpd %5, %%xmm1 \n\t"
                         "mulpd %5, %%xmm2"
                         :
                         :
                         "m" ((*X).c1),
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("mulpd %3, %%xmm3\n\t"
                         "mulpd %3, %%xmm4\n\t"
                         "mulpd %3, %%xmm5\n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c12),
                         "=m" ((*v).c13)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm4", "xmm5");

   __asm__ __volatile__ ("movupd %0, %%xmm0 \n\t"
                         "xorpd %%xmm1, %%xmm1\n\t"
                         "movhpd %2, %%xmm1 \n\t"
                         "movupd %3, %%xmm2 \n\t"
                         "mulpd %5, %%xmm2"
                         :
                         :
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c2),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("mulpd %3, %%xmm3\n\t"
                         "mulpd %3, %%xmm4\n\t"
                         "mulpd %3, %%xmm5\n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c21),
                         "=m" ((*v).c22),
                         "=m" ((*v).c23)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm4", "xmm5");

   __asm__ __volatile__ ("movupd %0, %%xmm0 \n\t"
                         "movupd %2, %%xmm1 \n\t"
                         "xorpd %%xmm2, %%xmm2\n\t"
                         "movhpd %4, %%xmm2"
                         :
                         :
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" ((*X).c3)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3dagxsu3vec(u);
   __asm__ __volatile__ ("mulpd %3, %%xmm3\n\t"
                         "mulpd %3, %%xmm4\n\t"
                         "mulpd %3, %%xmm5\n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c31),
                         "=m" ((*v).c32),
                         "=m" ((*v).c33)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm4", "xmm5");
}


void u3algxsu3dag(u3_alg_dble *X,su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("xorpd %%xmm0, %%xmm0\n\t"
                         "movhpd %0, %%xmm0 \n\t"
                         "movupd %1, %%xmm1 \n\t"
                         "movupd %3, %%xmm2 \n\t"
                         "mulpd %5, %%xmm1 \n\t"
                         "mulpd %5, %%xmm2"
                         :
                         :
                         "m" ((*X).c1),
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("mulpd %3, %%xmm3\n\t"
                         "mulpd %3, %%xmm4\n\t"
                         "mulpd %3, %%xmm5\n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c12),
                         "=m" ((*v).c13)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm4", "xmm5");

   __asm__ __volatile__ ("movupd %0, %%xmm0 \n\t"
                         "xorpd %%xmm1, %%xmm1\n\t"
                         "movhpd %2, %%xmm1 \n\t"
                         "movupd %3, %%xmm2 \n\t"
                         "mulpd %5, %%xmm2"
                         :
                         :
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c2),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("mulpd %3, %%xmm3\n\t"
                         "mulpd %3, %%xmm4\n\t"
                         "mulpd %3, %%xmm5\n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c21),
                         "=m" ((*v).c22),
                         "=m" ((*v).c23)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm4", "xmm5");

   __asm__ __volatile__ ("movupd %0, %%xmm0 \n\t"
                         "movupd %2, %%xmm1 \n\t"
                         "xorpd %%xmm2, %%xmm2\n\t"
                         "movhpd %4, %%xmm2"
                         :
                         :
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" ((*X).c9),
                         "m" ((*X).c3)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("mulpd %3, %%xmm3\n\t"
                         "mulpd %3, %%xmm4\n\t"
                         "mulpd %3, %%xmm5\n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c31),
                         "=m" ((*v).c32),
                         "=m" ((*v).c33)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm4", "xmm5");
}


double prod2su3alg(su3_dble *u,su3_dble *v,su3_alg_dble *X)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c21),
                         "m" ((*v).c31)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movsd %%xmm3, %%xmm15 \n\t"
                         "mulpd %6, %%xmm4 \n\t"
                         "mulpd %6, %%xmm5 \n\t"
                         "movhpd %%xmm3, %0 \n\t"
                         "movhpd %%xmm3, %1 \n\t"
                         "movapd %%xmm4, %2 \n\t"
                         "movapd %%xmm5, %4"
                         :
                         "=m" ((*X).c1),
                         "=m" ((*X).c2),
                         "=m" ((*X).c3),
                         "=m" ((*X).c4),
                         "=m" ((*X).c5),
                         "=m" ((*X).c6)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm4", "xmm5", "xmm15");

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c12),
                         "m" ((*v).c22),
                         "m" ((*v).c32)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("addsd %%xmm15, %%xmm4\n\t"
                         "movsd %%xmm4, %%xmm15 \n\t"
                         "addpd %0, %%xmm3 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "mulpd %2, %%xmm3 \n\t"
                         "subsd %3, %%xmm4 \n\t"
                         "mulpd %4, %%xmm5 \n\t"
                         "mulsd %5, %%xmm4"
                         :
                         :
                         "m" ((*X).c3),
                         "m" ((*X).c4),
                         "m" (c0),
                         "m" ((*X).c1),
                         "m" (_sse_sgn1_dble),
                         "m" (c1)
                         :
                         "xmm3", "xmm4", "xmm5", "xmm15");

   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movlpd %%xmm4, %2 \n\t"
                         "movapd %%xmm5, %3"
                         :
                         "=m" ((*X).c3),
                         "=m" ((*X).c4),
                         "=m" ((*X).c1),
                         "=m" ((*X).c7),
                         "=m" ((*X).c8));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c13),
                         "m" ((*v).c23),
                         "m" ((*v).c33)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("addsd %%xmm15, %%xmm5 \n\t"
                         "movsd %%xmm5, %%xmm15 \n\t"
                         "addpd %0, %%xmm3 \n\t"
                         "addpd %2, %%xmm4 \n\t"
                         "shufpd $0x1, %%xmm5, %%xmm5 \n\t"
                         "mulpd %4, %%xmm3 \n\t"
                         "subsd %5, %%xmm5 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "mulsd %6, %%xmm5"
                         :
                         :
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" (c0),
                         "m" ((*X).c2),
                         "m" (c1)
                         :
                         "xmm3", "xmm4", "xmm5", "xmm15");

   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %2 \n\t"
                         "movlpd %%xmm5, %4 \n\t"
                         "movlpd %%xmm15, %5"
                         :
                         "=m" ((*X).c5),
                         "=m" ((*X).c6),
                         "=m" ((*X).c7),
                         "=m" ((*X).c8),
                         "=m" ((*X).c2),
                         "=m" (tr));

   return tr;
}


void prod2u3alg(su3_dble *u,su3_dble *v,u3_alg_dble *X)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c21),
                         "m" ((*v).c31)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                         "mulpd %0, %%xmm4 \n\t"
                         "mulpd %0, %%xmm5 \n\t"
                         "addsd %%xmm3, %%xmm3"
                         :
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm4", "xmm5");

   __asm__ __volatile__ ("movupd %%xmm4, %0 \n\t"
                         "movupd %%xmm5, %2 \n\t"
                         "movlpd %%xmm3, %4"
                         :
                         "=m" ((*X).c4),
                         "=m" ((*X).c5),
                         "=m" ((*X).c6),
                         "=m" ((*X).c7),
                         "=m" ((*X).c1));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c12),
                         "m" ((*v).c22),
                         "m" ((*v).c32)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movupd %0, %%xmm6 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "mulpd %2, %%xmm5 \n\t"
                         "addpd %%xmm6, %%xmm3 \n\t"
                         "addsd %%xmm4, %%xmm4"
                         :
                         :
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm4", "xmm5",
                         "xmm6");

   __asm__ __volatile__ ("movupd %%xmm5, %0 \n\t"
                         "movupd %%xmm3, %2 \n\t"
                         "movlpd %%xmm4, %4"
                         :
                         "=m" ((*X).c8),
                         "=m" ((*X).c9),
                         "=m" ((*X).c4),
                         "=m" ((*X).c5),
                         "=m" ((*X).c2));

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2"
                         :
                         :
                         "m" ((*v).c13),
                         "m" ((*v).c23),
                         "m" ((*v).c33)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movupd %0, %%xmm6 \n\t"
                         "movupd %2, %%xmm7 \n\t"
                         "shufpd $0x1, %%xmm5, %%xmm5 \n\t"
                         "addpd %%xmm6, %%xmm3 \n\t"
                         "addpd %%xmm7, %%xmm4 \n\t"
                         "addsd %%xmm5, %%xmm5"
                         :
                         :
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8),
                         "m" ((*X).c9)
                         :
                         "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7");

   __asm__ __volatile__ ("movupd %%xmm3, %0 \n\t"
                         "movupd %%xmm4, %2 \n\t"
                         "movlpd %%xmm5, %4"
                         :
                         "=m" ((*X).c6),
                         "=m" ((*X).c7),
                         "=m" ((*X).c8),
                         "=m" ((*X).c9),
                         "=m" ((*X).c3));
}


void rotate_su3alg(su3_dble *u,su3_alg_dble *X)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %2, %%xmm3 \n\t"
                         "movapd %4, %%xmm5 \n\t"
                         "movapd %6, %%xmm7 \n\t"
                         "movapd %%xmm0, %%xmm1 \n\t"
                         "movsd %%xmm0, %%xmm2 \n\t"
                         :
                         :
                         "m" ((*X).c1),
                         "m" ((*X).c2),
                         "m" ((*X).c3),
                         "m" ((*X).c4),
                         "m" ((*X).c5),
                         "m" ((*X).c6),
                         "m" ((*X).c7),
                         "m" ((*X).c8)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm5", "xmm7");

   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm5, %1 \n\t"
                         "movapd %%xmm7, %2 \n\t"
                         "mulpd %6, %%xmm3 \n\t"
                         "mulpd %6, %%xmm5 \n\t"
                         "mulpd %6, %%xmm7 \n\t"
                         "movapd %%xmm3, %3 \n\t"
                         "movapd %%xmm5, %4 \n\t"
                         "movapd %%xmm7, %5"
                         :
                         "=m" (uX.c12),
                         "=m" (uX.c13),
                         "=m" (uX.c23),
                         "=m" (uX.c21),
                         "=m" (uX.c31),
                         "=m" (uX.c32)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm5", "xmm7");

   __asm__ __volatile__ ("shufpd $0x1, %%xmm0, %%xmm0 \n\t"
                         "addpd %%xmm1, %%xmm1 \n\t"
                         "addsd %%xmm0, %%xmm2 \n\t"
                         "subpd %%xmm1, %%xmm0 \n\t"
                         "xorpd %%xmm3, %%xmm3 \n\t"
                         "movlpd %%xmm2, %0 \n\t"
                         "movlpd %%xmm0, %1 \n\t"
                         "movhpd %%xmm0, %2 \n\t"
                         "movlpd %%xmm3, %3 \n\t"
                         "movlpd %%xmm3, %4 \n\t"
                         "movlpd %%xmm3, %5"
                         :
                         "=m" (uX.c11.im),
                         "=m" (uX.c22.im),
                         "=m" (uX.c33.im),
                         "=m" (uX.c11.re),
                         "=m" (uX.c22.re),
                         "=m" (uX.c33.re)
                         :
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3");

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm0 \n\t"
                         "mulpd %3, %%xmm1 \n\t"
                         "mulpd %3, %%xmm2"
                         :
                         :
                         "m" ((*u).c11),
                         "m" ((*u).c12),
                         "m" ((*u).c13),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(&uX);
   __asm__ __volatile__ ("movapd %%xmm3, %%xmm0 \n\t"
                         "movapd %%xmm4, %%xmm1 \n\t"
                         "movapd %%xmm5, %%xmm2"
                         :
                         :
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movhpd %%xmm3, %0 \n\t"
                         "mulpd %4, %%xmm5 \n\t"
                         "addpd %%xmm3, %%xmm3 \n\t"
                         "movapd %%xmm5, %1 \n\t"
                         "movhpd %%xmm3, %3"
                         :
                         "=m" ((*X).c1),
                         "=m" ((*X).c5),
                         "=m" ((*X).c6),
                         "=m" ((*X).c2)
                         :
                         "m" (_sse_sgn1_dble)
                         :
                         "xmm3", "xmm5");

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm0 \n\t"
                         "mulpd %3, %%xmm1 \n\t"
                         "mulpd %3, %%xmm2"
                         :
                         :
                         "m" ((*u).c21),
                         "m" ((*u).c22),
                         "m" ((*u).c23),
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(&uX);
   __asm__ __volatile__ ("movapd %%xmm3, %%xmm0 \n\t"
                         "movapd %%xmm4, %%xmm1 \n\t"
                         "movapd %%xmm5, %%xmm2"
                         :
                         :
                         :
                         "xmm0", "xmm1", "xmm2");
   su3xsu3vec(u);
   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "unpckhpd %%xmm4, %%xmm4 \n\t"
                         "mulpd %2, %%xmm5 \n\t"
                         "mulpd %3, %%xmm4 \n\t"
                         "subpd %4, %%xmm4"
                         :
                         "=m" ((*X).c3),
                         "=m" ((*X).c4)
                         :
                         "m" (_sse_sgn1_dble),
                         "m" (_sse_sgn2_dble),
                         "m" ((*X).c1),
                         "m" ((*X).c2)
                         :
                         "xmm4", "xmm5");

   __asm__ __volatile__ ("movapd %%xmm5, %0 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "movapd %%xmm4, %2"
                         :
                         "=m" ((*X).c7),
                         "=m" ((*X).c8),
                         "=m" ((*X).c1),
                         "=m" ((*X).c2)
                         :
                         "m" (c1)
                         :
                         "xmm4");
}

#else

static su3_vector_dble psi,chi;
static su3_dble uX;


static void su3xsu3vec(su3_dble *u)
{
   _su3_multiply(chi,*u,psi);
}


static void su3dagxsu3vec(su3_dble *u)
{
   _su3_inverse_multiply(chi,*u,psi);
}


void su3xsu3(su3_dble *u,su3_dble *v,su3_dble *w)
{
   psi.c1=(*v).c11;
   psi.c2=(*v).c21;
   psi.c3=(*v).c31;
   su3xsu3vec(u);
   (*w).c11=chi.c1;
   (*w).c21=chi.c2;
   (*w).c31=chi.c3;

   psi.c1=(*v).c12;
   psi.c2=(*v).c22;
   psi.c3=(*v).c32;
   su3xsu3vec(u);
   (*w).c12=chi.c1;
   (*w).c22=chi.c2;
   (*w).c32=chi.c3;

   psi.c1=(*v).c13;
   psi.c2=(*v).c23;
   psi.c3=(*v).c33;
   su3xsu3vec(u);
   (*w).c13=chi.c1;
   (*w).c23=chi.c2;
   (*w).c33=chi.c3;
}


void su3dagxsu3(su3_dble *u,su3_dble *v,su3_dble *w)
{
   psi.c1=(*v).c11;
   psi.c2=(*v).c21;
   psi.c3=(*v).c31;
   su3dagxsu3vec(u);
   (*w).c11=chi.c1;
   (*w).c21=chi.c2;
   (*w).c31=chi.c3;

   psi.c1=(*v).c12;
   psi.c2=(*v).c22;
   psi.c3=(*v).c32;
   su3dagxsu3vec(u);
   (*w).c12=chi.c1;
   (*w).c22=chi.c2;
   (*w).c32=chi.c3;

   psi.c1=(*v).c13;
   psi.c2=(*v).c23;
   psi.c3=(*v).c33;
   su3dagxsu3vec(u);
   (*w).c13=chi.c1;
   (*w).c23=chi.c2;
   (*w).c33=chi.c3;
}


void su3xsu3dag(su3_dble *u,su3_dble *v,su3_dble *w)
{
   psi.c1.re= (*v).c11.re;
   psi.c1.im=-(*v).c11.im;
   psi.c2.re= (*v).c12.re;
   psi.c2.im=-(*v).c12.im;
   psi.c3.re= (*v).c13.re;
   psi.c3.im=-(*v).c13.im;
   su3xsu3vec(u);
   (*w).c11=chi.c1;
   (*w).c21=chi.c2;
   (*w).c31=chi.c3;

   psi.c1.re= (*v).c21.re;
   psi.c1.im=-(*v).c21.im;
   psi.c2.re= (*v).c22.re;
   psi.c2.im=-(*v).c22.im;
   psi.c3.re= (*v).c23.re;
   psi.c3.im=-(*v).c23.im;
   su3xsu3vec(u);
   (*w).c12=chi.c1;
   (*w).c22=chi.c2;
   (*w).c32=chi.c3;

   psi.c1.re= (*v).c31.re;
   psi.c1.im=-(*v).c31.im;
   psi.c2.re= (*v).c32.re;
   psi.c2.im=-(*v).c32.im;
   psi.c3.re= (*v).c33.re;
   psi.c3.im=-(*v).c33.im;
   su3xsu3vec(u);
   (*w).c13=chi.c1;
   (*w).c23=chi.c2;
   (*w).c33=chi.c3;
}


void su3dagxsu3dag(su3_dble *u,su3_dble *v,su3_dble *w)
{
   psi.c1.re= (*v).c11.re;
   psi.c1.im=-(*v).c11.im;
   psi.c2.re= (*v).c12.re;
   psi.c2.im=-(*v).c12.im;
   psi.c3.re= (*v).c13.re;
   psi.c3.im=-(*v).c13.im;
   su3dagxsu3vec(u);
   (*w).c11=chi.c1;
   (*w).c21=chi.c2;
   (*w).c31=chi.c3;

   psi.c1.re= (*v).c21.re;
   psi.c1.im=-(*v).c21.im;
   psi.c2.re= (*v).c22.re;
   psi.c2.im=-(*v).c22.im;
   psi.c3.re= (*v).c23.re;
   psi.c3.im=-(*v).c23.im;
   su3dagxsu3vec(u);
   (*w).c12=chi.c1;
   (*w).c22=chi.c2;
   (*w).c32=chi.c3;

   psi.c1.re= (*v).c31.re;
   psi.c1.im=-(*v).c31.im;
   psi.c2.re= (*v).c32.re;
   psi.c2.im=-(*v).c32.im;
   psi.c3.re= (*v).c33.re;
   psi.c3.im=-(*v).c33.im;
   su3dagxsu3vec(u);
   (*w).c13=chi.c1;
   (*w).c23=chi.c2;
   (*w).c33=chi.c3;
}


void su3xu3alg(su3_dble *u,u3_alg_dble *X,su3_dble *v)
{
   psi.c1.re=0.0;
   psi.c1.im= (*X).c1;
   psi.c2.re=-(*X).c4;
   psi.c2.im= (*X).c5;
   psi.c3.re=-(*X).c6;
   psi.c3.im= (*X).c7;
   su3xsu3vec(u);
   (*v).c11=chi.c1;
   (*v).c21=chi.c2;
   (*v).c31=chi.c3;

   psi.c1.re= (*X).c4;
   psi.c1.im= (*X).c5;
   psi.c2.re=0.0;
   psi.c2.im= (*X).c2;
   psi.c3.re=-(*X).c8;
   psi.c3.im= (*X).c9;
   su3xsu3vec(u);
   (*v).c12=chi.c1;
   (*v).c22=chi.c2;
   (*v).c32=chi.c3;

   psi.c1.re= (*X).c6;
   psi.c1.im= (*X).c7;
   psi.c2.re= (*X).c8;
   psi.c2.im= (*X).c9;
   psi.c3.re=0.0;
   psi.c3.im= (*X).c3;
   su3xsu3vec(u);
   (*v).c13=chi.c1;
   (*v).c23=chi.c2;
   (*v).c33=chi.c3;
}


void su3dagxu3alg(su3_dble *u,u3_alg_dble *X,su3_dble *v)
{
   psi.c1.re=0.0;
   psi.c1.im= (*X).c1;
   psi.c2.re=-(*X).c4;
   psi.c2.im= (*X).c5;
   psi.c3.re=-(*X).c6;
   psi.c3.im= (*X).c7;
   su3dagxsu3vec(u);
   (*v).c11=chi.c1;
   (*v).c21=chi.c2;
   (*v).c31=chi.c3;

   psi.c1.re= (*X).c4;
   psi.c1.im= (*X).c5;
   psi.c2.re=0.0;
   psi.c2.im= (*X).c2;
   psi.c3.re=-(*X).c8;
   psi.c3.im= (*X).c9;
   su3dagxsu3vec(u);
   (*v).c12=chi.c1;
   (*v).c22=chi.c2;
   (*v).c32=chi.c3;

   psi.c1.re= (*X).c6;
   psi.c1.im= (*X).c7;
   psi.c2.re= (*X).c8;
   psi.c2.im= (*X).c9;
   psi.c3.re=0.0;
   psi.c3.im= (*X).c3;
   su3dagxsu3vec(u);
   (*v).c13=chi.c1;
   (*v).c23=chi.c2;
   (*v).c33=chi.c3;
}


void u3algxsu3(u3_alg_dble *X,su3_dble *u,su3_dble *v)
{
   psi.c1.re=0.0;
   psi.c1.im= (*X).c1;
   psi.c2.re=-(*X).c4;
   psi.c2.im= (*X).c5;
   psi.c3.re=-(*X).c6;
   psi.c3.im= (*X).c7;
   su3dagxsu3vec(u);
   (*v).c11.re=-chi.c1.re;
   (*v).c11.im= chi.c1.im;
   (*v).c12.re=-chi.c2.re;
   (*v).c12.im= chi.c2.im;
   (*v).c13.re=-chi.c3.re;
   (*v).c13.im= chi.c3.im;

   psi.c1.re= (*X).c4;
   psi.c1.im= (*X).c5;
   psi.c2.re=0.0;
   psi.c2.im= (*X).c2;
   psi.c3.re=-(*X).c8;
   psi.c3.im= (*X).c9;
   su3dagxsu3vec(u);
   (*v).c21.re=-chi.c1.re;
   (*v).c21.im= chi.c1.im;
   (*v).c22.re=-chi.c2.re;
   (*v).c22.im= chi.c2.im;
   (*v).c23.re=-chi.c3.re;
   (*v).c23.im= chi.c3.im;

   psi.c1.re= (*X).c6;
   psi.c1.im= (*X).c7;
   psi.c2.re= (*X).c8;
   psi.c2.im= (*X).c9;
   psi.c3.re=0.0;
   psi.c3.im= (*X).c3;
   su3dagxsu3vec(u);
   (*v).c31.re=-chi.c1.re;
   (*v).c31.im= chi.c1.im;
   (*v).c32.re=-chi.c2.re;
   (*v).c32.im= chi.c2.im;
   (*v).c33.re=-chi.c3.re;
   (*v).c33.im= chi.c3.im;
}


void u3algxsu3dag(u3_alg_dble *X,su3_dble *u,su3_dble *v)
{
   psi.c1.re=0.0;
   psi.c1.im= (*X).c1;
   psi.c2.re=-(*X).c4;
   psi.c2.im= (*X).c5;
   psi.c3.re=-(*X).c6;
   psi.c3.im= (*X).c7;
   su3xsu3vec(u);
   (*v).c11.re=-chi.c1.re;
   (*v).c11.im= chi.c1.im;
   (*v).c12.re=-chi.c2.re;
   (*v).c12.im= chi.c2.im;
   (*v).c13.re=-chi.c3.re;
   (*v).c13.im= chi.c3.im;

   psi.c1.re= (*X).c4;
   psi.c1.im= (*X).c5;
   psi.c2.re=0.0;
   psi.c2.im= (*X).c2;
   psi.c3.re=-(*X).c8;
   psi.c3.im= (*X).c9;
   su3xsu3vec(u);
   (*v).c21.re=-chi.c1.re;
   (*v).c21.im= chi.c1.im;
   (*v).c22.re=-chi.c2.re;
   (*v).c22.im= chi.c2.im;
   (*v).c23.re=-chi.c3.re;
   (*v).c23.im= chi.c3.im;

   psi.c1.re= (*X).c6;
   psi.c1.im= (*X).c7;
   psi.c2.re= (*X).c8;
   psi.c2.im= (*X).c9;
   psi.c3.re=0.0;
   psi.c3.im= (*X).c3;
   su3xsu3vec(u);
   (*v).c31.re=-chi.c1.re;
   (*v).c31.im= chi.c1.im;
   (*v).c32.re=-chi.c2.re;
   (*v).c32.im= chi.c2.im;
   (*v).c33.re=-chi.c3.re;
   (*v).c33.im= chi.c3.im;
}


double prod2su3alg(su3_dble *u,su3_dble *v,su3_alg_dble *X)
{
   double tr;

   psi.c1=(*v).c11;
   psi.c2=(*v).c21;
   psi.c3=(*v).c31;
   su3xsu3vec(u);
   tr=chi.c1.re;
   (*X).c1 = chi.c1.im;
   (*X).c2 = chi.c1.im;
   (*X).c3 =-chi.c2.re;
   (*X).c4 = chi.c2.im;
   (*X).c5 =-chi.c3.re;
   (*X).c6 = chi.c3.im;

   psi.c1=(*v).c12;
   psi.c2=(*v).c22;
   psi.c3=(*v).c32;
   su3xsu3vec(u);
   tr+=chi.c2.re;
   (*X).c3+= chi.c1.re;
   (*X).c4+= chi.c1.im;
   (*X).c1-= chi.c2.im;
   (*X).c7 =-chi.c3.re;
   (*X).c8 = chi.c3.im;

   psi.c1=(*v).c13;
   psi.c2=(*v).c23;
   psi.c3=(*v).c33;
   su3xsu3vec(u);
   tr+=chi.c3.re;
   (*X).c5+= chi.c1.re;
   (*X).c6+= chi.c1.im;
   (*X).c7+= chi.c2.re;
   (*X).c8+= chi.c2.im;
   (*X).c2-= chi.c3.im;

   (*X).c1*=(1.0/3.0);
   (*X).c2*=(1.0/3.0);
   (*X).c3*=0.5;
   (*X).c4*=0.5;
   (*X).c5*=0.5;
   (*X).c6*=0.5;
   (*X).c7*=0.5;
   (*X).c8*=0.5;

   return tr;
}


void prod2u3alg(su3_dble *u,su3_dble *v,u3_alg_dble *X)
{
   psi.c1=(*v).c11;
   psi.c2=(*v).c21;
   psi.c3=(*v).c31;
   su3xsu3vec(u);
   (*X).c1=chi.c1.im+chi.c1.im;
   (*X).c4=-chi.c2.re;
   (*X).c5=chi.c2.im;
   (*X).c6=-chi.c3.re;
   (*X).c7=chi.c3.im;

   psi.c1=(*v).c12;
   psi.c2=(*v).c22;
   psi.c3=(*v).c32;
   su3xsu3vec(u);
   (*X).c4+=chi.c1.re;
   (*X).c5+=chi.c1.im;
   (*X).c2=chi.c2.im+chi.c2.im;
   (*X).c8=-chi.c3.re;
   (*X).c9=chi.c3.im;

   psi.c1=(*v).c13;
   psi.c2=(*v).c23;
   psi.c3=(*v).c33;
   su3xsu3vec(u);
   (*X).c6+=chi.c1.re;
   (*X).c7+=chi.c1.im;
   (*X).c8+=chi.c2.re;
   (*X).c9+=chi.c2.im;
   (*X).c3=chi.c3.im+chi.c3.im;
}


void rotate_su3alg(su3_dble *u,su3_alg_dble *X)
{
   uX.c11.re=0.0;
   uX.c11.im=(*X).c1+(*X).c2;
   uX.c22.re=0.0;
   uX.c22.im=(*X).c2-(*X).c1-(*X).c1;
   uX.c33.re=0.0;
   uX.c33.im=(*X).c1-(*X).c2-(*X).c2;

   uX.c12.re= (*X).c3;
   uX.c12.im= (*X).c4;
   uX.c21.re=-(*X).c3;
   uX.c21.im= (*X).c4;

   uX.c13.re= (*X).c5;
   uX.c13.im= (*X).c6;
   uX.c31.re=-(*X).c5;
   uX.c31.im= (*X).c6;

   uX.c23.re= (*X).c7;
   uX.c23.im= (*X).c8;
   uX.c32.re=-(*X).c7;
   uX.c32.im= (*X).c8;

   psi.c1.re= (*u).c11.re;
   psi.c1.im=-(*u).c11.im;
   psi.c2.re= (*u).c12.re;
   psi.c2.im=-(*u).c12.im;
   psi.c3.re= (*u).c13.re;
   psi.c3.im=-(*u).c13.im;
   su3xsu3vec(&uX);
   psi=chi;
   su3xsu3vec(u);
   (*X).c1= chi.c1.im;
   (*X).c2= chi.c1.im+chi.c1.im;
   (*X).c5=-chi.c3.re;
   (*X).c6= chi.c3.im;

   psi.c1.re= (*u).c21.re;
   psi.c1.im=-(*u).c21.im;
   psi.c2.re= (*u).c22.re;
   psi.c2.im=-(*u).c22.im;
   psi.c3.re= (*u).c23.re;
   psi.c3.im=-(*u).c23.im;
   su3xsu3vec(&uX);
   psi=chi;
   su3xsu3vec(u);
   (*X).c3= chi.c1.re;
   (*X).c4= chi.c1.im;
   (*X).c1-=chi.c2.im;
   (*X).c2+=chi.c2.im;
   (*X).c7=-chi.c3.re;
   (*X).c8= chi.c3.im;

   (*X).c1*=(1.0/3.0);
   (*X).c2*=(1.0/3.0);
}

#endif
