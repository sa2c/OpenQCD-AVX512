/*******************************************************************************
*
* File pauli.c
*
* Copyright (C) 2005, 2009, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Basic functions for single-precision Hermitian 6x6 matrices.
*
* The externally accessible functions are
*
*   void mul_pauli(float mu,pauli *m,weyl *s,weyl *r)
*     Multiplies the Weyl spinor s by the matrix m+i*mu and assigns
*     the result to the Weyl spinor r. The source spinor is overwritten
*     if r=s and otherwise left unchanged.
*
*   void mul_pauli2(float mu,pauli *m,spinor *s,spinor *r)
*     Multiplies the spinor s by the matrix m+i*mu*gamma_5 and assigns
*     the result to the spinor r. The source spinor is overwritten
*     if r=s and otherwise left unchanged.
*
*   void assign_pauli(int vol,pauli_dble *md,pauli *m)
*     Assigns the field md[vol] of double-precision matrices to the field
*     m[vol] of single-precision matrices.
*
*   void apply_sw(int vol,float mu,pauli *m,spinor *s,spinor *r)
*     Applies the matrix field m[2*vol]+i*mu*gamma_5 to the spinor field
*     s[vol] and assigns the result to the field r[vol]. The source field
*     is overwritten if r=s and otherwise left unchanged (the arrays may
*     not overlap in this case).
*
* Notes:
*
* The storage format for Hermitian 6x6 matrices is described in the notes
* "Implementation of the lattice Dirac operator" (file doc/dirac.pdf).
*
* The programs perform no communications and can be called locally. If SSE
* or AVX instructions are used, the Pauli matrices, Weyl and Dirac spinors
* must be aligned to a 16 byte boundary.
*
*******************************************************************************/

#define PAULI_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "sw_term.h"

typedef union
{
   spinor s;
   weyl w[2];
} spin_t;

#if (defined x64)
#include "sse2.h"

void mul_pauli(float mu,pauli *m,weyl *s,weyl *r)
{
   m+=4;
   _prefetch_pauli(m);
   m-=4;

   __asm__ __volatile__ ("movss %0, %%xmm14 \n\t"
                         "movss %1, %%xmm2 \n\t"
                         "movss %2, %%xmm3 \n\t"
                         "movsd %3, %%xmm4 \n\t"
                         "shufps $0xb1, %%xmm14, %%xmm14"
                         :
                         :
                         "m" (mu),
                         "m" ((*m).u[0]),
                         "m" ((*m).u[1]),
                         "m" ((*m).u[8]),
                         "m" ((*m).u[9])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm14");

   __asm__ __volatile__ ("movhps %0, %%xmm2 \n\t"
                         "movhps %0, %%xmm3 \n\t"
                         "movhps %2, %%xmm4 \n\t"
                         "movsldup %4, %%xmm0 \n\t"
                         "movshdup %4, %%xmm1 \n\t"
                         "addps %%xmm14, %%xmm2 \n\t"
                         "subps %%xmm14, %%xmm3 \n\t"
                         "movaps %%xmm4, %%xmm10 \n\t"
                         "movaps %%xmm2, %%xmm8 \n\t"
                         "movaps %%xmm3, %%xmm9 \n\t"
                         "shufps $0x4e, %%xmm3, %%xmm3 \n\t"
                         "shufps $0xb1, %%xmm10, %%xmm10 \n\t"
                         "shufps $0xb1, %%xmm8, %%xmm8 \n\t"
                         "shufps $0x1b, %%xmm9, %%xmm9"
                         :
                         :
                         "m" ((*m).u[6]),
                         "m" ((*m).u[7]),
                         "m" ((*m).u[16]),
                         "m" ((*m).u[17]),
                         "m" ((*s).c1.c1),
                         "m" ((*s).c1.c2)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm8", "xmm9", "xmm10");

   __asm__ __volatile__ ("mulps %%xmm0, %%xmm2 \n\t"
                         "mulps %%xmm1, %%xmm3 \n\t"
                         "mulps %%xmm1, %%xmm4 \n\t"
                         "movsd %0, %%xmm5 \n\t"
                         "movsd %2, %%xmm6 \n\t"
                         "movsd %4, %%xmm7"
                         :
                         :
                         "m" ((*m).u[10]),
                         "m" ((*m).u[11]),
                         "m" ((*m).u[12]),
                         "m" ((*m).u[13]),
                         "m" ((*m).u[14]),
                         "m" ((*m).u[15])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7");

   s+=4;
   _prefetch_weyl(s);
   s-=4;

   __asm__ __volatile__ ("mulps %%xmm1, %%xmm8 \n\t"
                         "mulps %%xmm0, %%xmm9 \n\t"
                         "mulps %%xmm0, %%xmm10 \n\t"
                         "movhps %0, %%xmm5 \n\t"
                         "movhps %2, %%xmm6 \n\t"
                         "movhps %4, %%xmm7 \n\t"
                         "addsubps %%xmm8, %%xmm2 \n\t"
                         "addsubps %%xmm9, %%xmm3 \n\t"
                         "addsubps %%xmm10, %%xmm4"
                         :
                         :
                         "m" ((*m).u[18]),
                         "m" ((*m).u[19]),
                         "m" ((*m).u[20]),
                         "m" ((*m).u[21]),
                         "m" ((*m).u[22]),
                         "m" ((*m).u[23])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7", "xmm8", "xmm9",
                         "xmm10");

   __asm__ __volatile__ ("movaps %%xmm5, %%xmm8 \n\t"
                         "movaps %%xmm6, %%xmm9 \n\t"
                         "movaps %%xmm7, %%xmm10 \n\t"
                         "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                         "shufps $0xb1, %%xmm4, %%xmm4 \n\t"
                         "mulps %%xmm1, %%xmm5 \n\t"
                         "mulps %%xmm1, %%xmm6 \n\t"
                         "mulps %%xmm1, %%xmm7 \n\t"
                         "shufps $0xb1, %%xmm8, %%xmm8 \n\t"
                         "shufps $0xb1, %%xmm9, %%xmm9 \n\t"
                         "shufps $0xb1, %%xmm10, %%xmm10 \n\t"
                         "mulps %%xmm0, %%xmm8 \n\t"
                         "mulps %%xmm0, %%xmm9 \n\t"
                         "mulps %%xmm0, %%xmm10"
                         :
                         :
                         :
                         "xmm3", "xmm4", "xmm5", "xmm6",
                         "xmm7", "xmm8", "xmm9", "xmm10");

   __asm__ __volatile__ ("movss %0, %%xmm13 \n\t"
                         "movaps %1, %%xmm11"
                         :
                         :
                         "m" ((*m).u[2]),
                         "m" ((*m).u[8]),
                         "m" ((*m).u[9]),
                         "m" ((*m).u[10]),
                         "m" ((*m).u[11])
                         :
                         "xmm11", "xmm13");

   __asm__ __volatile__ ("movaps %0, %%xmm12 \n\t"
                         "addsubps %%xmm8, %%xmm5 \n\t"
                         "addsubps %%xmm9, %%xmm6 \n\t"
                         "addsubps %%xmm10, %%xmm7 \n\t"
                         "movhps %4, %%xmm13 \n\t"
                         "movsldup %6, %%xmm0 \n\t"
                         "movshdup %6, %%xmm1 \n\t"
                         "addps %%xmm14, %%xmm13 \n\t"
                         "movaps %%xmm11, %%xmm8 \n\t"
                         "movaps %%xmm12, %%xmm9 \n\t"
                         "movaps %%xmm13, %%xmm10 \n\t"
                         "shufps $0xb1, %%xmm8, %%xmm8 \n\t"
                         "shufps $0xb1, %%xmm9, %%xmm9 \n\t"
                         "shufps $0xb1, %%xmm10, %%xmm10"
                         :
                         :
                         "m" ((*m).u[16]),
                         "m" ((*m).u[17]),
                         "m" ((*m).u[18]),
                         "m" ((*m).u[19]),
                         "m" ((*m).u[24]),
                         "m" ((*m).u[25]),
                         "m" ((*s).c1.c3),
                         "m" ((*s).c2.c1)
                         :
                         "xmm0", "xmm1", "xmm5", "xmm6",
                         "xmm7", "xmm8", "xmm9", "xmm10",
                         "xmm12", "xmm13");

   __asm__ __volatile__ ("mulps %%xmm0, %%xmm11 \n\t"
                         "mulps %%xmm0, %%xmm12 \n\t"
                         "mulps %%xmm0, %%xmm13 \n\t"
                         "addps %%xmm11, %%xmm2 \n\t"
                         "addps %%xmm12, %%xmm3 \n\t"
                         "addps %%xmm13, %%xmm4 \n\t"
                         "movss %0, %%xmm11 \n\t"
                         "movsd %1, %%xmm12 \n\t"
                         "movsd %3, %%xmm13 \n\t"
                         "mulps %%xmm1, %%xmm8 \n\t"
                         "mulps %%xmm1, %%xmm9 \n\t"
                         "mulps %%xmm1, %%xmm10"
                         :
                         :
                         "m" ((*m).u[3]),
                         "m" ((*m).u[26]),
                         "m" ((*m).u[27]),
                         "m" ((*m).u[28]),
                         "m" ((*m).u[29])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm8",
                         "xmm9", "xmm10", "xmm11", "xmm12",
                         "xmm13");

   __asm__ __volatile__ ("movhps %0, %%xmm11 \n\t"
                         "movhps %2, %%xmm12 \n\t"
                         "movhps %4, %%xmm13 \n\t"
                         "subps %%xmm14, %%xmm11 \n\t"
                         "addsubps %%xmm8, %%xmm2 \n\t"
                         "addsubps %%xmm9, %%xmm3 \n\t"
                         "addsubps %%xmm10, %%xmm4"
                         :
                         :
                         "m" ((*m).u[24]),
                         "m" ((*m).u[25]),
                         "m" ((*m).u[30]),
                         "m" ((*m).u[31]),
                         "m" ((*m).u[32]),
                         "m" ((*m).u[33])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm11",
                         "xmm12", "xmm13");

   __asm__ __volatile__ ("movaps %%xmm11, %%xmm8 \n\t"
                         "movaps %%xmm12, %%xmm9 \n\t"
                         "movaps %%xmm13, %%xmm10 \n\t"
                         "shufps $0x4e, %%xmm11, %%xmm11 \n\t"
                         "shufps $0x1b, %%xmm8, %%xmm8 \n\t"
                         "shufps $0xb1, %%xmm9, %%xmm9 \n\t"
                         "shufps $0xb1, %%xmm10, %%xmm10 \n\t"
                         "mulps %%xmm1, %%xmm11 \n\t"
                         "mulps %%xmm1, %%xmm12 \n\t"
                         "mulps %%xmm1, %%xmm13 \n\t"
                         "addps %%xmm11, %%xmm5 \n\t"
                         "addps %%xmm12, %%xmm6 \n\t"
                         "addps %%xmm13, %%xmm7 \n\t"
                         "mulps %%xmm0, %%xmm8 \n\t"
                         "mulps %%xmm0, %%xmm9 \n\t"
                         "mulps %%xmm0, %%xmm10"
                         :
                         :
                         :
                         "xmm5", "xmm6","xmm7", "xmm8",
                         "xmm9", "xmm10", "xmm11", "xmm12",
                         "xmm13");

   __asm__ __volatile__ ("movaps %0, %%xmm11 \n\t"
                         "movaps %4, %%xmm12"
                         :
                         :
                         "m" ((*m).u[12]),
                         "m" ((*m).u[13]),
                         "m" ((*m).u[14]),
                         "m" ((*m).u[15]),
                         "m" ((*m).u[20]),
                         "m" ((*m).u[21]),
                         "m" ((*m).u[22]),
                         "m" ((*m).u[23])
                         :
                         "xmm11", "xmm12");

   __asm__ __volatile__ ("movups %0, %%xmm13 \n\t"
                         "addsubps %%xmm8, %%xmm5 \n\t"
                         "addsubps %%xmm9, %%xmm6 \n\t"
                         "addsubps %%xmm10, %%xmm7 \n\t"
                         "movaps %%xmm11, %%xmm8 \n\t"
                         "movaps %%xmm12, %%xmm9 \n\t"
                         "movaps %%xmm13, %%xmm10 \n\t"
                         "movsldup %4, %%xmm0 \n\t"
                         "movshdup %4, %%xmm1 \n\t"
                         "shufps $0xb1, %%xmm8, %%xmm8 \n\t"
                         "shufps $0xb1, %%xmm9, %%xmm9 \n\t"
                         "shufps $0xb1, %%xmm10, %%xmm10"
                         :
                         :
                         "m" ((*m).u[26]),
                         "m" ((*m).u[27]),
                         "m" ((*m).u[28]),
                         "m" ((*m).u[29]),
                         "m" ((*s).c2.c2),
                         "m" ((*s).c2.c3)
                         :
                         "xmm0", "xmm1", "xmm5", "xmm6",
                         "xmm7", "xmm8", "xmm9", "xmm10",
                         "xmm13");

   __asm__ __volatile__ ("mulps %%xmm0, %%xmm11 \n\t"
                         "mulps %%xmm0, %%xmm12 \n\t"
                         "mulps %%xmm0, %%xmm13 \n\t"
                         "addps %%xmm11, %%xmm2 \n\t"
                         "addps %%xmm12, %%xmm3 \n\t"
                         "addps %%xmm13, %%xmm4 \n\t"
                         "movups %0, %%xmm11 \n\t"
                         "movss %4, %%xmm12 \n\t"
                         "movss %5, %%xmm13 \n\t"
                         "mulps %%xmm1, %%xmm8 \n\t"
                         "mulps %%xmm1, %%xmm9 \n\t"
                         "mulps %%xmm1, %%xmm10"
                         :
                         :
                         "m" ((*m).u[30]),
                         "m" ((*m).u[31]),
                         "m" ((*m).u[32]),
                         "m" ((*m).u[33]),
                         "m" ((*m).u[4]),
                         "m" ((*m).u[5])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm8",
                         "xmm9", "xmm10", "xmm11", "xmm12",
                         "xmm13");

   __asm__ __volatile__ ("movhps %0, %%xmm12 \n\t"
                         "movhps %0, %%xmm13 \n\t"
                         "addsubps %%xmm8, %%xmm2 \n\t"
                         "addsubps %%xmm9, %%xmm3 \n\t"
                         "addps %%xmm14, %%xmm12 \n\t"
                         "subps %%xmm14, %%xmm13 \n\t"
                         "addsubps %%xmm10, %%xmm4"
                         :
                         :
                         "m" ((*m).u[34]),
                         "m" ((*m).u[35])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm12",
                         "xmm13");

   __asm__ __volatile__ ("movaps %%xmm11, %%xmm8 \n\t"
                         "movaps %%xmm12, %%xmm9 \n\t"
                         "movaps %%xmm13, %%xmm10 \n\t"
                         "shufps $0xb1, %%xmm8, %%xmm8 \n\t"
                         "shufps $0xb1, %%xmm9, %%xmm9 \n\t"
                         "shufps $0x1b, %%xmm10, %%xmm10 \n\t"
                         "mulps %%xmm1, %%xmm8 \n\t"
                         "mulps %%xmm1, %%xmm9 \n\t"
                         "mulps %%xmm0, %%xmm10 \n\t"
                         "shufps $0x4e, %%xmm13, %%xmm13 \n\t"
                         "shufps $0xb1, %%xmm5, %%xmm5 \n\t"
                         "shufps $0xb1, %%xmm6, %%xmm6 \n\t"
                         "mulps %%xmm0, %%xmm11 \n\t"
                         "mulps %%xmm0, %%xmm12 \n\t"
                         "mulps %%xmm1, %%xmm13"
                         :
                         :
                         :
                         "xmm5", "xmm6", "xmm8", "xmm9",
                         "xmm10", "xmm11", "xmm12", "xmm13");

   __asm__ __volatile__ ("addsubps %%xmm8, %%xmm5 \n\t"
                         "addsubps %%xmm9, %%xmm6 \n\t"
                         "addsubps %%xmm10, %%xmm7 \n\t"
                         "shufps $0xd8, %%xmm2, %%xmm2 \n\t"
                         "shufps $0xd8, %%xmm3, %%xmm3 \n\t"
                         "shufps $0xd8, %%xmm4, %%xmm4 \n\t"
                         "addps %%xmm11, %%xmm5 \n\t"
                         "addps %%xmm12, %%xmm6 \n\t"
                         "addps %%xmm13, %%xmm7 \n\t"
                         "shufps $0xd8, %%xmm5, %%xmm5 \n\t"
                         "shufps $0xd8, %%xmm6, %%xmm6 \n\t"
                         "shufps $0x8d, %%xmm7, %%xmm7 \n\t"
                         "haddps %%xmm3, %%xmm2 \n\t"
                         "haddps %%xmm5, %%xmm4 \n\t"
                         "haddps %%xmm7, %%xmm6 \n\t"
                         "movaps %%xmm2, %0 \n\t"
                         "movaps %%xmm4, %2 \n\t"
                         "movaps %%xmm6, %4"
                         :
                         "=m" ((*r).c1.c1),
                         "=m" ((*r).c1.c2),
                         "=m" ((*r).c1.c3),
                         "=m" ((*r).c2.c1),
                         "=m" ((*r).c2.c2),
                         "=m" ((*r).c2.c3)
                         :
                         :
                         "xmm2", "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7");
}

#if (defined AVX)
#include "avx.h"

void mul_pauli2(float mu,pauli *m,spinor *s,spinor *r)
{
   m+=4;
   _prefetch_pauli_dble(m);
   m-=4;

   __asm__ __volatile__ ("vmovss %0, %%xmm14 \n\t"
                         "vmovss %1, %%xmm2 \n\t"
                         "vmovss %2, %%xmm3 \n\t"
                         "vmovsd %3, %%xmm4 \n\t"
                         "vpermilps $0xb1, %%xmm14, %%xmm14 \n\t"
                         "vxorps %%xmm15, %%xmm15, %%xmm15 \n\t"
                         "vmovss %5, %%xmm8 \n\t"
                         "vmovss %6, %%xmm9 \n\t"
                         "vmovsd %7, %%xmm10 \n\t"
                         "vsubps %%xmm14, %%xmm15, %%xmm15"
                         :
                         :
                         "m" (mu),
                         "m" (m[0].u[0]),
                         "m" (m[0].u[1]),
                         "m" (m[0].u[8]),
                         "m" (m[0].u[9]),
                         "m" (m[1].u[0]),
                         "m" (m[1].u[1]),
                         "m" (m[1].u[8]),
                         "m" (m[1].u[9])
                         :
                         "xmm2", "xmm3", "xmm4",  "xmm8",
                         "xmm9", "xmm10", "xmm14", "xmm15");

   __asm__ __volatile__ ("vmovhps %0, %%xmm2, %%xmm2 \n\t"
                         "vmovhps %0, %%xmm3, %%xmm3 \n\t"
                         "vmovhps %2, %%xmm4, %%xmm4 \n\t"
                         "vmovaps %4, %%xmm0"
                         :
                         :
                         "m" (m[0].u[6]),
                         "m" (m[0].u[7]),
                         "m" (m[0].u[16]),
                         "m" (m[0].u[17]),
                         "m" ((*s).c1.c1),
                         "m" ((*s).c1.c2)
                         :
                         "xmm0", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("vmovhps %0, %%xmm8, %%xmm8 \n\t"
                         "vmovhps %0, %%xmm9, %%xmm9 \n\t"
                         "vmovhps %2, %%xmm10, %%xmm10 \n\t"
                         "vinsertf128 $0x1, %4, %%ymm0, %%ymm0 \n\t"
                         "vinsertf128 $0x1, %%xmm15, %%ymm14, %%ymm14 \n\t"
                         "vinsertf128 $0x1, %%xmm8, %%ymm2, %%ymm2 \n\t"
                         "vinsertf128 $0x1, %%xmm9, %%ymm3, %%ymm3 \n\t"
                         "vinsertf128 $0x1, %%xmm10, %%ymm4, %%ymm4 \n\t"
                         "vmovshdup %%ymm0, %%ymm1 \n\t"
                         "vmovsldup %%ymm0, %%ymm0"
                         :
                         :
                         "m" (m[1].u[6]),
                         "m" (m[1].u[7]),
                         "m" (m[1].u[16]),
                         "m" (m[1].u[17]),
                         "m" ((*s).c3.c1),
                         "m" ((*s).c3.c2)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm8", "xmm9", "xmm10",
                         "xmm14");

   __asm__ __volatile__ ("vaddps %%ymm14, %%ymm2, %%ymm2 \n\t"
                         "vsubps %%ymm14, %%ymm3, %%ymm3 \n\t"
                         "vpermilps $0xb1, %%ymm4, %%ymm10 \n\t"
                         "vpermilps $0x4e, %%ymm3, %%ymm3 \n\t"
                         "vpermilps $0xb1, %%ymm2, %%ymm8 \n\t"
                         "vpermilps $0xb1, %%ymm3, %%ymm9"
                         :
                         :
                         :
                         "xmm2", "xmm3", "xmm8", "xmm9",
                         "xmm10");

   __asm__ __volatile__ ("vmulps %%ymm0, %%ymm2, %%ymm2 \n\t"
                         "vmulps %%ymm1, %%ymm3, %%ymm3 \n\t"
                         "vmulps %%ymm1, %%ymm4, %%ymm4 \n\t"
                         "vmovsd %0, %%xmm5 \n\t"
                         "vmovsd %2, %%xmm6 \n\t"
                         "vmovsd %4, %%xmm7"
                         :
                         :
                         "m" (m[0].u[10]),
                         "m" (m[0].u[11]),
                         "m" (m[0].u[12]),
                         "m" (m[0].u[13]),
                         "m" (m[0].u[14]),
                         "m" (m[0].u[15])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7");

   __asm__ __volatile__ ("vmovsd %0, %%xmm11 \n\t"
                         "vmovsd %2, %%xmm12 \n\t"
                         "vmovsd %4, %%xmm13 \n\t"
                         "vmulps %%ymm1, %%ymm8, %%ymm8 \n\t"
                         "vmulps %%ymm0, %%ymm9, %%ymm9 \n\t"
                         "vmulps %%ymm0, %%ymm10, %%ymm10"
                         :
                         :
                         "m" (m[1].u[10]),
                         "m" (m[1].u[11]),
                         "m" (m[1].u[12]),
                         "m" (m[1].u[13]),
                         "m" (m[1].u[14]),
                         "m" (m[1].u[15])
                         :
                         "xmm8", "xmm9", "xmm10", "xmm11",
                         "xmm12", "xmm13");

   __asm__ __volatile__ ("vmovhps %0, %%xmm5, %%xmm5 \n\t"
                         "vmovhps %2, %%xmm6, %%xmm6 \n\t"
                         "vmovhps %4, %%xmm7, %%xmm7 \n\t"
                         "vaddsubps %%ymm8, %%ymm2, %%ymm2 \n\t"
                         "vaddsubps %%ymm9, %%ymm3, %%ymm3 \n\t"
                         "vaddsubps %%ymm10, %%ymm4, %%ymm4"
                         :
                         :
                         "m" (m[0].u[18]),
                         "m" (m[0].u[19]),
                         "m" (m[0].u[20]),
                         "m" (m[0].u[21]),
                         "m" (m[0].u[22]),
                         "m" (m[0].u[23])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7");

   __asm__ __volatile__ ("vmovhps %0, %%xmm11, %%xmm11 \n\t"
                         "vmovhps %2, %%xmm12, %%xmm12 \n\t"
                         "vmovhps %4, %%xmm13, %%xmm13 \n\t"
                         "vinsertf128 $0x1, %%xmm11, %%ymm5, %%ymm5 \n\t"
                         "vinsertf128 $0x1, %%xmm12, %%ymm6, %%ymm6 \n\t"
                         "vinsertf128 $0x1, %%xmm13, %%ymm7, %%ymm7"
                         :
                         :
                         "m" (m[1].u[18]),
                         "m" (m[1].u[19]),
                         "m" (m[1].u[20]),
                         "m" (m[1].u[21]),
                         "m" (m[1].u[22]),
                         "m" (m[1].u[23])
                         :
                         "xmm5", "xmm6", "xmm7", "xmm11",
                         "xmm12", "xmm13");

   s+=4;
   _prefetch_spinor(s);
   s-=4;

   __asm__ __volatile__ ("vpermilps $0xb1, %%ymm5, %%ymm8 \n\t"
                         "vpermilps $0xb1, %%ymm6, %%ymm9 \n\t"
                         "vpermilps $0xb1, %%ymm7, %%ymm10 \n\t"
                         "vpermilps $0xb1, %%ymm3, %%ymm3 \n\t"
                         "vpermilps $0xb1, %%ymm4, %%ymm4 \n\t"
                         "vmulps %%ymm1, %%ymm5, %%ymm5 \n\t"
                         "vmulps %%ymm1, %%ymm6, %%ymm6 \n\t"
                         "vmulps %%ymm1, %%ymm7, %%ymm7 \n\t"
                         "vmulps %%ymm0, %%ymm8, %%ymm8 \n\t"
                         "vmulps %%ymm0, %%ymm9, %%ymm9 \n\t"
                         "vmulps %%ymm0, %%ymm10, %%ymm10"
                         :
                         :
                         :
                         "xmm3", "xmm4", "xmm5", "xmm6",
                         "xmm7", "xmm8", "xmm9", "xmm10");

   __asm__ __volatile__ ("vmovaps %0, %%xmm11 \n\t"
                         "vmovaps %4, %%xmm12 \n\t"
                         "vmovss %8, %%xmm13 \n\t"
                         "vaddsubps %%ymm8, %%ymm5, %%ymm5 \n\t"
                         "vaddsubps %%ymm9, %%ymm6, %%ymm6 \n\t"
                         "vaddsubps %%ymm10, %%ymm7, %%ymm7"
                         :
                         :
                         "m" (m[0].u[8]),
                         "m" (m[0].u[9]),
                         "m" (m[0].u[10]),
                         "m" (m[0].u[11]),
                         "m" (m[0].u[16]),
                         "m" (m[0].u[17]),
                         "m" (m[0].u[18]),
                         "m" (m[0].u[19]),
                         "m" (m[0].u[2])
                         :
                         "xmm5", "xmm6", "xmm7", "xmm11",
                         "xmm12", "xmm13");

   __asm__ __volatile__ ("vmovaps %0, %%xmm8 \n\t"
                         "vmovaps %4, %%xmm9 \n\t"
                         "vmovss %8, %%xmm10"
                         :
                         :
                         "m" (m[1].u[8]),
                         "m" (m[1].u[9]),
                         "m" (m[1].u[10]),
                         "m" (m[1].u[11]),
                         "m" (m[1].u[16]),
                         "m" (m[1].u[17]),
                         "m" (m[1].u[18]),
                         "m" (m[1].u[19]),
                         "m" (m[1].u[2])
                         :
                         "xmm8", "xmm9", "xmm10");

   __asm__ __volatile__ ("vmovaps %0, %%xmm0 \n\t"
                         "vmovhps %2, %%xmm13, %%xmm13 \n\t"
                         "vmovhps %4, %%xmm10, %%xmm10 \n\t"
                         "vinsertf128 $0x1, %6, %%ymm0, %%ymm0 \n\t"
                         "vinsertf128 $0x1, %%xmm8, %%ymm11, %%ymm11 \n\t"
                         "vinsertf128 $0x1, %%xmm9, %%ymm12, %%ymm12 \n\t"
                         "vinsertf128 $0x1, %%xmm10, %%ymm13, %%ymm13 \n\t"
                         "vmovshdup %%ymm0, %%ymm1 \n\t"
                         "vmovsldup %%ymm0, %%ymm0 \n\t"
                         "vaddps %%ymm14, %%ymm13, %%ymm13 \n\t"
                         "vpermilps $0xb1, %%ymm11, %%ymm8 \n\t"
                         "vpermilps $0xb1, %%ymm12, %%ymm9 \n\t"
                         "vpermilps $0xb1, %%ymm13, %%ymm10"
                         :
                         :
                         "m" ((*s).c1.c3),
                         "m" ((*s).c2.c1),
                         "m" (m[0].u[24]),
                         "m" (m[0].u[25]),
                         "m" (m[1].u[24]),
                         "m" (m[1].u[25]),
                         "m" ((*s).c3.c3),
                         "m" ((*s).c4.c1)
                         :
                         "xmm0", "xmm1", "xmm8", "xmm9",
                         "xmm10", "xmm11", "xmm12", "xmm13");

   __asm__ __volatile__ ("vmulps %%ymm0, %%ymm11, %%ymm11 \n\t"
                         "vmulps %%ymm0, %%ymm12, %%ymm12 \n\t"
                         "vmulps %%ymm0, %%ymm13, %%ymm13 \n\t"
                         "vaddps %%ymm11, %%ymm2, %%ymm2 \n\t"
                         "vaddps %%ymm12, %%ymm3, %%ymm3 \n\t"
                         "vaddps %%ymm13, %%ymm4, %%ymm4 \n\t"
                         "vmovss %0, %%xmm11 \n\t"
                         "vmovsd %1, %%xmm12 \n\t"
                         "vmovsd %3, %%xmm13 \n\t"
                         "vmulps %%ymm1, %%ymm8, %%ymm8 \n\t"
                         "vmulps %%ymm1, %%ymm9, %%ymm9 \n\t"
                         "vmulps %%ymm1, %%ymm10, %%ymm10"
                         :
                         :
                         "m" (m[0].u[3]),
                         "m" (m[0].u[26]),
                         "m" (m[0].u[27]),
                         "m" (m[0].u[28]),
                         "m" (m[0].u[29])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm8",
                         "xmm9", "xmm10", "xmm11", "xmm12",
                         "xmm13");

   __asm__ __volatile__ ("vmovhps %0, %%xmm11, %%xmm11 \n\t"
                         "vmovhps %2, %%xmm12, %%xmm12 \n\t"
                         "vmovhps %4, %%xmm13, %%xmm13 \n\t"
                         "vaddsubps %%ymm8, %%ymm2, %%ymm2 \n\t"
                         "vaddsubps %%ymm9, %%ymm3, %%ymm3 \n\t"
                         "vaddsubps %%ymm10, %%ymm4, %%ymm4"
                         :
                         :
                         "m" (m[0].u[24]),
                         "m" (m[0].u[25]),
                         "m" (m[0].u[30]),
                         "m" (m[0].u[31]),
                         "m" (m[0].u[32]),
                         "m" (m[0].u[33])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm11",
                         "xmm12", "xmm13");

   __asm__ __volatile__ ("vmovss %0, %%xmm8 \n\t"
                         "vmovsd %1, %%xmm9 \n\t"
                         "vmovsd %3, %%xmm10"
                         :
                         :
                         "m" (m[1].u[3]),
                         "m" (m[1].u[26]),
                         "m" (m[1].u[27]),
                         "m" (m[1].u[28]),
                         "m" (m[1].u[29])
                         :
                         "xmm8", "xmm9", "xmm10");

   __asm__ __volatile__ ("vmovhps %0, %%xmm8, %%xmm8 \n\t"
                         "vmovhps %2, %%xmm9, %%xmm9 \n\t"
                         "vmovhps %4, %%xmm10, %%xmm10 \n\t"
                         "vinsertf128 $0x1, %%xmm8, %%ymm11, %%ymm11 \n\t"
                         "vinsertf128 $0x1, %%xmm9, %%ymm12, %%ymm12 \n\t"
                         "vinsertf128 $0x1, %%xmm10, %%ymm13, %%ymm13 \n\t"
                         "vsubps %%ymm14, %%ymm11, %%ymm11"
                         :
                         :
                         "m" (m[1].u[24]),
                         "m" (m[1].u[25]),
                         "m" (m[1].u[30]),
                         "m" (m[1].u[31]),
                         "m" (m[1].u[32]),
                         "m" (m[1].u[33])
                         :
                         "xmm8", "xmm9", "xmm10", "xmm11",
                         "xmm12", "xmm13");

   __asm__ __volatile__ ("vpermilps $0x4e, %%ymm11, %%ymm11 \n\t"
                         "vpermilps $0xb1, %%ymm11, %%ymm8 \n\t"
                         "vpermilps $0xb1, %%ymm12, %%ymm9 \n\t"
                         "vpermilps $0xb1, %%ymm13, %%ymm10 \n\t"
                         "vmulps %%ymm1, %%ymm11, %%ymm11 \n\t"
                         "vmulps %%ymm1, %%ymm12, %%ymm12 \n\t"
                         "vmulps %%ymm1, %%ymm13, %%ymm13 \n\t"
                         "vaddps %%ymm11, %%ymm5, %%ymm5 \n\t"
                         "vaddps %%ymm12, %%ymm6, %%ymm6 \n\t"
                         "vaddps %%ymm13, %%ymm7, %%ymm7 \n\t"
                         "vmulps %%ymm0, %%ymm8, %%ymm8 \n\t"
                         "vmulps %%ymm0, %%ymm9, %%ymm9 \n\t"
                         "vmulps %%ymm0, %%ymm10, %%ymm10"
                         :
                         :
                         :
                         "xmm5", "xmm6", "xmm7", "xmm8",
                         "xmm9", "xmm10", "xmm11", "xmm12",
                         "xmm13");

   __asm__ __volatile__ ("vmovaps %0, %%xmm11 \n\t"
                         "vmovaps %4, %%xmm12 \n\t"
                         "vaddsubps %%ymm8, %%ymm5, %%ymm5 \n\t"
                         "vaddsubps %%ymm9, %%ymm6, %%ymm6 \n\t"
                         "vaddsubps %%ymm10, %%ymm7, %%ymm7"
                         :
                         :
                         "m" (m[0].u[12]),
                         "m" (m[0].u[13]),
                         "m" (m[0].u[14]),
                         "m" (m[0].u[15]),
                         "m" (m[0].u[20]),
                         "m" (m[0].u[21]),
                         "m" (m[0].u[22]),
                         "m" (m[0].u[23])
                         :
                         "xmm5", "xmm6", "xmm7", "xmm11",
                         "xmm12");

   __asm__ __volatile__ ("vmovups %0, %%xmm13 \n\t"
                         "vmovaps %4, %%xmm0"
                         :
                         :
                         "m" (m[0].u[26]),
                         "m" (m[0].u[27]),
                         "m" (m[0].u[28]),
                         "m" (m[0].u[29]),
                         "m" ((*s).c2.c2),
                         "m" ((*s).c2.c3)
                         :
                         "xmm0", "xmm13");

   __asm__ __volatile__ ("vmovaps %0, %%xmm8 \n\t"
                         "vmovaps %4, %%xmm9"
                         :
                         :
                         "m" (m[1].u[12]),
                         "m" (m[1].u[13]),
                         "m" (m[1].u[14]),
                         "m" (m[1].u[15]),
                         "m" (m[1].u[20]),
                         "m" (m[1].u[21]),
                         "m" (m[1].u[22]),
                         "m" (m[1].u[23])
                         :
                         "xmm8", "xmm9");

   __asm__ __volatile__ ("vmovups %0, %%xmm10 \n\t"
                         "vinsertf128 $0x1, %4, %%ymm0, %%ymm0 \n\t"
                         "vinsertf128 $0x1, %%xmm8, %%ymm11, %%ymm11 \n\t"
                         "vinsertf128 $0x1, %%xmm9, %%ymm12, %%ymm12 \n\t"
                         "vinsertf128 $0x1, %%xmm10, %%ymm13, %%ymm13 \n\t"
                         "vmovshdup %%ymm0, %%ymm1 \n\t"
                         "vmovsldup %%ymm0, %%ymm0 \n\t"
                         "vpermilps $0xb1, %%ymm11, %%ymm8 \n\t"
                         "vpermilps $0xb1, %%ymm12, %%ymm9 \n\t"
                         "vpermilps $0xb1, %%ymm13, %%ymm10"
                         :
                         :
                         "m" (m[1].u[26]),
                         "m" (m[1].u[27]),
                         "m" (m[1].u[28]),
                         "m" (m[1].u[29]),
                         "m" ((*s).c4.c2),
                         "m" ((*s).c4.c3)
                         :
                         "xmm0", "xmm1", "xmm8", "xmm9",
                         "xmm10", "xmm11", "xmm12", "xmm13");

   __asm__ __volatile__ ("vmulps %%ymm0, %%ymm11, %%ymm11 \n\t"
                         "vmulps %%ymm0, %%ymm12, %%ymm12 \n\t"
                         "vmulps %%ymm0, %%ymm13, %%ymm13 \n\t"
                         "vaddps %%ymm11, %%ymm2, %%ymm2 \n\t"
                         "vaddps %%ymm12, %%ymm3, %%ymm3 \n\t"
                         "vaddps %%ymm13, %%ymm4, %%ymm4 \n\t"
                         "vmovups %0, %%xmm11 \n\t"
                         "vmovss %4, %%xmm12 \n\t"
                         "vmovss %5, %%xmm13 \n\t"
                         "vmulps %%ymm1, %%ymm8, %%ymm8 \n\t"
                         "vmulps %%ymm1, %%ymm9, %%ymm9 \n\t"
                         "vmulps %%ymm1, %%ymm10, %%ymm10"
                         :
                         :
                         "m" (m[0].u[30]),
                         "m" (m[0].u[31]),
                         "m" (m[0].u[32]),
                         "m" (m[0].u[33]),
                         "m" (m[0].u[4]),
                         "m" (m[0].u[5])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm8",
                         "xmm9", "xmm10", "xmm11", "xmm12",
                         "xmm13");

   __asm__ __volatile__ ("vmovhps %0, %%xmm12, %%xmm12 \n\t"
                         "vmovhps %0, %%xmm13, %%xmm13 \n\t"
                         "vaddsubps %%ymm8, %%ymm2, %%ymm2 \n\t"
                         "vaddsubps %%ymm9, %%ymm3, %%ymm3 \n\t"
                         "vaddsubps %%ymm10, %%ymm4, %%ymm4"
                         :
                         :
                         "m" (m[0].u[34]),
                         "m" (m[0].u[35])
                         :
                         "xmm2", "xmm3", "xmm4", "xmm12",
                         "xmm13");

   __asm__ __volatile__ ("vmovss %0, %%xmm9 \n\t"
                         "vmovss %1, %%xmm10 \n\t"
                         "vmovups %2, %%xmm8 \n\t"
                         "vmovhps %6, %%xmm9, %%xmm9 \n\t"
                         "vmovhps %6, %%xmm10, %%xmm10 \n\t"
                         "vinsertf128 $0x1, %%xmm8, %%ymm11, %%ymm11 \n\t"
                         "vinsertf128 $0x1, %%xmm9, %%ymm12, %%ymm12 \n\t"
                         "vinsertf128 $0x1, %%xmm10, %%ymm13, %%ymm13 \n\t"
                         "vaddps %%ymm14, %%ymm12, %%ymm12 \n\t"
                         "vsubps %%ymm14, %%ymm13, %%ymm13"
                         :
                         :
                         "m" (m[1].u[4]),
                         "m" (m[1].u[5]),
                         "m" (m[1].u[30]),
                         "m" (m[1].u[31]),
                         "m" (m[1].u[32]),
                         "m" (m[1].u[33]),
                         "m" (m[1].u[34]),
                         "m" (m[1].u[35])
                         :
                         "xmm8", "xmm9", "xmm10", "xmm11",
                         "xmm12", "xmm13");

   __asm__ __volatile__ ("vpermilps $0xb1, %%ymm5, %%ymm5 \n\t"
                         "vpermilps $0xb1, %%ymm6, %%ymm6 \n\t"
                         "vpermilps $0x4e, %%ymm13, %%ymm13 \n\t"
                         "vpermilps $0xb1, %%ymm11, %%ymm8 \n\t"
                         "vpermilps $0xb1, %%ymm12, %%ymm9 \n\t"
                         "vpermilps $0xb1, %%ymm13, %%ymm10 \n\t"
                         "vmulps %%ymm1, %%ymm8, %%ymm8 \n\t"
                         "vmulps %%ymm1, %%ymm9, %%ymm9 \n\t"
                         "vmulps %%ymm0, %%ymm10, %%ymm10 \n\t"
                         "vmulps %%ymm0, %%ymm11, %%ymm11 \n\t"
                         "vmulps %%ymm0, %%ymm12, %%ymm12 \n\t"
                         "vmulps %%ymm1, %%ymm13, %%ymm13"
                         :
                         :
                         :
                         "xmm5", "xmm6", "xmm8", "xmm9",
                         "xmm10", "xmm11", "xmm12", "xmm13");

   __asm__ __volatile__ ("vaddsubps %%ymm8, %%ymm5, %%ymm5 \n\t"
                         "vaddsubps %%ymm9, %%ymm6, %%ymm6 \n\t"
                         "vaddsubps %%ymm10, %%ymm7, %%ymm7 \n\t"
                         "vpermilps $0xd8, %%ymm2, %%ymm2 \n\t"
                         "vpermilps $0xd8, %%ymm3, %%ymm3 \n\t"
                         "vpermilps $0xd8, %%ymm4, %%ymm4 \n\t"
                         "vaddps %%ymm11, %%ymm5, %%ymm5 \n\t"
                         "vaddps %%ymm12, %%ymm6, %%ymm6 \n\t"
                         "vaddps %%ymm13, %%ymm7, %%ymm7 \n\t"
                         "vpermilps $0xd8, %%ymm5, %%ymm5 \n\t"
                         "vpermilps $0xd8, %%ymm6, %%ymm6 \n\t"
                         "vpermilps $0x8d, %%ymm7, %%ymm7 \n\t"
                         "vhaddps %%ymm3, %%ymm2, %%ymm2 \n\t"
                         "vhaddps %%ymm5, %%ymm4, %%ymm4 \n\t"
                         "vhaddps %%ymm7, %%ymm6, %%ymm6 \n\t"
                         "vmovaps %%xmm2, %0 \n\t"
                         "vmovaps %%xmm4, %2 \n\t"
                         "vmovaps %%xmm6, %4"
                         :
                         "=m" ((*r).c1.c1),
                         "=m" ((*r).c1.c2),
                         "=m" ((*r).c1.c3),
                         "=m" ((*r).c2.c1),
                         "=m" ((*r).c2.c2),
                         "=m" ((*r).c2.c3)
                         :
                         :
                         "xmm2", "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7");

   __asm__ __volatile__ ("vextractf128 $0x1, %%ymm2, %0 \n\t"
                         "vextractf128 $0x1, %%ymm4, %2 \n\t"
                         "vextractf128 $0x1, %%ymm6, %4"
                         :
                         "=m" ((*r).c3.c1),
                         "=m" ((*r).c3.c2),
                         "=m" ((*r).c3.c3),
                         "=m" ((*r).c4.c1),
                         "=m" ((*r).c4.c2),
                         "=m" ((*r).c4.c3));

   _avx_zeroupper();
}

#else

void mul_pauli2(float mu,pauli *m,spinor *s,spinor *r)
{
   spin_t *ps,*pr;

   ps=(spin_t*)(s);
   pr=(spin_t*)(r);

   mul_pauli(mu,m,(*ps).w,(*pr).w);
   mul_pauli(-mu,m+1,(*ps).w+1,(*pr).w+1);
}

#endif
#else

static weyl rs;


void mul_pauli(float mu,pauli *m,weyl *s,weyl *r)
{
   float *u;

   u=(*m).u;

   rs.c1.c1.re=
      u[ 0]*(*s).c1.c1.re-   mu*(*s).c1.c1.im+
      u[ 6]*(*s).c1.c2.re-u[ 7]*(*s).c1.c2.im+
      u[ 8]*(*s).c1.c3.re-u[ 9]*(*s).c1.c3.im+
      u[10]*(*s).c2.c1.re-u[11]*(*s).c2.c1.im+
      u[12]*(*s).c2.c2.re-u[13]*(*s).c2.c2.im+
      u[14]*(*s).c2.c3.re-u[15]*(*s).c2.c3.im;

   rs.c1.c1.im=
      u[ 0]*(*s).c1.c1.im+   mu*(*s).c1.c1.re+
      u[ 6]*(*s).c1.c2.im+u[ 7]*(*s).c1.c2.re+
      u[ 8]*(*s).c1.c3.im+u[ 9]*(*s).c1.c3.re+
      u[10]*(*s).c2.c1.im+u[11]*(*s).c2.c1.re+
      u[12]*(*s).c2.c2.im+u[13]*(*s).c2.c2.re+
      u[14]*(*s).c2.c3.im+u[15]*(*s).c2.c3.re;

   rs.c1.c2.re=
      u[ 6]*(*s).c1.c1.re+u[ 7]*(*s).c1.c1.im+
      u[ 1]*(*s).c1.c2.re-   mu*(*s).c1.c2.im+
      u[16]*(*s).c1.c3.re-u[17]*(*s).c1.c3.im+
      u[18]*(*s).c2.c1.re-u[19]*(*s).c2.c1.im+
      u[20]*(*s).c2.c2.re-u[21]*(*s).c2.c2.im+
      u[22]*(*s).c2.c3.re-u[23]*(*s).c2.c3.im;

   rs.c1.c2.im=
      u[ 6]*(*s).c1.c1.im-u[ 7]*(*s).c1.c1.re+
      u[ 1]*(*s).c1.c2.im+   mu*(*s).c1.c2.re+
      u[16]*(*s).c1.c3.im+u[17]*(*s).c1.c3.re+
      u[18]*(*s).c2.c1.im+u[19]*(*s).c2.c1.re+
      u[20]*(*s).c2.c2.im+u[21]*(*s).c2.c2.re+
      u[22]*(*s).c2.c3.im+u[23]*(*s).c2.c3.re;

   rs.c1.c3.re=
      u[ 8]*(*s).c1.c1.re+u[ 9]*(*s).c1.c1.im+
      u[16]*(*s).c1.c2.re+u[17]*(*s).c1.c2.im+
      u[ 2]*(*s).c1.c3.re-   mu*(*s).c1.c3.im+
      u[24]*(*s).c2.c1.re-u[25]*(*s).c2.c1.im+
      u[26]*(*s).c2.c2.re-u[27]*(*s).c2.c2.im+
      u[28]*(*s).c2.c3.re-u[29]*(*s).c2.c3.im;

   rs.c1.c3.im=
      u[ 8]*(*s).c1.c1.im-u[ 9]*(*s).c1.c1.re+
      u[16]*(*s).c1.c2.im-u[17]*(*s).c1.c2.re+
      u[ 2]*(*s).c1.c3.im+   mu*(*s).c1.c3.re+
      u[24]*(*s).c2.c1.im+u[25]*(*s).c2.c1.re+
      u[26]*(*s).c2.c2.im+u[27]*(*s).c2.c2.re+
      u[28]*(*s).c2.c3.im+u[29]*(*s).c2.c3.re;

   rs.c2.c1.re=
      u[10]*(*s).c1.c1.re+u[11]*(*s).c1.c1.im+
      u[18]*(*s).c1.c2.re+u[19]*(*s).c1.c2.im+
      u[24]*(*s).c1.c3.re+u[25]*(*s).c1.c3.im+
      u[ 3]*(*s).c2.c1.re-   mu*(*s).c2.c1.im+
      u[30]*(*s).c2.c2.re-u[31]*(*s).c2.c2.im+
      u[32]*(*s).c2.c3.re-u[33]*(*s).c2.c3.im;

   rs.c2.c1.im=
      u[10]*(*s).c1.c1.im-u[11]*(*s).c1.c1.re+
      u[18]*(*s).c1.c2.im-u[19]*(*s).c1.c2.re+
      u[24]*(*s).c1.c3.im-u[25]*(*s).c1.c3.re+
      u[ 3]*(*s).c2.c1.im+   mu*(*s).c2.c1.re+
      u[30]*(*s).c2.c2.im+u[31]*(*s).c2.c2.re+
      u[32]*(*s).c2.c3.im+u[33]*(*s).c2.c3.re;

   rs.c2.c2.re=
      u[12]*(*s).c1.c1.re+u[13]*(*s).c1.c1.im+
      u[20]*(*s).c1.c2.re+u[21]*(*s).c1.c2.im+
      u[26]*(*s).c1.c3.re+u[27]*(*s).c1.c3.im+
      u[30]*(*s).c2.c1.re+u[31]*(*s).c2.c1.im+
      u[ 4]*(*s).c2.c2.re-   mu*(*s).c2.c2.im+
      u[34]*(*s).c2.c3.re-u[35]*(*s).c2.c3.im;

   rs.c2.c2.im=
      u[12]*(*s).c1.c1.im-u[13]*(*s).c1.c1.re+
      u[20]*(*s).c1.c2.im-u[21]*(*s).c1.c2.re+
      u[26]*(*s).c1.c3.im-u[27]*(*s).c1.c3.re+
      u[30]*(*s).c2.c1.im-u[31]*(*s).c2.c1.re+
      u[ 4]*(*s).c2.c2.im+   mu*(*s).c2.c2.re+
      u[34]*(*s).c2.c3.im+u[35]*(*s).c2.c3.re;

   rs.c2.c3.re=
      u[14]*(*s).c1.c1.re+u[15]*(*s).c1.c1.im+
      u[22]*(*s).c1.c2.re+u[23]*(*s).c1.c2.im+
      u[28]*(*s).c1.c3.re+u[29]*(*s).c1.c3.im+
      u[32]*(*s).c2.c1.re+u[33]*(*s).c2.c1.im+
      u[34]*(*s).c2.c2.re+u[35]*(*s).c2.c2.im+
      u[ 5]*(*s).c2.c3.re-   mu*(*s).c2.c3.im;

   rs.c2.c3.im=
      u[14]*(*s).c1.c1.im-u[15]*(*s).c1.c1.re+
      u[22]*(*s).c1.c2.im-u[23]*(*s).c1.c2.re+
      u[28]*(*s).c1.c3.im-u[29]*(*s).c1.c3.re+
      u[32]*(*s).c2.c1.im-u[33]*(*s).c2.c1.re+
      u[34]*(*s).c2.c2.im-u[35]*(*s).c2.c2.re+
      u[ 5]*(*s).c2.c3.im+   mu*(*s).c2.c3.re;

   (*r)=rs;
}


void mul_pauli2(float mu,pauli *m,spinor *s,spinor *r)
{
   spin_t *ps,*pr;

   ps=(spin_t*)(s);
   pr=(spin_t*)(r);

   mul_pauli(mu,m,(*ps).w,(*pr).w);
   mul_pauli(-mu,m+1,(*ps).w+1,(*pr).w+1);
}

#endif

void assign_pauli(int vol,pauli_dble *md,pauli *m)
{
   float *u;
   double *ud,*um;
   pauli_dble *mm;

   mm=md+vol;

   for (;md<mm;md++)
   {
      u=(*m).u;
      ud=(*md).u;
      um=ud+36;

      for (;ud<um;ud+=9)
      {
         u[0]=(float)(ud[0]);
         u[1]=(float)(ud[1]);
         u[2]=(float)(ud[2]);
         u[3]=(float)(ud[3]);
         u[4]=(float)(ud[4]);
         u[5]=(float)(ud[5]);
         u[6]=(float)(ud[6]);
         u[7]=(float)(ud[7]);
         u[8]=(float)(ud[8]);

         u+=9;
      }

      m+=1;
   }
}


void apply_sw(int vol,float mu,pauli *m,spinor *s,spinor *r)
{
   spinor *sm;

   sm=s+vol;

   for (;s<sm;s++)
   {
      mul_pauli2(mu,m,s,r);
      m+=2;
      r+=1;
   }
}
