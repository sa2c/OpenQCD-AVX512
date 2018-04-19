
/*******************************************************************************
*
* File cmatrix.c
*
* Copyright (C) 2007, 2009, 2011, 2013,  Martin Luescher, Isabel Campos
*               2016
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Complex matrix algebra (single-precision version).
*
* The externally accessible functions are
*
*   void cmat_vec(int n,complex *a,complex *v,complex *w)
*     Computes w=a*v, where v and w are n-vectors and a an nxn matrix.
*
*   void cmat_vec_assign(int n,complex *a,complex *v,complex *w)
*     Adds a*v to w, where v and w are n-vectors and a an nxn matrix.
*
*   void cmat_add(int n,complex *a,complex *b,complex *c)
*     Computes the sum c=a+b of two nxn matrices a and b.
*
*   void cmat_sub(int n,complex *a,complex *b,complex *c)
*     Computes the difference c=a-b of two nxn matrices a and b.
*
*   void cmat_mul(int n,complex *a,complex *b,complex *c)
*     Computes the product c=a*b of two nxn matrices a and b.
*
*   void cmat_dag(int n,complex *a,complex *b)
*     Assigns the hermitian conjugate of a to b.
*
* Notes:
*
* All of these programs can be called locally. Complex nxn matrices with
* matrix elements A_{ij} are represented by linear arrays a of complex
* numbers such that
*
*   A_{ij} = a[i*n+j]
*
* where i,j=0,1,..,n-1. It is assumed that the input and output arrays do
* not overlap in memory (the results are otherwise unpredictable).
*
* If SSE or AVX instructions are to be used, and if n is even, it is taken
* for granted that the arrays are aligned to a 16 byte boundary.
*
*******************************************************************************/

#define CMATRIX_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "utils.h"
#include "linalg.h"

#if (defined AVX)
#include "avx.h"

#if (defined FMA3)

void cmat_vec(int n,complex *a,complex *v,complex *w)
{
   complex *b[4],*vv,*vm,*vmm,*wm;

   wm=w+n;

   if ((n&0x3)==0x0)
   {
      vm=v+n;
      b[3]=a;

      for (;w<wm;w+=4)
      {
         b[0]=b[3];
         b[1]=b[0]+n;
         b[2]=b[1]+n;
         b[3]=b[2]+n;

         __asm__ __volatile__ ("vxorps %%ymm0, %%ymm0, %%ymm0 \n\t"
                               "vxorps %%ymm1, %%ymm1, %%ymm1 \n\t"
                               "vxorps %%ymm2, %%ymm2, %%ymm2 \n\t"
                               "vxorps %%ymm3, %%ymm3, %%ymm3"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovshdup %0, %%ymm4 \n\t"
                                  "vmovsldup %0, %%ymm5"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3])
                                  :
                                  "xmm4", "xmm5");

            __asm__ __volatile__ ("vmovups %0, %%ymm6 \n\t"
                                  "vmovups %4, %%ymm7"
                                  :
                                  :
                                  "m" (b[0][0]),
                                  "m" (b[0][1]),
                                  "m" (b[0][2]),
                                  "m" (b[0][3]),
                                  "m" (b[1][0]),
                                  "m" (b[1][1]),
                                  "m" (b[1][2]),
                                  "m" (b[1][3])
                                  :
                                  "xmm6", "xmm7");

            __asm__ __volatile__ ("vmovups %0, %%ymm8 \n\t"
                                  "vmovups %4, %%ymm9"
                                  :
                                  :
                                  "m" (b[2][0]),
                                  "m" (b[2][1]),
                                  "m" (b[2][2]),
                                  "m" (b[2][3]),
                                  "m" (b[3][0]),
                                  "m" (b[3][1]),
                                  "m" (b[3][2]),
                                  "m" (b[3][3])
                                  :
                                  "xmm8", "xmm9");

            __asm__ __volatile__ ("vpermilps $0xb1, %%ymm6, %%ymm10 \n\t"
                                  "vpermilps $0xb1, %%ymm7, %%ymm11 \n\t"
                                  "vpermilps $0xb1, %%ymm8, %%ymm12 \n\t"
                                  "vpermilps $0xb1, %%ymm9, %%ymm13 \n\t"
                                  "vfmaddsub231ps %%ymm4, %%ymm10, %%ymm0 \n\t"
                                  "vfmaddsub231ps %%ymm4, %%ymm11, %%ymm1 \n\t"
                                  "vfmaddsub231ps %%ymm4, %%ymm12, %%ymm2 \n\t"
                                  "vfmaddsub231ps %%ymm4, %%ymm13, %%ymm3 \n\t"
                                  "vfmaddsub231ps %%ymm5, %%ymm6, %%ymm0 \n\t"
                                  "vfmaddsub231ps %%ymm5, %%ymm7, %%ymm1 \n\t"
                                  "vfmaddsub231ps %%ymm5, %%ymm8, %%ymm2 \n\t"
                                  "vfmaddsub231ps %%ymm5, %%ymm9, %%ymm3"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3",
                                  "xmm10", "xmm11", "xmm12", "xmm13");

            b[0]+=4;
            b[1]+=4;
            b[2]+=4;
            b[3]+=4;
         }

         __asm__ __volatile__ ("vextractf128 $0x1, %%ymm0, %%xmm10 \n\t"
                               "vextractf128 $0x1, %%ymm1, %%xmm11 \n\t"
                               "vextractf128 $0x1, %%ymm2, %%xmm12 \n\t"
                               "vextractf128 $0x1, %%ymm3, %%xmm13 \n\t"
                               "vaddps %%xmm10, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm11, %%xmm1, %%xmm1 \n\t"
                               "vaddps %%xmm12, %%xmm2, %%xmm2 \n\t"
                               "vaddps %%xmm13, %%xmm3, %%xmm3 \n\t"
                               "vinsertf128 $0x1, %%xmm2, %%ymm0, %%ymm0 \n\t"
                               "vinsertf128 $0x1, %%xmm3, %%ymm1, %%ymm1 \n\t"
                               "vpermilps $0xd8, %%ymm0, %%ymm0 \n\t"
                               "vpermilps $0xd8, %%ymm1, %%ymm1 \n\t"
                               "vhaddps %%ymm1, %%ymm0, %%ymm0 \n\t"
                               "vmovups %%ymm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1]),
                               "=m" (w[2]),
                               "=m" (w[3])
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm10", "xmm11", "xmm12", "xmm13");
      }
   }
   else if ((n&0x1)==0x0)
   {
      vm=v+n-2;
      b[1]=a;

      for (;w<wm;w+=2)
      {
         b[0]=b[1];
         b[1]=b[0]+n;

         __asm__ __volatile__ ("vxorps %%ymm0, %%ymm0, %%ymm0 \n\t"
                               "vxorps %%ymm1, %%ymm1, %%ymm1 \n\t"
                               "vxorps %%ymm2, %%ymm2, %%ymm2 \n\t"
                               "vxorps %%ymm3, %%ymm3, %%ymm3"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovsldup %0, %%ymm4 \n\t"
                                  "vmovshdup %0, %%ymm5"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3])
                                  :
                                  "xmm4", "xmm5");

            __asm__ __volatile__ ("vmovups %0, %%ymm6 \n\t"
                                  "vmovups %4, %%ymm7"
                                  :
                                  :
                                  "m" (b[0][0]),
                                  "m" (b[0][1]),
                                  "m" (b[0][2]),
                                  "m" (b[0][3]),
                                  "m" (b[1][0]),
                                  "m" (b[1][1]),
                                  "m" (b[1][2]),
                                  "m" (b[1][3])
                                  :
                                  "xmm6", "xmm7");

            __asm__ __volatile__ ("vfmadd231ps %%ymm4, %%ymm6, %%ymm0 \n\t"
                                  "vfmadd231ps %%ymm4, %%ymm7, %%ymm1 \n\t"
                                  "vfmadd231ps %%ymm5, %%ymm6, %%ymm2 \n\t"
                                  "vfmadd231ps %%ymm5, %%ymm7, %%ymm3"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3");

            b[0]+=4;
            b[1]+=4;
         }

         __asm__ __volatile__ ("vextractf128 $0x1, %%ymm0, %%xmm8 \n\t"
                               "vextractf128 $0x1, %%ymm1, %%xmm9 \n\t"
                               "vextractf128 $0x1, %%ymm2, %%xmm10 \n\t"
                               "vextractf128 $0x1, %%ymm3, %%xmm11\n\t"
                               "vaddps %%xmm8, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm9, %%xmm1, %%xmm1 \n\t"
                               "vaddps %%xmm10, %%xmm2, %%xmm2 \n\t"
                               "vaddps %%xmm11, %%xmm3, %%xmm3"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm8", "xmm9", "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovsldup %0, %%xmm4 \n\t"
                               "vmovshdup %0, %%xmm5 \n\t"
                               "vmovaps %2, %%xmm6 \n\t"
                               "vmovaps %4, %%xmm7 \n\t"
                               "vfmadd231ps %%xmm4, %%xmm6, %%xmm0 \n\t"
                               "vfmadd231ps %%xmm4, %%xmm7, %%xmm1 \n\t"
                               "vfmadd231ps %%xmm5, %%xmm6, %%xmm2 \n\t"
                               "vfmadd231ps %%xmm5, %%xmm7, %%xmm3"
                               :
                               :
                               "m" (vv[0]),
                               "m" (vv[1]),
                               "m" (b[0][0]),
                               "m" (b[0][1]),
                               "m" (b[1][0]),
                               "m" (b[1][1])
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm4", "xmm5", "xmm6", "xmm7");

         __asm__ __volatile__ ("vpermilps $0xb1, %%xmm2, %%xmm2 \n\t"
                               "vpermilps $0xb1, %%xmm3, %%xmm3 \n\t"
                               "vaddsubps %%xmm2, %%xmm0, %%xmm0 \n\t"
                               "vaddsubps %%xmm3, %%xmm1, %%xmm1"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         __asm__ __volatile__ ("vpermilps $0xd8, %%xmm0, %%xmm0 \n\t"
                               "vpermilps $0xd8, %%xmm1, %%xmm1 \n\t"
                               "vhaddps %%xmm1, %%xmm0, %%xmm0 \n\t"
                               "vmovaps %%xmm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1])
                               :
                               :
                               "xmm0", "xmm1");

         b[1]+=2;
      }
   }
   else
   {
      vmm=v+n;
      vm=vmm-(n&0x3);
      b[0]=a;

      for (;w<wm;w++)
      {
         __asm__ __volatile__ ("vxorps %%ymm0, %%ymm0, %%ymm0 \n\t"
                               "vxorps %%ymm1, %%ymm1, %%ymm1"
                               :
                               :
                               :
                               "xmm0", "xmm1");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovsldup %0, %%ymm4 \n\t"
                                  "vmovshdup %0, %%ymm5"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3])
                                  :
                                  "xmm4", "xmm5");

            __asm__ __volatile__ ("vmovups %0, %%ymm6 \n\t"
                                  "vfmadd231ps %%ymm4, %%ymm6, %%ymm0 \n\t"
                                  "vfmadd231ps %%ymm5, %%ymm6, %%ymm1"
                                  :
                                  :
                                  "m" (b[0][0]),
                                  "m" (b[0][1]),
                                  "m" (b[0][2]),
                                  "m" (b[0][3])
                                  :
                                  "xmm0", "xmm1", "xmm6");

            b[0]+=4;
         }

         __asm__ __volatile__ ("vextractf128 $0x1, %%ymm0, %%xmm7 \n\t"
                               "vextractf128 $0x1, %%ymm1, %%xmm8 \n\t"
                               "vaddps %%xmm7, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm8, %%xmm1, %%xmm1"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm7", "xmm8");

         for (;vv<vmm;vv++)
         {
            __asm__ __volatile__ ("vmovsd %0, %%xmm4 \n\t"
                                  "vmovsd %1, %%xmm6 \n\t"
                                  "vpermilps $0xe5, %%xmm4, %%xmm5 \n\t"
                                  "vpermilps $0xe0, %%xmm4, %%xmm4 \n\t"
                                  "vfmadd231ps %%xmm5, %%xmm6, %%xmm1 \n\t"
                                  "vfmadd231ps %%xmm4, %%xmm6, %%xmm0"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (b[0][0])
                                  :
                                  "xmm0", "xmm1", "xmm4", "xmm5",
                                  "xmm6");

            b[0]+=1;
         }

         __asm__ __volatile__ ("vpermilps $0xb1, %%xmm1, %%xmm1 \n\t"
                               "vaddsubps %%xmm1, %%xmm0, %%xmm2 \n\t"
                               "vpermilps $0x4e, %%xmm2, %%xmm3 \n\t"
                               "vaddps %%xmm2, %%xmm3, %%xmm0 \n\t"
                               "vmovsd %%xmm0, %0"
                               :
                               "=m" (w[0])
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");
      }
   }

   _avx_zeroupper();
}


void cmat_vec_assign(int n,complex *a,complex *v,complex *w)
{
   complex *b[4],*vv,*vm,*vmm,*wm;

   wm=w+n;

   if ((n&0x3)==0x0)
   {
      vm=v+n;
      b[3]=a;

      for (;w<wm;w+=4)
      {
         b[0]=b[3];
         b[1]=b[0]+n;
         b[2]=b[1]+n;
         b[3]=b[2]+n;

         __asm__ __volatile__ ("vmovsd %0, %%xmm0 \n\t"
                               "vmovsd %1, %%xmm1 \n\t"
                               "vmovsd %2, %%xmm2 \n\t"
                               "vmovsd %3, %%xmm3"
                               :
                               :
                               "m" (w[0]),
                               "m" (w[1]),
                               "m" (w[2]),
                               "m" (w[3])
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovshdup %0, %%ymm4 \n\t"
                                  "vmovsldup %0, %%ymm5"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3])
                                  :
                                  "xmm4", "xmm5");

            __asm__ __volatile__ ("vmovups %0, %%ymm6 \n\t"
                                  "vmovups %4, %%ymm7"
                                  :
                                  :
                                  "m" (b[0][0]),
                                  "m" (b[0][1]),
                                  "m" (b[0][2]),
                                  "m" (b[0][3]),
                                  "m" (b[1][0]),
                                  "m" (b[1][1]),
                                  "m" (b[1][2]),
                                  "m" (b[1][3])
                                  :
                                  "xmm6", "xmm7");

            __asm__ __volatile__ ("vmovups %0, %%ymm8 \n\t"
                                  "vmovups %4, %%ymm9"
                                  :
                                  :
                                  "m" (b[2][0]),
                                  "m" (b[2][1]),
                                  "m" (b[2][2]),
                                  "m" (b[2][3]),
                                  "m" (b[3][0]),
                                  "m" (b[3][1]),
                                  "m" (b[3][2]),
                                  "m" (b[3][3])
                                  :
                                  "xmm8", "xmm9");

            __asm__ __volatile__ ("vpermilps $0xb1, %%ymm6, %%ymm10 \n\t"
                                  "vpermilps $0xb1, %%ymm7, %%ymm11 \n\t"
                                  "vpermilps $0xb1, %%ymm8, %%ymm12 \n\t"
                                  "vpermilps $0xb1, %%ymm9, %%ymm13 \n\t"
                                  "vfmaddsub231ps %%ymm4, %%ymm10, %%ymm0 \n\t"
                                  "vfmaddsub231ps %%ymm4, %%ymm11, %%ymm1 \n\t"
                                  "vfmaddsub231ps %%ymm4, %%ymm12, %%ymm2 \n\t"
                                  "vfmaddsub231ps %%ymm4, %%ymm13, %%ymm3 \n\t"
                                  "vfmaddsub231ps %%ymm5, %%ymm6, %%ymm0 \n\t"
                                  "vfmaddsub231ps %%ymm5, %%ymm7, %%ymm1 \n\t"
                                  "vfmaddsub231ps %%ymm5, %%ymm8, %%ymm2 \n\t"
                                  "vfmaddsub231ps %%ymm5, %%ymm9, %%ymm3"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3",
                                  "xmm10", "xmm11", "xmm12", "xmm13");

            b[0]+=4;
            b[1]+=4;
            b[2]+=4;
            b[3]+=4;
         }

         __asm__ __volatile__ ("vextractf128 $0x1, %%ymm0, %%xmm10 \n\t"
                               "vextractf128 $0x1, %%ymm1, %%xmm11 \n\t"
                               "vextractf128 $0x1, %%ymm2, %%xmm12 \n\t"
                               "vextractf128 $0x1, %%ymm3, %%xmm13 \n\t"
                               "vaddps %%xmm10, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm11, %%xmm1, %%xmm1 \n\t"
                               "vaddps %%xmm12, %%xmm2, %%xmm2 \n\t"
                               "vaddps %%xmm13, %%xmm3, %%xmm3 \n\t"
                               "vinsertf128 $0x1, %%xmm2, %%ymm0, %%ymm0 \n\t"
                               "vinsertf128 $0x1, %%xmm3, %%ymm1, %%ymm1 \n\t"
                               "vpermilps $0xd8, %%ymm0, %%ymm0 \n\t"
                               "vpermilps $0xd8, %%ymm1, %%ymm1 \n\t"
                               "vhaddps %%ymm1, %%ymm0, %%ymm0 \n\t"
                               "vmovups %%ymm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1]),
                               "=m" (w[2]),
                               "=m" (w[3])
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm10", "xmm11", "xmm12", "xmm13");
      }
   }
   else if ((n&0x1)==0x0)
   {
      vm=v+n-2;
      b[1]=a;

      for (;w<wm;w+=2)
      {
         b[0]=b[1];
         b[1]=b[0]+n;

         __asm__ __volatile__ ("vmovsd %0, %%xmm0  \n\t"
                               "vmovsd %1, %%xmm1 \n\t"
                               "vxorps %%ymm2, %%ymm2, %%ymm2 \n\t"
                               "vxorps %%ymm3, %%ymm3, %%ymm3"
                               :
                               :
                               "m" (w[0]),
                               "m" (w[1])
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovsldup %0, %%ymm4 \n\t"
                                  "vmovshdup %0, %%ymm5"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3])
                                  :
                                  "xmm4", "xmm5");

            __asm__ __volatile__ ("vmovups %0, %%ymm6 \n\t"
                                  "vmovups %4, %%ymm7"
                                  :
                                  :
                                  "m" (b[0][0]),
                                  "m" (b[0][1]),
                                  "m" (b[0][2]),
                                  "m" (b[0][3]),
                                  "m" (b[1][0]),
                                  "m" (b[1][1]),
                                  "m" (b[1][2]),
                                  "m" (b[1][3])
                                  :
                                  "xmm6", "xmm7");

            __asm__ __volatile__ ("vfmadd231ps %%ymm4, %%ymm6, %%ymm0 \n\t"
                                  "vfmadd231ps %%ymm4, %%ymm7, %%ymm1 \n\t"
                                  "vfmadd231ps %%ymm5, %%ymm6, %%ymm2 \n\t"
                                  "vfmadd231ps %%ymm5, %%ymm7, %%ymm3"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3");

            b[0]+=4;
            b[1]+=4;
         }

         __asm__ __volatile__ ("vextractf128 $0x1, %%ymm0, %%xmm8 \n\t"
                               "vextractf128 $0x1, %%ymm1, %%xmm9 \n\t"
                               "vextractf128 $0x1, %%ymm2, %%xmm10 \n\t"
                               "vextractf128 $0x1, %%ymm3, %%xmm11\n\t"
                               "vaddps %%xmm8, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm9, %%xmm1, %%xmm1 \n\t"
                               "vaddps %%xmm10, %%xmm2, %%xmm2 \n\t"
                               "vaddps %%xmm11, %%xmm3, %%xmm3"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm8", "xmm9", "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovsldup %0, %%xmm4 \n\t"
                               "vmovshdup %0, %%xmm5 \n\t"
                               "vmovaps %2, %%xmm6 \n\t"
                               "vmovaps %4, %%xmm7 \n\t"
                               "vfmadd231ps %%xmm4, %%xmm6, %%xmm0 \n\t"
                               "vfmadd231ps %%xmm4, %%xmm7, %%xmm1 \n\t"
                               "vfmadd231ps %%xmm5, %%xmm6, %%xmm2 \n\t"
                               "vfmadd231ps %%xmm5, %%xmm7, %%xmm3"
                               :
                               :
                               "m" (vv[0]),
                               "m" (vv[1]),
                               "m" (b[0][0]),
                               "m" (b[0][1]),
                               "m" (b[1][0]),
                               "m" (b[1][1])
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm4", "xmm5", "xmm6", "xmm7");

         __asm__ __volatile__ ("vpermilps $0xb1, %%xmm2, %%xmm2 \n\t"
                               "vpermilps $0xb1, %%xmm3, %%xmm3 \n\t"
                               "vaddsubps %%xmm2, %%xmm0, %%xmm0 \n\t"
                               "vaddsubps %%xmm3, %%xmm1, %%xmm1"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         __asm__ __volatile__ ("vpermilps $0xd8, %%xmm0, %%xmm0 \n\t"
                               "vpermilps $0xd8, %%xmm1, %%xmm1 \n\t"
                               "vhaddps %%xmm1, %%xmm0, %%xmm0 \n\t"
                               "vmovaps %%xmm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1])
                               :
                               :
                               "xmm0", "xmm1");

         b[1]+=2;
      }
   }
   else
   {
      vmm=v+n;
      vm=vmm-(n&0x3);
      b[0]=a;

      for (;w<wm;w+=1)
      {
         __asm__ __volatile__ ("vmovsd %0, %%xmm0 \n\t"
                               "vxorps %%ymm1, %%ymm1, %%ymm1"
                               :
                               :
                               "m" (w[0])
                               :
                               "xmm0", "xmm1");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovsldup %0, %%ymm4 \n\t"
                                  "vmovshdup %0, %%ymm5"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3])
                                  :
                                  "xmm4", "xmm5");

            __asm__ __volatile__ ("vmovups %0, %%ymm6 \n\t"
                                  "vfmadd231ps %%ymm4, %%ymm6, %%ymm0 \n\t"
                                  "vfmadd231ps %%ymm5, %%ymm6, %%ymm1"
                                  :
                                  :
                                  "m" (b[0][0]),
                                  "m" (b[0][1]),
                                  "m" (b[0][2]),
                                  "m" (b[0][3])
                                  :
                                  "xmm0", "xmm1", "xmm6");

            b[0]+=4;
         }

         __asm__ __volatile__ ("vextractf128 $0x1, %%ymm0, %%xmm7 \n\t"
                               "vextractf128 $0x1, %%ymm1, %%xmm8 \n\t"
                               "vaddps %%xmm7, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm8, %%xmm1, %%xmm1"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm7", "xmm8");

         for (;vv<vmm;vv++)
         {
            __asm__ __volatile__ ("vmovsd %0, %%xmm4 \n\t"
                                  "vmovsd %1, %%xmm6 \n\t"
                                  "vpermilps $0xe5, %%xmm4, %%xmm5 \n\t"
                                  "vpermilps $0xe0, %%xmm4, %%xmm4 \n\t"
                                  "vfmadd231ps %%xmm5, %%xmm6, %%xmm1 \n\t"
                                  "vfmadd231ps %%xmm4, %%xmm6, %%xmm0"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (b[0][0])
                                  :
                                  "xmm0", "xmm1", "xmm4", "xmm5",
                                  "xmm6");

            b[0]+=1;
         }

         __asm__ __volatile__ ("vpermilps $0xb1, %%xmm1, %%xmm1 \n\t"
                               "vaddsubps %%xmm1, %%xmm0, %%xmm2 \n\t"
                               "vpermilps $0x4e, %%xmm2, %%xmm3 \n\t"
                               "vaddps %%xmm2, %%xmm3, %%xmm0 \n\t"
                               "vmovsd %%xmm0, %0"
                               :
                               "=m" (w[0])
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");
      }
   }

   _avx_zeroupper();
}

#else

void cmat_vec(int n,complex *a,complex *v,complex *w)
{
   complex *b,*vv,*vm,*wm;

   if ((n&0x3)==0x0)
   {
      vm=v+n;
      wm=w+n;
      b=a;

      for (;w<wm;w+=2)
      {
         a=b;
         b=a+n;

         __asm__ __volatile__ ("vxorps %%ymm0, %%ymm0, %%ymm0 \n\t"
                               "vxorps %%ymm1, %%ymm1, %%ymm1 \n\t"
                               "vxorps %%ymm2, %%ymm2, %%ymm2 \n\t"
                               "vxorps %%ymm3, %%ymm3, %%ymm3"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovaps %0, %%xmm8 \n\t"
                                  "vmovsldup %2, %%ymm4 \n\t"
                                  "vmovshdup %2, %%ymm5"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm8");

            __asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm8, %%ymm8 \n\t"
                                  "vmovsldup %2, %%ymm6 \n\t"
                                  "vmovshdup %2, %%ymm7 \n\t"
                                  "vpermilps $0xb1, %%ymm8, %%ymm9"
                                  :
                                  :
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (b[0]),
                                  "m" (b[1]),
                                  "m" (b[2]),
                                  "m" (b[3])
                                  :
                                  "xmm6", "xmm7", "xmm8", "xmm9");

            __asm__ __volatile__ ("vmulps %%ymm8, %%ymm4, %%ymm4 \n\t"
                                  "vmulps %%ymm9, %%ymm5, %%ymm5 \n\t"
                                  "vmulps %%ymm8, %%ymm6, %%ymm6 \n\t"
                                  "vmulps %%ymm9, %%ymm7, %%ymm7 \n\t"
                                  "vaddps %%ymm4, %%ymm0, %%ymm0 \n\t"
                                  "vaddps %%ymm5, %%ymm1, %%ymm1 \n\t"
                                  "vaddps %%ymm6, %%ymm2, %%ymm2 \n\t"
                                  "vaddps %%ymm7, %%ymm3, %%ymm3"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3",
                                  "xmm4", "xmm5", "xmm6", "xmm7");

            a+=4;
            b+=4;
         }

         __asm__ __volatile__ ("vaddsubps %%ymm1, %%ymm0, %%ymm0 \n\t"
                               "vaddsubps %%ymm3, %%ymm2, %%ymm2 \n\t"
                               "vextractf128 $0x1, %%ymm0, %%xmm8 \n\t"
                               "vextractf128 $0x1, %%ymm2, %%xmm9 \n\t"
                               "vaddps %%xmm8, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm9, %%xmm2, %%xmm2 \n\t"
                               "vpermilps $0xd8, %%xmm0, %%xmm0 \n\t"
                               "vpermilps $0xd8, %%xmm2, %%xmm2 \n\t"
                               "vhaddps %%xmm2, %%xmm0, %%xmm0 \n\t"
                               "vmovaps %%xmm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1])
                               :
                               :
                               "xmm0", "xmm2", "xmm8", "xmm9");
      }
   }
   else if ((n&0x1)==0x0)
   {
      vm=v+n-2;
      wm=w+n;
      b=a;

      for (;w<wm;w+=2)
      {
         a=b;
         b=a+n;

         __asm__ __volatile__ ("vxorps %%ymm0, %%ymm0, %%ymm0 \n\t"
                               "vxorps %%ymm1, %%ymm1, %%ymm1 \n\t"
                               "vxorps %%ymm2, %%ymm2, %%ymm2 \n\t"
                               "vxorps %%ymm3, %%ymm3, %%ymm3"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovaps %0, %%xmm8 \n\t"
                                  "vmovsldup %2, %%ymm4 \n\t"
                                  "vmovshdup %2, %%ymm5"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm8");

            __asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm8, %%ymm8 \n\t"
                                  "vmovsldup %2, %%ymm6 \n\t"
                                  "vmovshdup %2, %%ymm7 \n\t"
                                  "vpermilps $0xb1, %%ymm8, %%ymm9"
                                  :
                                  :
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (b[0]),
                                  "m" (b[1]),
                                  "m" (b[2]),
                                  "m" (b[3])
                                  :
                                  "xmm6", "xmm7", "xmm8", "xmm9");

            __asm__ __volatile__ ("vmulps %%ymm8, %%ymm4, %%ymm4 \n\t"
                                  "vmulps %%ymm9, %%ymm5, %%ymm5 \n\t"
                                  "vmulps %%ymm8, %%ymm6, %%ymm6 \n\t"
                                  "vmulps %%ymm9, %%ymm7, %%ymm7 \n\t"
                                  "vaddps %%ymm4, %%ymm0, %%ymm0 \n\t"
                                  "vaddps %%ymm5, %%ymm1, %%ymm1 \n\t"
                                  "vaddps %%ymm6, %%ymm2, %%ymm2 \n\t"
                                  "vaddps %%ymm7, %%ymm3, %%ymm3"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3",
                                  "xmm4", "xmm5", "xmm6", "xmm7");

            a+=4;
            b+=4;
         }

         __asm__ __volatile__ ("vaddsubps %%ymm1, %%ymm0, %%ymm0 \n\t"
                               "vaddsubps %%ymm3, %%ymm2, %%ymm2 \n\t"
                               "vextractf128 $0x1, %%ymm0, %%xmm8 \n\t"
                               "vextractf128 $0x1, %%ymm2, %%xmm9 \n\t"
                               "vaddps %%xmm8, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm9, %%xmm2, %%xmm2"
                               :
                               :
                               :
                               "xmm0", "xmm2", "xmm8", "xmm9");

         __asm__ __volatile__ ("vmovaps %0, %%xmm4 \n\t"
                               "vmovsldup %2, %%xmm6 \n\t"
                               "vmovshdup %2, %%xmm7 \n\t"
                               "vpermilps $0xb1, %%xmm4, %%xmm5 \n\t"
                               "vmovsldup %4, %%xmm8 \n\t"
                               "vmovshdup %4, %%xmm9 \n\t"
                               "vmulps %%xmm4, %%xmm6, %%xmm6 \n\t"
                               "vmulps %%xmm5, %%xmm7, %%xmm7 \n\t"
                               "vmulps %%xmm4, %%xmm8, %%xmm8 \n\t"
                               "vmulps %%xmm5, %%xmm9, %%xmm9 \n\t"
                               "vaddsubps %%xmm7, %%xmm6, %%xmm6 \n\t"
                               "vaddsubps %%xmm9, %%xmm8, %%xmm8 \n\t"
                               "vaddps %%xmm6, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm8, %%xmm2, %%xmm2"
                               :
                               :
                               "m" (vv[0]),
                               "m" (vv[1]),
                               "m" (a[0]),
                               "m" (a[1]),
                               "m" (b[0]),
                               "m" (b[1])
                               :
                               "xmm0", "xmm2", "xmm4", "xmm5",
                               "xmm6", "xmm7", "xmm8", "xmm9");

         a+=2;
         b+=2;

         __asm__ __volatile__ ("vpermilps $0xd8, %%xmm0, %%xmm0 \n\t"
                               "vpermilps $0xd8, %%xmm2, %%xmm2 \n\t"
                               "vhaddps %%xmm2, %%xmm0, %%xmm0 \n\t"
                               "vmovaps %%xmm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1])
                               :
                               :
                               "xmm0", "xmm2");
      }
   }
   else
   {
      vm=v+n-(n&0x3);
      wm=w+n;

      __asm__ __volatile__ ("vxorps %%ymm9, %%ymm9, %%ymm9"
                            :
                            :
                            :
                            "xmm9");

      for (;w<wm;w+=1)
      {
         __asm__ __volatile__ ("vxorps %%ymm0, %%ymm0, %%ymm0 \n\t"
                               "vxorps %%ymm1, %%ymm1, %%ymm1"
                               :
                               :
                               :
                               "xmm0", "xmm1");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovups %0, %%ymm6 \n\t"
                                  "vmovsldup %4, %%ymm4 \n\t"
                                  "vmovshdup %4, %%ymm5 \n\t"
                                  "vpermilps $0xb1, %%ymm6, %%ymm7"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7");

            __asm__ __volatile__ ("vmulps %%ymm6, %%ymm4, %%ymm4 \n\t"
                                  "vmulps %%ymm7, %%ymm5, %%ymm5 \n\t"
                                  "vaddps %%ymm4, %%ymm0, %%ymm0 \n\t"
                                  "vaddps %%ymm5, %%ymm1, %%ymm1"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm4", "xmm5");

            a+=4;
         }

         __asm__ __volatile__ ("vaddsubps %%ymm1, %%ymm0, %%ymm2 \n\t"
                               "vextractf128 $0x1, %%ymm2, %%xmm3"
                               :
                               :
                               :
                               "xmm2", "xmm3");

         if ((n&0x2)!=0x0)
         {
            __asm__ __volatile__ ("vmovups %0, %%xmm4 \n\t"
                                  "vmovsldup %2, %%xmm6 \n\t"
                                  "vmovshdup %2, %%xmm7 \n\t"
                                  "vpermilps $0xb1, %%xmm4, %%xmm5 \n\t"
                                  "vmulps %%xmm4, %%xmm6, %%xmm6 \n\t"
                                  "vmulps %%xmm5, %%xmm7, %%xmm7 \n\t"
                                  "vaddps %%xmm6, %%xmm2, %%xmm2 \n\t"
                                  "vaddsubps %%xmm7, %%xmm3, %%xmm3"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (a[0]),
                                  "m" (a[1])
                                  :
                                  "xmm2", "xmm3", "xmm4", "xmm5",
                                  "xmm6", "xmm7");

            a+=2;
            vv+=2;
         }

         __asm__ __volatile__ ("vmovlps %0, %%xmm9, %%xmm4 \n\t"
                               "vmovlps %1, %%xmm9, %%xmm8 \n\t"
                               "vmovsldup %%xmm8, %%xmm6 \n\t"
                               "vmovshdup %%xmm8, %%xmm7 \n\t"
                               "vpermilps $0xb1, %%xmm4, %%xmm5 \n\t"
                               "vmulps %%xmm4, %%xmm6, %%xmm6 \n\t"
                               "vmulps %%xmm5, %%xmm7, %%xmm7 \n\t"
                               "vaddps %%xmm6, %%xmm2, %%xmm2 \n\t"
                               "vaddsubps %%xmm7, %%xmm3, %%xmm3"
                               :
                               :
                               "m" (vv[0]),
                               "m" (a[0])
                               :
                               "xmm2", "xmm3", "xmm4", "xmm5",
                               "xmm6", "xmm7", "xmm8");

         a+=1;

         __asm__ __volatile__ ("vaddps %%xmm3, %%xmm2, %%xmm2 \n\t"
                               "vpermilps $0xd8, %%xmm2, %%xmm0 \n\t"
                               "vhaddps %%xmm0, %%xmm0, %%xmm0 \n\t"
                               "vmovlps %%xmm0, %0"
                               :
                               "=m" (w[0])
                               :
                               :
                               "xmm0", "xmm2");
      }
   }

   _avx_zeroupper();
}


void cmat_vec_assign(int n,complex *a,complex *v,complex *w)
{
   complex *b,*vv,*vm,*wm;

   if ((n&0x3)==0x0)
   {
      vm=v+n;
      wm=w+n;
      b=a;

      for (;w<wm;w+=2)
      {
         a=b;
         b=a+n;

         __asm__ __volatile__ ("vxorps %%ymm1, %%ymm1, %%ymm1 \n\t"
                               "vxorps %%ymm3, %%ymm3, %%ymm3 \n\t"
                               "vmovlps %0, %%xmm1, %%xmm0 \n\t"
                               "vmovlps %1, %%xmm3, %%xmm2"
                               :
                               :
                               "m" (w[0]),
                               "m" (w[1])
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovaps %0, %%xmm8 \n\t"
                                  "vmovsldup %2, %%ymm4 \n\t"
                                  "vmovshdup %2, %%ymm5"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm8");

            __asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm8, %%ymm8 \n\t"
                                  "vmovsldup %2, %%ymm6 \n\t"
                                  "vmovshdup %2, %%ymm7 \n\t"
                                  "vpermilps $0xb1, %%ymm8, %%ymm9"
                                  :
                                  :
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (b[0]),
                                  "m" (b[1]),
                                  "m" (b[2]),
                                  "m" (b[3])
                                  :
                                  "xmm6", "xmm7", "xmm8", "xmm9");

            __asm__ __volatile__ ("vmulps %%ymm8, %%ymm4, %%ymm4 \n\t"
                                  "vmulps %%ymm9, %%ymm5, %%ymm5 \n\t"
                                  "vmulps %%ymm8, %%ymm6, %%ymm6 \n\t"
                                  "vmulps %%ymm9, %%ymm7, %%ymm7 \n\t"
                                  "vaddps %%ymm4, %%ymm0, %%ymm0 \n\t"
                                  "vaddps %%ymm5, %%ymm1, %%ymm1 \n\t"
                                  "vaddps %%ymm6, %%ymm2, %%ymm2 \n\t"
                                  "vaddps %%ymm7, %%ymm3, %%ymm3"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3",
                                  "xmm4", "xmm5", "xmm6", "xmm7");

            a+=4;
            b+=4;
         }

         __asm__ __volatile__ ("vaddsubps %%ymm1, %%ymm0, %%ymm0 \n\t"
                               "vaddsubps %%ymm3, %%ymm2, %%ymm2 \n\t"
                               "vextractf128 $0x1, %%ymm0, %%xmm8 \n\t"
                               "vextractf128 $0x1, %%ymm2, %%xmm9 \n\t"
                               "vaddps %%xmm8, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm9, %%xmm2, %%xmm2 \n\t"
                               "vpermilps $0xd8, %%xmm0, %%xmm0 \n\t"
                               "vpermilps $0xd8, %%xmm2, %%xmm2 \n\t"
                               "vhaddps %%xmm2, %%xmm0, %%xmm0 \n\t"
                               "vmovaps %%xmm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1])
                               :
                               :
                               "xmm0", "xmm2", "xmm8", "xmm9");
      }
   }
   else if ((n&0x1)==0x0)
   {
      vm=v+n-2;
      wm=w+n;
      b=a;

      for (;w<wm;w+=2)
      {
         a=b;
         b=a+n;

         __asm__ __volatile__ ("vxorps %%ymm1, %%ymm1, %%ymm1 \n\t"
                               "vxorps %%ymm3, %%ymm3, %%ymm3 \n\t"
                               "vmovlps %0, %%xmm1, %%xmm0 \n\t"
                               "vmovlps %1, %%xmm3, %%xmm2"
                               :
                               :
                               "m" (w[0]),
                               "m" (w[1])
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovaps %0, %%xmm8 \n\t"
                                  "vmovsldup %2, %%ymm4 \n\t"
                                  "vmovshdup %2, %%ymm5"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm8");

            __asm__ __volatile__ ("vinsertf128 $0x1, %0, %%ymm8, %%ymm8 \n\t"
                                  "vmovsldup %2, %%ymm6 \n\t"
                                  "vmovshdup %2, %%ymm7 \n\t"
                                  "vpermilps $0xb1, %%ymm8, %%ymm9"
                                  :
                                  :
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (b[0]),
                                  "m" (b[1]),
                                  "m" (b[2]),
                                  "m" (b[3])
                                  :
                                  "xmm6", "xmm7", "xmm8", "xmm9");

            __asm__ __volatile__ ("vmulps %%ymm8, %%ymm4, %%ymm4 \n\t"
                                  "vmulps %%ymm9, %%ymm5, %%ymm5 \n\t"
                                  "vmulps %%ymm8, %%ymm6, %%ymm6 \n\t"
                                  "vmulps %%ymm9, %%ymm7, %%ymm7 \n\t"
                                  "vaddps %%ymm4, %%ymm0, %%ymm0 \n\t"
                                  "vaddps %%ymm5, %%ymm1, %%ymm1 \n\t"
                                  "vaddps %%ymm6, %%ymm2, %%ymm2 \n\t"
                                  "vaddps %%ymm7, %%ymm3, %%ymm3"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3",
                                  "xmm4", "xmm5", "xmm6", "xmm7");

            a+=4;
            b+=4;
         }

         __asm__ __volatile__ ("vaddsubps %%ymm1, %%ymm0, %%ymm0 \n\t"
                               "vaddsubps %%ymm3, %%ymm2, %%ymm2 \n\t"
                               "vextractf128 $0x1, %%ymm0, %%xmm8 \n\t"
                               "vextractf128 $0x1, %%ymm2, %%xmm9 \n\t"
                               "vaddps %%xmm8, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm9, %%xmm2, %%xmm2"
                               :
                               :
                               :
                               "xmm0", "xmm2", "xmm8", "xmm9");

         __asm__ __volatile__ ("vmovaps %0, %%xmm4 \n\t"
                               "vmovsldup %2, %%xmm6 \n\t"
                               "vmovshdup %2, %%xmm7 \n\t"
                               "vpermilps $0xb1, %%xmm4, %%xmm5 \n\t"
                               "vmovsldup %4, %%xmm8 \n\t"
                               "vmovshdup %4, %%xmm9 \n\t"
                               "vmulps %%xmm4, %%xmm6, %%xmm6 \n\t"
                               "vmulps %%xmm5, %%xmm7, %%xmm7 \n\t"
                               "vmulps %%xmm4, %%xmm8, %%xmm8 \n\t"
                               "vmulps %%xmm5, %%xmm9, %%xmm9 \n\t"
                               "vaddsubps %%xmm7, %%xmm6, %%xmm6 \n\t"
                               "vaddsubps %%xmm9, %%xmm8, %%xmm8 \n\t"
                               "vaddps %%xmm6, %%xmm0, %%xmm0 \n\t"
                               "vaddps %%xmm8, %%xmm2, %%xmm2"
                               :
                               :
                               "m" (vv[0]),
                               "m" (vv[1]),
                               "m" (a[0]),
                               "m" (a[1]),
                               "m" (b[0]),
                               "m" (b[1])
                               :
                               "xmm0", "xmm2", "xmm4", "xmm5",
                               "xmm6", "xmm7", "xmm8", "xmm9");

         a+=2;
         b+=2;

         __asm__ __volatile__ ("vpermilps $0xd8, %%xmm0, %%xmm0 \n\t"
                               "vpermilps $0xd8, %%xmm2, %%xmm2 \n\t"
                               "vhaddps %%xmm2, %%xmm0, %%xmm0 \n\t"
                               "vmovaps %%xmm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1])
                               :
                               :
                               "xmm0", "xmm2");
      }
   }
   else
   {
      vm=v+n-(n&0x3);
      wm=w+n;

      __asm__ __volatile__ ("vxorps %%ymm9, %%ymm9, %%ymm9"
                            :
                            :
                            :
                            "xmm9");

      for (;w<wm;w+=1)
      {
         __asm__ __volatile__ ("vmovlps %0, %%xmm9, %%xmm0 \n\t"
                               "vxorps %%ymm1, %%ymm1, %%ymm1"
                               :
                               :
                               "m" (w[0])
                               :
                               "xmm0", "xmm1");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("vmovups %0, %%ymm6 \n\t"
                                  "vmovsldup %4, %%ymm4 \n\t"
                                  "vmovshdup %4, %%ymm5 \n\t"
                                  "vpermilps $0xb1, %%ymm6, %%ymm7"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7");

            __asm__ __volatile__ ("vmulps %%ymm6, %%ymm4, %%ymm4 \n\t"
                                  "vmulps %%ymm7, %%ymm5, %%ymm5 \n\t"
                                  "vaddps %%ymm4, %%ymm0, %%ymm0 \n\t"
                                  "vaddps %%ymm5, %%ymm1, %%ymm1"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm4", "xmm5");

            a+=4;
         }

         __asm__ __volatile__ ("vaddsubps %%ymm1, %%ymm0, %%ymm2 \n\t"
                               "vextractf128 $0x1, %%ymm2, %%xmm3"
                               :
                               :
                               :
                               "xmm2", "xmm3");

         if ((n&0x2)!=0x0)
         {
            __asm__ __volatile__ ("vmovups %0, %%xmm4 \n\t"
                                  "vmovsldup %2, %%xmm6 \n\t"
                                  "vmovshdup %2, %%xmm7 \n\t"
                                  "vpermilps $0xb1, %%xmm4, %%xmm5 \n\t"
                                  "vmulps %%xmm4, %%xmm6, %%xmm6 \n\t"
                                  "vmulps %%xmm5, %%xmm7, %%xmm7 \n\t"
                                  "vaddps %%xmm6, %%xmm2, %%xmm2 \n\t"
                                  "vaddsubps %%xmm7, %%xmm3, %%xmm3"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (a[0]),
                                  "m" (a[1])
                                  :
                                  "xmm2", "xmm3", "xmm4", "xmm5",
                                  "xmm6", "xmm7");

            a+=2;
            vv+=2;
         }

         __asm__ __volatile__ ("vmovlps %0, %%xmm9, %%xmm4 \n\t"
                               "vmovlps %1, %%xmm9, %%xmm8 \n\t"
                               "vmovsldup %%xmm8, %%xmm6 \n\t"
                               "vmovshdup %%xmm8, %%xmm7 \n\t"
                               "vpermilps $0xb1, %%xmm4, %%xmm5 \n\t"
                               "vmulps %%xmm4, %%xmm6, %%xmm6 \n\t"
                               "vmulps %%xmm5, %%xmm7, %%xmm7 \n\t"
                               "vaddps %%xmm6, %%xmm2, %%xmm2 \n\t"
                               "vaddsubps %%xmm7, %%xmm3, %%xmm3"
                               :
                               :
                               "m" (vv[0]),
                               "m" (a[0])
                               :
                               "xmm2", "xmm3", "xmm4", "xmm5",
                               "xmm6", "xmm7", "xmm8");

         a+=1;

         __asm__ __volatile__ ("vaddps %%xmm3, %%xmm2, %%xmm2 \n\t"
                               "vpermilps $0xd8, %%xmm2, %%xmm0 \n\t"
                               "vhaddps %%xmm0, %%xmm0, %%xmm0 \n\t"
                               "vmovlps %%xmm0, %0"
                               :
                               "=m" (w[0])
                               :
                               :
                               "xmm0", "xmm2");
      }
   }

   _avx_zeroupper();
}

#endif

#elif (defined x64)
#include "sse2.h"

void cmat_vec(int n,complex *a,complex *v,complex *w)
{
   complex *vv,*vm,*wm;

   if ((n&0x3)==0x0)
   {
      vm=v+n;
      wm=w+n;

      for (;w<wm;w+=2)
      {
         __asm__ __volatile__ ("xorps %%xmm0, %%xmm0 \n\t"
                               "xorps %%xmm1, %%xmm1 \n\t"
                               "xorps %%xmm2, %%xmm2 \n\t"
                               "xorps %%xmm3, %%xmm3"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                                  "movaps %2, %%xmm5 \n\t"
                                  "movsldup %4, %%xmm6 \n\t"
                                  "movshdup %4, %%xmm7 \n\t"
                                  "movsldup %6, %%xmm8 \n\t"
                                  "movshdup %6, %%xmm9"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            __asm__ __volatile__ ("mulps %%xmm4, %%xmm6 \n\t"
                                  "mulps %%xmm4, %%xmm7 \n\t"
                                  "mulps %%xmm5, %%xmm8 \n\t"
                                  "mulps %%xmm5, %%xmm9 \n\t"
                                  "addps %%xmm6, %%xmm0 \n\t"
                                  "addps %%xmm7, %%xmm1 \n\t"
                                  "addps %%xmm8, %%xmm0 \n\t"
                                  "addps %%xmm9, %%xmm1"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            a+=4;
         }

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                                  "movaps %2, %%xmm5 \n\t"
                                  "movsldup %4, %%xmm6 \n\t"
                                  "movshdup %4, %%xmm7 \n\t"
                                  "movsldup %6, %%xmm8 \n\t"
                                  "movshdup %6, %%xmm9"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            __asm__ __volatile__ ("mulps %%xmm4, %%xmm6 \n\t"
                                  "mulps %%xmm4, %%xmm7 \n\t"
                                  "mulps %%xmm5, %%xmm8 \n\t"
                                  "mulps %%xmm5, %%xmm9 \n\t"
                                  "addps %%xmm6, %%xmm2 \n\t"
                                  "addps %%xmm7, %%xmm3 \n\t"
                                  "addps %%xmm8, %%xmm2 \n\t"
                                  "addps %%xmm9, %%xmm3"
                                  :
                                  :
                                  :
                                  "xmm2", "xmm3", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            a+=4;
         }

         __asm__ __volatile__ ("shufps $0xb1, %%xmm1, %%xmm1 \n\t"
                               "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                               "addsubps %%xmm1, %%xmm0 \n\t"
                               "addsubps %%xmm3, %%xmm2 \n\t"
                               "shufps $0xd8, %%xmm0, %%xmm0 \n\t"
                               "shufps $0xd8, %%xmm2, %%xmm2 \n\t"
                               "haddps %%xmm2, %%xmm0 \n\t"
                               "movaps %%xmm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1])
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");
      }
   }
   else if ((n&0x1)==0x0)
   {
      vm=v+n-2;
      wm=w+n;

      for (;w<wm;w+=2)
      {
         __asm__ __volatile__ ("xorps %%xmm0, %%xmm0 \n\t"
                               "xorps %%xmm1, %%xmm1 \n\t"
                               "xorps %%xmm2, %%xmm2 \n\t"
                               "xorps %%xmm3, %%xmm3"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                                  "movaps %2, %%xmm5 \n\t"
                                  "movsldup %4, %%xmm6 \n\t"
                                  "movshdup %4, %%xmm7 \n\t"
                                  "movsldup %6, %%xmm8 \n\t"
                                  "movshdup %6, %%xmm9"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            __asm__ __volatile__ ("mulps %%xmm4, %%xmm6 \n\t"
                                  "mulps %%xmm4, %%xmm7 \n\t"
                                  "mulps %%xmm5, %%xmm8 \n\t"
                                  "mulps %%xmm5, %%xmm9 \n\t"
                                  "addps %%xmm6, %%xmm0 \n\t"
                                  "addps %%xmm7, %%xmm1 \n\t"
                                  "addps %%xmm8, %%xmm0 \n\t"
                                  "addps %%xmm9, %%xmm1"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            a+=4;
         }

         __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                               "movsldup %2, %%xmm6 \n\t"
                               "movshdup %2, %%xmm7 \n\t"
                               "mulps %%xmm4, %%xmm6 \n\t"
                               "mulps %%xmm4, %%xmm7 \n\t"
                               "addps %%xmm6, %%xmm0 \n\t"
                               "addps %%xmm7, %%xmm1"
                               :
                               :
                               "m" (vv[0]),
                               "m" (vv[1]),
                               "m" (a[0]),
                               "m" (a[1])
                               :
                               "xmm0", "xmm1", "xmm4", "xmm6",
                               "xmm7");

         a+=2;

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                                  "movaps %2, %%xmm5 \n\t"
                                  "movsldup %4, %%xmm6 \n\t"
                                  "movshdup %4, %%xmm7 \n\t"
                                  "movsldup %6, %%xmm8 \n\t"
                                  "movshdup %6, %%xmm9"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            __asm__ __volatile__ ("mulps %%xmm4, %%xmm6 \n\t"
                                  "mulps %%xmm4, %%xmm7 \n\t"
                                  "mulps %%xmm5, %%xmm8 \n\t"
                                  "mulps %%xmm5, %%xmm9 \n\t"
                                  "addps %%xmm6, %%xmm2 \n\t"
                                  "addps %%xmm7, %%xmm3 \n\t"
                                  "addps %%xmm8, %%xmm2 \n\t"
                                  "addps %%xmm9, %%xmm3"
                                  :
                                  :
                                  :
                                  "xmm2", "xmm3", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            a+=4;
         }

         __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                               "movsldup %2, %%xmm6 \n\t"
                               "movshdup %2, %%xmm7 \n\t"
                               "mulps %%xmm4, %%xmm6 \n\t"
                               "mulps %%xmm4, %%xmm7 \n\t"
                               "addps %%xmm6, %%xmm2 \n\t"
                               "addps %%xmm7, %%xmm3"
                               :
                               :
                               "m" (vv[0]),
                               "m" (vv[1]),
                               "m" (a[0]),
                               "m" (a[1])
                               :
                               "xmm2", "xmm3", "xmm4", "xmm6",
                               "xmm7");

         __asm__ __volatile__ ("shufps $0xb1, %%xmm1, %%xmm1 \n\t"
                               "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                               "addsubps %%xmm1, %%xmm0 \n\t"
                               "addsubps %%xmm3, %%xmm2 \n\t"
                               "shufps $0xd8, %%xmm0, %%xmm0 \n\t"
                               "shufps $0xd8, %%xmm2, %%xmm2 \n\t"
                               "haddps %%xmm2, %%xmm0 \n\t"
                               "movaps %%xmm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1])
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         a+=2;
      }
   }
   else
   {
      vm=v+n-1;
      wm=w+n;

      for (;w<wm;w+=1)
      {
         __asm__ __volatile__ ("xorps %%xmm0, %%xmm0 \n\t"
                               "xorps %%xmm1, %%xmm1"
                               :
                               :
                               :
                               "xmm0", "xmm1");

         for (vv=v;vv<vm;vv+=2)
         {
            __asm__ __volatile__ ("movups %0, %%xmm3 \n\t"
                                  "movups %2, %%xmm4 \n\t"
                                  "movsldup %%xmm3, %%xmm2 \n\t"
                                  "movshdup %%xmm3, %%xmm3 \n\t"
                                  "mulps %%xmm4, %%xmm2 \n\t"
                                  "mulps %%xmm4, %%xmm3 \n\t"
                                  "addps %%xmm2, %%xmm0 \n\t"
                                  "addps %%xmm3, %%xmm1"
                                  :
                                  :
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (vv[0]),
                                  "m" (vv[1])
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3",
                                  "xmm4");

            a+=2;
         }

         __asm__ __volatile__ ("movsd %1, %%xmm3 \n\t"
                               "movsd %2, %%xmm4 \n\t"
                               "movsldup %%xmm3, %%xmm2 \n\t"
                               "movshdup %%xmm3, %%xmm3 \n\t"
                               "mulps %%xmm4, %%xmm2 \n\t"
                               "mulps %%xmm4, %%xmm3 \n\t"
                               "addps %%xmm2, %%xmm0 \n\t"
                               "addps %%xmm3, %%xmm1 \n\t"
                               "xorps %%xmm4, %%xmm4 \n\t"
                               "shufps $0xb1, %%xmm1, %%xmm1 \n\t"
                               "addsubps %%xmm1, %%xmm0 \n\t"
                               "movhlps %%xmm0, %%xmm4\n\t"
                               "addps %%xmm4, %%xmm0 \n\t"
                               "movsd %%xmm0, %0"
                               :
                               "=m" (w[0])
                               :
                               "m" (a[0]),
                               "m" (vv[0])
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm4");

         a+=1;
      }
   }
}


void cmat_vec_assign(int n,complex *a,complex *v,complex *w)
{
   complex *vv,*vm,*wm;

   if ((n&0x3)==0x0)
   {
      vm=v+n;
      wm=w+n;

      for (;w<wm;w+=2)
      {
         __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                               "xorps %%xmm1, %%xmm1 \n\t"
                               "movsd %1, %%xmm2 \n\t"
                               "xorps %%xmm3, %%xmm3"
                               :
                               :
                               "m" (w[0]),
                               "m" (w[1])
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                                  "movaps %2, %%xmm5 \n\t"
                                  "movsldup %4, %%xmm6 \n\t"
                                  "movshdup %4, %%xmm7 \n\t"
                                  "movsldup %6, %%xmm8 \n\t"
                                  "movshdup %6, %%xmm9"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            __asm__ __volatile__ ("mulps %%xmm4, %%xmm6 \n\t"
                                  "mulps %%xmm4, %%xmm7 \n\t"
                                  "mulps %%xmm5, %%xmm8 \n\t"
                                  "mulps %%xmm5, %%xmm9 \n\t"
                                  "addps %%xmm6, %%xmm0 \n\t"
                                  "addps %%xmm7, %%xmm1 \n\t"
                                  "addps %%xmm8, %%xmm0 \n\t"
                                  "addps %%xmm9, %%xmm1"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            a+=4;
         }

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                                  "movaps %2, %%xmm5 \n\t"
                                  "movsldup %4, %%xmm6 \n\t"
                                  "movshdup %4, %%xmm7 \n\t"
                                  "movsldup %6, %%xmm8 \n\t"
                                  "movshdup %6, %%xmm9"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            __asm__ __volatile__ ("mulps %%xmm4, %%xmm6 \n\t"
                                  "mulps %%xmm4, %%xmm7 \n\t"
                                  "mulps %%xmm5, %%xmm8 \n\t"
                                  "mulps %%xmm5, %%xmm9 \n\t"
                                  "addps %%xmm6, %%xmm2 \n\t"
                                  "addps %%xmm7, %%xmm3 \n\t"
                                  "addps %%xmm8, %%xmm2 \n\t"
                                  "addps %%xmm9, %%xmm3"
                                  :
                                  :
                                  :
                                  "xmm2", "xmm3", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            a+=4;
         }

         __asm__ __volatile__ ("shufps $0xb1, %%xmm1, %%xmm1 \n\t"
                               "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                               "addsubps %%xmm1, %%xmm0 \n\t"
                               "addsubps %%xmm3, %%xmm2 \n\t"
                               "shufps $0xd8, %%xmm0, %%xmm0 \n\t"
                               "shufps $0xd8, %%xmm2, %%xmm2 \n\t"
                               "haddps %%xmm2, %%xmm0 \n\t"
                               "movaps %%xmm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1])
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");
      }
   }
   else if ((n&0x1)==0x0)
   {
      vm=v+n-2;
      wm=w+n;

      for (;w<wm;w+=2)
      {
         __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                               "xorps %%xmm1, %%xmm1 \n\t"
                               "movsd %1, %%xmm2 \n\t"
                               "xorps %%xmm3, %%xmm3"
                               :
                               :
                               "m" (w[0]),
                               "m" (w[1])
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                                  "movaps %2, %%xmm5 \n\t"
                                  "movsldup %4, %%xmm6 \n\t"
                                  "movshdup %4, %%xmm7 \n\t"
                                  "movsldup %6, %%xmm8 \n\t"
                                  "movshdup %6, %%xmm9"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            __asm__ __volatile__ ("mulps %%xmm4, %%xmm6 \n\t"
                                  "mulps %%xmm4, %%xmm7 \n\t"
                                  "mulps %%xmm5, %%xmm8 \n\t"
                                  "mulps %%xmm5, %%xmm9 \n\t"
                                  "addps %%xmm6, %%xmm0 \n\t"
                                  "addps %%xmm7, %%xmm1 \n\t"
                                  "addps %%xmm8, %%xmm0 \n\t"
                                  "addps %%xmm9, %%xmm1"
                                  :
                                  :
                                  :
                                  "xmm0", "xmm1", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            a+=4;
         }

         __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                               "movsldup %2, %%xmm6 \n\t"
                               "movshdup %2, %%xmm7 \n\t"
                               "mulps %%xmm4, %%xmm6 \n\t"
                               "mulps %%xmm4, %%xmm7 \n\t"
                               "addps %%xmm6, %%xmm0 \n\t"
                               "addps %%xmm7, %%xmm1"
                               :
                               :
                               "m" (vv[0]),
                               "m" (vv[1]),
                               "m" (a[0]),
                               "m" (a[1])
                               :
                               "xmm0", "xmm1", "xmm4", "xmm6",
                               "xmm7");

         a+=2;

         for (vv=v;vv<vm;vv+=4)
         {
            __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                                  "movaps %2, %%xmm5 \n\t"
                                  "movsldup %4, %%xmm6 \n\t"
                                  "movshdup %4, %%xmm7 \n\t"
                                  "movsldup %6, %%xmm8 \n\t"
                                  "movshdup %6, %%xmm9"
                                  :
                                  :
                                  "m" (vv[0]),
                                  "m" (vv[1]),
                                  "m" (vv[2]),
                                  "m" (vv[3]),
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (a[2]),
                                  "m" (a[3])
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            __asm__ __volatile__ ("mulps %%xmm4, %%xmm6 \n\t"
                                  "mulps %%xmm4, %%xmm7 \n\t"
                                  "mulps %%xmm5, %%xmm8 \n\t"
                                  "mulps %%xmm5, %%xmm9 \n\t"
                                  "addps %%xmm6, %%xmm2 \n\t"
                                  "addps %%xmm7, %%xmm3 \n\t"
                                  "addps %%xmm8, %%xmm2 \n\t"
                                  "addps %%xmm9, %%xmm3"
                                  :
                                  :
                                  :
                                  "xmm2", "xmm3", "xmm6", "xmm7",
                                  "xmm8", "xmm9");

            a+=4;
         }

         __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                               "movsldup %2, %%xmm6 \n\t"
                               "movshdup %2, %%xmm7 \n\t"
                               "mulps %%xmm4, %%xmm6 \n\t"
                               "mulps %%xmm4, %%xmm7 \n\t"
                               "addps %%xmm6, %%xmm2 \n\t"
                               "addps %%xmm7, %%xmm3"
                               :
                               :
                               "m" (vv[0]),
                               "m" (vv[1]),
                               "m" (a[0]),
                               "m" (a[1])
                               :
                               "xmm2", "xmm3", "xmm4", "xmm6",
                               "xmm7");

         __asm__ __volatile__ ("shufps $0xb1, %%xmm1, %%xmm1 \n\t"
                               "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                               "addsubps %%xmm1, %%xmm0 \n\t"
                               "addsubps %%xmm3, %%xmm2 \n\t"
                               "shufps $0xd8, %%xmm0, %%xmm0 \n\t"
                               "shufps $0xd8, %%xmm2, %%xmm2 \n\t"
                               "haddps %%xmm2, %%xmm0 \n\t"
                               "movaps %%xmm0, %0"
                               :
                               "=m" (w[0]),
                               "=m" (w[1])
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");

         a+=2;
      }
   }
   else
   {
      vm=v+n-1;
      wm=w+n;

      for (;w<wm;w+=1)
      {
         __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                               "xorps %%xmm1, %%xmm1"
                               :
                               :
                               "m" (w[0])
                               :
                               "xmm0", "xmm1");

         for (vv=v;vv<vm;vv+=2)
         {
            __asm__ __volatile__ ("movups %0, %%xmm3 \n\t"
                                  "movups %2, %%xmm4 \n\t"
                                  "movsldup %%xmm3, %%xmm2 \n\t"
                                  "movshdup %%xmm3, %%xmm3 \n\t"
                                  "mulps %%xmm4, %%xmm2 \n\t"
                                  "mulps %%xmm4, %%xmm3 \n\t"
                                  "addps %%xmm2, %%xmm0 \n\t"
                                  "addps %%xmm3, %%xmm1"
                                  :
                                  :
                                  "m" (a[0]),
                                  "m" (a[1]),
                                  "m" (vv[0]),
                                  "m" (vv[1])
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3",
                                  "xmm4");

            a+=2;
         }

         __asm__ __volatile__ ("movsd %1, %%xmm3 \n\t"
                               "movsd %2, %%xmm4 \n\t"
                               "movsldup %%xmm3, %%xmm2 \n\t"
                               "movshdup %%xmm3, %%xmm3 \n\t"
                               "mulps %%xmm4, %%xmm2 \n\t"
                               "mulps %%xmm4, %%xmm3 \n\t"
                               "addps %%xmm2, %%xmm0 \n\t"
                               "addps %%xmm3, %%xmm1 \n\t"
                               "xorps %%xmm4, %%xmm4 \n\t"
                               "shufps $0xb1, %%xmm1, %%xmm1 \n\t"
                               "addsubps %%xmm1, %%xmm0 \n\t"
                               "movhlps %%xmm0, %%xmm4\n\t"
                               "addps %%xmm4, %%xmm0 \n\t"
                               "movsd %%xmm0, %0"
                               :
                               "=m" (w[0])
                               :
                               "m" (a[0]),
                               "m" (vv[0])
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm4");

         a+=1;
      }
   }
}

#else

void cmat_vec(int n,complex *a,complex *v,complex *w)
{
   complex *vv,*vm,*wm;

   vm=v+n;
   wm=w+n;

   for (;w<wm;w++)
   {
      (*w).re=0.0f;
      (*w).im=0.0f;

      for (vv=v;vv<vm;vv++)
      {
         (*w).re+=((*a).re*(*vv).re-(*a).im*(*vv).im);
         (*w).im+=((*a).re*(*vv).im+(*a).im*(*vv).re);
         a+=1;
      }
   }
}


void cmat_vec_assign(int n,complex *a,complex *v,complex *w)
{
   complex *vv,*vm,*wm;

   vm=v+n;
   wm=w+n;

   for (;w<wm;w++)
   {
      for (vv=v;vv<vm;vv++)
      {
         (*w).re+=((*a).re*(*vv).re-(*a).im*(*vv).im);
         (*w).im+=((*a).re*(*vv).im+(*a).im*(*vv).re);
         a+=1;
      }
   }
}

#endif

void cmat_add(int n,complex *a,complex *b,complex *c)
{
   complex *am;

   am=a+n*n;

   for (;a<am;a++)
   {
      (*c).re=(*a).re+(*b).re;
      (*c).im=(*a).im+(*b).im;
      b+=1;
      c+=1;
   }
}


void cmat_sub(int n,complex *a,complex *b,complex *c)
{
   complex *am;

   am=a+n*n;

   for (;a<am;a++)
   {
      (*c).re=(*a).re-(*b).re;
      (*c).im=(*a).im-(*b).im;
      b+=1;
      c+=1;
   }
}


void cmat_mul(int n,complex *a,complex *b,complex *c)
{
   complex *aa,*bb,*am,*bm,*bbm;

   am=a+n*n;
   bm=b+n;
   bbm=b+n*n;

   for (;a<am;a+=n)
   {
      for (;b<bm;b++)
      {
         (*c).re=0.0f;
         (*c).im=0.0f;
         aa=a;

         for (bb=b;bb<bbm;bb+=n)
         {
            (*c).re+=((*aa).re*(*bb).re-(*aa).im*(*bb).im);
            (*c).im+=((*aa).re*(*bb).im+(*aa).im*(*bb).re);
            aa+=1;
         }

         c+=1;
      }

      b-=n;
   }
}


void cmat_dag(int n,complex *a,complex *b)
{
   complex *bb,*am,*bbm;

   am=a+n*n;
   bbm=b+n*n;

   for (;a<am;)
   {
      for (bb=b;bb<bbm;bb+=n)
      {
         (*bb).re=(*a).re;
         (*bb).im=-(*a).im;
         a+=1;
      }

      b+=1;
   }
}
