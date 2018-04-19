
/*******************************************************************************
*
* File blk_solv.c
*
* Copyright (C) 2005, 2011-2013, 2016 Martin Luescher, Isabel Campos
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Solution of the Dirac equation on the blocks of the SAP_BLOCKS grid.
*
* The externally accessible functions are
*
*   void blk_mres(int n,float mu,int nmr)
*     Depending on whether the twisted-mass flag is set or not, this
*     program approximately solves (Dw+i*mu*gamma_5*1e)*b.s[0]=b.s[1] or
*     (Dw+i*mu*gamma_5)*b.s[0]=b.s[1] on the n'th block b of the SAP_BLOCKS
*     grid. The solution is obtained by applying nmr minimal residual steps,
*     using b.s[2] as workspace. On exit, the approximate solution and its
*     residue are in b.s[0] and b.s[1], respectively.
*
*   void blk_eo_mres(int n,float mu,int nmr)
*     Approximate solution of (Dwhat+i*mu*gamma_5)*b.s[0]=b.s[1] for given
*     b.s[1] on the n'th block b of the SAP_BLOCKS grid. The solution is
*     obtained by applying nmr minimal residual steps, using b.s[2] as
*     workspace. On exit, the approximate solution and its residue are in
*     b.s[0] and b.s[1], respectively, while b.s[0],b.s[1] and b.s[2] are
*     unchanged on the odd points.
*
* Notes:
*
* The twisted-mass flag is retrieved from the parameter data base (see
* flags/lat_parms.c). These programs do not perform any communications and
* can be called locally. It is taken for granted that the SAP_BLOCKS grid
* is allocated and that the gauge field and the SW term on the blocks are
* in the proper condition.
*
*******************************************************************************/

#define BLK_SOLV_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "su3.h"
#include "utils.h"
#include "sflds.h"
#include "linalg.h"
#include "block.h"
#include "dirac.h"
#include "sap.h"

static int vol;
static spinor **s;

#if (defined x64)
#include "sse2.h"

#if (defined AVX)
#include "avx.h"

static float unity=1.0f;

#if (defined FMA3)

static void scalar_prods(float *r,complex *z)
{
   spinor *s1,*s2,*sm,*smb;

   s1=s[1];
   s2=s[2];
   sm=s1+vol;

   __asm__ __volatile__ ("vxorpd %%ymm12, %%ymm12, %%ymm12 \n\t"
                         "vxorpd %%ymm13, %%ymm13, %%ymm13 \n\t"
                         "vxorpd %%ymm14, %%ymm14, %%ymm14 \n\t"
                         "vxorpd %%ymm15, %%ymm15, %%ymm15"
                         :
                         :
                         :
                         "xmm12", "xmm13", "xmm14", "xmm15");

   while (s1<sm)
   {
      __asm__ __volatile__ ("vmovaps %%ymm15, %%ymm6 \n\t"
                            "vmovaps %%ymm15, %%ymm7 \n\t"
                            "vmovaps %%ymm15, %%ymm8 \n\t"
                            "vmovaps %%ymm15, %%ymm9 \n\t"
                            "vmovaps %%ymm15, %%ymm10 \n\t"
                            "vmovaps %%ymm15, %%ymm11"
                            :
                            :
                            :
                            "xmm6", "xmm7", "xmm8", "xmm9",
                            "xmm10", "xmm11");

      smb=s1+8;

      while (s1<smb)
      {
         __asm__ __volatile__ ("vmovaps %0, %%ymm0 \n\t"
                               "vmovaps %4, %%ymm1"
                               :
                               :
                               "m" ((*s2).c1.c1),
                               "m" ((*s2).c1.c2),
                               "m" ((*s2).c1.c3),
                               "m" ((*s2).c2.c1),
                               "m" ((*s2).c2.c2),
                               "m" ((*s2).c2.c3),
                               "m" ((*s2).c3.c1),
                               "m" ((*s2).c3.c2)
                               :
                               "xmm0", "xmm1");

         __asm__ __volatile__ ("vmovsldup %0, %%ymm2 \n\t"
                               "vmovsldup %4, %%ymm3 \n\t"
                               "vmovshdup %0, %%ymm4 \n\t"
                               "vmovshdup %4, %%ymm5"
                               :
                               :
                               "m" ((*s1).c1.c1),
                               "m" ((*s1).c1.c2),
                               "m" ((*s1).c1.c3),
                               "m" ((*s1).c2.c1),
                               "m" ((*s1).c2.c2),
                               "m" ((*s1).c2.c3),
                               "m" ((*s1).c3.c1),
                               "m" ((*s1).c3.c2)
                               :
                               "xmm2", "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("vfmadd231ps %%ymm0, %%ymm0, %%ymm6 \n\t"
                               "vfmadd231ps %%ymm1, %%ymm1, %%ymm7 \n\t"
                               "vfmadd231ps %%ymm0, %%ymm2, %%ymm8 \n\t"
                               "vfmadd231ps %%ymm1, %%ymm3, %%ymm9 \n\t"
                               "vfmadd231ps %%ymm0, %%ymm4, %%ymm10 \n\t"
                               "vfmadd231ps %%ymm1, %%ymm5, %%ymm11"
                               :
                               :
                               :
                               "xmm6", "xmm7", "xmm8", "xmm9",
                               "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovaps %0, %%ymm0 \n\t"
                               "vmovaps %4, %%ymm1"
                               :
                               :
                               "m" ((*s2).c3.c3),
                               "m" ((*s2).c4.c1),
                               "m" ((*s2).c4.c2),
                               "m" ((*s2).c4.c3),
                               "m" ((*(s2+1)).c1.c1),
                               "m" ((*(s2+1)).c1.c2),
                               "m" ((*(s2+1)).c1.c3),
                               "m" ((*(s2+1)).c2.c1)
                               :
                               "xmm0", "xmm1");

         __asm__ __volatile__ ("vmovsldup %0, %%ymm2 \n\t"
                               "vmovsldup %4, %%ymm3 \n\t"
                               "vmovshdup %0, %%ymm4 \n\t"
                               "vmovshdup %4, %%ymm5"
                               :
                               :
                               "m" ((*s1).c3.c3),
                               "m" ((*s1).c4.c1),
                               "m" ((*s1).c4.c2),
                               "m" ((*s1).c4.c3),
                               "m" ((*(s1+1)).c1.c1),
                               "m" ((*(s1+1)).c1.c2),
                               "m" ((*(s1+1)).c1.c3),
                               "m" ((*(s1+1)).c2.c1)
                               :
                               "xmm2", "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("vfmadd231ps %%ymm0, %%ymm0, %%ymm6 \n\t"
                               "vfmadd231ps %%ymm1, %%ymm1, %%ymm7 \n\t"
                               "vfmadd231ps %%ymm0, %%ymm2, %%ymm8 \n\t"
                               "vfmadd231ps %%ymm1, %%ymm3, %%ymm9 \n\t"
                               "vfmadd231ps %%ymm0, %%ymm4, %%ymm10 \n\t"
                               "vfmadd231ps %%ymm1, %%ymm5, %%ymm11"
                               :
                               :
                               :
                               "xmm6", "xmm7", "xmm8", "xmm9",
                               "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovaps %0, %%ymm0 \n\t"
                               "vmovaps %4, %%ymm1"
                               :
                               :
                               "m" ((*(s2+1)).c2.c2),
                               "m" ((*(s2+1)).c2.c3),
                               "m" ((*(s2+1)).c3.c1),
                               "m" ((*(s2+1)).c3.c2),
                               "m" ((*(s2+1)).c3.c3),
                               "m" ((*(s2+1)).c4.c1),
                               "m" ((*(s2+1)).c4.c2),
                               "m" ((*(s2+1)).c4.c3)
                               :
                               "xmm0", "xmm1");

         __asm__ __volatile__ ("vmovsldup %0, %%ymm2 \n\t"
                               "vmovsldup %4, %%ymm3 \n\t"
                               "vmovshdup %0, %%ymm4 \n\t"
                               "vmovshdup %4, %%ymm5"
                               :
                               :
                               "m" ((*(s1+1)).c2.c2),
                               "m" ((*(s1+1)).c2.c3),
                               "m" ((*(s1+1)).c3.c1),
                               "m" ((*(s1+1)).c3.c2),
                               "m" ((*(s1+1)).c3.c3),
                               "m" ((*(s1+1)).c4.c1),
                               "m" ((*(s1+1)).c4.c2),
                               "m" ((*(s1+1)).c4.c3)
                               :
                               "xmm2", "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("vfmadd231ps %%ymm0, %%ymm0, %%ymm6 \n\t"
                               "vfmadd231ps %%ymm1, %%ymm1, %%ymm7 \n\t"
                               "vfmadd231ps %%ymm0, %%ymm2, %%ymm8 \n\t"
                               "vfmadd231ps %%ymm1, %%ymm3, %%ymm9 \n\t"
                               "vfmadd231ps %%ymm0, %%ymm4, %%ymm10 \n\t"
                               "vfmadd231ps %%ymm1, %%ymm5, %%ymm11"
                               :
                               :
                               :
                               "xmm6", "xmm7", "xmm8", "xmm9",
                               "xmm10", "xmm11");

         s1+=2;
         s2+=2;
      }

      __asm__ __volatile__ ("vaddps %%ymm6, %%ymm7, %%ymm2 \n\t"
                            "vaddps %%ymm8, %%ymm9, %%ymm3 \n\t"
                            "vaddps %%ymm10, %%ymm11, %%ymm4 \n\t"
                            "vextractf128 $0x1, %%ymm2, %%xmm5 \n\t"
                            "vextractf128 $0x1, %%ymm3, %%xmm0 \n\t"
                            "vextractf128 $0x1, %%ymm4, %%xmm1 \n\t"
                            "vaddps %%xmm2, %%xmm5, %%xmm6 \n\t"
                            "vaddps %%xmm3, %%xmm0, %%xmm7 \n\t"
                            "vaddps %%xmm4, %%xmm1, %%xmm8 \n\t"
                            "vcvtps2pd %%xmm6, %%ymm9 \n\t"
                            "vcvtps2pd %%xmm7, %%ymm10 \n\t"
                            "vcvtps2pd %%xmm8, %%ymm11 \n\t"
                            "vaddpd %%ymm9, %%ymm12, %%ymm12 \n\t"
                            "vaddpd %%ymm10, %%ymm13, %%ymm13 \n\t"
                            "vaddpd %%ymm11, %%ymm14, %%ymm14"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5", "xmm6", "xmm7",
                            "xmm8", "xmm9", "xmm10", "xmm11",
                            "xmm12", "xmm13", "xmm14");
   }

   __asm__ __volatile__ ("vextractf128 $0x1, %%ymm12, %%xmm0 \n\t"
                         "vextractf128 $0x1, %%ymm13, %%xmm1 \n\t"
                         "vextractf128 $0x1, %%ymm14, %%xmm2 \n\t"
                         "vaddpd %%xmm12, %%xmm0, %%xmm3 \n\t"
                         "vaddpd %%xmm13, %%xmm1, %%xmm4 \n\t"
                         "vaddpd %%xmm14, %%xmm2, %%xmm5 \n\t"
                         "vhaddpd %%xmm3, %%xmm3, %%xmm6 \n\t"
                         "vpermilpd $0x1, %%xmm4, %%xmm7 \n\t"
                         "vaddsubpd %%xmm7, %%xmm5, %%xmm8 \n\t"
                         "vpermilpd $0x1, %%xmm8, %%xmm9 \n\t"
                         "vcvtpd2ps %%xmm6, %%xmm10 \n\t"
                         "vcvtpd2ps %%xmm9, %%xmm11 \n\t"
                         "vmovss %%xmm10, %0 \n\t"
                         "vmovlps %%xmm11, %1"
                         :
                         "=m" (*r),
                         "=m" (*z)
                         :
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7",
                         "xmm8", "xmm9", "xmm10", "xmm11");

   _avx_zeroupper();
}


static void linear_cmbs(float *r,complex *z)
{
   spinor *s0,*s1,*s2,*sm;

   s0=s[0];
   s1=s[1];
   s2=s[2];
   sm=s0+vol;

   __asm__ __volatile__ ("vmovss %0, %%xmm7 \n\t"
                         "vmovss %1, %%xmm8 \n\t"
                         "vdivss %%xmm7, %%xmm8, %%xmm9 \n\t"
                         "vinsertf128 $0x1, %%xmm9, %%ymm9, %%ymm9 \n\t"
                         "vpermilps $0x0, %%ymm9, %%ymm10 \n\t"
                         "vxorps %%ymm13, %%ymm13, %%ymm13 \n\t"
                         "vbroadcastss %2, %%ymm11 \n\t"
                         "vbroadcastss %3, %%ymm12 \n\t"
                         "vaddsubps %%ymm11, %%ymm13, %%ymm13 \n\t"
                         "vmulps %%ymm10, %%ymm12, %%ymm12 \n\t"
                         "vmulps %%ymm10, %%ymm13, %%ymm13"
                         :
                         :
                         "m" (*r),
                         "m" (unity),
                         "m" ((*z).im),
                         "m" ((*z).re)
                         :
                         "xmm7", "xmm8", "xmm9", "xmm10",
                         "xmm11", "xmm12", "xmm13");

   while (s0<sm)
   {
      _avx_spinor_load(*s1);
      _avx_spinor_load_up(*s0);

      __asm__ __volatile__ ("vpermilps $0xb1, %%ymm0, %%ymm6 \n\t"
                            "vpermilps $0xb1, %%ymm1, %%ymm7 \n\t"
                            "vpermilps $0xb1, %%ymm2, %%ymm8"
                            :
                            :
                            :
                            "xmm6", "xmm7", "xmm8");

      __asm__ __volatile__ ("vmovaps %0, %%ymm9 \n\t"
                            "vmovaps %4, %%ymm10"
                            :
                            :
                            "m" ((*s2).c1.c1),
                            "m" ((*s2).c1.c2),
                            "m" ((*s2).c1.c3),
                            "m" ((*s2).c2.c1),
                            "m" ((*s2).c2.c2),
                            "m" ((*s2).c2.c3),
                            "m" ((*s2).c3.c1),
                            "m" ((*s2).c3.c2)
                            :
                            "xmm9", "xmm10");

      __asm__ __volatile__ ("vmovaps %0, %%ymm11"
                            :
                            :
                            "m" ((*s2).c3.c3),
                            "m" ((*s2).c4.c1),
                            "m" ((*s2).c4.c2),
                            "m" ((*s2).c4.c3)
                            :
                            "xmm11");

      __asm__ __volatile__ ("vfmadd231ps %%ymm12, %%ymm0, %%ymm3 \n\t"
                            "vfmadd231ps %%ymm12, %%ymm1, %%ymm4 \n\t"
                            "vfmadd231ps %%ymm12, %%ymm2, %%ymm5 \n\t"
                            "vfnmadd231ps %%ymm12, %%ymm9, %%ymm0 \n\t"
                            "vfnmadd231ps %%ymm12, %%ymm10, %%ymm1 \n\t"
                            "vfnmadd231ps %%ymm12, %%ymm11, %%ymm2 \n\t"
                            "vpermilps $0xb1, %%ymm9, %%ymm9 \n\t"
                            "vpermilps $0xb1, %%ymm10, %%ymm10 \n\t"
                            "vpermilps $0xb1, %%ymm11, %%ymm11 \n\t"
                            "vfmadd231ps %%ymm13, %%ymm6, %%ymm3 \n\t"
                            "vfmadd231ps %%ymm13, %%ymm7, %%ymm4 \n\t"
                            "vfmadd231ps %%ymm13, %%ymm8, %%ymm5 \n\t"
                            "vfnmadd231ps %%ymm13, %%ymm9, %%ymm0 \n\t"
                            "vfnmadd231ps %%ymm13, %%ymm10, %%ymm1 \n\t"
                            "vfnmadd231ps %%ymm13, %%ymm11, %%ymm2"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5", "xmm9", "xmm10",
                            "xmm11");

      _avx_spinor_store_up(*s0);
      _avx_spinor_store(*s1);

      s0+=1;
      s1+=1;
      s2+=1;
   }

   _avx_zeroupper();
}

#else

static void scalar_prods(float *r,complex *z)
{
   spinor *s1,*s2,*sm;

   __asm__ __volatile__ ("vxorpd %%ymm12, %%ymm12, %%ymm12 \n\t"
                         "vxorpd %%ymm13, %%ymm13, %%ymm13 \n\t"
                         "vxorpd %%ymm14, %%ymm14, %%ymm14"
                         :
                         :
                         :
                         "xmm12", "xmm13", "xmm14");

   s1=s[1];
   s2=s[2];
   sm=s1+vol;

   for (;s1<sm;s1++)
   {
      _avx_spinor_load(*s2);

      __asm__ __volatile__ ("vmovsldup %0, %%ymm3 \n\t"
                            "vmovsldup %4, %%ymm4"
                            :
                            :
                            "m" ((*s1).c1.c1),
                            "m" ((*s1).c1.c2),
                            "m" ((*s1).c1.c3),
                            "m" ((*s1).c2.c1),
                            "m" ((*s1).c2.c2),
                            "m" ((*s1).c2.c3),
                            "m" ((*s1).c3.c1),
                            "m" ((*s1).c3.c2)
                            :
                            "xmm3", "xmm4");

      __asm__ __volatile__ ("vmovsldup %0, %%ymm5 \n\t"
                            "vmovshdup %4, %%ymm6"
                            :
                            :
                            "m" ((*s1).c3.c3),
                            "m" ((*s1).c4.c1),
                            "m" ((*s1).c4.c2),
                            "m" ((*s1).c4.c3),
                            "m" ((*s1).c1.c1),
                            "m" ((*s1).c1.c2),
                            "m" ((*s1).c1.c3),
                            "m" ((*s1).c2.c1)
                            :
                            "xmm5", "xmm6");

      __asm__ __volatile__ ("vmovshdup %0, %%ymm7 \n\t"
                            "vmovshdup %4, %%ymm8"
                            :
                            :
                            "m" ((*s1).c2.c2),
                            "m" ((*s1).c2.c3),
                            "m" ((*s1).c3.c1),
                            "m" ((*s1).c3.c2),
                            "m" ((*s1).c3.c3),
                            "m" ((*s1).c4.c1),
                            "m" ((*s1).c4.c2),
                            "m" ((*s1).c4.c3)
                            :
                            "xmm7", "xmm8");

      __asm__ __volatile__ ("vmulps %%ymm0, %%ymm3, %%ymm3 \n\t"
                            "vmulps %%ymm1, %%ymm4, %%ymm4 \n\t"
                            "vmulps %%ymm2, %%ymm5, %%ymm5 \n\t"
                            "vmulps %%ymm0, %%ymm6, %%ymm6 \n\t"
                            "vmulps %%ymm1, %%ymm7, %%ymm7 \n\t"
                            "vmulps %%ymm2, %%ymm8, %%ymm8 \n\t"
                            "vmulps %%ymm0, %%ymm0, %%ymm0 \n\t"
                            "vmulps %%ymm1, %%ymm1, %%ymm1 \n\t"
                            "vmulps %%ymm2, %%ymm2, %%ymm2 \n\t"
                            "vaddps %%ymm4, %%ymm3, %%ymm3 \n\t"
                            "vaddps %%ymm7, %%ymm6, %%ymm6 \n\t"
                            "vaddps %%ymm1, %%ymm0, %%ymm0 \n\t"
                            "vaddps %%ymm5, %%ymm3, %%ymm3 \n\t"
                            "vaddps %%ymm8, %%ymm6, %%ymm6 \n\t"
                            "vaddps %%ymm2, %%ymm0, %%ymm0"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5", "xmm6", "xmm7",
                            "xmm8");

      __asm__ __volatile__ ("vextractf128 $0x1, %%ymm3, %%xmm9 \n\t"
                            "vextractf128 $0x1, %%ymm6, %%xmm10 \n\t"
                            "vextractf128 $0x1, %%ymm0, %%xmm11 \n\t"
                            "vaddps %%xmm9, %%xmm3, %%xmm3 \n\t"
                            "vaddps %%xmm10, %%xmm6, %%xmm6 \n\t"
                            "vaddps %%xmm11, %%xmm0, %%xmm0 \n\t"
                            "vcvtps2pd %%xmm3, %%ymm4 \n\t"
                            "vcvtps2pd %%xmm6, %%ymm7 \n\t"
                            "vcvtps2pd %%xmm0, %%ymm1 \n\t"
                            "vaddpd %%ymm4, %%ymm12, %%ymm12 \n\t"
                            "vaddpd %%ymm7, %%ymm13, %%ymm13 \n\t"
                            "vaddpd %%ymm1, %%ymm14, %%ymm14"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm3", "xmm4",
                            "xmm6", "xmm7", "xmm9", "xmm10",
                            "xmm11", "xmm12", "xmm13", "xmm14");

      s2+=1;
   }

   __asm__ __volatile__ ("vextractf128 $0x1, %%ymm12, %%xmm3 \n\t"
                         "vextractf128 $0x1, %%ymm13, %%xmm4 \n\t"
                         "vextractf128 $0x1, %%ymm14, %%xmm5 \n\t"
                         "vaddpd %%xmm3, %%xmm12, %%xmm12 \n\t"
                         "vaddpd %%xmm4, %%xmm13, %%xmm13 \n\t"
                         "vaddpd %%xmm5, %%xmm14, %%xmm14"
                         :
                         :
                         :
                         "xmm3", "xmm4", "xmm5", "xmm12",
                         "xmm13", "xmm14");

   _avx_zeroupper();

   __asm__ __volatile__ ("shufpd $0x1, %%xmm12, %%xmm12 \n\t"
                         "haddpd %%xmm14, %%xmm14 \n\t"
                         "addsubpd %%xmm12, %%xmm13 \n\t"
                         "cvtsd2ss %%xmm14, %%xmm0 \n\t"
                         "shufpd $0x1, %%xmm13, %%xmm13 \n\t"
                         "movss %%xmm0, %0 \n\t"
                         "cvtpd2ps %%xmm13, %%xmm1 \n\t"
                         "movlps %%xmm1, %1"
                         :
                         "=m" (*r),
                         "=m" (*z)
                         :
                         :
                         "xmm0", "xmm1", "xmm12", "xmm13",
                         "xmm14");
}


static void linear_cmbs(float *r,complex *z)
{
   spinor *s0,*s1,*s2,*sm;

   s0=s[0];
   s1=s[1];
   s2=s[2];
   sm=s0+vol;

   __asm__ __volatile__ ("vmovss %0, %%xmm7 \n\t"
                         "vmovss %1, %%xmm8 \n\t"
                         "vdivss %%xmm7, %%xmm8, %%xmm9 \n\t"
                         "vinsertf128 $0x1, %%xmm9, %%ymm9, %%ymm9 \n\t"
                         "vpermilps $0x0, %%ymm9, %%ymm10 \n\t"
                         "vxorps %%ymm13, %%ymm13, %%ymm13 \n\t"
                         "vbroadcastss %2, %%ymm11 \n\t"
                         "vbroadcastss %3, %%ymm12 \n\t"
                         "vaddsubps %%ymm11, %%ymm13, %%ymm13 \n\t"
                         "vmulps %%ymm10, %%ymm12, %%ymm12 \n\t"
                         "vmulps %%ymm10, %%ymm13, %%ymm13"
                         :
                         :
                         "m" (*r),
                         "m" (unity),
                         "m" ((*z).im),
                         "m" ((*z).re)
                         :
                         "xmm7", "xmm8", "xmm9", "xmm10",
                         "xmm11", "xmm12", "xmm13");

   for (;s0<sm;s0++)
   {
      _avx_mulc_spinor(*s2);
      _avx_spinor_load_up(*s1);

      __asm__ __volatile__ ("vsubps %%ymm0, %%ymm3, %%ymm0 \n\t"
                            "vsubps %%ymm1, %%ymm4, %%ymm1 \n\t"
                            "vsubps %%ymm2, %%ymm5, %%ymm2"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2");

      _avx_spinor_store(*s1);

      __asm__ __volatile__ ("vpermilps $0xb1, %%ymm3, %%ymm6 \n\t"
                            "vpermilps $0xb1, %%ymm4, %%ymm7 \n\t"
                            "vpermilps $0xb1, %%ymm5, %%ymm8 \n\t"
                            "vmulps %%ymm12, %%ymm3, %%ymm3 \n\t"
                            "vmulps %%ymm13, %%ymm6, %%ymm6 \n\t"
                            "vmulps %%ymm12, %%ymm4, %%ymm4 \n\t"
                            "vmulps %%ymm13, %%ymm7, %%ymm7 \n\t"
                            "vmulps %%ymm12, %%ymm5, %%ymm5 \n\t"
                            "vmulps %%ymm13, %%ymm8, %%ymm8"
                            :
                            :
                            :
                            "xmm3", "xmm4", "xmm5",
                            "xmm6", "xmm7", "xmm8");

      __asm__ __volatile__ ("vaddps %0, %%ymm3, %%ymm3"
                            :
                            :
                            "m" ((*s0).c1.c1),
                            "m" ((*s0).c1.c2),
                            "m" ((*s0).c1.c3),
                            "m" ((*s0).c2.c1)
                            :
                            "xmm3");

      __asm__ __volatile__ ("vaddps %0, %%ymm4, %%ymm4"
                            :
                            :
                            "m" ((*s0).c2.c2),
                            "m" ((*s0).c2.c3),
                            "m" ((*s0).c3.c1),
                            "m" ((*s0).c3.c2)
                            :
                            "xmm4");

      __asm__ __volatile__ ("vaddps %0, %%ymm5, %%ymm5"
                            :
                            :
                            "m" ((*s0).c3.c3),
                            "m" ((*s0).c4.c1),
                            "m" ((*s0).c4.c2),
                            "m" ((*s0).c4.c3)
                            :
                            "xmm5");

      __asm__ __volatile__ ("vaddps %%ymm6, %%ymm3, %%ymm3 \n\t"
                            "vaddps %%ymm7, %%ymm4, %%ymm4 \n\t"
                            "vaddps %%ymm8, %%ymm5, %%ymm5"
                            :
                            :
                            :
                            "xmm3", "xmm4", "xmm5");

      _avx_spinor_store_up(*s0);

      s1+=1;
      s2+=1;
   }

   _avx_zeroupper();
}

#endif
#else

static void scalar_prods(float *r,complex *z)
{
   spinor *s1,*s2,*sm;

   __asm__ __volatile__ ("xorpd %%xmm12, %%xmm12 \n\t"
                         "xorpd %%xmm13, %%xmm13 \n\t"
                         "xorpd %%xmm14, %%xmm14"
                         :
                         :
                         :
                         "xmm12", "xmm13", "xmm14");

   s1=s[1];
   s2=s[2];
   sm=s1+vol;

   for (;s1<sm;s1++)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %2, %%xmm1 \n\t"
                            "movaps %4, %%xmm2"
                            :
                            :
                            "m" ((*s2).c1.c1),
                            "m" ((*s2).c1.c2),
                            "m" ((*s2).c1.c3),
                            "m" ((*s2).c2.c1),
                            "m" ((*s2).c2.c2),
                            "m" ((*s2).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("movsldup %0, %%xmm4 \n\t"
                            "movsldup %2, %%xmm5 \n\t"
                            "movsldup %4, %%xmm6 \n\t"
                            "movshdup %0, %%xmm8 \n\t"
                            "movshdup %2, %%xmm9 \n\t"
                            "movshdup %4, %%xmm10"
                            :
                            :
                            "m" ((*s1).c1.c1),
                            "m" ((*s1).c1.c2),
                            "m" ((*s1).c1.c3),
                            "m" ((*s1).c2.c1),
                            "m" ((*s1).c2.c2),
                            "m" ((*s1).c2.c3)
                            :
                            "xmm4", "xmm5", "xmm6", "xmm8",
                            "xmm9", "xmm10");

      __asm__ __volatile__ ("mulps %%xmm0, %%xmm4 \n\t"
                            "mulps %%xmm1, %%xmm5 \n\t"
                            "mulps %%xmm2, %%xmm6 \n\t"
                            "mulps %%xmm0, %%xmm8 \n\t"
                            "mulps %%xmm1, %%xmm9 \n\t"
                            "mulps %%xmm2, %%xmm10 \n\t"
                            "mulps %%xmm0, %%xmm0 \n\t"
                            "mulps %%xmm1, %%xmm1 \n\t"
                            "mulps %%xmm2, %%xmm2 \n\t"
                            "addps %%xmm5, %%xmm4 \n\t"
                            "addps %%xmm9, %%xmm8 \n\t"
                            "addps %%xmm1, %%xmm0 \n\t"
                            "addps %%xmm6, %%xmm4 \n\t"
                            "addps %%xmm10, %%xmm8 \n\t"
                            "addps %%xmm2, %%xmm0"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm4",
                            "xmm5", "xmm6", "xmm8", "xmm9",
                            "xmm10");

      __asm__ __volatile__ ("movaps %0, %%xmm1 \n\t"
                            "movaps %2, %%xmm2 \n\t"
                            "movaps %4, %%xmm3"
                            :
                            :
                            "m" ((*s2).c3.c1),
                            "m" ((*s2).c3.c2),
                            "m" ((*s2).c3.c3),
                            "m" ((*s2).c4.c1),
                            "m" ((*s2).c4.c2),
                            "m" ((*s2).c4.c3)
                            :
                            "xmm1", "xmm2", "xmm3");

      __asm__ __volatile__ ("movsldup %0, %%xmm5 \n\t"
                            "movsldup %2, %%xmm6 \n\t"
                            "movsldup %4, %%xmm7 \n\t"
                            "movshdup %0, %%xmm9 \n\t"
                            "movshdup %2, %%xmm10 \n\t"
                            "movshdup %4, %%xmm11"
                            :
                            :
                            "m" ((*s1).c3.c1),
                            "m" ((*s1).c3.c2),
                            "m" ((*s1).c3.c3),
                            "m" ((*s1).c4.c1),
                            "m" ((*s1).c4.c2),
                            "m" ((*s1).c4.c3)
                            :
                            "xmm5", "xmm6", "xmm7","xmm9",
                            "xmm10", "xmm11");

      __asm__ __volatile__ ("mulps %%xmm1, %%xmm5 \n\t"
                            "mulps %%xmm2, %%xmm6 \n\t"
                            "mulps %%xmm3, %%xmm7 \n\t"
                            "mulps %%xmm1, %%xmm9 \n\t"
                            "mulps %%xmm2, %%xmm10 \n\t"
                            "mulps %%xmm3, %%xmm11 \n\t"
                            "mulps %%xmm1, %%xmm1 \n\t"
                            "mulps %%xmm2, %%xmm2 \n\t"
                            "mulps %%xmm3, %%xmm3 \n\t"
                            "addps %%xmm5, %%xmm4 \n\t"
                            "addps %%xmm9, %%xmm8 \n\t"
                            "addps %%xmm1, %%xmm0 \n\t"
                            "addps %%xmm6, %%xmm4 \n\t"
                            "addps %%xmm10, %%xmm8 \n\t"
                            "addps %%xmm2, %%xmm0 \n\t"
                            "addps %%xmm7, %%xmm4 \n\t"
                            "addps %%xmm11, %%xmm8 \n\t"
                            "addps %%xmm3, %%xmm0"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5", "xmm6", "xmm7",
                            "xmm8", "xmm9", "xmm10", "xmm11");

      __asm__ __volatile__ ("movhlps %%xmm4, %%xmm5 \n\t"
                            "movhlps %%xmm8, %%xmm9 \n\t"
                            "movhlps %%xmm0, %%xmm1 \n\t"
                            "addps %%xmm5, %%xmm4 \n\t"
                            "addps %%xmm9, %%xmm8 \n\t"
                            "addps %%xmm1, %%xmm0 \n\t"
                            "cvtps2pd %%xmm4, %%xmm6 \n\t"
                            "cvtps2pd %%xmm8, %%xmm10 \n\t"
                            "cvtps2pd %%xmm0, %%xmm2 \n\t"
                            "addpd %%xmm6, %%xmm12 \n\t"
                            "addpd %%xmm10, %%xmm13 \n\t"
                            "addpd %%xmm2, %%xmm14"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm4",
                            "xmm5", "xmm6", "xmm8", "xmm9",
                            "xmm10", "xmm12", "xmm13", "xmm14");

      s2+=1;
   }

   __asm__ __volatile__ ("shufpd $0x1, %%xmm12, %%xmm12 \n\t"
                         "haddpd %%xmm14, %%xmm14 \n\t"
                         "addsubpd %%xmm12, %%xmm13 \n\t"
                         "cvtsd2ss %%xmm14, %%xmm0 \n\t"
                         "shufpd $0x1, %%xmm13, %%xmm13 \n\t"
                         "movss %%xmm0, %0 \n\t"
                         "cvtpd2ps %%xmm13, %%xmm1 \n\t"
                         "movlps %%xmm1, %1"
                         :
                         "=m" (*r),
                         "=m" (*z)
                         :
                         :
                         "xmm0", "xmm1", "xmm12", "xmm13",
                         "xmm14");
}


static void linear_cmbs(float *r,complex *z)
{
   spinor *s0,*s1,*s2,*sm;

   s0=s[0];
   s1=s[1];
   s2=s[2];
   sm=s0+vol;

   __asm__ __volatile__ ("movss %0, %%xmm4 \n\t"
                         "movss %1, %%xmm5 \n\t"
                         "movss %2, %%xmm6 \n\t"
                         "movss %3, %%xmm7 \n\t"
                         "divss %%xmm4, %%xmm5 \n\t"
                         "mulss %%xmm5, %%xmm6 \n\t"
                         "mulss %%xmm5, %%xmm7 \n\t"
                         "shufps $0x0, %%xmm6, %%xmm6 \n\t"
                         "shufps $0x0, %%xmm7, %%xmm7 \n\t"
                         "mulps %1, %%xmm7"
                         :
                         :
                         "m" (*r),
                         "m" (_sse_sgn13),
                         "m" ((*z).re),
                         "m" ((*z).im)
                         :
                         "xmm4", "xmm5", "xmm6", "xmm7");

   for (;s0<sm;s0++)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %2, %%xmm1 \n\t"
                            "movaps %4, %%xmm2"
                            :
                            :
                            "m" ((*s2).c1.c1),
                            "m" ((*s2).c1.c2),
                            "m" ((*s2).c1.c3),
                            "m" ((*s2).c2.c1),
                            "m" ((*s2).c2.c2),
                            "m" ((*s2).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("movaps %%xmm0, %%xmm8 \n\t"
                            "movaps %%xmm1, %%xmm9 \n\t"
                            "movaps %%xmm2, %%xmm10 \n\t"
                            "mulps %%xmm6, %%xmm0 \n\t"
                            "mulps %%xmm6, %%xmm1 \n\t"
                            "mulps %%xmm6, %%xmm2 \n\t"
                            "shufps $0xb1, %%xmm8, %%xmm8 \n\t"
                            "shufps $0xb1, %%xmm9, %%xmm9 \n\t"
                            "shufps $0xb1, %%xmm10, %%xmm10 \n\t"
                            "mulps %%xmm7, %%xmm8 \n\t"
                            "mulps %%xmm7, %%xmm9 \n\t"
                            "mulps %%xmm7, %%xmm10 \n\t"
                            "addps %%xmm8, %%xmm0 \n\t"
                            "addps %%xmm9, %%xmm1 \n\t"
                            "addps %%xmm10, %%xmm2"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm8",
                            "xmm9", "xmm10");

      __asm__ __volatile__ ("movaps %0, %%xmm3 \n\t"
                            "movaps %2, %%xmm4 \n\t"
                            "movaps %4, %%xmm5 \n\t"
                            "addps %%xmm3, %%xmm0 \n\t"
                            "addps %%xmm4, %%xmm1 \n\t"
                            "addps %%xmm5, %%xmm2"
                            :
                            :
                            "m" ((*s1).c1.c1),
                            "m" ((*s1).c1.c2),
                            "m" ((*s1).c1.c3),
                            "m" ((*s1).c2.c1),
                            "m" ((*s1).c2.c2),
                            "m" ((*s1).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movaps %%xmm3, %%xmm8 \n\t"
                            "movaps %%xmm4, %%xmm9 \n\t"
                            "movaps %%xmm5, %%xmm10 \n\t"
                            "mulps %%xmm6, %%xmm3 \n\t"
                            "mulps %%xmm6, %%xmm4 \n\t"
                            "mulps %%xmm6, %%xmm5 \n\t"
                            "shufps $0xb1, %%xmm8, %%xmm8 \n\t"
                            "shufps $0xb1, %%xmm9, %%xmm9 \n\t"
                            "shufps $0xb1, %%xmm10, %%xmm10 \n\t"
                            "movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %2 \n\t"
                            "movaps %%xmm2, %4 \n\t"
                            "mulps %%xmm7, %%xmm8 \n\t"
                            "mulps %%xmm7, %%xmm9 \n\t"
                            "mulps %%xmm7, %%xmm10"
                            :
                            "=m" ((*s1).c1.c1),
                            "=m" ((*s1).c1.c2),
                            "=m" ((*s1).c1.c3),
                            "=m" ((*s1).c2.c1),
                            "=m" ((*s1).c2.c2),
                            "=m" ((*s1).c2.c3)
                            :
                            :
                            "xmm3", "xmm4", "xmm5", "xmm8",
                            "xmm9", "xmm10");

      __asm__ __volatile__ ("movaps %0, %%xmm11 \n\t"
                            "movaps %2, %%xmm12 \n\t"
                            "movaps %4, %%xmm13 \n\t"
                            "addps %%xmm8, %%xmm3 \n\t"
                            "addps %%xmm9, %%xmm4 \n\t"
                            "addps %%xmm10, %%xmm5 \n\t"
                            "subps %%xmm3, %%xmm11 \n\t"
                            "subps %%xmm4, %%xmm12 \n\t"
                            "subps %%xmm5, %%xmm13"
                            :
                            :
                            "m" ((*s0).c1.c1),
                            "m" ((*s0).c1.c2),
                            "m" ((*s0).c1.c3),
                            "m" ((*s0).c2.c1),
                            "m" ((*s0).c2.c2),
                            "m" ((*s0).c2.c3)
                            :
                            "xmm3", "xmm4", "xmm5", "xmm11",
                            "xmm12", "xmm13");

      __asm__ __volatile__ ("movaps %%xmm11, %0 \n\t"
                            "movaps %%xmm12, %2 \n\t"
                            "movaps %%xmm13, %4"
                            :
                            "=m" ((*s0).c1.c1),
                            "=m" ((*s0).c1.c2),
                            "=m" ((*s0).c1.c3),
                            "=m" ((*s0).c2.c1),
                            "=m" ((*s0).c2.c2),
                            "=m" ((*s0).c2.c3));

      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %2, %%xmm1 \n\t"
                            "movaps %4, %%xmm2"
                            :
                            :
                            "m" ((*s2).c3.c1),
                            "m" ((*s2).c3.c2),
                            "m" ((*s2).c3.c3),
                            "m" ((*s2).c4.c1),
                            "m" ((*s2).c4.c2),
                            "m" ((*s2).c4.c3)
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("movaps %%xmm0, %%xmm8 \n\t"
                            "movaps %%xmm1, %%xmm9 \n\t"
                            "movaps %%xmm2, %%xmm10 \n\t"
                            "mulps %%xmm6, %%xmm0 \n\t"
                            "mulps %%xmm6, %%xmm1 \n\t"
                            "mulps %%xmm6, %%xmm2 \n\t"
                            "shufps $0xb1, %%xmm8, %%xmm8 \n\t"
                            "shufps $0xb1, %%xmm9, %%xmm9 \n\t"
                            "shufps $0xb1, %%xmm10, %%xmm10 \n\t"
                            "mulps %%xmm7, %%xmm8 \n\t"
                            "mulps %%xmm7, %%xmm9 \n\t"
                            "mulps %%xmm7, %%xmm10 \n\t"
                            "addps %%xmm8, %%xmm0 \n\t"
                            "addps %%xmm9, %%xmm1 \n\t"
                            "addps %%xmm10, %%xmm2"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm8",
                            "xmm9", "xmm10");

      __asm__ __volatile__ ("movaps %0, %%xmm3 \n\t"
                            "movaps %2, %%xmm4 \n\t"
                            "movaps %4, %%xmm5 \n\t"
                            "addps %%xmm3, %%xmm0 \n\t"
                            "addps %%xmm4, %%xmm1 \n\t"
                            "addps %%xmm5, %%xmm2"
                            :
                            :
                            "m" ((*s1).c3.c1),
                            "m" ((*s1).c3.c2),
                            "m" ((*s1).c3.c3),
                            "m" ((*s1).c4.c1),
                            "m" ((*s1).c4.c2),
                            "m" ((*s1).c4.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movaps %%xmm3, %%xmm8 \n\t"
                            "movaps %%xmm4, %%xmm9 \n\t"
                            "movaps %%xmm5, %%xmm10 \n\t"
                            "mulps %%xmm6, %%xmm3 \n\t"
                            "mulps %%xmm6, %%xmm4 \n\t"
                            "mulps %%xmm6, %%xmm5 \n\t"
                            "shufps $0xb1, %%xmm8, %%xmm8 \n\t"
                            "shufps $0xb1, %%xmm9, %%xmm9 \n\t"
                            "shufps $0xb1, %%xmm10, %%xmm10 \n\t"
                            "movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %2 \n\t"
                            "movaps %%xmm2, %4 \n\t"
                            "mulps %%xmm7, %%xmm8 \n\t"
                            "mulps %%xmm7, %%xmm9 \n\t"
                            "mulps %%xmm7, %%xmm10"
                            :
                            "=m" ((*s1).c3.c1),
                            "=m" ((*s1).c3.c2),
                            "=m" ((*s1).c3.c3),
                            "=m" ((*s1).c4.c1),
                            "=m" ((*s1).c4.c2),
                            "=m" ((*s1).c4.c3)
                            :
                            :
                            "xmm3", "xmm4", "xmm5", "xmm8",
                            "xmm9", "xmm10");

      __asm__ __volatile__ ("movaps %0, %%xmm11 \n\t"
                            "movaps %2, %%xmm12 \n\t"
                            "movaps %4, %%xmm13 \n\t"
                            "addps %%xmm8, %%xmm3 \n\t"
                            "addps %%xmm9, %%xmm4 \n\t"
                            "addps %%xmm10, %%xmm5 \n\t"
                            "subps %%xmm3, %%xmm11 \n\t"
                            "subps %%xmm4, %%xmm12 \n\t"
                            "subps %%xmm5, %%xmm13"
                            :
                            :
                            "m" ((*s0).c3.c1),
                            "m" ((*s0).c3.c2),
                            "m" ((*s0).c3.c3),
                            "m" ((*s0).c4.c1),
                            "m" ((*s0).c4.c2),
                            "m" ((*s0).c4.c3)
                            :
                            "xmm3", "xmm4", "xmm5", "xmm11",
                            "xmm12", "xmm13");

      __asm__ __volatile__ ("movaps %%xmm11, %0 \n\t"
                            "movaps %%xmm12, %2 \n\t"
                            "movaps %%xmm13, %4"
                            :
                            "=m" ((*s0).c3.c1),
                            "=m" ((*s0).c3.c2),
                            "=m" ((*s0).c3.c3),
                            "=m" ((*s0).c4.c1),
                            "=m" ((*s0).c4.c2),
                            "=m" ((*s0).c4.c3));

      s1+=1;
      s2+=1;
   }
}

#endif

void blk_mres(int n,float mu,int nmr)
{
   int nb,isw,imr;
   float r;
   complex z;
   block_t *b;

   b=blk_list(SAP_BLOCKS,&nb,&isw);
   vol=(*b).vol;
   s=(*b).s;

   set_s2zero(vol,s[0]);

   for (imr=0;imr<nmr;imr++)
   {
      Dw_blk(SAP_BLOCKS,n,mu,1,2);
      scalar_prods(&r,&z);

      if (r<(2.0f*FLT_MIN))
         return;

      linear_cmbs(&r,&z);
   }
}


void blk_eo_mres(int n,float mu,int nmr)
{
   int nb,isw,imr;
   float r;
   complex z;
   block_t *b;

   b=blk_list(SAP_BLOCKS,&nb,&isw);
   vol=(*b).vol/2;
   s=(*b).s;

   set_s2zero(vol,s[0]);

   for (imr=0;imr<nmr;imr++)
   {
      Dwhat_blk(SAP_BLOCKS,n,mu,1,2);
      scalar_prods(&r,&z);

      if (r<(2.0f*FLT_MIN))
         return;

      linear_cmbs(&r,&z);
   }
}

#else

void blk_mres(int n,float mu,int nmr)
{
   int nb,isw,imr;
   float r;
   complex z;
   block_t *b;

   b=blk_list(SAP_BLOCKS,&nb,&isw);
   vol=(*b).vol;
   s=(*b).s;

   set_s2zero(vol,s[0]);

   for (imr=0;imr<nmr;imr++)
   {
      Dw_blk(SAP_BLOCKS,n,mu,1,2);
      r=norm_square(vol,0,s[2]);

      if (r<(2.0f*FLT_MIN))
         return;

      z=spinor_prod(vol,0,s[2],s[1]);

      r=1.0f/r;
      z.re*=r;
      z.im*=r;
      mulc_spinor_add(vol,s[0],s[1],z);

      z.re=-z.re;
      z.im=-z.im;
      mulc_spinor_add(vol,s[1],s[2],z);
   }
}


void blk_eo_mres(int n,float mu,int nmr)
{
   int nb,isw,imr;
   float r;
   complex z;
   block_t *b;

   b=blk_list(SAP_BLOCKS,&nb,&isw);
   vol=(*b).vol/2;
   s=(*b).s;

   set_s2zero(vol,s[0]);

   for (imr=0;imr<nmr;imr++)
   {
      Dwhat_blk(SAP_BLOCKS,n,mu,1,2);
      r=norm_square(vol,0,s[2]);

      if (r<(2.0f*FLT_MIN))
         return;

      z=spinor_prod(vol,0,s[2],s[1]);

      r=1.0f/r;
      z.re*=r;
      z.im*=r;
      mulc_spinor_add(vol,s[0],s[1],z);

      z.re=-z.re;
      z.im=-z.im;
      mulc_spinor_add(vol,s[1],s[2],z);
   }
}

#endif
