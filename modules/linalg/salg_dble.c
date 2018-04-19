
/*******************************************************************************
*
* File salg_dble.c
*
* Copyright (C) 2005, 2007, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic linear algebra routines for double-precision Dirac fields.
*
* The externally accessible functions are
*
*   complex_dble spinor_prod_dble(int vol,int icom,spinor_dble *s,
*                                 spinor_dble *r)
*     Computes the scalar product of the fields s and r.
*
*   double spinor_prod_re_dble(int vol,int icom,spinor_dble *s,
*                              spinor_dble *r)
*     Computes the real part of the scalar product of the fields
*     s and r.
*
*   complex_dble spinor_prod5_dble(int vol,int icom,spinor_dble *s,
*                                  spinor_dble *r)
*     Computes the scalar product of the fields s and gamma_5*r.
*
*   double norm_square_dble(int vol,int icom,spinor_dble *s)
*     Computes the square of the norm of the field s.
*
*   void mulc_spinor_add_dble(int vol,spinor_dble *s,spinor_dble *r,
*                             complex_dble z)
*     Replaces the field s by s+z*r.
*
*   void mulr_spinor_add_dble(int vol,spinor_dble *s,spinor_dble *r,
*                             double c)
*     Replaces the field s by s+c*r.
*
*   void combine_spinor_dble(int vol,spinor_dble *s,spinor_dble *r,
*                            double cs,double cr)
*     Replaces the field s by cs*s+cr*r.
*
*   void project_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
*     Replaces the field s by s-(r,s)*r.
*
*   void scale_dble(int vol,double c,spinor_dble *s)
*     Replaces the field s by c*s.
*
*   double normalize_dble(int vol,int icom,spinor_dble *s)
*     Replaces the field s by s/||s|| and returns the norm ||s||.
*
*   void rotate_dble(int vol,int n,spinor_dble **ppk,complex_dble *v)
*     Replaces the fields pk by sum_j pj*v[n*j+k] where 0<=k,j<n and
*     pk=ppk[k].
*
*   void mulg5_dble(int vol,spinor_dble *s)
*     Multiplies the field s with gamma_5.
*
*   void mulmg5_dble(int vol,spinor_dble *s)
*     Multiplies the field s with -gamma_5.
*
* Notes:
*
* All these programs act on arrays of spinor fields whose base address
* is passed through the arguments. The length of the arrays is specified
* by the parameter vol. Scalar products are globally summed if the
* parameter icom is equal to 1. In this case, the calculated values are
* guaranteed to be exactly the same on all processes.
*
* The programs perform no communications except in the case of the scalar
* products if these are globally summed. If SSE (AVX) instructions are used,
* the spinor fields must be aligned to a 16 (32) byte boundary.
*
*******************************************************************************/

#define SALG_DBLE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "sflds.h"
#include "linalg.h"
#include "global.h"

static int nrot=0;
static int isx,isz,init=0;
static double smx ALIGNED8;
static complex_dble smz ALIGNED16;
static spinor_dble *psi;


static void alloc_wrotate(int n)
{
   if (nrot>0)
      afree(psi);

   psi=amalloc(n*sizeof(*psi),ALIGN);
   error_loc(psi==NULL,1,"alloc_wrotate [salg_dble.c]",
             "Unable to allocate workspace");
   set_sd2zero(n,psi);
   nrot=n;
}

#if (defined AVX)
#include "avx.h"

#if (defined FMA3)

complex_dble spinor_prod_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isz);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("vxorpd %%ymm0, %%ymm0, %%ymm0 \n\t"
                            "vxorpd %%ymm1, %%ymm1, %%ymm1 \n\t"
                            "vxorpd %%ymm2, %%ymm2, %%ymm2 \n\t"
                            "vxorpd %%ymm3, %%ymm3, %%ymm3 \n\t"
                            "vxorpd %%ymm4, %%ymm4, %%ymm4 \n\t"
                            "vxorpd %%ymm5, %%ymm5, %%ymm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      for (;s<smb;s++)
      {
         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmovddup %0, %%ymm9 \n\t"
                               "vmovddup %2, %%ymm10 \n\t"
                               "vmovddup %4, %%ymm11"
                               :
                               :
                               "m" ((*s).c1.c1.re),
                               "m" ((*s).c1.c2.re),
                               "m" ((*s).c1.c3.re),
                               "m" ((*s).c2.c1.re),
                               "m" ((*s).c2.c2.re),
                               "m" ((*s).c2.c3.re)
                               :
                               "xmm9", "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovddup %0, %%ymm12 \n\t"
                               "vmovddup %2, %%ymm13 \n\t"
                               "vmovddup %4, %%ymm14"
                               :
                               :
                               "m" ((*s).c1.c1.im),
                               "m" ((*s).c1.c2.im),
                               "m" ((*s).c1.c3.im),
                               "m" ((*s).c2.c1.im),
                               "m" ((*s).c2.c2.im),
                               "m" ((*s).c2.c3.im)
                               :
                               "xmm12", "xmm13", "xmm14");

         __asm__ __volatile__ ("vfmadd231pd %%ymm6, %%ymm9, %%ymm0 \n\t"
                               "vfmadd231pd %%ymm7, %%ymm10, %%ymm1 \n\t"
                               "vfmadd231pd %%ymm8, %%ymm11, %%ymm2 \n\t"
                               "vfnmadd231pd %%ymm6, %%ymm12, %%ymm3 \n\t"
                               "vfnmadd231pd %%ymm7, %%ymm13, %%ymm4 \n\t"
                               "vfnmadd231pd %%ymm8, %%ymm14, %%ymm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm4", "xmm5");

         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmovddup %0, %%ymm9 \n\t"
                               "vmovddup %2, %%ymm10 \n\t"
                               "vmovddup %4, %%ymm11"
                               :
                               :
                               "m" ((*s).c3.c1.re),
                               "m" ((*s).c3.c2.re),
                               "m" ((*s).c3.c3.re),
                               "m" ((*s).c4.c1.re),
                               "m" ((*s).c4.c2.re),
                               "m" ((*s).c4.c3.re)
                               :
                               "xmm9", "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovddup %0, %%ymm12 \n\t"
                               "vmovddup %2, %%ymm13 \n\t"
                               "vmovddup %4, %%ymm14"
                               :
                               :
                               "m" ((*s).c3.c1.im),
                               "m" ((*s).c3.c2.im),
                               "m" ((*s).c3.c3.im),
                               "m" ((*s).c4.c1.im),
                               "m" ((*s).c4.c2.im),
                               "m" ((*s).c4.c3.im)
                               :
                               "xmm12", "xmm13", "xmm14");

         __asm__ __volatile__ ("vfmadd231pd %%ymm6, %%ymm9, %%ymm0 \n\t"
                               "vfmadd231pd %%ymm7, %%ymm10, %%ymm1 \n\t"
                               "vfmadd231pd %%ymm8, %%ymm11, %%ymm2 \n\t"
                               "vfnmadd231pd %%ymm6, %%ymm12, %%ymm3 \n\t"
                               "vfnmadd231pd %%ymm7, %%ymm13, %%ymm4 \n\t"
                               "vfnmadd231pd %%ymm8, %%ymm14, %%ymm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm4", "xmm5");

         r+=1;
      }

      __asm__ __volatile__ ("vaddpd %%ymm1, %%ymm0, %%ymm0 \n\t"
                            "vaddpd %%ymm4, %%ymm3, %%ymm3 \n\t"
                            "vaddpd %%ymm2, %%ymm0, %%ymm0 \n\t"
                            "vaddpd %%ymm5, %%ymm3, %%ymm3 \n\t"
                            "vpermilpd $0x5, %%ymm3, %%ymm3 \n\t"
                            "vaddsubpd %%ymm3, %%ymm0, %%ymm0 \n\t"
                            "vextractf128 $0x1, %%ymm0, %%xmm1 \n\t"
                            "vaddpd %%xmm1, %%xmm0, %%xmm0 \n\t"
                            "vmovapd %%xmm0, %0"
                            :
                            "=m" (smz)
                            :
                            :
                            "xmm0", "xmm1", "xmm3");

      _avx_zeroupper();

      add_to_hsum(isz,(double*)(&smz));
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isz,(double*)(&smz));
   else
      local_hsum(isz,(double*)(&smz));

   return smz;
}


double spinor_prod_re_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isx);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("vxorpd %%ymm0, %%ymm0, %%ymm0 \n\t"
                            "vxorpd %%ymm1, %%ymm1, %%ymm1 \n\t"
                            "vxorpd %%ymm2, %%ymm2, %%ymm2"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2");

      for (;s<smb;s++)
      {
         __asm__ __volatile__ ("vmovapd %0, %%ymm3 \n\t"
                               "vmovapd %2, %%ymm4 \n\t"
                               "vmovapd %4, %%ymm5"
                               :
                               :
                               "m" ((*s).c1.c1),
                               "m" ((*s).c1.c2),
                               "m" ((*s).c1.c3),
                               "m" ((*s).c2.c1),
                               "m" ((*s).c2.c2),
                               "m" ((*s).c2.c3)
                               :
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vfmadd231pd %%ymm3, %%ymm6, %%ymm0 \n\t"
                               "vfmadd231pd %%ymm4, %%ymm7, %%ymm1 \n\t"
                               "vfmadd231pd %%ymm5, %%ymm8, %%ymm2"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2");

         __asm__ __volatile__ ("vmovapd %0, %%ymm9 \n\t"
                               "vmovapd %2, %%ymm10 \n\t"
                               "vmovapd %4, %%ymm11"
                               :
                               :
                               "m" ((*s).c3.c1),
                               "m" ((*s).c3.c2),
                               "m" ((*s).c3.c3),
                               "m" ((*s).c4.c1),
                               "m" ((*s).c4.c2),
                               "m" ((*s).c4.c3)
                               :
                               "xmm9", "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovapd %0, %%ymm12 \n\t"
                               "vmovapd %2, %%ymm13 \n\t"
                               "vmovapd %4, %%ymm14"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm12", "xmm13", "xmm14");

         __asm__ __volatile__ ("vfmadd231pd %%ymm9, %%ymm12, %%ymm0 \n\t"
                               "vfmadd231pd %%ymm10, %%ymm13, %%ymm1 \n\t"
                               "vfmadd231pd %%ymm11, %%ymm14, %%ymm2"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2");

         r+=1;
      }

      __asm__ __volatile__ ("vaddpd %%ymm1, %%ymm0, %%ymm0 \n\t"
                            "vaddpd %%ymm2, %%ymm0, %%ymm0 \n\t"
                            "vextractf128 $0x1, %%ymm0, %%xmm1 \n\t"
                            "vaddpd %%xmm1, %%xmm0, %%xmm0 \n\t"
                            "vhaddpd %%xmm0, %%xmm0, %%xmm2 \n\t"
                            "vmovsd %%xmm2, %0"
                            :
                            "=m" (smx)
                            :
                            :
                            "xmm0", "xmm1", "xmm2");

      _avx_zeroupper();

      add_to_hsum(isx,&smx);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isx,&smx);
   else
      local_hsum(isx,&smx);

   return smx;
}


complex_dble spinor_prod5_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isz);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("vxorpd %%ymm0, %%ymm0, %%ymm0 \n\t"
                            "vxorpd %%ymm1, %%ymm1, %%ymm1 \n\t"
                            "vxorpd %%ymm2, %%ymm2, %%ymm2 \n\t"
                            "vxorpd %%ymm3, %%ymm3, %%ymm3 \n\t"
                            "vxorpd %%ymm4, %%ymm4, %%ymm4 \n\t"
                            "vxorpd %%ymm5, %%ymm5, %%ymm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      for (;s<smb;s++)
      {
         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmovddup %0, %%ymm9 \n\t"
                               "vmovddup %2, %%ymm10 \n\t"
                               "vmovddup %4, %%ymm11"
                               :
                               :
                               "m" ((*s).c1.c1.re),
                               "m" ((*s).c1.c2.re),
                               "m" ((*s).c1.c3.re),
                               "m" ((*s).c2.c1.re),
                               "m" ((*s).c2.c2.re),
                               "m" ((*s).c2.c3.re)
                               :
                               "xmm9", "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovddup %0, %%ymm12 \n\t"
                               "vmovddup %2, %%ymm13 \n\t"
                               "vmovddup %4, %%ymm14"
                               :
                               :
                               "m" ((*s).c1.c1.im),
                               "m" ((*s).c1.c2.im),
                               "m" ((*s).c1.c3.im),
                               "m" ((*s).c2.c1.im),
                               "m" ((*s).c2.c2.im),
                               "m" ((*s).c2.c3.im)
                               :
                               "xmm12", "xmm13", "xmm14");

         __asm__ __volatile__ ("vfmadd231pd %%ymm6, %%ymm9, %%ymm0 \n\t"
                               "vfmadd231pd %%ymm7, %%ymm10, %%ymm1 \n\t"
                               "vfmadd231pd %%ymm8, %%ymm11, %%ymm2 \n\t"
                               "vfnmadd231pd %%ymm6, %%ymm12, %%ymm3 \n\t"
                               "vfnmadd231pd %%ymm7, %%ymm13, %%ymm4 \n\t"
                               "vfnmadd231pd %%ymm8, %%ymm14, %%ymm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm4", "xmm5");

         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmovddup %0, %%ymm9 \n\t"
                               "vmovddup %2, %%ymm10 \n\t"
                               "vmovddup %4, %%ymm11"
                               :
                               :
                               "m" ((*s).c3.c1.re),
                               "m" ((*s).c3.c2.re),
                               "m" ((*s).c3.c3.re),
                               "m" ((*s).c4.c1.re),
                               "m" ((*s).c4.c2.re),
                               "m" ((*s).c4.c3.re)
                               :
                               "xmm9", "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovddup %0, %%ymm12 \n\t"
                               "vmovddup %2, %%ymm13 \n\t"
                               "vmovddup %4, %%ymm14"
                               :
                               :
                               "m" ((*s).c3.c1.im),
                               "m" ((*s).c3.c2.im),
                               "m" ((*s).c3.c3.im),
                               "m" ((*s).c4.c1.im),
                               "m" ((*s).c4.c2.im),
                               "m" ((*s).c4.c3.im)
                               :
                               "xmm12", "xmm13", "xmm14");

         __asm__ __volatile__ ("vfnmadd231pd %%ymm6, %%ymm9, %%ymm0 \n\t"
                               "vfnmadd231pd %%ymm7, %%ymm10, %%ymm1 \n\t"
                               "vfnmadd231pd %%ymm8, %%ymm11, %%ymm2 \n\t"
                               "vfmadd231pd %%ymm6, %%ymm12, %%ymm3 \n\t"
                               "vfmadd231pd %%ymm7, %%ymm13, %%ymm4 \n\t"
                               "vfmadd231pd %%ymm8, %%ymm14, %%ymm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3",
                               "xmm4", "xmm5");

         r+=1;
      }

      __asm__ __volatile__ ("vaddpd %%ymm1, %%ymm0, %%ymm0 \n\t"
                            "vaddpd %%ymm4, %%ymm3, %%ymm3 \n\t"
                            "vaddpd %%ymm2, %%ymm0, %%ymm0 \n\t"
                            "vaddpd %%ymm5, %%ymm3, %%ymm3 \n\t"
                            "vpermilpd $0x5, %%ymm3, %%ymm3 \n\t"
                            "vaddsubpd %%ymm3, %%ymm0, %%ymm0 \n\t"
                            "vextractf128 $0x1, %%ymm0, %%xmm1 \n\t"
                            "vaddpd %%xmm1, %%xmm0, %%xmm0 \n\t"
                            "vmovapd %%xmm0, %0"
                            :
                            "=m" (smz)
                            :
                            :
                            "xmm0", "xmm1", "xmm3");

      _avx_zeroupper();

      add_to_hsum(isz,(double*)(&smz));
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isz,(double*)(&smz));
   else
      local_hsum(isz,(double*)(&smz));

   return smz;
}


double norm_square_dble(int vol,int icom,spinor_dble *s)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isx);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("vxorpd %%ymm6, %%ymm6, %%ymm6 \n\t"
                            "vxorpd %%ymm7, %%ymm7, %%ymm7 \n\t"
                            "vxorpd %%ymm8, %%ymm8, %%ymm8 \n\t"
                            "vxorpd %%ymm9, %%ymm9, %%ymm9 \n\t"
                            "vxorpd %%ymm10, %%ymm10, %%ymm10 \n\t"
                            "vxorpd %%ymm11, %%ymm11, %%ymm11"
                            :
                            :
                            :
                            "xmm6", "xmm7", "xmm8", "xmm9",
                            "xmm10", "xmm11");

      for (;s<smb;s++)
      {
         _avx_spinor_load_dble(*s);

         __asm__ __volatile__ ("vfmadd231pd %%ymm0, %%ymm0, %%ymm6 \n\t"
                               "vfmadd231pd %%ymm1, %%ymm1, %%ymm7 \n\t"
                               "vfmadd231pd %%ymm2, %%ymm2, %%ymm8 \n\t"
                               "vfmadd231pd %%ymm3, %%ymm3, %%ymm9 \n\t"
                               "vfmadd231pd %%ymm4, %%ymm4, %%ymm10 \n\t"
                               "vfmadd231pd %%ymm5, %%ymm5, %%ymm11 \n\t"
                               :
                               :
                               :
                               "xmm6", "xmm7", "xmm8", "xmm9",
                               "xmm10", "xmm11");

      }

      __asm__ __volatile__ ("vaddpd %%ymm6, %%ymm7, %%ymm7 \n\t"
                            "vaddpd %%ymm8, %%ymm9, %%ymm9 \n\t"
                            "vaddpd %%ymm10, %%ymm11, %%ymm11 \n\t"
                            "vaddpd %%ymm7, %%ymm9, %%ymm9 \n\t"
                            "vaddpd %%ymm9, %%ymm11, %%ymm11 \n\t"
                            "vextractf128 $0x1, %%ymm11, %%xmm12 \n\t"
                            "vaddpd %%xmm11, %%xmm12, %%xmm12 \n\t"
                            "vhaddpd %%xmm12, %%xmm12, %%xmm13 \n\t"
                            "vmovsd %%xmm13, %0 \n\t"
                            :
                            "=m" (smx)
                            :
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm7", "xmm9", "xmm11");

      _avx_zeroupper();

      add_to_hsum(isx,&smx);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isx,&smx);
   else
      local_hsum(isx,&smx);

   return smx;
}

#else

complex_dble spinor_prod_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isz);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("vxorpd %%ymm0, %%ymm0, %%ymm0 \n\t"
                            "vxorpd %%ymm1, %%ymm1, %%ymm1 \n\t"
                            "vxorpd %%ymm2, %%ymm2, %%ymm2 \n\t"
                            "vxorpd %%ymm3, %%ymm3, %%ymm3 \n\t"
                            "vxorpd %%ymm4, %%ymm4, %%ymm4 \n\t"
                            "vxorpd %%ymm5, %%ymm5, %%ymm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm3", "xmm4", "xmm5");

      for (;s<smb;s++)
      {
         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*s).c1.c1),
                               "m" ((*s).c1.c2),
                               "m" ((*s).c1.c3),
                               "m" ((*s).c2.c1),
                               "m" ((*s).c2.c2),
                               "m" ((*s).c2.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmovapd %0, %%ymm12 \n\t"
                               "vmovapd %2, %%ymm13 \n\t"
                               "vmovapd %4, %%ymm14"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm12", "xmm13", "xmm14");

         __asm__ __volatile__ ("vpermilpd $0x5, %%ymm6, %%ymm9 \n\t"
                               "vpermilpd $0x5, %%ymm7, %%ymm10 \n\t"
                               "vpermilpd $0x5, %%ymm8, %%ymm11 \n\t"
                               "vmulpd %%ymm12, %%ymm6, %%ymm6 \n\t"
                               "vmulpd %%ymm13, %%ymm7, %%ymm7 \n\t"
                               "vmulpd %%ymm14, %%ymm8, %%ymm8 \n\t"
                               "vmulpd %%ymm12, %%ymm9, %%ymm9 \n\t"
                               "vmulpd %%ymm13, %%ymm10, %%ymm10 \n\t"
                               "vmulpd %%ymm14, %%ymm11, %%ymm11 \n\t"
                               "vaddpd %%ymm6, %%ymm0, %%ymm0 \n\t"
                               "vaddpd %%ymm7, %%ymm1, %%ymm1 \n\t"
                               "vaddpd %%ymm8, %%ymm2, %%ymm2 \n\t"
                               "vaddsubpd %%ymm9, %%ymm3, %%ymm3 \n\t"
                               "vaddsubpd %%ymm10, %%ymm4, %%ymm4 \n\t"
                               "vaddsubpd %%ymm11, %%ymm5, %%ymm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5",
                               "xmm6", "xmm7", "xmm8",
                               "xmm9", "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*s).c3.c1),
                               "m" ((*s).c3.c2),
                               "m" ((*s).c3.c3),
                               "m" ((*s).c4.c1),
                               "m" ((*s).c4.c2),
                               "m" ((*s).c4.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmovapd %0, %%ymm12 \n\t"
                               "vmovapd %2, %%ymm13 \n\t"
                               "vmovapd %4, %%ymm14"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm12", "xmm13", "xmm14");

         __asm__ __volatile__ ("vpermilpd $0x5, %%ymm6, %%ymm9 \n\t"
                               "vpermilpd $0x5, %%ymm7, %%ymm10 \n\t"
                               "vpermilpd $0x5, %%ymm8, %%ymm11 \n\t"
                               "vmulpd %%ymm12, %%ymm6, %%ymm6 \n\t"
                               "vmulpd %%ymm13, %%ymm7, %%ymm7 \n\t"
                               "vmulpd %%ymm14, %%ymm8, %%ymm8 \n\t"
                               "vmulpd %%ymm12, %%ymm9, %%ymm9 \n\t"
                               "vmulpd %%ymm13, %%ymm10, %%ymm10 \n\t"
                               "vmulpd %%ymm14, %%ymm11, %%ymm11 \n\t"
                               "vaddpd %%ymm6, %%ymm0, %%ymm0 \n\t"
                               "vaddpd %%ymm7, %%ymm1, %%ymm1 \n\t"
                               "vaddpd %%ymm8, %%ymm2, %%ymm2 \n\t"
                               "vaddsubpd %%ymm9, %%ymm3, %%ymm3 \n\t"
                               "vaddsubpd %%ymm10, %%ymm4, %%ymm4 \n\t"
                               "vaddsubpd %%ymm11, %%ymm5, %%ymm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5",
                               "xmm6", "xmm7", "xmm8",
                               "xmm9", "xmm10", "xmm11");

         r+=1;
      }

      __asm__ __volatile__ ("vaddpd %%ymm0, %%ymm2, %%ymm2 \n\t"
                            "vaddpd %%ymm3, %%ymm5, %%ymm5 \n\t"
                            "vaddpd %%ymm1, %%ymm2, %%ymm2 \n\t"
                            "vaddpd %%ymm4, %%ymm5, %%ymm5 \n\t"
                            "vhaddpd %%ymm5, %%ymm2, %%ymm0 \n\t"
                            "vextractf128 $0x1, %%ymm0, %%xmm1 \n\t"
                            "vaddpd %%xmm0, %%xmm1, %%xmm2 \n\t"
                            "vmovapd %%xmm2, %0 \n\t"
                            "vzeroupper"
                            :
                            "=m" (smz)
                            :
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm5");

      add_to_hsum(isz,(double*)(&smz));
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isz,(double*)(&smz));
   else
      local_hsum(isz,(double*)(&smz));

   return smz;
}


double spinor_prod_re_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isx);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("vxorpd %%ymm0, %%ymm0, %%ymm0 \n\t"
                            "vxorpd %%ymm1, %%ymm1, %%ymm1 \n\t"
                            "vxorpd %%ymm2, %%ymm2, %%ymm2 \n\t"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2");

      for (;s<smb;s++)
      {
         __asm__ __volatile__ ("vmovapd %0, %%ymm3 \n\t"
                               "vmovapd %2, %%ymm4 \n\t"
                               "vmovapd %4, %%ymm5"
                               :
                               :
                               "m" ((*s).c1.c1),
                               "m" ((*s).c1.c2),
                               "m" ((*s).c1.c3),
                               "m" ((*s).c2.c1),
                               "m" ((*s).c2.c2),
                               "m" ((*s).c2.c3)
                               :
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmulpd %%ymm3, %%ymm6, %%ymm6 \n\t"
                               "vmulpd %%ymm4, %%ymm7, %%ymm7 \n\t"
                               "vmulpd %%ymm5, %%ymm8, %%ymm8 \n\t"
                               "vaddpd %%ymm6, %%ymm0, %%ymm0 \n\t"
                               "vaddpd %%ymm7, %%ymm1, %%ymm1 \n\t"
                               "vaddpd %%ymm8, %%ymm2, %%ymm2 \n\t"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmovapd %0, %%ymm3 \n\t"
                               "vmovapd %2, %%ymm4 \n\t"
                               "vmovapd %4, %%ymm5"
                               :
                               :
                               "m" ((*s).c3.c1),
                               "m" ((*s).c3.c2),
                               "m" ((*s).c3.c3),
                               "m" ((*s).c4.c1),
                               "m" ((*s).c4.c2),
                               "m" ((*s).c4.c3)
                               :
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmulpd %%ymm3, %%ymm6, %%ymm6 \n\t"
                               "vmulpd %%ymm4, %%ymm7, %%ymm7 \n\t"
                               "vmulpd %%ymm5, %%ymm8, %%ymm8 \n\t"
                               "vaddpd %%ymm6, %%ymm0, %%ymm0 \n\t"
                               "vaddpd %%ymm7, %%ymm1, %%ymm1 \n\t"
                               "vaddpd %%ymm8, %%ymm2, %%ymm2 \n\t"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm6", "xmm7", "xmm8");

         r+=1;
      }

      __asm__ __volatile__ ("vaddpd %%ymm0, %%ymm2, %%ymm2 \n\t"
                            "vaddpd %%ymm1, %%ymm2, %%ymm2 \n\t"
                            "vextractf128 $0x1, %%ymm2, %%xmm1 \n\t"
                            "vaddpd %%xmm1, %%xmm2, %%xmm2 \n\t"
                            "vhaddpd %%xmm2, %%xmm2, %%xmm0 \n\t"
                            "vmovsd %%xmm0, %0 \n\t"
                            "vzeroupper"
                            :
                            "=m" (smx)
                            :
                            :
                            "xmm0", "xmm1", "xmm2");

      add_to_hsum(isx,&smx);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isx,&smx);
   else
      local_hsum(isx,&smx);

   return smx;
}


complex_dble spinor_prod5_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isz);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("vxorpd %%ymm0, %%ymm0, %%ymm0 \n\t"
                            "vxorpd %%ymm1, %%ymm1, %%ymm1 \n\t"
                            "vxorpd %%ymm2, %%ymm2, %%ymm2 \n\t"
                            "vxorpd %%ymm3, %%ymm3, %%ymm3 \n\t"
                            "vxorpd %%ymm4, %%ymm4, %%ymm4 \n\t"
                            "vxorpd %%ymm5, %%ymm5, %%ymm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm3", "xmm4", "xmm5");

      for (;s<smb;s++)
      {
         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*s).c1.c1),
                               "m" ((*s).c1.c2),
                               "m" ((*s).c1.c3),
                               "m" ((*s).c2.c1),
                               "m" ((*s).c2.c2),
                               "m" ((*s).c2.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmovapd %0, %%ymm12 \n\t"
                               "vmovapd %2, %%ymm13 \n\t"
                               "vmovapd %4, %%ymm14"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm12", "xmm13", "xmm14");

         __asm__ __volatile__ ("vpermilpd $0x5, %%ymm6, %%ymm9 \n\t"
                               "vpermilpd $0x5, %%ymm7, %%ymm10 \n\t"
                               "vpermilpd $0x5, %%ymm8, %%ymm11 \n\t"
                               "vmulpd %%ymm12, %%ymm6, %%ymm6 \n\t"
                               "vmulpd %%ymm13, %%ymm7, %%ymm7 \n\t"
                               "vmulpd %%ymm14, %%ymm8, %%ymm8 \n\t"
                               "vmulpd %%ymm12, %%ymm9, %%ymm9 \n\t"
                               "vmulpd %%ymm13, %%ymm10, %%ymm10 \n\t"
                               "vmulpd %%ymm14, %%ymm11, %%ymm11 \n\t"
                               "vaddpd %%ymm6, %%ymm0, %%ymm0 \n\t"
                               "vaddpd %%ymm7, %%ymm1, %%ymm1 \n\t"
                               "vaddpd %%ymm8, %%ymm2, %%ymm2 \n\t"
                               "vaddsubpd %%ymm9, %%ymm3, %%ymm3 \n\t"
                               "vaddsubpd %%ymm10, %%ymm4, %%ymm4 \n\t"
                               "vaddsubpd %%ymm11, %%ymm5, %%ymm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5",
                               "xmm6", "xmm7", "xmm8",
                               "xmm9", "xmm10", "xmm11");

         __asm__ __volatile__ ("vmovapd %0, %%ymm6 \n\t"
                               "vmovapd %2, %%ymm7 \n\t"
                               "vmovapd %4, %%ymm8"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm6", "xmm7", "xmm8");

         __asm__ __volatile__ ("vmovapd %0, %%ymm12 \n\t"
                               "vmovapd %2, %%ymm13 \n\t"
                               "vmovapd %4, %%ymm14"
                               :
                               :
                               "m" ((*s).c3.c1),
                               "m" ((*s).c3.c2),
                               "m" ((*s).c3.c3),
                               "m" ((*s).c4.c1),
                               "m" ((*s).c4.c2),
                               "m" ((*s).c4.c3)
                               :
                               "xmm12", "xmm13", "xmm14");

         __asm__ __volatile__ ("vpermilpd $0x5, %%ymm6, %%ymm9 \n\t"
                               "vpermilpd $0x5, %%ymm7, %%ymm10 \n\t"
                               "vpermilpd $0x5, %%ymm8, %%ymm11 \n\t"
                               "vmulpd %%ymm12, %%ymm6, %%ymm6 \n\t"
                               "vmulpd %%ymm13, %%ymm7, %%ymm7 \n\t"
                               "vmulpd %%ymm14, %%ymm8, %%ymm8 \n\t"
                               "vmulpd %%ymm12, %%ymm9, %%ymm9 \n\t"
                               "vmulpd %%ymm13, %%ymm10, %%ymm10 \n\t"
                               "vmulpd %%ymm14, %%ymm11, %%ymm11 \n\t"
                               "vsubpd %%ymm6, %%ymm0, %%ymm0 \n\t"
                               "vsubpd %%ymm7, %%ymm1, %%ymm1 \n\t"
                               "vsubpd %%ymm8, %%ymm2, %%ymm2 \n\t"
                               "vaddsubpd %%ymm9, %%ymm3, %%ymm3 \n\t"
                               "vaddsubpd %%ymm10, %%ymm4, %%ymm4 \n\t"
                               "vaddsubpd %%ymm11, %%ymm5, %%ymm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5",
                               "xmm6", "xmm7", "xmm8",
                               "xmm9", "xmm10", "xmm11");

         r+=1;
      }

      __asm__ __volatile__ ("vaddpd %%ymm0, %%ymm2, %%ymm2 \n\t"
                            "vaddpd %%ymm3, %%ymm5, %%ymm5 \n\t"
                            "vaddpd %%ymm1, %%ymm2, %%ymm2 \n\t"
                            "vaddpd %%ymm4, %%ymm5, %%ymm5 \n\t"
                            "vhaddpd %%ymm5, %%ymm2, %%ymm0 \n\t"
                            "vextractf128 $0x1, %%ymm0, %%xmm1 \n\t"
                            "vaddpd %%xmm0, %%xmm1, %%xmm2 \n\t"
                            "vmovapd %%xmm2, %0 \n\t"
                            "vzeroupper"
                            :
                            "=m" (smz)
                            :
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm5");

      add_to_hsum(isz,(double*)(&smz));
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isz,(double*)(&smz));
   else
      local_hsum(isz,(double*)(&smz));

   return smz;
}


double norm_square_dble(int vol,int icom,spinor_dble *s)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isx);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("vxorpd %%ymm6, %%ymm6, %%ymm6 \n\t"
                            "vxorpd %%ymm7, %%ymm7, %%ymm7 \n\t"
                            "vxorpd %%ymm8, %%ymm8, %%ymm8 \n\t"
                            "vxorpd %%ymm9, %%ymm9, %%ymm9 \n\t"
                            "vxorpd %%ymm10, %%ymm10, %%ymm10 \n\t"
                            "vxorpd %%ymm11, %%ymm11, %%ymm11"
                            :
                            :
                            :
                            "xmm6", "xmm7", "xmm8",
                            "xmm9", "xmm10", "xmm11");

      for (;s<smb;s++)
      {
         _avx_spinor_load_dble(*s);

         __asm__ __volatile__ ("vmulpd %%ymm0, %%ymm0, %%ymm0 \n\t"
                               "vmulpd %%ymm1, %%ymm1, %%ymm1 \n\t"
                               "vmulpd %%ymm2, %%ymm2, %%ymm2 \n\t"
                               "vmulpd %%ymm3, %%ymm3, %%ymm3 \n\t"
                               "vmulpd %%ymm4, %%ymm4, %%ymm4 \n\t"
                               "vmulpd %%ymm5, %%ymm5, %%ymm5 \n\t"
                               "vaddpd %%ymm0, %%ymm6, %%ymm6 \n\t"
                               "vaddpd %%ymm1, %%ymm7, %%ymm7 \n\t"
                               "vaddpd %%ymm2, %%ymm8, %%ymm8 \n\t"
                               "vaddpd %%ymm3, %%ymm9, %%ymm9 \n\t"
                               "vaddpd %%ymm4, %%ymm10, %%ymm10 \n\t"
                               "vaddpd %%ymm5, %%ymm11, %%ymm11"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5",
                               "xmm6", "xmm7", "xmm8",
                               "xmm9", "xmm10", "xmm11");

      }

      __asm__ __volatile__ ("vaddpd %%ymm6, %%ymm7, %%ymm7 \n\t"
                            "vaddpd %%ymm8, %%ymm9, %%ymm9 \n\t"
                            "vaddpd %%ymm10, %%ymm11, %%ymm11 \n\t"
                            "vaddpd %%ymm7, %%ymm9, %%ymm9 \n\t"
                            "vaddpd %%ymm9, %%ymm11, %%ymm11 \n\t"
                            "vextractf128 $0x1, %%ymm11, %%xmm0 \n\t"
                            "vaddpd %%xmm11, %%xmm0, %%xmm1 \n\t"
                            "vhaddpd %%xmm1, %%xmm1, %%xmm2 \n\t"
                            "vmovsd %%xmm2, %0 \n\t"
                            "vzeroupper"
                            :
                            "=m" (smx)
                            :
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm7", "xmm9", "xmm11");

      add_to_hsum(isx,&smx);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isx,&smx);
   else
      local_hsum(isx,&smx);

   return smx;
}

#endif

void mulc_spinor_add_dble(int vol,spinor_dble *s,spinor_dble *r,
                          complex_dble z)
{
   spinor_dble *sm;

   _avx_load_cmplx_up_dble(z);
   sm=s+vol;

   for (;s<sm;s++)
   {
      _avx_spinor_load_dble(*s);
      _avx_mulc_spinor_add_dble(*r);
      _avx_spinor_store_dble(*s);

      r+=1;
   }

   _avx_zeroupper();
}


void mulr_spinor_add_dble(int vol,spinor_dble *s,spinor_dble *r,
                          double c)
{
   spinor_dble *sm;

   _avx_load_real_up_dble(c);
   sm=s+vol;

   for (;s<sm;s++)
   {
      _avx_spinor_load_dble(*s);
      _avx_mulr_spinor_add_dble(*r);
      _avx_spinor_store_dble(*s);

      r+=1;
   }

   _avx_zeroupper();
}


void combine_spinor_dble(int vol,spinor_dble *s,spinor_dble *r,
                         double cs,double cr)
{
   spinor_dble *sm;

   _avx_load_real_dble(cs);
   _avx_load_real_up_dble(cr);
   sm=s+vol;

   for (;s<sm;s++)
   {
      _avx_mulr_spinor_dble(*s);
      _avx_mulr_spinor_add_dble(*r);
      _avx_spinor_store_dble(*s);

      r+=1;
   }

   _avx_zeroupper();
}


void scale_dble(int vol,double c,spinor_dble *s)
{
   spinor_dble *sm;

   _avx_load_real_dble(c);
   sm=s+vol;

   for (;s<sm;s++)
   {
      _avx_mulr_spinor_dble(*s);
      _avx_spinor_store_dble(*s);
   }

   _avx_zeroupper();
}


void rotate_dble(int vol,int n,spinor_dble **ppk,complex_dble *v)
{
   int k,j,ix;
   complex_dble *z;
   spinor_dble *pk,*pj;

   if (n>nrot)
      alloc_wrotate(n);

   for (ix=0;ix<vol;ix++)
   {
      for (k=0;k<n;k++)
      {
         pj=ppk[0]+ix;
         z=v+k;

         _avx_load_cmplx_dble(*z);
         _avx_mulc_spinor_dble(*pj);

         for (j=1;j<n;j++)
         {
            pj=ppk[j]+ix;
            z+=n;
            _avx_load_cmplx_up_dble(*z);
            _avx_mulc_spinor_add_dble(*pj);
         }

         pk=psi+k;
         _avx_spinor_store_dble(*pk);
      }

      for (k=0;k<n;k++)
      {
         pk=psi+k;
         pj=ppk[k]+ix;

         _avx_spinor_load_dble(*pk);
         _avx_spinor_store_dble(*pj);
      }
   }

   _avx_zeroupper();
}


void mulg5_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;

   __asm__ __volatile__ ("vxorpd %%ymm3, %%ymm3, %%ymm3 \n\t"
                         "vxorpd %%ymm4, %%ymm4, %%ymm4 \n\t"
                         "vxorpd %%ymm5, %%ymm5, %%ymm5"
                         :
                         :
                         :
                         "xmm3", "xmm4", "xmm5");

   sm=s+vol;

   for (;s<sm;s++)
   {
      __asm__ __volatile__ ("vmovapd %0, %%ymm0 \n\t"
                            "vmovapd %2, %%ymm1 \n\t"
                            "vmovapd %4, %%ymm2 \n\t"
                            "vsubpd %%ymm0, %%ymm3, %%ymm0 \n\t"
                            "vsubpd %%ymm1, %%ymm4, %%ymm1 \n\t"
                            "vsubpd %%ymm2, %%ymm5, %%ymm2"
                            :
                            :
                            "m" ((*s).c3.c1),
                            "m" ((*s).c3.c2),
                            "m" ((*s).c3.c3),
                            "m" ((*s).c4.c1),
                            "m" ((*s).c4.c2),
                            "m" ((*s).c4.c3)
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("vmovapd %%ymm0, %0 \n\t"
                            "vmovapd %%ymm1, %2 \n\t"
                            "vmovapd %%ymm2, %4"
                            :
                            "=m" ((*s).c3.c1),
                            "=m" ((*s).c3.c2),
                            "=m" ((*s).c3.c3),
                            "=m" ((*s).c4.c1),
                            "=m" ((*s).c4.c2),
                            "=m" ((*s).c4.c3));
   }

   _avx_zeroupper();
}


void mulmg5_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;

   __asm__ __volatile__ ("vxorpd %%ymm3, %%ymm3, %%ymm3 \n\t"
                         "vxorpd %%ymm4, %%ymm4, %%ymm4 \n\t"
                         "vxorpd %%ymm5, %%ymm5, %%ymm5"
                         :
                         :
                         :
                         "xmm3", "xmm4", "xmm5");

   sm=s+vol;

   for (;s<sm;s++)
   {
      __asm__ __volatile__ ("vmovapd %0, %%ymm0 \n\t"
                            "vmovapd %2, %%ymm1 \n\t"
                            "vmovapd %4, %%ymm2 \n\t"
                            "vsubpd %%ymm0, %%ymm3, %%ymm0 \n\t"
                            "vsubpd %%ymm1, %%ymm4, %%ymm1 \n\t"
                            "vsubpd %%ymm2, %%ymm5, %%ymm2"
                            :
                            :
                            "m" ((*s).c1.c1),
                            "m" ((*s).c1.c2),
                            "m" ((*s).c1.c3),
                            "m" ((*s).c2.c1),
                            "m" ((*s).c2.c2),
                            "m" ((*s).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("vmovapd %%ymm0, %0 \n\t"
                            "vmovapd %%ymm1, %2 \n\t"
                            "vmovapd %%ymm2, %4"
                            :
                            "=m" ((*s).c1.c1),
                            "=m" ((*s).c1.c2),
                            "=m" ((*s).c1.c3),
                            "=m" ((*s).c2.c1),
                            "=m" ((*s).c2.c2),
                            "=m" ((*s).c2.c3));
   }

   _avx_zeroupper();
}

#elif (defined x64)
#include "sse2.h"

complex_dble spinor_prod_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isz);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("xorpd %%xmm6, %%xmm6 \n\t"
                            "xorpd %%xmm7, %%xmm7 \n\t"
                            "xorpd %%xmm8, %%xmm8 \n\t"
                            "xorpd %%xmm9, %%xmm9 \n\t"
                            "xorpd %%xmm10, %%xmm10 \n\t"
                            "xorpd %%xmm11, %%xmm11"
                            :
                            :
                            :
                            "xmm6", "xmm7", "xmm8",
                            "xmm9", "xmm10", "xmm11");

      for (;s<smb;s++)
      {
         _sse_load_dble((*s).c1);
         _sse_load_up_dble((*s).c2);

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm6 \n\t"
                               "addpd %%xmm3, %%xmm7 \n\t"
                               "addpd %%xmm5, %%xmm8"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm6", "xmm7", "xmm8");

         _sse_load_dble((*s).c3);
         _sse_load_up_dble((*s).c4);

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm6 \n\t"
                               "addpd %%xmm3, %%xmm7 \n\t"
                               "addpd %%xmm5, %%xmm8"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm6", "xmm7", "xmm8");

         _sse_load_dble((*s).c1);
         _sse_load_up_dble((*s).c2);

         __asm__ __volatile__ ("shufpd $0x1, %%xmm0, %%xmm0 \n\t"
                               "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                               "shufpd $0x1, %%xmm2, %%xmm2 \n\t"
                               "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                               "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                               "shufpd $0x1, %%xmm5, %%xmm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm9 \n\t"
                               "addpd %%xmm3, %%xmm10 \n\t"
                               "addpd %%xmm5, %%xmm11"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm9", "xmm10", "xmm11");

         _sse_load_dble((*s).c3);
         _sse_load_up_dble((*s).c4);

         __asm__ __volatile__ ("shufpd $0x1, %%xmm0, %%xmm0 \n\t"
                               "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                               "shufpd $0x1, %%xmm2, %%xmm2 \n\t"
                               "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                               "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                               "shufpd $0x1, %%xmm5, %%xmm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm9 \n\t"
                               "addpd %%xmm3, %%xmm10 \n\t"
                               "addpd %%xmm5, %%xmm11"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm9", "xmm10", "xmm11");

         r+=1;
      }

      __asm__ __volatile__ ("addpd %%xmm6, %%xmm8 \n\t"
                            "addpd %%xmm9, %%xmm11 \n\t"
                            "addpd %%xmm7, %%xmm8 \n\t"
                            "addpd %%xmm10, %%xmm11 \n\t"
                            "haddpd %%xmm8, %%xmm8 \n\t"
                            "hsubpd %%xmm11, %%xmm11 \n\t"
                            "movsd %%xmm8, %0 \n\t"
                            "movsd %%xmm11, %1"
                            :
                            "=m" (smz.re),
                            "=m" (smz.im)
                            :
                            :
                            "xmm8", "xmm11");

      add_to_hsum(isz,(double*)(&smz));
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isz,(double*)(&smz));
   else
      local_hsum(isz,(double*)(&smz));

   smz.im=-smz.im;

   return smz;
}


double spinor_prod_re_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isx);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("xorpd %%xmm6, %%xmm6 \n\t"
                            "xorpd %%xmm7, %%xmm7 \n\t"
                            "xorpd %%xmm8, %%xmm8"
                            :
                            :
                            :
                            "xmm6", "xmm7", "xmm8");

      for (;s<smb;s++)
      {
         _sse_load_dble((*s).c1);
         _sse_load_up_dble((*s).c2);

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm6 \n\t"
                               "addpd %%xmm3, %%xmm7 \n\t"
                               "addpd %%xmm5, %%xmm8"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm6", "xmm7", "xmm8");

         _sse_load_dble((*s).c3);
         _sse_load_up_dble((*s).c4);

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm6 \n\t"
                               "addpd %%xmm3, %%xmm7 \n\t"
                               "addpd %%xmm5, %%xmm8"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm6", "xmm7", "xmm8");

         r+=1;
      }

      __asm__ __volatile__ ("addpd %%xmm6, %%xmm8 \n\t"
                            "addpd %%xmm7, %%xmm8 \n\t"
                            "haddpd %%xmm8, %%xmm8 \n\t"
                            "movsd %%xmm8, %0"
                            :
                            "=m" (smx)
                            :
                            :
                            "xmm8");

      add_to_hsum(isx,&smx);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isx,&smx);
   else
      local_hsum(isx,&smx);

   return smx;
}


complex_dble spinor_prod5_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isz);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("xorpd %%xmm6, %%xmm6 \n\t"
                            "xorpd %%xmm7, %%xmm7 \n\t"
                            "xorpd %%xmm8, %%xmm8 \n\t"
                            "xorpd %%xmm9, %%xmm9 \n\t"
                            "xorpd %%xmm10, %%xmm10 \n\t"
                            "xorpd %%xmm11, %%xmm11"
                            :
                            :
                            :
                            "xmm6", "xmm7", "xmm8",
                            "xmm9", "xmm10", "xmm11");

      for (;s<smb;s++)
      {
         _sse_load_dble((*s).c1);
         _sse_load_up_dble((*s).c2);

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm6 \n\t"
                               "addpd %%xmm3, %%xmm7 \n\t"
                               "addpd %%xmm5, %%xmm8"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm6", "xmm7", "xmm8");

         _sse_load_dble((*s).c3);
         _sse_load_up_dble((*s).c4);

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "subpd %%xmm1, %%xmm6 \n\t"
                               "subpd %%xmm3, %%xmm7 \n\t"
                               "subpd %%xmm5, %%xmm8"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm6", "xmm7", "xmm8");

         _sse_load_dble((*s).c1);
         _sse_load_up_dble((*s).c2);

         __asm__ __volatile__ ("shufpd $0x1, %%xmm0, %%xmm0 \n\t"
                               "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                               "shufpd $0x1, %%xmm2, %%xmm2 \n\t"
                               "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                               "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                               "shufpd $0x1, %%xmm5, %%xmm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*r).c1.c1),
                               "m" ((*r).c1.c2),
                               "m" ((*r).c1.c3),
                               "m" ((*r).c2.c1),
                               "m" ((*r).c2.c2),
                               "m" ((*r).c2.c3)
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm9 \n\t"
                               "addpd %%xmm3, %%xmm10 \n\t"
                               "addpd %%xmm5, %%xmm11"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm9", "xmm10", "xmm11");

         _sse_load_dble((*s).c3);
         _sse_load_up_dble((*s).c4);

         __asm__ __volatile__ ("shufpd $0x1, %%xmm0, %%xmm0 \n\t"
                               "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                               "shufpd $0x1, %%xmm2, %%xmm2 \n\t"
                               "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                               "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                               "shufpd $0x1, %%xmm5, %%xmm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*r).c3.c1),
                               "m" ((*r).c3.c2),
                               "m" ((*r).c3.c3),
                               "m" ((*r).c4.c1),
                               "m" ((*r).c4.c2),
                               "m" ((*r).c4.c3)
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "subpd %%xmm1, %%xmm9 \n\t"
                               "subpd %%xmm3, %%xmm10 \n\t"
                               "subpd %%xmm5, %%xmm11"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm9", "xmm10", "xmm11");

         r+=1;
      }

      __asm__ __volatile__ ("addpd %%xmm6, %%xmm8 \n\t"
                            "addpd %%xmm9, %%xmm11 \n\t"
                            "addpd %%xmm7, %%xmm8 \n\t"
                            "addpd %%xmm10, %%xmm11 \n\t"
                            "haddpd %%xmm8, %%xmm8 \n\t"
                            "hsubpd %%xmm11, %%xmm11 \n\t"
                            "movsd %%xmm8, %0 \n\t"
                            "movsd %%xmm11, %1"
                            :
                            "=m" (smz.re),
                            "=m" (smz.im)
                            :
                            :
                            "xmm8", "xmm11");

      add_to_hsum(isz,(double*)(&smz));
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isz,(double*)(&smz));
   else
      local_hsum(isz,(double*)(&smz));

   smz.im=-smz.im;

   return smz;
}


double norm_square_dble(int vol,int icom,spinor_dble *s)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isx);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      __asm__ __volatile__ ("xorpd %%xmm6, %%xmm6 \n\t"
                            "xorpd %%xmm7, %%xmm7 \n\t"
                            "xorpd %%xmm8, %%xmm8"
                            :
                            :
                            :
                            "xmm6", "xmm7", "xmm8");

      for (;s<smb;s++)
      {
         _sse_load_dble((*s).c1);
         _sse_load_up_dble((*s).c2);

         __asm__ __volatile__ ("mulpd %%xmm0, %%xmm0 \n\t"
                               "mulpd %%xmm1, %%xmm1 \n\t"
                               "mulpd %%xmm2, %%xmm2 \n\t"
                               "mulpd %%xmm3, %%xmm3 \n\t"
                               "mulpd %%xmm4, %%xmm4 \n\t"
                               "mulpd %%xmm5, %%xmm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm6 \n\t"
                               "addpd %%xmm3, %%xmm7 \n\t"
                               "addpd %%xmm5, %%xmm8"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm6", "xmm7", "xmm8");

         _sse_load_dble((*s).c3);
         _sse_load_up_dble((*s).c4);

         __asm__ __volatile__ ("mulpd %%xmm0, %%xmm0 \n\t"
                               "mulpd %%xmm1, %%xmm1 \n\t"
                               "mulpd %%xmm2, %%xmm2 \n\t"
                               "mulpd %%xmm3, %%xmm3 \n\t"
                               "mulpd %%xmm4, %%xmm4 \n\t"
                               "mulpd %%xmm5, %%xmm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2",
                               "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm6 \n\t"
                               "addpd %%xmm3, %%xmm7 \n\t"
                               "addpd %%xmm5, %%xmm8"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5",
                               "xmm6", "xmm7", "xmm8");
      }

      __asm__ __volatile__ ("addpd %%xmm6, %%xmm8 \n\t"
                            "addpd %%xmm7, %%xmm8 \n\t"
                            "haddpd %%xmm8, %%xmm8 \n\t"
                            "movsd %%xmm8, %0"
                            :
                            "=m" (smx)
                            :
                            :
                            "xmm8");

      add_to_hsum(isx,&smx);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isx,&smx);
   else
      local_hsum(isx,&smx);

   return smx;
}


void mulc_spinor_add_dble(int vol,spinor_dble *s,spinor_dble *r,
                          complex_dble z)
{
   spinor_dble *sm;

   _sse_load_cmplx_dble(z);
   sm=s+vol;

   for (;s<sm;s++)
   {
      _sse_load_dble((*s).c1);
      _sse_load_up_dble((*s).c2);
      _sse_mulc_vector_add_dble((*r).c1);
      _sse_mulc_vector_add_up_dble((*r).c2);
      _sse_store_dble((*s).c1);
      _sse_store_up_dble((*s).c2);

      _sse_load_dble((*s).c3);
      _sse_load_up_dble((*s).c4);
      _sse_mulc_vector_add_dble((*r).c3);
      _sse_mulc_vector_add_up_dble((*r).c4);
      _sse_store_dble((*s).c3);
      _sse_store_up_dble((*s).c4);

      r+=1;
   }
}


void mulr_spinor_add_dble(int vol,spinor_dble *s,spinor_dble *r,
                          double c)
{
   spinor_dble *sm;

   _sse_load_real_dble(c);

   sm=s+vol;

   for (;s<sm;s++)
   {
      _sse_load_dble((*r).c1);
      _sse_mulr_vector_add_dble((*s).c1);
      _sse_load_up_dble((*r).c2);
      _sse_mulr_vector_add_up_dble((*s).c2);
      _sse_store_dble((*s).c1);
      _sse_store_up_dble((*s).c2);

      _sse_load_dble((*r).c3);
      _sse_mulr_vector_add_dble((*s).c3);
      _sse_load_up_dble((*r).c4);
      _sse_mulr_vector_add_up_dble((*s).c4);
      _sse_store_dble((*s).c3);
      _sse_store_up_dble((*s).c4);

      r+=1;
   }
}


void combine_spinor_dble(int vol,spinor_dble *s,spinor_dble *r,
                         double cs,double cr)
{
   spinor_dble *sm;

   __asm__ __volatile__ ("movddup %0, %%xmm12 \n\t"
                         "movddup %1, %%xmm14 \n\t"
                         "movapd %%xmm12, %%xmm13 \n\t"
                         "movapd %%xmm14, %%xmm15 \n\t"
                         :
                         :
                         "m" (cs),
                         "m" (cr)
                         :
                         "xmm12", "xmm13", "xmm14", "xmm15");

   sm=s+vol;

   for (;s<sm;s++)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5 \n\t"
                            "mulpd %%xmm12, %%xmm0 \n\t"
                            "mulpd %%xmm13, %%xmm1 \n\t"
                            "mulpd %%xmm12, %%xmm2 \n\t"
                            "mulpd %%xmm13, %%xmm3 \n\t"
                            "mulpd %%xmm12, %%xmm4 \n\t"
                            "mulpd %%xmm13, %%xmm5"
                            :
                            :
                            "m" ((*s).c1.c1),
                            "m" ((*s).c1.c2),
                            "m" ((*s).c1.c3),
                            "m" ((*s).c2.c1),
                            "m" ((*s).c2.c2),
                            "m" ((*s).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                            "movapd %1, %%xmm7 \n\t"
                            "movapd %2, %%xmm8 \n\t"
                            "movapd %3, %%xmm9 \n\t"
                            "movapd %4, %%xmm10 \n\t"
                            "movapd %5, %%xmm11 \n\t"
                            "mulpd %%xmm14, %%xmm6 \n\t"
                            "mulpd %%xmm15, %%xmm7 \n\t"
                            "mulpd %%xmm14, %%xmm8 \n\t"
                            "mulpd %%xmm15, %%xmm9 \n\t"
                            "mulpd %%xmm14, %%xmm10 \n\t"
                            "mulpd %%xmm15, %%xmm11"
                            :
                            :
                            "m" ((*r).c1.c1),
                            "m" ((*r).c1.c2),
                            "m" ((*r).c1.c3),
                            "m" ((*r).c2.c1),
                            "m" ((*r).c2.c2),
                            "m" ((*r).c2.c3)
                            :
                            "xmm6", "xmm7", "xmm8", "xmm9",
                            "xmm10", "xmm11");

      __asm__ __volatile__ ("addpd %%xmm6, %%xmm0 \n\t"
                            "addpd %%xmm7, %%xmm1 \n\t"
                            "addpd %%xmm8, %%xmm2 \n\t"
                            "addpd %%xmm9, %%xmm3 \n\t"
                            "addpd %%xmm10, %%xmm4 \n\t"
                            "addpd %%xmm11, %%xmm5 \n\t"
                            "movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"
                            :
                            "=m" ((*s).c1.c1),
                            "=m" ((*s).c1.c2),
                            "=m" ((*s).c1.c3),
                            "=m" ((*s).c2.c1),
                            "=m" ((*s).c2.c2),
                            "=m" ((*s).c2.c3)
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5 \n\t"
                            "mulpd %%xmm12, %%xmm0 \n\t"
                            "mulpd %%xmm13, %%xmm1 \n\t"
                            "mulpd %%xmm12, %%xmm2 \n\t"
                            "mulpd %%xmm13, %%xmm3 \n\t"
                            "mulpd %%xmm12, %%xmm4 \n\t"
                            "mulpd %%xmm13, %%xmm5"
                            :
                            :
                            "m" ((*s).c3.c1),
                            "m" ((*s).c3.c2),
                            "m" ((*s).c3.c3),
                            "m" ((*s).c4.c1),
                            "m" ((*s).c4.c2),
                            "m" ((*s).c4.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                            "movapd %1, %%xmm7 \n\t"
                            "movapd %2, %%xmm8 \n\t"
                            "movapd %3, %%xmm9 \n\t"
                            "movapd %4, %%xmm10 \n\t"
                            "movapd %5, %%xmm11 \n\t"
                            "mulpd %%xmm14, %%xmm6 \n\t"
                            "mulpd %%xmm15, %%xmm7 \n\t"
                            "mulpd %%xmm14, %%xmm8 \n\t"
                            "mulpd %%xmm15, %%xmm9 \n\t"
                            "mulpd %%xmm14, %%xmm10 \n\t"
                            "mulpd %%xmm15, %%xmm11"
                            :
                            :
                            "m" ((*r).c3.c1),
                            "m" ((*r).c3.c2),
                            "m" ((*r).c3.c3),
                            "m" ((*r).c4.c1),
                            "m" ((*r).c4.c2),
                            "m" ((*r).c4.c3)
                            :
                            "xmm6", "xmm7", "xmm8", "xmm9",
                            "xmm10", "xmm11");

      __asm__ __volatile__ ("addpd %%xmm6, %%xmm0 \n\t"
                            "addpd %%xmm7, %%xmm1 \n\t"
                            "addpd %%xmm8, %%xmm2 \n\t"
                            "addpd %%xmm9, %%xmm3 \n\t"
                            "addpd %%xmm10, %%xmm4 \n\t"
                            "addpd %%xmm11, %%xmm5 \n\t"
                            "movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"
                            :
                            "=m" ((*s).c3.c1),
                            "=m" ((*s).c3.c2),
                            "=m" ((*s).c3.c3),
                            "=m" ((*s).c4.c1),
                            "=m" ((*s).c4.c2),
                            "=m" ((*s).c4.c3)
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

         r+=1;
   }
}


void scale_dble(int vol,double c,spinor_dble *s)
{
   spinor_dble *sm;

   _sse_load_real_dble(c);

   sm=s+vol;

   for (;s<sm;s++)
   {
      _sse_load_dble((*s).c1);
      _sse_mulr_vector_dble();
      _sse_load_up_dble((*s).c2);
      _sse_mulr_vector_up_dble();
      _sse_store_dble((*s).c1);
      _sse_store_up_dble((*s).c2);

      _sse_load_dble((*s).c3);
      _sse_mulr_vector_dble();
      _sse_load_up_dble((*s).c4);
      _sse_mulr_vector_up_dble();
      _sse_store_dble((*s).c3);
      _sse_store_up_dble((*s).c4);
   }
}


void rotate_dble(int vol,int n,spinor_dble **ppk,complex_dble *v)
{
   int k,j,ix;
   complex_dble *z;
   spinor_dble *pk,*pj;

   if (n>nrot)
      alloc_wrotate(n);

   for (ix=0;ix<vol;ix++)
   {
      for (k=0;k<n;k++)
      {
         pj=ppk[0]+ix;
         z=v+k;

         _sse_load_cmplx_dble(*z);
         _sse_load_dble((*pj).c1);
         _sse_load_up_dble((*pj).c2);
         _sse_mulc_vector_dble();
         _sse_mulc_vector_up_dble();

         for (j=1;j<n;j++)
         {
            pj=ppk[j]+ix;
            z+=n;
            _sse_load_cmplx_dble(*z);
            _sse_mulc_vector_add_dble((*pj).c1);
            _sse_mulc_vector_add_up_dble((*pj).c2);
         }

         pk=psi+k;
         _sse_store_dble((*pk).c1);
         _sse_store_up_dble((*pk).c2);
      }

      for (k=0;k<n;k++)
      {
         pj=ppk[0]+ix;
         z=v+k;

         _sse_load_cmplx_dble(*z);
         _sse_load_dble((*pj).c3);
         _sse_load_up_dble((*pj).c4);
         _sse_mulc_vector_dble();
         _sse_mulc_vector_up_dble();

         for (j=1;j<n;j++)
         {
            pj=ppk[j]+ix;
            z+=n;
            _sse_load_cmplx_dble(*z);
            _sse_mulc_vector_add_dble((*pj).c3);
            _sse_mulc_vector_add_up_dble((*pj).c4);
         }

         pk=psi+k;
         _sse_store_dble((*pk).c3);
         _sse_store_up_dble((*pk).c4);
      }

      for (k=0;k<n;k++)
      {
         pk=psi+k;
         pj=ppk[k]+ix;

         _sse_load_dble((*pk).c1);
         _sse_load_up_dble((*pk).c2);
         _sse_store_dble((*pj).c1);
         _sse_store_up_dble((*pj).c2);

         _sse_load_dble((*pk).c3);
         _sse_load_up_dble((*pk).c4);
         _sse_store_dble((*pj).c3);
         _sse_store_up_dble((*pj).c4);
      }
   }
}


void mulg5_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;


   __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                         "movapd %%xmm6, %%xmm7"
                         :
                         :
                         "m" (_sse_sgn_dble)
                         :
                         "xmm6", "xmm7");

   sm=s+vol;

   for (;s<sm;s++)
   {
      _sse_load_dble((*s).c3);
      _sse_mulr_vector_dble();
      _sse_load_up_dble((*s).c4);
      _sse_mulr_vector_up_dble();
      _sse_store_dble((*s).c3);
      _sse_store_up_dble((*s).c4);
   }
}


void mulmg5_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;

   __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                         "movapd %%xmm6, %%xmm7"
                         :
                         :
                         "m" (_sse_sgn_dble)
                         :
                         "xmm6", "xmm7");

   sm=s+vol;

   for (;s<sm;s++)
   {
      _sse_load_dble((*s).c1);
      _sse_mulr_vector_dble();
      _sse_load_up_dble((*s).c2);
      _sse_mulr_vector_up_dble();
      _sse_store_dble((*s).c1);
      _sse_store_up_dble((*s).c2);
   }
}

#else

complex_dble spinor_prod_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isz);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      smz.re=0.0;
      smz.im=0.0;

      for (;s<smb;s++)
      {
         smz.re+=(_vector_prod_re((*s).c1,(*r).c1)+
                  _vector_prod_re((*s).c2,(*r).c2)+
                  _vector_prod_re((*s).c3,(*r).c3)+
                  _vector_prod_re((*s).c4,(*r).c4));

         smz.im+=(_vector_prod_im((*s).c1,(*r).c1)+
                  _vector_prod_im((*s).c2,(*r).c2)+
                  _vector_prod_im((*s).c3,(*r).c3)+
                  _vector_prod_im((*s).c4,(*r).c4));

         r+=1;
      }

      add_to_hsum(isz,(double*)(&smz));
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isz,(double*)(&smz));
   else
      local_hsum(isz,(double*)(&smz));

   return smz;
}


double spinor_prod_re_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isx);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      smx=0.0;

      for (;s<smb;s++)
      {
         smx+=(_vector_prod_re((*s).c1,(*r).c1)+
               _vector_prod_re((*s).c2,(*r).c2)+
               _vector_prod_re((*s).c3,(*r).c3)+
               _vector_prod_re((*s).c4,(*r).c4));

         r+=1;
      }

      add_to_hsum(isx,&smx);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isx,&smx);
   else
      local_hsum(isx,&smx);

   return smx;
}


complex_dble spinor_prod5_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isz);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      smz.re=0.0;
      smz.im=0.0;

      for (;s<smb;s++)
      {
         smz.re+=(_vector_prod_re((*s).c1,(*r).c1)+
                  _vector_prod_re((*s).c2,(*r).c2));

         smz.re-=(_vector_prod_re((*s).c3,(*r).c3)+
                  _vector_prod_re((*s).c4,(*r).c4));

         smz.im+=(_vector_prod_im((*s).c1,(*r).c1)+
                  _vector_prod_im((*s).c2,(*r).c2));

         smz.im-=(_vector_prod_im((*s).c3,(*r).c3)+
                  _vector_prod_im((*s).c4,(*r).c4));

         r+=1;
      }

      add_to_hsum(isz,(double*)(&smz));
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isz,(double*)(&smz));
   else
      local_hsum(isz,(double*)(&smz));

   return smz;
}


double norm_square_dble(int vol,int icom,spinor_dble *s)
{
   spinor_dble *sm,*smb;

   if (init==0)
   {
      isx=init_hsum(1);
      isz=init_hsum(2);
      init=1;
   }

   reset_hsum(isx);
   sm=s+vol;

   while (s<sm)
   {
      smb=s+8;
      if (smb>sm)
         smb=sm;

      smx=0.0;

      for (;s<smb;s++)
      {
         smx+=(_vector_prod_re((*s).c1,(*s).c1)+
               _vector_prod_re((*s).c2,(*s).c2)+
               _vector_prod_re((*s).c3,(*s).c3)+
               _vector_prod_re((*s).c4,(*s).c4));
      }

      add_to_hsum(isx,&smx);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(isx,&smx);
   else
      local_hsum(isx,&smx);

   return smx;
}


void mulc_spinor_add_dble(int vol,spinor_dble *s,spinor_dble *r,
                          complex_dble z)
{
   spinor_dble *sm;

   sm=s+vol;

   for (;s<sm;s++)
   {
      _vector_mulc_assign((*s).c1,z,(*r).c1);
      _vector_mulc_assign((*s).c2,z,(*r).c2);
      _vector_mulc_assign((*s).c3,z,(*r).c3);
      _vector_mulc_assign((*s).c4,z,(*r).c4);

      r+=1;
   }
}


void mulr_spinor_add_dble(int vol,spinor_dble *s,spinor_dble *r,
                          double c)
{
   spinor_dble *sm;

   sm=s+vol;

   for (;s<sm;s++)
   {
      _vector_mulr_assign((*s).c1,c,(*r).c1);
      _vector_mulr_assign((*s).c2,c,(*r).c2);
      _vector_mulr_assign((*s).c3,c,(*r).c3);
      _vector_mulr_assign((*s).c4,c,(*r).c4);

      r+=1;
   }
}


void combine_spinor_dble(int vol,spinor_dble *s,spinor_dble *r,
                         double cs,double cr)
{
   spinor_dble *sm;

   sm=s+vol;

   for (;s<sm;s++)
   {
      _vector_combine((*s).c1,(*r).c1,cs,cr);
      _vector_combine((*s).c2,(*r).c2,cs,cr);
      _vector_combine((*s).c3,(*r).c3,cs,cr);
      _vector_combine((*s).c4,(*r).c4,cs,cr);

      r+=1;
   }
}


void scale_dble(int vol,double c,spinor_dble *s)
{
   spinor_dble *sm;

   sm=s+vol;

   for (;s<sm;s++)
   {
      _vector_mul((*s).c1,c,(*s).c1);
      _vector_mul((*s).c2,c,(*s).c2);
      _vector_mul((*s).c3,c,(*s).c3);
      _vector_mul((*s).c4,c,(*s).c4);
   }
}


void rotate_dble(int vol,int n,spinor_dble **ppk,complex_dble *v)
{
   int k,j,ix;
   complex_dble *z;
   spinor_dble *pk,*pj;

   if (n>nrot)
      alloc_wrotate(n);

   for (ix=0;ix<vol;ix++)
   {
      for (k=0;k<n;k++)
      {
         pk=psi+k;
         pj=ppk[0]+ix;
         z=v+k;

         _vector_mulc((*pk).c1,*z,(*pj).c1);
         _vector_mulc((*pk).c2,*z,(*pj).c2);
         _vector_mulc((*pk).c3,*z,(*pj).c3);
         _vector_mulc((*pk).c4,*z,(*pj).c4);

         for (j=1;j<n;j++)
         {
            pj=ppk[j]+ix;
            z+=n;

            _vector_mulc_assign((*pk).c1,*z,(*pj).c1);
            _vector_mulc_assign((*pk).c2,*z,(*pj).c2);
            _vector_mulc_assign((*pk).c3,*z,(*pj).c3);
            _vector_mulc_assign((*pk).c4,*z,(*pj).c4);
         }
      }

      for (k=0;k<n;k++)
         *(ppk[k]+ix)=psi[k];
   }
}


void mulg5_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;

   sm=s+vol;

   for (;s<sm;s++)
   {
      (*s).c3.c1.re=-(*s).c3.c1.re;
      (*s).c3.c1.im=-(*s).c3.c1.im;
      (*s).c3.c2.re=-(*s).c3.c2.re;
      (*s).c3.c2.im=-(*s).c3.c2.im;
      (*s).c3.c3.re=-(*s).c3.c3.re;
      (*s).c3.c3.im=-(*s).c3.c3.im;
      (*s).c4.c1.re=-(*s).c4.c1.re;
      (*s).c4.c1.im=-(*s).c4.c1.im;
      (*s).c4.c2.re=-(*s).c4.c2.re;
      (*s).c4.c2.im=-(*s).c4.c2.im;
      (*s).c4.c3.re=-(*s).c4.c3.re;
      (*s).c4.c3.im=-(*s).c4.c3.im;
   }
}


void mulmg5_dble(int vol,spinor_dble *s)
{
   spinor_dble *sm;

   sm=s+vol;

   for (;s<sm;s++)
   {
      (*s).c1.c1.re=-(*s).c1.c1.re;
      (*s).c1.c1.im=-(*s).c1.c1.im;
      (*s).c1.c2.re=-(*s).c1.c2.re;
      (*s).c1.c2.im=-(*s).c1.c2.im;
      (*s).c1.c3.re=-(*s).c1.c3.re;
      (*s).c1.c3.im=-(*s).c1.c3.im;
      (*s).c2.c1.re=-(*s).c2.c1.re;
      (*s).c2.c1.im=-(*s).c2.c1.im;
      (*s).c2.c2.re=-(*s).c2.c2.re;
      (*s).c2.c2.im=-(*s).c2.c2.im;
      (*s).c2.c3.re=-(*s).c2.c3.re;
      (*s).c2.c3.im=-(*s).c2.c3.im;
   }
}

#endif

void project_dble(int vol,int icom,spinor_dble *s,spinor_dble *r)
{
   complex_dble z;

   z=spinor_prod_dble(vol,icom,r,s);
   z.re=-z.re;
   z.im=-z.im;
   mulc_spinor_add_dble(vol,s,r,z);
}


double normalize_dble(int vol,int icom,spinor_dble *s)
{
   double r;

   r=norm_square_dble(vol,icom,s);
   r=sqrt(r);

   if (r!=0.0)
      scale_dble(vol,1.0/r,s);
   else
      error_loc(1,1,"normalize_dble [salg_dble.c]",
                "Vector has vanishing norm");

   return r;
}
