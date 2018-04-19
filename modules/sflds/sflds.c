
/*******************************************************************************
*
* File sflds.c
*
* Copyright (C) 2005, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic assignment and initialization programs for single- and double-
* precision spinor fields
*
* The externally accessible functions are
*
*   void set_s2zero(int vol,spinor *s)
*     Sets the single-precision spinor field s to zero.
*
*   void set_sd2zero(int vol,spinor_dble *sd)
*     Sets the double-precision spinor field sd to zero.
*
*   void random_s(int vol,spinor *s,float sigma)
*     Initializes the components of the single-precision field s
*     to (complex) random values z with distribution proportional
*     to exp{-|z|^2/sigma^2}.
*
*   void random_sd(int vol,spinor_dble *sd,double sigma)
*     Initializes the components of the double-precision field sd
*     to (complex) random values z with distribution proportional
*     to exp{-|z|^2/sigma^2}.
*
*   void assign_s2s(int vol,spinor *s,spinor *r)
*     Assigns the single-precision field s to the single-precision
*     field r.
*
*   void assign_s2sd(int vol,spinor *s,spinor_dble *rd)
*     Assigns the single-precision field s to the double-precision
*     field rd.
*
*   void assign_sd2s(int vol,spinor_dble *sd,spinor *r)
*     Assigns the double-precision field sd to the single-precision
*     field r.
*
*   void assign_sd2sd(int vol,spinor_dble *sd,spinor_dble *rd)
*     Assigns the double-precision field sd to the double-precision
*     field rd.
*
*   void diff_s2s(int vol,spinor *s,spinor *r)
*     Assigns the difference s-r of the single-precision fields s and
*     r to r.
*
*   void add_s2sd(int vol,spinor *s,spinor_dble *rd)
*     Adds the single-precision field s to the double-precision field
*     rd.
*
*   void diff_sd2s(int vol,spinor_dble *sd,spinor_dble *rd,spinor *r)
*     Assigns the difference sd-rd of the double-precision fields sd 
*     and rd to the single-precision field r.
* 
* Notes:
*
* All these programs operate on arrays of spinor fields, whose base address
* is passed through the arguments. The length of the arrays is specified by
* the parameter vol. 
*
* Since no communications are performed, all programs in this file can be
* called locally. If SSE instructions are used, the fields must be aligned
* to a 16 byte boundary.
*
*******************************************************************************/

#define SFLDS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "random.h"
#include "sflds.h"

#if (defined x64)
#include "sse2.h"

void set_s2zero(int vol,spinor *s)
{
   spinor *sm;

   __asm__ __volatile__ ("xorps %%xmm0, %%xmm0 \n\t"
                         "xorps %%xmm1, %%xmm1 \n\t"
                         "xorps %%xmm2, %%xmm2"
                         :
                         :
                         :
                         "xmm0", "xmm1", "xmm2");
   
   sm=s+vol;
   
   for (;s<sm;s++)
   {
      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %2 \n\t"
                            "movaps %%xmm2, %4"
                            :
                            "=m" ((*s).c1.c1),
                            "=m" ((*s).c1.c2),
                            "=m" ((*s).c1.c3),
                            "=m" ((*s).c2.c1),
                            "=m" ((*s).c2.c2),
                            "=m" ((*s).c2.c3));

      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %2 \n\t"
                            "movaps %%xmm2, %4"
                            :
                            "=m" ((*s).c3.c1),
                            "=m" ((*s).c3.c2),
                            "=m" ((*s).c3.c3),
                            "=m" ((*s).c4.c1),
                            "=m" ((*s).c4.c2),
                            "=m" ((*s).c4.c3));      
   }
}


void set_sd2zero(int vol,spinor_dble *sd)
{
   spinor_dble *sm;

   __asm__ __volatile__ ("xorpd %%xmm0, %%xmm0 \n\t"
                         "xorpd %%xmm1, %%xmm1 \n\t"
                         "xorpd %%xmm2, %%xmm2"
                         :
                         :
                         :
                         "xmm0", "xmm1", "xmm2");
   
   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm0, %3 \n\t"
                            "movapd %%xmm1, %4 \n\t"
                            "movapd %%xmm2, %5"                            
                            :
                            "=m" ((*sd).c1.c1),
                            "=m" ((*sd).c1.c2),
                            "=m" ((*sd).c1.c3),
                            "=m" ((*sd).c2.c1),
                            "=m" ((*sd).c2.c2),
                            "=m" ((*sd).c2.c3));

      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm0, %3 \n\t"
                            "movapd %%xmm1, %4 \n\t"
                            "movapd %%xmm2, %5"                            
                            :
                            "=m" ((*sd).c3.c1),
                            "=m" ((*sd).c3.c2),
                            "=m" ((*sd).c3.c3),
                            "=m" ((*sd).c4.c1),
                            "=m" ((*sd).c4.c2),
                            "=m" ((*sd).c4.c3));
   }
}


void assign_s2s(int vol,spinor *s,spinor *r)
{
   spinor *rm;

   rm=r+vol;
   
   for (;r<rm;r++)
   {
      _sse_spinor_load(*s);
      s+=4;
      _prefetch_spinor(s);
      s-=3;
      _sse_spinor_store(*r);
   }
}


void assign_s2sd(int vol,spinor *s,spinor_dble *rd)
{
   spinor_dble *rm;

   rm=rd+vol;
   
   for (;rd<rm;rd++)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %2, %%xmm2 \n\t"
                            "movaps %4, %%xmm4"
                            :
                            :
                            "m" ((*s).c1.c1),
                            "m" ((*s).c1.c2),
                            "m" ((*s).c1.c3),
                            "m" ((*s).c2.c1),
                            "m" ((*s).c2.c2),
                            "m" ((*s).c2.c3)
                            :
                            "xmm0", "xmm2", "xmm4");

      __asm__ __volatile__ ("movaps %0, %%xmm6 \n\t"
                            "movaps %2, %%xmm8 \n\t"
                            "movaps %4, %%xmm10"
                            :
                            :
                            "m" ((*s).c3.c1),
                            "m" ((*s).c3.c2),
                            "m" ((*s).c3.c3),
                            "m" ((*s).c4.c1),
                            "m" ((*s).c4.c2),
                            "m" ((*s).c4.c3)
                            :
                            "xmm6", "xmm8", "xmm10");      

      s+=4;
      _prefetch_spinor(s);
      s-=3;
      
      __asm__ __volatile__ ("movhlps %%xmm0, %%xmm1 \n\t"
                            "movhlps %%xmm2, %%xmm3 \n\t"
                            "movhlps %%xmm4, %%xmm5 \n\t"
                            "movhlps %%xmm6, %%xmm7 \n\t"
                            "movhlps %%xmm8, %%xmm9 \n\t"
                            "movhlps %%xmm10, %%xmm11"
                            :
                            :
                            :
                            "xmm1", "xmm3", "xmm5", "xmm7",
                            "xmm9", "xmm10");

      __asm__ __volatile__ ("cvtps2pd %%xmm0, %%xmm0 \n\t"
                            "cvtps2pd %%xmm1, %%xmm1 \n\t"
                            "cvtps2pd %%xmm2, %%xmm2 \n\t"
                            "cvtps2pd %%xmm3, %%xmm3 \n\t"
                            "cvtps2pd %%xmm4, %%xmm4 \n\t"
                            "cvtps2pd %%xmm5, %%xmm5 \n\t"
                            "cvtps2pd %%xmm6, %%xmm6 \n\t"
                            "cvtps2pd %%xmm7, %%xmm7 \n\t"
                            "cvtps2pd %%xmm8, %%xmm8 \n\t"
                            "cvtps2pd %%xmm9, %%xmm9 \n\t"
                            "cvtps2pd %%xmm10, %%xmm10 \n\t"
                            "cvtps2pd %%xmm11, %%xmm11"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5", "xmm6", "xmm7",
                            "xmm8", "xmm9", "xmm10", "xmm11");
      
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"                            
                            :
                            "=m" ((*rd).c1.c1),
                            "=m" ((*rd).c1.c2),
                            "=m" ((*rd).c1.c3),
                            "=m" ((*rd).c2.c1),
                            "=m" ((*rd).c2.c2),
                            "=m" ((*rd).c2.c3));

      __asm__ __volatile__ ("movapd %%xmm6, %0 \n\t"
                            "movapd %%xmm7, %1 \n\t"
                            "movapd %%xmm8, %2 \n\t"
                            "movapd %%xmm9, %3 \n\t"
                            "movapd %%xmm10, %4 \n\t"
                            "movapd %%xmm11, %5"                            
                            :
                            "=m" ((*rd).c3.c1),
                            "=m" ((*rd).c3.c2),
                            "=m" ((*rd).c3.c3),
                            "=m" ((*rd).c4.c1),
                            "=m" ((*rd).c4.c2),
                            "=m" ((*rd).c4.c3));      
   }
}


void assign_sd2s(int vol,spinor_dble *sd,spinor *r)
{
   spinor *rm;

   rm=r+vol;
   
   for (;r<rm;r++)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"                            
                            :
                            :
                            "m" ((*sd).c1.c1),
                            "m" ((*sd).c1.c2),
                            "m" ((*sd).c1.c3),
                            "m" ((*sd).c2.c1),
                            "m" ((*sd).c2.c2),
                            "m" ((*sd).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                            "movapd %1, %%xmm7 \n\t"
                            "movapd %2, %%xmm8 \n\t"
                            "movapd %3, %%xmm9 \n\t"
                            "movapd %4, %%xmm10 \n\t"
                            "movapd %5, %%xmm11"                            
                            :
                            :
                            "m" ((*sd).c3.c1),
                            "m" ((*sd).c3.c2),
                            "m" ((*sd).c3.c3),
                            "m" ((*sd).c4.c1),
                            "m" ((*sd).c4.c2),
                            "m" ((*sd).c4.c3)
                            :
                            "xmm6", "xmm7", "xmm8",
                            "xmm9", "xmm10", "xmm11");      

      sd+=4;
      _prefetch_spinor_dble(sd);
      sd-=3;

      __asm__ __volatile__ ("cvtpd2ps %%xmm0, %%xmm0 \n\t"
                            "cvtpd2ps %%xmm1, %%xmm1 \n\t"
                            "cvtpd2ps %%xmm2, %%xmm2 \n\t"
                            "cvtpd2ps %%xmm3, %%xmm3 \n\t"
                            "cvtpd2ps %%xmm4, %%xmm4 \n\t"
                            "cvtpd2ps %%xmm5, %%xmm5 \n\t"
                            "cvtpd2ps %%xmm6, %%xmm6 \n\t"
                            "cvtpd2ps %%xmm7, %%xmm7 \n\t"
                            "cvtpd2ps %%xmm8, %%xmm8 \n\t"
                            "cvtpd2ps %%xmm9, %%xmm9 \n\t"
                            "cvtpd2ps %%xmm10, %%xmm10 \n\t"
                            "cvtpd2ps %%xmm11, %%xmm11"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5", "xmm6", "xmm7",
                            "xmm8", "xmm9", "xmm10", "xmm11");

      __asm__ __volatile__ ("movlhps %%xmm1, %%xmm0 \n\t"
                            "movlhps %%xmm3, %%xmm2 \n\t"
                            "movlhps %%xmm5, %%xmm4 \n\t"
                            "movlhps %%xmm7, %%xmm6 \n\t"
                            "movlhps %%xmm9, %%xmm8 \n\t"
                            "movlhps %%xmm11, %%xmm10"
                            :
                            :
                            :
                            "xmm0", "xmm2", "xmm4", "xmm6",
                            "xmm8", "xmm10");

      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm2, %2 \n\t"
                            "movaps %%xmm4, %4"
                            :
                            "=m" ((*r).c1.c1),
                            "=m" ((*r).c1.c2),
                            "=m" ((*r).c1.c3),
                            "=m" ((*r).c2.c1),
                            "=m" ((*r).c2.c2),
                            "=m" ((*r).c2.c3));

      __asm__ __volatile__ ("movaps %%xmm6, %0 \n\t"
                            "movaps %%xmm8, %2 \n\t"
                            "movaps %%xmm10, %4"
                            :
                            "=m" ((*r).c3.c1),
                            "=m" ((*r).c3.c2),
                            "=m" ((*r).c3.c3),
                            "=m" ((*r).c4.c1),
                            "=m" ((*r).c4.c2),
                            "=m" ((*r).c4.c3));
   }
}


void assign_sd2sd(int vol,spinor_dble *sd,spinor_dble *rd)
{
   spinor_dble *rm;

   rm=rd+vol;
   
   for (;rd<rm;rd++)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"                            
                            :
                            :
                            "m" ((*sd).c1.c1),
                            "m" ((*sd).c1.c2),
                            "m" ((*sd).c1.c3),
                            "m" ((*sd).c2.c1),
                            "m" ((*sd).c2.c2),
                            "m" ((*sd).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                            "movapd %1, %%xmm7 \n\t"
                            "movapd %2, %%xmm8 \n\t"
                            "movapd %3, %%xmm9 \n\t"
                            "movapd %4, %%xmm10 \n\t"
                            "movapd %5, %%xmm11"                            
                            :
                            :
                            "m" ((*sd).c3.c1),
                            "m" ((*sd).c3.c2),
                            "m" ((*sd).c3.c3),
                            "m" ((*sd).c4.c1),
                            "m" ((*sd).c4.c2),
                            "m" ((*sd).c4.c3)
                            :
                            "xmm6", "xmm7", "xmm8",
                            "xmm9", "xmm10", "xmm11");      

      sd+=4;
      _prefetch_spinor_dble(sd);
      sd-=3;
      
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"                            
                            :
                            "=m" ((*rd).c1.c1),
                            "=m" ((*rd).c1.c2),
                            "=m" ((*rd).c1.c3),
                            "=m" ((*rd).c2.c1),
                            "=m" ((*rd).c2.c2),
                            "=m" ((*rd).c2.c3));

      __asm__ __volatile__ ("movapd %%xmm6, %0 \n\t"
                            "movapd %%xmm7, %1 \n\t"
                            "movapd %%xmm8, %2 \n\t"
                            "movapd %%xmm9, %3 \n\t"
                            "movapd %%xmm10, %4 \n\t"
                            "movapd %%xmm11, %5"                            
                            :
                            "=m" ((*rd).c3.c1),
                            "=m" ((*rd).c3.c2),
                            "=m" ((*rd).c3.c3),
                            "=m" ((*rd).c4.c1),
                            "=m" ((*rd).c4.c2),
                            "=m" ((*rd).c4.c3));      
   }
}


void diff_s2s(int vol,spinor *s,spinor *r)
{
   spinor *sm;

   sm=s+vol;
   
   while (s<sm)
   {
      _sse_spinor_load(*s);

      s+=4;
      _prefetch_spinor(s);
      s-=3;
      
      __asm__ __volatile__ ("subps %0, %%xmm0 \n\t"
                            "subps %2, %%xmm1 \n\t"
                            "subps %4, %%xmm2"
                            :
                            :
                            "m" ((*r).c1.c1),
                            "m" ((*r).c1.c2),
                            "m" ((*r).c1.c3),
                            "m" ((*r).c2.c1),
                            "m" ((*r).c2.c2),
                            "m" ((*r).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("subps %0, %%xmm3 \n\t"
                            "subps %2, %%xmm4 \n\t"
                            "subps %4, %%xmm5"
                            :
                            :
                            "m" ((*r).c3.c1),
                            "m" ((*r).c3.c2),
                            "m" ((*r).c3.c3),
                            "m" ((*r).c4.c1),
                            "m" ((*r).c4.c2),
                            "m" ((*r).c4.c3)
                            :
                            "xmm3", "xmm4", "xmm5");      

      _sse_spinor_store(*r);
      
      r+=4;
      _prefetch_spinor(r);
      r-=3;
   }
}


void add_s2sd(int vol,spinor *s,spinor_dble *rd)
{
   spinor_dble *rm;

   rm=rd+vol;
   
   for (;rd<rm;rd++)
   {
      rd+=4;
      _prefetch_spinor_dble(rd);
      rd-=4;
      
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %2, %%xmm2 \n\t"
                            "movaps %4, %%xmm4"
                            :
                            :
                            "m" ((*s).c1.c1),
                            "m" ((*s).c1.c2),
                            "m" ((*s).c1.c3),
                            "m" ((*s).c2.c1),
                            "m" ((*s).c2.c2),
                            "m" ((*s).c2.c3)
                            :
                            "xmm0", "xmm2", "xmm4");

      __asm__ __volatile__ ("movaps %0, %%xmm6 \n\t"
                            "movaps %2, %%xmm8 \n\t"
                            "movaps %4, %%xmm10"
                            :
                            :
                            "m" ((*s).c3.c1),
                            "m" ((*s).c3.c2),
                            "m" ((*s).c3.c3),
                            "m" ((*s).c4.c1),
                            "m" ((*s).c4.c2),
                            "m" ((*s).c4.c3)
                            :
                            "xmm6", "xmm8", "xmm10");      

      s+=4;
      _prefetch_spinor(s);
      s-=3;
      
      __asm__ __volatile__ ("movhlps %%xmm0, %%xmm1 \n\t"
                            "movhlps %%xmm2, %%xmm3 \n\t"
                            "movhlps %%xmm4, %%xmm5 \n\t"
                            "movhlps %%xmm6, %%xmm7 \n\t"
                            "movhlps %%xmm8, %%xmm9 \n\t"
                            "movhlps %%xmm10, %%xmm11"
                            :
                            :
                            :
                            "xmm1", "xmm3", "xmm5", "xmm7",
                            "xmm9", "xmm10");

      __asm__ __volatile__ ("cvtps2pd %%xmm0, %%xmm0 \n\t"
                            "cvtps2pd %%xmm1, %%xmm1 \n\t"
                            "cvtps2pd %%xmm2, %%xmm2 \n\t"
                            "cvtps2pd %%xmm3, %%xmm3 \n\t"
                            "cvtps2pd %%xmm4, %%xmm4 \n\t"
                            "cvtps2pd %%xmm5, %%xmm5 \n\t"
                            "cvtps2pd %%xmm6, %%xmm6 \n\t"
                            "cvtps2pd %%xmm7, %%xmm7 \n\t"
                            "cvtps2pd %%xmm8, %%xmm8 \n\t"
                            "cvtps2pd %%xmm9, %%xmm9 \n\t"
                            "cvtps2pd %%xmm10, %%xmm10 \n\t"
                            "cvtps2pd %%xmm11, %%xmm11"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5", "xmm6", "xmm7",
                            "xmm8", "xmm9", "xmm10", "xmm11");

      __asm__ __volatile__ ("addpd %0, %%xmm0 \n\t"
                            "addpd %1, %%xmm1 \n\t"
                            "addpd %2, %%xmm2 \n\t"
                            "addpd %3, %%xmm3 \n\t"
                            "addpd %4, %%xmm4 \n\t"
                            "addpd %5, %%xmm5"                            
                            :
                            :
                            "m" ((*rd).c1.c1),
                            "m" ((*rd).c1.c2),
                            "m" ((*rd).c1.c3),
                            "m" ((*rd).c2.c1),
                            "m" ((*rd).c2.c2),
                            "m" ((*rd).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("addpd %0, %%xmm6 \n\t"
                            "addpd %1, %%xmm7 \n\t"
                            "addpd %2, %%xmm8 \n\t"
                            "addpd %3, %%xmm9 \n\t"
                            "addpd %4, %%xmm10 \n\t"
                            "addpd %5, %%xmm11"                            
                            :
                            :
                            "m" ((*rd).c3.c1),
                            "m" ((*rd).c3.c2),
                            "m" ((*rd).c3.c3),
                            "m" ((*rd).c4.c1),
                            "m" ((*rd).c4.c2),
                            "m" ((*rd).c4.c3)
                            :
                            "xmm6", "xmm7", "xmm8",
                            "xmm9", "xmm10", "xmm11");      
      
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"                            
                            :
                            "=m" ((*rd).c1.c1),
                            "=m" ((*rd).c1.c2),
                            "=m" ((*rd).c1.c3),
                            "=m" ((*rd).c2.c1),
                            "=m" ((*rd).c2.c2),
                            "=m" ((*rd).c2.c3));

      __asm__ __volatile__ ("movapd %%xmm6, %0 \n\t"
                            "movapd %%xmm7, %1 \n\t"
                            "movapd %%xmm8, %2 \n\t"
                            "movapd %%xmm9, %3 \n\t"
                            "movapd %%xmm10, %4 \n\t"
                            "movapd %%xmm11, %5"                            
                            :
                            "=m" ((*rd).c3.c1),
                            "=m" ((*rd).c3.c2),
                            "=m" ((*rd).c3.c3),
                            "=m" ((*rd).c4.c1),
                            "=m" ((*rd).c4.c2),
                            "=m" ((*rd).c4.c3));      
   }
}

void diff_sd2s(int vol,spinor_dble *sd,spinor_dble *rd,spinor *r)
{
   spinor *rm;

   rm=r+vol;
   
   for (;r<rm;r++)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"                            
                            :
                            :
                            "m" ((*sd).c1.c1),
                            "m" ((*sd).c1.c2),
                            "m" ((*sd).c1.c3),
                            "m" ((*sd).c2.c1),
                            "m" ((*sd).c2.c2),
                            "m" ((*sd).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                            "movapd %1, %%xmm7 \n\t"
                            "movapd %2, %%xmm8 \n\t"
                            "movapd %3, %%xmm9 \n\t"
                            "movapd %4, %%xmm10 \n\t"
                            "movapd %5, %%xmm11"                            
                            :
                            :
                            "m" ((*sd).c3.c1),
                            "m" ((*sd).c3.c2),
                            "m" ((*sd).c3.c3),
                            "m" ((*sd).c4.c1),
                            "m" ((*sd).c4.c2),
                            "m" ((*sd).c4.c3)
                            :
                            "xmm6", "xmm7", "xmm8",
                            "xmm9", "xmm10", "xmm11");      

      sd+=4;
      _prefetch_spinor_dble(sd);
      sd-=3;
      
      __asm__ __volatile__ ("subpd %0, %%xmm0 \n\t"
                            "subpd %1, %%xmm1 \n\t"
                            "subpd %2, %%xmm2 \n\t"
                            "subpd %3, %%xmm3 \n\t"
                            "subpd %4, %%xmm4 \n\t"
                            "subpd %5, %%xmm5"                            
                            :
                            :
                            "m" ((*rd).c1.c1),
                            "m" ((*rd).c1.c2),
                            "m" ((*rd).c1.c3),
                            "m" ((*rd).c2.c1),
                            "m" ((*rd).c2.c2),
                            "m" ((*rd).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2",
                            "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("subpd %0, %%xmm6 \n\t"
                            "subpd %1, %%xmm7 \n\t"
                            "subpd %2, %%xmm8 \n\t"
                            "subpd %3, %%xmm9 \n\t"
                            "subpd %4, %%xmm10 \n\t"
                            "subpd %5, %%xmm11"                            
                            :
                            :
                            "m" ((*rd).c3.c1),
                            "m" ((*rd).c3.c2),
                            "m" ((*rd).c3.c3),
                            "m" ((*rd).c4.c1),
                            "m" ((*rd).c4.c2),
                            "m" ((*rd).c4.c3)
                            :
                            "xmm6", "xmm7", "xmm8",
                            "xmm9", "xmm10", "xmm11");      

      rd+=4;
      _prefetch_spinor_dble(rd);
      rd-=3;

      __asm__ __volatile__ ("cvtpd2ps %%xmm0, %%xmm0 \n\t"
                            "cvtpd2ps %%xmm1, %%xmm1 \n\t"
                            "cvtpd2ps %%xmm2, %%xmm2 \n\t"
                            "cvtpd2ps %%xmm3, %%xmm3 \n\t"
                            "cvtpd2ps %%xmm4, %%xmm4 \n\t"
                            "cvtpd2ps %%xmm5, %%xmm5 \n\t"
                            "cvtpd2ps %%xmm6, %%xmm6 \n\t"
                            "cvtpd2ps %%xmm7, %%xmm7 \n\t"
                            "cvtpd2ps %%xmm8, %%xmm8 \n\t"
                            "cvtpd2ps %%xmm9, %%xmm9 \n\t"
                            "cvtpd2ps %%xmm10, %%xmm10 \n\t"
                            "cvtpd2ps %%xmm11, %%xmm11"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5", "xmm6", "xmm7",
                            "xmm8", "xmm9", "xmm10", "xmm11");

      __asm__ __volatile__ ("movlhps %%xmm1, %%xmm0 \n\t"
                            "movlhps %%xmm3, %%xmm2 \n\t"
                            "movlhps %%xmm5, %%xmm4 \n\t"
                            "movlhps %%xmm7, %%xmm6 \n\t"
                            "movlhps %%xmm9, %%xmm8 \n\t"
                            "movlhps %%xmm11, %%xmm10"
                            :
                            :
                            :
                            "xmm0", "xmm2", "xmm4", "xmm6",
                            "xmm8", "xmm10");

      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm2, %2 \n\t"
                            "movaps %%xmm4, %4"
                            :
                            "=m" ((*r).c1.c1),
                            "=m" ((*r).c1.c2),
                            "=m" ((*r).c1.c3),
                            "=m" ((*r).c2.c1),
                            "=m" ((*r).c2.c2),
                            "=m" ((*r).c2.c3));

      __asm__ __volatile__ ("movaps %%xmm6, %0 \n\t"
                            "movaps %%xmm8, %2 \n\t"
                            "movaps %%xmm10, %4"
                            :
                            "=m" ((*r).c3.c1),
                            "=m" ((*r).c3.c2),
                            "=m" ((*r).c3.c3),
                            "=m" ((*r).c4.c1),
                            "=m" ((*r).c4.c2),
                            "=m" ((*r).c4.c3));
   }
}

#else

static const spinor s0={{{0.0f}}};
static const spinor_dble sd0={{{0.0}}};


void set_s2zero(int vol,spinor *s)
{
   spinor *sm;

   sm=s+vol;

   for (;s<sm;s++)
      (*s)=s0;
}


void set_sd2zero(int vol,spinor_dble *sd)
{
   spinor_dble *sm;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
      (*sd)=sd0;
}


void assign_s2s(int vol,spinor *s,spinor *r)
{
   spinor *sm;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      (*r)=(*s);
      r+=1;
   }
}


void assign_s2sd(int vol,spinor *s,spinor_dble *rd)
{
   spinor *sm;

   sm=s+vol;
   
   for (;s<sm;s++)
   {   
      (*rd).c1.c1.re=(double)((*s).c1.c1.re);
      (*rd).c1.c1.im=(double)((*s).c1.c1.im);
      (*rd).c1.c2.re=(double)((*s).c1.c2.re);
      (*rd).c1.c2.im=(double)((*s).c1.c2.im);
      (*rd).c1.c3.re=(double)((*s).c1.c3.re);
      (*rd).c1.c3.im=(double)((*s).c1.c3.im);

      (*rd).c2.c1.re=(double)((*s).c2.c1.re);
      (*rd).c2.c1.im=(double)((*s).c2.c1.im);
      (*rd).c2.c2.re=(double)((*s).c2.c2.re);
      (*rd).c2.c2.im=(double)((*s).c2.c2.im);
      (*rd).c2.c3.re=(double)((*s).c2.c3.re);
      (*rd).c2.c3.im=(double)((*s).c2.c3.im);

      (*rd).c3.c1.re=(double)((*s).c3.c1.re);
      (*rd).c3.c1.im=(double)((*s).c3.c1.im);
      (*rd).c3.c2.re=(double)((*s).c3.c2.re);
      (*rd).c3.c2.im=(double)((*s).c3.c2.im);
      (*rd).c3.c3.re=(double)((*s).c3.c3.re);
      (*rd).c3.c3.im=(double)((*s).c3.c3.im);

      (*rd).c4.c1.re=(double)((*s).c4.c1.re);
      (*rd).c4.c1.im=(double)((*s).c4.c1.im);
      (*rd).c4.c2.re=(double)((*s).c4.c2.re);
      (*rd).c4.c2.im=(double)((*s).c4.c2.im);
      (*rd).c4.c3.re=(double)((*s).c4.c3.re);
      (*rd).c4.c3.im=(double)((*s).c4.c3.im);      

      rd+=1;
   }
}


void assign_sd2s(int vol,spinor_dble *sd,spinor *r)
{
   spinor_dble *sm;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {   
      (*r).c1.c1.re=(float)((*sd).c1.c1.re);
      (*r).c1.c1.im=(float)((*sd).c1.c1.im);
      (*r).c1.c2.re=(float)((*sd).c1.c2.re);
      (*r).c1.c2.im=(float)((*sd).c1.c2.im);
      (*r).c1.c3.re=(float)((*sd).c1.c3.re);
      (*r).c1.c3.im=(float)((*sd).c1.c3.im);

      (*r).c2.c1.re=(float)((*sd).c2.c1.re);
      (*r).c2.c1.im=(float)((*sd).c2.c1.im);
      (*r).c2.c2.re=(float)((*sd).c2.c2.re);
      (*r).c2.c2.im=(float)((*sd).c2.c2.im);
      (*r).c2.c3.re=(float)((*sd).c2.c3.re);
      (*r).c2.c3.im=(float)((*sd).c2.c3.im);

      (*r).c3.c1.re=(float)((*sd).c3.c1.re);
      (*r).c3.c1.im=(float)((*sd).c3.c1.im);
      (*r).c3.c2.re=(float)((*sd).c3.c2.re);
      (*r).c3.c2.im=(float)((*sd).c3.c2.im);
      (*r).c3.c3.re=(float)((*sd).c3.c3.re);
      (*r).c3.c3.im=(float)((*sd).c3.c3.im);

      (*r).c4.c1.re=(float)((*sd).c4.c1.re);
      (*r).c4.c1.im=(float)((*sd).c4.c1.im);
      (*r).c4.c2.re=(float)((*sd).c4.c2.re);
      (*r).c4.c2.im=(float)((*sd).c4.c2.im);
      (*r).c4.c3.re=(float)((*sd).c4.c3.re);
      (*r).c4.c3.im=(float)((*sd).c4.c3.im);      
      
      r+=1;
   }
}


void assign_sd2sd(int vol,spinor_dble *sd,spinor_dble *rd)
{
   spinor_dble *sm;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      (*rd)=(*sd);
      rd+=1;
   }
}


void diff_s2s(int vol,spinor *s,spinor *r)
{
   spinor *sm;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      _vector_sub((*r).c1,(*s).c1,(*r).c1);
      _vector_sub((*r).c2,(*s).c2,(*r).c2);
      _vector_sub((*r).c3,(*s).c3,(*r).c3);
      _vector_sub((*r).c4,(*s).c4,(*r).c4);

      r+=1;
   }
}


void add_s2sd(int vol,spinor *s,spinor_dble *rd)
{
   spinor *sm;

   sm=s+vol;
   
   for (;s<sm;s++)
   {   
      (*rd).c1.c1.re+=(double)((*s).c1.c1.re);
      (*rd).c1.c1.im+=(double)((*s).c1.c1.im);
      (*rd).c1.c2.re+=(double)((*s).c1.c2.re);
      (*rd).c1.c2.im+=(double)((*s).c1.c2.im);
      (*rd).c1.c3.re+=(double)((*s).c1.c3.re);
      (*rd).c1.c3.im+=(double)((*s).c1.c3.im);

      (*rd).c2.c1.re+=(double)((*s).c2.c1.re);
      (*rd).c2.c1.im+=(double)((*s).c2.c1.im);
      (*rd).c2.c2.re+=(double)((*s).c2.c2.re);
      (*rd).c2.c2.im+=(double)((*s).c2.c2.im);
      (*rd).c2.c3.re+=(double)((*s).c2.c3.re);
      (*rd).c2.c3.im+=(double)((*s).c2.c3.im);

      (*rd).c3.c1.re+=(double)((*s).c3.c1.re);
      (*rd).c3.c1.im+=(double)((*s).c3.c1.im);
      (*rd).c3.c2.re+=(double)((*s).c3.c2.re);
      (*rd).c3.c2.im+=(double)((*s).c3.c2.im);
      (*rd).c3.c3.re+=(double)((*s).c3.c3.re);
      (*rd).c3.c3.im+=(double)((*s).c3.c3.im);

      (*rd).c4.c1.re+=(double)((*s).c4.c1.re);
      (*rd).c4.c1.im+=(double)((*s).c4.c1.im);
      (*rd).c4.c2.re+=(double)((*s).c4.c2.re);
      (*rd).c4.c2.im+=(double)((*s).c4.c2.im);
      (*rd).c4.c3.re+=(double)((*s).c4.c3.re);
      (*rd).c4.c3.im+=(double)((*s).c4.c3.im);      

      rd+=1;
   }
}


void diff_sd2s(int vol,spinor_dble *sd,spinor_dble *rd,spinor *r)
{
   spinor_dble *sm;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {   
      (*r).c1.c1.re=(float)((*sd).c1.c1.re-(*rd).c1.c1.re);
      (*r).c1.c1.im=(float)((*sd).c1.c1.im-(*rd).c1.c1.im);
      (*r).c1.c2.re=(float)((*sd).c1.c2.re-(*rd).c1.c2.re);
      (*r).c1.c2.im=(float)((*sd).c1.c2.im-(*rd).c1.c2.im);
      (*r).c1.c3.re=(float)((*sd).c1.c3.re-(*rd).c1.c3.re);
      (*r).c1.c3.im=(float)((*sd).c1.c3.im-(*rd).c1.c3.im);

      (*r).c2.c1.re=(float)((*sd).c2.c1.re-(*rd).c2.c1.re);
      (*r).c2.c1.im=(float)((*sd).c2.c1.im-(*rd).c2.c1.im);
      (*r).c2.c2.re=(float)((*sd).c2.c2.re-(*rd).c2.c2.re);
      (*r).c2.c2.im=(float)((*sd).c2.c2.im-(*rd).c2.c2.im);
      (*r).c2.c3.re=(float)((*sd).c2.c3.re-(*rd).c2.c3.re);
      (*r).c2.c3.im=(float)((*sd).c2.c3.im-(*rd).c2.c3.im);

      (*r).c3.c1.re=(float)((*sd).c3.c1.re-(*rd).c3.c1.re);
      (*r).c3.c1.im=(float)((*sd).c3.c1.im-(*rd).c3.c1.im);
      (*r).c3.c2.re=(float)((*sd).c3.c2.re-(*rd).c3.c2.re);
      (*r).c3.c2.im=(float)((*sd).c3.c2.im-(*rd).c3.c2.im);
      (*r).c3.c3.re=(float)((*sd).c3.c3.re-(*rd).c3.c3.re);
      (*r).c3.c3.im=(float)((*sd).c3.c3.im-(*rd).c3.c3.im);

      (*r).c4.c1.re=(float)((*sd).c4.c1.re-(*rd).c4.c1.re);
      (*r).c4.c1.im=(float)((*sd).c4.c1.im-(*rd).c4.c1.im);
      (*r).c4.c2.re=(float)((*sd).c4.c2.re-(*rd).c4.c2.re);
      (*r).c4.c2.im=(float)((*sd).c4.c2.im-(*rd).c4.c2.im);
      (*r).c4.c3.re=(float)((*sd).c4.c3.re-(*rd).c4.c3.re);
      (*r).c4.c3.im=(float)((*sd).c4.c3.im-(*rd).c4.c3.im);      
      
      r+=1;
      rd+=1;
   }
}

#endif

void random_s(int vol,spinor *s,float sigma)
{
   float r[24];
   spinor *sm;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      gauss(r,24);

      (*s).c1.c1.re=sigma*r[ 0];
      (*s).c1.c1.im=sigma*r[ 1];
      (*s).c1.c2.re=sigma*r[ 2];
      (*s).c1.c2.im=sigma*r[ 3];
      (*s).c1.c3.re=sigma*r[ 4];
      (*s).c1.c3.im=sigma*r[ 5];

      (*s).c2.c1.re=sigma*r[ 6];
      (*s).c2.c1.im=sigma*r[ 7];
      (*s).c2.c2.re=sigma*r[ 8];
      (*s).c2.c2.im=sigma*r[ 9];
      (*s).c2.c3.re=sigma*r[10];
      (*s).c2.c3.im=sigma*r[11];

      (*s).c3.c1.re=sigma*r[12];
      (*s).c3.c1.im=sigma*r[13];
      (*s).c3.c2.re=sigma*r[14];
      (*s).c3.c2.im=sigma*r[15];
      (*s).c3.c3.re=sigma*r[16];
      (*s).c3.c3.im=sigma*r[17];      
      
      (*s).c4.c1.re=sigma*r[18];
      (*s).c4.c1.im=sigma*r[19];
      (*s).c4.c2.re=sigma*r[20];
      (*s).c4.c2.im=sigma*r[21];
      (*s).c4.c3.re=sigma*r[22];
      (*s).c4.c3.im=sigma*r[23]; 
   }
}


void random_sd(int vol,spinor_dble *sd,double sigma)
{
   double r[24];
   spinor_dble *sm;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      gauss_dble(r,24);

      (*sd).c1.c1.re=sigma*r[ 0];
      (*sd).c1.c1.im=sigma*r[ 1];
      (*sd).c1.c2.re=sigma*r[ 2];
      (*sd).c1.c2.im=sigma*r[ 3];
      (*sd).c1.c3.re=sigma*r[ 4];
      (*sd).c1.c3.im=sigma*r[ 5];

      (*sd).c2.c1.re=sigma*r[ 6];
      (*sd).c2.c1.im=sigma*r[ 7];
      (*sd).c2.c2.re=sigma*r[ 8];
      (*sd).c2.c2.im=sigma*r[ 9];
      (*sd).c2.c3.re=sigma*r[10];
      (*sd).c2.c3.im=sigma*r[11];

      (*sd).c3.c1.re=sigma*r[12];
      (*sd).c3.c1.im=sigma*r[13];
      (*sd).c3.c2.re=sigma*r[14];
      (*sd).c3.c2.im=sigma*r[15];
      (*sd).c3.c3.re=sigma*r[16];
      (*sd).c3.c3.im=sigma*r[17];      
      
      (*sd).c4.c1.re=sigma*r[18];
      (*sd).c4.c1.im=sigma*r[19];
      (*sd).c4.c2.re=sigma*r[20];
      (*sd).c4.c2.im=sigma*r[21];
      (*sd).c4.c3.re=sigma*r[22];
      (*sd).c4.c3.im=sigma*r[23]; 
   }
}
