
/*******************************************************************************
*
* File cm3x3.c
*
* Copyright (C) 2009, 2010, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Complex 3x3 matrix operations
*
* The externally accessible functions are
*
*   void cm3x3_zero(int vol,su3_dble *u)
*     Sets the elements of the array u[] to zero
*
*   void cm3x3_unity(int vol,su3_dble *u)
*     Sets the elements of the array u[] to the unit matrix
*
*   void cm3x3_assign(int vol,su3_dble *u,su3_dble *v)
*     Assigns the elements of the array u[] to those of the array v[]
*
*   void cm3x3_swap(int vol,su3_dble *u,su3_dble *v)
*     Swaps the elements of the array u[] with those of the array v[]
*
*   void cm3x3_dagger(su3_dble *u,su3_dble *v)
*     Assigns the hermitian conjugate of (*u) to (*v)
*
*   void cm3x3_tr(su3_dble *u,su3_dble *v,complex_dble *tr)
*     Assigns the trace of (*u)*(*v) to (*tr)
*
*   void cm3x3_retr(su3_dble *u,su3_dble *v,double *tr)
*     Assigns the real part of the trace of (*u)*(*v) to (*tr)
*
*   void cm3x3_imtr(su3_dble *u,su3_dble *v,double *tr)
*     Assigns the imaginary part of the trace of (*u)*(*v) to (*tr)
*
*   void cm3x3_add(su3_dble *u,su3_dble *v)
*     Adds (*u) to (*v). The input matrix is unchanged unless u=v
*
*   void cm3x3_mul_add(su3_dble *u,su3_dble *v,su3_dble *w)
*     Adds (*u)*(*v) to (*w) assuming that w!=u. The input matrix (*u)
*     is unchanged and also (*v) unless v=w
*
*   void cm3x3_mulr(double *r,su3_dble *u,su3_dble *v)
*     Assigns (*r)*(*u) to (*v). The input matrix is unchanged
*     unless u=v 
*
*   void cm3x3_mulr_add(double *r,su3_dble *u,su3_dble *v)
*     Adds (*r)*(*u) to (*v). The input matrix is unchanged 
*     unless u=v
*
*   void cm3x3_mulc(complex_dble *c,su3_dble *u,su3_dble *v)
*     Assigns (*c)*(*u) to (*v). The input matrix is unchanged 
*     unless u=v 
*
*   void cm3x3_mulc_add(complex_dble *c,su3_dble *u,su3_dble *v)
*     Adds (*c)*(*u) to (*v). The input matrix is unchanged 
*     unless u=v 
*
*   void cm3x3_lc1(complex_dble *c,su3_dble *u,su3_dble *v)
*     Assigns c[0]+c[1]*(*u) to (*v). The input matrix is unchanged
*     unless u=v
*
*   void cm3x3_lc2(complex_dble *c,su3_dble *u,su3_dble *v)
*     Assigns c[0]+c[1]*u[0]+c[2]*u[1] to (*v) assuming v!=u+1. The
*     input matrix u[1] is unchanged and also u[0] unless u=v
*
* Notes:
*
* The programs in this module do not perform any communications and can be 
* called locally. The parameter vol specifies the number of elements of
* the arrays in the argument list.
*
* If SSE2 instructions are used, it is assumed that the matrices and complex
* coefficients are aligned to a 16 byte boundary.
*
*******************************************************************************/

#define CM3X3_C

#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "su3fcts.h"

#if (defined x64)
#include "sse2.h"

static const double one=1.0;


void cm3x3_zero(int vol,su3_dble *u)
{
   su3_dble *um;
   
   __asm__ __volatile__ ("xorpd %%xmm0, %%xmm0 \n\t"
                         "xorpd %%xmm1, %%xmm1 \n\t"
                         "xorpd %%xmm2, %%xmm2"
                         :
                         :
                         :
                         "xmm0", "xmm1", "xmm2");

   um=u+vol;

   for (;u<um;u++)
   {
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2"                         
                            :
                            "=m" ((*u).c11),
                            "=m" ((*u).c12),
                            "=m" ((*u).c13));

      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2"                         
                            :
                            "=m" ((*u).c21),
                            "=m" ((*u).c22),
                            "=m" ((*u).c23));

      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2"                         
                            :
                            "=m" ((*u).c31),
                            "=m" ((*u).c32),
                            "=m" ((*u).c33));
   }
}


void cm3x3_unity(int vol,su3_dble *u)
{
   su3_dble *um;

   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "xorpd %%xmm1, %%xmm1 \n\t"
                         "xorpd %%xmm2, %%xmm2"
                         :
                         :
                         "m" (one)
                         :
                         "xmm0", "xmm1", "xmm2");

   um=u+vol;

   for (;u<um;u++)
   {
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2"                         
                            :
                            "=m" ((*u).c11),
                            "=m" ((*u).c12),
                            "=m" ((*u).c13));

      __asm__ __volatile__ ("movapd %%xmm1, %0 \n\t"
                            "movapd %%xmm0, %1 \n\t"
                            "movapd %%xmm2, %2"                         
                            :
                            "=m" ((*u).c21),
                            "=m" ((*u).c22),
                            "=m" ((*u).c23));

      __asm__ __volatile__ ("movapd %%xmm1, %0 \n\t"
                            "movapd %%xmm2, %1 \n\t"
                            "movapd %%xmm0, %2"                         
                            :
                            "=m" ((*u).c31),
                            "=m" ((*u).c32),
                            "=m" ((*u).c33));
   }
}


void cm3x3_assign(int vol,su3_dble *u,su3_dble *v)
{
   su3_dble *um;

   um=u+vol;

   for (;u<um;u++)
   {
      __asm__ __volatile__ ("movapd %3, %%xmm0 \n\t"
                            "movapd %4, %%xmm1 \n\t"
                            "movapd %5, %%xmm2 \n\t"
                            "movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2"
                            :
                            "=m" ((*v).c11),
                            "=m" ((*v).c12),
                            "=m" ((*v).c13)
                            :
                            "m" ((*u).c11),
                            "m" ((*u).c12),
                            "m" ((*u).c13)
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5 \n\t"
                            "movapd %%xmm3, %0 \n\t"
                            "movapd %%xmm4, %1 \n\t"
                            "movapd %%xmm5, %2"
                            :
                            "=m" ((*v).c21),
                            "=m" ((*v).c22),
                            "=m" ((*v).c23)
                            :
                            "m" ((*u).c21),
                            "m" ((*u).c22),
                            "m" ((*u).c23)
                            :
                            "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %3, %%xmm6 \n\t"
                            "movapd %4, %%xmm7 \n\t"
                            "movapd %5, %%xmm8 \n\t"
                            "movapd %%xmm6, %0 \n\t"
                            "movapd %%xmm7, %1 \n\t"
                            "movapd %%xmm8, %2"
                            :
                            "=m" ((*v).c31),
                            "=m" ((*v).c32),
                            "=m" ((*v).c33)
                            :
                            "m" ((*u).c31),
                            "m" ((*u).c32),
                            "m" ((*u).c33)
                            :
                            "xmm6", "xmm7", "xmm8");

      v+=1;
   }
}


void cm3x3_swap(int vol,su3_dble *u,su3_dble *v)
{
   su3_dble *um;

   um=u+vol;

   for (;u<um;u++)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"                            
                            :
                            :
                            "m" ((*u).c11),
                            "m" ((*u).c12),
                            "m" ((*u).c13),
                            "m" ((*v).c11),
                            "m" ((*v).c12),
                            "m" ((*v).c13)                            
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"                            
                            :
                            "=m" ((*v).c11),
                            "=m" ((*v).c12),
                            "=m" ((*v).c13),
                            "=m" ((*u).c11),
                            "=m" ((*u).c12),
                            "=m" ((*u).c13));

      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"                            
                            :
                            :
                            "m" ((*u).c21),
                            "m" ((*u).c22),
                            "m" ((*u).c23),
                            "m" ((*v).c21),
                            "m" ((*v).c22),
                            "m" ((*v).c23)                            
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"                            
                            :
                            "=m" ((*v).c21),
                            "=m" ((*v).c22),
                            "=m" ((*v).c23),
                            "=m" ((*u).c21),
                            "=m" ((*u).c22),
                            "=m" ((*u).c23));

      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"                            
                            :
                            :
                            "m" ((*u).c31),
                            "m" ((*u).c32),
                            "m" ((*u).c33),
                            "m" ((*v).c31),
                            "m" ((*v).c32),
                            "m" ((*v).c33)                            
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"                            
                            :
                            "=m" ((*v).c31),
                            "=m" ((*v).c32),
                            "=m" ((*v).c33),
                            "=m" ((*u).c31),
                            "=m" ((*u).c32),
                            "=m" ((*u).c33));

      v+=1;
   }
}


void cm3x3_dagger(su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %0, %%xmm1"
                         :
                         :
                         "m" (_sse_sgn2_dble)
                         :
                         "xmm0", "xmm1");

   __asm__ __volatile__ ("movapd %3, %%xmm2 \n\t"
                         "movapd %4, %%xmm3 \n\t"
                         "movapd %5, %%xmm4 \n\t"
                         "mulpd %%xmm0, %%xmm2 \n\t"
                         "mulpd %%xmm1, %%xmm3 \n\t"
                         "mulpd %%xmm0, %%xmm4 \n\t"                         
                         "movapd %%xmm2, %0 \n\t"
                         "movapd %%xmm3, %1 \n\t"
                         "movapd %%xmm4, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c21),
                         "=m" ((*v).c12)
                         :
                         "m" ((*u).c11),
                         "m" ((*u).c12),
                         "m" ((*u).c21)
                         :
                         "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movapd %3, %%xmm5 \n\t"
                         "movapd %4, %%xmm6 \n\t"
                         "movapd %5, %%xmm7 \n\t"
                         "mulpd %%xmm1, %%xmm5 \n\t"
                         "mulpd %%xmm0, %%xmm6 \n\t"
                         "mulpd %%xmm1, %%xmm7 \n\t"                         
                         "movapd %%xmm5, %0 \n\t"
                         "movapd %%xmm6, %1 \n\t"
                         "movapd %%xmm7, %2"
                         :
                         "=m" ((*v).c22),
                         "=m" ((*v).c31),
                         "=m" ((*v).c13)
                         :
                         "m" ((*u).c22),
                         "m" ((*u).c13),
                         "m" ((*u).c31)
                         :
                         "xmm5", "xmm6", "xmm7");    

   __asm__ __volatile__ ("movapd %3, %%xmm8 \n\t"
                         "movapd %4, %%xmm9 \n\t"
                         "movapd %5, %%xmm10 \n\t"
                         "mulpd %%xmm0, %%xmm8 \n\t"
                         "mulpd %%xmm1, %%xmm9 \n\t"
                         "mulpd %%xmm0, %%xmm10 \n\t"                         
                         "movapd %%xmm8, %0 \n\t"
                         "movapd %%xmm9, %1 \n\t"
                         "movapd %%xmm10, %2"
                         :
                         "=m" ((*v).c33),
                         "=m" ((*v).c32),
                         "=m" ((*v).c23)
                         :
                         "m" ((*u).c33),
                         "m" ((*u).c23),
                         "m" ((*u).c32)
                         :
                         "xmm8", "xmm9", "xmm10");
}


void cm3x3_tr(su3_dble *u,su3_dble *v,complex_dble *tr)
{
   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm3 \n\t"
                         "movddup %4, %%xmm4 \n\t"
                         "movddup %5, %%xmm5"                         
                         :
                         :
                         "m" ((*u).c11.im),
                         "m" ((*u).c12.im),
                         "m" ((*u).c13.im),
                         "m" ((*u).c11.re),
                         "m" ((*u).c12.re),
                         "m" ((*u).c13.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %0, %%xmm3 \n\t"
                         "mulpd %1, %%xmm4 \n\t"
                         "mulpd %2, %%xmm5"                         
                         :
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c21),
                         "m" ((*v).c31)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5");   

   __asm__ __volatile__ ("movddup %0, %%xmm6 \n\t"
                         "movddup %1, %%xmm7 \n\t"
                         "movddup %2, %%xmm8 \n\t"
                         "movddup %3, %%xmm9 \n\t"
                         "movddup %4, %%xmm10 \n\t"
                         "movddup %5, %%xmm11"                         
                         :
                         :
                         "m" ((*u).c21.im),
                         "m" ((*u).c22.im),
                         "m" ((*u).c23.im),
                         "m" ((*u).c21.re),
                         "m" ((*u).c22.re),
                         "m" ((*u).c23.re)
                         :
                         "xmm6", "xmm7", "xmm8", "xmm9",
                         "xmm10", "xmm11");

   __asm__ __volatile__ ("mulpd %0, %%xmm6 \n\t"
                         "mulpd %1, %%xmm7 \n\t"
                         "mulpd %2, %%xmm8 \n\t"
                         "mulpd %0, %%xmm9 \n\t"
                         "mulpd %1, %%xmm10 \n\t"
                         "mulpd %2, %%xmm11 \n\t"
                         "addpd %%xmm6, %%xmm0 \n\t"
                         "addpd %%xmm7, %%xmm1 \n\t"
                         "addpd %%xmm8, %%xmm2 \n\t"
                         "addpd %%xmm9, %%xmm3 \n\t"
                         "addpd %%xmm10, %%xmm4 \n\t"
                         "addpd %%xmm11, %%xmm5"                         
                         :
                         :
                         "m" ((*v).c12),
                         "m" ((*v).c22),
                         "m" ((*v).c32)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7",
                         "xmm8", "xmm9", "xmm10", "xmm11");   

   __asm__ __volatile__ ("movddup %0, %%xmm6 \n\t"
                         "movddup %1, %%xmm7 \n\t"
                         "movddup %2, %%xmm8 \n\t"
                         "movddup %3, %%xmm9 \n\t"
                         "movddup %4, %%xmm10 \n\t"
                         "movddup %5, %%xmm11"                         
                         :
                         :
                         "m" ((*u).c31.im),
                         "m" ((*u).c32.im),
                         "m" ((*u).c33.im),
                         "m" ((*u).c31.re),
                         "m" ((*u).c32.re),
                         "m" ((*u).c33.re)
                         :
                         "xmm6", "xmm7", "xmm8", "xmm9",
                         "xmm10", "xmm11");

   __asm__ __volatile__ ("mulpd %0, %%xmm6 \n\t"
                         "mulpd %1, %%xmm7 \n\t"
                         "mulpd %2, %%xmm8 \n\t"
                         "mulpd %0, %%xmm9 \n\t"
                         "mulpd %1, %%xmm10 \n\t"
                         "mulpd %2, %%xmm11 \n\t"
                         "addpd %%xmm6, %%xmm0 \n\t"
                         "addpd %%xmm7, %%xmm1 \n\t"
                         "addpd %%xmm8, %%xmm2 \n\t"
                         "addpd %%xmm9, %%xmm3 \n\t"
                         "addpd %%xmm10, %%xmm4 \n\t"
                         "addpd %%xmm11, %%xmm5"                         
                         :
                         :
                         "m" ((*v).c13),
                         "m" ((*v).c23),
                         "m" ((*v).c33)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7",
                         "xmm8", "xmm9", "xmm10", "xmm11");   
   
   __asm__ __volatile__ ("addpd %%xmm1, %%xmm0 \n\t"
                         "addpd %%xmm4, %%xmm3 \n\t"
                         "addpd %%xmm2, %%xmm0 \n\t"
                         "addpd %%xmm5, %%xmm3 \n\t"
                         "shufpd $0x1, %%xmm0, %%xmm0 \n\t"
                         "addsubpd %%xmm0, %%xmm3 \n\t"
                         "movapd %%xmm3, %0"
                         :
                         "=m" (*tr)
                         :
                         :
                         "xmm0", "xmm3");
}


void cm3x3_retr(su3_dble *u,su3_dble *v,double *tr)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm0 \n\t"
                         "mulpd %4, %%xmm1 \n\t"
                         "mulpd %5, %%xmm2"
                         :
                         :
                         "m" ((*u).c11),
                         "m" ((*u).c21),
                         "m" ((*u).c31),
                         "m" ((*v).c11),
                         "m" ((*v).c12),
                         "m" ((*v).c13)                         
                         :
                         "xmm0", "xmm1", "xmm2");

   __asm__ __volatile__ ("movapd %0, %%xmm3 \n\t"
                         "movapd %1, %%xmm4 \n\t"
                         "movapd %2, %%xmm5 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "mulpd %5, %%xmm5 \n\t"
                         "addpd %%xmm3, %%xmm0 \n\t"
                         "addpd %%xmm4, %%xmm1 \n\t"
                         "addpd %%xmm5, %%xmm2"
                         :
                         :
                         "m" ((*u).c12),
                         "m" ((*u).c22),
                         "m" ((*u).c32),
                         "m" ((*v).c21),
                         "m" ((*v).c22),
                         "m" ((*v).c23)                         
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5");

   __asm__ __volatile__ ("movapd %0, %%xmm3 \n\t"
                         "movapd %1, %%xmm4 \n\t"
                         "movapd %2, %%xmm5 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "mulpd %5, %%xmm5 \n\t"
                         "addpd %%xmm3, %%xmm0 \n\t"
                         "addpd %%xmm4, %%xmm1 \n\t"
                         "addpd %%xmm5, %%xmm2"
                         :
                         :
                         "m" ((*u).c13),
                         "m" ((*u).c23),
                         "m" ((*u).c33),
                         "m" ((*v).c31),
                         "m" ((*v).c32),
                         "m" ((*v).c33)                         
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5");
   
   __asm__ __volatile__ ("addpd %%xmm1, %%xmm0 \n\t"
                         "addpd %%xmm2, %%xmm0 \n\t"
                         "movapd %%xmm0, %%xmm1 \n\t"
                         "shufpd $0x1, %%xmm0, %%xmm0 \n\t"
                         "subsd %%xmm0, %%xmm1 \n\t"
                         "movsd %%xmm1, %0"
                         :
                         "=m" (*tr)
                         :
                         :
                         "xmm0", "xmm1");   
}


void cm3x3_imtr(su3_dble *u,su3_dble *v,double *tr)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "shufpd $0x1, %%xmm0, %%xmm0 \n\t"
                         "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                         "shufpd $0x1, %%xmm2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm0 \n\t"
                         "mulpd %4, %%xmm1 \n\t"
                         "mulpd %5, %%xmm2"
                         :
                         :
                         "m" ((*u).c11),
                         "m" ((*u).c21),
                         "m" ((*u).c31),
                         "m" ((*v).c11),
                         "m" ((*v).c12),
                         "m" ((*v).c13)                         
                         :
                         "xmm0", "xmm1", "xmm2");

   __asm__ __volatile__ ("movapd %0, %%xmm3 \n\t"
                         "movapd %1, %%xmm4 \n\t"
                         "movapd %2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "shufpd $0x1, %%xmm5, %%xmm5 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "mulpd %5, %%xmm5 \n\t"
                         "addpd %%xmm3, %%xmm0 \n\t"
                         "addpd %%xmm4, %%xmm1 \n\t"
                         "addpd %%xmm5, %%xmm2"
                         :
                         :
                         "m" ((*u).c12),
                         "m" ((*u).c22),
                         "m" ((*u).c32),
                         "m" ((*v).c21),
                         "m" ((*v).c22),
                         "m" ((*v).c23)                         
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5");

   __asm__ __volatile__ ("movapd %0, %%xmm3 \n\t"
                         "movapd %1, %%xmm4 \n\t"
                         "movapd %2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "shufpd $0x1, %%xmm5, %%xmm5 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "mulpd %5, %%xmm5 \n\t"
                         "addpd %%xmm3, %%xmm0 \n\t"
                         "addpd %%xmm4, %%xmm1 \n\t"
                         "addpd %%xmm5, %%xmm2"
                         :
                         :
                         "m" ((*u).c13),
                         "m" ((*u).c23),
                         "m" ((*u).c33),
                         "m" ((*v).c31),
                         "m" ((*v).c32),
                         "m" ((*v).c33)                         
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5");
   
   __asm__ __volatile__ ("addpd %%xmm1, %%xmm0 \n\t"
                         "addpd %%xmm2, %%xmm0 \n\t"
                         "movapd %%xmm0, %%xmm1 \n\t"
                         "shufpd $0x1, %%xmm0, %%xmm0 \n\t"
                         "addsd %%xmm0, %%xmm1 \n\t"
                         "movsd %%xmm1, %0"
                         :
                         "=m" (*tr)
                         :
                         :
                         "xmm0", "xmm1");   
}


void cm3x3_add(su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "addpd %3, %%xmm0 \n\t"
                         "addpd %4, %%xmm1 \n\t"
                         "addpd %5, %%xmm2"
                         :
                         :
                         "m" ((*u).c11),
                         "m" ((*u).c12),
                         "m" ((*u).c13),
                         "m" ((*v).c11),
                         "m" ((*v).c12),
                         "m" ((*v).c13)                         
                         :
                         "xmm0", "xmm1", "xmm2");

   __asm__ __volatile__ ("movapd %0, %%xmm3 \n\t"
                         "movapd %1, %%xmm4 \n\t"
                         "movapd %2, %%xmm5 \n\t"
                         "addpd %3, %%xmm3 \n\t"
                         "addpd %4, %%xmm4 \n\t"
                         "addpd %5, %%xmm5"
                         :
                         :
                         "m" ((*u).c21),
                         "m" ((*u).c22),
                         "m" ((*u).c23),
                         "m" ((*v).c21),
                         "m" ((*v).c22),
                         "m" ((*v).c23)                         
                         :
                         "xmm3", "xmm4", "xmm5");   

   __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                         "movapd %1, %%xmm7 \n\t"
                         "movapd %2, %%xmm8 \n\t"
                         "addpd %3, %%xmm6 \n\t"
                         "addpd %4, %%xmm7 \n\t"
                         "addpd %5, %%xmm8"
                         :
                         :
                         "m" ((*u).c31),
                         "m" ((*u).c32),
                         "m" ((*u).c33),
                         "m" ((*v).c31),
                         "m" ((*v).c32),
                         "m" ((*v).c33)                         
                         :
                         "xmm6", "xmm7", "xmm8");

   __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                         "movapd %%xmm1, %1 \n\t"
                         "movapd %%xmm2, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c12),
                         "=m" ((*v).c13));

   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c21),
                         "=m" ((*v).c22),
                         "=m" ((*v).c23));    

   __asm__ __volatile__ ("movapd %%xmm6, %0 \n\t"
                         "movapd %%xmm7, %1 \n\t"
                         "movapd %%xmm8, %2"
                         :
                         "=m" ((*v).c31),
                         "=m" ((*v).c32),
                         "=m" ((*v).c33));
}


void cm3x3_mul_add(su3_dble *u,su3_dble *v,su3_dble *w)
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

   _sse_su3_multiply_dble(*u); 

   __asm__ __volatile__ ("addpd %3, %%xmm3 \n\t"
                         "addpd %4, %%xmm4 \n\t"
                         "addpd %5, %%xmm5 \n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c11),
                         "=m" ((*w).c21),
                         "=m" ((*w).c31)                         
                         :
                         "m" ((*w).c11),
                         "m" ((*w).c21),
                         "m" ((*w).c31)                         
                         :
                         "xmm3", "xmm4", "xmm5");

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

   _sse_su3_multiply_dble(*u); 

   __asm__ __volatile__ ("addpd %3, %%xmm3 \n\t"
                         "addpd %4, %%xmm4 \n\t"
                         "addpd %5, %%xmm5 \n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c12),
                         "=m" ((*w).c22),
                         "=m" ((*w).c32)                         
                         :
                         "m" ((*w).c12),
                         "m" ((*w).c22),
                         "m" ((*w).c32)                         
                         :
                         "xmm3", "xmm4", "xmm5");

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

   _sse_su3_multiply_dble(*u); 

   __asm__ __volatile__ ("addpd %3, %%xmm3 \n\t"
                         "addpd %4, %%xmm4 \n\t"
                         "addpd %5, %%xmm5 \n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*w).c13),
                         "=m" ((*w).c23),
                         "=m" ((*w).c33)                         
                         :
                         "m" ((*w).c13),
                         "m" ((*w).c23),
                         "m" ((*w).c33)                         
                         :
                         "xmm3", "xmm4", "xmm5");
}


void cm3x3_mulr(double *r,su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %0, %%xmm1 \n\t"
                         "movddup %0, %%xmm2"
                         :
                         :
                         "m" (*r)
                         :
                         "xmm0", "xmm1", "xmm2");   

   __asm__ __volatile__ ("movapd %0, %%xmm3 \n\t"
                         "movapd %1, %%xmm4 \n\t"
                         "movapd %2, %%xmm5 \n\t"
                         "mulpd %%xmm0, %%xmm3 \n\t"
                         "mulpd %%xmm1, %%xmm4 \n\t"
                         "mulpd %%xmm2, %%xmm5"
                         :
                         :
                         "m" ((*u).c11),
                         "m" ((*u).c12),
                         "m" ((*u).c13)
                         :
                         "xmm3", "xmm4", "xmm5");

   __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c12),
                         "=m" ((*v).c13));
   
   __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                         "movapd %1, %%xmm7 \n\t"
                         "movapd %2, %%xmm8 \n\t"
                         "mulpd %%xmm0, %%xmm6 \n\t"
                         "mulpd %%xmm1, %%xmm7 \n\t"
                         "mulpd %%xmm2, %%xmm8"
                         :
                         :
                         "m" ((*u).c21),
                         "m" ((*u).c22),
                         "m" ((*u).c23)
                         :
                         "xmm6", "xmm7", "xmm8");

   __asm__ __volatile__ ("movapd %%xmm6, %0 \n\t"
                         "movapd %%xmm7, %1 \n\t"
                         "movapd %%xmm8, %2"
                         :
                         "=m" ((*v).c21),
                         "=m" ((*v).c22),
                         "=m" ((*v).c23));    

   __asm__ __volatile__ ("movapd %0, %%xmm9 \n\t"
                         "movapd %1, %%xmm10 \n\t"
                         "movapd %2, %%xmm11 \n\t"
                         "mulpd %%xmm0, %%xmm9 \n\t"
                         "mulpd %%xmm1, %%xmm10 \n\t"
                         "mulpd %%xmm2, %%xmm11"
                         :
                         :
                         "m" ((*u).c31),
                         "m" ((*u).c32),
                         "m" ((*u).c33)
                         :
                         "xmm9", "xmm10", "xmm11");

   __asm__ __volatile__ ("movapd %%xmm9, %0 \n\t"
                         "movapd %%xmm10, %1 \n\t"
                         "movapd %%xmm11, %2"
                         :
                         "=m" ((*v).c31),
                         "=m" ((*v).c32),
                         "=m" ((*v).c33));
}


void cm3x3_mulr_add(double *r,su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %0, %%xmm1 \n\t"
                         "movddup %0, %%xmm2"
                         :
                         :
                         "m" (*r)
                         :
                         "xmm0", "xmm1", "xmm2");   

   __asm__ __volatile__ ("movapd %0, %%xmm3 \n\t"
                         "movapd %1, %%xmm4 \n\t"
                         "movapd %2, %%xmm5 \n\t"
                         "mulpd %%xmm0, %%xmm3 \n\t"
                         "mulpd %%xmm1, %%xmm4 \n\t"
                         "mulpd %%xmm2, %%xmm5"
                         :
                         :
                         "m" ((*u).c11),
                         "m" ((*u).c12),
                         "m" ((*u).c13)
                         :
                         "xmm3", "xmm4", "xmm5");

   __asm__ __volatile__ ("addpd %3, %%xmm3 \n\t"
                         "addpd %4, %%xmm4 \n\t"
                         "addpd %5, %%xmm5 \n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm5, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c12),
                         "=m" ((*v).c13)                         
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c12),
                         "m" ((*v).c13)                         
                         :
                         "xmm3", "xmm4", "xmm5");
   
   __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                         "movapd %1, %%xmm7 \n\t"
                         "movapd %2, %%xmm8 \n\t"
                         "mulpd %%xmm0, %%xmm6 \n\t"
                         "mulpd %%xmm1, %%xmm7 \n\t"
                         "mulpd %%xmm2, %%xmm8"
                         :
                         :
                         "m" ((*u).c21),
                         "m" ((*u).c22),
                         "m" ((*u).c23)
                         :
                         "xmm6", "xmm7", "xmm8");

   __asm__ __volatile__ ("addpd %3, %%xmm6 \n\t"
                         "addpd %4, %%xmm7 \n\t"
                         "addpd %5, %%xmm8 \n\t"
                         "movapd %%xmm6, %0 \n\t"
                         "movapd %%xmm7, %1 \n\t"
                         "movapd %%xmm8, %2"
                         :
                         "=m" ((*v).c21),
                         "=m" ((*v).c22),
                         "=m" ((*v).c23)                         
                         :
                         "m" ((*v).c21),
                         "m" ((*v).c22),
                         "m" ((*v).c23)                         
                         :
                         "xmm6", "xmm7", "xmm8");
   
   __asm__ __volatile__ ("movapd %0, %%xmm9 \n\t"
                         "movapd %1, %%xmm10 \n\t"
                         "movapd %2, %%xmm11 \n\t"
                         "mulpd %%xmm0, %%xmm9 \n\t"
                         "mulpd %%xmm1, %%xmm10 \n\t"
                         "mulpd %%xmm2, %%xmm11"
                         :
                         :
                         "m" ((*u).c31),
                         "m" ((*u).c32),
                         "m" ((*u).c33)
                         :
                         "xmm9", "xmm10", "xmm11");

   __asm__ __volatile__ ("addpd %3, %%xmm9 \n\t"
                         "addpd %4, %%xmm10 \n\t"
                         "addpd %5, %%xmm11 \n\t"
                         "movapd %%xmm9, %0 \n\t"
                         "movapd %%xmm10, %1 \n\t"
                         "movapd %%xmm11, %2"
                         :
                         "=m" ((*v).c31),
                         "=m" ((*v).c32),
                         "=m" ((*v).c33)                         
                         :
                         "m" ((*v).c31),
                         "m" ((*v).c32),
                         "m" ((*v).c33)                         
                         :
                         "xmm9", "xmm10", "xmm11");
}


void cm3x3_mulc(complex_dble *c,su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %0, %%xmm1 \n\t"
                         "shufpd $0x1, %%xmm1, %%xmm1"
                         :
                         :
                         "m" (*c)
                         :
                         "xmm0", "xmm1");

   __asm__ __volatile__ ("movddup %0, %%xmm2 \n\t"
                         "movddup %1, %%xmm3 \n\t"
                         "movddup %2, %%xmm4 \n\t"
                         "movddup %3, %%xmm5 \n\t"
                         "movddup %4, %%xmm6 \n\t"
                         "movddup %5, %%xmm7 \n\t"
                         "mulpd %%xmm0, %%xmm2 \n\t"
                         "mulpd %%xmm1, %%xmm3 \n\t"
                         "mulpd %%xmm0, %%xmm4 \n\t"
                         "mulpd %%xmm1, %%xmm5 \n\t"
                         "mulpd %%xmm0, %%xmm6 \n\t"
                         "mulpd %%xmm1, %%xmm7"                         
                         :
                         :
                         "m" ((*u).c11.re),
                         "m" ((*u).c11.im),
                         "m" ((*u).c12.re),
                         "m" ((*u).c12.im),
                         "m" ((*u).c13.re),
                         "m" ((*u).c13.im)
                         :
                         "xmm2", "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7");

   __asm__ __volatile__ ("addsubpd %%xmm3, %%xmm2 \n\t"
                         "addsubpd %%xmm5, %%xmm4 \n\t"
                         "addsubpd %%xmm7, %%xmm6 \n\t"
                         "movapd %%xmm2, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm6, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c12),
                         "=m" ((*v).c13)
                         :
                         :
                         "xmm2", "xmm4", "xmm6");

   __asm__ __volatile__ ("movddup %0, %%xmm8 \n\t"
                         "movddup %1, %%xmm9 \n\t"
                         "movddup %2, %%xmm10 \n\t"
                         "movddup %3, %%xmm11 \n\t"
                         "movddup %4, %%xmm12 \n\t"
                         "movddup %5, %%xmm13 \n\t"
                         "mulpd %%xmm0, %%xmm8 \n\t"
                         "mulpd %%xmm1, %%xmm9 \n\t"
                         "mulpd %%xmm0, %%xmm10 \n\t"
                         "mulpd %%xmm1, %%xmm11 \n\t"
                         "mulpd %%xmm0, %%xmm12 \n\t"
                         "mulpd %%xmm1, %%xmm13"                         
                         :
                         :
                         "m" ((*u).c21.re),
                         "m" ((*u).c21.im),
                         "m" ((*u).c22.re),
                         "m" ((*u).c22.im),
                         "m" ((*u).c23.re),
                         "m" ((*u).c23.im)
                         :
                         "xmm8", "xmm9", "xmm10", "xmm11",
                         "xmm12", "xmm13");

   __asm__ __volatile__ ("addsubpd %%xmm9, %%xmm8 \n\t"
                         "addsubpd %%xmm11, %%xmm10 \n\t"
                         "addsubpd %%xmm13, %%xmm12 \n\t"
                         "movapd %%xmm8, %0 \n\t"
                         "movapd %%xmm10, %1 \n\t"
                         "movapd %%xmm12, %2"
                         :
                         "=m" ((*v).c21),
                         "=m" ((*v).c22),
                         "=m" ((*v).c23)                         
                         :
                         :
                         "xmm8", "xmm10", "xmm12");

   __asm__ __volatile__ ("movddup %0, %%xmm2 \n\t"
                         "movddup %1, %%xmm3 \n\t"
                         "movddup %2, %%xmm4 \n\t"
                         "movddup %3, %%xmm5 \n\t"
                         "movddup %4, %%xmm6 \n\t"
                         "movddup %5, %%xmm7 \n\t"
                         "mulpd %%xmm0, %%xmm2 \n\t"
                         "mulpd %%xmm1, %%xmm3 \n\t"
                         "mulpd %%xmm0, %%xmm4 \n\t"
                         "mulpd %%xmm1, %%xmm5 \n\t"
                         "mulpd %%xmm0, %%xmm6 \n\t"
                         "mulpd %%xmm1, %%xmm7"                         
                         :
                         :
                         "m" ((*u).c31.re),
                         "m" ((*u).c31.im),
                         "m" ((*u).c32.re),
                         "m" ((*u).c32.im),
                         "m" ((*u).c33.re),
                         "m" ((*u).c33.im)
                         :
                         "xmm2", "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7");

   __asm__ __volatile__ ("addsubpd %%xmm3, %%xmm2 \n\t"
                         "addsubpd %%xmm5, %%xmm4 \n\t"
                         "addsubpd %%xmm7, %%xmm6 \n\t"
                         "movapd %%xmm2, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm6, %2"
                         :
                         "=m" ((*v).c31),
                         "=m" ((*v).c32),
                         "=m" ((*v).c33)                         
                         :
                         :
                         "xmm2", "xmm4", "xmm6");
}


void cm3x3_mulc_add(complex_dble *c,su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %0, %%xmm1 \n\t"
                         "shufpd $0x1, %%xmm1, %%xmm1"
                         :
                         :
                         "m" (*c)
                         :
                         "xmm0", "xmm1");

   __asm__ __volatile__ ("movddup %0, %%xmm2 \n\t"
                         "movddup %1, %%xmm3 \n\t"
                         "movddup %2, %%xmm4 \n\t"
                         "movddup %3, %%xmm5 \n\t"
                         "movddup %4, %%xmm6 \n\t"
                         "movddup %5, %%xmm7 \n\t"
                         "mulpd %%xmm0, %%xmm2 \n\t"
                         "mulpd %%xmm1, %%xmm3 \n\t"
                         "mulpd %%xmm0, %%xmm4 \n\t"
                         "mulpd %%xmm1, %%xmm5 \n\t"
                         "mulpd %%xmm0, %%xmm6 \n\t"
                         "mulpd %%xmm1, %%xmm7"                         
                         :
                         :
                         "m" ((*u).c11.re),
                         "m" ((*u).c11.im),
                         "m" ((*u).c12.re),
                         "m" ((*u).c12.im),
                         "m" ((*u).c13.re),
                         "m" ((*u).c13.im)
                         :
                         "xmm2", "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7");

   __asm__ __volatile__ ("addsubpd %%xmm3, %%xmm2 \n\t"
                         "addsubpd %%xmm5, %%xmm4 \n\t"
                         "addsubpd %%xmm7, %%xmm6 \n\t"
                         "addpd %3, %%xmm2 \n\t"
                         "addpd %4, %%xmm4 \n\t"
                         "addpd %5, %%xmm6 \n\t"                         
                         "movapd %%xmm2, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm6, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c12),
                         "=m" ((*v).c13)                         
                         :
                         "m" ((*v).c11),
                         "m" ((*v).c12),
                         "m" ((*v).c13)                          
                         :
                         "xmm2", "xmm4", "xmm6");

   __asm__ __volatile__ ("movddup %0, %%xmm8 \n\t"
                         "movddup %1, %%xmm9 \n\t"
                         "movddup %2, %%xmm10 \n\t"
                         "movddup %3, %%xmm11 \n\t"
                         "movddup %4, %%xmm12 \n\t"
                         "movddup %5, %%xmm13 \n\t"
                         "mulpd %%xmm0, %%xmm8 \n\t"
                         "mulpd %%xmm1, %%xmm9 \n\t"
                         "mulpd %%xmm0, %%xmm10 \n\t"
                         "mulpd %%xmm1, %%xmm11 \n\t"
                         "mulpd %%xmm0, %%xmm12 \n\t"
                         "mulpd %%xmm1, %%xmm13"                         
                         :
                         :
                         "m" ((*u).c21.re),
                         "m" ((*u).c21.im),
                         "m" ((*u).c22.re),
                         "m" ((*u).c22.im),
                         "m" ((*u).c23.re),
                         "m" ((*u).c23.im)
                         :
                         "xmm8", "xmm9", "xmm10", "xmm11",
                         "xmm12", "xmm13");

   __asm__ __volatile__ ("addsubpd %%xmm9, %%xmm8 \n\t"
                         "addsubpd %%xmm11, %%xmm10 \n\t"
                         "addsubpd %%xmm13, %%xmm12 \n\t"
                         "addpd %3, %%xmm8 \n\t"
                         "addpd %4, %%xmm10 \n\t"
                         "addpd %5, %%xmm12 \n\t"                            
                         "movapd %%xmm8, %0 \n\t"
                         "movapd %%xmm10, %1 \n\t"
                         "movapd %%xmm12, %2"
                         :
                         "=m" ((*v).c21),
                         "=m" ((*v).c22),
                         "=m" ((*v).c23)                         
                         :
                         "m" ((*v).c21),
                         "m" ((*v).c22),
                         "m" ((*v).c23)
                         :
                         "xmm8", "xmm10", "xmm12");

   __asm__ __volatile__ ("movddup %0, %%xmm2 \n\t"
                         "movddup %1, %%xmm3 \n\t"
                         "movddup %2, %%xmm4 \n\t"
                         "movddup %3, %%xmm5 \n\t"
                         "movddup %4, %%xmm6 \n\t"
                         "movddup %5, %%xmm7 \n\t"
                         "mulpd %%xmm0, %%xmm2 \n\t"
                         "mulpd %%xmm1, %%xmm3 \n\t"
                         "mulpd %%xmm0, %%xmm4 \n\t"
                         "mulpd %%xmm1, %%xmm5 \n\t"
                         "mulpd %%xmm0, %%xmm6 \n\t"
                         "mulpd %%xmm1, %%xmm7"                         
                         :
                         :
                         "m" ((*u).c31.re),
                         "m" ((*u).c31.im),
                         "m" ((*u).c32.re),
                         "m" ((*u).c32.im),
                         "m" ((*u).c33.re),
                         "m" ((*u).c33.im)
                         :
                         "xmm2", "xmm3", "xmm4", "xmm5",
                         "xmm6", "xmm7");

   __asm__ __volatile__ ("addsubpd %%xmm3, %%xmm2 \n\t"
                         "addsubpd %%xmm5, %%xmm4 \n\t"
                         "addsubpd %%xmm7, %%xmm6 \n\t"
                         "addpd %3, %%xmm2 \n\t"
                         "addpd %4, %%xmm4 \n\t"
                         "addpd %5, %%xmm6 \n\t"
                         "movapd %%xmm2, %0 \n\t"
                         "movapd %%xmm4, %1 \n\t"
                         "movapd %%xmm6, %2"
                         :
                         "=m" ((*v).c31),
                         "=m" ((*v).c32),
                         "=m" ((*v).c33)                         
                         :
                         "m" ((*v).c31),
                         "m" ((*v).c32),
                         "m" ((*v).c33) 
                         :
                         "xmm2", "xmm4", "xmm6");
}


void cm3x3_lc1(complex_dble *c,su3_dble *u,su3_dble *v)
{
   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %1, %%xmm2 \n\t"                         
                         "shufpd $0x1, %%xmm2, %%xmm2"
                         :
                         :
                         "m" (c[0]),
                         "m" (c[1])
                         :
                         "xmm0", "xmm1", "xmm2");

   __asm__ __volatile__ ("movddup %0, %%xmm3 \n\t"
                         "movddup %1, %%xmm4 \n\t"
                         "movddup %2, %%xmm5 \n\t"
                         "movddup %3, %%xmm6 \n\t"
                         "movddup %4, %%xmm7 \n\t"
                         "movddup %5, %%xmm8 \n\t"
                         "mulpd %%xmm1, %%xmm3 \n\t"
                         "mulpd %%xmm2, %%xmm4 \n\t"
                         "mulpd %%xmm1, %%xmm5 \n\t"
                         "mulpd %%xmm2, %%xmm6 \n\t"
                         "addpd %%xmm0, %%xmm3 \n\t"
                         "mulpd %%xmm1, %%xmm7 \n\t"
                         "mulpd %%xmm2, %%xmm8"                         
                         :
                         :
                         "m" ((*u).c11.re),
                         "m" ((*u).c11.im),
                         "m" ((*u).c12.re),
                         "m" ((*u).c12.im),
                         "m" ((*u).c13.re),
                         "m" ((*u).c13.im)
                         :
                         "xmm3", "xmm4", "xmm5", "xmm6",
                         "xmm7", "xmm8");

   __asm__ __volatile__ ("addsubpd %%xmm4, %%xmm3 \n\t"
                         "addsubpd %%xmm6, %%xmm5 \n\t"
                         "addsubpd %%xmm8, %%xmm7 \n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm5, %1 \n\t"
                         "movapd %%xmm7, %2"
                         :
                         "=m" ((*v).c11),
                         "=m" ((*v).c12),
                         "=m" ((*v).c13)
                         :
                         :
                         "xmm3", "xmm5", "xmm7");

   __asm__ __volatile__ ("movddup %0, %%xmm9 \n\t"
                         "movddup %1, %%xmm10 \n\t"
                         "movddup %2, %%xmm11 \n\t"
                         "movddup %3, %%xmm12 \n\t"
                         "movddup %4, %%xmm13 \n\t"
                         "movddup %5, %%xmm14 \n\t"
                         "mulpd %%xmm1, %%xmm9 \n\t"
                         "mulpd %%xmm2, %%xmm10 \n\t"
                         "mulpd %%xmm1, %%xmm11 \n\t"
                         "mulpd %%xmm2, %%xmm12 \n\t"
                         "mulpd %%xmm1, %%xmm13 \n\t"
                         "addpd %%xmm0, %%xmm11 \n\t"
                         "mulpd %%xmm2, %%xmm14"                         
                         :
                         :
                         "m" ((*u).c21.re),
                         "m" ((*u).c21.im),
                         "m" ((*u).c22.re),
                         "m" ((*u).c22.im),
                         "m" ((*u).c23.re),
                         "m" ((*u).c23.im)
                         :
                         "xmm9", "xmm10", "xmm11", "xmm12",
                         "xmm13", "xmm14");

   __asm__ __volatile__ ("addsubpd %%xmm10, %%xmm9 \n\t"
                         "addsubpd %%xmm12, %%xmm11 \n\t"
                         "addsubpd %%xmm14, %%xmm13 \n\t"
                         "movapd %%xmm9, %0 \n\t"
                         "movapd %%xmm11, %1 \n\t"
                         "movapd %%xmm13, %2"
                         :
                         "=m" ((*v).c21),
                         "=m" ((*v).c22),
                         "=m" ((*v).c23)                         
                         :
                         :
                         "xmm9", "xmm11", "xmm13");

   __asm__ __volatile__ ("movddup %0, %%xmm3 \n\t"
                         "movddup %1, %%xmm4 \n\t"
                         "movddup %2, %%xmm5 \n\t"
                         "movddup %3, %%xmm6 \n\t"
                         "movddup %4, %%xmm7 \n\t"
                         "movddup %5, %%xmm8 \n\t"
                         "mulpd %%xmm1, %%xmm3 \n\t"
                         "mulpd %%xmm2, %%xmm4 \n\t"
                         "mulpd %%xmm1, %%xmm5 \n\t"
                         "mulpd %%xmm2, %%xmm6 \n\t"
                         "mulpd %%xmm1, %%xmm7 \n\t"
                         "mulpd %%xmm2, %%xmm8 \n\t"
                         "addpd %%xmm0, %%xmm7"                         
                         :
                         :
                         "m" ((*u).c31.re),
                         "m" ((*u).c31.im),
                         "m" ((*u).c32.re),
                         "m" ((*u).c32.im),
                         "m" ((*u).c33.re),
                         "m" ((*u).c33.im)
                         :
                         "xmm3", "xmm4", "xmm5", "xmm6",
                         "xmm7", "xmm8");

   __asm__ __volatile__ ("addsubpd %%xmm4, %%xmm3 \n\t"
                         "addsubpd %%xmm6, %%xmm5 \n\t"
                         "addsubpd %%xmm8, %%xmm7 \n\t"
                         "movapd %%xmm3, %0 \n\t"
                         "movapd %%xmm5, %1 \n\t"
                         "movapd %%xmm7, %2"
                         :
                         "=m" ((*v).c31),
                         "=m" ((*v).c32),
                         "=m" ((*v).c33)
                         :
                         :
                         "xmm3", "xmm5", "xmm7");
}

#else

static su3_dble ww;


void cm3x3_zero(int vol,su3_dble *u)
{
   su3_dble *um;

   um=u+vol;

   for (;u<um;u++)
   {
      (*u).c11.re=0.0;
      (*u).c11.im=0.0;
      (*u).c12.re=0.0;
      (*u).c12.im=0.0;
      (*u).c13.re=0.0;
      (*u).c13.im=0.0;   

      (*u).c21.re=0.0;
      (*u).c21.im=0.0;
      (*u).c22.re=0.0;
      (*u).c22.im=0.0;
      (*u).c23.re=0.0;
      (*u).c23.im=0.0;

      (*u).c31.re=0.0;
      (*u).c31.im=0.0;
      (*u).c32.re=0.0;
      (*u).c32.im=0.0;
      (*u).c33.re=0.0;
      (*u).c33.im=0.0;
   }
}


void cm3x3_unity(int vol,su3_dble *u)
{
   su3_dble *um;

   um=u+vol;

   for (;u<um;u++)
   {   
      (*u).c11.re=1.0;
      (*u).c11.im=0.0;
      (*u).c12.re=0.0;
      (*u).c12.im=0.0;
      (*u).c13.re=0.0;
      (*u).c13.im=0.0;   

      (*u).c21.re=0.0;
      (*u).c21.im=0.0;
      (*u).c22.re=1.0;
      (*u).c22.im=0.0;
      (*u).c23.re=0.0;
      (*u).c23.im=0.0;

      (*u).c31.re=0.0;
      (*u).c31.im=0.0;
      (*u).c32.re=0.0;
      (*u).c32.im=0.0;
      (*u).c33.re=1.0;
      (*u).c33.im=0.0;
   }
}


void cm3x3_assign(int vol,su3_dble *u,su3_dble *v)
{
   su3_dble *um;

   um=u+vol;

   for (;u<um;u++)
   {    
      (*v).c11.re=(*u).c11.re;
      (*v).c11.im=(*u).c11.im;
      (*v).c12.re=(*u).c12.re;
      (*v).c12.im=(*u).c12.im;
      (*v).c13.re=(*u).c13.re;
      (*v).c13.im=(*u).c13.im;

      (*v).c21.re=(*u).c21.re;
      (*v).c21.im=(*u).c21.im;
      (*v).c22.re=(*u).c22.re;
      (*v).c22.im=(*u).c22.im;
      (*v).c23.re=(*u).c23.re;
      (*v).c23.im=(*u).c23.im;    

      (*v).c31.re=(*u).c31.re;
      (*v).c31.im=(*u).c31.im;
      (*v).c32.re=(*u).c32.re;
      (*v).c32.im=(*u).c32.im;
      (*v).c33.re=(*u).c33.re;
      (*v).c33.im=(*u).c33.im;

      v+=1;
   }
}


void cm3x3_swap(int vol,su3_dble *u,su3_dble *v)
{
   complex_dble z;
   su3_dble *um;

   um=u+vol;

   for (;u<um;u++)
   {
      z.re=(*v).c11.re;
      z.im=(*v).c11.im;
      (*v).c11.re=(*u).c11.re;
      (*v).c11.im=(*u).c11.im;
      (*u).c11.re=z.re;
      (*u).c11.im=z.im;

      z.re=(*v).c12.re;
      z.im=(*v).c12.im;
      (*v).c12.re=(*u).c12.re;
      (*v).c12.im=(*u).c12.im;
      (*u).c12.re=z.re;
      (*u).c12.im=z.im;

      z.re=(*v).c13.re;
      z.im=(*v).c13.im;
      (*v).c13.re=(*u).c13.re;
      (*v).c13.im=(*u).c13.im;
      (*u).c13.re=z.re;
      (*u).c13.im=z.im;
      
      z.re=(*v).c21.re;
      z.im=(*v).c21.im;
      (*v).c21.re=(*u).c21.re;
      (*v).c21.im=(*u).c21.im;
      (*u).c21.re=z.re;
      (*u).c21.im=z.im;

      z.re=(*v).c22.re;
      z.im=(*v).c22.im;
      (*v).c22.re=(*u).c22.re;
      (*v).c22.im=(*u).c22.im;
      (*u).c22.re=z.re;
      (*u).c22.im=z.im;

      z.re=(*v).c23.re;
      z.im=(*v).c23.im;
      (*v).c23.re=(*u).c23.re;
      (*v).c23.im=(*u).c23.im;
      (*u).c23.re=z.re;
      (*u).c23.im=z.im;

      z.re=(*v).c31.re;
      z.im=(*v).c31.im;
      (*v).c31.re=(*u).c31.re;
      (*v).c31.im=(*u).c31.im;
      (*u).c31.re=z.re;
      (*u).c31.im=z.im;

      z.re=(*v).c32.re;
      z.im=(*v).c32.im;
      (*v).c32.re=(*u).c32.re;
      (*v).c32.im=(*u).c32.im;
      (*u).c32.re=z.re;
      (*u).c32.im=z.im;

      z.re=(*v).c33.re;
      z.im=(*v).c33.im;
      (*v).c33.re=(*u).c33.re;
      (*v).c33.im=(*u).c33.im;
      (*u).c33.re=z.re;
      (*u).c33.im=z.im;
      
      v+=1;
   }
}


void cm3x3_dagger(su3_dble *u,su3_dble *v)
{
   complex_dble z;
   
   (*v).c11.re= (*u).c11.re;
   (*v).c11.im=-(*u).c11.im;
   (*v).c22.re= (*u).c22.re;
   (*v).c22.im=-(*u).c22.im;
   (*v).c33.re= (*u).c33.re;
   (*v).c33.im=-(*u).c33.im;
   
   z.re= (*u).c12.re;
   z.im=-(*u).c12.im;
   (*v).c12.re= (*u).c21.re;
   (*v).c12.im=-(*u).c21.im;
   (*v).c21.re=z.re;
   (*v).c21.im=z.im;

   z.re= (*u).c13.re;
   z.im=-(*u).c13.im;
   (*v).c13.re= (*u).c31.re;
   (*v).c13.im=-(*u).c31.im;
   (*v).c31.re=z.re;
   (*v).c31.im=z.im;   

   z.re= (*u).c23.re;
   z.im=-(*u).c23.im;
   (*v).c23.re= (*u).c32.re;
   (*v).c23.im=-(*u).c32.im;
   (*v).c32.re=z.re;
   (*v).c32.im=z.im;   
}


void cm3x3_tr(su3_dble *u,su3_dble *v,complex_dble *tr)
{
   complex_dble z;

   z.re =(*u).c11.re*(*v).c11.re-(*u).c11.im*(*v).c11.im;
   z.im =(*u).c11.re*(*v).c11.im+(*u).c11.im*(*v).c11.re;
   z.re+=(*u).c12.re*(*v).c21.re-(*u).c12.im*(*v).c21.im;
   z.im+=(*u).c12.re*(*v).c21.im+(*u).c12.im*(*v).c21.re;
   z.re+=(*u).c13.re*(*v).c31.re-(*u).c13.im*(*v).c31.im;
   z.im+=(*u).c13.re*(*v).c31.im+(*u).c13.im*(*v).c31.re;

   z.re+=(*u).c21.re*(*v).c12.re-(*u).c21.im*(*v).c12.im;
   z.im+=(*u).c21.re*(*v).c12.im+(*u).c21.im*(*v).c12.re;
   z.re+=(*u).c22.re*(*v).c22.re-(*u).c22.im*(*v).c22.im;
   z.im+=(*u).c22.re*(*v).c22.im+(*u).c22.im*(*v).c22.re;
   z.re+=(*u).c23.re*(*v).c32.re-(*u).c23.im*(*v).c32.im;
   z.im+=(*u).c23.re*(*v).c32.im+(*u).c23.im*(*v).c32.re;

   z.re+=(*u).c31.re*(*v).c13.re-(*u).c31.im*(*v).c13.im;
   z.im+=(*u).c31.re*(*v).c13.im+(*u).c31.im*(*v).c13.re;
   z.re+=(*u).c32.re*(*v).c23.re-(*u).c32.im*(*v).c23.im;
   z.im+=(*u).c32.re*(*v).c23.im+(*u).c32.im*(*v).c23.re;
   z.re+=(*u).c33.re*(*v).c33.re-(*u).c33.im*(*v).c33.im;
   z.im+=(*u).c33.re*(*v).c33.im+(*u).c33.im*(*v).c33.re;     
   
   (*tr).re=z.re;
   (*tr).im=z.im;
}


void cm3x3_retr(su3_dble *u,su3_dble *v,double *tr)
{
   double r;

   r =(*u).c11.re*(*v).c11.re-(*u).c11.im*(*v).c11.im;
   r+=(*u).c12.re*(*v).c21.re-(*u).c12.im*(*v).c21.im;
   r+=(*u).c13.re*(*v).c31.re-(*u).c13.im*(*v).c31.im;

   r+=(*u).c21.re*(*v).c12.re-(*u).c21.im*(*v).c12.im;
   r+=(*u).c22.re*(*v).c22.re-(*u).c22.im*(*v).c22.im;
   r+=(*u).c23.re*(*v).c32.re-(*u).c23.im*(*v).c32.im;

   r+=(*u).c31.re*(*v).c13.re-(*u).c31.im*(*v).c13.im;
   r+=(*u).c32.re*(*v).c23.re-(*u).c32.im*(*v).c23.im;
   r+=(*u).c33.re*(*v).c33.re-(*u).c33.im*(*v).c33.im;
   
   (*tr)=r;
}


void cm3x3_imtr(su3_dble *u,su3_dble *v,double *tr)
{
   double r;

   r =(*u).c11.re*(*v).c11.im+(*u).c11.im*(*v).c11.re;
   r+=(*u).c12.re*(*v).c21.im+(*u).c12.im*(*v).c21.re;
   r+=(*u).c13.re*(*v).c31.im+(*u).c13.im*(*v).c31.re;

   r+=(*u).c21.re*(*v).c12.im+(*u).c21.im*(*v).c12.re;
   r+=(*u).c22.re*(*v).c22.im+(*u).c22.im*(*v).c22.re;
   r+=(*u).c23.re*(*v).c32.im+(*u).c23.im*(*v).c32.re;

   r+=(*u).c31.re*(*v).c13.im+(*u).c31.im*(*v).c13.re;
   r+=(*u).c32.re*(*v).c23.im+(*u).c32.im*(*v).c23.re;
   r+=(*u).c33.re*(*v).c33.im+(*u).c33.im*(*v).c33.re;
   
   (*tr)=r;
}


void cm3x3_add(su3_dble *u,su3_dble *v)
{
   (*v).c11.re+=(*u).c11.re;
   (*v).c11.im+=(*u).c11.im;
   (*v).c12.re+=(*u).c12.re;
   (*v).c12.im+=(*u).c12.im;
   (*v).c13.re+=(*u).c13.re;
   (*v).c13.im+=(*u).c13.im;

   (*v).c21.re+=(*u).c21.re;
   (*v).c21.im+=(*u).c21.im;
   (*v).c22.re+=(*u).c22.re;
   (*v).c22.im+=(*u).c22.im;
   (*v).c23.re+=(*u).c23.re;
   (*v).c23.im+=(*u).c23.im;    

   (*v).c31.re+=(*u).c31.re;
   (*v).c31.im+=(*u).c31.im;
   (*v).c32.re+=(*u).c32.re;
   (*v).c32.im+=(*u).c32.im;
   (*v).c33.re+=(*u).c33.re;
   (*v).c33.im+=(*u).c33.im; 
}


void cm3x3_mul_add(su3_dble *u,su3_dble *v,su3_dble *w)
{
   su3xsu3(u,v,&ww);
   cm3x3_add(&ww,w);
}


void cm3x3_mulr(double *r,su3_dble *u,su3_dble *v)
{
   double rr;

   rr=(*r);

   (*v).c11.re=rr*(*u).c11.re;
   (*v).c11.im=rr*(*u).c11.im;
   (*v).c12.re=rr*(*u).c12.re;
   (*v).c12.im=rr*(*u).c12.im;
   (*v).c13.re=rr*(*u).c13.re;
   (*v).c13.im=rr*(*u).c13.im;

   (*v).c21.re=rr*(*u).c21.re;
   (*v).c21.im=rr*(*u).c21.im;
   (*v).c22.re=rr*(*u).c22.re;
   (*v).c22.im=rr*(*u).c22.im;
   (*v).c23.re=rr*(*u).c23.re;
   (*v).c23.im=rr*(*u).c23.im;

   (*v).c31.re=rr*(*u).c31.re;
   (*v).c31.im=rr*(*u).c31.im;
   (*v).c32.re=rr*(*u).c32.re;
   (*v).c32.im=rr*(*u).c32.im;
   (*v).c33.re=rr*(*u).c33.re;
   (*v).c33.im=rr*(*u).c33.im;
}


void cm3x3_mulr_add(double *r,su3_dble *u,su3_dble *v)
{
   double rr;

   rr=(*r);

   (*v).c11.re+=rr*(*u).c11.re;
   (*v).c11.im+=rr*(*u).c11.im;
   (*v).c12.re+=rr*(*u).c12.re;
   (*v).c12.im+=rr*(*u).c12.im;
   (*v).c13.re+=rr*(*u).c13.re;
   (*v).c13.im+=rr*(*u).c13.im;

   (*v).c21.re+=rr*(*u).c21.re;
   (*v).c21.im+=rr*(*u).c21.im;
   (*v).c22.re+=rr*(*u).c22.re;
   (*v).c22.im+=rr*(*u).c22.im;
   (*v).c23.re+=rr*(*u).c23.re;
   (*v).c23.im+=rr*(*u).c23.im;

   (*v).c31.re+=rr*(*u).c31.re;
   (*v).c31.im+=rr*(*u).c31.im;
   (*v).c32.re+=rr*(*u).c32.re;
   (*v).c32.im+=rr*(*u).c32.im;
   (*v).c33.re+=rr*(*u).c33.re;
   (*v).c33.im+=rr*(*u).c33.im;
}


void cm3x3_mulc(complex_dble *c,su3_dble *u,su3_dble *v)
{
   double rr;
   complex_dble cc;

   cc.re=(*c).re;
   cc.im=(*c).im;

   rr=          cc.re*(*u).c11.im+cc.im*(*u).c11.re;   
   (*v).c11.re=(cc.re*(*u).c11.re-cc.im*(*u).c11.im);
   (*v).c11.im=rr;
   rr=          cc.re*(*u).c12.im+cc.im*(*u).c12.re;   
   (*v).c12.re=(cc.re*(*u).c12.re-cc.im*(*u).c12.im);   
   (*v).c12.im=rr;
   rr=          cc.re*(*u).c13.im+cc.im*(*u).c13.re;   
   (*v).c13.re=(cc.re*(*u).c13.re-cc.im*(*u).c13.im);   
   (*v).c13.im=rr;

   rr=          cc.re*(*u).c21.im+cc.im*(*u).c21.re;   
   (*v).c21.re=(cc.re*(*u).c21.re-cc.im*(*u).c21.im);
   (*v).c21.im=rr;
   rr=          cc.re*(*u).c22.im+cc.im*(*u).c22.re;   
   (*v).c22.re=(cc.re*(*u).c22.re-cc.im*(*u).c22.im);   
   (*v).c22.im=rr;
   rr=          cc.re*(*u).c23.im+cc.im*(*u).c23.re;   
   (*v).c23.re=(cc.re*(*u).c23.re-cc.im*(*u).c23.im);   
   (*v).c23.im=rr;

   rr=          cc.re*(*u).c31.im+cc.im*(*u).c31.re;   
   (*v).c31.re=(cc.re*(*u).c31.re-cc.im*(*u).c31.im);
   (*v).c31.im=rr;
   rr=          cc.re*(*u).c32.im+cc.im*(*u).c32.re;   
   (*v).c32.re=(cc.re*(*u).c32.re-cc.im*(*u).c32.im);   
   (*v).c32.im=rr;
   rr=          cc.re*(*u).c33.im+cc.im*(*u).c33.re;   
   (*v).c33.re=(cc.re*(*u).c33.re-cc.im*(*u).c33.im);   
   (*v).c33.im=rr;
}


void cm3x3_mulc_add(complex_dble *c,su3_dble *u,su3_dble *v)
{
   double rr;
   complex_dble cc;

   cc.re=(*c).re;
   cc.im=(*c).im;
   
   rr=           cc.re*(*u).c11.im+cc.im*(*u).c11.re;   
   (*v).c11.re+=(cc.re*(*u).c11.re-cc.im*(*u).c11.im);
   (*v).c11.im+=rr;
   rr=           cc.re*(*u).c12.im+cc.im*(*u).c12.re;   
   (*v).c12.re+=(cc.re*(*u).c12.re-cc.im*(*u).c12.im);   
   (*v).c12.im+=rr;
   rr=           cc.re*(*u).c13.im+cc.im*(*u).c13.re;   
   (*v).c13.re+=(cc.re*(*u).c13.re-cc.im*(*u).c13.im);   
   (*v).c13.im+=rr;

   rr=           cc.re*(*u).c21.im+cc.im*(*u).c21.re;   
   (*v).c21.re+=(cc.re*(*u).c21.re-cc.im*(*u).c21.im);
   (*v).c21.im+=rr;
   rr=           cc.re*(*u).c22.im+cc.im*(*u).c22.re;   
   (*v).c22.re+=(cc.re*(*u).c22.re-cc.im*(*u).c22.im);   
   (*v).c22.im+=rr;
   rr=           cc.re*(*u).c23.im+cc.im*(*u).c23.re;   
   (*v).c23.re+=(cc.re*(*u).c23.re-cc.im*(*u).c23.im);   
   (*v).c23.im+=rr;

   rr=           cc.re*(*u).c31.im+cc.im*(*u).c31.re;   
   (*v).c31.re+=(cc.re*(*u).c31.re-cc.im*(*u).c31.im);
   (*v).c31.im+=rr;
   rr=           cc.re*(*u).c32.im+cc.im*(*u).c32.re;   
   (*v).c32.re+=(cc.re*(*u).c32.re-cc.im*(*u).c32.im);   
   (*v).c32.im+=rr;
   rr=           cc.re*(*u).c33.im+cc.im*(*u).c33.re;   
   (*v).c33.re+=(cc.re*(*u).c33.re-cc.im*(*u).c33.im);   
   (*v).c33.im+=rr;
}


void cm3x3_lc1(complex_dble *c,su3_dble *u,su3_dble *v)
{
   complex_dble cc;
   
   cm3x3_mulc(c+1,u,v);

   cc.re=(*c).re;
   cc.im=(*c).im;
   
   (*v).c11.re+=cc.re;
   (*v).c11.im+=cc.im;
   (*v).c22.re+=cc.re;
   (*v).c22.im+=cc.im;
   (*v).c33.re+=cc.re;
   (*v).c33.im+=cc.im;   
}

#endif

void cm3x3_lc2(complex_dble *c,su3_dble *u,su3_dble *v)
{
   cm3x3_lc1(c,u,v);
   cm3x3_mulc_add(c+2,u+1,v);
}
