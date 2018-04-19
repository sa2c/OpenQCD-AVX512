
/*******************************************************************************
*
* File Aw_gen.c
*
* Copyright (C) 2007, 2008, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic programs needed for the computation of the little Dirac operator
*
* The externally accessible functions are
*
*   void gather_ud(int vol,int *imb,su3_dble *ud,su3_dble *vd)
*     Assigns the 3x3 matrices ud[imb[i]] to vd[i] (i=0,..,vol-1).
*
*   void gather_sd(int vol,int *imb,spinor_dble *sd,spinor_dble *rd)
*     Assigns the spinors sd[imb[i]] to rd[i] (i=0,..,vol-1).
*  
*   void apply_u2sd(int vol,int *imb,su3_dble *ud,spinor_dble *sd,
*                   spinor_dble *rd)
*     Multiplies the spinors sd[imb[i]] by the 3x3 matrices ud[i] and 
*     assigns the result to rd[i] (i=0,..,vol-1). 
*
*   void apply_udag2sd(int vol,int *imb,su3_dble *ud,spinor_dble *sd,
*                      spinor_dble *rd)
*     Multiplies the spinors sd[imb[i]] by the adjoint of the 3x3 matrices
*     ud[i] and assigns the result to rd[i] (i=0,..,vol-1).
*
* The following is an array of functions indexed by the direction mu=0,..,3:
*
*   void (*spinor_prod_gamma[])(int vol,spinor_dble *sd,spinor_dble *rd,
*                               complex_dble *sp)
*      Computes the scalar products (sd,rd) and (sd,gamma_mu*rd), where
*      gamma_mu denotes the Dirac matrix with index mu and the spinor
*      fields are assumed to have vol elements. On exit the calculated
*      products are assigned to sp[0] and sp[1], respectively.
*
* Notes:
*
* The representation of the Dirac matrices is specified in the notes
* "Implementation of the lattice Dirac operator" (file doc/dirac.pdf).
* The input and output fields may not overlap in the case of the programs
* gather_ud(), gather_sd(), apply_u2sd() and apply_udag2sd().
*
* All these programs can be called locally. If SSE inline-assembly is used
* (i.e. if x64 is set), it is taken for granted  that the field arrays are
* aligned to 16 byte boundaries.
*
*******************************************************************************/

#define AW_GEN_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "little.h"

#define MAX_LEVELS 8
#define BLK_LENGTH 8

static int cnt[MAX_LEVELS];
static complex_dble sm0[MAX_LEVELS] ALIGNED16;
static complex_dble sm1[MAX_LEVELS] ALIGNED16;


static void init_sm(void)
{
   int n;

   for (n=0;n<MAX_LEVELS;n++)
   {
      cnt[n]=0;
      sm0[n].re=0.0;
      sm0[n].im=0.0;
      sm1[n].re=0.0;
      sm1[n].im=0.0;
   }   
}


static void acc_sm(void) 
{
   int n;

   cnt[0]+=1;
   
   for (n=1;(cnt[n-1]>=BLK_LENGTH)&&(n<MAX_LEVELS);n++)
   {
      cnt[n]+=1;
      sm0[n].re+=sm0[n-1].re;
      sm0[n].im+=sm0[n-1].im;
      sm1[n].re+=sm1[n-1].re;
      sm1[n].im+=sm1[n-1].im;      

      cnt[n-1]=0;
      sm0[n-1].re=0.0;
      sm0[n-1].im=0.0;
      sm1[n-1].re=0.0;
      sm1[n-1].im=0.0;
   }
}


static void sum_sm(void)
{
   int n;

   for (n=1;n<MAX_LEVELS;n++)
   {
      sm0[0].re+=sm0[n].re;
      sm0[0].im+=sm0[n].im;      
      sm1[0].re+=sm1[n].re;
      sm1[0].im+=sm1[n].im;      
   }
}


void gather_ud(int vol,int *imb,su3_dble *ud,su3_dble *vd)
{
   int *imm;

   imm=imb+vol;

   for (;imb<imm;imb++)
   {
      (*vd)=ud[*imb];
      vd+=1;
   }
}


void gather_sd(int vol,int *imb,spinor_dble *sd,spinor_dble *rd)
{
   int *imm;

   imm=imb+vol;

   for (;imb<imm;imb++)
   {
      (*rd)=sd[*imb];
      rd+=1;
   }
}

#if (defined x64)
#include "sse2.h"

#define _start_sm() \
__asm__ __volatile__ ("xorpd %%xmm12, %%xmm12 \n\t" \
                      "xorpd %%xmm13, %%xmm13 \n\t" \
                      "xorpd %%xmm14, %%xmm14 \n\t" \
                      "xorpd %%xmm15, %%xmm15" \
                      : \
                      : \
                      : \
                      "xmm12", "xmm13", "xmm14", "xmm15")

#define _store_sm() \
__asm__ __volatile__ ("shufpd $0x1, %%xmm12, %%xmm12 \n\t" \
                      "shufpd $0x1, %%xmm14, %%xmm14 \n\t" \
                      "addsubpd %%xmm12, %%xmm13 \n\t" \
                      "addsubpd %%xmm14, %%xmm15 \n\t" \
                      "shufpd $0x1, %%xmm13, %%xmm13 \n\t" \
                      "shufpd $0x1, %%xmm15, %%xmm15 \n\t" \
                      "addpd %2, %%xmm13 \n\t" \
                      "addpd %3, %%xmm15 \n\t" \
                      "movapd %%xmm13, %0 \n\t" \
                      "movapd %%xmm15, %1 \n\t" \
                      : \
                      "=m" (sm0[0]), \
                      "=m" (sm1[0]) \
                      : \
                      "m" (sm0[0]), \
                      "m" (sm1[0]) \
                      : \
                      "xmm12", "xmm13", "xmm14", "xmm15")

#define _load_chi(s) \
__asm__ __volatile__ ("movddup %0, %%xmm0 \n\t" \
                      "movddup %1, %%xmm1 \n\t" \
                      "movddup %2, %%xmm2 \n\t" \
                      "movddup %3, %%xmm3 \n\t" \
                      "movddup %4, %%xmm4 \n\t" \
                      "movddup %5, %%xmm5" \
                      : \
                      : \
                      "m" ((s).c1.re), \
                      "m" ((s).c2.re), \
                      "m" ((s).c3.re), \
                      "m" ((s).c1.im), \
                      "m" ((s).c2.im), \
                      "m" ((s).c3.im) \
                      : \
                      "xmm0", "xmm1", "xmm2", "xmm3", \
                      "xmm4", "xmm5")

#define _load_psi0(s) \
__asm__ __volatile__ ("movapd %0, %%xmm6 \n\t" \
                      "movapd %1, %%xmm7 \n\t" \
                      "movapd %2, %%xmm8 \n\t" \
                      "movapd %0, %%xmm9 \n\t" \
                      "movapd %1, %%xmm10 \n\t" \
                      "movapd %2, %%xmm11 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm8 \n\t" \
                      "mulpd %%xmm3, %%xmm9 \n\t" \
                      "mulpd %%xmm4, %%xmm10 \n\t" \
                      "mulpd %%xmm5, %%xmm11" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm6", "xmm7", "xmm8", "xmm9", \
                      "xmm10", "xmm11")

#define _load_psi1_add(s) \
__asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t" \
                      "mulpd %1, %%xmm1 \n\t" \
                      "mulpd %2, %%xmm2 \n\t" \
                      "mulpd %0, %%xmm3 \n\t" \
                      "mulpd %1, %%xmm4 \n\t" \
                      "mulpd %2, %%xmm5 \n\t" \
                      "addpd %%xmm6, %%xmm12 \n\t" \
                      "addpd %%xmm9, %%xmm13 \n\t" \
                      "addpd %%xmm0, %%xmm14 \n\t" \
                      "addpd %%xmm3, %%xmm15 \n\t" \
                      "addpd %%xmm7, %%xmm12 \n\t" \
                      "addpd %%xmm10, %%xmm13 \n\t" \
                      "addpd %%xmm1, %%xmm14 \n\t" \
                      "addpd %%xmm4, %%xmm15 \n\t" \
                      "addpd %%xmm8, %%xmm12 \n\t" \
                      "addpd %%xmm11, %%xmm13 \n\t" \
                      "addpd %%xmm2, %%xmm14 \n\t" \
                      "addpd %%xmm5, %%xmm15" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm0", "xmm1", "xmm2", "xmm3", \
                      "xmm4", "xmm5", "xmm12", "xmm13", \
                      "xmm14", "xmm15")

#define _load_psi1_sub(s) \
__asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t" \
                      "mulpd %1, %%xmm1 \n\t" \
                      "mulpd %2, %%xmm2 \n\t" \
                      "mulpd %0, %%xmm3 \n\t" \
                      "mulpd %1, %%xmm4 \n\t" \
                      "mulpd %2, %%xmm5 \n\t" \
                      "addpd %%xmm6, %%xmm12 \n\t" \
                      "addpd %%xmm9, %%xmm13 \n\t" \
                      "subpd %%xmm0, %%xmm14 \n\t" \
                      "subpd %%xmm3, %%xmm15 \n\t" \
                      "addpd %%xmm7, %%xmm12 \n\t" \
                      "addpd %%xmm10, %%xmm13 \n\t" \
                      "subpd %%xmm1, %%xmm14 \n\t" \
                      "subpd %%xmm4, %%xmm15 \n\t" \
                      "addpd %%xmm8, %%xmm12 \n\t" \
                      "addpd %%xmm11, %%xmm13 \n\t" \
                      "subpd %%xmm2, %%xmm14 \n\t" \
                      "subpd %%xmm5, %%xmm15" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm0", "xmm1", "xmm2", "xmm3", \
                      "xmm4", "xmm5", "xmm12", "xmm13", \
                      "xmm14", "xmm15")


void apply_u2sd(int vol,int *imb,su3_dble *ud,spinor_dble *sd,
                spinor_dble *rd)
{
   int *imm;
   spinor_dble *si;

   imm=imb+vol;

   while (imb<imm)
   {
      si=sd+(*imb);
      imb+=1;

      if (imb<imm)
         _prefetch_spinor_dble(sd+(*imb));

      _sse_load_dble((*si).c1);
      _sse_su3_multiply_dble(*ud);
      _sse_store_up_dble((*rd).c1);

      _sse_load_dble((*si).c2);
      _sse_su3_multiply_dble(*ud);
      _sse_store_up_dble((*rd).c2);

      _sse_load_dble((*si).c3);
      _sse_su3_multiply_dble(*ud);
      _sse_store_up_dble((*rd).c3);

      _sse_load_dble((*si).c4);
      _sse_su3_multiply_dble(*ud);
      _sse_store_up_dble((*rd).c4);      

      ud+=1;
      rd+=1;
   }
}


void apply_udag2sd(int vol,int *imb,su3_dble *ud,spinor_dble *sd,
                   spinor_dble *rd)
{
   int *imm;
   spinor_dble *si;

   imm=imb+vol;

   while (imb<imm)
   {
      si=sd+(*imb);
      imb+=1;

      if (imb<imm)
         _prefetch_spinor_dble(sd+(*imb));

      _sse_load_dble((*si).c1);
      _sse_su3_inverse_multiply_dble(*ud);
      _sse_store_up_dble((*rd).c1);

      _sse_load_dble((*si).c2);
      _sse_su3_inverse_multiply_dble(*ud);
      _sse_store_up_dble((*rd).c2);

      _sse_load_dble((*si).c3);
      _sse_su3_inverse_multiply_dble(*ud);
      _sse_store_up_dble((*rd).c3);

      _sse_load_dble((*si).c4);
      _sse_su3_inverse_multiply_dble(*ud);
      _sse_store_up_dble((*rd).c4);      

      ud+=1;
      rd+=1;
   }
}


static void spinor_prod_gamma0(int vol,spinor_dble *sd,spinor_dble *rd,
                               complex_dble *sp)
{
   spinor_dble *rt,*rm;

   init_sm();
   rt=rd+vol;
   rm=rd;

   while (rm<rt)
   {
      rm+=BLK_LENGTH;
      if (rm>rt)
         rm=rt;

      _start_sm();      

      for (;rd<rm;rd++)
      {
         _load_chi((*rd).c1);
         _load_psi0((*sd).c1);
         _load_psi1_add((*sd).c3);

         _load_chi((*rd).c2);
         _load_psi0((*sd).c2);
         _load_psi1_add((*sd).c4);         

         _load_chi((*rd).c3);
         _load_psi0((*sd).c3);
         _load_psi1_add((*sd).c1);

         _load_chi((*rd).c4);
         _load_psi0((*sd).c4);
         _load_psi1_add((*sd).c2);         
         
         sd+=1;
      }

      _store_sm();
      acc_sm();
   }

   sum_sm();
   sp[0].re=sm0[0].re;
   sp[0].im=sm0[0].im;
   sp[1].re=-sm1[0].re;
   sp[1].im=-sm1[0].im;
}


static void spinor_prod_gamma1(int vol,spinor_dble *sd,spinor_dble *rd,
                               complex_dble *sp)
{
   spinor_dble *rt,*rm;

   init_sm();
   rt=rd+vol;
   rm=rd;

   while (rm<rt)
   {
      rm+=BLK_LENGTH;
      if (rm>rt)
         rm=rt;

      _start_sm();      

      for (;rd<rm;rd++)
      {
         _load_chi((*rd).c1);
         _load_psi0((*sd).c1);
         _load_psi1_sub((*sd).c4);

         _load_chi((*rd).c2);
         _load_psi0((*sd).c2);
         _load_psi1_sub((*sd).c3);         

         _load_chi((*rd).c3);
         _load_psi0((*sd).c3);
         _load_psi1_add((*sd).c2);

         _load_chi((*rd).c4);
         _load_psi0((*sd).c4);
         _load_psi1_add((*sd).c1);         
         
         sd+=1;
      }

      _store_sm();
      acc_sm();
   }

   sum_sm();
   sp[0].re=sm0[0].re;
   sp[0].im=sm0[0].im;
   sp[1].re=sm1[0].im;
   sp[1].im=-sm1[0].re;
}


static void spinor_prod_gamma2(int vol,spinor_dble *sd,spinor_dble *rd,
                               complex_dble *sp)
{
   spinor_dble *rt,*rm;

   init_sm();
   rt=rd+vol;
   rm=rd;

   while (rm<rt)
   {
      rm+=BLK_LENGTH;
      if (rm>rt)
         rm=rt;

      _start_sm();      

      for (;rd<rm;rd++)
      {
         _load_chi((*rd).c1);
         _load_psi0((*sd).c1);
         _load_psi1_add((*sd).c4);

         _load_chi((*rd).c2);
         _load_psi0((*sd).c2);
         _load_psi1_sub((*sd).c3);         

         _load_chi((*rd).c3);
         _load_psi0((*sd).c3);
         _load_psi1_sub((*sd).c2);

         _load_chi((*rd).c4);
         _load_psi0((*sd).c4);
         _load_psi1_add((*sd).c1);         
         
         sd+=1;
      }

      _store_sm();
      acc_sm();
   }

   sum_sm();
   sp[0].re=sm0[0].re;
   sp[0].im=sm0[0].im;
   sp[1].re=-sm1[0].re;
   sp[1].im=-sm1[0].im;
}


static void spinor_prod_gamma3(int vol,spinor_dble *sd,spinor_dble *rd,
                               complex_dble *sp)
{
   spinor_dble *rt,*rm;

   init_sm();
   rt=rd+vol;
   rm=rd;

   while (rm<rt)
   {
      rm+=BLK_LENGTH;
      if (rm>rt)
         rm=rt;

      _start_sm();      

      for (;rd<rm;rd++)
      {
         _load_chi((*rd).c1);
         _load_psi0((*sd).c1);
         _load_psi1_sub((*sd).c3);

         _load_chi((*rd).c2);
         _load_psi0((*sd).c2);
         _load_psi1_add((*sd).c4);         

         _load_chi((*rd).c3);
         _load_psi0((*sd).c3);
         _load_psi1_add((*sd).c1);

         _load_chi((*rd).c4);
         _load_psi0((*sd).c4);
         _load_psi1_sub((*sd).c2);         
         
         sd+=1;
      }

      _store_sm();
      acc_sm();
   }

   sum_sm();
   sp[0].re=sm0[0].re;
   sp[0].im=sm0[0].im;
   sp[1].re=sm1[0].im;
   sp[1].im=-sm1[0].re;
}

#else

void apply_u2sd(int vol,int *imb,su3_dble *ud,spinor_dble *sd,
                spinor_dble *rd)
{
   int *imm;
   spinor_dble *si;

   imm=imb+vol;

   for (;imb<imm;imb++)
   {
      si=sd+(*imb);

      _su3_multiply((*rd).c1,(*ud),(*si).c1);
      _su3_multiply((*rd).c2,(*ud),(*si).c2);
      _su3_multiply((*rd).c3,(*ud),(*si).c3);
      _su3_multiply((*rd).c4,(*ud),(*si).c4);      

      ud+=1;
      rd+=1;
   }
}


void apply_udag2sd(int vol,int *imb,su3_dble *ud,spinor_dble *sd,
                   spinor_dble *rd)
{
   int *imm;
   spinor_dble *si;

   imm=imb+vol;

   for (;imb<imm;imb++)
   {
      si=sd+(*imb);

      _su3_inverse_multiply((*rd).c1,(*ud),(*si).c1);
      _su3_inverse_multiply((*rd).c2,(*ud),(*si).c2);
      _su3_inverse_multiply((*rd).c3,(*ud),(*si).c3);
      _su3_inverse_multiply((*rd).c4,(*ud),(*si).c4);      

      ud+=1;
      rd+=1;
   }
}  


static void spinor_prod_gamma0(int vol,spinor_dble *sd,spinor_dble *rd,
                               complex_dble *sp)
{
   complex_dble z0,z1;
   spinor_dble *rt,*rm;

   init_sm();
   rt=rd+vol;
   rm=rd;

   while (rm<rt)
   {
      rm+=BLK_LENGTH;
      if (rm>rt)
         rm=rt;

      z0.re=0.0;
      z0.im=0.0;
      z1.re=0.0;
      z1.im=0.0;
      
      for (;rd<rm;rd++)
      {
         z0.re+=(_vector_prod_re((*sd).c1,(*rd).c1)+
                 _vector_prod_re((*sd).c2,(*rd).c2)+
                 _vector_prod_re((*sd).c3,(*rd).c3)+
                 _vector_prod_re((*sd).c4,(*rd).c4));

         z0.im+=(_vector_prod_im((*sd).c1,(*rd).c1)+
                 _vector_prod_im((*sd).c2,(*rd).c2)+
                 _vector_prod_im((*sd).c3,(*rd).c3)+
                 _vector_prod_im((*sd).c4,(*rd).c4));      

         z1.re+=(_vector_prod_re((*sd).c1,(*rd).c3)+
                 _vector_prod_re((*sd).c2,(*rd).c4)+
                 _vector_prod_re((*sd).c3,(*rd).c1)+
                 _vector_prod_re((*sd).c4,(*rd).c2));               

         z1.im+=(_vector_prod_im((*sd).c1,(*rd).c3)+
                 _vector_prod_im((*sd).c2,(*rd).c4)+
                 _vector_prod_im((*sd).c3,(*rd).c1)+
                 _vector_prod_im((*sd).c4,(*rd).c2));               

         sd+=1;
      }

      sm0[0].re+=z0.re;
      sm0[0].im+=z0.im;
      sm1[0].re+=z1.re;
      sm1[0].im+=z1.im;      
      acc_sm();
   }

   sum_sm();
   sp[0].re=sm0[0].re;
   sp[0].im=sm0[0].im;
   sp[1].re=-sm1[0].re;
   sp[1].im=-sm1[0].im;
}


static void spinor_prod_gamma1(int vol,spinor_dble *sd,spinor_dble *rd,
                               complex_dble *sp)
{
   complex_dble z0,z1;
   spinor_dble *rt,*rm;

   init_sm();
   rt=rd+vol;
   rm=rd;

   while (rm<rt)
   {
      rm+=BLK_LENGTH;
      if (rm>rt)
         rm=rt;

      z0.re=0.0;
      z0.im=0.0;
      z1.re=0.0;
      z1.im=0.0;      
      
      for (;rd<rm;rd++)
      {
         z0.re+=(_vector_prod_re((*sd).c1,(*rd).c1)+
                 _vector_prod_re((*sd).c2,(*rd).c2)+
                 _vector_prod_re((*sd).c3,(*rd).c3)+
                 _vector_prod_re((*sd).c4,(*rd).c4));

         z0.im+=(_vector_prod_im((*sd).c1,(*rd).c1)+
                 _vector_prod_im((*sd).c2,(*rd).c2)+
                 _vector_prod_im((*sd).c3,(*rd).c3)+
                 _vector_prod_im((*sd).c4,(*rd).c4));      

         z1.re+=(_vector_prod_re((*sd).c1,(*rd).c4)+
                 _vector_prod_re((*sd).c2,(*rd).c3));

         z1.re-=(_vector_prod_re((*sd).c3,(*rd).c2)+
                 _vector_prod_re((*sd).c4,(*rd).c1));               

         z1.im+=(_vector_prod_im((*sd).c1,(*rd).c4)+
                 _vector_prod_im((*sd).c2,(*rd).c3));

         z1.im-=(_vector_prod_im((*sd).c3,(*rd).c2)+
                 _vector_prod_im((*sd).c4,(*rd).c1));               

         sd+=1;
      }

      sm0[0].re+=z0.re;
      sm0[0].im+=z0.im;
      sm1[0].re+=z1.re;
      sm1[0].im+=z1.im;      
      acc_sm();
   }

   sum_sm();
   sp[0].re=sm0[0].re;
   sp[0].im=sm0[0].im;
   sp[1].re=sm1[0].im;
   sp[1].im=-sm1[0].re;
}


static void spinor_prod_gamma2(int vol,spinor_dble *sd,spinor_dble *rd,
                               complex_dble *sp)
{
   complex_dble z0,z1;
   spinor_dble *rt,*rm;

   init_sm();
   rt=rd+vol;
   rm=rd;

   while (rm<rt)
   {
      rm+=BLK_LENGTH;
      if (rm>rt)
         rm=rt;

      z0.re=0.0;
      z0.im=0.0;
      z1.re=0.0;
      z1.im=0.0;      
      
      for (;rd<rm;rd++)
      {
         z0.re+=(_vector_prod_re((*sd).c1,(*rd).c1)+
                 _vector_prod_re((*sd).c2,(*rd).c2)+
                 _vector_prod_re((*sd).c3,(*rd).c3)+
                 _vector_prod_re((*sd).c4,(*rd).c4));

         z0.im+=(_vector_prod_im((*sd).c1,(*rd).c1)+
                 _vector_prod_im((*sd).c2,(*rd).c2)+
                 _vector_prod_im((*sd).c3,(*rd).c3)+
                 _vector_prod_im((*sd).c4,(*rd).c4));      

         z1.re+=(_vector_prod_re((*sd).c1,(*rd).c4)+
                 _vector_prod_re((*sd).c4,(*rd).c1));

         z1.re-=(_vector_prod_re((*sd).c2,(*rd).c3)+
                 _vector_prod_re((*sd).c3,(*rd).c2));               

         z1.im+=(_vector_prod_im((*sd).c1,(*rd).c4)+
                 _vector_prod_im((*sd).c4,(*rd).c1));

         z1.im-=(_vector_prod_im((*sd).c2,(*rd).c3)+
                 _vector_prod_im((*sd).c3,(*rd).c2));               

         sd+=1;
      }

      sm0[0].re+=z0.re;
      sm0[0].im+=z0.im;
      sm1[0].re+=z1.re;
      sm1[0].im+=z1.im;
      acc_sm();
   }

   sum_sm();
   sp[0].re=sm0[0].re;
   sp[0].im=sm0[0].im;
   sp[1].re=-sm1[0].re;
   sp[1].im=-sm1[0].im;
}


static void spinor_prod_gamma3(int vol,spinor_dble *sd,spinor_dble *rd,
                               complex_dble *sp)
{
   complex_dble z0,z1;
   spinor_dble *rt,*rm;

   init_sm();
   rt=rd+vol;
   rm=rd;

   while (rm<rt)
   {
      rm+=BLK_LENGTH;
      if (rm>rt)
         rm=rt;

      z0.re=0.0;
      z0.im=0.0;
      z1.re=0.0;
      z1.im=0.0;
      
      for (;rd<rm;rd++)
      {
         z0.re+=(_vector_prod_re((*sd).c1,(*rd).c1)+
                 _vector_prod_re((*sd).c2,(*rd).c2)+
                 _vector_prod_re((*sd).c3,(*rd).c3)+
                 _vector_prod_re((*sd).c4,(*rd).c4));

         z0.im+=(_vector_prod_im((*sd).c1,(*rd).c1)+
                 _vector_prod_im((*sd).c2,(*rd).c2)+
                 _vector_prod_im((*sd).c3,(*rd).c3)+
                 _vector_prod_im((*sd).c4,(*rd).c4));      

         z1.re+=(_vector_prod_re((*sd).c1,(*rd).c3)+
                 _vector_prod_re((*sd).c4,(*rd).c2));

         z1.re-=(_vector_prod_re((*sd).c2,(*rd).c4)+
                 _vector_prod_re((*sd).c3,(*rd).c1));               

         z1.im+=(_vector_prod_im((*sd).c1,(*rd).c3)+
                 _vector_prod_im((*sd).c4,(*rd).c2));

         z1.im-=(_vector_prod_im((*sd).c2,(*rd).c4)+
                 _vector_prod_im((*sd).c3,(*rd).c1));               

         sd+=1;
      }

      sm0[0].re+=z0.re;
      sm0[0].im+=z0.im;
      sm1[0].re+=z1.re;
      sm1[0].im+=z1.im;      
      acc_sm();
   }

   sum_sm();
   sp[0].re=sm0[0].re;
   sp[0].im=sm0[0].im;
   sp[1].re=sm1[0].im;
   sp[1].im=-sm1[0].re;
}

#endif

void (*spinor_prod_gamma[4])
(int vol,spinor_dble *sd,spinor_dble *rd,complex_dble *sp)=
{spinor_prod_gamma0,spinor_prod_gamma1,spinor_prod_gamma2,spinor_prod_gamma3};

