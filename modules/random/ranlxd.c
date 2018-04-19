
/*******************************************************************************
*
* File ranlxd.c
*
* Copyright (C) 2005, 2008, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Random number generator "ranlxd" version 3.0. See the notes
*
*   "User's guide for ranlxs and ranlxd v3.0" (May 2001)
*
*   "Algorithms used in ranlux v3.0" (May 2001)
*
* for a detailed description.
*
* The externally accessible functions are
*
*   void ranlxd(double r[],int n)
*     Computes the next n double-precision random numbers and
*     assigns them to the elements r[0],...,r[n-1] of the array r[]
*
*   void rlxd_init(int level,int seed)
*     Initialization of the generator
*
*   int rlxd_size(void)
*     Returns the number of integers required to save the state of
*     the generator
*
*   void rlxd_get(int state[])
*     Extracts the current state of the generator and stores the
*     information in the array state[N] where N>=rlxd_size()
*
*   void rlxd_reset(int state[])
*     Resets the generator to the state defined by the array state[N]
*
*******************************************************************************/

#define RANLXD_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "random.h"

#if ((defined AVX)||(defined x64))

typedef struct __attribute__ ((aligned (16)))
{
   float c1,c2,c3,c4;
} vec_t;

typedef struct
{
   vec_t c1,c2;
} dble_vec_t;

static int init=0,pr,prm,ir,jr,is,is_old,next[96];
static vec_t one,one_bit,carry;

static union
{
   dble_vec_t vec[12];
   float num[96];
} x;

#define STEP(pi,pj) \
  __asm__ __volatile__ ("movaps %4, %%xmm4 \n\t" \
                        "movaps %%xmm2, %%xmm3 \n\t" \
                        "subps %2, %%xmm4 \n\t" \
                        "movaps %%xmm1, %%xmm5 \n\t" \
                        "cmpps $0x6, %%xmm4, %%xmm2 \n\t" \
                        "andps %%xmm2, %%xmm5 \n\t" \
                        "subps %%xmm3, %%xmm4 \n\t" \
                        "andps %%xmm0, %%xmm2 \n\t" \
                        "addps %%xmm4, %%xmm5 \n\t" \
                        "movaps %%xmm5, %0 \n\t" \
                        "movaps %5, %%xmm6 \n\t" \
                        "movaps %%xmm2, %%xmm3 \n\t" \
                        "subps %3, %%xmm6 \n\t" \
                        "movaps %%xmm1, %%xmm7 \n\t" \
                        "cmpps $0x6, %%xmm6, %%xmm2 \n\t" \
                        "andps %%xmm2, %%xmm7 \n\t" \
                        "subps %%xmm3, %%xmm6 \n\t" \
                        "andps %%xmm0, %%xmm2 \n\t" \
                        "addps %%xmm6, %%xmm7 \n\t" \
                        "movaps %%xmm7, %1" \
                        : \
                        "=m" ((*pi).c1), \
                        "=m" ((*pi).c2) \
                        : \
                        "m" ((*pi).c1), \
                        "m" ((*pi).c2), \
                        "m" ((*pj).c1), \
                        "m" ((*pj).c2) \
                        : \
                        "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7")


static void update(void)
{
   int k,kmax;
   dble_vec_t *pmin,*pmax,*pi,*pj;

   kmax=pr;
   pmin=&x.vec[0];
   pmax=pmin+12;
   pi=&x.vec[ir];
   pj=&x.vec[jr];

   __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                         "movaps %1, %%xmm1 \n\t"
                         "movaps %2, %%xmm2"
                         :
                         :
                         "m" (one_bit),
                         "m" (one),
                         "m" (carry)
                         :
                         "xmm0", "xmm1", "xmm2");

   for (k=0;k<kmax;k++)
   {
      STEP(pi,pj);
      pi+=1;
      pj+=1;
      if (pi==pmax)
         pi=pmin;
      if (pj==pmax)
         pj=pmin;
   }

   __asm__ __volatile__ ("movaps %%xmm2, %0"
                         :
                         "=m" (carry));

   ir+=prm;
   jr+=prm;
   if (ir>=12)
      ir-=12;
   if (jr>=12)
      jr-=12;
   is=8*ir;
   is_old=is;
}


static void define_constants(void)
{
   int k;
   float b;

   one.c1=1.0f;
   one.c2=1.0f;
   one.c3=1.0f;
   one.c4=1.0f;

   b=(float)(ldexp(1.0,-24));
   one_bit.c1=b;
   one_bit.c2=b;
   one_bit.c3=b;
   one_bit.c4=b;

   for (k=0;k<96;k++)
   {
      next[k]=(k+1)%96;
      if ((k%4)==3)
         next[k]=(k+5)%96;
   }
}


void rlxd_init(int level,int seed)
{
   int i,k,l;
   int ibit,jbit,xbit[31];
   int ix,iy;

   define_constants();

   error_loc((level<1)||(level>2),1,"rlxd_init [ranlxd.c]",
             "Bad choice of luxury level (should be 1 or 2)");

   if (level==1)
      pr=202;
   else if (level==2)
      pr=397;

   i=seed;

   for (k=0;k<31;k++)
   {
      xbit[k]=i%2;
      i/=2;
   }

   error_loc((seed<=0)||(i!=0),1,"rlxd_init [ranlxd.c]",
             "Bad choice of seed (should be between 1 and 2^31-1)");

   ibit=0;
   jbit=18;

   for (i=0;i<4;i++)
   {
      for (k=0;k<24;k++)
      {
         ix=0;

         for (l=0;l<24;l++)
         {
            iy=xbit[ibit];
            ix=2*ix+iy;

            xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
            ibit=(ibit+1)%31;
            jbit=(jbit+1)%31;
         }

         if ((k%4)!=i)
            ix=16777215-ix;

         x.num[4*k+i]=(float)(ldexp((double)(ix),-24));
      }
   }

   carry.c1=0.0f;
   carry.c2=0.0f;
   carry.c3=0.0f;
   carry.c4=0.0f;

   ir=0;
   jr=7;
   is=91;
   is_old=0;
   prm=pr%12;
   init=1;
}


void ranlxd(double r[],int n)
{
   int k;

   if (init==0)
      rlxd_init(1,1);

   for (k=0;k<n;k++)
   {
      is=next[is];
      if (is==is_old)
         update();
      r[k]=(double)(x.num[is+4])+(double)(one_bit.c1*x.num[is]);
   }
}


int rlxd_size(void)
{
   return(105);
}


void rlxd_get(int state[])
{
   int k;
   float base;

   error_loc(init==0,1,"rlxd_get [ranlxd.c]",
             "Undefined state (ranlxd is not initialized)");

   base=(float)(ldexp(1.0,24));
   state[0]=rlxd_size();

   for (k=0;k<96;k++)
      state[k+1]=(int)(base*x.num[k]);

   state[97]=(int)(base*carry.c1);
   state[98]=(int)(base*carry.c2);
   state[99]=(int)(base*carry.c3);
   state[100]=(int)(base*carry.c4);

   state[101]=pr;
   state[102]=ir;
   state[103]=jr;
   state[104]=is;
}


void rlxd_reset(int state[])
{
   int k;

   define_constants();

   error_loc(state[0]!=rlxd_size(),1,"rlxd_reset [ranlxd.c]",
             "Unexpected input data");

   for (k=0;k<96;k++)
   {
      error_loc((state[k+1]<0)||(state[k+1]>=167777216),1,
                "rlxd_reset [ranlxd.c]","Unexpected input data");

      x.num[k]=(float)(ldexp((double)(state[k+1]),-24));
   }

   error_loc(((state[97]!=0)&&(state[97]!=1))||
             ((state[98]!=0)&&(state[98]!=1))||
             ((state[99]!=0)&&(state[99]!=1))||
             ((state[100]!=0)&&(state[100]!=1)),1,
             "rlxd_reset [ranlxd.c]","Unexpected input data");

   carry.c1=(float)(ldexp((double)(state[97]),-24));
   carry.c2=(float)(ldexp((double)(state[98]),-24));
   carry.c3=(float)(ldexp((double)(state[99]),-24));
   carry.c4=(float)(ldexp((double)(state[100]),-24));

   pr=state[101];
   ir=state[102];
   jr=state[103];
   is=state[104];
   is_old=8*ir;
   prm=pr%12;
   init=1;

   error_loc(((pr!=202)&&(pr!=397))||
             (ir<0)||(ir>11)||(jr<0)||(jr>11)||(jr!=((ir+7)%12))||
             (is<0)||(is>91),1,
             "rlxd_reset [ranlxd.c]","Unexpected input data");
}

#else

#define BASE 0x1000000
#define MASK 0xffffff

typedef struct
{
   int c1,c2,c3,c4;
} vec_t;

typedef struct
{
   vec_t c1,c2;
} dble_vec_t;

static int init=0,pr,prm,ir,jr,is,is_old,next[96];
static double one_bit;
static vec_t carry;

static union
{
   dble_vec_t vec[12];
   int num[96];
} x;

#define STEP(pi,pj) \
      d=(*pj).c1.c1-(*pi).c1.c1-carry.c1; \
      (*pi).c2.c1+=(d<0); \
      d+=BASE; \
      (*pi).c1.c1=d&MASK; \
      d=(*pj).c1.c2-(*pi).c1.c2-carry.c2; \
      (*pi).c2.c2+=(d<0); \
      d+=BASE; \
      (*pi).c1.c2=d&MASK; \
      d=(*pj).c1.c3-(*pi).c1.c3-carry.c3; \
      (*pi).c2.c3+=(d<0); \
      d+=BASE; \
      (*pi).c1.c3=d&MASK; \
      d=(*pj).c1.c4-(*pi).c1.c4-carry.c4; \
      (*pi).c2.c4+=(d<0); \
      d+=BASE; \
      (*pi).c1.c4=d&MASK; \
      d=(*pj).c2.c1-(*pi).c2.c1; \
      carry.c1=(d<0); \
      d+=BASE; \
      (*pi).c2.c1=d&MASK; \
      d=(*pj).c2.c2-(*pi).c2.c2; \
      carry.c2=(d<0); \
      d+=BASE; \
      (*pi).c2.c2=d&MASK; \
      d=(*pj).c2.c3-(*pi).c2.c3; \
      carry.c3=(d<0); \
      d+=BASE; \
      (*pi).c2.c3=d&MASK; \
      d=(*pj).c2.c4-(*pi).c2.c4; \
      carry.c4=(d<0); \
      d+=BASE; \
      (*pi).c2.c4=d&MASK


static void update(void)
{
   int k,kmax,d;
   dble_vec_t *pmin,*pmax,*pi,*pj;

   kmax=pr;
   pmin=&x.vec[0];
   pmax=pmin+12;
   pi=&x.vec[ir];
   pj=&x.vec[jr];

   for (k=0;k<kmax;k++)
   {
      STEP(pi,pj);
      pi+=1;
      pj+=1;
      if (pi==pmax)
         pi=pmin;
      if (pj==pmax)
         pj=pmin;
   }

   ir+=prm;
   jr+=prm;
   if (ir>=12)
      ir-=12;
   if (jr>=12)
      jr-=12;
   is=8*ir;
   is_old=is;
}


static void define_constants(void)
{
   int k;

   one_bit=ldexp(1.0,-24);

   for (k=0;k<96;k++)
   {
      next[k]=(k+1)%96;
      if ((k%4)==3)
         next[k]=(k+5)%96;
   }
}


void rlxd_init(int level,int seed)
{
   int i,k,l;
   int ibit,jbit,xbit[31];
   int ix,iy;

   error_loc((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24)||
             (DBL_MANT_DIG<48),1,"rlxd_init [ranlxd.c]",
             "Arithmetic on this machine is not suitable for ranlxd");

   define_constants();

   error_loc((level<1)||(level>2),1,"rlxd_init [ranlxd.c]",
             "Bad choice of luxury level (should be 1 or 2)");

   if (level==1)
      pr=202;
   else if (level==2)
      pr=397;

   i=seed;

   for (k=0;k<31;k++)
   {
      xbit[k]=i%2;
      i/=2;
   }

   error_loc((seed<=0)||(i!=0),1,"rlxd_init [ranlxd.c]",
             "Bad choice of seed (should be between 1 and 2^31-1)");

   ibit=0;
   jbit=18;

   for (i=0;i<4;i++)
   {
      for (k=0;k<24;k++)
      {
         ix=0;

         for (l=0;l<24;l++)
         {
            iy=xbit[ibit];
            ix=2*ix+iy;

            xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
            ibit=(ibit+1)%31;
            jbit=(jbit+1)%31;
         }

         if ((k%4)!=i)
            ix=16777215-ix;

         x.num[4*k+i]=ix;
      }
   }

   carry.c1=0;
   carry.c2=0;
   carry.c3=0;
   carry.c4=0;

   ir=0;
   jr=7;
   is=91;
   is_old=0;
   prm=pr%12;
   init=1;
}


void ranlxd(double r[],int n)
{
   int k;

   if (init==0)
      rlxd_init(1,1);

   for (k=0;k<n;k++)
   {
      is=next[is];
      if (is==is_old)
         update();
      r[k]=one_bit*((double)(x.num[is+4])+one_bit*(double)(x.num[is]));
   }
}


int rlxd_size(void)
{
   return(105);
}


void rlxd_get(int state[])
{
   int k;

   error_loc(init==0,1,"rlxd_get [ranlxd.c]",
             "Undefined state (ranlxd is not initialized)");

   state[0]=rlxd_size();

   for (k=0;k<96;k++)
      state[k+1]=x.num[k];

   state[97]=carry.c1;
   state[98]=carry.c2;
   state[99]=carry.c3;
   state[100]=carry.c4;

   state[101]=pr;
   state[102]=ir;
   state[103]=jr;
   state[104]=is;
}


void rlxd_reset(int state[])
{
   int k;

   error_loc((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24)||
             (DBL_MANT_DIG<48),1,"rlxd_init [ranlxd.c]",
             "Arithmetic on this machine is not suitable for ranlxd");


   define_constants();

   error_loc(state[0]!=rlxd_size(),1,"rlxd_reset [ranlxd.c]",
             "Unexpected input data");

   for (k=0;k<96;k++)
   {
      error_loc((state[k+1]<0)||(state[k+1]>=167777216),1,
                "rlxd_reset [ranlxd.c]","Unexpected input data");

      x.num[k]=state[k+1];
   }

   error_loc(((state[97]!=0)&&(state[97]!=1))||
             ((state[98]!=0)&&(state[98]!=1))||
             ((state[99]!=0)&&(state[99]!=1))||
             ((state[100]!=0)&&(state[100]!=1)),1,
             "rlxd_reset [ranlxd.c]","Unexpected input data");

   carry.c1=state[97];
   carry.c2=state[98];
   carry.c3=state[99];
   carry.c4=state[100];

   pr=state[101];
   ir=state[102];
   jr=state[103];
   is=state[104];
   is_old=8*ir;
   prm=pr%12;
   init=1;

   error_loc(((pr!=202)&&(pr!=397))||
             (ir<0)||(ir>11)||(jr<0)||(jr>11)||(jr!=((ir+7)%12))||
             (is<0)||(is>91),1,
             "rlxd_reset [ranlxd.c]","Unexpected input data");
}

#endif
