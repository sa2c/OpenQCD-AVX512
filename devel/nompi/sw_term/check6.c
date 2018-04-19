
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2005, 2009, 2010, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Accuracy of inv_pauli_dble()
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "utils.h"
#include "random.h"
#include "linalg.h"
#include "sw_term.h"

#define NM 10000

typedef union
{
   weyl w;
   complex c[6];
} spin_t;

typedef union
{
   weyl_dble w;
   complex_dble c[6];
} spin_dble_t;

static spin_t vs ALIGNED16;
static spin_dble_t vd ALIGNED16;
static const weyl vs0={{{0.0f,0.0f},{0.0f,0.0f},{0.0f,0.0f}},
                       {{0.0f,0.0f},{0.0f,0.0f},{0.0f,0.0f}}};
static const weyl_dble vd0={{{0.0,0.0},{0.0,0.0},{0.0,0.0}},
                            {{0.0,0.0},{0.0,0.0},{0.0,0.0}}};


int main(void)
{
   int n,k,l,itot,*is;
   double mu,fact,d,dmax;
   pauli *ms,*ims,*msb,*imsb;
   pauli_dble *md,*imd,*mdb,*imdb;

   printf("\n");
   printf("Accuracy of inv_pauli_dble()\n");
   printf("----------------------------\n\n");

#if (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   is=amalloc(NM*sizeof(*is),3);
   msb=amalloc(3*NM*sizeof(*msb),4);
   mdb=amalloc(3*NM*sizeof(*mdb),4);
   error((is==NULL)||(msb==NULL)||(mdb==NULL),1,
         "main [check6.c]","Unable to allocate auxiliary arrays");

   imsb=msb+NM;
   imdb=mdb+NM;

   rlxd_init(1,1234);
   mu=0.0123;
   md=mdb;
   imd=imdb;
   itot=0;
   dmax=0.0;
   fact=sqrt(2.0);

   for (n=0;n<NM;n++)
   {
      gauss_dble((*md).u,36);

      for (k=0;k<6;k++)
         (*md).u[k]*=fact;

      is[n]=(inv_pauli_dble(0.0,md,imd)|inv_pauli_dble(mu,md,imd+NM));

      if (is[n]==0)
      {
         for (k=0;k<6;k++)
         {
            vd.w=vd0;
            vd.c[k].re=1.0;

            mul_pauli_dble(mu,md,&(vd.w),&(vd.w));
            mul_pauli_dble(0.0,imd,&(vd.w),&(vd.w));
            mul_pauli_dble(-mu,md,&(vd.w),&(vd.w));
            mul_pauli_dble(0.0,imd+NM,&(vd.w),&(vd.w));

            vd.c[k].re-=1.0;

            for (l=0;l<6;l++)
            {
               d=vd.c[l].re*vd.c[l].re+vd.c[l].im*vd.c[l].im;
               if (d>dmax)
                  dmax=d;
            }
         }
      }
      else
         itot+=1;

      md+=1;
      imd+=1;
   }

   printf("Double-precision program, mu=%.4f:\n",mu);
   printf("%d Gaussian random matrices, %d inversion failures\n",NM,itot);
   printf("Maximal relative deviation = %.1e ",sqrt(dmax));
   printf("(safe cases only)\n\n");

   assign_pauli(NM,mdb,msb);
   assign_pauli(2*NM,imdb,imsb);

   ms=msb;
   ims=imsb;
   dmax=0.0;

   for (n=0;n<NM;n++)
   {
      if (is[n]==0)
      {
         for (k=0;k<6;k++)
         {
            vs.w=vs0;
            vs.c[k].re=1.0f;

            mul_pauli((float)(mu),ms,&(vs.w),&(vs.w));
            mul_pauli(0.0f,ims,&(vs.w),&(vs.w));
            mul_pauli((float)(-mu),ms,&(vs.w),&(vs.w));
            mul_pauli(0.0f,ims+NM,&(vs.w),&(vs.w));

            vs.c[k].re-=1.0f;

            for (l=0;l<6;l++)
            {
               d=(double)(vs.c[l].re*vs.c[l].re+vs.c[l].im*vs.c[l].im);
               if (d>dmax)
                  dmax=d;
            }
         }
      }

      ms+=1;
      ims+=1;
   }

   printf("After assignment to single-precision matrices:\n");
   printf("Maximal relative deviation = %.1e ",sqrt(dmax));
   printf("(safe cases only)\n\n");
   exit(0);
}
