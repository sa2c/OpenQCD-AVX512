
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of assign_pauli() and apply_sw()
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

#define NM 131

typedef union
{
   spinor s;
   complex c[12];
   float r[24];
} spin_t;

static pauli m[2*NM] ALIGNED16;
static pauli_dble md[2*NM] ALIGNED16;
static spin_t sp1[NM] ALIGNED16;
static spin_t sp2[NM] ALIGNED16;
static spin_t rp1[NM] ALIGNED16;
static spin_t rp2[NM] ALIGNED16;
static complex mv[36] ALIGNED16;


static void random_pauli_dble(void)
{
   int i;
   double *u;

   for (i=0;i<(2*NM);i++)
   {
      u=md[i].u;
      gauss_dble(u,36);
   }
}


static float diff_pauli(void)
{
   int i,j;
   float d,dmax,*u;
   double *ud;

   dmax=0.0f;

   for (i=0;i<(2*NM);i++)
   {
      u=m[i].u;
      ud=md[i].u;

      for (j=0;j<36;j++)
      {
         d=u[j]-(float)(ud[j]);
         if (d<0.0f)
            d=-d;
         if (d>dmax)
            dmax=d;
      }
   }

   return dmax;
}


static void random_spin(void)
{
   int i;

   for (i=0;i<NM;i++)
      gauss(sp1[i].r,24);
}


static void cp_spin(spin_t *sp,spin_t *rp)
{
   int i,j;

   for (i=0;i<NM;i++)
   {
      for (j=0;j<24;j++)
         rp[i].r[j]=sp[i].r[j];
   }
}


static float diff_spin(spin_t *sp,spin_t *rp)
{
   int i,j;
   float d,dmax;

   dmax=0.0f;

   for (i=0;i<NM;i++)
   {
      for (j=0;j<24;j++)
      {
         d=sp[i].r[j]-rp[i].r[j];
         if (d<0.0f)
            d=-d;
         if (d>dmax)
            dmax=d;
      }
   }

   return dmax;
}


static void pauli2mv(float mu,pauli *mp)
{
   int i,j,k;
   float *u;

   u=(*mp).u;
   k=6;

   for (i=0;i<6;i++)
   {
      mv[6*i+i].re=u[i];
      mv[6*i+i].im=mu;

      for (j=i+1;j<6;j++)
      {
         mv[6*i+j].re=u[k];
         mv[6*j+i].re=u[k];
         k+=1;
         mv[6*i+j].im=u[k];
         mv[6*j+i].im=-u[k];
         k+=1;
      }
   }
}


int main(void)
{
   int i;
   float mu;
   spinor *s1,*r1;

   printf("\n");
   printf("Check of assign_pauli() and apply_sw()\n");
   printf("--------------------------------------\n\n");

#if (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   rlxs_init(0,3898);
   random_pauli_dble();
   assign_pauli(2*NM,md,m);

   printf("Check of assign_pauli():\nAbsolute deviation = %.1e\n\n",
          diff_pauli());

   random_spin();
   cp_spin(sp1,sp2);
   mu=0.1234f;

   s1=(spinor*)(sp1);
   r1=(spinor*)(rp1);
   apply_sw(NM,mu,m,s1,r1);

   error(diff_spin(sp1,sp2)!=0.0f,1,"main [check3.c]",
         "apply_sw() does not preserve the input spinor field");

   for (i=0;i<NM;i++)
   {
      pauli2mv(mu,m+2*i);
      cmat_vec(6,mv,sp1[i].c,rp2[i].c);

      pauli2mv(-mu,m+2*i+1);
      cmat_vec(6,mv,sp1[i].c+6,rp2[i].c+6);
   }

   printf("Check of apply_sw():\nAbsolute deviation = %.1e\n",
          diff_spin(rp1,rp2));
   printf("Input spinor field preserved\n");

   apply_sw(NM,mu,m,s1,s1);
   error(diff_spin(sp1,rp1)!=0.0f,1,"main [check3.c]",
         "apply_sw() does not work correctly if r=s");
   printf("Works correctly if input and output fields coincide\n\n");

   exit(0);
}
