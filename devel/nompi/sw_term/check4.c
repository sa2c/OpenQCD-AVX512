
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of apply_sw_dble() and apply_swinv_dble()
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

#define NM 1001

typedef union
{
   spinor_dble s;
   complex_dble c[12];
   double r[24];
} spin_t;

static pauli_dble m[2*NM] ALIGNED16;
static spin_t sp1[NM] ALIGNED16;
static spin_t sp2[NM] ALIGNED16;
static spin_t rp1[NM] ALIGNED16;
static spin_t rp2[NM] ALIGNED16;
static complex_dble mv[36] ALIGNED16;


static void random_pauli_dble(void)
{
   int i,j;
   double *u;

   for (i=0;i<(2*NM);i++)
   {
      u=m[i].u;
      gauss_dble(u,36);

      for (j=0;j<6;j++)
         u[j]+=10.0;
   }
}


static void random_spin(void)
{
   int i;

   for (i=0;i<NM;i++)
      gauss_dble(sp1[i].r,24);
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


static double diff_spin(spin_t *sp,spin_t *rp)
{
   int i,j;
   float d,dmax;

   dmax=0.0;

   for (i=0;i<NM;i++)
   {
      for (j=0;j<24;j++)
      {
         d=sp[i].r[j]-rp[i].r[j];
         if (d<0.0)
            d=-d;
         if (d>dmax)
            dmax=d;
      }
   }

   return dmax;
}


static void pauli2mv(double mu,pauli_dble *mp)
{
   int i,j,k;
   double *u;

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
   int i,ie;
   double mu;
   spinor_dble *s1,*r1;

   printf("\n");
   printf("Check of apply_sw_dble() and apply_swinv_dble()\n");
   printf("-----------------------------------------------\n\n");

#if (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   rlxd_init(1,3898);
   s1=(spinor_dble*)(sp1);
   r1=(spinor_dble*)(rp1);
   mu=0.0123;

   random_pauli_dble();
   random_spin();
   cp_spin(sp1,sp2);
   apply_sw_dble(NM,mu,m,s1,r1);

   error(diff_spin(sp1,sp2)!=0.0,1,"main [check4.c]",
         "apply_sw_dble() does not preserve the input spinor field");

   for (i=0;i<NM;i++)
   {
      pauli2mv(mu,m+2*i);
      cmat_vec_dble(6,mv,sp1[i].c,rp2[i].c);

      pauli2mv(-mu,m+2*i+1);
      cmat_vec_dble(6,mv,sp1[i].c+6,rp2[i].c+6);
   }

   printf("Check of apply_sw_dble():\nAbsolute deviation = %.1e\n",
          diff_spin(rp1,rp2));
   printf("Input spinor field preserved\n");

   apply_sw_dble(NM,mu,m,s1,s1);
   error(diff_spin(sp1,rp1)!=0.0f,1,"main [check4.c]",
         "apply_sw_dble() does not work correctly if r=s");
   printf("Works correctly if input and output fields coincide\n\n");

   random_pauli_dble();
   random_spin();
   cp_spin(sp1,sp2);
   ie=apply_swinv_dble(NM,mu,m,s1,r1);

   error(diff_spin(sp1,sp2)!=0.0,1,"main [check4.c]",
         "apply_swinv_dble() does not preserve the input spinor field");

   for (i=0;i<NM;i++)
   {
      pauli2mv(mu,m+2*i);
      cmat_vec_dble(6,mv,rp1[i].c,rp2[i].c);

      pauli2mv(-mu,m+2*i+1);
      cmat_vec_dble(6,mv,rp1[i].c+6,rp2[i].c+6);
   }

   error(ie!=0,1,"main [check4.c]",
         "apply_swinv_dble(): not all matrix inversions were safe");
   printf("Check of apply_swinv_dble() (safe inversions only):\n");
   printf("Absolute deviation = %.1e\n",diff_spin(sp1,rp2));
   printf("Input spinor field preserved\n");

   apply_swinv_dble(NM,mu,m,s1,s1);
   error(diff_spin(sp1,rp1)!=0.0f,1,"main [check4.c]",
         "apply_swinv_dble() does not work correctly if r=s");
   printf("Works correctly if input and output fields coincide\n\n");

   exit(0);
}
