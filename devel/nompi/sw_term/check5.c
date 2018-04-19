
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of inv_pauli_dble()
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

static pauli_dble ma[3] ALIGNED16;
static complex_dble aa[4][36] ALIGNED16;


static void random_pauli(pauli_dble *m)
{
   int i;
   double *u;
 
   u=(*m).u;
   gauss_dble(u,36);

   for (i=0;i<6;i++)
      (*m).u[i]+=10.0;
}


static void pauli2mat(pauli_dble *m,complex_dble *a)
{
   int i,j,k;
   double *u;

   u=(*m).u;
   k=6;

   for (i=0;i<6;i++)
   {
      a[6*i+i].re=u[i];
      a[6*i+i].im=0.0;

      for (j=i+1;j<6;j++)
      {
         a[6*i+j].re=u[k];
         a[6*j+i].re=u[k];         
         k+=1;
         a[6*i+j].im=u[k];
         a[6*j+i].im=-u[k];         
         k+=1;
      }
   }
}


int main(void)
{
   int i,j,ie;
   double mu,d,dmax;

   printf("\n");
   printf("Check of inv_pauli_dble()\n");
   printf("-------------------------\n\n");

#if (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   rlxd_init(1,3898);
   mu=0.1234;
   ie=1;

   while (ie)
   {
      random_pauli(ma);
      ie=inv_pauli_dble(mu,ma,ma+1);
   }
   
   pauli2mat(ma,aa[0]);
   pauli2mat(ma+1,aa[1]);   
   cmat_mul_dble(6,aa[0],aa[0],aa[2]);

   for (i=0;i<6;i++)
      aa[2][6*i+i].re+=mu*mu;

   cmat_mul_dble(6,aa[1],aa[2],aa[3]);
   cmat_sub_dble(6,aa[3],aa[0],aa[2]);
   dmax=0.0;

   for (i=0;i<6;i++)
   {
      for (j=0;j<6;j++)
      {
         d=aa[2][6*i+j].re;

         if (d<0.0)
            d=-d;
         if (d>dmax)
            dmax=d;

         d=aa[2][6*i+j].im;

         if (d<0.0)
            d=-d;
         if (d>dmax)
            dmax=d;        
      }
   }

   printf("Maximal absolute deviation = %.1e\n",dmax);

   inv_pauli_dble(mu,ma,ma);
   dmax=0.0;

   for (i=0;i<36;i++)
   {
      d=ma[0].u[i]-ma[1].u[i];

      if (d<0.0)
         d=-d;
      if (d>dmax)
         dmax=d;
   }
   
   error(dmax!=0.0,1,"main [check5.c]",
         "inv_pauli_dble() is incorrect when m=im");    
   printf("Works correctly if input and output matrices coincide\n\n");
   
   exit(0);
}
