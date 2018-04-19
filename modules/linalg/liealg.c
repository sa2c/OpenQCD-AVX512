
/*******************************************************************************
*
* File liealg.c
*
* Copyright (C) 2005, 2009-2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Basic functions for fields with values in the Lie algebra of SU(3)
*
* The externally accessible functions are
*
*   void random_alg(int vol,su3_alg_dble *X)
*     Initializes the Lie algebra elements X to random values
*     with distribution proportional to exp{tr[X^2]}.
*
*   double norm_square_alg(int vol,int icom,su3_alg_dble *X)
*     Computes the square of the norm of the norm squared of the field X.
*
*   double scalar_prod_alg(int vol,int icom,su3_alg_dble *X,su3_alg_dble *Y)
*     Computes the scalar product of the fields X and Y.
*
*   void set_alg2zero(int vol,su3_alg_dble *X)
*     Sets the array elements X to zero.
*
*   void set_ualg2zero(int vol,u3_alg_dble *X)
*     Sets the array elements X to zero.
*
*   void assign_alg2alg(int vol,su3_alg_dble *X,su3_alg_dble *Y)
*     Assigns the field X to the field Y.
*
*   void swap_alg(int vol,su3_alg_dble *X,su3_alg_dble *Y)
*     Swaps the fields X and Y.
*
*   void muladd_assign_alg(int vol,double r,su3_alg_dble *X,su3_alg_dble *Y)
*     Adds r*X to Y.
*
* Notes:
*
* Lie algebra elements X are traceless antihermitian 3x3 matrices that
* are represented by structures with real elements x1,...,x8 through
*
*  X_11=i*(x1+x2), X_22=i*(x2-2*x1), X_33=i*(x1-2*x2),
*
*  X_12=x3+i*x4, X_13=x5+i*x6, X_23=x7+i*x8
*
* The scalar product (X,Y) of any two elements of the Lie algebra is
*
*  (X,Y)=-2*tr{XY}
*
* and the norm of X is (X,X)^(1/2).
*
* All programs in this module operate on arrays of Lie algebra elements whose
* base address is passed through the arguments. The length of the array is
* specified by the parameter vol. Scalar products etc. are globally summed if
* the parameter icom is equal to 1. In this case the calculated values are
* guaranteed to be exactly the same on all processes.
*
*******************************************************************************/

#define LIEALG_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "random.h"
#include "linalg.h"
#include "global.h"

static int ism,init=0;
static double c1=0.0,c2,c3,rb[8];


void random_alg(int vol,su3_alg_dble *X)
{
   su3_alg_dble *Xm;

   if (c1==0.0)
   {
      c1=(sqrt(3.0)+1.0)/6.0;
      c2=(sqrt(3.0)-1.0)/6.0;
      c3=1.0/sqrt(2.0);
   }

   Xm=X+vol;

   for (;X<Xm;X++)
   {
      gauss_dble(rb,8);

      (*X).c1=c1*rb[0]+c2*rb[1];
      (*X).c2=c1*rb[1]+c2*rb[0];
      (*X).c3=c3*rb[2];
      (*X).c4=c3*rb[3];
      (*X).c5=c3*rb[4];
      (*X).c6=c3*rb[5];
      (*X).c7=c3*rb[6];
      (*X).c8=c3*rb[7];
   }
}


double norm_square_alg(int vol,int icom,su3_alg_dble *X)
{
   double sm;
   su3_alg_dble *Xm;

   if (init==0)
   {
      ism=init_hsum(1);
      init=1;
   }

   reset_hsum(ism);
   Xm=X+vol;

   for (;X<Xm;X++)
   {
      sm=3.0*((*X).c1*(*X).c1+(*X).c2*(*X).c2-(*X).c1*(*X).c2)+
         (*X).c3*(*X).c3+(*X).c4*(*X).c4+(*X).c5*(*X).c5+
         (*X).c6*(*X).c6+(*X).c7*(*X).c7+(*X).c8*(*X).c8;

      add_to_hsum(ism,&sm);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(ism,&sm);
   else
      local_hsum(ism,&sm);

   return 4.0*sm;
}


double scalar_prod_alg(int vol,int icom,su3_alg_dble *X,su3_alg_dble *Y)
{
   double sm;
   su3_alg_dble *Xm;

   if (init==0)
   {
      ism=init_hsum(1);
      init=1;
   }

   reset_hsum(ism);
   Xm=X+vol;

   for (;X<Xm;X++)
   {
      sm=12.0*((*X).c1*(*Y).c1+(*X).c2*(*Y).c2)
         -6.0*((*X).c1*(*Y).c2+(*X).c2*(*Y).c1)
         +4.0*((*X).c3*(*Y).c3+(*X).c4*(*Y).c4+(*X).c5*(*Y).c5+
               (*X).c6*(*Y).c6+(*X).c7*(*Y).c7+(*X).c8*(*Y).c8);

      Y+=1;
      add_to_hsum(ism,&sm);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(ism,&sm);
   else
      local_hsum(ism,&sm);

   return sm;
}


void set_alg2zero(int vol,su3_alg_dble *X)
{
   su3_alg_dble *Xm;

   Xm=X+vol;

   for (;X<Xm;X++)
   {
      (*X).c1=0.0;
      (*X).c2=0.0;
      (*X).c3=0.0;
      (*X).c4=0.0;
      (*X).c5=0.0;
      (*X).c6=0.0;
      (*X).c7=0.0;
      (*X).c8=0.0;
   }
}


void set_ualg2zero(int vol,u3_alg_dble *X)
{
   u3_alg_dble *Xm;

   Xm=X+vol;

   for (;X<Xm;X++)
   {
      (*X).c1=0.0;
      (*X).c2=0.0;
      (*X).c3=0.0;
      (*X).c4=0.0;
      (*X).c5=0.0;
      (*X).c6=0.0;
      (*X).c7=0.0;
      (*X).c8=0.0;
      (*X).c9=0.0;
   }
}


void assign_alg2alg(int vol,su3_alg_dble *X,su3_alg_dble *Y)
{
   su3_alg_dble *Xm;

   Xm=X+vol;

   for (;X<Xm;X++)
   {
      (*Y).c1=(*X).c1;
      (*Y).c2=(*X).c2;
      (*Y).c3=(*X).c3;
      (*Y).c4=(*X).c4;
      (*Y).c5=(*X).c5;
      (*Y).c6=(*X).c6;
      (*Y).c7=(*X).c7;
      (*Y).c8=(*X).c8;

      Y+=1;
   }
}


void swap_alg(int vol,su3_alg_dble *X,su3_alg_dble *Y)
{
   double r;
   su3_alg_dble *Xm;

   Xm=X+vol;

   for (;X<Xm;X++)
   {
      r=(*Y).c1;
      (*Y).c1=(*X).c1;
      (*X).c1=r;

      r=(*Y).c2;
      (*Y).c2=(*X).c2;
      (*X).c2=r;

      r=(*Y).c3;
      (*Y).c3=(*X).c3;
      (*X).c3=r;

      r=(*Y).c4;
      (*Y).c4=(*X).c4;
      (*X).c4=r;

      r=(*Y).c5;
      (*Y).c5=(*X).c5;
      (*X).c5=r;

      r=(*Y).c6;
      (*Y).c6=(*X).c6;
      (*X).c6=r;

      r=(*Y).c7;
      (*Y).c7=(*X).c7;
      (*X).c7=r;

      r=(*Y).c8;
      (*Y).c8=(*X).c8;
      (*X).c8=r;

      Y+=1;
   }
}


void muladd_assign_alg(int vol,double r,su3_alg_dble *X,su3_alg_dble *Y)
{
   su3_alg_dble *Xm;

   Xm=X+vol;

   for (;X<Xm;X++)
   {
      (*Y).c1+=r*(*X).c1;
      (*Y).c2+=r*(*X).c2;
      (*Y).c3+=r*(*X).c3;
      (*Y).c4+=r*(*X).c4;
      (*Y).c5+=r*(*X).c5;
      (*Y).c6+=r*(*X).c6;
      (*Y).c7+=r*(*X).c7;
      (*Y).c8+=r*(*X).c8;

      Y+=1;
   }
}
