
/*******************************************************************************
*
* File random_su3.c
*
* Copyright (C) 2004, 2009 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generation of uniformly distributed single- and double-precision
* SU(3) matrices
*
* The externally accessible function is
*
*   void random_su3(su3 *u)
*     Generates a random single-precision SU(3) matrix and assigns it to *u
*
*   void random_su3_dble(su3_dble *u)
*     Generates a random double-precision SU(3) matrix and assigns it to *u
*
* Notes:
*
* The random matrices are uniformly distributed over SU(3) to a precision
* given by the number of significant bits of the random numbers returned by
* ranlxs and ranlxd respectively. Rougly speaking one can expect the matrices
* to be uniformly random up to systematic deviations from 1 at the level of 
* 10^(-7) and 10^(-14) in the single- and double-precision programs.
*
*******************************************************************************/

#define RANDOM_SU3_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "su3.h"
#include "random.h"
#include "su3fcts.h"

static float rs[6];
static double rd[6];
static su3_vector vs1,vs2,vs3;
static su3_vector_dble vd1,vd2,vd3;


static void random_su3_vector(su3_vector *v)
{
   float norm,fact;

   norm=0.0f;

   while (norm<=0.1f)
   {
      gauss(rs,6);
      norm=rs[0]*rs[0]+rs[1]*rs[1]+rs[2]*rs[2]+
           rs[3]*rs[3]+rs[4]*rs[4]+rs[5]*rs[5];
   }

   fact=1.0f/(float)sqrt((double)(norm));

   (*v).c1.re=fact*rs[0];
   (*v).c1.im=fact*rs[1];   
   (*v).c2.re=fact*rs[2];
   (*v).c2.im=fact*rs[3];
   (*v).c3.re=fact*rs[4];
   (*v).c3.im=fact*rs[5];
}


void random_su3(su3 *u)
{
   float norm,fact;

   random_su3_vector(&vs1);
   norm=0.0f;

   while (norm<=0.1f)
   {
      random_su3_vector(&vs2);
      _vector_cross_prod(vs3,vs1,vs2);
      norm=_vector_prod_re(vs3,vs3);
   }        

   fact=1.0f/(float)sqrt((double)(norm));

   vs3.c1.re*=fact;
   vs3.c1.im*=fact;
   vs3.c2.re*=fact;
   vs3.c2.im*=fact;
   vs3.c3.re*=fact;
   vs3.c3.im*=fact;

   _vector_cross_prod(vs2,vs3,vs1);

   (*u).c11.re=vs1.c1.re;
   (*u).c11.im=vs1.c1.im;
   (*u).c12.re=vs1.c2.re;
   (*u).c12.im=vs1.c2.im;
   (*u).c13.re=vs1.c3.re;
   (*u).c13.im=vs1.c3.im;   

   (*u).c21.re=vs2.c1.re;
   (*u).c21.im=vs2.c1.im;
   (*u).c22.re=vs2.c2.re;
   (*u).c22.im=vs2.c2.im;
   (*u).c23.re=vs2.c3.re;
   (*u).c23.im=vs2.c3.im;

   (*u).c31.re=vs3.c1.re;
   (*u).c31.im=vs3.c1.im;
   (*u).c32.re=vs3.c2.re;
   (*u).c32.im=vs3.c2.im;
   (*u).c33.re=vs3.c3.re;
   (*u).c33.im=vs3.c3.im;   
}


static void random_su3_vector_dble(su3_vector_dble *v)
{
   double norm,fact;

   norm=0.0;

   while (norm<=0.1)
   {
      gauss_dble(rd,6);
      norm=rd[0]*rd[0]+rd[1]*rd[1]+rd[2]*rd[2]+
           rd[3]*rd[3]+rd[4]*rd[4]+rd[5]*rd[5];
   }

   fact=1.0/sqrt(norm);

   (*v).c1.re=fact*rd[0];
   (*v).c1.im=fact*rd[1];   
   (*v).c2.re=fact*rd[2];
   (*v).c2.im=fact*rd[3];
   (*v).c3.re=fact*rd[4];
   (*v).c3.im=fact*rd[5];
}


void random_su3_dble(su3_dble *u)
{
   double norm,fact;

   random_su3_vector_dble(&vd1);
   norm=0.0;

   while (norm<=0.1)
   {
      random_su3_vector_dble(&vd2);
      _vector_cross_prod(vd3,vd1,vd2);      
      norm=_vector_prod_re(vd3,vd3);
   }        

   fact=1.0/sqrt(norm);

   vd3.c1.re*=fact;
   vd3.c1.im*=fact;
   vd3.c2.re*=fact;
   vd3.c2.im*=fact;
   vd3.c3.re*=fact;
   vd3.c3.im*=fact;

   _vector_cross_prod(vd2,vd3,vd1);

   (*u).c11.re=vd1.c1.re;
   (*u).c11.im=vd1.c1.im;
   (*u).c12.re=vd1.c2.re;
   (*u).c12.im=vd1.c2.im;
   (*u).c13.re=vd1.c3.re;
   (*u).c13.im=vd1.c3.im;   

   (*u).c21.re=vd2.c1.re;
   (*u).c21.im=vd2.c1.im;
   (*u).c22.re=vd2.c2.re;
   (*u).c22.im=vd2.c2.im;
   (*u).c23.re=vd2.c3.re;
   (*u).c23.im=vd2.c3.im;

   (*u).c31.re=vd3.c1.re;
   (*u).c31.im=vd3.c1.im;
   (*u).c32.re=vd3.c2.re;
   (*u).c32.im=vd3.c2.im;
   (*u).c33.re=vd3.c3.re;
   (*u).c33.im=vd3.c3.im;   
}
