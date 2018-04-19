
/*******************************************************************************
*
* File su3ren.c
*
* Copyright (C) 2005, 2009, 2010, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Renormalization of SU(3) matrices
*
* The externally accessible function are
*
*   void project_to_su3(su3 *u)
*     Projects an approximate single-precision SU(3) matrix back to SU(3).
*     No action is performed if the matrix is degenerate
*
*   void project_to_su3_dble(su3_dble *u)
*     Projects an approximate double-precision SU(3) matrix back to SU(3).
*     No action is performed if the matrix is degenerate
*
* Notes:
*
* The programs in this module do not perform any communications and can be 
* called locally. A matrix is considered to be degenerate if the first
* column vector or the cross-product of the first and the second vector is
* exactly equal to zero.
*
*******************************************************************************/

#define SU3REN_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "su3fcts.h"


static int normalize(su3_vector *v)
{
   float r;

   r=_vector_prod_re((*v),(*v));
   r=(float)sqrt((double)(r));
   
   if (r==0.0f)
      return 1;
   else
   {
      r=1.0f/r;
      _vector_mul((*v),r,(*v));
      return 0;
   }
}


static int normalize_dble(su3_vector_dble *v)
{
   double r;

   r=_vector_prod_re((*v),(*v));
   r=sqrt(r);

   if (r==0.0)
      return 1;
   else
   {
      r=1.0/r;
      _vector_mul((*v),r,(*v));
      return 0;
   }
}


void project_to_su3(su3 *u)
{
   int it;
   su3_vector v1,v2,v3;

   v1.c1.re=(*u).c11.re;
   v1.c1.im=(*u).c11.im;   
   v1.c2.re=(*u).c12.re;
   v1.c2.im=(*u).c12.im;
   v1.c3.re=(*u).c13.re;
   v1.c3.im=(*u).c13.im;

   v2.c1.re=(*u).c21.re;
   v2.c1.im=(*u).c21.im;   
   v2.c2.re=(*u).c22.re;
   v2.c2.im=(*u).c22.im;
   v2.c3.re=(*u).c23.re;
   v2.c3.im=(*u).c23.im;

   v3.c1.re=(*u).c31.re;
   v3.c1.im=(*u).c31.im;   
   v3.c2.re=(*u).c32.re;
   v3.c2.im=(*u).c32.im;
   v3.c3.re=(*u).c33.re;
   v3.c3.im=(*u).c33.im;    
   
   it=normalize(&v1);
   _vector_cross_prod(v3,v1,v2);
   it|=normalize(&v3);
   _vector_cross_prod(v2,v3,v1);   

   if (it==0)
   {
      (*u).c11.re=v1.c1.re;
      (*u).c11.im=v1.c1.im;
      (*u).c12.re=v1.c2.re;
      (*u).c12.im=v1.c2.im;
      (*u).c13.re=v1.c3.re;
      (*u).c13.im=v1.c3.im;

      (*u).c21.re=v2.c1.re;
      (*u).c21.im=v2.c1.im;
      (*u).c22.re=v2.c2.re;
      (*u).c22.im=v2.c2.im;
      (*u).c23.re=v2.c3.re;
      (*u).c23.im=v2.c3.im;

      (*u).c31.re=v3.c1.re;
      (*u).c31.im=v3.c1.im;
      (*u).c32.re=v3.c2.re;
      (*u).c32.im=v3.c2.im;
      (*u).c33.re=v3.c3.re;
      (*u).c33.im=v3.c3.im;
   }
}


void project_to_su3_dble(su3_dble *u)
{
   int it;
   su3_vector_dble v1,v2,v3;

   v1.c1.re=(*u).c11.re;
   v1.c1.im=(*u).c11.im;   
   v1.c2.re=(*u).c12.re;
   v1.c2.im=(*u).c12.im;
   v1.c3.re=(*u).c13.re;
   v1.c3.im=(*u).c13.im;

   v2.c1.re=(*u).c21.re;
   v2.c1.im=(*u).c21.im;   
   v2.c2.re=(*u).c22.re;
   v2.c2.im=(*u).c22.im;
   v2.c3.re=(*u).c23.re;
   v2.c3.im=(*u).c23.im;

   v3.c1.re=(*u).c31.re;
   v3.c1.im=(*u).c31.im;   
   v3.c2.re=(*u).c32.re;
   v3.c2.im=(*u).c32.im;
   v3.c3.re=(*u).c33.re;
   v3.c3.im=(*u).c33.im;    
   
   it=normalize_dble(&v1);
   _vector_cross_prod(v3,v1,v2);
   it|=normalize_dble(&v3);
   _vector_cross_prod(v2,v3,v1);   

   if (it==0)
   {
      (*u).c11.re=v1.c1.re;
      (*u).c11.im=v1.c1.im;
      (*u).c12.re=v1.c2.re;
      (*u).c12.im=v1.c2.im;
      (*u).c13.re=v1.c3.re;
      (*u).c13.im=v1.c3.im;

      (*u).c21.re=v2.c1.re;
      (*u).c21.im=v2.c1.im;
      (*u).c22.re=v2.c2.re;
      (*u).c22.im=v2.c2.im;
      (*u).c23.re=v2.c3.re;
      (*u).c23.im=v2.c3.im;

      (*u).c31.re=v3.c1.re;
      (*u).c31.im=v3.c1.im;
      (*u).c32.re=v3.c2.re;
      (*u).c32.im=v3.c2.im;
      (*u).c33.re=v3.c3.re;
      (*u).c33.im=v3.c3.im;
   }
}
