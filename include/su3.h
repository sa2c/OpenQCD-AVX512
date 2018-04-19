
/*******************************************************************************
*
* File su3.h
*
* Copyright (C) 2005, 2009, 2011, 2013 Martin Luescher, Filippo Palombi
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Type definitions and macros for SU(3) matrices, SU(3) vectors and Dirac
* spinors
*
*******************************************************************************/

#ifndef SU3_H
#define SU3_H

#if ((defined AVX)&&(!(defined x64)))
#define x64
#endif

#if (defined x64)
#define ALIGNED8 __attribute__ ((aligned (8)))
#define ALIGNED16 __attribute__ ((aligned (16)))
#define ALIGNED32 __attribute__ ((aligned (32)))
#else
#define ALIGNED8
#define ALIGNED16
#define ALIGNED32
#endif

typedef struct
{
   float re,im;
} complex;

typedef struct
{
   complex c1,c2,c3;
} su3_vector;

typedef struct
{
   complex c11,c12,c13,c21,c22,c23,c31,c32,c33;
} su3;

typedef struct
{
   float c1,c2,c3,c4,c5,c6,c7,c8;
} su3_alg;

typedef struct
{
   su3_vector c1,c2;
} weyl;

typedef struct
{
   su3_vector c1,c2,c3,c4;
} spinor;

typedef struct
{
   float u[36];
} pauli;

typedef struct
{
   float c1,c2,c3,c4,c5,c6,c7,c8,c9;
} u3_alg;

typedef struct
{
   double re,im;
} complex_dble;

typedef struct
{
   complex_dble c1,c2,c3;
} su3_vector_dble;

typedef struct
{
   complex_dble c11,c12,c13,c21,c22,c23,c31,c32,c33;
} su3_dble;

typedef struct
{
   double c1,c2,c3,c4,c5,c6,c7,c8;
} su3_alg_dble;

typedef struct
{
   su3_vector_dble c1,c2;
} weyl_dble;

typedef struct
{
   su3_vector_dble c1,c2,c3,c4;
} spinor_dble;

typedef struct
{
   double u[36];
} pauli_dble;

typedef struct
{
   double c1,c2,c3,c4,c5,c6,c7,c8,c9;
} u3_alg_dble;

/*******************************************************************************
*
* The following macros are the same for single and double precision types
*
* Depending on the macro, arguments are variables of type su3_vector and su3
* (or su3_vector_dble and su3_dble)
*
*******************************************************************************/

/*
* r.c1=c*s.c1 (c real)
* r.c2=c*s.c2
* r.c3=c*s.c3
*/

#define _vector_mul(r,c,s) \
   (r).c1.re=(c)*(s).c1.re; \
   (r).c1.im=(c)*(s).c1.im; \
   (r).c2.re=(c)*(s).c2.re; \
   (r).c2.im=(c)*(s).c2.im; \
   (r).c3.re=(c)*(s).c3.re; \
   (r).c3.im=(c)*(s).c3.im

/*
* r.c1=i*c*s.c1 (c real)
* r.c2=i*c*s.c2
* r.c3=i*c*s.c3
*/

#define _vector_imul(r,c,s) \
   (r).c1.re=-(c)*(s).c1.im; \
   (r).c1.im= (c)*(s).c1.re; \
   (r).c2.re=-(c)*(s).c2.im; \
   (r).c2.im= (c)*(s).c2.re; \
   (r).c3.re=-(c)*(s).c3.im; \
   (r).c3.im= (c)*(s).c3.re

/*
* r.c1=c*s.c1 (c complex)
* r.c2=c*s.c2
* r.c3=c*s.c3
*/

#define _vector_mulc(r,c,s) \
   (r).c1.re=(c).re*(s).c1.re-(c).im*(s).c1.im; \
   (r).c1.im=(c).re*(s).c1.im+(c).im*(s).c1.re; \
   (r).c2.re=(c).re*(s).c2.re-(c).im*(s).c2.im; \
   (r).c2.im=(c).re*(s).c2.im+(c).im*(s).c2.re; \
   (r).c3.re=(c).re*(s).c3.re-(c).im*(s).c3.im; \
   (r).c3.im=(c).re*(s).c3.im+(c).im*(s).c3.re

/*
* r.c1=s1.c1+s2.c1
* r.c2=s1.c2+s2.c2
* r.c3=s1.c3+s2.c3
*/

#define _vector_add(r,s1,s2) \
   (r).c1.re=(s1).c1.re+(s2).c1.re; \
   (r).c1.im=(s1).c1.im+(s2).c1.im; \
   (r).c2.re=(s1).c2.re+(s2).c2.re; \
   (r).c2.im=(s1).c2.im+(s2).c2.im; \
   (r).c3.re=(s1).c3.re+(s2).c3.re; \
   (r).c3.im=(s1).c3.im+(s2).c3.im

/*
* r.c1=s1.c1-s2.c1
* r.c2=s1.c2-s2.c2
* r.c3=s1.c3-s2.c3
*/

#define _vector_sub(r,s1,s2) \
   (r).c1.re=(s1).c1.re-(s2).c1.re; \
   (r).c1.im=(s1).c1.im-(s2).c1.im; \
   (r).c2.re=(s1).c2.re-(s2).c2.re; \
   (r).c2.im=(s1).c2.im-(s2).c2.im; \
   (r).c3.re=(s1).c3.re-(s2).c3.re; \
   (r).c3.im=(s1).c3.im-(s2).c3.im

/*
* r.c1=s1.c1+i*s2.c1
* r.c2=s1.c2+i*s2.c2
* r.c3=s1.c3+i*s2.c3
*/

#define _vector_i_add(r,s1,s2) \
   (r).c1.re=(s1).c1.re-(s2).c1.im; \
   (r).c1.im=(s1).c1.im+(s2).c1.re; \
   (r).c2.re=(s1).c2.re-(s2).c2.im; \
   (r).c2.im=(s1).c2.im+(s2).c2.re; \
   (r).c3.re=(s1).c3.re-(s2).c3.im; \
   (r).c3.im=(s1).c3.im+(s2).c3.re

/*
* r.c1=s1.c1+i*s2.c1
* r.c2=s1.c2+i*s2.c2
* r.c3=s1.c3+i*s2.c3
*/

#define _vector_i_sub(r,s1,s2) \
   (r).c1.re=(s1).c1.re+(s2).c1.im; \
   (r).c1.im=(s1).c1.im-(s2).c1.re; \
   (r).c2.re=(s1).c2.re+(s2).c2.im; \
   (r).c2.im=(s1).c2.im-(s2).c2.re; \
   (r).c3.re=(s1).c3.re+(s2).c3.im; \
   (r).c3.im=(s1).c3.im-(s2).c3.re

/*
* r.c1+=s.c1
* r.c2+=s.c2
* r.c3+=s.c3
*/

#define _vector_add_assign(r,s) \
   (r).c1.re+=(s).c1.re; \
   (r).c1.im+=(s).c1.im; \
   (r).c2.re+=(s).c2.re; \
   (r).c2.im+=(s).c2.im; \
   (r).c3.re+=(s).c3.re; \
   (r).c3.im+=(s).c3.im

/*
* r.c1-=s.c1
* r.c2-=s.c2
* r.c3-=s.c3
*/

#define _vector_sub_assign(r,s) \
   (r).c1.re-=(s).c1.re; \
   (r).c1.im-=(s).c1.im; \
   (r).c2.re-=(s).c2.re; \
   (r).c2.im-=(s).c2.im; \
   (r).c3.re-=(s).c3.re; \
   (r).c3.im-=(s).c3.im

/*
* r.c1+=i*s.c1
* r.c2+=i*s.c2
* r.c3+=i*s.c3
*/

#define _vector_i_add_assign(r,s) \
   (r).c1.re-=(s).c1.im; \
   (r).c1.im+=(s).c1.re; \
   (r).c2.re-=(s).c2.im; \
   (r).c2.im+=(s).c2.re; \
   (r).c3.re-=(s).c3.im; \
   (r).c3.im+=(s).c3.re

/*
* r.c1-=i*s.c1
* r.c2-=i*s.c2
* r.c3-=i*s.c3
*/

#define _vector_i_sub_assign(r,s) \
   (r).c1.re+=(s).c1.im; \
   (r).c1.im-=(s).c1.re; \
   (r).c2.re+=(s).c2.im; \
   (r).c2.im-=(s).c2.re; \
   (r).c3.re+=(s).c3.im; \
   (r).c3.im-=(s).c3.re

/*
* Real part of the scalar product (r,s)
*/

#define _vector_prod_re(r,s) \
   (r).c1.re*(s).c1.re+(r).c1.im*(s).c1.im+ \
   (r).c2.re*(s).c2.re+(r).c2.im*(s).c2.im+ \
   (r).c3.re*(s).c3.re+(r).c3.im*(s).c3.im

/*
* Imaginary part of the scalar product (r,s)
*/

#define _vector_prod_im(r,s) \
   (r).c1.re*(s).c1.im-(r).c1.im*(s).c1.re+ \
   (r).c2.re*(s).c2.im-(r).c2.im*(s).c2.re+ \
   (r).c3.re*(s).c3.im-(r).c3.im*(s).c3.re

/*
* r.c1+=c*s.c1 (c real)
* r.c2+=c*s.c2
* r.c3+=c*s.c3
*/

#define _vector_mulr_assign(r,c,s) \
   (r).c1.re+=(c)*(s).c1.re; \
   (r).c1.im+=(c)*(s).c1.im; \
   (r).c2.re+=(c)*(s).c2.re; \
   (r).c2.im+=(c)*(s).c2.im; \
   (r).c3.re+=(c)*(s).c3.re; \
   (r).c3.im+=(c)*(s).c3.im

/*
* r.c1+=i*c*s.c1 (c real)
* r.c2+=i*c*s.c2
* r.c3+=i*c*s.c3
*/

#define _vector_mulir_assign(r,c,s) \
   (r).c1.re-=(c)*(s).c1.im; \
   (r).c1.im+=(c)*(s).c1.re; \
   (r).c2.re-=(c)*(s).c2.im; \
   (r).c2.im+=(c)*(s).c2.re; \
   (r).c3.re-=(c)*(s).c3.im; \
   (r).c3.im+=(c)*(s).c3.re

/*
* r.c1+=z*s.c1 (z of type complex)
* r.c2+=z*s.c2
* r.c3+=z*s.c3
*/

#define _vector_mulc_assign(r,z,s) \
   (r).c1.re+=((z).re*(s).c1.re-(z).im*(s).c1.im); \
   (r).c1.im+=((z).re*(s).c1.im+(z).im*(s).c1.re); \
   (r).c2.re+=((z).re*(s).c2.re-(z).im*(s).c2.im); \
   (r).c2.im+=((z).re*(s).c2.im+(z).im*(s).c2.re); \
   (r).c3.re+=((z).re*(s).c3.re-(z).im*(s).c3.im); \
   (r).c3.im+=((z).re*(s).c3.im+(z).im*(s).c3.re)

/*
* r.c1-=z*s.c1 (z of type complex)
* r.c2-=z*s.c2
* r.c3-=z*s.c3
*/

#define _vector_project(r,z,s) \
   (r).c1.re-=((z).re*(s).c1.re-(z).im*(s).c1.im); \
   (r).c1.im-=((z).re*(s).c1.im+(z).im*(s).c1.re); \
   (r).c2.re-=((z).re*(s).c2.re-(z).im*(s).c2.im); \
   (r).c2.im-=((z).re*(s).c2.im+(z).im*(s).c2.re); \
   (r).c3.re-=((z).re*(s).c3.re-(z).im*(s).c3.im); \
   (r).c3.im-=((z).re*(s).c3.im+(z).im*(s).c3.re)

/*
* r.c1=c*r.c1+s.c1 (c real)
* r.c2=c*r.c2+s.c2
* r.c3=c*r.c3+s.c3
*/

#define _vector_mulr_add(r,c,s) \
   (r).c1.re=(c)*(r).c1.re+(s).c1.re; \
   (r).c1.im=(c)*(r).c1.im+(s).c1.im; \
   (r).c2.re=(c)*(r).c2.re+(s).c2.re; \
   (r).c2.im=(c)*(r).c2.im+(s).c2.im; \
   (r).c3.re=(c)*(r).c3.re+(s).c3.re; \
   (r).c3.im=(c)*(r).c3.im+(s).c3.im

/*
* r.c1=cr*r.c1+cs*s.c1 (cr,cs real)
* r.c2=cr*r.c2+cs*s.c2
* r.c3=cr*r.c3+cs*s.c3
*/

#define _vector_combine(r,s,cr,cs) \
   (r).c1.re=(cr)*(r).c1.re+(cs)*(s).c1.re; \
   (r).c1.im=(cr)*(r).c1.im+(cs)*(s).c1.im; \
   (r).c2.re=(cr)*(r).c2.re+(cs)*(s).c2.re; \
   (r).c2.im=(cr)*(r).c2.im+(cs)*(s).c2.im; \
   (r).c3.re=(cr)*(r).c3.re+(cs)*(s).c3.re; \
   (r).c3.im=(cr)*(r).c3.im+(cs)*(s).c3.im

/*
* v.c1=(w.c2*z.c3-w.c3*z.c2)^*
* v.c2=(w.c3*z.c1-w.c1*z.c3)^*
* v.c3=(w.c1*z.c2-w.c2*z.c1)^*
*/

#define _vector_cross_prod(v,w,z) \
   (v).c1.re= (w).c2.re*(z).c3.re-(w).c2.im*(z).c3.im  \
             -(w).c3.re*(z).c2.re+(w).c3.im*(z).c2.im; \
   (v).c1.im= (w).c3.re*(z).c2.im+(w).c3.im*(z).c2.re  \
             -(w).c2.re*(z).c3.im-(w).c2.im*(z).c3.re; \
   (v).c2.re= (w).c3.re*(z).c1.re-(w).c3.im*(z).c1.im  \
             -(w).c1.re*(z).c3.re+(w).c1.im*(z).c3.im; \
   (v).c2.im= (w).c1.re*(z).c3.im+(w).c1.im*(z).c3.re  \
             -(w).c3.re*(z).c1.im-(w).c3.im*(z).c1.re; \
   (v).c3.re= (w).c1.re*(z).c2.re-(w).c1.im*(z).c2.im  \
             -(w).c2.re*(z).c1.re+(w).c2.im*(z).c1.im; \
   (v).c3.im= (w).c2.re*(z).c1.im+(w).c2.im*(z).c1.re  \
             -(w).c1.re*(z).c2.im-(w).c1.im*(z).c2.re

/*
* SU(3) matrix u times SU(3) vector s
*
* r.c1=(u*s).c1
* r.c2=(u*s).c2
* r.c3=(u*s).c3
*/

#define _su3_multiply(r,u,s) \
   (r).c1.re= (u).c11.re*(s).c1.re-(u).c11.im*(s).c1.im  \
             +(u).c12.re*(s).c2.re-(u).c12.im*(s).c2.im  \
             +(u).c13.re*(s).c3.re-(u).c13.im*(s).c3.im; \
   (r).c1.im= (u).c11.re*(s).c1.im+(u).c11.im*(s).c1.re  \
             +(u).c12.re*(s).c2.im+(u).c12.im*(s).c2.re  \
             +(u).c13.re*(s).c3.im+(u).c13.im*(s).c3.re; \
   (r).c2.re= (u).c21.re*(s).c1.re-(u).c21.im*(s).c1.im  \
             +(u).c22.re*(s).c2.re-(u).c22.im*(s).c2.im  \
             +(u).c23.re*(s).c3.re-(u).c23.im*(s).c3.im; \
   (r).c2.im= (u).c21.re*(s).c1.im+(u).c21.im*(s).c1.re  \
             +(u).c22.re*(s).c2.im+(u).c22.im*(s).c2.re  \
             +(u).c23.re*(s).c3.im+(u).c23.im*(s).c3.re; \
   (r).c3.re= (u).c31.re*(s).c1.re-(u).c31.im*(s).c1.im  \
             +(u).c32.re*(s).c2.re-(u).c32.im*(s).c2.im  \
             +(u).c33.re*(s).c3.re-(u).c33.im*(s).c3.im; \
   (r).c3.im= (u).c31.re*(s).c1.im+(u).c31.im*(s).c1.re  \
             +(u).c32.re*(s).c2.im+(u).c32.im*(s).c2.re  \
             +(u).c33.re*(s).c3.im+(u).c33.im*(s).c3.re

/*
* SU(3) matrix u^dagger times SU(3) vector s
*
* r.c1=(u^dagger*s).c1
* r.c2=(u^dagger*s).c2
* r.c3=(u^dagger*s).c3
*/

#define _su3_inverse_multiply(r,u,s) \
   (r).c1.re= (u).c11.re*(s).c1.re+(u).c11.im*(s).c1.im  \
             +(u).c21.re*(s).c2.re+(u).c21.im*(s).c2.im  \
             +(u).c31.re*(s).c3.re+(u).c31.im*(s).c3.im; \
   (r).c1.im= (u).c11.re*(s).c1.im-(u).c11.im*(s).c1.re  \
             +(u).c21.re*(s).c2.im-(u).c21.im*(s).c2.re  \
             +(u).c31.re*(s).c3.im-(u).c31.im*(s).c3.re; \
   (r).c2.re= (u).c12.re*(s).c1.re+(u).c12.im*(s).c1.im  \
             +(u).c22.re*(s).c2.re+(u).c22.im*(s).c2.im  \
             +(u).c32.re*(s).c3.re+(u).c32.im*(s).c3.im; \
   (r).c2.im= (u).c12.re*(s).c1.im-(u).c12.im*(s).c1.re  \
             +(u).c22.re*(s).c2.im-(u).c22.im*(s).c2.re  \
             +(u).c32.re*(s).c3.im-(u).c32.im*(s).c3.re; \
   (r).c3.re= (u).c13.re*(s).c1.re+(u).c13.im*(s).c1.im  \
             +(u).c23.re*(s).c2.re+(u).c23.im*(s).c2.im  \
             +(u).c33.re*(s).c3.re+(u).c33.im*(s).c3.im; \
   (r).c3.im= (u).c13.re*(s).c1.im-(u).c13.im*(s).c1.re  \
             +(u).c23.re*(s).c2.im-(u).c23.im*(s).c2.re  \
             +(u).c33.re*(s).c3.im-(u).c33.im*(s).c3.re

/*******************************************************************************
*
* Macros for SU(3) matrices
*
* Arguments are variables of type su3
*
*******************************************************************************/

/*
* u=v^dagger
*/

#define _su3_dagger(u,v) \
   (u).c11.re= (v).c11.re; \
   (u).c11.im=-(v).c11.im; \
   (u).c12.re= (v).c21.re; \
   (u).c12.im=-(v).c21.im; \
   (u).c13.re= (v).c31.re; \
   (u).c13.im=-(v).c31.im; \
   (u).c21.re= (v).c12.re; \
   (u).c21.im=-(v).c12.im; \
   (u).c22.re= (v).c22.re; \
   (u).c22.im=-(v).c22.im; \
   (u).c23.re= (v).c32.re; \
   (u).c23.im=-(v).c32.im; \
   (u).c31.re= (v).c13.re; \
   (u).c31.im=-(v).c13.im; \
   (u).c32.re= (v).c23.re; \
   (u).c32.im=-(v).c23.im; \
   (u).c33.re= (v).c33.re; \
   (u).c33.im=-(v).c33.im

/*
* u=v*w
*/

#define _su3_times_su3(u,v,w) \
   (u).c11.re= (v).c11.re*(w).c11.re-(v).c11.im*(w).c11.im  \
              +(v).c12.re*(w).c21.re-(v).c12.im*(w).c21.im  \
              +(v).c13.re*(w).c31.re-(v).c13.im*(w).c31.im; \
   (u).c11.im= (v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
              +(v).c12.re*(w).c21.im+(v).c12.im*(w).c21.re  \
              +(v).c13.re*(w).c31.im+(v).c13.im*(w).c31.re; \
   (u).c12.re= (v).c11.re*(w).c12.re-(v).c11.im*(w).c12.im  \
              +(v).c12.re*(w).c22.re-(v).c12.im*(w).c22.im  \
              +(v).c13.re*(w).c32.re-(v).c13.im*(w).c32.im; \
   (u).c12.im= (v).c11.re*(w).c12.im+(v).c11.im*(w).c12.re  \
              +(v).c12.re*(w).c22.im+(v).c12.im*(w).c22.re  \
              +(v).c13.re*(w).c32.im+(v).c13.im*(w).c32.re; \
   (u).c13.re= (v).c11.re*(w).c13.re-(v).c11.im*(w).c13.im  \
              +(v).c12.re*(w).c23.re-(v).c12.im*(w).c23.im  \
              +(v).c13.re*(w).c33.re-(v).c13.im*(w).c33.im; \
   (u).c13.im= (v).c11.re*(w).c13.im+(v).c11.im*(w).c13.re  \
              +(v).c12.re*(w).c23.im+(v).c12.im*(w).c23.re  \
              +(v).c13.re*(w).c33.im+(v).c13.im*(w).c33.re; \
   (u).c21.re= (v).c21.re*(w).c11.re-(v).c21.im*(w).c11.im  \
              +(v).c22.re*(w).c21.re-(v).c22.im*(w).c21.im  \
              +(v).c23.re*(w).c31.re-(v).c23.im*(w).c31.im; \
   (u).c21.im= (v).c21.re*(w).c11.im+(v).c21.im*(w).c11.re  \
              +(v).c22.re*(w).c21.im+(v).c22.im*(w).c21.re  \
              +(v).c23.re*(w).c31.im+(v).c23.im*(w).c31.re; \
   (u).c22.re= (v).c21.re*(w).c12.re-(v).c21.im*(w).c12.im  \
              +(v).c22.re*(w).c22.re-(v).c22.im*(w).c22.im  \
              +(v).c23.re*(w).c32.re-(v).c23.im*(w).c32.im; \
   (u).c22.im= (v).c21.re*(w).c12.im+(v).c21.im*(w).c12.re  \
              +(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re  \
              +(v).c23.re*(w).c32.im+(v).c23.im*(w).c32.re; \
   (u).c23.re= (v).c21.re*(w).c13.re-(v).c21.im*(w).c13.im  \
              +(v).c22.re*(w).c23.re-(v).c22.im*(w).c23.im  \
              +(v).c23.re*(w).c33.re-(v).c23.im*(w).c33.im; \
   (u).c23.im= (v).c21.re*(w).c13.im+(v).c21.im*(w).c13.re  \
              +(v).c22.re*(w).c23.im+(v).c22.im*(w).c23.re  \
              +(v).c23.re*(w).c33.im+(v).c23.im*(w).c33.re; \
   (u).c31.re= (v).c31.re*(w).c11.re-(v).c31.im*(w).c11.im  \
              +(v).c32.re*(w).c21.re-(v).c32.im*(w).c21.im  \
              +(v).c33.re*(w).c31.re-(v).c33.im*(w).c31.im; \
   (u).c31.im= (v).c31.re*(w).c11.im+(v).c31.im*(w).c11.re  \
              +(v).c32.re*(w).c21.im+(v).c32.im*(w).c21.re  \
              +(v).c33.re*(w).c31.im+(v).c33.im*(w).c31.re; \
   (u).c32.re= (v).c31.re*(w).c12.re-(v).c31.im*(w).c12.im  \
              +(v).c32.re*(w).c22.re-(v).c32.im*(w).c22.im  \
              +(v).c33.re*(w).c32.re-(v).c33.im*(w).c32.im; \
   (u).c32.im= (v).c31.re*(w).c12.im+(v).c31.im*(w).c12.re  \
              +(v).c32.re*(w).c22.im+(v).c32.im*(w).c22.re  \
              +(v).c33.re*(w).c32.im+(v).c33.im*(w).c32.re; \
   (u).c33.re= (v).c31.re*(w).c13.re-(v).c31.im*(w).c13.im  \
              +(v).c32.re*(w).c23.re-(v).c32.im*(w).c23.im  \
              +(v).c33.re*(w).c33.re-(v).c33.im*(w).c33.im; \
   (u).c33.im= (v).c31.re*(w).c13.im+(v).c31.im*(w).c13.re  \
              +(v).c32.re*(w).c23.im+(v).c32.im*(w).c23.re  \
              +(v).c33.re*(w).c33.im+(v).c33.im*(w).c33.re

/*******************************************************************************
*
* Macros for variables of type su3_alg
*
*******************************************************************************/

/*
* r+=s
*/

#define _su3_alg_add_assign(r,s) \
   (r).c1+=(s).c1; \
   (r).c2+=(s).c2; \
   (r).c3+=(s).c3; \
   (r).c4+=(s).c4; \
   (r).c5+=(s).c5; \
   (r).c6+=(s).c6; \
   (r).c7+=(s).c7; \
   (r).c8+=(s).c8

/*
* r-=s
*/

#define _su3_alg_sub_assign(r,s) \
   (r).c1-=(s).c1; \
   (r).c2-=(s).c2; \
   (r).c3-=(s).c3; \
   (r).c4-=(s).c4; \
   (r).c5-=(s).c5; \
   (r).c6-=(s).c6; \
   (r).c7-=(s).c7; \
   (r).c8-=(s).c8

/*
* s*=c, c real
*/

#define _su3_alg_mul_assign(s,c) \
   (s).c1*=(c); \
   (s).c2*=(c); \
   (s).c3*=(c); \
   (s).c4*=(c); \
   (s).c5*=(c); \
   (s).c6*=(c); \
   (s).c7*=(c); \
   (s).c8*=(c)

/*
* r+=c*s, c real
*/

#define _su3_alg_mul_add_assign(r,c,s) \
   (r).c1+=(c)*(s).c1; \
   (r).c2+=(c)*(s).c2; \
   (r).c3+=(c)*(s).c3; \
   (r).c4+=(c)*(s).c4; \
   (r).c5+=(c)*(s).c5; \
   (r).c6+=(c)*(s).c6; \
   (r).c7+=(c)*(s).c7; \
   (r).c8+=(c)*(s).c8

/*
* r-=c*s, c real
*/

#define _su3_alg_mul_sub_assign(r,c,s) \
   (r).c1-=(c)*(s).c1; \
   (r).c2-=(c)*(s).c2; \
   (r).c3-=(c)*(s).c3; \
   (r).c4-=(c)*(s).c4; \
   (r).c5-=(c)*(s).c5; \
   (r).c6-=(c)*(s).c6; \
   (r).c7-=(c)*(s).c7; \
   (r).c8-=(c)*(s).c8

/*******************************************************************************
*
* Macros for variables of type u3_alg
*
*******************************************************************************/

/*
* r=c*(u+v)
*/

#define _u3_alg_mul_add(r,c,u,v) \
   (r).c1=(c)*((u).c1+(v).c1); \
   (r).c2=(c)*((u).c2+(v).c2); \
   (r).c3=(c)*((u).c3+(v).c3); \
   (r).c4=(c)*((u).c4+(v).c4); \
   (r).c5=(c)*((u).c5+(v).c5); \
   (r).c6=(c)*((u).c6+(v).c6); \
   (r).c7=(c)*((u).c7+(v).c7); \
   (r).c8=(c)*((u).c8+(v).c8); \
   (r).c9=(c)*((u).c9+(v).c9)

/*
* r=c*(u-v)
*/

#define _u3_alg_mul_sub(r,c,u,v) \
   (r).c1=(c)*((u).c1-(v).c1); \
   (r).c2=(c)*((u).c2-(v).c2); \
   (r).c3=(c)*((u).c3-(v).c3); \
   (r).c4=(c)*((u).c4-(v).c4); \
   (r).c5=(c)*((u).c5-(v).c5); \
   (r).c6=(c)*((u).c6-(v).c6); \
   (r).c7=(c)*((u).c7-(v).c7); \
   (r).c8=(c)*((u).c8-(v).c8); \
   (r).c9=(c)*((u).c9-(v).c9)

#endif
