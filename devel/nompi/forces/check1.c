
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of det2xt and prod2xt
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "random.h"
#include "utils.h"
#include "su3fcts.h"
#include "sw_term.h"
#include "forces.h"

typedef union
{
   spinor_dble s;
   weyl_dble w[2];
   double r[24];
} spin_t;

typedef union
{
   su3_vector_dble v;
   double r[6];
} vec_t;

static int pln[6][2]={{0,1},{0,2},{0,3},{2,3},{3,1},{1,2}};
static const su3_vector_dble vd0={{0.0,0.0},{0.0,0.0},{0.0,0.0}};
static const spinor_dble sd0={{{0.0,0.0},{0.0,0.0},{0.0,0.0}},
                              {{0.0,0.0},{0.0,0.0},{0.0,0.0}},
                              {{0.0,0.0},{0.0,0.0},{0.0,0.0}},
                              {{0.0,0.0},{0.0,0.0},{0.0,0.0}}};
static su3_dble Q ALIGNED16;
static spin_t s1 ALIGNED16;
static spin_t s2 ALIGNED16;
static spin_t s3 ALIGNED16;
static spin_t s4 ALIGNED16;
static pauli_dble m[2] ALIGNED16;


static su3_vector_dble mul_cplx(complex_dble z,su3_vector_dble s)
{
   su3_vector_dble r;

   r.c1.re=z.re*s.c1.re-z.im*s.c1.im;
   r.c1.im=z.im*s.c1.re+z.re*s.c1.im;
   r.c2.re=z.re*s.c2.re-z.im*s.c2.im;
   r.c2.im=z.im*s.c2.re+z.re*s.c2.im;
   r.c3.re=z.re*s.c3.re-z.im*s.c3.im;
   r.c3.im=z.im*s.c3.re+z.re*s.c3.im;

   return r;
}


static spinor_dble mul_gamma(int mu,spinor_dble s)
{
   spinor_dble r;
   complex_dble i,m_i,m_1;

   i.re=0.0;
   i.im=1.0;

   m_i.re=0.0;
   m_i.im=-1.0;

   m_1.re=-1.0;
   m_1.im=0.0;

   if (mu==0)
   {
      r.c1=mul_cplx(m_1,s.c3);
      r.c2=mul_cplx(m_1,s.c4);
      r.c3=mul_cplx(m_1,s.c1);
      r.c4=mul_cplx(m_1,s.c2);
   }
   else if (mu==1)
   {
      r.c1=mul_cplx(m_i,s.c4);
      r.c2=mul_cplx(m_i,s.c3);
      r.c3=mul_cplx(i,s.c2);
      r.c4=mul_cplx(i,s.c1);
   }
   else if (mu==2)
   {
      r.c1=mul_cplx(m_1,s.c4);
      r.c2=s.c3;
      r.c3=s.c2;
      r.c4=mul_cplx(m_1,s.c1);
   }
   else if (mu==3)
   {
      r.c1=mul_cplx(m_i,s.c3);
      r.c2=mul_cplx(i,s.c4);
      r.c3=mul_cplx(i,s.c1);
      r.c4=mul_cplx(m_i,s.c2);
   }
   else
   {
      r.c1=s.c1;
      r.c2=s.c2;
      r.c3=mul_cplx(m_1,s.c3);
      r.c4=mul_cplx(m_1,s.c4);
   }

   return r;
}


static spinor_dble mul_sigma(int mu,int nu,spinor_dble s)
{
   complex_dble z;
   spinor_dble r1,r2;

   r1=mul_gamma(nu,s);
   r1=mul_gamma(mu,r1);

   r2=mul_gamma(mu,s);
   r2=mul_gamma(nu,r2);

   _vector_sub_assign(r1.c1,r2.c1);
   _vector_sub_assign(r1.c2,r2.c2);
   _vector_sub_assign(r1.c3,r2.c3);
   _vector_sub_assign(r1.c4,r2.c4);

   z.re=0.0;
   z.im=0.5;
   _vector_mulc(r2.c1,z,r1.c1);
   _vector_mulc(r2.c2,z,r1.c2);
   _vector_mulc(r2.c3,z,r1.c3);
   _vector_mulc(r2.c4,z,r1.c4);

   return r2;
}


static spinor_dble mul_Fhat(su3_dble Q,spinor_dble s)
{
   su3_dble F;
   spinor_dble r;

   F.c11.re=0.0;
   F.c11.im=0.25*Q.c11.im;
   F.c22.re=0.0;
   F.c22.im=0.25*Q.c22.im;
   F.c33.re=0.0;
   F.c33.im=0.25*Q.c33.im;

   F.c12.re=0.125*(Q.c12.re-Q.c21.re);
   F.c12.im=0.125*(Q.c12.im+Q.c21.im);
   F.c21.re=-F.c12.re;
   F.c21.im=F.c12.im;

   F.c13.re=0.125*(Q.c13.re-Q.c31.re);
   F.c13.im=0.125*(Q.c13.im+Q.c31.im);
   F.c31.re=-F.c13.re;
   F.c31.im=F.c13.im;

   F.c23.re=0.125*(Q.c23.re-Q.c32.re);
   F.c23.im=0.125*(Q.c23.im+Q.c32.im);
   F.c32.re=-F.c23.re;
   F.c32.im=F.c23.im;

   _su3_multiply(r.c1,F,s.c1);
   _su3_multiply(r.c2,F,s.c2);
   _su3_multiply(r.c3,F,s.c3);
   _su3_multiply(r.c4,F,s.c4);

   return r;
}


static su3_vector_dble mul_X(u3_alg_dble X,su3_vector_dble s)
{
   su3_dble M;
   su3_vector_dble r;

   M.c11.re=0.0;
   M.c11.im=X.c1;
   M.c22.re=0.0;
   M.c22.im=X.c2;
   M.c33.re=0.0;
   M.c33.im=X.c3;

   M.c12.re=X.c4;
   M.c12.im=X.c5;
   M.c21.re=-X.c4;
   M.c21.im=X.c5;

   M.c13.re=X.c6;
   M.c13.im=X.c7;
   M.c31.re=-X.c6;
   M.c31.im=X.c7;

   M.c23.re=X.c8;
   M.c23.im=X.c9;
   M.c32.re=-X.c8;
   M.c32.im=X.c9;

   _su3_multiply(r,M,s);

   return r;
}


int main(void)
{
   int n,mu,nu,i;
   complex_dble z;
   vec_t v1,v2,v3;
   u3_alg_dble X[6];

   printf("\n");
   printf("Check of det2xt and prod2xt\n");
   printf("---------------------------\n\n");

   rlxd_init(1,23456);

   ranlxd(v1.r,6);
   ranlxd(v2.r,6);
   ranlxd(v3.r,6);

   ranlxd(s1.r,24);
   ranlxd(s2.r,24);
   ranlxd(s3.r,24);
   ranlxd(s4.r,24);

   ranlxd(m[0].u,36);
   ranlxd(m[1].u,36);

   det2xt(m,X);

   printf("det2xt:\n");

   for (n=0;n<6;n++)
   {
      mu=pln[n][0];
      nu=pln[n][1];

      random_su3_dble(&Q);
      z.im=0.0;

      for (i=0;i<12;i++)
      {
         s1.s=sd0;
         s1.r[2*i]=1.0;

         mul_pauli_dble(0.0,m,s1.w,s2.w);
         mul_pauli_dble(0.0,m+1,s1.w+1,s2.w+1);
         s1.s=mul_sigma(mu,nu,s2.s);
         s2.s=mul_Fhat(Q,s1.s);

         z.im-=s2.r[2*i+1];
      }

      z.re=0.0;

      for (i=0;i<3;i++)
      {
         v1.v=vd0;
         v1.r[2*i]=1.0;

         v2.v=mul_X(X[n],v1.v);
         _su3_multiply(v3.v,Q,v2.v);

         z.re+=v3.r[2*i];
      }

      printf("mu,nu = %d,%d: %.2e\n",
             mu,nu,fabs(2.0*z.re-8.0*z.im));
   }

   ranlxd(s1.r,24);
   ranlxd(s2.r,24);

   prod2xt(&s1.s,&s2.s,X);

   printf("\n");
   printf("prod2xt:\n");

   for (n=0;n<6;n++)
   {
      mu=pln[n][0];
      nu=pln[n][1];

      random_su3_dble(&Q);
      z.im=0.0;

      s3.s=mul_sigma(mu,nu,s2.s);
      s4.s=mul_gamma(5,s3.s);
      s3.s=mul_Fhat(Q,s4.s);

      z.im =_vector_prod_im(s1.s.c1,s3.s.c1);
      z.im+=_vector_prod_im(s1.s.c2,s3.s.c2);
      z.im+=_vector_prod_im(s1.s.c3,s3.s.c3);
      z.im+=_vector_prod_im(s1.s.c4,s3.s.c4);

      z.re=0.0;

      for (i=0;i<3;i++)
      {
         v1.v=vd0;
         v1.r[2*i]=1.0;

         v2.v=mul_X(X[n],v1.v);
         _su3_multiply(v3.v,Q,v2.v);

         z.re+=v3.r[2*i];
      }

      printf("mu,nu = %d,%d: %.2e\n",
             mu,nu,fabs(2.0*z.re+16.0*z.im));
   }

   printf("\n");
   exit(0);
}
