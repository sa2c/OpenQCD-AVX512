
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of prod2xv
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

#define _re(z,w) ((z).re*(w).re+(z).im*(w).im)
#define _im(z,w) ((z).im*(w).re-(z).re*(w).im)

typedef union
{
   su3_dble u;
   complex_dble c[9];
} umat_t;

static su3_dble u ALIGNED16;
static su3_dble v ALIGNED16;
static spinor_dble rx ALIGNED16;
static spinor_dble ry ALIGNED16;
static spinor_dble sx ALIGNED16;
static spinor_dble sy ALIGNED16;
static spinor_dble sw ALIGNED16;


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


static void add_tensor(su3_vector_dble *r,su3_vector_dble *s,su3_dble *p)
{
   (*p).c11.re+=_re((*r).c1,(*s).c1);
   (*p).c11.im+=_im((*r).c1,(*s).c1);
   (*p).c12.re+=_re((*r).c1,(*s).c2);
   (*p).c12.im+=_im((*r).c1,(*s).c2);
   (*p).c13.re+=_re((*r).c1,(*s).c3);
   (*p).c13.im+=_im((*r).c1,(*s).c3);

   (*p).c21.re+=_re((*r).c2,(*s).c1);
   (*p).c21.im+=_im((*r).c2,(*s).c1);
   (*p).c22.re+=_re((*r).c2,(*s).c2);
   (*p).c22.im+=_im((*r).c2,(*s).c2);
   (*p).c23.re+=_re((*r).c2,(*s).c3);
   (*p).c23.im+=_im((*r).c2,(*s).c3);

   (*p).c31.re+=_re((*r).c3,(*s).c1);
   (*p).c31.im+=_im((*r).c3,(*s).c1);
   (*p).c32.re+=_re((*r).c3,(*s).c2);
   (*p).c32.im+=_im((*r).c3,(*s).c2);
   (*p).c33.re+=_re((*r).c3,(*s).c3);
   (*p).c33.im+=_im((*r).c3,(*s).c3);
}


static double max_dev(su3_dble *u,su3_dble *v)
{
   int i;
   double nrm,dev;
   umat_t uu,uv;

   uu.u=(*u);
   uv.u=(*v);

   nrm=0.0;
   dev=0.0;

   for (i=0;i<9;i++)
   {
      nrm+=uu.c[i].re*uu.c[i].re+uu.c[i].im*uu.c[i].im;

      dev+=(uu.c[i].re-uv.c[i].re)*(uu.c[i].re-uv.c[i].re)+
           (uu.c[i].im-uv.c[i].im)*(uu.c[i].im-uv.c[i].im);
   }

   return sqrt(dev/nrm);
}


int main(void)
{
   int mu;

   printf("\n");
   printf("Check of prod2xv\n");
   printf("-----------------\n\n");

   rlxd_init(1,567);

   gauss_dble((double*)(&rx),24);
   gauss_dble((double*)(&ry),24);
   gauss_dble((double*)(&sx),24);
   gauss_dble((double*)(&sy),24);

   for (mu=0;mu<4;mu++)
   {
      prod2xv[mu](&rx,&ry,&sx,&sy,&u);
      cm3x3_zero(1,&v);

      sw=mul_gamma(mu,ry);
      _vector_sub(sw.c1,ry.c1,sw.c1);
      _vector_sub(sw.c2,ry.c2,sw.c2);
      _vector_sub(sw.c3,ry.c3,sw.c3);
      _vector_sub(sw.c4,ry.c4,sw.c4);
      sw=mul_gamma(5,sw);

      add_tensor(&sw.c1,&sx.c1,&v);
      add_tensor(&sw.c2,&sx.c2,&v);
      add_tensor(&sw.c3,&sx.c3,&v);
      add_tensor(&sw.c4,&sx.c4,&v);

      sw=mul_gamma(mu,sy);
      _vector_sub(sw.c1,sy.c1,sw.c1);
      _vector_sub(sw.c2,sy.c2,sw.c2);
      _vector_sub(sw.c3,sy.c3,sw.c3);
      _vector_sub(sw.c4,sy.c4,sw.c4);
      sw=mul_gamma(5,sw);

      add_tensor(&sw.c1,&rx.c1,&v);
      add_tensor(&sw.c2,&rx.c2,&v);
      add_tensor(&sw.c3,&rx.c3,&v);
      add_tensor(&sw.c4,&rx.c4,&v);

      printf("mu = %d: %.2e\n",mu,max_dev(&u,&v));
   }

   printf("\n");
   exit(0);
}
