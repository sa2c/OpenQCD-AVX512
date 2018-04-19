
/*******************************************************************************
*
* File frcfcts.c
*
* Copyright (C) 2005, 2011, 2012, 2016 Martin Luescher, Stefan Schaefer
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic functions used for the force calculation.
*
* The externally accessible functions are
*
*   void det2xt(pauli_dble *m,u3_alg_dble *X)
*     Computes the matrices X[0],..,X[5] associated to the SW term on a
*     given lattice point (see the notes). The program expects that m[0]
*     and m[1] contain the hermitian part of the inverse of the SW term
*     at the chosen point.
*
*   void prod2xt(spinor_dble *r,spinor_dble *s,u3_alg_dble *X)
*     Computes the matrices X[0],..,X[5] associated to a pair of spinors
*     r and s at a given lattice point (see the notes).
*
* The following is an array of functions indexed by the direction mu=0,..,3:
*
*   void (*prod2xv[])(spinor_dble *rx,spinor_dble *ry,
*                      spinor_dble *sx,spinor_dble *sy,su3_dble *u)
*     Computes the complex 3x3 matrix
*
*       u=tr{gamma_5*(1-gamma_mu)*[(sy x rx^dag)+(ry x sx^dag)]}
*
*     where ..x.. denotes the tensor product in spinor space and the trace
*     is taken over the Dirac indices.
*
* Notes:
*
* As discussed in the notes
*
*  M. Luescher: "Molecular-dynamics quark forces" (January 2012)
*
* the programs in this module serve to compute the spin part of the quark
* forces. The data type u3_alg_dble is described at the top of the module
* su3fcts/su3prod.c.
*
* The matrices computed by the program det2xt() are
*
*  X[n]=i*tr{sigma_{mu,nu}*diag(m[0],m[1])}
*
* where (mu,nu)=(0,1),(0,2),(0,3),(2,3),(3,1),(1,2) for n=0,..,5. Similarly,
* the program prod2xt() computes
*
*  X[n]=i*tr{(gamma_5*sigma_{mu,nu}*s) x (r^dag)+(s<->r)}
*
* where ..x.. denotes the tensor product in spinor space. In both cases,
* the trace is taken over the Dirac indices only.
*
*******************************************************************************/

#define FRCFCTS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "forces.h"

#define _re(z,w) ((z).re*(w).re+(z).im*(w).im)
#define _im(z,w) ((z).im*(w).re-(z).re*(w).im)

typedef union
{
   spinor_dble s;
   weyl_dble w[2];
} spin_t;

static su3_vector_dble psi1 ALIGNED16;
static su3_vector_dble psi2 ALIGNED16;
static su3_vector_dble chi1 ALIGNED16;
static su3_vector_dble chi2 ALIGNED16;
static pauli_dble ms[2] ALIGNED16;


void det2xt(pauli_dble *m,u3_alg_dble *X)
{
   double x,*up,*um;
   u3_alg_dble *X0,*X1;

   up=m[0].u;
   um=m[1].u;

   X0=X;
   X1=X+3;

   x=up[10]+up[10];
   (*X0).c1=x;
   (*X1).c1=-x;
   x=um[10]+um[10];
   (*X0).c1-=x;
   (*X1).c1-=x;

   x=up[20]+up[20];
   (*X0).c2=x;
   (*X1).c2=-x;
   x=um[20]+um[20];
   (*X0).c2-=x;
   (*X1).c2-=x;

   x=up[28]+up[28];
   (*X0).c3=x;
   (*X1).c3=-x;
   x=um[28]+um[28];
   (*X0).c3-=x;
   (*X1).c3-=x;

   x=up[19]-up[13];
   (*X0).c4=x;
   (*X1).c4=-x;
   x=um[19]-um[13];
   (*X0).c4-=x;
   (*X1).c4-=x;

   x=up[12]+up[18];
   (*X0).c5=x;
   (*X1).c5=-x;
   x=um[12]+um[18];
   (*X0).c5-=x;
   (*X1).c5-=x;

   x=up[25]-up[15];
   (*X0).c6=x;
   (*X1).c6=-x;
   x=um[25]-um[15];
   (*X0).c6-=x;
   (*X1).c6-=x;

   x=up[14]+up[24];
   (*X0).c7=x;
   (*X1).c7=-x;
   x=um[14]+um[24];
   (*X0).c7-=x;
   (*X1).c7-=x;

   x=up[27]-up[23];
   (*X0).c8=x;
   (*X1).c8=-x;
   x=um[27]-um[23];
   (*X0).c8-=x;
   (*X1).c8-=x;

   x=up[22]+up[26];
   (*X0).c9=x;
   (*X1).c9=-x;
   x=um[22]+um[26];
   (*X0).c9-=x;
   (*X1).c9-=x;

   X0=X+1;
   X1=X+4;

   x=up[11]+up[11];
   (*X0).c1=-x;
   (*X1).c1=x;
   x=um[11]+um[11];
   (*X0).c1+=x;
   (*X1).c1+=x;

   x=up[21]+up[21];
   (*X0).c2=-x;
   (*X1).c2=x;
   x=um[21]+um[21];
   (*X0).c2+=x;
   (*X1).c2+=x;

   x=up[29]+up[29];
   (*X0).c3=-x;
   (*X1).c3=x;
   x=um[29]+um[29];
   (*X0).c3+=x;
   (*X1).c3+=x;

   x=up[18]-up[12];
   (*X0).c4=x;
   (*X1).c4=-x;
   x=um[18]-um[12];
   (*X0).c4-=x;
   (*X1).c4-=x;

   x=up[13]+up[19];
   (*X0).c5=-x;
   (*X1).c5=x;
   x=um[13]+um[19];
   (*X0).c5+=x;
   (*X1).c5+=x;

   x=up[24]-up[14];
   (*X0).c6=x;
   (*X1).c6=-x;
   x=um[24]-um[14];
   (*X0).c6-=x;
   (*X1).c6-=x;

   x=up[25]+up[15];
   (*X0).c7=-x;
   (*X1).c7=x;
   x=um[25]+um[15];
   (*X0).c7+=x;
   (*X1).c7+=x;

   x=up[26]-up[22];
   (*X0).c8=x;
   (*X1).c8=-x;
   x=um[26]-um[22];
   (*X0).c8-=x;
   (*X1).c8-=x;

   x=up[27]+up[23];
   (*X0).c9=-x;
   (*X1).c9=x;
   x=um[27]+um[23];
   (*X0).c9+=x;
   (*X1).c9+=x;

   X0=X+2;
   X1=X+5;

   x=up[0]-up[3];
   (*X0).c1=x;
   (*X1).c1=-x;
   x=um[0]-um[3];
   (*X0).c1-=x;
   (*X1).c1-=x;

   x=up[1]-up[4];
   (*X0).c2=x;
   (*X1).c2=-x;
   x=um[1]-um[4];
   (*X0).c2-=x;
   (*X1).c2-=x;

   x=up[2]-up[5];
   (*X0).c3=x;
   (*X1).c3=-x;
   x=um[2]-um[5];
   (*X0).c3-=x;
   (*X1).c3-=x;

   x=up[31]-up[7];
   (*X0).c4=x;
   (*X1).c4=-x;
   x=um[31]-um[7];
   (*X0).c4-=x;
   (*X1).c4-=x;

   x=up[6]-up[30];
   (*X0).c5=x;
   (*X1).c5=-x;
   x=um[6]-um[30];
   (*X0).c5-=x;
   (*X1).c5-=x;

   x=up[33]-up[9];
   (*X0).c6=x;
   (*X1).c6=-x;
   x=um[33]-um[9];
   (*X0).c6-=x;
   (*X1).c6-=x;

   x=up[8]-up[32];
   (*X0).c7=x;
   (*X1).c7=-x;
   x=um[8]-um[32];
   (*X0).c7-=x;
   (*X1).c7-=x;

   x=up[35]-up[17];
   (*X0).c8=x;
   (*X1).c8=-x;
   x=um[35]-um[17];
   (*X0).c8-=x;
   (*X1).c8-=x;

   x=up[16]-up[34];
   (*X0).c9=x;
   (*X1).c9=-x;
   x=um[16]-um[34];
   (*X0).c9-=x;
   (*X1).c9-=x;
}


static void det2xt5(pauli_dble *m,u3_alg_dble *X)
{
   double x,*up,*um;
   u3_alg_dble *X0,*X1;

   up=m[0].u;
   um=m[1].u;

   X0=X;
   X1=X+3;

   x=up[10]+up[10];
   (*X0).c1=x;
   (*X1).c1=-x;
   x=um[10]+um[10];
   (*X0).c1+=x;
   (*X1).c1+=x;

   x=up[20]+up[20];
   (*X0).c2=x;
   (*X1).c2=-x;
   x=um[20]+um[20];
   (*X0).c2+=x;
   (*X1).c2+=x;

   x=up[28]+up[28];
   (*X0).c3=x;
   (*X1).c3=-x;
   x=um[28]+um[28];
   (*X0).c3+=x;
   (*X1).c3+=x;

   x=up[19]-up[13];
   (*X0).c4=x;
   (*X1).c4=-x;
   x=um[19]-um[13];
   (*X0).c4+=x;
   (*X1).c4+=x;

   x=up[12]+up[18];
   (*X0).c5=x;
   (*X1).c5=-x;
   x=um[12]+um[18];
   (*X0).c5+=x;
   (*X1).c5+=x;

   x=up[25]-up[15];
   (*X0).c6=x;
   (*X1).c6=-x;
   x=um[25]-um[15];
   (*X0).c6+=x;
   (*X1).c6+=x;

   x=up[14]+up[24];
   (*X0).c7=x;
   (*X1).c7=-x;
   x=um[14]+um[24];
   (*X0).c7+=x;
   (*X1).c7+=x;

   x=up[27]-up[23];
   (*X0).c8=x;
   (*X1).c8=-x;
   x=um[27]-um[23];
   (*X0).c8+=x;
   (*X1).c8+=x;

   x=up[22]+up[26];
   (*X0).c9=x;
   (*X1).c9=-x;
   x=um[22]+um[26];
   (*X0).c9+=x;
   (*X1).c9+=x;

   X0=X+1;
   X1=X+4;

   x=up[11]+up[11];
   (*X0).c1=-x;
   (*X1).c1=x;
   x=um[11]+um[11];
   (*X0).c1-=x;
   (*X1).c1-=x;

   x=up[21]+up[21];
   (*X0).c2=-x;
   (*X1).c2=x;
   x=um[21]+um[21];
   (*X0).c2-=x;
   (*X1).c2-=x;

   x=up[29]+up[29];
   (*X0).c3=-x;
   (*X1).c3=x;
   x=um[29]+um[29];
   (*X0).c3-=x;
   (*X1).c3-=x;

   x=up[18]-up[12];
   (*X0).c4=x;
   (*X1).c4=-x;
   x=um[18]-um[12];
   (*X0).c4+=x;
   (*X1).c4+=x;

   x=up[13]+up[19];
   (*X0).c5=-x;
   (*X1).c5=x;
   x=um[13]+um[19];
   (*X0).c5-=x;
   (*X1).c5-=x;

   x=up[24]-up[14];
   (*X0).c6=x;
   (*X1).c6=-x;
   x=um[24]-um[14];
   (*X0).c6+=x;
   (*X1).c6+=x;

   x=up[25]+up[15];
   (*X0).c7=-x;
   (*X1).c7=x;
   x=um[25]+um[15];
   (*X0).c7-=x;
   (*X1).c7-=x;

   x=up[26]-up[22];
   (*X0).c8=x;
   (*X1).c8=-x;
   x=um[26]-um[22];
   (*X0).c8+=x;
   (*X1).c8+=x;

   x=up[27]+up[23];
   (*X0).c9=-x;
   (*X1).c9=x;
   x=um[27]+um[23];
   (*X0).c9-=x;
   (*X1).c9-=x;

   X0=X+2;
   X1=X+5;

   x=up[0]-up[3];
   (*X0).c1=x;
   (*X1).c1=-x;
   x=um[0]-um[3];
   (*X0).c1+=x;
   (*X1).c1+=x;

   x=up[1]-up[4];
   (*X0).c2=x;
   (*X1).c2=-x;
   x=um[1]-um[4];
   (*X0).c2+=x;
   (*X1).c2+=x;

   x=up[2]-up[5];
   (*X0).c3=x;
   (*X1).c3=-x;
   x=um[2]-um[5];
   (*X0).c3+=x;
   (*X1).c3+=x;

   x=up[31]-up[7];
   (*X0).c4=x;
   (*X1).c4=-x;
   x=um[31]-um[7];
   (*X0).c4+=x;
   (*X1).c4+=x;

   x=up[6]-up[30];
   (*X0).c5=x;
   (*X1).c5=-x;
   x=um[6]-um[30];
   (*X0).c5+=x;
   (*X1).c5+=x;

   x=up[33]-up[9];
   (*X0).c6=x;
   (*X1).c6=-x;
   x=um[33]-um[9];
   (*X0).c6+=x;
   (*X1).c6+=x;

   x=up[8]-up[32];
   (*X0).c7=x;
   (*X1).c7=-x;
   x=um[8]-um[32];
   (*X0).c7+=x;
   (*X1).c7+=x;

   x=up[35]-up[17];
   (*X0).c8=x;
   (*X1).c8=-x;
   x=um[35]-um[17];
   (*X0).c8+=x;
   (*X1).c8+=x;

   x=up[16]-up[34];
   (*X0).c9=x;
   (*X1).c9=-x;
   x=um[16]-um[34];
   (*X0).c9+=x;
   (*X1).c9+=x;
}


static void vec2pauli(weyl_dble *r,weyl_dble *s,pauli_dble *m)
{
   double *u;
   su3_vector_dble *r1,*r2,*s1,*s2;

   u=(*m).u;
   r1=&((*r).c1);
   r2=&((*r).c2);
   s1=&((*s).c1);
   s2=&((*s).c2);

   u[ 0]=_re((*s1).c1,(*r1).c1)+_re((*s1).c1,(*r1).c1);
   u[ 1]=_re((*s1).c2,(*r1).c2)+_re((*s1).c2,(*r1).c2);
   u[ 2]=_re((*s1).c3,(*r1).c3)+_re((*s1).c3,(*r1).c3);

   u[ 3]=_re((*s2).c1,(*r2).c1)+_re((*s2).c1,(*r2).c1);
   u[ 4]=_re((*s2).c2,(*r2).c2)+_re((*s2).c2,(*r2).c2);
   u[ 5]=_re((*s2).c3,(*r2).c3)+_re((*s2).c3,(*r2).c3);

   u[ 6]=_re((*s1).c1,(*r1).c2)+_re((*r1).c1,(*s1).c2);
   u[ 7]=_im((*s1).c1,(*r1).c2)+_im((*r1).c1,(*s1).c2);
   u[ 8]=_re((*s1).c1,(*r1).c3)+_re((*r1).c1,(*s1).c3);
   u[ 9]=_im((*s1).c1,(*r1).c3)+_im((*r1).c1,(*s1).c3);

   u[10]=_re((*s1).c1,(*r2).c1)+_re((*r1).c1,(*s2).c1);
   u[11]=_im((*s1).c1,(*r2).c1)+_im((*r1).c1,(*s2).c1);
   u[12]=_re((*s1).c1,(*r2).c2)+_re((*r1).c1,(*s2).c2);
   u[13]=_im((*s1).c1,(*r2).c2)+_im((*r1).c1,(*s2).c2);
   u[14]=_re((*s1).c1,(*r2).c3)+_re((*r1).c1,(*s2).c3);
   u[15]=_im((*s1).c1,(*r2).c3)+_im((*r1).c1,(*s2).c3);

   u[16]=_re((*s1).c2,(*r1).c3)+_re((*r1).c2,(*s1).c3);
   u[17]=_im((*s1).c2,(*r1).c3)+_im((*r1).c2,(*s1).c3);

   u[18]=_re((*s1).c2,(*r2).c1)+_re((*r1).c2,(*s2).c1);
   u[19]=_im((*s1).c2,(*r2).c1)+_im((*r1).c2,(*s2).c1);
   u[20]=_re((*s1).c2,(*r2).c2)+_re((*r1).c2,(*s2).c2);
   u[21]=_im((*s1).c2,(*r2).c2)+_im((*r1).c2,(*s2).c2);
   u[22]=_re((*s1).c2,(*r2).c3)+_re((*r1).c2,(*s2).c3);
   u[23]=_im((*s1).c2,(*r2).c3)+_im((*r1).c2,(*s2).c3);

   u[24]=_re((*s1).c3,(*r2).c1)+_re((*r1).c3,(*s2).c1);
   u[25]=_im((*s1).c3,(*r2).c1)+_im((*r1).c3,(*s2).c1);
   u[26]=_re((*s1).c3,(*r2).c2)+_re((*r1).c3,(*s2).c2);
   u[27]=_im((*s1).c3,(*r2).c2)+_im((*r1).c3,(*s2).c2);
   u[28]=_re((*s1).c3,(*r2).c3)+_re((*r1).c3,(*s2).c3);
   u[29]=_im((*s1).c3,(*r2).c3)+_im((*r1).c3,(*s2).c3);

   u[30]=_re((*s2).c1,(*r2).c2)+_re((*r2).c1,(*s2).c2);
   u[31]=_im((*s2).c1,(*r2).c2)+_im((*r2).c1,(*s2).c2);
   u[32]=_re((*s2).c1,(*r2).c3)+_re((*r2).c1,(*s2).c3);
   u[33]=_im((*s2).c1,(*r2).c3)+_im((*r2).c1,(*s2).c3);
   u[34]=_re((*s2).c2,(*r2).c3)+_re((*r2).c2,(*s2).c3);
   u[35]=_im((*s2).c2,(*r2).c3)+_im((*r2).c2,(*s2).c3);
}


void prod2xt(spinor_dble *r,spinor_dble *s,u3_alg_dble *X)
{
   spin_t *spr,*sps;

   spr=(spin_t*)(r);
   sps=(spin_t*)(s);

   vec2pauli((*spr).w,(*sps).w,ms);
   vec2pauli((*spr).w+1,(*sps).w+1,ms+1);

   det2xt5(ms,X);
}


static void set2mat(su3_dble *u)
{
   (*u).c11.re=_re(psi1.c1,chi1.c1)+_re(psi2.c1,chi2.c1);
   (*u).c11.im=_im(psi1.c1,chi1.c1)+_im(psi2.c1,chi2.c1);
   (*u).c12.re=_re(psi1.c1,chi1.c2)+_re(psi2.c1,chi2.c2);
   (*u).c12.im=_im(psi1.c1,chi1.c2)+_im(psi2.c1,chi2.c2);
   (*u).c13.re=_re(psi1.c1,chi1.c3)+_re(psi2.c1,chi2.c3);
   (*u).c13.im=_im(psi1.c1,chi1.c3)+_im(psi2.c1,chi2.c3);

   (*u).c21.re=_re(psi1.c2,chi1.c1)+_re(psi2.c2,chi2.c1);
   (*u).c21.im=_im(psi1.c2,chi1.c1)+_im(psi2.c2,chi2.c1);
   (*u).c22.re=_re(psi1.c2,chi1.c2)+_re(psi2.c2,chi2.c2);
   (*u).c22.im=_im(psi1.c2,chi1.c2)+_im(psi2.c2,chi2.c2);
   (*u).c23.re=_re(psi1.c2,chi1.c3)+_re(psi2.c2,chi2.c3);
   (*u).c23.im=_im(psi1.c2,chi1.c3)+_im(psi2.c2,chi2.c3);

   (*u).c31.re=_re(psi1.c3,chi1.c1)+_re(psi2.c3,chi2.c1);
   (*u).c31.im=_im(psi1.c3,chi1.c1)+_im(psi2.c3,chi2.c1);
   (*u).c32.re=_re(psi1.c3,chi1.c2)+_re(psi2.c3,chi2.c2);
   (*u).c32.im=_im(psi1.c3,chi1.c2)+_im(psi2.c3,chi2.c2);
   (*u).c33.re=_re(psi1.c3,chi1.c3)+_re(psi2.c3,chi2.c3);
   (*u).c33.im=_im(psi1.c3,chi1.c3)+_im(psi2.c3,chi2.c3);
}


static void add2mat(su3_dble *u)
{
   (*u).c11.re+=_re(psi1.c1,chi1.c1)+_re(psi2.c1,chi2.c1);
   (*u).c11.im+=_im(psi1.c1,chi1.c1)+_im(psi2.c1,chi2.c1);
   (*u).c12.re+=_re(psi1.c1,chi1.c2)+_re(psi2.c1,chi2.c2);
   (*u).c12.im+=_im(psi1.c1,chi1.c2)+_im(psi2.c1,chi2.c2);
   (*u).c13.re+=_re(psi1.c1,chi1.c3)+_re(psi2.c1,chi2.c3);
   (*u).c13.im+=_im(psi1.c1,chi1.c3)+_im(psi2.c1,chi2.c3);

   (*u).c21.re+=_re(psi1.c2,chi1.c1)+_re(psi2.c2,chi2.c1);
   (*u).c21.im+=_im(psi1.c2,chi1.c1)+_im(psi2.c2,chi2.c1);
   (*u).c22.re+=_re(psi1.c2,chi1.c2)+_re(psi2.c2,chi2.c2);
   (*u).c22.im+=_im(psi1.c2,chi1.c2)+_im(psi2.c2,chi2.c2);
   (*u).c23.re+=_re(psi1.c2,chi1.c3)+_re(psi2.c2,chi2.c3);
   (*u).c23.im+=_im(psi1.c2,chi1.c3)+_im(psi2.c2,chi2.c3);

   (*u).c31.re+=_re(psi1.c3,chi1.c1)+_re(psi2.c3,chi2.c1);
   (*u).c31.im+=_im(psi1.c3,chi1.c1)+_im(psi2.c3,chi2.c1);
   (*u).c32.re+=_re(psi1.c3,chi1.c2)+_re(psi2.c3,chi2.c2);
   (*u).c32.im+=_im(psi1.c3,chi1.c2)+_im(psi2.c3,chi2.c2);
   (*u).c33.re+=_re(psi1.c3,chi1.c3)+_re(psi2.c3,chi2.c3);
   (*u).c33.im+=_im(psi1.c3,chi1.c3)+_im(psi2.c3,chi2.c3);
}


static void prod2xv0(spinor_dble *rx,spinor_dble *ry,
                      spinor_dble *sx,spinor_dble *sy,su3_dble *u)
{
   _vector_add(psi1,(*ry).c1,(*ry).c3);
   _vector_add(psi2,(*ry).c2,(*ry).c4);
   _vector_sub(chi1,(*sx).c1,(*sx).c3);
   _vector_sub(chi2,(*sx).c2,(*sx).c4);
   set2mat(u);

   _vector_add(psi1,(*sy).c1,(*sy).c3);
   _vector_add(psi2,(*sy).c2,(*sy).c4);
   _vector_sub(chi1,(*rx).c1,(*rx).c3);
   _vector_sub(chi2,(*rx).c2,(*rx).c4);
   add2mat(u);
}


static void prod2xv1(spinor_dble *rx,spinor_dble *ry,
                      spinor_dble *sx,spinor_dble *sy,su3_dble *u)
{
   _vector_i_add(psi1,(*ry).c1,(*ry).c4);
   _vector_i_add(psi2,(*ry).c2,(*ry).c3);
   _vector_i_sub(chi1,(*sx).c1,(*sx).c4);
   _vector_i_sub(chi2,(*sx).c2,(*sx).c3);
   set2mat(u);

   _vector_i_add(psi1,(*sy).c1,(*sy).c4);
   _vector_i_add(psi2,(*sy).c2,(*sy).c3);
   _vector_i_sub(chi1,(*rx).c1,(*rx).c4);
   _vector_i_sub(chi2,(*rx).c2,(*rx).c3);
   add2mat(u);
}


static void prod2xv2(spinor_dble *rx,spinor_dble *ry,
                      spinor_dble *sx,spinor_dble *sy,su3_dble *u)
{
   _vector_add(psi1,(*ry).c1,(*ry).c4);
   _vector_sub(psi2,(*ry).c2,(*ry).c3);
   _vector_sub(chi1,(*sx).c1,(*sx).c4);
   _vector_add(chi2,(*sx).c2,(*sx).c3);
   set2mat(u);

   _vector_add(psi1,(*sy).c1,(*sy).c4);
   _vector_sub(psi2,(*sy).c2,(*sy).c3);
   _vector_sub(chi1,(*rx).c1,(*rx).c4);
   _vector_add(chi2,(*rx).c2,(*rx).c3);
   add2mat(u);
}


static void prod2xv3(spinor_dble *rx,spinor_dble *ry,
                      spinor_dble *sx,spinor_dble *sy,su3_dble *u)
{
   _vector_i_add(psi1,(*ry).c1,(*ry).c3);
   _vector_i_sub(psi2,(*ry).c2,(*ry).c4);
   _vector_i_sub(chi1,(*sx).c1,(*sx).c3);
   _vector_i_add(chi2,(*sx).c2,(*sx).c4);
   set2mat(u);

   _vector_i_add(psi1,(*sy).c1,(*sy).c3);
   _vector_i_sub(psi2,(*sy).c2,(*sy).c4);
   _vector_i_sub(chi1,(*rx).c1,(*rx).c3);
   _vector_i_add(chi2,(*rx).c2,(*rx).c4);
   add2mat(u);
}


void (*prod2xv[4])(spinor_dble *rx,spinor_dble *ry,
                    spinor_dble *sx,spinor_dble *sy,su3_dble *u)=
{prod2xv0,prod2xv1,prod2xv2,prod2xv3};
