
/*******************************************************************************
*
* File xtensor.c
*
* Copyright (C) 2011, 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Spin parts of the quark force.
*
* The externally accessible functions are
*
*   u3_alg_dble **xtensor(void)
*     Returns the pointers xt[0],..,xt[5] to the X tensor field components
*     with Lorentz indices (0,1),(0,2),(0,3),(2,3),(3,1),(1,2). The arrays
*     are automatically allocated and initialized to zero if they are not
*     already allocated.
*
*   void set_xt2zero(void)
*     Sets the X tensor field to zero.
*
*   int add_det2xt(double c,ptset_t set)
*     Computes the spin part of the SW force deriving from the action
*     -Tr{ln(D)}, where D=D_ee,D_oo,D_ee+D_oo or 1 when set=EVEN_PTS,
*     ODD_PTS,ALL_PTS or NO_PTS (see the notes). The calculated matrices
*     are then multiplied by c and are added to the X tensor field. When
*     needed, the program recomputes and inverts the SW term. The program
*     returns 0 if all inversions were safe and a non-zero value otherwise.
*
*   void add_prod2xt(double c,spinor_dble *r,spinor_dble *s)
*     Computes the spin part of the SW force deriving from the "action"
*     -2*Re(r,gamma_5*Dw*s), where Dw denotes the lattice Dirac operator
*     (see the notes). The calculated matrices are then multiplied by c
*     and are added to the X tensor field.
*
*   su3_dble *xvector(void)
*     Returns the pointer xv to the X vector field. The components of
*     field are stored in memory in the same order as the link variables.
*     The array automatically allocated and initialized to zero if it is
*     not already allocated.
*
*   void set_xv2zero(void)
*     Sets the X vector field to zero.
*
*   void add_prod2xv(double c,spinor_dble *r,spinor_dble *s)
*     Computes the spin part of the force deriving from the hopping terms
*     in the "action" -2*Re(r,gamma_5*Dw*s), where Dw denotes the lattice
*     Dirac operator (see the notes). The calculated matrices are then
*     multiplied by c and are added to the X vector field.
*
* Notes:
*
* The computation of the quark forces is described in "Molecular-dynamics
* quark forces" (doc/forces.pdf). For unexplained notation concerning the
* SW term see "Implementation of the lattice Dirac operator" (doc/dirac.pdf).
*
* The SW contribution to the n'th component of the X tensor at the point x
* is given by
*
*  X[n]=i*tr{sigma_{mu,nu}*M(x)^(-1)}
*
* where M(x)=is the 12x12 matrix representing the SW term at x and n labels
* the (mu,nu)=(0,1),(0,2),(0,3),(2,3),(3,1),(1,2) index pairs. Similarly, for
* given spinor fields r and s, the associated X tensor at x is defined by
*
*  X[n]=i*tr{[gamma_5*sigma_{mu,nu}*s(x) x r^dag(x)]+(s<->r)}
*
* The contribution of the fields r,s to the X vector component on the link
* (x,x+mu) is given by
*
*  X=tr{[gamma_5*(1-gamma_mu)*s(x+mu) x r^dag(x)]+(s<->r)}
*
* In all cases, the trace is taken over the Dirac indices only.
*
* The components of the X tensor field are of type u3_alg_dble. As in the
* case of symmetric gauge-field tensor, the field array includes additional
* space for the field components on the boundaries of the local lattice
* (see tcharge/ftensor.c and lattice/README.ftidx). The type u3_alg_dble
* is explained in the module su3fcts/su3prod.c.
*
* The programs in this module may perform global operations and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define XTENSOR_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "utils.h"
#include "lattice.h"
#include "sw_term.h"
#include "sflds.h"
#include "linalg.h"
#include "forces.h"
#include "global.h"

#define _u3_alg_mul_add_assign(r,c,s) \
   (r).c1+=(c)*(s).c1; \
   (r).c2+=(c)*(s).c2; \
   (r).c3+=(c)*(s).c3; \
   (r).c4+=(c)*(s).c4; \
   (r).c5+=(c)*(s).c5; \
   (r).c6+=(c)*(s).c6; \
   (r).c7+=(c)*(s).c7; \
   (r).c8+=(c)*(s).c8; \
   (r).c9+=(c)*(s).c9

#define _su3_mul_add_assign(r,c,s) \
   (r).c11.re+=(c)*(s).c11.re; \
   (r).c11.im+=(c)*(s).c11.im; \
   (r).c12.re+=(c)*(s).c12.re; \
   (r).c12.im+=(c)*(s).c12.im; \
   (r).c13.re+=(c)*(s).c13.re; \
   (r).c13.im+=(c)*(s).c13.im; \
   (r).c21.re+=(c)*(s).c21.re; \
   (r).c21.im+=(c)*(s).c21.im; \
   (r).c22.re+=(c)*(s).c22.re; \
   (r).c22.im+=(c)*(s).c22.im; \
   (r).c23.re+=(c)*(s).c23.re; \
   (r).c23.im+=(c)*(s).c23.im; \
   (r).c31.re+=(c)*(s).c31.re; \
   (r).c31.im+=(c)*(s).c31.im; \
   (r).c32.re+=(c)*(s).c32.re; \
   (r).c32.im+=(c)*(s).c32.im; \
   (r).c33.re+=(c)*(s).c33.re; \
   (r).c33.im+=(c)*(s).c33.im

static u3_alg_dble X[6];
static u3_alg_dble **xts=NULL,**xt;
static const su3_dble ud0={{0.0,0.0},{0.0,0.0},{0.0,0.0},
                           {0.0,0.0},{0.0,0.0},{0.0,0.0},
                           {0.0,0.0},{0.0,0.0},{0.0,0.0}};
static su3_dble w ALIGNED16;
static su3_dble *xvs=NULL,*xv;


static void alloc_xts(void)
{
   int n,nt,nxt[6];
   u3_alg_dble **pp,*p;
   ftidx_t *idx;

   idx=ftidx();
   nt=0;

   for (n=0;n<6;n++)
   {
      nxt[n]=VOLUME+idx[n].nft[0]+idx[n].nft[1];
      nt+=nxt[n];
   }

   pp=malloc(12*sizeof(*pp));
   p=amalloc(nt*sizeof(*p),ALIGN);
   error((pp==NULL)||(p==NULL),1,"alloc_xts [xtensor.c]",
         "Unable to allocate field arrays");

   set_ualg2zero(nt,p);
   xts=pp;
   xt=pp+6;

   for (n=0;n<6;n++)
   {
      (*pp)=p;
      pp+=1;
      p+=nxt[n];
   }
}


u3_alg_dble **xtensor(void)
{
   int n;

   if (xts==NULL)
      alloc_xts();

   for (n=0;n<6;n++)
      xt[n]=xts[n];

   return xt;
}


void set_xt2zero(void)
{
   int n;

   if (xts==NULL)
      alloc_xts();
   else
   {
      for (n=0;n<6;n++)
         set_ualg2zero(VOLUME,xts[n]);
   }
}


int add_det2xt(double c,ptset_t set)
{
   int n,ifail;
   pauli_dble *m,*mm;

   if (set==NO_PTS)
      return 0;

   ifail=sw_term(set);

   if (ifail!=0)
      return ifail;

   if (xts==NULL)
      alloc_xts();

   if (set==ODD_PTS)
   {
      for (n=0;n<6;n++)
         xt[n]=xts[n]+(VOLUME/2);

      m=swdfld()+VOLUME;
   }
   else
   {
      for (n=0;n<6;n++)
         xt[n]=xts[n];

      m=swdfld();
   }

   if (set==ALL_PTS)
      mm=m+(2*VOLUME);
   else
      mm=m+VOLUME;

   for (;m<mm;m+=2)
   {
      det2xt(m,X);

      for (n=0;n<6;n++)
      {
         _u3_alg_mul_add_assign(xt[n][0],c,X[n]);
         xt[n]+=1;
      }
   }

   return 0;
}


void add_prod2xt(double c,spinor_dble *r,spinor_dble *s)
{
   int n;
   spinor_dble *rm;

   if (xts==NULL)
      alloc_xts();

   for (n=0;n<6;n++)
      xt[n]=xts[n];

   rm=r+VOLUME;

   for (;r<rm;r++)
   {
      prod2xt(r,s,X);

      for (n=0;n<6;n++)
      {
         _u3_alg_mul_add_assign(xt[n][0],c,X[n]);
         xt[n]+=1;
      }

      s+=1;
   }
}


static void alloc_xvs(void)
{
   int ix;

   xvs=amalloc(4*VOLUME*sizeof(*xv),ALIGN);
   error(xvs==NULL,1,"alloc_xvs [xtensor.c]",
         "Unable to allocate field array");

   for (ix=0;ix<(4*VOLUME);ix++)
      xvs[ix]=ud0;
}


su3_dble *xvector(void)
{
   if (xvs==NULL)
      alloc_xvs();

   return xvs;
}


void set_xv2zero(void)
{
   int ix;

   if (xvs==NULL)
      alloc_xvs();
   else
   {
      for (ix=0;ix<(4*VOLUME);ix++)
         xvs[ix]=ud0;
   }
}


void add_prod2xv(double c,spinor_dble *r,spinor_dble *s)
{
   int mu,*piup,*pidn;
   su3_dble *xvm;
   spinor_dble *ro,*so;

   if (xvs==NULL)
      alloc_xvs();

   cpsd_int_bnd(0x1,r);
   cpsd_int_bnd(0x1,s);

   piup=iup[VOLUME/2];
   pidn=idn[VOLUME/2];

   ro=r+(VOLUME/2);
   so=s+(VOLUME/2);

   xv=xvs;
   xvm=xv+4*VOLUME;

   while (xv<xvm)
   {
      for (mu=0;mu<4;mu++)
      {
         prod2xv[mu](ro,r+(*piup),so,s+(*piup),&w);
         _su3_mul_add_assign(*xv,c,w);

         piup+=1;
         xv+=1;

         prod2xv[mu](r+(*pidn),ro,s+(*pidn),so,&w);
         _su3_mul_add_assign(*xv,c,w);

         pidn+=1;
         xv+=1;
      }

      ro+=1;
      so+=1;
   }
}
