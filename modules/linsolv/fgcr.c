
/*******************************************************************************
*
* File fgcr.c
*
* Copyright (C) 2005, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic flexible GCR solver program for the lattice Dirac equation.
*
* The externally accessible function is
*
*   double fgcr(int vol,int icom,
*               void (*Dop)(spinor_dble *s,spinor_dble *r),
*               void (*Mop)(int k,spinor *rho,spinor *phi,spinor *chi),
*               spinor **ws,spinor_dble **wsd,int nkv,int nmx,double res,
*               spinor_dble *eta,spinor_dble *psi,int *status)
*     Solution of the Dirac equation D*psi=eta for given source eta, using
*     the preconditioned GCR algorithm. See the notes for the explanation
*     of the parameters of the program.
*
* Notes:
*
* This program uses single-precision arithmetic to reduce the execution
* time, but obtains the solution with double-precision accuracy.
*
* The programs Dop() and Mop() for the operator D and the preconditioner M
* are assumed to have the following properties:
*
*   void Dop(spinor_dble *s,spinor_dble *r)
*     Application of the operator D to the Dirac field s and assignment of
*     the result to r. On exit s may be changed but must satisfy D*s=r.
*
*   void Mop(int k,spinor *rho,spinor *phi,spinor *chi)
*     Approximate solution of the equation D*phi=rho in the k'th step of
*     the GCR algorithm. On exit rho is unchanged and chi=D*phi.
*
* Mop() is not required to be a linear operator and may involve an iterative
* procedure with a dynamical stopping criterion, for example. The field phi
* merely defines the next search direction and can in principle be chosen
* arbitrarily.
*
* The other parameters of the program fgcr() are:
*
*   vol     Number of spinors in the Dirac fields.
*
*   icom    Indicates whether the equation to be solved is a local
*           equation (icom=0) or a global one (icom=1). Scalar products
*           are summed over all MPI processes if icom=1, while no
*           communications are performed if icom=0.
*
*   nkv     Maximal number of Krylov vectors generated before the GCR
*           algorithm is restarted.
*
*   nmx     Maximal total number of Krylov vectors that may be generated.
*
*   res     Desired maximal relative residue |eta-D*psi|/|eta| of the
*           calculated solution.
*
*   ws      Array of at least 2*nkv+1 single-precision spinor fields
*           (used as work space).
*
*   wsd     Array of at least 1 double-precision spinor field (used
*           as work space).
*
*   eta     Source field (unchanged on exit).
*
*   psi     Calculated approximate solution of the Dirac equation
*           D*psi=eta.
*
*   status  On exit, this parameter reports the total number of Krylov
*           vectors that were generated, or a negative value if the
*           program failed.
*
* Independently of whether the program succeeds in solving the Dirac equation
* to the desired accuracy, the program returns the norm of the residue of
* the field psi.
*
* Some debugging output is printed to stdout on process 0 if FGCR_DBG is
* defined at compilation time.
*
*******************************************************************************/

#define FGCR_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "sflds.h"
#include "linalg.h"
#include "linsolv.h"
#include "global.h"

#define PRECISION_LIMIT ((double)(100.0f*FLT_EPSILON))

static int nkm=0;
static float *b;
static complex *a,*c;
static double rn;
static spinor **phi,**chi,*rho;
static spinor_dble *wrk;


static int alloc_arrays(int nkv)
{
   if (nkm>0)
   {
      afree(a);
      afree(b);
   }

   a=amalloc(nkv*(nkv+1)*sizeof(*a),ALIGN);
   b=amalloc(nkv*sizeof(*b),ALIGN);

   if ((a==NULL)||(b==NULL))
      return 1;

   c=a+nkv*nkv;
   nkm=nkv;

   return 0;
}


static void gcr_init(int vol,int icom,int nkv,spinor **ws,spinor_dble **wsd,
                     spinor_dble *eta,spinor_dble *psi)
{
   phi=ws;
   rho=ws[nkv];
   chi=ws+nkv+1;
   wrk=wsd[0];

   set_sd2zero(vol,psi);
   assign_sd2s(vol,eta,rho);

   rn=(double)(norm_square(vol,icom,rho));
   rn=sqrt(rn);
}


static void gcr_step(int vol,int icom,int k,int nkv,
                     void (*Mop)(int k,spinor *rho,spinor *phi,spinor *chi))
{
   int l;
   complex z;

   (*Mop)(k,rho,phi[k],chi[k]);

   for (l=0;l<k;l++)
   {
      a[nkv*l+k]=spinor_prod(vol,icom,chi[l],chi[k]);
      z.re=-a[nkv*l+k].re;
      z.im=-a[nkv*l+k].im;
      mulc_spinor_add(vol,chi[k],chi[l],z);
   }

   b[k]=normalize(vol,icom,chi[k]);
   c[k]=spinor_prod(vol,icom,chi[k],rho);
   z.re=-c[k].re;
   z.im=-c[k].im;
   mulc_spinor_add(vol,rho,chi[k],z);

   rn=(double)(norm_square(vol,icom,rho));
   rn=sqrt(rn);
}


static void update_psi(int vol,int icom,int k,int nkv,
                       spinor_dble *eta,spinor_dble *psi,
                       void (*Dop)(spinor_dble *s,spinor_dble *r))
{
   int l,i;
   float r;
   complex z;

   for (l=k;l>=0;l--)
   {
      z.re=c[l].re;
      z.im=c[l].im;

      for (i=(l+1);i<=k;i++)
      {
         z.re-=(a[l*nkv+i].re*c[i].re-a[l*nkv+i].im*c[i].im);
         z.im-=(a[l*nkv+i].re*c[i].im+a[l*nkv+i].im*c[i].re);
      }

      r=1.0f/b[l];
      c[l].re=z.re*r;
      c[l].im=z.im*r;
   }

   set_s2zero(vol,rho);

   for (l=k;l>=0;l--)
      mulc_spinor_add(vol,rho,phi[l],c[l]);

   add_s2sd(vol,rho,psi);
   (*Dop)(psi,wrk);
   diff_sd2s(vol,eta,wrk,rho);

   rn=(double)(norm_square(vol,icom,rho));
   rn=sqrt(rn);
}


double fgcr(int vol,int icom,
            void (*Dop)(spinor_dble *s,spinor_dble *r),
            void (*Mop)(int k,spinor *eta,spinor *psi,spinor *chi),
            spinor **ws,spinor_dble **wsd,int nkv,int nmx,double res,
            spinor_dble *eta,spinor_dble *psi,int *status)
{
   int ie,k,iprms[3];
   double rn_old,tol,dprms[1];

   if ((icom==1)&&(NPROC>1))
   {
      iprms[0]=vol;
      iprms[1]=nkv;
      iprms[2]=nmx;
      dprms[0]=res;

      MPI_Bcast(iprms,3,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      error((iprms[0]!=vol)||(iprms[1]!=nkv)||(iprms[2]!=nmx)||
            (dprms[0]!=res),1,"fgcr [fgcr.c]","Parameters are not global");

      error_root((vol<=0)||(nkv<1)||(nmx<1)||(res<=DBL_EPSILON),1,
                 "fgcr [fgcr.c]","Improper choice of vol,nkv,nmx or res");

      if (nkv>nkm)
      {
         ie=alloc_arrays(nkv);
         error(ie,1,"fgcr [fgcr.c]","Unable to allocate auxiliary arrays");
      }
   }
   else
   {
      if ((vol<=0)||(nkv<1)||(nmx<1)||(res<=DBL_EPSILON))
      {
         error_loc(1,1,"fgcr [fgcr.c]",
                   "Improper choice of vol,nkv,nmx or res");
         (*status)=0;
         return 1.0;
      }

      if (nkv>nkm)
      {
         ie=alloc_arrays(nkv);

         if (ie)
         {
            error_loc(1,1,"fgcr [fgcr.c]",
                      "Unable to allocate auxiliary arrays");
            (*status)=0;
            return 1.0;
         }
      }
   }

   gcr_init(vol,icom,nkv,ws,wsd,eta,psi);
   tol=res*rn;
   (*status)=0;

   while (rn>tol)
   {
#ifdef FGCR_DBG
      message("[fgcr]: rn_old = %.2e\n",rn);
#endif
      rn_old=rn;

      for (k=0;;k++)
      {
         gcr_step(vol,icom,k,nkv,Mop);
         (*status)+=1;
#ifdef FGCR_DBG
         message("[fgcr]: k = %d, rn = %.2e\n",k,rn);
#endif
         if ((rn<=tol)||(rn<(PRECISION_LIMIT*rn_old))||
             ((k+1)==nkv)||((*status)==nmx))
            break;
      }

      update_psi(vol,icom,k,nkv,eta,psi,Dop);

      if (((*status)==nmx)&&(rn>tol))
      {
         (*status)=-1;
         return rn;
      }
   }

   return rn;
}
