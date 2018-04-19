
/*******************************************************************************
*
* File fgcr4vd.c
*
* Copyright (C) 2007, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic flexible GCR solver program for the little Dirac equation.
*
* The externally accessible function is
*
*   double fgcr4vd(int vol,int icom,
*                  void (*Dop)(complex_dble *v,complex_dble *w),
*                  void (*Mop)(int k,complex *rho,complex *phi,complex *chi),
*                  complex *wv[],complex_dble *wvd[],int nkv,int nmx,double res,
*                  complex_dble *eta,complex_dble *psi,int *status)
*     Solution of the little equation D*psi=eta for given source eta, using
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
*   void Dop(complex_dble *v,complex_dble *w)
*     Application of the operator D to the complex field v and assignment
*     of the result to w. On exit v may be changed but must satisfy D*v=w.
*
*   void Mop(int k,complex *rho,complex *phi,complex *chi)
*     Approximate solution of the equation D*phi=rho in the k'th step of
*     the GCR algorithm. On exit rho is unchanged and chi=D*phi.
*
* Mop() is not required to be a linear operator and may involve an iterative
* procedure with a dynamical stopping criterion, for example. The field phi
* merely defines the next search direction and can in principle be chosen
* arbitrarily.
*
* The other parameters of the program fgcr4vd() are:
*
*   vol     Number of complex components of the fields on which the
*           operator D acts.
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
*   wv      Array of at least 2*nkv+1 single-precision complex fields
*           (used as work space).
*
*   wvd     Array of at least 1 double-precision complex field (used
*           as work space).
*
*   eta     Source field (unchanged on exit).
*
*   psi     Calculated approximate solution of the little equation
*           D*psi=eta.
*
*   status  On exit, this parameter reports the total number of Krylov
*           vectors that were generated or -1 if the algorithm did not
*           converge.
*
* Independently of whether the program succeeds in solving the little equation
* to the desired accuracy, the program returns the norm of the residue of
* the field psi.
*
* Some debugging output is printed to stdout on process 0 if FGCR4VD_DBG is
* defined at compilation time.
*
*******************************************************************************/

#define FGCR4VD_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "vflds.h"
#include "linalg.h"
#include "linsolv.h"
#include "global.h"

#define PRECISION_LIMIT ((double)(100.0f*FLT_EPSILON))

static int nkm=0;
static float *b;
static complex *a,*c;
static double rn;
static complex **phi,**chi,*rho;
static complex_dble *wrk,*cs1,*cs2;


static int alloc_arrays(int nkv)
{
   if (nkm>0)
   {
      afree(a);
      afree(b);
      afree(cs1);
   }

   a=amalloc(nkv*(nkv+1)*sizeof(*a),ALIGN);
   b=amalloc(nkv*sizeof(*b),ALIGN);
   cs1=amalloc(2*(nkv+2)*sizeof(*cs1),ALIGN);

   if ((a==NULL)||(b==NULL)||(cs1==NULL))
      return 1;

   c=a+nkv*nkv;
   cs2=cs1+nkv+2;
   nkm=nkv;

   return 0;
}


static void gcr_init(int vol,int icom,int nkv,complex **wv,complex_dble **wvd,
                     complex_dble *eta,complex_dble *psi)
{
   phi=wv;
   rho=wv[nkv];
   chi=wv+nkv+1;
   wrk=wvd[0];

   set_vd2zero(vol,psi);
   assign_vd2v(vol,eta,rho);

   rn=(double)(vnorm_square(vol,icom,rho));
   rn=sqrt(rn);
}


static void sum_vprod(int icom,int n)
{
   int i;

   if ((icom==1)&&(NPROC>1))
   {
      MPI_Reduce((double*)(cs1),(double*)(cs2),2*n,MPI_DOUBLE,
                 MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast((double*)(cs2),2*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      for (i=0;i<n;i++)
      {
         cs2[i].re=cs1[i].re;
         cs2[i].im=cs1[i].im;
      }
   }
}


static void gcr_step(int vol,int icom,int k,int nkv,
                     void (*Mop)(int k,complex *rho,complex *phi,complex *chi))
{
   int l;
   float r;
   complex z;

   (*Mop)(k,rho,phi[k],chi[k]);

   for (l=0;l<k;l++)
   {
      z=vprod(vol,0,chi[l],chi[k]);
      cs1[l].re=(double)(z.re);
      cs1[l].im=(double)(z.im);
   }

   sum_vprod(icom,k);

   for (l=0;l<k;l++)
   {
      a[nkv*l+k].re=(float)(cs2[l].re);
      a[nkv*l+k].im=(float)(cs2[l].im);
      z.re=-a[nkv*l+k].re;
      z.im=-a[nkv*l+k].im;
      mulc_vadd(vol,chi[k],chi[l],z);
   }

   r=vnorm_square(vol,0,chi[k]);
   cs1[0].re=(double)(r);
   cs1[0].im=0.0;
   z=vprod(vol,0,chi[k],rho);
   cs1[1].re=(double)(z.re);
   cs1[1].im=(double)(z.im);
   sum_vprod(icom,2);

   b[k]=(float)(sqrt(cs2[0].re));

   if (b[k]==0.0f)
      return;

   r=1.0f/b[k];
   vscale(vol,r,chi[k]);
   c[k].re=r*(float)(cs2[1].re);
   c[k].im=r*(float)(cs2[1].im);
   z.re=-c[k].re;
   z.im=-c[k].im;
   mulc_vadd(vol,rho,chi[k],z);

   rn=(double)(vnorm_square(vol,icom,rho));
   rn=sqrt(rn);
}


static void update_psi(int vol,int icom,int k,int nkv,
                       complex_dble *eta,complex_dble *psi,
                       void (*Dop)(complex_dble *v,complex_dble *w))
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

   set_v2zero(vol,rho);

   for (l=k;l>=0;l--)
      mulc_vadd(vol,rho,phi[l],c[l]);

   add_v2vd(vol,rho,psi);
   (*Dop)(psi,wrk);
   diff_vd2v(vol,eta,wrk,rho);

   rn=(double)(vnorm_square(vol,icom,rho));
   rn=sqrt(rn);
}


double fgcr4vd(int vol,int icom,
               void (*Dop)(complex_dble *v,complex_dble *w),
               void (*Mop)(int k,complex *eta,complex *psi,complex *chi),
               complex **wv,complex_dble **wvd,int nkv,int nmx,double res,
               complex_dble *eta,complex_dble *psi,int *status)
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
            (dprms[0]!=res),1,"fgcr4vd [fgcr4vd.c]",
            "Parameters are not global");

      error_root((vol<=0)||(nkv<1)||(nmx<1)||(res<=DBL_EPSILON),1,
                 "fgcr4vd [fgcr4vd.c]",
                 "Improper choice of vol,nkv,nmx or res");

      if (nkv>nkm)
      {
         ie=alloc_arrays(nkv);
         error(ie,1,"fgcr4vd [fgcr4vd.c]",
               "Unable to allocate auxiliary arrays");
      }
   }
   else
   {
      if ((vol<=0)||(nkv<1)||(nmx<1)||(res<=DBL_EPSILON))
      {
         error_loc(1,1,"fgcr4vd [fgcrvvd.c]",
                   "Improper choice of vol,nkv,nmx or res");
         (*status)=0;
         return 1.0;
      }

      if (nkv>nkm)
      {
         ie=alloc_arrays(nkv);

         if (ie)
         {
            error_loc(1,1,"fgcr4vd [fgcr4vd.c]",
                      "Unable to allocate auxiliary arrays");
            (*status)=0;
            return 1.0;
         }
      }
   }

   gcr_init(vol,icom,nkv,wv,wvd,eta,psi);
   tol=res*rn;
   (*status)=0;

   while (rn>tol)
   {
#ifdef FGCR4VD_DBG
      message("[fgcr4vd]: rn_old = %.2e\n",rn);
#endif
      rn_old=rn;

      for (k=0;;k++)
      {
         gcr_step(vol,icom,k,nkv,Mop);
         (*status)+=1;
#ifdef FGCR4VD_DBG
         message("[fgcr4vd]: k = %d, rn = %.2e\n",k,rn);
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
