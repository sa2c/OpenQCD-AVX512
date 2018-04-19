
/*******************************************************************************
*
* File mscg.c
*
* Copyright (C) 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic multi-shift CG solver program for the lattice Dirac equation
*
* The externally accessible function is
*
*   void mscg(int vol,int icom,int nmu,double *mu,
*             void (*Dop_dble)(double mu,spinor_dble *s,spinor_dble *r),
*             spinor_dble **wsd,int nmx,double *res,
*             spinor_dble *eta,spinor_dble **psi,int *status)
*     Solution of the Dirac equation (D^dag*D+mu^2)*psi=eta for a given
*     source eta and one or more values of mu using the multi-shift CG
*     algorithm. See the notes for the explanation of the parameters of
*     the program.
*
* Notes:
*
* The algorithm implemented in this module is described in the notes
* "Multi-shift conjugate gradient algorithm" (file doc/mscg.pdf).
*  
* The program Dop_dble() for the Dirac operator is assumed to have the 
* following properties:
*
*   void Dop_dble(double mu,spinor_dble *s,spinor_dble *r)
*     Application of an operator Op or its hermitian conjugate Op^dag
*     to the double-precision Dirac field s and assignment of the result
*     to r (where r is different from s). The operator must be such that
*     the identity Op^dag*Op=D^dag*D+mu^2 holds. Op and Op^dag are applied
*     alternatingly, i.e. the first call of the program applies Op, the
*     next call Op^dag, then Op again and so on. In all cases, the source
*     field s remains unchanged.
*
* The other parameters of the program mscg() are:
*
*   vol     Number of spinors in the Dirac fields.         
*
*   icom    Indicates whether the equation to be solved is a local
*           equation (icom=0) or a global one (icom=1). Scalar products
*           are summed over all MPI processes if icom=1, while no
*           communications are performed if icom=0.
*
*   nmu     Number of shifts mu.
*
*   mu      Array of the shifts mu (nmu elements).
*
*   nmx     Maximal number of CG iterations that may be applied.
*
*   res     Array of the desired maximal relative residues of the
*           calculated solutions (nmu elements). 
*
*   wsd     Array of at least 3+nmu (5 if nmu=1) double-precision spinor
*           fields (used as work space).
*
*   eta     Source field (unchanged on exit).
*
*   psi     Array of the calculated approximate solutions of the Dirac 
*           equations (D^dag*D+mu^2)*psi=eta (nmu elements).
*
*   status  On exit, this parameter reports the number of CG iterations
*           that were required, or a negative value if the program failed.
*
* The spinor fields must have at least vol elements and must be such that
* the program Dop_dble() acts correctly on them. Some debugging output is
* printed to stdout on process 0 if the macro MSCG_DBG is defined.
*
*******************************************************************************/

#define MSCG_C

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

typedef struct
{
   int stop;
   double s,tol;
   double ah,bh,gh,rh;
   spinor_dble *xh,*ph;
} cgsh_t;

typedef struct
{
   int k,stop;
   double mu,tol;
   double a,b;
   double rn0,rn,rnsq;
   spinor_dble *x,*r,*p,*ap,*w;
} cgs_t;

static int ns=0;
static double *dprms;
static cgs_t cgs;
static cgsh_t *cgsh;


static int alloc_cgs(int nmu,double *mu,double *res,spinor_dble **wsd,
                     spinor_dble **psi)
{
   int k,l,k0;
   
   if (nmu>ns)
   {
      if (ns>0)
         free(dprms);
      if (ns>1)
         free(cgsh);
      
      dprms=malloc(2*nmu*sizeof(*dprms));
      if (dprms==NULL)
         return 1;

      if (nmu>1)
      {
         cgsh=malloc((nmu-1)*sizeof(*cgsh));
         if (cgsh==NULL)
            return 1;
      }
      
      ns=nmu;
   }

   k0=0;

   for (k=1;k<nmu;k++)
   {
      if (fabs(mu[k])<fabs(mu[k0]))
         k0=k;
   }

   cgs.k=k0;
   cgs.stop=0;
   cgs.mu=mu[k0];
   cgs.a=1.0;
   cgs.b=0.0;
   cgs.tol=res[k0];
   
   cgs.x=psi[k0];
   cgs.r=wsd[0];
   cgs.p=wsd[1];
   cgs.ap=wsd[2];
   cgs.w=wsd[3];

   l=0;

   for (k=0;k<nmu;k++)
   {
      if (k!=k0)
      {
         cgsh[l].stop=0;
         cgsh[l].s=mu[k]*mu[k]-mu[k0]*mu[k0];
         cgsh[l].tol=res[k];
         cgsh[l].ah=1.0;
         cgsh[l].bh=0.0;
         cgsh[l].gh=1.0;
         cgsh[l].rh=1.0;

         cgsh[l].xh=psi[k];
         cgsh[l].ph=wsd[4+l];

         l+=1;
      }
   }

   return 0;
}


static void cg_init(int vol,int icom,int nmu,spinor_dble *eta)
{
   int k;
   
   set_sd2zero(vol,cgs.x);
   assign_sd2sd(vol,eta,cgs.r);
   assign_sd2sd(vol,eta,cgs.p);

   cgs.rnsq=norm_square_dble(vol,icom,eta);
   cgs.rn=sqrt(cgs.rnsq);
   cgs.rn0=cgs.rn;
   cgs.tol*=cgs.rn0;

   for (k=0;k<(nmu-1);k++)
   {
      set_sd2zero(vol,cgsh[k].xh);
      assign_sd2sd(vol,eta,cgsh[k].ph);

      cgsh[k].tol*=cgs.rn0;
   }
}


static void cg_step1(int vol,int icom,int nmu,
                     void (*Dop_dble)(double mu,spinor_dble *s,spinor_dble *r))
{
   int k;
   double om;
   
   Dop_dble(cgs.mu,cgs.p,cgs.w);
   Dop_dble(cgs.mu,cgs.w,cgs.ap);

   om=cgs.b/cgs.a;
   cgs.a=cgs.rnsq/norm_square_dble(vol,icom,cgs.w);
   om*=cgs.a;

   for (k=0;k<(nmu-1);k++)
   {
      if (cgsh[k].stop==0)
      {
         cgsh[k].rh=1.0/(1.0+cgsh[k].s*cgs.a+(1.0-cgsh[k].rh)*om);
         cgsh[k].ah=cgsh[k].rh*cgs.a;
      }
   }
}


static void cg_step2(int vol,int nmu)
{
   int k;
   
   mulr_spinor_add_dble(vol,cgs.x,cgs.p,cgs.a);
   mulr_spinor_add_dble(vol,cgs.r,cgs.ap,-cgs.a);

   for (k=0;k<(nmu-1);k++)
   {
      if (cgsh[k].stop==0)
         mulr_spinor_add_dble(vol,cgsh[k].xh,cgsh[k].ph,cgsh[k].ah);
   }
}


static void cg_step3(int vol,int icom,int nmu)
{
   int k;
   double rnsq,rh;

   rnsq=norm_square_dble(vol,icom,cgs.r);
   cgs.b=rnsq/cgs.rnsq;
   cgs.rnsq=rnsq;
   cgs.rn=sqrt(rnsq);

   for (k=0;k<(nmu-1);k++)
   {
      if (cgsh[k].stop==0)
      {
         rh=cgsh[k].rh;
         cgsh[k].bh=rh*rh*cgs.b;
         cgsh[k].gh*=rh;
      }
   }
}


static void cg_step4(int vol,int nmu)
{
   int k;

   combine_spinor_dble(vol,cgs.p,cgs.r,cgs.b,1.0);

   for (k=0;k<(nmu-1);k++)
   {
      if (cgsh[k].stop==0)
         combine_spinor_dble(vol,cgsh[k].ph,cgs.r,cgsh[k].bh,cgsh[k].gh);
   }
}


static int set_stop_flag(int nmu)
{
   int nstop,k;
   double rn;

   rn=cgs.rn;
   cgs.stop|=(rn<(0.99*cgs.tol));
   nstop=cgs.stop;

   for (k=0;k<(nmu-1);k++)
   {
      cgsh[k].stop|=((cgsh[k].gh*rn)<(0.99*cgsh[k].tol));
      nstop+=cgsh[k].stop;
   }
   
   return nstop;
}


static int check_res(int vol,int icom,double mu,
                     void (*Dop_dble)(double mu,spinor_dble *s,spinor_dble *r),
                     spinor_dble **wsd,int nmx,double res,
                     spinor_dble *eta,spinor_dble *psi,int *ncg)
{
   double tol,rnsq,rn,a,b;
   spinor_dble *x,*r,*p,*ap,*w;

   tol=res*cgs.rn0;
   r=wsd[0];
   w=wsd[1];
   
   Dop_dble(mu,psi,w);
   Dop_dble(mu,w,r);

   mulr_spinor_add_dble(vol,r,eta,-1.0);
   rnsq=norm_square_dble(vol,icom,r);
   rn=sqrt(rnsq);

   if (rn<=tol)
      return 0;
   else if ((*ncg)>=nmx)
      return 1;

   x=wsd[2];
   p=wsd[3];
   ap=wsd[4];

   set_sd2zero(vol,x);
   assign_sd2sd(vol,r,p);
   
   while ((rn>tol)&&((*ncg)<nmx))
   {
      Dop_dble(mu,p,w);
      Dop_dble(mu,w,ap);

      a=rnsq/norm_square_dble(vol,icom,w);
      mulr_spinor_add_dble(vol,x,p,a);
      mulr_spinor_add_dble(vol,r,ap,-a);

      rn=norm_square_dble(vol,icom,r);
      b=rn/rnsq;
      rnsq=rn;
      rn=sqrt(rnsq);
         
      combine_spinor_dble(vol,p,r,b,1.0);
      (*ncg)+=1;
   }

   mulr_spinor_add_dble(vol,psi,x,-1.0);

   if (rn<=tol)
      return 0;
   else
      return 1;
}


void mscg(int vol,int icom,int nmu,double *mu,
          void (*Dop_dble)(double mu,spinor_dble *s,spinor_dble *r),
          spinor_dble **wsd,int nmx,double *res,
          spinor_dble *eta,spinor_dble **psi,int *status)
{
   int ncg,nstop,k,ie;
   int iprms[3];

   if ((icom==1)&&(NPROC>1))
   {
      iprms[0]=vol;
      iprms[1]=nmu;
      iprms[2]=nmx;

      MPI_Bcast(iprms,3,MPI_INT,0,MPI_COMM_WORLD);
      error((iprms[0]!=vol)||(iprms[1]!=nmu)||(iprms[2]!=nmx),1,
            "mscg [mscg.c]","Integer parameters are not global");
      error_root((vol<1)||(nmu<1)||(nmx<1),1,"mscg [mscg.c]",
                 "Improper choice of vol,nmu or nmx");
      
      ie=alloc_cgs(nmu,mu,res,wsd,psi);
      error(ie!=0,1,"mscg [mscg.c]","Unable to allocate auxiliary arrays");

      for (k=0;k<nmu;k++)
      {
         dprms[k]=mu[k];
         dprms[nmu+k]=res[k];
      }

      MPI_Bcast(dprms,2*nmu,MPI_DOUBLE,0,MPI_COMM_WORLD);

      for (k=0;k<nmu;k++)
      {
         ie|=(dprms[k]!=mu[k]);
         ie|=(dprms[nmu+k]!=res[k]);
      }      

      error(ie!=0,1,"mscg [mscg.c]","Shifts or residues are not global");

      for (k=0;k<nmu;k++)
         ie|=(res[k]<=DBL_EPSILON);

      error_root(ie!=0,1,"mscg [mscg.c]","Improper choice of residues");
   }
   else
   {
      if ((vol<1)||(nmu<1)||(nmx<1))
      {
         error_loc(1,1,"mscg [mscg.c]",
                   "Improper choice of vol,nmu or nmx");
         (*status)=0;
         return;
      }

      ie=alloc_cgs(nmu,mu,res,wsd,psi);

      if (ie!=0)
      {
         error_loc(1,1,"mscg [mscg.c]",
                   "Unable to allocate auxiliary arrays");
         (*status)=0;
         return;
      }

      for (k=0;k<nmu;k++)
         ie|=(res[k]<=DBL_EPSILON);
      
      if (ie!=0)
      {
         error_loc(1,1,"mscg [mscg.c]",
                   "Improper choice of residues");
         (*status)=0;
         return;
      }      
   }

#ifdef MSCG_DBG
   message("[mscg]: nmu = %d, mu = %.2e",nmu,cgs.mu);

   for (k=0;k<nmu;k++)
   {
      if (k!=cgs.k)
         message(", %.2e",mu[k]);
   }
#endif

   cg_init(vol,icom,nmu,eta);
   nstop=set_stop_flag(nmu);
   ncg=0;
   
#ifdef MSCG_DBG   
   message("\n[mscg]: tol = %.2e",cgs.tol);

   for (k=0;k<(nmu-1);k++)
      message(", %.2e",cgsh[k].tol);

   message("\n");
#endif
   
   while ((ncg<nmx)&&(nstop<nmu))
   {
#ifdef MSCG_DBG
      if (cgs.stop)
         message("[mscg]: res =         ");
      else
         message("[mscg]: res = %.2e",cgs.rn);

      for (k=0;k<(nmu-1);k++)
      {
         if (cgsh[k].stop)
            message(",         ");
         else
            message(", %.2e",cgsh[k].gh*cgs.rn);
      }

      message("\n");
#endif

      cg_step1(vol,icom,nmu,Dop_dble);
      cg_step2(vol,nmu);
      cg_step3(vol,icom,nmu);
      cg_step4(vol,nmu);
      
      nstop=set_stop_flag(nmu);
      ncg+=1;
   }

#ifdef MSCG_DBG
   message("[mscg]: ncg = %d, nstop = %d\n",ncg,nstop);
#endif
   
   if ((ncg==nmx)&&(nstop<nmu))
   {
      (*status)=-1;
      return;
   }

   k=cgs.k;
   ie=check_res(vol,icom,mu[k],Dop_dble,wsd,nmx,res[k],eta,psi[k],&ncg);

#ifdef MSCG_DBG
   message("[mscg]: ie,ncg = %d,%d",ie,ncg);
#endif

   for (k=0;k<nmu;k++)
   {
      if (k!=cgs.k)
      {
         ie+=check_res(vol,icom,mu[k],Dop_dble,wsd,nmx,res[k],eta,psi[k],&ncg);

#ifdef MSCG_DBG
         message("; %d,%d",ie,ncg);
#endif
      }
   }

#ifdef MSCG_DBG   
   message("\n");
#endif   
   
   if (ie!=0)
      (*status)=-1;
   else
      (*status)=ncg;

   return;
}
