
/*******************************************************************************
*
* File rwrat.c
*
* Copyright (C) 2012-2014 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Rational function reweighting factor.
*
* The externally accessible function is
*
*   double rwrat(int irp,int n,int *np,int *isp,double *sqn,int **status)
*     Generates a random pseudo-fermion field with normal distribution,
*     assigns its square norm to sqn and returns -ln(r) (see the notes).
*
* Notes:
*
* The computation of the reweighting factor needed to correct for the inexact
* pseudo-fermion action used for the charm and the strange quark is discussed
* in the notes "Charm and strange quark in openQCD simulations" (doc/rhmc.pdf).
*
* As explained there, an unbiased stochastic estimate r of the reweighting
* factor is obtained by choosing a pseudo-fermion field eta randomly with
* distribution proportional to exp{-(eta,eta)} and by calculating
*
*  r=exp{-(eta,[(1+Z)^(-1/2)-1]*eta)},
*
* where
*
*  Z=Dwhat^dag*Dwhat*R^2-1
*
* involves the Zolotarev rational function R of Dwhat^dag*Dwhat employed
* in the simulations and Dwhat denotes the even-odd preconditioned, O(a)-
* improved Wilson-Dirac operator.
*
* The computation requires R to be applied a number of times to a given
* spinor field. For this calculation, R is divided into n parts according
* to
*
*  R=A*{1+R_0+R_1+..+R_{n-1}},
*
*  R_k=Sum{rmu[j]/(Dwhat^dag*Dwhat+mu[j]^2),j=l[k]..l[k]+np[k]-1}
*
*  l[k+1]=l[k]+np[k], l[0]=0,
*
* mu[j] and rmu[j] being the poles and associated residues of R. The constant
* A is such that Z is of order delta, the approximation error of the Zolotarev
* rational function (see ratfcts/ratfcts.c).
*
* The arguments of the program rwrat() are:
*
*  irp      Index of the Zolotarev rational function R in the parameter
*           data base.
*
*  n        Number of parts R_k of R.
*
*  np       Array of the numbers np[k] of poles of the parts R_k. The
*           poles and zeros of R are ordered such that larger values come
*           first. R_0 includes the first np[0] poles and zeros, R_1 the
*           next np[1] poles and zeros, and so on.
*
*  isp      Array of the indices isp[k] of the solvers to be used for the
*           computation of the action of R_k on a given spinor field (the
*           supported solvers are MSCG, SAP_GCR and DFL_SAP_GCR).
*
*  sqn      Square norm of the generated random pseudo-fermion field.
*
*  status   Array of the average of the status variables returned by the
*           solvers. The array status[k] refers to the part R_k and must
*           have as many elements as are returned by the solver with index
*           isp[k] (1 for MSCG and SAP_GCR and 3 for DFL_SAP_GCR). In the
*           case of the DFL_SAP_GCR solver, status[k][2] reports the
*           number of subspace regenerations that were required.
*
* It is taken for granted that the solver parameters have been set by
* set_solver_parms() [flags/sparms.c]. The bare quark mass is taken to be
* the one last set by set_sw_parms() [flags/lat_parms.c].
*
* The computation of -ln(r) involves a series expansion of (1+Z)^(-1/2),
* which is stopped when the remainder of the series is estimated to be
* less than PRECISION_LIMIT (a macro defined below). The true accuracy
* of -ln(r) however also depends on the chosen solver residues.
*
* The program requires a workspace of double-precision spinor fields, whose
* size depends on the chosen solvers and division of R into parts R_k. For
* a given k, the required workspace is 3+np[k] (MSCG) and 5 (SAP_GCR and
* DFL_SAP_GCR). Note that these figures do not include the workspace for the
* solver programs (see forces/tmcgm.c, sap/sap_gcr.c and dfl/dfl_sap_gcr.c).
* The numbers of fields to be allocated must be greater or equal to those
* required for any k=0,..,n-1.
*
* The program in this module performs global communications and must be
* called simultaneously on all MPI processes. Some debugging information
* is printed to stdout if the macro RWRAT_DBG is defined.
*
*******************************************************************************/

#define RWRAT_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "sflds.h"
#include "sw_term.h"
#include "dirac.h"
#include "linalg.h"
#include "sap.h"
#include "dfl.h"
#include "ratfcts.h"
#include "forces.h"
#include "update.h"
#include "global.h"

#define PRECISION_LIMIT 1.0e-10

static int nps=0,ns=0,*nsps;
static double *rs;
static double cfs[4]={-0.5,0.375,-0.3125,0.2734375};


static void set_nsps(int n,int *np,int *isp)
{
   int k;

   if (n>ns)
   {
      if (ns>0)
         free(nsps);

      nsps=malloc(2*n*sizeof(*nsps));
      error(nsps==NULL,1,"set_nsps [rwrat.c]",
            "Unable to allocate auxiliary array");
      ns=n;
   }

   for (k=0;k<n;k++)
   {
      nsps[k]=np[k];
      nsps[n+k]=isp[k];
   }
}


static double set_eta(spinor_dble *eta)
{
   random_sd(VOLUME/2,eta,1.0);
   set_sd2zero(VOLUME/2,eta+(VOLUME/2));
   bnd_sd2zero(EVEN_PTS,eta);

   return norm_square_dble(VOLUME/2,1,eta);
}


static void set_res(int np,double res)
{
   int k;

   if (np>nps)
   {
      if (nps>0)
         free(rs);

      rs=malloc(np*sizeof(*rs));
      error(rs==NULL,1,"set_res [rwrat.c]",
            "Unable to allocate auxiliary array");
      nps=np;
   }

   for (k=0;k<np;k++)
      rs[k]=res;
}


static void apply_Rk(int np,int isp,double *mu,double *rmu,
                     spinor_dble *eta,spinor_dble *psi,int *status)
{
   int k,l,stat[6];
   spinor_dble **rsd;
   solver_parms_t sp;
   sap_parms_t sap;

   sp=solver_parms(isp);

   if (sp.solver==MSCG)
   {
      rsd=reserve_wsd(np);

      set_res(np,sp.res);
      tmcgm(sp.nmx,rs,np,mu,eta,rsd,status);

      error_root(status[0]<0,1,"apply_Rk [rwrat.c]","MSCG solver failed "
                 "(isp=%d, status=%d)",isp,status[0]);

      for (k=0;k<np;k++)
         mulr_spinor_add_dble(VOLUME/2,psi,rsd[k],rmu[k]);

      release_wsd();
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      rsd=reserve_wsd(2);
      mulg5_dble(VOLUME/2,eta);
      set_sd2zero(VOLUME/2,eta+(VOLUME/2));
      status[0]=0;

      for (k=0;k<np;k++)
      {
         sap_gcr(sp.nkv,sp.nmx,sp.res,mu[k],eta,rsd[0],stat);
         mulg5_dble(VOLUME/2,rsd[0]);
         set_sd2zero(VOLUME/2,rsd[0]+(VOLUME/2));
         sap_gcr(sp.nkv,sp.nmx,sp.res,-mu[k],rsd[0],rsd[1],stat+1);

         error_root((stat[0]<0)||(stat[1]<0),1,"apply_Rk [rwrat.c]",
                    "SAP_GCR solver failed (isp=%d, status=%d;%d)",
                    isp,stat[0],stat[1]);

         status[0]+=stat[0];
         status[0]+=stat[1];

         mulr_spinor_add_dble(VOLUME/2,psi,rsd[1],rmu[k]);
      }

      status[0]=(status[0]+np)/(2*np);
      mulg5_dble(VOLUME/2,eta);
      release_wsd();
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      rsd=reserve_wsd(2);
      mulg5_dble(VOLUME/2,eta);
      set_sd2zero(VOLUME/2,eta+(VOLUME/2));

      for (l=0;l<3;l++)
         status[l]=0;

      for (k=0;k<np;k++)
      {
         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu[k],eta,rsd[0],stat);
         mulg5_dble(VOLUME/2,rsd[0]);
         set_sd2zero(VOLUME/2,rsd[0]+(VOLUME/2));
         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,-mu[k],rsd[0],rsd[1],stat+3);

         error_root((stat[0]<0)||(stat[1]<0)||(stat[3]<0)||(stat[4]<0),1,
                    "apply_Rk [rwrat.c]","DFL_SAP_GCR solver failed (isp=%d, "
                    "status=%d,%d,%d;%d,%d,%d)",isp,stat[0],stat[1],
                    stat[2],stat[3],stat[4],stat[5]);

         for (l=0;l<3;l++)
         {
            status[l]+=stat[l];
            status[l]+=stat[l+3];
         }

         status[2]+=(stat[2]!=0);
         status[2]+=(stat[5]!=0);

         mulr_spinor_add_dble(VOLUME/2,psi,rsd[1],rmu[k]);
      }

      for (l=0;l<2;l++)
         status[l]=(status[l]+np)/(2*np);

      mulg5_dble(VOLUME/2,eta);
      release_wsd();
   }
   else
      error_root(1,1,"apply_Rk [rwrat.c]","Unknown or unsupported solver");
}


static void apply_QR(int n,int *np,int *isp,ratfct_t *rf,
                     spinor_dble *eta,spinor_dble *psi,int **status)
{
   int k,l,stat[3];
   double *mu,*rmu;
   spinor_dble **wsd;
   solver_parms_t sp;

   mu=(*rf).mu;
   rmu=(*rf).rmu;
   assign_sd2sd(VOLUME/2,eta,psi);

   for (k=0;k<n;k++)
   {
      apply_Rk(np[k],isp[k],mu,rmu,eta,psi,stat);
      sp=solver_parms(isp[k]);

      if (sp.solver==DFL_SAP_GCR)
      {
         for (l=0;l<2;l++)
            status[k][l]+=stat[l];

         status[k][2]+=(stat[2]!=0);
      }
      else
         status[k][0]+=stat[0];

      mu+=np[k];
      rmu+=np[k];
   }

   wsd=reserve_wsd(1);

   sw_term(ODD_PTS);
   Dwhat_dble(0.0,psi,wsd[0]);
   mulg5_dble(VOLUME/2,wsd[0]);
   assign_sd2sd(VOLUME/2,wsd[0],psi);
   scale_dble(VOLUME/2,(*rf).A,psi);

   release_wsd();
}


static void apply_Z(int n,int *np,int *isp,ratfct_t *rf,
                     spinor_dble *eta,spinor_dble *psi,int **status)
{
   spinor_dble **wsd;

   wsd=reserve_wsd(1);

   apply_QR(n,np,isp,rf,eta,wsd[0],status);
   apply_QR(n,np,isp,rf,wsd[0],psi,status);
   mulr_spinor_add_dble(VOLUME/2,psi,eta,-1.0);

   release_wsd();
}


static void init_stat(int n,int *isp,int **status)
{
   int k,l;
   solver_parms_t sp;

   for (k=0;k<n;k++)
   {
      sp=solver_parms(isp[k]);

      if (sp.solver==DFL_SAP_GCR)
      {
         for (l=0;l<3;l++)
            status[k][l]=0;
      }
      else
         status[k][0]=0;
   }
}


static void avg_stat(int nz,int n,int *isp,int **status)
{
   int k,l;
   solver_parms_t sp;

   for (k=0;k<n;k++)
   {
      sp=solver_parms(isp[k]);

      if (sp.solver==DFL_SAP_GCR)
      {
         for (l=0;l<2;l++)
            status[k][l]=(status[k][l]+nz)/(2*nz);

#ifdef RWRAT_DBG
         message("[rwrat]: status[%d] = %d,%d,%d\n",
                 k,status[k][0],status[k][1],status[k][2]);
#endif
      }
      else
      {
         status[k][0]=(status[k][0]+nz)/(2*nz);

#ifdef RWRAT_DBG
         message("[rwrat]: status[%d] = %d\n",k,status[k][0]);
#endif
      }
   }
}


double rwrat(int irp,int n,int *np,int *isp,double *sqn,int **status)
{
   int k,l,ie,irat[3],iprms[2];
   double lnr,delta,r[2];
   spinor_dble **wsd;
   ratfct_t rf;
   tm_parms_t tm;
   rat_parms_t rp;

   if (NPROC>1)
   {
      iprms[0]=irp;
      iprms[1]=n;

      MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

      error((iprms[0]!=irp)||(iprms[1]!=n),1,
            "rwrat [rwrat.c]","Parameter irp or n is not global");
   }

   rp=rat_parms(irp);
   error_root((rp.degree==0)||(n<1),1,"rwrat [rwrat.c]",
              "Undefined rational function or improper choice of n");

   if (NPROC>1)
   {
      set_nsps(n,np,isp);
      MPI_Bcast(nsps,2*n,MPI_INT,0,MPI_COMM_WORLD);
      ie=0;

      for (k=0;k<n;k++)
      {
         ie|=(nsps[k]!=np[k]);
         ie|=(nsps[n+k]!=isp[k]);
      }

      error(ie!=0,1,"rwrat [rwrat.c]","Parameters np or isp are not global");
   }

   ie=0;
   l=0;

   for (k=0;k<n;k++)
   {
      ie|=(np[k]<1);
      l+=np[k];
   }

   error_root((ie!=0)||(l!=rp.degree),1,"rwrat [rwrat.c]",
              "Improper choice of np");

   tm=tm_parms();
   if (tm.eoflg!=1)
      set_tm_parms(1);

   irat[0]=irp;
   irat[1]=0;
   irat[2]=rp.degree-1;
   rf=ratfct(irat);
   delta=rf.delta;

   init_stat(n,isp,status);
   wsd=reserve_wsd(2);
   (*sqn)=set_eta(wsd[0]);

   k=1;
   apply_Z(n,np,isp,&rf,wsd[0],wsd[1],status);
   r[0]=spinor_prod_re_dble(VOLUME/2,1,wsd[0],wsd[1]);
   r[1]=norm_square_dble(VOLUME/2,1,wsd[1]);
   lnr=cfs[0]*r[0]+cfs[1]*r[1];

#ifdef RWRAT_DBG
   message("[rwrat]: irp = %d, delta = %.1e, n = %d, precision limit = %.1e\n",
           irp,delta,n,PRECISION_LIMIT);
   message("[rwrat]: sqn = %.4e, <Z> = %.4e, <Z^2> = %.4e",
           (*sqn),r[0],r[1]);
#endif

   if ((delta*r[1])>PRECISION_LIMIT)
   {
      k=2;
      apply_Z(n,np,isp,&rf,wsd[1],wsd[0],status);
      r[0]=spinor_prod_re_dble(VOLUME/2,1,wsd[1],wsd[0]);
      r[1]=norm_square_dble(VOLUME/2,1,wsd[0]);
      lnr+=(cfs[2]*r[0]+cfs[3]*r[1]);

#ifdef RWRAT_DBG
      message(", <Z^3> = %.4e, <Z^4> = %.4e",r[0],r[1]);
#endif

      error_root((delta*r[1])>PRECISION_LIMIT,1,"rwrat [rwrat.c]",
                 "Unable to reach the required precision");
   }

#ifdef RWRAT_DBG
   message("\n");
#endif

   avg_stat(k,n,isp,status);
   release_wsd();

#ifdef RWRAT_DBG
   message("[rwrat]: -ln(r) = %.4e\n",lnr);
#endif

   return lnr;
}
