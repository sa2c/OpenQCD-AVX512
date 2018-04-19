
/*******************************************************************************
*
* File rwtm.c
*
* Copyright (C) 2012-2014 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Twisted-mass reweighting factors.
*
* The externally accessible functions are
*
*   double rwtm1(double mu1,double mu2,int isp,double *sqn,int *status)
*     Generates a random pseudo-fermion field with normal distribution,
*     assigns its square norm to sqn and returns -ln(r1) (see the notes).
*     The twisted-mass Dirac equation is solved using the solver specified
*     by the parameter set number isp.
*      The argument status must be an array of at least 1,1 and 3 elements,
*     respectively, in the case of the CGNE, SAP_GCR and DFL_SAP_GCR solver.
*     On exit the array elements contain the status values returned by the
*     solver program (when the DFL_SAP_GCR solver is used, status[2] reports
*     the number of deflation subspace regenerations that were required).
*
*   double rwtm2(double mu1,double mu2,int isp,double *sqn,int *status)
*     Generates a random pseudo-fermion field with normal distribution,
*     assigns its square norm to sqn and returns -ln(r2) (see the notes).
*     The twisted-mass Dirac equation is solved using the solver specified
*     by the parameter set number isp.
*      The argument status must be an array of at least 1,1 and 3 elements,
*     respectively, in the case of the CGNE, SAP_GCR and DFL_SAP_GCR solver.
*     On exit the array elements contain the average of the status values
*     returned by the solver program (when the DFL_SAP_GCR solver is used,
*     status[2] reports the number of deflation subspace regenerations that
*     were required).
*
* Notes:
*
* Twisted-mass reweighting of the quark determinant was introduced in
*
*  M. Luescher, F. Palombi: "Fluctuations and reweighting of the quark
*  determinant on large lattices", PoS LATTICE2008 (2008) 049.
*
* The values returned by the programs in this module are stochastic estimates
* of the factors in a product decomposition of the reweighting factors. See
* section 6 of the notes
*
*  M. Luescher: "Parameters of the openQCD main programs" [doc/parms.pdf].
*
* For a given random pseudo-fermion field eta with distribution proportional
* to exp{-(eta,eta)}, the factors r1 and r2 are defined by
*
*  r1=exp{-(eta,[R_1-1]*eta)},   r2=exp{-(eta,[R_2-1]*eta)},
*
*  R1=(X+mu2^2)*(X+mu1^2)^(-1),
*
*  R2=R1^2*(X+2*mu1^2)*(X+2*mu2^2)^(-1),  X=Dw^dag*Dw,
*
* where Dw denotes the massive O(a)-improved Wilson-Dirac operator. In
* both cases, the twisted masses must satisfy
*
*  0<=mu1<mu2.
*
* The bare quark mass is taken to be the one last set by set_sw_parms()
* [flags/lat_parms.c] and it is assumed that the chosen solver parameters
* have been set by set_solver_parms() [flags/solver_parms.c].
*
* The programs in this module perform global communications and must be
* called simultaneously on all MPI processes. They require a workspace of
* 2 double-precision spinor fields in addition to the workspace needed for
* the solver program.
*
*******************************************************************************/

#define RWTM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "sflds.h"
#include "linalg.h"
#include "sap.h"
#include "dfl.h"
#include "forces.h"
#include "update.h"
#include "global.h"


static void check_parms(double mu1,double mu2,int isp)
{
   int iprms[1];
   double dprms[2];

   if (NPROC>1)
   {
      iprms[0]=isp;
      dprms[0]=mu1;
      dprms[1]=mu2;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

      error((iprms[0]!=isp)||(dprms[0]!=mu1)||(dprms[1]!=mu2),1,
            "check_parms [rwtm.c]","Parameters are not global");
   }

   error_root((mu1<0.0)||(mu2<=mu1),1,"check_parms [rwtm.c]",
              "Twisted masses mu1,mu2 are out of range");
}


static double set_eta(spinor_dble *eta)
{
   random_sd(VOLUME,eta,1.0);
   bnd_sd2zero(ALL_PTS,eta);

   return norm_square_dble(VOLUME,1,eta);
}


double rwtm1(double mu1,double mu2,int isp,double *sqn,int *status)
{
   double lnr;
   spinor_dble *eta,*phi,**wsd;
   solver_parms_t sp;
   sap_parms_t sap;

   check_parms(mu1,mu2,isp);
   wsd=reserve_wsd(2);
   eta=wsd[0];
   phi=wsd[1];
   (*sqn)=set_eta(eta);
   sp=solver_parms(isp);

   if (sp.solver==CGNE)
   {
      tmcg(sp.nmx,sp.res,mu1,eta,phi,status);

      error_root(status[0]<0,1,"rwtm1 [rwtm.c]",
                 "CGNE solver failed (mu = %.2e, parameter set no %d, "
                 "status = %d)",mu1,isp,status[0]);

      lnr=spinor_prod_re_dble(VOLUME,1,eta,phi);
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);
      mulg5_dble(VOLUME,eta);
      sap_gcr(sp.nkv,sp.nmx,sp.res,mu1,eta,phi,status);

      error_root(status[0]<0,1,"rwtm1 [rwtm.c]",
                 "SAP_GCR solver failed (mu = %.2e, parameter set no %d, "
                 "status = %d)",mu1,isp,status[0]);

      lnr=norm_square_dble(VOLUME,1,phi);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);
      mulg5_dble(VOLUME,eta);
      dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu1,eta,phi,status);

      error_root((status[0]<0)||(status[1]<0),1,
                 "rwtm1 [rwtm.c]","DFL_SAP_GCR solver failed "
                 "(mu = %.2e, parameter set no %d, status = (%d,%d,%d))",
                 mu1,isp,status[0],status[1],status[2]);

      status[2]=(status[2]!=0);
      lnr=norm_square_dble(VOLUME,1,phi);
   }
   else
   {
      lnr=0.0;
      error_root(1,1,"rwtm1 [rwtm.c]","Unknown solver");
   }

   release_wsd();

   return (mu2*mu2-mu1*mu1)*lnr;
}


double rwtm2(double mu1,double mu2,int isp,double *sqn,int *status)
{
   int stat[3];
   double lnr1,lnr2;
   spinor_dble *eta,*phi,**wsd;
   solver_parms_t sp;
   sap_parms_t sap;

   check_parms(mu1,mu2,isp);
   wsd=reserve_wsd(2);
   eta=wsd[0];
   phi=wsd[1];
   (*sqn)=set_eta(eta);
   sp=solver_parms(isp);

   if (sp.solver==CGNE)
   {
      tmcg(sp.nmx,sp.res,mu1,eta,phi,status);

      error_root(status[0]<0,1,"rwtm2 [rwtm.c]",
                 "CGNE solver failed (mu = %.2e, parameter set no %d, "
                 "status = %d)",mu1,isp,status[0]);

      tmcg(sp.nmx,sp.res,sqrt(2.0)*mu2,eta,eta,stat);

      error_root(stat[0]<0,1,"rwtm2 [rwtm.c]",
                 "CGNE solver failed (mu = %.2e, parameter set no %d, "
                 "status = %d)",sqrt(2.0)*mu2,isp,stat[0]);
      status[0]=(status[0]+stat[0]+1)/2;

      if (mu1>0.0)
         lnr1=norm_square_dble(VOLUME,1,phi);
      else
         lnr1=0.0;

      lnr2=spinor_prod_re_dble(VOLUME,1,eta,phi);
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);
      mulg5_dble(VOLUME,eta);
      sap_gcr(sp.nkv,sp.nmx,sp.res,mu1,eta,phi,status);

      error_root(status[0]<0,1,"rwtm2 [rwtm.c]",
                 "SAP_GCR solver failed (mu = %.2e, parameter set no %d, "
                 "status = %d)",mu1,isp,status[0]);

      mulg5_dble(VOLUME,phi);
      sap_gcr(sp.nkv,sp.nmx,sp.res,sqrt(2.0)*mu2,phi,eta,stat);

      error_root(stat[0]<0,2,"rwtm2 [rwtm.c]",
                 "SAP_GCR solver failed (mu = %.2e, parameter set no %d, "
                 "status = %d)",sqrt(2.0)*mu2,isp,stat[0]);
      status[0]+=stat[0];

      if (mu1>0.0)
      {
         sap_gcr(sp.nkv,sp.nmx,sp.res,mu1,phi,phi,stat);
         error_root(stat[0]<0,3,"rwtm2 [rwtm.c]",
                    "SAP_GCR solver failed (mu = %.2e, parameter set no %d, "
                    "status = %d)",mu1,isp,stat[0]);
         status[0]=(status[0]+stat[0]+1)/3;

         lnr1=norm_square_dble(VOLUME,1,phi);
      }
      else
      {
         status[0]=(status[0]+1)/2;
         lnr1=0.0;
      }

      lnr2=norm_square_dble(VOLUME,1,eta);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);
      mulg5_dble(VOLUME,eta);
      dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu1,eta,phi,status);

      error_root((status[0]<0)||(status[1]<0),1,
                 "rwtm2 [rwtm.c]","DFL_SAP_GCR solver failed "
                 "(mu = %.2e, parameter set no %d, status = (%d,%d,%d))",
                 mu1,isp,status[0],status[1],status[2]);
      status[2]=(status[2]!=0);

      mulg5_dble(VOLUME,phi);
      dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,sqrt(2.0)*mu2,phi,eta,stat);

      error_root((stat[0]<0)||(stat[1]<0),2,
                 "rwtm2 [rwtm.c]","DFL_SAP_GCR solver failed "
                 "(mu = %.2e, parameter set no %d, status = (%d,%d,%d)",
                 sqrt(2.0)*mu2,isp,stat[0],stat[1],stat[2]);
      status[0]+=stat[0];
      status[1]+=stat[1];
      status[2]+=(stat[2]!=0);

      if (mu1>0.0)
      {
         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu1,phi,phi,stat);

         error_root((stat[0]<0)||(stat[1]<0),3,
                    "rwtm2 [rwtm.c]","DFL_SAP_GCR solver failed "
                    "(mu = %.2e, parameter set no %d, status = (%d,%d,%d)",
                    mu1,isp,stat[0],stat[1],stat[2]);

         status[0]=(status[0]+stat[0]+1)/3;
         status[1]=(status[1]+stat[1]+1)/3;
         status[2]+=(stat[2]!=0);

         lnr1=norm_square_dble(VOLUME,1,phi);
      }
      else
      {
         status[0]=(status[0]+1)/2;
         status[1]=(status[1]+1)/2;
         lnr1=0.0;
      }

      lnr2=norm_square_dble(VOLUME,1,eta);
   }
   else
   {
      lnr1=0.0;
      lnr2=0.0;
      error_root(1,1,"rwtm2 [rwtm.c]","Unknown solver");
   }

   release_wsd();

   mu1=mu1*mu1;
   mu2=mu2*mu2;

   return ((mu2-mu1)/(2.0*mu2-mu1))*(mu1*(mu2-mu1)*lnr1+2.0*mu2*mu2*lnr2);
}
