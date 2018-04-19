
/*******************************************************************************
*
* File force1.c
*
* Copyright (C) 2011-2013 Stefan Schaefer, Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Twisted mass pseudo-fermion action and force.
*
* The externally accessible functions are
*
*   double setpf1(double mu,int ipf,int icom)
*     Generates a pseudo-fermion field phi with probability proportional
*     to exp(-Spf) and returns the action Spf (see the notes).
*
*   void force1(double mu,int ipf,int isp,int icr,double c,int *status)
*     Computes the force deriving from the action Spf (see the notes).
*     The calculated force is multiplied by c and added to the molecular-
*     dynamics force field.
*
*   double action1(double mu,int ipf,int isp,int icom,int *status)
*     Returns the action Spf (see the notes).
*
* Notes:
*
* The pseudo-fermion action Spf is given by
*
*   Spf=(phi,(Dw^dag*Dw+mu^2)^(-1)*phi),
*
* where Dw denotes the (improved) Wilson-Dirac operator and phi the pseudo-
* fermion field.
*
* The common parameters of the programs in this module are:
*
*   mu            Twisted mass parameter in Spf.
*
*   ipf           Index of the pseudo-fermion field phi in the
*                 structure returned by mdflds() [mdflds.c].
*
*   isp           Index of the solver parameter set that describes
*                 the solver to be used for the solution of the
*                 Dirac equation.
*
*   icom          The action returned by the programs setpf3() and
*                 action3() is summed over all MPI processes if icom=1.
*                 Otherwise the local part of the action is returned.
*
*   status        Status values returned by the solver used for the
*                 solution of the Dirac equation.
*
* The supported solvers are CGNE, SAP_GCR and DFL_SAP_GCR. Depending
* on the program and the solver, the number of status variables varies
* and is given by:
*
*                  CGNE         SAP_GCR       DFL_SAP_GCR
*   force1()         1             2               6
*   action1()        1             1               3
*
* Note that, in force1(), the GCR solvers solve the Dirac equations twice.
* In these cases, the program writes the status values one after the other
* to the array. The bare quark mass m0 is the one last set by sw_parms()
* [flags/lat_parms.c] and it is taken for granted that the parameters of
* the solver have been set by set_solver_parms() [flags/solver_parms.c].
*
* The program force1() attempts to propagate the solutions of the Dirac
* equation along the molecular-dynamics trajectories, using the field
* stack number icr (no fields are propagated if icr=0). If this feature
* is used, the program setup_chrono() [update/chrono.c] must be called
* before force1() is called for the first time.
*
* The required workspaces of double-precision spinor fields are
*
*                  CGNE         SAP_GCR       DFL_SAP_GCR
*   setpf1()         1             1               1
*   force1()     2+(icr>0)     2+2*(icr>0)     2+2*(icr>0)
*   action1()        1             1               1
*
* (these figures do not include the workspace required by the solvers).
*
* The programs in this module perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define FORCE1_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "mdflds.h"
#include "sw_term.h"
#include "sflds.h"
#include "dirac.h"
#include "linalg.h"
#include "sap.h"
#include "dfl.h"
#include "update.h"
#include "forces.h"
#include "global.h"


double setpf1(double mu,int ipf,int icom)
{
   double act;
   spinor_dble **wsd,*phi;
   mdflds_t *mdfs;
   tm_parms_t tm;

   tm=tm_parms();
   if (tm.eoflg==1)
      set_tm_parms(0);

   wsd=reserve_wsd(1);
   random_sd(VOLUME,wsd[0],1.0);
   bnd_sd2zero(ALL_PTS,wsd[0]);
   act=norm_square_dble(VOLUME,icom,wsd[0]);

   sw_term(NO_PTS);

   mdfs=mdflds();
   phi=(*mdfs).pf[ipf];
   Dw_dble(mu,wsd[0],phi);
   mulg5_dble(VOLUME,phi);
   release_wsd();

   return act;
}


void force1(double mu,int ipf,int isp,int icr,double c,int *status)
{
   int l;
   double res0,res1;
   spinor_dble *phi,*chi,*psi,**wsd;
   spinor_dble *rho,*eta,**rsd;
   mdflds_t *mdfs;
   solver_parms_t sp;
   sap_parms_t sap;
   tm_parms_t tm;

   tm=tm_parms();
   if (tm.eoflg==1)
      set_tm_parms(0);

   mdfs=mdflds();
   sp=solver_parms(isp);
   sw_term(NO_PTS);

   wsd=reserve_wsd(2);
   phi=(*mdfs).pf[ipf];
   psi=wsd[0];
   chi=wsd[1];

   if (sp.solver==CGNE)
   {
      if (get_chrono(icr,chi))
      {
         rsd=reserve_wsd(1);
         rho=rsd[0];

         Dw_dble(-mu,chi,psi);
         mulg5_dble(VOLUME,psi);
         Dw_dble(mu,psi,rho);
         mulg5_dble(VOLUME,rho);
         mulr_spinor_add_dble(VOLUME,rho,phi,-1.0);

         res0=norm_square_dble(VOLUME,1,phi);
         res1=norm_square_dble(VOLUME,1,rho);
         res1=sqrt(res1/res0);

         if (res1<1.0)
         {
            if (res1>sp.res)
            {
               tmcg(sp.nmx,sp.res/res1,mu,rho,psi,status);
               mulr_spinor_add_dble(VOLUME,chi,psi,-1.0);
            }
            else
               status[0]=0;
         }
         else
            tmcg(sp.nmx,sp.res,mu,phi,chi,status);

         release_wsd();
      }
      else
         tmcg(sp.nmx,sp.res,mu,phi,chi,status);

      error_root(status[0]<0,1,"force1 [force1.c]",
                 "CGNE solver failed (mu = %.4e, parameter set no %d, "
                 "status = %d)",mu,isp,status[0]);
      if (icr)
         add_chrono(icr,chi);
      Dw_dble(-mu,chi,psi);
      mulg5_dble(VOLUME,psi);
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      if (get_chrono(icr,chi))
      {
         rsd=reserve_wsd(2);
         rho=rsd[0];
         eta=rsd[1];

         Dw_dble(-mu,chi,psi);
         mulg5_dble(VOLUME,psi);
         Dw_dble(mu,psi,rho);
         mulg5_dble(VOLUME,rho);
         mulr_spinor_add_dble(VOLUME,rho,phi,-1.0);

         res0=norm_square_dble(VOLUME,1,phi);
         res1=norm_square_dble(VOLUME,1,rho);
         res1=sqrt(res1/res0);

         if (res1<1.0)
         {
            if (res1>sp.res)
            {
               mulg5_dble(VOLUME,rho);
               sap_gcr(sp.nkv,sp.nmx,sp.res/res1,mu,rho,eta,status);
               mulr_spinor_add_dble(VOLUME,psi,eta,-1.0);

               res0=norm_square_dble(VOLUME,1,psi);
               res1=norm_square_dble(VOLUME,1,eta);
               res1=sqrt(res1/res0);

               if (res1<1.0)
               {
                  if (res1>sp.res)
                  {
                     mulg5_dble(VOLUME,eta);
                     sap_gcr(sp.nkv,sp.nmx,sp.res/res1,-mu,eta,rho,status+1);
                     mulr_spinor_add_dble(VOLUME,chi,rho,-1.0);
                  }
                  else
                     status[1]=0;
               }
               else
               {
                  mulg5_dble(VOLUME,psi);
                  sap_gcr(sp.nkv,sp.nmx,sp.res,-mu,psi,chi,status+1);
                  mulg5_dble(VOLUME,psi);
               }
            }
            else
            {
               status[0]=0;
               status[1]=0;
            }
         }
         else
         {
            mulg5_dble(VOLUME,phi);
            sap_gcr(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);
            mulg5_dble(VOLUME,phi);
            mulg5_dble(VOLUME,psi);
            sap_gcr(sp.nkv,sp.nmx,sp.res,-mu,psi,chi,status+1);
            mulg5_dble(VOLUME,psi);
         }

         release_wsd();
      }
      else
      {
         mulg5_dble(VOLUME,phi);
         sap_gcr(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);
         mulg5_dble(VOLUME,phi);
         mulg5_dble(VOLUME,psi);
         sap_gcr(sp.nkv,sp.nmx,sp.res,-mu,psi,chi,status+1);
         mulg5_dble(VOLUME,psi);
      }

      error_root((status[0]<0)||(status[1]<0),1,"force1 [force1.c]",
                 "SAP_GCR solver failed (mu = %.4e, parameter set no %d, "
                 "status = %d;%d)",mu,isp,status[0],status[1]);
      if (icr)
         add_chrono(icr,chi);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      if (get_chrono(icr,chi))
      {
         rsd=reserve_wsd(2);
         rho=rsd[0];
         eta=rsd[1];

         Dw_dble(-mu,chi,psi);
         mulg5_dble(VOLUME,psi);
         Dw_dble(mu,psi,rho);
         mulg5_dble(VOLUME,rho);
         mulr_spinor_add_dble(VOLUME,rho,phi,-1.0);

         res0=norm_square_dble(VOLUME,1,phi);
         res1=norm_square_dble(VOLUME,1,rho);
         res1=sqrt(res1/res0);

         if (res1<1.0)
         {
            if (res1>sp.res)
            {
               mulg5_dble(VOLUME,rho);
               dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res/res1,mu,rho,eta,status);
               mulr_spinor_add_dble(VOLUME,psi,eta,-1.0);

               res0=norm_square_dble(VOLUME,1,psi);
               res1=norm_square_dble(VOLUME,1,eta);
               res1=sqrt(res1/res0);

               if (res1<1.0)
               {
                  if (res1>sp.res)
                  {
                     mulg5_dble(VOLUME,eta);
                     dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res/res1,-mu,eta,rho,
                                  status+3);
                     mulr_spinor_add_dble(VOLUME,chi,rho,-1.0);
                  }
                  else
                  {
                     for (l=3;l<6;l++)
                        status[l]=0;
                  }
               }
               else
               {
                  mulg5_dble(VOLUME,psi);
                  dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,-mu,psi,chi,status+3);
                  mulg5_dble(VOLUME,psi);
               }
            }
            else
            {
               for (l=0;l<6;l++)
                  status[l]=0;
            }
         }
         else
         {
            mulg5_dble(VOLUME,phi);
            dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);
            mulg5_dble(VOLUME,phi);
            mulg5_dble(VOLUME,psi);
            dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,-mu,psi,chi,status+3);
            mulg5_dble(VOLUME,psi);
         }

         release_wsd();
      }
      else
      {
         mulg5_dble(VOLUME,phi);
         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);
         mulg5_dble(VOLUME,phi);
         mulg5_dble(VOLUME,psi);
         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,-mu,psi,chi,status+3);
         mulg5_dble(VOLUME,psi);
      }

      error_root((status[0]<0)||(status[1]<0)||(status[3]<0)||(status[4]<0),1,
                 "force1 [force1.c]","DFL_SAP_GCR solver failed "
                 "(mu = %.4e, parameter set no %d, status = %d,%d,%d;%d,%d,%d)",
                 mu,isp,status[0],status[1],status[2],
                 status[3],status[4],status[5]);

      if (icr)
         add_chrono(icr,chi);
   }
   else
      error_root(1,1,"force1 [force1.c]","Unknown solver");

   set_xt2zero();
   add_prod2xt(1.0,chi,psi);
   sw_frc(c);

   set_xv2zero();
   add_prod2xv(1.0,chi,psi);
   hop_frc(c);

   release_wsd();
}


double action1(double mu,int ipf,int isp,int icom,int *status)
{
   double act;
   spinor_dble *phi,*psi,**wsd,**rsd;
   mdflds_t *mdfs;
   solver_parms_t sp;
   sap_parms_t sap;
   tm_parms_t tm;

   tm=tm_parms();
   if (tm.eoflg==1)
      set_tm_parms(0);

   mdfs=mdflds();
   sp=solver_parms(isp);

   wsd=reserve_wsd(1);
   psi=wsd[0];
   phi=(*mdfs).pf[ipf];

   if (sp.solver==CGNE)
   {
      tmcg(sp.nmx,sp.res,mu,phi,psi,status);

      error_root(status[0]<0,1,"action1 [force1.c]",
                 "CGNE solver failed (mu = %.4e, parameter set no %d, "
                 "status = %d)",mu,isp,status[0]);

      rsd=reserve_wsd(1);
      Dw_dble(-mu,psi,rsd[0]);
      act=norm_square_dble(VOLUME,icom,rsd[0]);
      release_wsd();
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      mulg5_dble(VOLUME,phi);
      sap_gcr(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);
      mulg5_dble(VOLUME,phi);

      error_root(status[0]<0,1,"action1 [force1.c]",
                 "SAP_GCR solver failed (mu = %.4e, parameter set no %d, "
                 "status = %d)",mu,isp,status[0]);

      act=norm_square_dble(VOLUME,icom,psi);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      mulg5_dble(VOLUME,phi);
      dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);
      mulg5_dble(VOLUME,phi);

      error_root((status[0]<0)||(status[1]<0),1,
                 "action1 [force1.c]","DFL_SAP_GCR solver failed "
                 "(mu = %.4e, parameter set no %d, status = %d,%d,%d)",
                 mu,isp,status[0],status[1],status[2]);

      act=norm_square_dble(VOLUME,icom,psi);
   }
   else
   {
      error_root(1,1,"action1 [force1.c]","Unknown solver");
      act=0.0;
   }

   release_wsd();

   return act;
}
