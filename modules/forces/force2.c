
/*******************************************************************************
*
* File force2.c
*
* Copyright (C) 2011-2013 Stefan Schaefer, Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Hasenbusch twisted_mass pseudo-fermion action and force.
*
* The externally accessible functions are
*
*   double setpf2(double mu0,double mu1,int ipf,int isp,int icom,
*                 int *status)
*     Generates a pseudo-fermion field phi with probability proportional
*     to exp(-Spf) and returns the action Spf-(phi,phi) (see the notes).
*
*   void force2(double mu0,int mu1,int ipf,int isp,int icr,double c,
*               int *status)
*     Computes the force deriving from the action Spf (see the notes).
*     The calculated force is multiplied by c and added to the molecular-
*     dynamics force field.
*
*   double action2(double mu0,double mu1,int ipf,int isp,int icom,
*                  int *status)
*     Returns the action Spf-(phi,phi) (see the notes).
*
* Notes:
*
* The pseudo-fermion action Spf is given by
*
*   Spf=(phi,(Dw^dag*Dw+mu1^2)(Dw^dag*Dw+mu0^2)^(-1)*phi)
*
*      =(phi,phi)+(mu1^2-mu0^2)*(phi,(Dw^dag*Dw+mu0^2)^(-1)*phi)
*
* where Dw denotes the (improved) Wilson-Dirac operator and phi the pseudo-
* fermion field.
*
* The common parameters of the programs in this module are:
*
*   mu0,mu1       Twisted mass parameters in Spf.
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
*   setpf2()         1             1               3
*   force2()         1             2               6
*   action2()        1             1               3
*
* The solver used in the case of setpf2() is for the Dirac equation with
* twisted mass mu1, while force2() and action2() use the solver for the
* equation with twisted mass mu0. Different solvers may be needed in the
* two cases if mu1>>mu0, for example.
*
* Note that, in force2(), the GCR solvers solve the Dirac equations twice.
* In these cases, the program writes the status values one after the other
* to the array. The bare quark mass m0 is the one last set by sw_parms()
* [flags/lat_parms.c] and it is taken for granted that the parameters of
* the solver have been set by set_solver_parms() [flags/solver_parms.c].
*
* The program force2() attempts to propagate the solutions of the Dirac
* equation along the molecular-dynamics trajectories, using the field
* stack number icr (no fields are propagated if icr=0). If this feature
* is used, the program setup_chrono() [update/chrono.c] must be called
* before force2() is called for the first time.
*
* The required workspaces of double-precision spinor fields are
*
*                  CGNE         SAP_GCR       DFL_SAP_GCR
*   setpf2()         1             1               1
*   force2()     2+(icr>0)     2+2*(icr>0)     2+2*(icr>0)
*   action2()        1             1               1
*
* (these figures do not include the workspace required by the solvers).
*
* The programs in this module perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define FORCE2_C

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
#include "forces.h"
#include "global.h"


double setpf2(double mu0,double mu1,int ipf,int isp,int icom,int *status)
{
   double act;
   complex_dble z;
   spinor_dble **wsd,**rsd;
   spinor_dble *phi,*psi,*chi;
   mdflds_t *mdfs;
   solver_parms_t sp;
   sap_parms_t sap;
   tm_parms_t tm;

   tm=tm_parms();
   if (tm.eoflg==1)
      set_tm_parms(0);

   mdfs=mdflds();
   phi=(*mdfs).pf[ipf];
   wsd=reserve_wsd(1);
   psi=wsd[0];

   random_sd(VOLUME,phi,1.0);
   bnd_sd2zero(ALL_PTS,phi);
   sp=solver_parms(isp);

   if (sp.solver==CGNE)
   {
      tmcg(sp.nmx,sp.res,mu1,phi,psi,status);

      error_root(status[0]<0,1,"setpf2 [force2.c]","CGNE solver failed "
                 "(mu = %.4e, parameter set no %d, status = %d)",
                 mu1,isp,status[0]);

      rsd=reserve_wsd(1);
      chi=rsd[0];
      assign_sd2sd(VOLUME,psi,chi);
      Dw_dble(-mu1,chi,psi);
      mulg5_dble(VOLUME,psi);
      release_wsd();
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      mulg5_dble(VOLUME,phi);
      sap_gcr(sp.nkv,sp.nmx,sp.res,mu1,phi,psi,status);
      mulg5_dble(VOLUME,phi);

      error_root(status[0]<0,1,"setpf2 [force2.c]","SAP_GCR solver failed "
                 "(mu = %.4e, parameter set no %d, status = %d)",
                 mu1,isp,status[0]);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      mulg5_dble(VOLUME,phi);
      dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu1,phi,psi,status);
      mulg5_dble(VOLUME,phi);

      error_root((status[0]<0)||(status[1]<0),1,
                 "setpf2 [force2.c]","DFL_SAP_GCR solver failed "
                 "(mu = %.4e, parameter set no %d, status = %d,%d,%d)",
                 mu1,isp,status[0],status[1],status[2]);
   }
   else
      error_root(1,1,"setpf2 [force2.c]","Unknown solver");

   z.re=0.0;
   z.im=mu0-mu1;
   mulc_spinor_add_dble(VOLUME,phi,psi,z);
   act=(mu1*mu1-mu0*mu0)*norm_square_dble(VOLUME,icom,psi);
   release_wsd();

   return act;
}


void force2(double mu0,double mu1,int ipf,int isp,int icr,
            double c,int *status)
{
   double dmu2;

   dmu2=mu1*mu1-mu0*mu0;

   force1(mu0,ipf,isp,icr,dmu2*c,status);
}


double action2(double mu0,double mu1,int ipf,int isp,int icom,int *status)
{
   double dmu2,act;

   dmu2=mu1*mu1-mu0*mu0;
   act=dmu2*action1(mu0,ipf,isp,icom,status);

   return act;
}
