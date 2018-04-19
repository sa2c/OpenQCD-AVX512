
/*******************************************************************************
*
* File sap_gcr.c
*
* Copyright (C) 2005, 2011-2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* SAP+GCR solver for the Wilson-Dirac equation.
*
* The externally accessible function is
*
*   double sap_gcr(int nkv,int nmx,double res,double mu,
*                  spinor_dble *eta,spinor_dble *psi,int *status)
*     Obtains an approximate solution psi of the Wilson-Dirac equation for
*     given source eta using the SAP-preconditioned GCR algorithm. See the
*     notes for the explanation of the parameters of the program.
*
* Notes:
*
* Depending on whether the twisted-mass flag is set or not, the program
* solves the equation
*
*   (Dw+i*mu*gamma_5*1e)*psi=eta  or  (Dw+i*mu*gamma_5)*psi=eta
*
* respectively. The twisted-mass flag is retrieved from the parameter data
* base (see flags/lat_parms.c).

* The program is based on the flexible GCR algorithm (see linsolv/fgcr.c).
* It assumes that the improvement coefficients, the quark mass in the SW
* term and the parameters of the SAP preconditioner have been set through
* set_lat_parms(), set_sw_parms() and set_sap_parms() (see the modules in
* the modules/flags directory).
*
* All other parameters are passed through the argument list:
*
*   nkv     Maximal number of Krylov vectors generated before the GCR
*           algorithm is restarted.
*
*   nmx     Maximal total number of Krylov vectors that may be generated.
*
*   res     Desired maximal relative residue |eta-D*psi|/|eta| of the
*           calculated solution.
*
*   mu      Value of the twisted mass in the Dirac equation.
*
*   eta     Source field. Note that source fields must vanish at global
*           time 0 and NPR0C0*L0-1, as has to be the case for physical
*           quark fields. eta is unchanged on exit unless psi=eta (which
*           is permissible).
*
*   psi     Calculated approximate solution of the Dirac equation. psi
*           vanishes at global time 0 and NPROC0*L0-1.
*
*   status  If the program is able to solve the Dirac equation to the
*           desired accuracy, status reports the total number of Krylov
*           vectors that were required for the solution. Negative values
*           indicate that the program failed (-1: the algorithm did not
*           converge, -2: the inversion of the SW term on the odd points
*           was not safe).
*
* The program returns the norm of the residue of the calculated approximate
* solution if status[0]>=-1. Otherwise the field psi is set to zero and the
* program returns the norm of the source eta.
*
* The SAP_BLOCKS blocks grid is automatically allocated and the SW term is
* recalculated when needed. The gauge and SW fields are then copied to the
* block grid if they are not in the proper condition.
*
* Evidently the SAP+GCR solver is a global program that must be called on
* all processes simultaneously. The required workspaces are
*
*  spinor              2*nkv+1
*  spinor_dble         2
*
* (see utils/wspace.c).
*
*******************************************************************************/

#define SAP_GCR_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "sflds.h"
#include "linalg.h"
#include "block.h"
#include "sw_term.h"
#include "dirac.h"
#include "linsolv.h"
#include "sap.h"
#include "global.h"

static double mud;
static sap_parms_t spr;


static void Dop(spinor_dble *s,spinor_dble *r)
{
   Dw_dble(mud,s,r);
}


static void Mop(int k,spinor *rho,spinor *phi,spinor *chi)
{
   int n;

   set_s2zero(VOLUME,phi);
   assign_s2s(VOLUME,rho,chi);

   for (n=0;n<spr.ncy;n++)
      sap((float)(mud),spr.isolv,spr.nmr,phi,chi);

   diff_s2s(VOLUME,rho,chi);
}


double sap_gcr(int nkv,int nmx,double res,double mu,
               spinor_dble *eta,spinor_dble *psi,int *status)
{
   int nb,isw,ifail;
   int swde,swdo,swu,swe,swo;
   double rho,rho0,fact;
   spinor **ws;
   spinor_dble **rsd,**wsd;

   spr=sap_parms();
   error_root(spr.ncy==0,1,"sap_gcr [sap_gcr.c]","SAP parameters are not set");

   blk_list(SAP_BLOCKS,&nb,&isw);

   if (nb==0)
      alloc_bgr(SAP_BLOCKS);

   if (query_grid_flags(SAP_BLOCKS,UBGR_MATCH_UD)!=1)
      assign_ud2ubgr(SAP_BLOCKS);

   if (query_flags(SWD_UP2DATE)!=1)
      sw_term(NO_PTS);

   swde=query_flags(SWD_E_INVERTED);
   swdo=query_flags(SWD_O_INVERTED);

   swu=query_grid_flags(SAP_BLOCKS,SW_UP2DATE);
   swe=query_grid_flags(SAP_BLOCKS,SW_E_INVERTED);
   swo=query_grid_flags(SAP_BLOCKS,SW_O_INVERTED);
   ifail=0;

   if (spr.isolv==0)
   {
      if ((swde==1)||(swdo==1))
         sw_term(NO_PTS);

      if ((swu!=1)||(swe==1)||(swo==1))
         assign_swd2swbgr(SAP_BLOCKS,NO_PTS);
   }
   else if (spr.isolv==1)
   {
      if ((swde!=1)&&(swdo==1))
      {
         if ((swu!=1)||(swe==1)||(swo!=1))
            assign_swd2swbgr(SAP_BLOCKS,NO_PTS);

         sw_term(NO_PTS);
      }
      else
      {
         if ((swde==1)||(swdo==1))
            sw_term(NO_PTS);

         if ((swu!=1)||(swe==1)||(swo!=1))
            ifail=assign_swd2swbgr(SAP_BLOCKS,ODD_PTS);
      }
   }
   else
      error_root(1,1,"sap_gcr [sap_gcr.c]","Unknown block solver");

   rho0=sqrt(norm_square_dble(VOLUME,1,eta));
   rho=rho0;
   status[0]=0;

   if (ifail)
      status[0]=-2;
   else
   {
      ws=reserve_ws(2*nkv+1);
      wsd=reserve_wsd(1);
      rsd=reserve_wsd(1);

      mud=mu;
      fact=rho0/sqrt((double)(VOLUME)*(double)(24*NPROC));

      if (fact!=0.0)
      {
         assign_sd2sd(VOLUME,eta,rsd[0]);
         scale_dble(VOLUME,1.0/fact,rsd[0]);

         rho=fgcr(VOLUME,1,Dop,Mop,ws,wsd,nkv,nmx,res,rsd[0],psi,status);

         scale_dble(VOLUME,fact,psi);
         rho*=fact;
      }
      else
      {
         rho=0.0;
         set_sd2zero(VOLUME,psi);
      }

      release_wsd();
      release_wsd();
      release_ws();
   }

   if (status[0]<-1)
   {
      rho=rho0;
      set_sd2zero(VOLUME,psi);
   }

   return rho;
}
