
/*******************************************************************************
*
* File tmcg.c
*
* Copyright (C) 2011-2013 Martin Luescher, Stefan Schaefer
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* CG solver for the normal Wilson-Dirac equation with a twisted-mass term.
*
* The externally accessible function is
*
*   double tmcg(int nmx,double res,double mu,
*               spinor_dble *eta,spinor_dble *psi,int *status)
*     Obtains an approximate solution psi of the normal Wilson-Dirac
*     equation for given source eta (see the notes).
*
*   double tmcgeo(int nmx,double res,double mu,
*                 spinor_dble *eta,spinor_dble *psi,int *status)
*     Obtains an approximate solution psi of the normal even-odd
*     preconditioned Wilson-Dirac equation for given source eta (see
*     the notes).
*
* Notes:
*
* The normal and the normal even-odd preconditioned Wilson-Dirac equations
* are
*
*   (Dw^dag*Dw+mu^2)*psi=eta
*
*   (Dwhat^dag*Dwhat+mu^2)*psi=eta
*
* respectively.
*
* The programs are based on the standard CG algorithm (see linsolv/cgne.c).
* They assume that the improvement coefficients and the quark mass in the
* SW term have been set through set_lat_parms() and set_sw_parms() (see
* flags/lat_parms.c).
*
* All other parameters are passed through the argument list:
*
*   nmx     Maximal total number of CG iterations that may be performed.
*
*   res     Desired maximal relative residue |eta-(Dw^dag*Dw+mu^2)*psi|/|eta|
*           of the calculated solution.
*
*   mu      Value of the twisted mass in the Dirac equation.
*
*   eta     Source field. Note that source fields must respect the chosen
*           boundary conditions at time 0 and NPR0C0*L0-1, as has to be the
*           the case for physical quark fields (see doc/dirac.pdf). On exit
*           eta is unchanged unless psi=eta (which is permissible).
*
*   psi     Calculated approximate solution of the Dirac equation.
*
*   status  If the program is able to solve the Dirac equation to the
*           desired accuracy, status reports the number of CG iterations
*           that were required for the solution. Negative values indicate
*           that the program failed (-1: the algorithm did not converge,
*           -2: the solution went out of range, -3: the inversion of the
*           SW term was not safe).
*
* The programs return the norm of the residue of the calculated approximate
* solution if status[0]>=-1. Otherwise the field psi is set to zero and the
* programs return the norm of the source eta.
*
* The SW term is recalculated when needed. Evidently the solver is a global
* program that must be called on all processes simultaneously. The required
* workspaces are
*
*  spinor              5
*  spinor_dble         3
*
* (see utils/wspace.c).
*
*******************************************************************************/

#define TMCG_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "linsolv.h"
#include "forces.h"
#include "global.h"

static float mus;
static double mud;


static void Dop(spinor *s,spinor *r)
{
   Dw(mus,s,r);
   mulg5(VOLUME,r);
   mus=-mus;
}


static void Dop_dble(spinor_dble *s,spinor_dble *r)
{
   Dw_dble(mud,s,r);
   mulg5_dble(VOLUME,r);
   mud=-mud;
}


double tmcg(int nmx,double res,double mu,
            spinor_dble *eta,spinor_dble *psi,int *status)
{
   double rho,rho0,fact;
   spinor **ws;
   spinor_dble **rsd,**wsd;
   tm_parms_t tm;

   tm=tm_parms();
   if (tm.eoflg==1)
      set_tm_parms(0);

   if (query_flags(U_MATCH_UD)!=1)
      assign_ud2u();

   sw_term(NO_PTS);

   if ((query_flags(SW_UP2DATE)!=1)||
       (query_flags(SW_E_INVERTED)!=0)||(query_flags(SW_O_INVERTED)!=0))
      assign_swd2sw();

   ws=reserve_ws(5);
   wsd=reserve_wsd(2);
   rsd=reserve_wsd(1);

   mus=(float)(mu);
   mud=mu;
   rho0=sqrt(norm_square_dble(VOLUME,1,eta));
   fact=rho0/sqrt((double)(VOLUME)*(double)(24*NPROC));

   if (fact!=0.0)
   {
      assign_sd2sd(VOLUME,eta,rsd[0]);
      scale_dble(VOLUME,1.0/fact,rsd[0]);

      rho=cgne(VOLUME,1,Dop,Dop_dble,ws,wsd,nmx,res,rsd[0],psi,status);

      scale_dble(VOLUME,fact,psi);
      rho*=fact;
   }
   else
   {
      status[0]=0;
      rho=0.0;
      set_sd2zero(VOLUME,psi);
   }

   release_wsd();
   release_wsd();
   release_ws();

   if (status[0]<-1)
   {
      rho=rho0;
      set_sd2zero(VOLUME,psi);
   }

   return rho;
}


static void Doph(spinor *s,spinor *r)
{
   Dwhat(mus,s,r);
   mulg5(VOLUME/2,r);
   mus=-mus;
}


static void Doph_dble(spinor_dble *s,spinor_dble *r)
{
   Dwhat_dble(mud,s,r);
   mulg5_dble(VOLUME/2,r);
   mud=-mud;
}


double tmcgeo(int nmx,double res,double mu,
              spinor_dble *eta,spinor_dble *psi,int *status)
{
   int ifail;
   double rho,rho0,fact;
   spinor **ws;
   spinor_dble **rsd,**wsd;

   rho0=sqrt(norm_square_dble(VOLUME/2,1,eta));
   ifail=sw_term(ODD_PTS);

   if (ifail)
   {
      status[0]=-2;
      rho=rho0;
   }
   else
   {
      if (query_flags(U_MATCH_UD)!=1)
         assign_ud2u();

      if ((query_flags(SW_UP2DATE)!=1)||
          (query_flags(SW_E_INVERTED)!=0)||(query_flags(SW_O_INVERTED)!=1))
         assign_swd2sw();

      ws=reserve_ws(5);
      wsd=reserve_wsd(2);
      rsd=reserve_wsd(1);

      mus=(float)(mu);
      mud=mu;
      fact=rho0/sqrt((double)(VOLUME/2)*(double)(24*NPROC));

      if (fact!=0.0)
      {
         assign_sd2sd(VOLUME/2,eta,rsd[0]);
         scale_dble(VOLUME/2,1.0/fact,rsd[0]);

         rho=cgne(VOLUME/2,1,Doph,Doph_dble,ws,wsd,nmx,res,rsd[0],psi,status);

         scale_dble(VOLUME/2,fact,psi);
         rho*=fact;
      }
      else
      {
         status[0]=0;
         rho=0.0;
         set_sd2zero(VOLUME/2,psi);
      }

      release_wsd();
      release_wsd();
      release_ws();
   }

   if (status[0]<-1)
   {
      rho=rho0;
      set_sd2zero(VOLUME/2,psi);
   }

   return rho;
}
