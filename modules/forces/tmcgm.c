
/*******************************************************************************
*
* File tmcgm.c
*
* Copyright (C) 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Multi-shift CG solver for the normal even-odd preconditioned Wilson-Dirac
* equation (Dwhat^dag*Dwhat+mu^2)*psi=eta with a twisted-mass term.
*
* The externally accessible function is
*
*   void tmcgm(int nmx,double *res,int nmu,double *mu,
*              spinor_dble *eta,spinor_dble **psi,int *status)
*     Obtains approximate solutions psi[0],..,psi[nmu-1] of the normal
*     even-odd preconditioned Wilson-Dirac equation for given source eta
*     and nmu values of the twisted-mass parameter mu. See the notes for
*     the explanation of the parameters of the program.
*
* Notes:
*
* The program is based on the multi-shift CG algorithm (see linsolv/mscg.c).
* It assumes that the improvement coefficients and the quark mass in the
* SW term have been set through set_lat_parms() and set_sw_parms() (see
* flags/lat_parms.c).
*
* All other parameters are passed through the argument list:
*
*   nmx     Maximal total number of CG iterations that may be performed.
*
*   res     Array of the desired maximal relative residues of the
*           calculated solutions (nmu elements)
*
*   nmu     Number of twisted masses mu.
*
*   mu      Array of the twisted masses (nmu elements)
*
*   eta     Source field. Note that source fields must respect the chosen
*           boundary conditions at time 0 and NPR0C0*L0-1, as has to be the
*           the case for physical quark fields (see doc/dirac.pdf).
*
*   psi     Array of the calculated approximate solutions of the Dirac
*           equations (Dwhat^dag*Dwhat+mu^2)*psi=eta (nmu elements).
*
*   status  If the program was able to solve the Dirac equations to the
*           desired accuracy, status[0] reports the total number of CG
*           iterations that were required. Negative values indicate that
*           the program failed (-1: the algorithm did not converge, -2:
*           the inversion of the SW term on the odd points was not safe).
*
* The source field eta must be different from psi[0],..,psi[nmu-1]. If
* status[0]>=-1 the calculated approximate solutions are returned. In
* all other cases, the fields are set to zero.
*
* The SW term is recalculated when needed. Evidently the solver is a global
* program that must be called on all processes simultaneously. The required
* workspace is
*
*  spinor_dble         3+nmu (5 if nmu=1)
*
* (see utils/wspace.c).
*
*******************************************************************************/

#define TMCGM_C

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

static int iop=0;


static void Dop_dble(double mu,spinor_dble *s,spinor_dble *r)
{
   if (iop==0)
      Dwhat_dble(mu,s,r);
   else
      Dwhat_dble(-mu,s,r);

   mulg5_dble(VOLUME/2,r);
   iop^=0x1;
}


void tmcgm(int nmx,double *res,int nmu,double *mu,
           spinor_dble *eta,spinor_dble **psi,int *status)
{
   int ifail,k;
   spinor_dble **wsd;

   ifail=sw_term(ODD_PTS);

   if (ifail)
   {
      status[0]=-2;

      for (k=0;k<nmu;k++)
         set_sd2zero(VOLUME/2,psi[k]);
   }
   else
   {
      if (nmu==1)
         wsd=reserve_wsd(5);
      else
         wsd=reserve_wsd(3+nmu);

      mscg(VOLUME/2,1,nmu,mu,Dop_dble,wsd,nmx,res,eta,psi,status);
      release_wsd();
   }
}
