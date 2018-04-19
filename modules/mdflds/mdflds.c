
/*******************************************************************************
*
* File mdflds.c
*
* Copyright (C) 2011, 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation and initialization of the MD auxiliary fields.
*
* The externally accessible functions are
*
*   mdflds_t *mdflds(void)
*     Returns the pointer to a mdflds_t structure containing the force and
*     momentum field. The fields are automatically allocated if needed.
*
*   void set_frc2zero(void)
*     Sets all force variables, including those on the boundary, to zero.
*
*   void bnd_mom2zero(void)
*     Sets the components of the momentum field on the static links
*     to zero (see the notes).
*
*   void random_mom(void)
*     Sets the elements X of the momentum field on the active links to
*     random values with distribution proportional to exp(tr{X^2}). On
*     the static links the field is set to zero (see the notes).
*
*   double momentum_action(int icom)
*     Returns the action of the momentum field. The action is summed
*     over all MPI processes if (and only if) icom=1.
*
* Notes:
*
* The arrays *.mom and *.frc in the structure returned by mflds() are the
* molecular-dynamics momentum and force fields. Their elements are ordered
* in the same way as the link variables (see main/README.global). Moreover,
* the force field includes space for 7*(BNDRY/4) additional links as do the
* gauge fields (see lattice/README.uidx).
*
* Before the momentum and force fields are allocated, the geometry arrays
* must be set. The sets of static and active links depend on the chosen
* boundary conditions. Only the field variables on the active links are
* updated in the simulations.
*
* The number npf of pseudo-fermion fields is retrieved from the parameter
* data base (see flags/hmc_parms.c). It is thus assumed that npf has been
* set when the programs in this module are called for the first time (the
* field array is otherwise set to NULL).
*
* Pseudo-fermion fields are of the same size NSPIN as other quark fields.
* In the structure returned by mdflds(), the address of the pseudo-fermion
* field with index ipf is *.pf[ipf].
*
* The programs potentially perform global operations and must be called
* simultaneously on all MPI processes.
*
*******************************************************************************/

#define MDFLDS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "sflds.h"
#include "linalg.h"
#include "mdflds.h"
#include "global.h"

static const su3_alg_dble md0={0.0};
static mdflds_t *mdfs=NULL;


static void alloc_mdflds(void)
{
   int npf,ipf;
   su3_alg_dble *mom;
   spinor_dble **pp,*p;
   hmc_parms_t hmc;

   error_root(sizeof(su3_alg_dble)!=(8*sizeof(double)),1,
              "alloc_mdflds [mdflds.c]",
              "The su3_alg_dble structures are not properly packed");

   error(iup[0][0]==0,1,"alloc_mdflds [mdflds.c]",
         "The geometry arrays are not set");

   mdfs=malloc(sizeof(*mdfs));
   mom=amalloc((8*VOLUME+7*(BNDRY/4))*sizeof(*mom),ALIGN);
   error((mdfs==NULL)||(mom==NULL),1,"alloc_mdflds [mdflds.c]",
         "Unable to allocate momentum and force fields");

   set_alg2zero(8*VOLUME+7*(BNDRY/4),mom);
   (*mdfs).mom=mom;
   (*mdfs).frc=mom+4*VOLUME;

   hmc=hmc_parms();
   npf=hmc.npf;

   if (npf>0)
   {
      pp=malloc(npf*sizeof(*pp));
      p=amalloc(npf*NSPIN*sizeof(*p),ALIGN);
      error((pp==NULL)||(p==NULL),1,"alloc_mdflds [mdflds.c]",
            "Unable to allocate pseudo-fermion fields");
      set_sd2zero(npf*NSPIN,p);

      for (ipf=0;ipf<npf;ipf++)
      {
         pp[ipf]=p;
         p+=NSPIN;
      }

      (*mdfs).npf=npf;
      (*mdfs).pf=pp;
   }
   else
   {
      (*mdfs).npf=0;
      (*mdfs).pf=NULL;
   }
}


mdflds_t *mdflds(void)
{
   if (mdfs==NULL)
      alloc_mdflds();

   return mdfs;
}


void set_frc2zero(void)
{
   if (mdfs==NULL)
      alloc_mdflds();
   else
      set_alg2zero(4*VOLUME+7*(BNDRY/4),(*mdfs).frc);
}


void bnd_mom2zero(void)
{
   int bc,ifc;
   int nlks,*lks,*lkm,npts,*pts,*ptm;
   su3_alg_dble *mom,*m;

   bc=bc_type();

   if ((bc==0)||(bc==1))
   {
      if (mdfs==NULL)
         alloc_mdflds();
      mom=(*mdfs).mom;

      if (bc==0)
      {
         lks=bnd_lks(&nlks);
         lkm=lks+nlks;

         for (;lks<lkm;lks++)
            mom[*lks]=md0;
      }
      else if (bc==1)
      {
         pts=bnd_pts(&npts);
         ptm=pts+npts;
         pts+=(npts/2);

         for (;pts<ptm;pts++)
         {
            m=mom+8*(pts[0]-(VOLUME/2));

            for (ifc=2;ifc<8;ifc++)
               m[ifc]=md0;
         }
      }
   }
}


void random_mom(void)
{
   if (mdfs==NULL)
      alloc_mdflds();
   random_alg(4*VOLUME,(*mdfs).mom);
   bnd_mom2zero();
}


double momentum_action(int icom)
{
   return 0.5*norm_square_alg(4*VOLUME,icom,(*mdfs).mom);
}
