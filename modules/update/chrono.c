
/*******************************************************************************
*
* File chrono.c
*
* Copyright (C) 2007, 2011, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs needed for the propagation of solutions of the Dirac equation
* along the molecular-dynamics trajectories
*
* The externally accessible functions are
*
*   void setup_chrono(void)
*     Allocates the required memory space for the stacks of previous
*     solutions to be used in the course of the molecular-dynamics
*     trajectories. The number and size of the stacks is inferred from
*     the parameter data base.
*
*   double mdtime(void)
*     Returns the current molecular-dynamics time.
*
*   void step_mdtime(double dt)
*     Advances the molecular-dynamics time by dt.
*
*   void add_chrono(int icr,spinor_dble *psi)
*     Adds the solution psi obtained at the current molecular-dynamics 
*     time to the stack number icr of previously calculated solutions. 
*
*   int get_chrono(int icr,spinor_dble *psi)
*     Extrapolates the solutions stored in the stack number icr to the  
*     current molecular-dynamics time. The program returns 0 and leaves
*     psi unchanged if the stack does not contain any previous solutions.
*     Otherwise the program assigns the extrapolated solution to psi and
*     returns 1.
*
*   void reset_chrono(void)
*     Sets the molecular-dynamics time and all counters of previously
*     computed solutions to zero.
*
* Notes:
*
* The propagation of the solutions of the Dirac equation was proposed by
*
*   R.C. Brower et al., "Chronological inversion method for the Dirac
*   matrix in Hybrid Monte Carlo", Nucl. Phys. B484 (1997) 353
*
* Here the solutions are propagated using a polynomial extrapolation. The
* maximal number of solutions to be kept in memory can be chosen for each
* solution stack separately.
*
* Each quark force specified in the parameter data base may have up to 4
* solution stacks associated with it (see flags/force_parms.c). In all
* cases, the chronological propagation of the solutions can be turned off
* by setting the maximal numbers of fields to be kept in memory to zero.
* Internally the stacks are labeled by an index icr>=0, where the empty
* stack has index icr=0 and all other stacks have index icr>0. The stack
* indices are included in the force parameter sets.
*
* The module includes a clock that serves to keep track of the molecular-
* dynamics times at which the Dirac equation is solved. The clock is
* advanced by the molecular-dynamics integrator (see update/mdint.c).
*
* All programs in this module should be called simultaneously on all MPI
* processes.
*
*******************************************************************************/

#define CHRONO_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "sflds.h"
#include "linalg.h"
#include "update.h"
#include "global.h"

typedef struct
{
   int ncr;
   int isd,nsd;
   double *ta;
   spinor_dble **sd;
} stack_t;

static int nst=0;
static double mdt=0.0;
static stack_t *st=NULL;


static void init_stacks(void)
{
   int ncr,icr,k;
   double *ta;
   spinor_dble **sd;

   for (icr=0;icr<nst;icr++)
   {
      st[icr].isd=0;
      st[icr].nsd=0;
      
      ncr=st[icr].ncr;
      ta=st[icr].ta;
      sd=st[icr].sd;

      for (k=0;k<ncr;k++)
      {
         ta[k]=0.0;
         set_sd2zero(VOLUME,sd[k]);
      }
   }
}


static void free_stacks(void)
{
   if (nst>0)
   {
      free(st[1].ta);
      afree(st[1].sd[0]);
      free(st[1].sd);
      free(st);
      
      nst=0;
      st=NULL;
   }
}


static void alloc_stacks(void)
{
   int i,j,k;
   hmc_parms_t hmc;
   mdint_parms_t mdp;
   force_parms_t fp;

   hmc=hmc_parms();
   
   for (i=0;i<hmc.nlv;i++)
   {
      mdp=mdint_parms(i);

      for (j=0;j<mdp.nfr;j++)
      {
         fp=force_parms(mdp.ifr[j]);

         for (k=0;k<4;k++)
         {
            if (fp.icr[k]>nst)
               nst=fp.icr[k];
         }
      }
   }  

   if (nst>0)
   {
      nst+=1;
      st=malloc(nst*sizeof(*st));
      error(st==NULL,1,"alloc_stacks [chrono.c]",
            "Unable to allocate stack structures");

      for (i=0;i<nst;i++)
      {
         st[i].ncr=0;
         st[i].ta=NULL;
         st[i].sd=NULL;
      }
      
      for (i=0;i<hmc.nlv;i++)
      {
         mdp=mdint_parms(i);

         for (j=0;j<mdp.nfr;j++)
         {
            fp=force_parms(mdp.ifr[j]);

            for (k=0;k<4;k++)
               st[fp.icr[k]].ncr=fp.ncr[k];
         }
      }
   }
}


void setup_chrono(void)
{
   int ncr,icr,k;
   double *ta;
   spinor_dble **sd,*s;
   
   free_stacks();
   alloc_stacks();

   if (nst>0)
   {
      ncr=0;

      for (icr=0;icr<nst;icr++)
         ncr+=st[icr].ncr;

      if (ncr>0)
      {
         ta=malloc(ncr*sizeof(*ta));
         sd=malloc(ncr*sizeof(*sd));
         s=amalloc(ncr*VOLUME*sizeof(*s),ALIGN);

         error((ta==NULL)||(sd==NULL)||(s==NULL),1,"alloc_stacks [chrono.c]",
               "Unable to allocate field stacks");
      }
      else
      {
         ta=NULL;
         sd=NULL;
         s=NULL;
      }

      for (icr=1;icr<nst;icr++)
      {
         ncr=st[icr].ncr;

         if (ncr>0)
         {
            st[icr].ta=ta;
            st[icr].sd=sd;
            
            for (k=0;k<ncr;k++)
            {
               sd[k]=s;
               s+=VOLUME;
            }

            ta+=ncr;
            sd+=ncr;
         }
      }
      
      init_stacks();
   }

   mdt=0.0;
}


double mdtime(void)
{
   return mdt;
}


void step_mdtime(double dt)
{
   mdt+=dt;
}


void add_chrono(int icr,spinor_dble *psi)
{
   int ncr,isd,jsd,nsd;
   
   if ((icr>0)&&(icr<nst))
   {
      ncr=st[icr].ncr;
      isd=st[icr].isd;
      nsd=st[icr].nsd;

      if (nsd==ncr)
      {
         st[icr].ta[isd]=mdt;
         assign_sd2sd(VOLUME,psi,st[icr].sd[isd]);

         isd+=1;
         if (isd==ncr)
            isd=0;
         st[icr].isd=isd;
      }
      else
      {         
         jsd=isd+nsd;
         if (jsd>=ncr)
            jsd-=ncr;
         
         st[icr].ta[jsd]=mdt;
         assign_sd2sd(VOLUME,psi,st[icr].sd[jsd]);
         st[icr].nsd+=1;
      }
   }
   else
      error_loc(1,1,"add_chrono [chrono.c]","Unknown field stack");
}


int get_chrono(int icr,spinor_dble *psi)
{
   int ncr,nsd,isd;
   int k,l,ksd,lsd;
   double *ta,c;
   spinor_dble **sd;

   if ((icr>0)&&(icr<nst))
   {
      nsd=st[icr].nsd;

      if (nsd==0)
         return 0;

      ncr=st[icr].ncr;      
      isd=st[icr].isd;
      ta=st[icr].ta;
      sd=st[icr].sd;

      set_sd2zero(VOLUME,psi);

      for (k=0;k<nsd;k++)
      {
         ksd=isd+k;
         if (ksd>=ncr)
            ksd-=ncr;
         c=1.0;

         for (l=0;l<nsd;l++)
         {
            if (l!=k)
            {
               lsd=isd+l;
               if (lsd>=ncr)
                  lsd-=ncr;

               c*=((mdt-ta[lsd])/(ta[ksd]-ta[lsd]));
            }
         }

         mulr_spinor_add_dble(VOLUME,psi,sd[ksd],c);
      }

      return 1;
   }
   else if (icr==0)
      return 0;
   else
   {
      error_loc(1,1,"get_chrono [chrono.c]","Unknown field stack");
      return 0;
   }
}


void reset_chrono(void)
{
   int icr;

   for (icr=0;icr<nst;icr++)
   {
      st[icr].isd=0;
      st[icr].nsd=0;
   }

   mdt=0.0;
}
