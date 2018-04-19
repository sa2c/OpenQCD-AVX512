
/*******************************************************************************
*
* File force3.c
*
* Copyright (C) 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Rational function forces.
*
* The externally accessible functions are
*
*   double setpf3(int *irat,int ipf,int isw,int isp,int icom,int *status)
*     Generates a pseudo-fermion field phi with probability proportional
*     to exp(-Spf) and returns the action Spf+Sdet-(phi,phi) if isw=1 or
*     Spf-(phi,phi) if isw!=1 (see the notes).
*
*   void force3(int *irat,int ipf,int isw,int isp,double c,int *status)
*     Computes the force deriving from the action Spf+Sdet if isw=1 or
*     Spf if isw!=1 (see the notes). The calculated force is multiplied
*     by c and added to the molecular-dynamics force field.
*
*   double action3(int *irat,int ipf,int isw,int isp,int icom,int *status)
*     Returns the action Spf+Sdet-(phi,phi) if isw=1 or Spf-(phi,phi) if
*     isw!=1 (see the notes).
*
* Notes:
*
* Simulations including the charm and/or the strange quark are based on
* a version of the RHMC algorithm. See the notes "Charm and strange quark
* in openQCD simulations" (file doc/rhmc.pdf).
*
* The pseudo-fermion action Spf is given by
*
*   Spf=(phi,P_{k,l}*phi),
*
* where P_{k,l} is the fraction of a Zolotarev rational function, which
* is defined by the parameters:
*
*   irat[0]       Index of the Zolotarev rational function in the
*                 parameter data base.
*
*   irat[1]       Lower end k of the selected coefficient range.
*
*   irat[2]       Upper end l of the selected coefficient range.
*
* See ratfcts/ratfcts.c for further explanations. The inclusion of the
* "small quark determinant" amounts to adding the action
*
*   Sdet=-ln{det(1e+Doo)}+constant
*
* to the molecular-dynamics Hamilton function, where 1e is the projector
* to the quark fields that vanish on the odd lattice sites and Doo the
* odd-odd component of the Dirac operator (the constant is adjusted so
* as to reduce the significance losses when the action differences are
* computed at the end of the molecular-dynamics trajectories).
*
* The other parameters of the programs in this module are:
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
*   status        Array of the average status values returned by the
*                 solver used for the solution of the Dirac equation
*                 (in the case of the DFL_SAP_GCR solver, status[2]
*                 and status[5] are not averaged).
*
* The supported solvers are MSCG, SAP_GCR and DFL_SAP_GCR. Depending
* on the program and the solver, the number of status variables varies
* and is given by:
*
*                  MSCG         SAP_GCR       DFL_SAP_GCR
*   setpf3()         1             1               3
*   force3()         1             2               6
*   action3()        1             1               3
*
* Note that, in force3(), the GCR solvers solve the Dirac equations twice.
* In these cases, the program writes the status values one after the other
* to the array. The bare quark mass m0 is the one last set by sw_parms()
* [flags/lat_parms.c] and it is taken for granted that the parameters of
* the solver have been set by set_solver_parms() [flags/solver_parms.c].
*
* The required workspaces of double-precision spinor fields are
*
*                  MSCG         SAP_GCR       DFL_SAP_GCR
*   setpf3()        np             2               2
*   force3()        np             3               3
*   action3()       np             1               1
*
* where np is the number of poles of P_{k,l} (these figures do not include
* the workspace required by the solvers).
*
* The programs in this module perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define FORCE3_C

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
#include "ratfcts.h"
#include "forces.h"
#include "global.h"

#define N0 (NPROC0*L0)

static int isx,init=0,nps=0;
static double *rs;


static void set_res(int np,double res)
{
   int k;

   if (np>nps)
   {
      if (nps>0)
         free(rs);

      rs=malloc(np*sizeof(*rs));
      error(rs==NULL,1,"set_res [force3.c]",
            "Unable to allocate auxiliary array");
   }

   for (k=0;k<np;k++)
      rs[k]=res;
}


static double sdet(void)
{
   int bc,ix,iy,t,ie;
   double c,p,smx;
   complex_dble z;
   pauli_dble *m;
   sw_parms_t swp;

   if (init==0)
   {
      isx=init_hsum(1);
      init=1;
   }

   swp=sw_parms();

   if ((4.0+swp.m0)>1.0)
      c=pow(4.0+swp.m0,-6.0);
   else
      c=1.0;

   sw_term(NO_PTS);
   reset_hsum(isx);
   m=swdfld()+VOLUME;
   bc=bc_type();
   ix=(VOLUME/2);
   ie=0;

   while (ix<VOLUME)
   {
      p=1.0;
      iy=ix+8;
      if (iy>VOLUME)
         iy=VOLUME;

      for (;ix<iy;ix++)
      {
         t=global_time(ix);

         if (((t>0)||(bc==3))&&((t<(N0-1))||(bc!=0)))
         {
            z=det_pauli_dble(0.0,m);

            if (z.re>0.0)
               p*=(c*z.re);
            else
               ie=1;

            z=det_pauli_dble(0.0,m+1);

            if (z.re>0.0)
               p*=(c*z.re);
            else
               ie=1;
         }

         m+=2;
      }

      if (p>0.0)
      {
         smx=-log(p);
         add_to_hsum(isx,&smx);
      }
      else
         ie=1;
   }

   error(ie!=0,1,"sdet [force3.c]",
         "SW term has negative or vanishing determinant");
   local_hsum(isx,&smx);

   return smx;
}


double setpf3(int *irat,int ipf,int isw,int isp,int icom,int *status)
{
   int np,k,l,stat[3];
   double r,act,*nu,*rnu;
   complex_dble z;
   spinor_dble *phi,**wsd,**rsd;
   mdflds_t *mdfs;
   ratfct_t rf;
   tm_parms_t tm;
   solver_parms_t sp;
   sap_parms_t sap;

   tm=tm_parms();
   if (tm.eoflg!=1)
      set_tm_parms(1);

   mdfs=mdflds();
   phi=(*mdfs).pf[ipf];
   random_sd(VOLUME/2,phi,1.0);
   bnd_sd2zero(EVEN_PTS,phi);
   set_sd2zero(VOLUME/2,phi+(VOLUME/2));

   rf=ratfct(irat);
   np=rf.np;
   nu=rf.nu;
   rnu=rf.rnu;

   sp=solver_parms(isp);

   if (isw==1)
      act=sdet();
   else
      act=0.0;

   if (sp.solver==MSCG)
   {
      rsd=reserve_wsd(np);

      set_res(np,sp.res);
      tmcgm(sp.nmx,rs,np,nu,phi,rsd,status);

      error_root(status[0]<0,1,"setpf3 [force3.c]","MSCG solver failed "
                 "(irat=%d,%d,%d, isp=%d, status=%d)",
                 irat[0],irat[1],irat[2],isp,status[0]);

      wsd=reserve_wsd(2);
      set_sd2zero(VOLUME/2,wsd[0]);

      for (k=0;k<np;k++)
      {
         Dwhat_dble(-nu[k],rsd[k],wsd[1]);
         mulg5_dble(VOLUME/2,wsd[1]);
         mulr_spinor_add_dble(VOLUME/2,wsd[0],wsd[1],rnu[k]);
         act-=2.0*nu[k]*rnu[k]*norm_square_dble(VOLUME/2,0,wsd[1]);
      }

      act-=norm_square_dble(VOLUME/2,0,wsd[0]);
      z.re=0.0;
      z.im=1.0;
      mulc_spinor_add_dble(VOLUME/2,phi,wsd[0],z);

      release_wsd();
      release_wsd();
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      rsd=reserve_wsd(1);
      wsd=reserve_wsd(1);
      set_sd2zero(VOLUME/2,wsd[0]);
      mulg5_dble(VOLUME/2,phi);
      status[0]=0;

      for (k=0;k<np;k++)
      {
         sap_gcr(sp.nkv,sp.nmx,sp.res,nu[k],phi,rsd[0],stat);

         error_root(stat[0]<0,1,"setpf3 [force3.c]","SAP_GCR solver failed "
                    "(irat=%d,%d,%d, isp=%d, status=%d)",
                    irat[0],irat[1],irat[2],isp,stat[0]);

         status[0]+=stat[0];
         mulr_spinor_add_dble(VOLUME/2,wsd[0],rsd[0],rnu[k]);
         act-=2.0*nu[k]*rnu[k]*norm_square_dble(VOLUME/2,0,rsd[0]);
      }

      status[0]=(status[0]+(np/2))/np;
      mulg5_dble(VOLUME/2,phi);
      act-=norm_square_dble(VOLUME/2,0,wsd[0]);
      z.re=0.0;
      z.im=1.0;
      mulc_spinor_add_dble(VOLUME/2,phi,wsd[0],z);

      release_wsd();
      release_wsd();
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      rsd=reserve_wsd(1);
      wsd=reserve_wsd(1);
      set_sd2zero(VOLUME/2,wsd[0]);
      mulg5_dble(VOLUME/2,phi);

      for (l=0;l<3;l++)
         status[l]=0;

      for (k=0;k<np;k++)
      {
         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,nu[k],phi,rsd[0],stat);

         error_root((stat[0]<0)||(stat[1]<0),1,
                    "setpf3 [force3.c]","DFL_SAP_GCR solver failed "
                    "(irat=%d,%d,%d, isp=%d, status=%d,%d,%d)",
                    irat[0],irat[1],irat[2],isp,stat[0],stat[1],stat[2]);

         for (l=0;l<3;l++)
            status[l]+=stat[l];

         mulr_spinor_add_dble(VOLUME/2,wsd[0],rsd[0],rnu[k]);
         act-=2.0*nu[k]*rnu[k]*norm_square_dble(VOLUME/2,0,rsd[0]);
      }

      for (l=0;l<2;l++)
         status[l]=(status[l]+(np/2))/np;

      mulg5_dble(VOLUME/2,phi);
      act-=norm_square_dble(VOLUME/2,0,wsd[0]);
      z.re=0.0;
      z.im=1.0;
      mulc_spinor_add_dble(VOLUME/2,phi,wsd[0],z);

      release_wsd();
      release_wsd();
   }
   else
      error_root(1,1,"setpf3 [force3.c]","Unknown solver");

   if (icom==1)
   {
      r=act;
      MPI_Reduce(&r,&act,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&act,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return act;
}


void force3(int *irat,int ipf,int isw,int isp,double c,int *status)
{
   int np,k,l,ifail,stat[6];
   double *mu,*rmu;
   spinor_dble *phi,**rsd,**wsd;
   mdflds_t *mdfs;
   ratfct_t rf;
   tm_parms_t tm;
   solver_parms_t sp;
   sap_parms_t sap;

   tm=tm_parms();
   if (tm.eoflg!=1)
      set_tm_parms(1);

   mdfs=mdflds();
   phi=(*mdfs).pf[ipf];

   rf=ratfct(irat);
   np=rf.np;
   mu=rf.mu;
   rmu=rf.rmu;

   sp=solver_parms(isp);
   set_xt2zero();
   set_xv2zero();

   if (sp.solver==MSCG)
   {
      rsd=reserve_wsd(np);

      set_res(np,sp.res);
      tmcgm(sp.nmx,rs,np,mu,phi,rsd,status);

      error_root(status[0]<0,1,"force3 [force3.c]","MSCG solver failed "
                 "(irat=%d,%d,%d, isp=%d, status=%d)",
                 irat[0],irat[1],irat[2],isp,status[0]);

      wsd=reserve_wsd(1);

      for (k=0;k<np;k++)
      {
         Dwoe_dble(rsd[k],rsd[k]);
         Dwoo_dble(0.0,rsd[k],rsd[k]);

         Dwhat_dble(-mu[k],rsd[k],wsd[0]);
         mulg5_dble(VOLUME/2,wsd[0]);
         Dwoe_dble(wsd[0],wsd[0]);
         Dwoo_dble(0.0,wsd[0],wsd[0]);

         add_prod2xt(rmu[k],rsd[k],wsd[0]);
         add_prod2xv(-rmu[k],rsd[k],wsd[0]);
      }

      release_wsd();
      release_wsd();
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      rsd=reserve_wsd(2);
      wsd=reserve_wsd(1);
      mulg5_dble(VOLUME/2,phi);
      set_sd2zero(VOLUME/2,phi+(VOLUME/2));
      set_sd2zero(VOLUME/2,wsd[0]+(VOLUME/2));

      for (l=0;l<2;l++)
         status[l]=0;

      for (k=0;k<np;k++)
      {
         sap_gcr(sp.nkv,sp.nmx,sp.res,mu[k],phi,rsd[0],stat);
         assign_sd2sd(VOLUME/2,rsd[0],wsd[0]);
         mulg5_dble(VOLUME/2,wsd[0]);
         sap_gcr(sp.nkv,sp.nmx,sp.res,-mu[k],wsd[0],rsd[1],stat+1);

         error_root((stat[0]<0)||(stat[1]<0),1,"force3 [force3.c]",
                    "SAP_GCR solver failed (irat=%d,%d,%d, isp=%d, "
                    "status=%d;%d)",irat[0],irat[1],irat[2],
                    isp,stat[0],stat[1]);

         for (l=0;l<2;l++)
            status[l]+=stat[l];

         add_prod2xt(rmu[k],rsd[1],rsd[0]);
         add_prod2xv(rmu[k],rsd[1],rsd[0]);
      }

      for (l=0;l<2;l++)
         status[l]=(status[l]+(np/2))/np;

      mulg5_dble(VOLUME/2,phi);
      release_wsd();
      release_wsd();
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      rsd=reserve_wsd(2);
      wsd=reserve_wsd(1);
      mulg5_dble(VOLUME/2,phi);
      set_sd2zero(VOLUME/2,phi+(VOLUME/2));
      set_sd2zero(VOLUME/2,wsd[0]+(VOLUME/2));

      for (l=0;l<6;l++)
         status[l]=0;

      for (k=0;k<np;k++)
      {
         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu[k],phi,rsd[0],stat);
         assign_sd2sd(VOLUME/2,rsd[0],wsd[0]);
         mulg5_dble(VOLUME/2,wsd[0]);
         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,-mu[k],wsd[0],rsd[1],stat+3);

         error_root((stat[0]<0)||(stat[1]<0)||(stat[3]<0)||(stat[4]<0),1,
                    "force3 [force3.c]","DFL_SAP_GCR solver failed "
                    "(irat=%d,%d,%d, isp=%d, ""status=%d,%d,%d;%d,%d,%d)",
                    irat[0],irat[1],irat[2],isp,
                    stat[0],stat[1],stat[2],stat[3],stat[4],stat[5]);

         for (l=0;l<6;l++)
            status[l]+=stat[l];

         add_prod2xt(rmu[k],rsd[1],rsd[0]);
         add_prod2xv(rmu[k],rsd[1],rsd[0]);
      }

      for (l=0;l<2;l++)
      {
         status[l]=(status[l]+(np/2))/np;
         status[l+3]=(status[l+3]+(np/2))/np;
      }

      mulg5_dble(VOLUME/2,phi);
      release_wsd();
      release_wsd();
   }
   else
      error_root(1,1,"force3 [force3.c]","Unknown solver");

   if (isw==1)
   {
      ifail=add_det2xt(1.0,ODD_PTS);

      error_root(ifail!=0,1,"force3 [force3.c]",
                 "Inversion of the SW term was not safe");
   }

   sw_frc(c);
   hop_frc(c);
}


double action3(int *irat,int ipf,int isw,int isp,int icom,int *status)
{
   int np,k,l,stat[3];
   double act,r,*mu,*rmu;
   spinor_dble *phi,**rsd,**wsd;
   mdflds_t *mdfs;
   ratfct_t rf;
   solver_parms_t sp;
   sap_parms_t sap;
   tm_parms_t tm;

   tm=tm_parms();
   if (tm.eoflg!=1)
      set_tm_parms(1);

   mdfs=mdflds();
   phi=(*mdfs).pf[ipf];

   rf=ratfct(irat);
   np=rf.np;
   mu=rf.mu;
   rmu=rf.rmu;

   sp=solver_parms(isp);

   if (isw==1)
      act=sdet();
   else
      act=0.0;

   if (sp.solver==MSCG)
   {
      rsd=reserve_wsd(np);

      set_res(np,sp.res);
      tmcgm(sp.nmx,rs,np,mu,phi,rsd,status);

      error_root(status[0]<0,1,"action3 [force3.c]","MSCG solver failed "
                 "(irat=%d,%d,%d, isp=%d, status=%d)",
                 irat[0],irat[1],irat[2],isp,status[0]);

      wsd=reserve_wsd(1);

      for (k=0;k<np;k++)
      {
         Dwhat_dble(-mu[k],rsd[k],wsd[0]);
         act+=rmu[k]*norm_square_dble(VOLUME/2,0,wsd[0]);
      }

      release_wsd();
      release_wsd();
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);
      rsd=reserve_wsd(1);
      mulg5_dble(VOLUME/2,phi);
      set_sd2zero(VOLUME/2,phi+(VOLUME/2));
      status[0]=0;

      for (k=0;k<np;k++)
      {
         sap_gcr(sp.nkv,sp.nmx,sp.res,mu[k],phi,rsd[0],stat);

         error_root(stat[0]<0,1,"action3 [force3.c]",
                    "SAP_GCR solver failed (irat=%d,%d,%d, isp=%d, "
                    "status=%d)",irat[0],irat[1],irat[2],
                    isp,stat[0]);

         status[0]+=stat[0];
         act+=rmu[k]*norm_square_dble(VOLUME/2,0,rsd[0]);
      }

      status[0]=(status[0]+(np/2))/np;

      mulg5_dble(VOLUME/2,phi);
      release_wsd();
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);
      rsd=reserve_wsd(1);
      mulg5_dble(VOLUME/2,phi);
      set_sd2zero(VOLUME/2,phi+(VOLUME/2));

      for (l=0;l<3;l++)
         status[l]=0;

      for (k=0;k<np;k++)
      {
         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu[k],phi,rsd[0],stat);

         error_root((stat[0]<0)||(stat[1]<0),1,"action3 [force3.c]",
                    "DFL_SAP_GCR solver failed (irat=%d,%d,%d, isp=%d, "
                    "status=%d,%d,%d)",irat[0],irat[1],irat[2],isp,
                    stat[0],stat[1],stat[2]);

         for (l=0;l<3;l++)
            status[l]+=stat[l];

         act+=rmu[k]*norm_square_dble(VOLUME/2,0,rsd[0]);
      }

      for (l=0;l<2;l++)
         status[l]=(status[l]+(np/2))/np;

      mulg5_dble(VOLUME/2,phi);
      release_wsd();
   }
   else
      error_root(1,1,"action3 [force3.c]","Unknown solver");

   if (icom==1)
   {
      r=act;
      MPI_Reduce(&r,&act,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&act,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return act;
}
