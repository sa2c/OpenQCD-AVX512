
/*******************************************************************************
*
* File force4.c
*
* Copyright (C) 2012, 2013 Martin Luescher, Stefan Schaefer
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Twisted mass pseudo-fermion action and force with even-odd preconditioning.
*
* The externally accessible functions are
*
*   double setpf4(double mu,int ipf,int isw,int icom)
*     Generates a pseudo-fermion field phi with probability proportional
*     to exp(-Spf) and returns the action Spf+Sdet if isw=1 or Spf if
*     isw!=1 (see the notes).
*
*   void force4(double mu,int ipf,int isw,int isp,int icr,double c,
*               int *status)
*     Computes the force deriving from the action Spf+Sdet if isw=1 or
*     Spf if isw!=1 (see the notes). The calculated force is multiplied
*     by c and added to the molecular-dynamics force field.
*
*   double action4(double mu,int ipf,int isw,int isp,int icom,
*                  int *status)
*     Returns the action Spf+Sdet if isw=1 or Spf if isw!=1 (see the
*     notes).
*
* Notes:
*
* The pseudo-fermion action Spf is given by
*
*   Spf=(phi,(Dwhat^dag*Dwhat+mu^2)^(-1)*phi),
*
* where Dwhat denotes the even-odd preconditioned (improved) Wilson-Dirac
* operator and phi the pseudo-fermion field. The latter vanishes on the
* odd lattice sites.
*
* The inclusion of the "small quark determinant" amounts to adding the
* action
*
*   Sdet=-2*ln{det(1e+Doo)}+constant
*
* to the molecular-dynamics Hamilton function, where 1e is the projector
* to the quark fields that vanish on the odd lattice sites and Doo the
* odd-odd component of the Dirac operator (the constant is adjusted so
* as to reduce the significance losses when the action differences are
* computed at the end of the molecular-dynamics trajectories).
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
*   icom          The action returned by the programs setpf4() and
*                 action4() is summed over all MPI processes if icom=1.
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
*   force4()         1             2               6
*   action4()        1             1               3
*
* Note that, in force4(), the GCR solvers solve the Dirac equations twice.
* In these cases, the program writes the status values one after the other
* to the array. The bare quark mass m0 is the one last set by sw_parms()
* [flags/lat_parms.c] and it is taken for granted that the parameters of
* the solver have been set by set_solver_parms() [flags/solver_parms.c].
*
* The program force4() attempts to propagate the solutions of the Dirac
* equation along the molecular-dynamics trajectories, using the field
* stack number icr (no fields are propagated if icr=0). If this feature
* is used, the program setup_chrono() [update/chrono.c] must be called
* before force4() is called for the first time.
*
* The required workspaces of double-precision spinor fields are
*
*                  CGNE         SAP_GCR       DFL_SAP_GCR
*   setpf4()         1             1               1
*   force4()     2+(icr>0)     2+2*(icr>0)     2+2*(icr>0)
*   action4()        1             1               1
*
* (these figures do not include the workspace required by the solvers).
*
* The programs in this module perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define FORCE4_C

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

#define N0 (NPROC0*L0)

static int isx,init=0;


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

      if (p!=0.0)
      {
         smx=-2.0*log(p);
         add_to_hsum(isx,&smx);
      }
      else
         ie=1;
   }

   error(ie!=0,1,"sdet [force4.c]",
         "SW term has negative or vanishing determinant");
   local_hsum(isx,&smx);

   return smx;
}


double setpf4(double mu,int ipf,int isw,int icom)
{
   int ifail;
   double act,r;
   spinor_dble *phi,**wsd;
   mdflds_t *mdfs;

   if (isw==1)
      act=sdet();
   else
      act=0.0;

   wsd=reserve_wsd(1);
   random_sd(VOLUME/2,wsd[0],1.0);
   bnd_sd2zero(EVEN_PTS,wsd[0]);
   act+=norm_square_dble(VOLUME/2,0,wsd[0]);

   ifail=sw_term(ODD_PTS);
   error_root(ifail!=0,1,"set_pf4 [force4.c]",
              "Inversion of the SW term was not safe");

   mdfs=mdflds();
   phi=(*mdfs).pf[ipf];
   Dwhat_dble(mu,wsd[0],phi);
   mulg5_dble(VOLUME/2,phi);
   set_sd2zero(VOLUME/2,phi+(VOLUME/2));
   release_wsd();

   if (icom==1)
   {
      r=act;
      MPI_Reduce(&r,&act,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&act,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return act;
}


void force4(double mu,int ipf,int isw,int isp,int icr,double c,int *status)
{
   int ifail,l;
   double res0,res1;
   spinor_dble *phi,*chi,*psi,**wsd;
   spinor_dble *rho,*eta,**rsd;
   mdflds_t *mdfs;
   solver_parms_t sp;
   sap_parms_t sap;
   tm_parms_t tm;

   tm=tm_parms();
   if (tm.eoflg!=1)
      set_tm_parms(1);

   sp=solver_parms(isp);

   mdfs=mdflds();
   phi=(*mdfs).pf[ipf];
   wsd=reserve_wsd(2);
   psi=wsd[0];
   chi=wsd[1];

   set_xt2zero();
   set_xv2zero();

   if (isw==1)
   {
      ifail=add_det2xt(2.0,ODD_PTS);
      error_root(ifail!=0,1,"force4 [force4.c]",
                 "Inversion of the SW term was not safe");
   }

   if (sp.solver==CGNE)
   {
      if (get_chrono(icr,chi))
      {
         rsd=reserve_wsd(1);
         rho=rsd[0];

         ifail=sw_term(ODD_PTS);
         error_root(ifail!=0,1,"force4 [force4.c]",
                    "Inversion of the SW term was not safe");

         Dwhat_dble(-mu,chi,psi);
         mulg5_dble(VOLUME/2,psi);
         Dwhat_dble(mu,psi,rho);
         mulg5_dble(VOLUME/2,rho);
         mulr_spinor_add_dble(VOLUME/2,rho,phi,-1.0);

         res0=norm_square_dble(VOLUME/2,1,phi);
         res1=norm_square_dble(VOLUME/2,1,rho);
         res1=sqrt(res1/res0);

         if (res1<1.0)
         {
            if (res1>sp.res)
            {
               tmcgeo(sp.nmx,sp.res/res1,mu,rho,psi,status);
               mulr_spinor_add_dble(VOLUME/2,chi,psi,-1.0);
            }
            else
               status[0]=0;
         }
         else
            tmcgeo(sp.nmx,sp.res,mu,phi,chi,status);

         release_wsd();
      }
      else
         tmcgeo(sp.nmx,sp.res,mu,phi,chi,status);

      error_root(status[0]<0,1,"force4 [force4.c]",
                 "CGNE solver failed (mu = %.4e, parameter set no %d, "
                 "status = %d)",mu,isp,status[0]);

      Dwoe_dble(chi,chi);
      Dwoo_dble(0.0,chi,chi);
      Dwhat_dble(-mu,chi,psi);
      mulg5_dble(VOLUME/2,psi);
      Dwoe_dble(psi,psi);
      Dwoo_dble(0.0,psi,psi);

      if (icr)
         add_chrono(icr,chi);

      add_prod2xt(1.0,chi,psi);
      add_prod2xv(-1.0,chi,psi);
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

         ifail=sw_term(ODD_PTS);
         error_root(ifail!=0,1,"force4 [force4.c]",
                    "Inversion of the SW term was not safe");

         Dwhat_dble(-mu,chi,psi);
         mulg5_dble(VOLUME/2,psi);
         Dwhat_dble(mu,psi,rho);
         mulg5_dble(VOLUME/2,rho);
         mulr_spinor_add_dble(VOLUME/2,rho,phi,-1.0);

         res0=norm_square_dble(VOLUME/2,1,phi);
         res1=norm_square_dble(VOLUME/2,1,rho);
         res1=sqrt(res1/res0);

         if (res1<1.0)
         {
            Dwoe_dble(chi,chi);
            Dwoo_dble(0.0,chi,chi);
            scale_dble(VOLUME/2,-1.0,chi+(VOLUME/2));

            Dwoe_dble(psi,psi);
            Dwoo_dble(0.0,psi,psi);
            scale_dble(VOLUME/2,-1.0,psi+(VOLUME/2));

            if (res1>sp.res)
            {
               mulg5_dble(VOLUME/2,rho);
               set_sd2zero(VOLUME/2,rho+(VOLUME/2));

               sap_gcr(sp.nkv,sp.nmx,sp.res/res1,mu,rho,eta,status);

               mulr_spinor_add_dble(VOLUME,psi,eta,-1.0);

               res0=norm_square_dble(VOLUME/2,1,psi);
               res1=norm_square_dble(VOLUME/2,1,eta);
               res1=sqrt(res1/res0);

               if (res1<1.0)
               {
                  if (res1>sp.res)
                  {
                     mulg5_dble(VOLUME/2,eta);
                     set_sd2zero(VOLUME/2,eta+(VOLUME/2));

                     sap_gcr(sp.nkv,sp.nmx,sp.res/res1,-mu,eta,rho,status+1);

                     mulr_spinor_add_dble(VOLUME,chi,rho,-1.0);
                  }
                  else
                     status[1]=0;
               }
               else
               {
                  assign_sd2sd(VOLUME/2,psi,eta);
                  mulg5_dble(VOLUME/2,eta);
                  set_sd2zero(VOLUME/2,eta+(VOLUME/2));

                  sap_gcr(sp.nkv,sp.nmx,sp.res,-mu,eta,chi,status+1);
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
            mulg5_dble(VOLUME/2,phi);
            set_sd2zero(VOLUME/2,phi+(VOLUME/2));

            sap_gcr(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);

            mulg5_dble(VOLUME/2,phi);
            assign_sd2sd(VOLUME/2,psi,eta);
            mulg5_dble(VOLUME/2,eta);
            set_sd2zero(VOLUME/2,eta+(VOLUME/2));

            sap_gcr(sp.nkv,sp.nmx,sp.res,-mu,eta,chi,status+1);
         }

         release_wsd();
      }
      else
      {
         rsd=reserve_wsd(1);
         eta=rsd[0];

         mulg5_dble(VOLUME/2,phi);
         set_sd2zero(VOLUME/2,phi+(VOLUME/2));

         sap_gcr(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);

         mulg5_dble(VOLUME/2,phi);
         assign_sd2sd(VOLUME/2,psi,eta);
         mulg5_dble(VOLUME/2,eta);
         set_sd2zero(VOLUME/2,eta+(VOLUME/2));

         sap_gcr(sp.nkv,sp.nmx,sp.res,-mu,eta,chi,status+1);

         release_wsd();
      }

      error_root((status[0]<0)||(status[1]<0),1,"force4 [force4.c]",
                 "SAP_GCR solver failed (mu = %.4e, parameter set no %d, "
                 "status = %d;%d)",mu,isp,status[0],status[1]);

      if (icr)
         add_chrono(icr,chi);

      add_prod2xt(1.0,chi,psi);
      add_prod2xv(1.0,chi,psi);
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

         ifail=sw_term(ODD_PTS);
         error_root(ifail!=0,1,"force4 [force4.c]",
                    "Inversion of the SW term was not safe");

         Dwhat_dble(-mu,chi,psi);
         mulg5_dble(VOLUME/2,psi);
         Dwhat_dble(mu,psi,rho);
         mulg5_dble(VOLUME/2,rho);
         mulr_spinor_add_dble(VOLUME/2,rho,phi,-1.0);

         res0=norm_square_dble(VOLUME/2,1,phi);
         res1=norm_square_dble(VOLUME/2,1,rho);
         res1=sqrt(res1/res0);

         if (res1<1.0)
         {
            Dwoe_dble(chi,chi);
            Dwoo_dble(0.0,chi,chi);
            scale_dble(VOLUME/2,-1.0,chi+(VOLUME/2));

            Dwoe_dble(psi,psi);
            Dwoo_dble(0.0,psi,psi);
            scale_dble(VOLUME/2,-1.0,psi+(VOLUME/2));

            if (res1>sp.res)
            {
               mulg5_dble(VOLUME/2,rho);
               set_sd2zero(VOLUME/2,rho+(VOLUME/2));

               dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res/res1,mu,rho,eta,
                            status);

               mulr_spinor_add_dble(VOLUME,psi,eta,-1.0);

               res0=norm_square_dble(VOLUME/2,1,psi);
               res1=norm_square_dble(VOLUME/2,1,eta);
               res1=sqrt(res1/res0);

               if (res1<1.0)
               {
                  if (res1>sp.res)
                  {
                     mulg5_dble(VOLUME/2,eta);
                     set_sd2zero(VOLUME/2,eta+(VOLUME/2));

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
                  assign_sd2sd(VOLUME/2,psi,eta);
                  mulg5_dble(VOLUME/2,eta);
                  set_sd2zero(VOLUME/2,eta+(VOLUME/2));

                  dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,-mu,eta,chi,status+3);
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
            mulg5_dble(VOLUME/2,phi);
            set_sd2zero(VOLUME/2,phi+(VOLUME/2));

            dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);

            mulg5_dble(VOLUME/2,phi);
            assign_sd2sd(VOLUME/2,psi,eta);
            mulg5_dble(VOLUME/2,eta);
            set_sd2zero(VOLUME/2,eta+(VOLUME/2));

            dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,-mu,eta,chi,status+3);
         }

         release_wsd();
      }
      else
      {
         rsd=reserve_wsd(1);
         eta=rsd[0];

         mulg5_dble(VOLUME/2,phi);
         set_sd2zero(VOLUME/2,phi+(VOLUME/2));

         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);

         mulg5_dble(VOLUME/2,phi);
         assign_sd2sd(VOLUME/2,psi,eta);
         mulg5_dble(VOLUME/2,eta);
         set_sd2zero(VOLUME/2,eta+(VOLUME/2));

         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,-mu,eta,chi,status+3);

         release_wsd();
      }

      error_root((status[0]<0)||(status[1]<0)||(status[3]<0)||(status[4]<0),1,
                 "force4 [force4.c]","DFL_SAP_GCR solver failed "
                 "(mu = %.4e, parameter set no %d, status = %d,%d,%d;%d,%d,%d)",
                 mu,isp,status[0],status[1],status[2],status[3],
                 status[4],status[5]);

      if (icr)
         add_chrono(icr,chi);

      add_prod2xt(1.0,chi,psi);
      add_prod2xv(1.0,chi,psi);
   }
   else
      error_root(1,1,"force4 [force4.c]","Unknown solver");

   sw_frc(c);
   hop_frc(c);

   release_wsd();
}


double action4(double mu,int ipf,int isw,int isp,int icom,int *status)
{
   double act,r;
   spinor_dble *phi,*chi,*psi;
   spinor_dble **rsd,**wsd;
   mdflds_t *mdfs;
   solver_parms_t sp;
   sap_parms_t sap;
   tm_parms_t tm;

   tm=tm_parms();
   if (tm.eoflg!=1)
      set_tm_parms(1);

   mdfs=mdflds();
   phi=(*mdfs).pf[ipf];
   sp=solver_parms(isp);

   if (isw==1)
      act=sdet();
   else
      act=0.0;

   if (sp.solver==CGNE)
   {
      rsd=reserve_wsd(1);
      chi=rsd[0];

      tmcgeo(sp.nmx,sp.res,mu,phi,chi,status);

      error_root(status[0]<0,1,"action4 [force4.c]",
                 "CGNE solver failed (mu = %.4e, parameter set no %d, "
                 "status = %d)",mu,isp,status[0]);

      wsd=reserve_wsd(1);
      psi=wsd[0];

      Dwhat_dble(-mu,chi,psi);
      act+=norm_square_dble(VOLUME/2,0,psi);

      release_wsd();
      release_wsd();
   }
   else if (sp.solver==SAP_GCR)
   {
      rsd=reserve_wsd(1);
      psi=rsd[0];

      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);
      mulg5_dble(VOLUME/2,phi);
      set_sd2zero(VOLUME/2,phi+(VOLUME/2));

      sap_gcr(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);

      error_root(status[0]<0,1,"action4 [force4.c]",
                 "SAP_GCR solver failed (mu = %.4e, parameter set no %d, "
                 "status = %d)",mu,isp,status[0]);

      mulg5_dble(VOLUME/2,phi);
      act+=norm_square_dble(VOLUME/2,0,psi);

      release_wsd();
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      rsd=reserve_wsd(1);
      psi=rsd[0];

      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);
      mulg5_dble(VOLUME/2,phi);
      set_sd2zero(VOLUME/2,phi+(VOLUME/2));

      dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mu,phi,psi,status);

      error_root((status[0]<0)||(status[1]<0),1,"action4 [force4.c]",
                 "DFL_SAP_GCR solver failed (mu = %.4e, parameter set "
                 "no %d, status = %d,%d,%d)",mu,isp,
                 status[0],status[1],status[2]);

      mulg5_dble(VOLUME/2,phi);
      act+=norm_square_dble(VOLUME/2,0,psi);

      release_wsd();
   }
   else
      error_root(1,1,"action4 [force4.c]","Unknown solver");

   if (icom==1)
   {
      r=act;
      MPI_Reduce(&r,&act,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&act,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return act;
}
