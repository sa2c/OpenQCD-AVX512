
/*******************************************************************************
*
* File mdint.c
*
* Copyright (C) 2011-2013 Stefan Schaefer, Martin Luescher, John Bulava
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Integration of the molecular-dynamics equations.
*
* The externally accessible functions are
*
*   void run_mdint(void)
*     Integrates the molecular-dynamics equations using the current
*     integrator (see the notes).
*
* Notes:
*
* The integrator used is the one defined by the array of elementary operations
* returned by mdsteps() (see update/mdsteps.c). It is assumed that the fields
* and the integrator have been properly initialized.
*
* In the course of the integration, the solver iteration numbers are added
* to the appropriate counters provided by the module update/counters.c.
*
* The program in this module performs global communications and must be
* called simultaneously on all MPI processes.
*
* Some debugging information is printed to stdout if the macro MDINT_DBG is
* defined. The norm of the forces printed is the norm per active link.
*
*******************************************************************************/

#define MDINT_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "mpi.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "mdflds.h"
#include "su3fcts.h"
#include "linalg.h"
#include "dfl.h"
#include "forces.h"
#include "update.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int nsm;
static double rtau,dtau;


static void chk_mode_regen(int isp,int *status)
{
   int i,is;
   solver_parms_t sp;

   sp=solver_parms(isp);

   if (sp.solver==DFL_SAP_GCR)
   {
      is=status[2];

      for (i=2;i<4;i++)
         status[i]=status[i+1];

      status[4]=is;

      if (status[4]>0)
         add2counter("modes",2,status+4);
      if (status[5]>0)
         add2counter("modes",2,status+5);
   }
}


static void update_mom(void)
{
   int bc,ix,t,ifc;
   su3_alg_dble *mom,*frc;
   mdflds_t *mdfs;

   bc=bc_type();
   mdfs=mdflds();
   mom=(*mdfs).mom;
   frc=(*mdfs).frc;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         _su3_alg_sub_assign(mom[0],frc[0]);
         mom+=1;
         frc+=1;

         if (bc!=0)
         {
            _su3_alg_sub_assign(mom[0],frc[0]);
         }

         mom+=1;
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
            {
               _su3_alg_sub_assign(mom[0],frc[0]);
            }

            mom+=1;
            frc+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
         {
            _su3_alg_sub_assign(mom[0],frc[0]);
         }

         mom+=1;
         frc+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            _su3_alg_sub_assign(mom[0],frc[0]);
            mom+=1;
            frc+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            _su3_alg_sub_assign(mom[0],frc[0]);
            mom+=1;
            frc+=1;
         }
      }
   }
}


static void update_ud(double eps)
{
   int bc,ix,t,ifc;
   su3_dble *u;
   su3_alg_dble *mom;
   mdflds_t *mdfs;

   bc=bc_type();
   mdfs=mdflds();
   mom=(*mdfs).mom;
   u=udfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         expXsu3(eps,mom,u);
         u+=1;
         mom+=1;

         if (bc!=0)
            expXsu3(eps,mom,u);
         u+=1;
         mom+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
               expXsu3(eps,mom,u);
            u+=1;
            mom+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            expXsu3(eps,mom,u);
         u+=1;
         mom+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            u+=1;
            mom+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            u+=1;
            mom+=1;
         }
      }
   }

   set_flags(UPDATED_UD);
}


static void start_dfl_upd(void)
{
   dfl_upd_parms_t dup;

   dup=dfl_upd_parms();
   dtau=dup.dtau;
   nsm=dup.nsm;
   rtau=0.0;
}


static void dfl_upd(int isp)
{
   int status[2];
   solver_parms_t sp;

   if ((nsm>0)&&(rtau>dtau))
   {
      sp=solver_parms(isp);

      if (sp.solver==DFL_SAP_GCR)
      {
         dfl_update2(nsm,status);
         error_root((status[1]<0)||((status[1]==0)&&(status[0]<0)),1,
                    "dfl_upd [mdint.c]","Deflation subspace update "
                    "failed (status = %d;%d)",status[0],status[1]);

         if (status[1]==0)
            add2counter("modes",1,status);
         else
            add2counter("modes",2,status+1);

         rtau=0.0;
      }
   }
}

#ifdef MDINT_DBG

void run_mdint(void)
{
   int my_rank,nop,itu;
   int iop,status[6];
   double *mu,eps,nlk,nrm;
   mdflds_t *mdfs;
   mdstep_t *s,*sm;
   hmc_parms_t hmc;
   force_parms_t fp;
   double wt1, wt2;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   mdfs=mdflds();
   hmc=hmc_parms();
   mu=hmc.mu;
   reset_chrono();
   start_dfl_upd();

   nlk=(double)(4*N0*N1)*(double)(N2*N3);
   if (bc_type()==0)
      nlk-=(double)(N1)*(double)(N2*N3);
   else if (bc_type()==1)
      nlk-=(double)(3*N1)*(double)(N2*N3);

   s=mdsteps(&nop,&itu);
   sm=s+nop;

   for (;s<sm;s++)
   {
      iop=(*s).iop;
      eps=(*s).eps;

      if (iop<itu)
      {
         fp=force_parms(iop);

         MPI_Barrier(MPI_COMM_WORLD);
         wt1=MPI_Wtime();

         if (fp.force==FRG)
            force0(eps);
         else
         {
            dfl_upd(fp.isp[0]);
            set_sw_parms(sea_quark_mass(fp.im0));
            set_frc2zero();
            status[2]=0;
            status[5]=0;

            if (fp.force==FRF_TM1)
               force1(mu[fp.imu[0]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM1_EO)
               force4(mu[fp.imu[0]],fp.ipf,0,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM1_EO_SDET)
               force4(mu[fp.imu[0]],fp.ipf,1,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM2)
               force2(mu[fp.imu[0]],mu[fp.imu[1]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM2_EO)
               force5(mu[fp.imu[0]],mu[fp.imu[1]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_RAT)
               force3(fp.irat,fp.ipf,0,fp.isp[0],
                      eps,status);
            else if (fp.force==FRF_RAT_SDET)
               force3(fp.irat,fp.ipf,1,fp.isp[0],
                      eps,status);

            chk_mode_regen(fp.isp[0],status);
            add2counter("force",iop,status);
         }

         MPI_Barrier(MPI_COMM_WORLD);
         wt2=MPI_Wtime();

         update_mom();
         nrm=norm_square_alg(4*VOLUME,1,(*mdfs).frc);
         nrm=sqrt(nrm/nlk);

         if (my_rank==0)
         {
            if (fp.force==FRG)
               printf("Force FRG:              ");
            else if (fp.force==FRF_TM1)
               printf("Force FRF_TM1:          ");
            else if (fp.force==FRF_TM1_EO)
               printf("Force FRF_TM1_EO:       ");
            else if (fp.force==FRF_TM1_EO_SDET)
               printf("Force FRF_TM1_EO_SDET:  ");
            else if (fp.force==FRF_TM2)
               printf("Force FRF_TM2:          ");
            else if (fp.force==FRF_TM2_EO)
               printf("Force FRF_TM2_EO:       ");
            else if (fp.force==FRF_RAT)
               printf("Force FRF_RAT:          ");
            else if (fp.force==FRF_RAT_SDET)
               printf("Force FRF_RAT_SDET:     ");

            printf("nrm = %.2e, eps = % .2e, nrm*|eps| = %.2e, "
                   "time = %.2e sec\n",nrm/fabs(eps),eps,nrm,wt2-wt1);
         }
      }
      else if (iop==itu)
      {
         update_ud(eps);
         step_mdtime(eps);
         rtau+=eps;
      }
   }
}

#else

void run_mdint(void)
{
   int nop,itu;
   int iop,status[6];
   double *mu,eps;
   mdstep_t *s,*sm;
   hmc_parms_t hmc;
   force_parms_t fp;

   hmc=hmc_parms();
   mu=hmc.mu;
   reset_chrono();
   start_dfl_upd();

   s=mdsteps(&nop,&itu);
   sm=s+nop;

   for (;s<sm;s++)
   {
      iop=(*s).iop;
      eps=(*s).eps;

      if (iop<itu)
      {
         fp=force_parms(iop);

         if (fp.force==FRG)
            force0(eps);
         else
         {
            dfl_upd(fp.isp[0]);
            set_sw_parms(sea_quark_mass(fp.im0));
            status[2]=0;
            status[5]=0;

            if (fp.force==FRF_TM1)
               force1(mu[fp.imu[0]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM1_EO)
               force4(mu[fp.imu[0]],fp.ipf,0,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM1_EO_SDET)
               force4(mu[fp.imu[0]],fp.ipf,1,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM2)
               force2(mu[fp.imu[0]],mu[fp.imu[1]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM2_EO)
               force5(mu[fp.imu[0]],mu[fp.imu[1]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_RAT)
               force3(fp.irat,fp.ipf,0,fp.isp[0],
                      eps,status);
            else if (fp.force==FRF_RAT_SDET)
               force3(fp.irat,fp.ipf,1,fp.isp[0],
                      eps,status);

            chk_mode_regen(fp.isp[0],status);
            add2counter("force",iop,status);
         }
      }
      else if (iop==itu)
      {
         update_mom();
         update_ud(eps);
         step_mdtime(eps);
         rtau+=eps;
      }
      else
         update_mom();
   }
}

#endif
