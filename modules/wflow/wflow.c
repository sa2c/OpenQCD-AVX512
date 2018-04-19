
/*******************************************************************************
*
* File wflow.c
*
* Copyright (C) 2009-2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Integration of the Wilson flow.
*
* The externally accessible functions are
*
*   void fwd_euler(int n,double eps)
*     Applies n forward Euler integration steps, with step size eps, to the
*     current gauge field.
*
*   void fwd_rk2(int n,double eps)
*     Applies n forward 2nd-order Runge-Kutta integration steps, with step
*     size eps, to the current gauge field.
*
*   void fwd_rk3(int n,double eps)
*     Applies n forward 3rd-order Runge-Kutta integration steps, with step
*     size eps, to the current gauge field.
*
* Notes:
*
* On lattices with periodic boundary conditions, the Wilson flow is defined
* through equations (1.3) and (1.4) in
*
*   M. Luescher: "Properties and uses of the Wilson flow in lattice QCD",
*   JHEP 1008 (2010) 071.
*
* The Runge-Kutta integrators used here are described in appendix C of this
* paper.
*
* In the case of open, SF and open-SF boundary conditions, the flow evolves
* the active link variables only, where the action on the right of the flow
* equation is taken to be the tree-level O(a)-improved plaquette action.
* O(a)-improvement moreover requires the derivative of the action on the
* space-like links at the boundaries with open boundary conditions to be
* multiplied by 2. See appendix B in
*
*   M. Luescher, S. Schaefer: "Lattice QCD with open boundary conditions
*   and twisted-mass reweighting", Comput.Phys.Commun. 184 (2013) 519,
*
* for the exact form of the flow equation in this case.
*
* The integration programs make use of the force field in the structure
* returned by mdflds [mdflds.c]. On exit the force field must therefore be
* expected to be changed.
*
* The Runge-Kutta integrators require a workspace of 1 force field. All
* programs in this module involve global communications and must be called
* on all MPI processes simultaneously with the same values of the parameters.
*
*******************************************************************************/

#define WFLOW_C

#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "flags.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "mdflds.h"
#include "linalg.h"
#include "forces.h"
#include "wflow.h"
#include "global.h"

#define N0 (NPROC0*L0)


static void update_ud(double eps,su3_alg_dble *frc)
{
   int bc,ix,t,ifc;
   su3_dble *u;

   bc=bc_type();
   u=udfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         expXsu3(eps,frc,u);
         frc+=1;
         u+=1;

         if (bc!=0)
            expXsu3(eps,frc,u);
         frc+=1;
         u+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if ((bc==0)||(bc==2))
               expXsu3(2.0*eps,frc,u);
            else if (bc==3)
               expXsu3(eps,frc,u);
            frc+=1;
            u+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            expXsu3(eps,frc,u);
         frc+=1;
         u+=1;

         expXsu3(eps,frc,u);
         frc+=1;
         u+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc==0)
               expXsu3(2.0*eps,frc,u);
            else
               expXsu3(eps,frc,u);
            frc+=1;
            u+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            expXsu3(eps,frc,u);
            frc+=1;
            u+=1;
         }
      }
   }

   set_flags(UPDATED_UD);
}


static void update_fro1(double c,su3_alg_dble *frc,su3_alg_dble *fro)
{
   su3_alg_dble *frm;

   frm=frc+4*VOLUME;

   for (;frc<frm;frc++)
   {
      (*fro).c1-=c*(*frc).c1;
      (*fro).c2-=c*(*frc).c2;
      (*fro).c3-=c*(*frc).c3;
      (*fro).c4-=c*(*frc).c4;
      (*fro).c5-=c*(*frc).c5;
      (*fro).c6-=c*(*frc).c6;
      (*fro).c7-=c*(*frc).c7;
      (*fro).c8-=c*(*frc).c8;

      fro+=1;
   }
}


static void update_fro2(double c,su3_alg_dble *frc,su3_alg_dble *fro)
{
   su3_alg_dble *frm;

   frm=frc+4*VOLUME;

   for (;frc<frm;frc++)
   {
      (*fro).c1=(*frc).c1+c*(*fro).c1;
      (*fro).c2=(*frc).c2+c*(*fro).c2;
      (*fro).c3=(*frc).c3+c*(*fro).c3;
      (*fro).c4=(*frc).c4+c*(*fro).c4;
      (*fro).c5=(*frc).c5+c*(*fro).c5;
      (*fro).c6=(*frc).c6+c*(*fro).c6;
      (*fro).c7=(*frc).c7+c*(*fro).c7;
      (*fro).c8=(*frc).c8+c*(*fro).c8;

      fro+=1;
   }
}


void fwd_euler(int n,double eps)
{
   int iprms[1],k;
   double dprms[1];
   su3_alg_dble *frc;
   mdflds_t *mdfs;

   if (NPROC>1)
   {
      iprms[0]=n;
      dprms[0]=eps;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      error((iprms[0]!=n)||(dprms[0]!=eps),1,
            "fwd_euler [wflow.c]","Parameters are not global");
   }

   if (n>0)
   {
      mdfs=mdflds();
      frc=(*mdfs).frc;

      for (k=0;k<n;k++)
      {
         plaq_frc();
         update_ud(-eps,frc);
      }
   }
}


void fwd_rk2(int n,double eps)
{
   int iprms[1],k;
   double dprms[1];
   su3_alg_dble *frc,*fro,**fsv;
   mdflds_t *mdfs;

   if (NPROC>1)
   {
      iprms[0]=n;
      dprms[0]=eps;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      error((iprms[0]!=n)||(dprms[0]!=eps),1,
            "fwd_rk2 [wflow.c]","Parameters are not global");
   }

   if (n>0)
   {
      mdfs=mdflds();
      frc=(*mdfs).frc;
      fsv=reserve_wfd(1);
      fro=fsv[0];

      for (k=0;k<n;k++)
      {
         plaq_frc();
         assign_alg2alg(4*VOLUME,frc,fro);
         update_ud(-0.5*eps,frc);

         plaq_frc();
         update_fro2(-0.5,frc,fro);
         update_ud(-eps,fro);
      }

      release_wfd();
   }
}


void fwd_rk3(int n,double eps)
{
   int iprms[1],k;
   double dprms[1];
   su3_alg_dble *frc,*fro,**fsv;
   mdflds_t *mdfs;

   if (NPROC>1)
   {
      iprms[0]=n;
      dprms[0]=eps;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      error((iprms[0]!=n)||(dprms[0]!=eps),1,
            "fwd_rk3 [wflow.c]","Parameters are not global");
   }

   if (n>0)
   {
      mdfs=mdflds();
      frc=(*mdfs).frc;
      fsv=reserve_wfd(1);
      fro=fsv[0];

      for (k=0;k<n;k++)
      {
         plaq_frc();
         assign_alg2alg(4*VOLUME,frc,fro);
         update_ud(-0.25*eps,frc);

         plaq_frc();
         update_fro1(32.0/17.0,frc,fro);
         update_ud((17.0/36.0)*eps,fro);

         plaq_frc();
         update_fro2(17.0/27.0,frc,fro);
         update_ud(-0.75*eps,fro);
      }

      release_wfd();
   }
}
