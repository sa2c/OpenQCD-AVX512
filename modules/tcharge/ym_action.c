
/*******************************************************************************
*
* File ym_action.c
*
* Copyright (C) 2010-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the Yang-Mills action using the symmetric field tensor.
*
* The externally accessible functions are
*
*   double ym_action(void)
*     Returns the Yang-Mills action S (w/o prefactor 1/g0^2) of the
*     double-precision gauge field, using a symmetric expression for
*     the gauge-field tensor.
*
*   double ym_action_slices(double *asl)
*     Computes the sum asl[t] of the Yang-Mills action density (w/o
*     prefactor 1/g0^2) of the double-precision gauge field at time
*     t=0,1,...,N0-1 (where N0=NPROC0*L0). The program returns the
*     total action.
*
* Notes:
*
* The Yang-Mills action density s(x) is defined by
*
*  s(x)=(1/4)*sum_{mu,nu} [F_{mu,nu}^a(x)]^2
*
* where
*
*  F_{mu,nu}^a(x)=-2*tr{F_{mu,nu}(x)*T^a}, a=1,..,8,
*
* are the SU(3) components of the symmetric field tensor returned by the
* program ftensor() [ftensor.c]. At the boundaries of the lattice (if any),
* the action density is set to zero. The total action S is the sum of s(x)
* over all points x with time component in the range
*
*  0<x0<NPROC0*L0-1        (open bc),
*
*  0<x0<NPROC0*L0          (SF and open-SF bc),
*
*  0<=x0<NPROC0*L0         (periodic bc).
*
* The programs in this module perform global operations and must be called
* simultaneously on all MPI processes.
*
*******************************************************************************/

#define YM_ACTION_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "linalg.h"
#include "tcharge.h"
#include "global.h"

#define N0 (NPROC0*L0)

static int isx[L0],init=0;
static double asl0[N0];
static u3_alg_dble **ft;


static double prodXX(u3_alg_dble *X)
{
   double sm;

   sm=(-2.0/3.0)*((*X).c1+(*X).c2+(*X).c3)*((*X).c1+(*X).c2+(*X).c3)+
      2.0*((*X).c1*(*X).c1+(*X).c2*(*X).c2+(*X).c3*(*X).c3)+
      4.0*((*X).c4*(*X).c4+(*X).c5*(*X).c5+(*X).c6*(*X).c6+
           (*X).c7*(*X).c7+(*X).c8*(*X).c8+(*X).c9*(*X).c9);

   return sm;
}


static double density(int ix)
{
   double sm;

   sm=prodXX(ft[0]+ix)+prodXX(ft[1]+ix)+prodXX(ft[2]+ix)+
      prodXX(ft[3]+ix)+prodXX(ft[4]+ix)+prodXX(ft[5]+ix);

   return sm;
}


double ym_action(void)
{
   int bc,ix,t,tmx;
   double S;

   if (init==0)
   {
      for (t=0;t<L0;t++)
         isx[t]=init_hsum(1);

      init=1;
   }

   ft=ftensor();
   bc=bc_type();
   if (bc==0)
      tmx=N0-1;
   else
      tmx=N0;
   reset_hsum(isx[0]);

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (((t>0)&&(t<tmx))||(bc==3))
      {
         S=density(ix);
         add_to_hsum(isx[0],&S);
      }
   }

   if (NPROC>1)
      global_hsum(isx[0],&S);
   else
      local_hsum(isx[0],&S);

   return 0.5*S;
}


double ym_action_slices(double *asl)
{
   int bc,ix,t,t0,tmx;
   double S;

   if (init==0)
   {
      for (t=0;t<L0;t++)
         isx[t]=init_hsum(1);

      init=1;
   }

   ft=ftensor();
   bc=bc_type();
   if (bc==0)
      tmx=N0-1;
   else
      tmx=N0;
   t0=cpr[0]*L0;

   for (t=0;t<L0;t++)
      reset_hsum(isx[t]);

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (((t>0)&&(t<tmx))||(bc==3))
      {
         t-=t0;
         S=density(ix);
         add_to_hsum(isx[t],&S);
      }
   }

   for (t=0;t<N0;t++)
      asl0[t]=0.0;

   for (t=0;t<L0;t++)
   {
      local_hsum(isx[t],&S);
      asl0[t+t0]=0.5*S;
   }

   if (NPROC>1)
   {
      MPI_Reduce(asl0,asl,N0,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(asl,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      for (t=0;t<N0;t++)
         asl[t]=asl0[t];
   }

   S=0.0;

   for (t=0;t<N0;t++)
      S+=asl[t];

   return S;
}
