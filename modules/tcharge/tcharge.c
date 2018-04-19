
/*******************************************************************************
*
* File tcharge.c
*
* Copyright (C) 2010-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the topological charge using the symmetric field tensor.
*
* The externally accessible functions are
*
*   double tcharge(void)
*     Returns the "field-theoretic" topological charge Q of the global
*     double-precision gauge field, using a symmetric expression for the
*     gauge-field tensor.
*
*   double tcharge_slices(double *qsl)
*     Computes the sum qsl[x0] of the "field-theoretic" topological charge
*     density of the double-precision gauge field at time x0=0,1,...,N0-1
*     (where N0=NPROC0*L0). The program returns the total charge.
*
* Notes:
*
* The topological charge density q(x) is defined by
*
*  q(x)=(8*Pi^2)^(-1)*{F_{01}^a(x)*F_{23}^a(x)+
*                      F_{02}^a(x)*F_{31}^a(x)+
*                      F_{03}^a(x)*F_{12}^a(x)},
*
* where
*
*  F_{mu,nu}^a(x)=-2*tr{F_{mu,nu}(x)*T^a}, a=1,..,8,
*
* are the SU(3) components of the symmetric field tensor returned by the
* program ftensor() [ftensor.c]. At the boundaries of the lattice (if any),
* the charge density is set to zero. The total charge Q is the sum of q(x)
* over all points x with time component in the range
*
*  0<x0<NPROC0*L0-1        (open bc),
*
*  0<x0<NPROC0*L0          (SF and open-SF bc),
*
*  0<=x0<NPROC0*L0         (periodic bc).
*
* The programs in this module perform global communications and must be
* called simultaneously on all processes.
*
*******************************************************************************/

#define TCHARGE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "tcharge.h"
#include "global.h"

#define N0 (NPROC0*L0)

static int isx[L0],init=0;
static double qsl0[N0];
static u3_alg_dble **ft;


static double prodXY(u3_alg_dble *X,u3_alg_dble *Y)
{
   double sm;

   sm=(-2.0/3.0)*((*X).c1+(*X).c2+(*X).c3)*((*Y).c1+(*Y).c2+(*Y).c3)+
      2.0*((*X).c1*(*Y).c1+(*X).c2*(*Y).c2+(*X).c3*(*Y).c3)+
      4.0*((*X).c4*(*Y).c4+(*X).c5*(*Y).c5+(*X).c6*(*Y).c6+
           (*X).c7*(*Y).c7+(*X).c8*(*Y).c8+(*X).c9*(*Y).c9);

   return sm;
}


static double density(int ix)
{
   double sm;

   sm=prodXY(ft[0]+ix,ft[3]+ix)+
      prodXY(ft[1]+ix,ft[4]+ix)+
      prodXY(ft[2]+ix,ft[5]+ix);

   return sm;
}


double tcharge(void)
{
   int bc,ix,t,tmx;
   double pi,Q;

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
         Q=density(ix);
         add_to_hsum(isx[0],&Q);
      }
   }

   if (NPROC>1)
      global_hsum(isx[0],&Q);
   else
      local_hsum(isx[0],&Q);

   pi=4.0*atan(1.0);

   return Q/(8.0*pi*pi);
}


double tcharge_slices(double *qsl)
{
   int bc,ix,t,t0,tmx;
   double pi,fact,Q;

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
         Q=density(ix);
         add_to_hsum(isx[t],&Q);
      }
   }

   for (t=0;t<N0;t++)
      qsl0[t]=0.0;

   pi=4.0*atan(1.0);
   fact=1.0/(8.0*pi*pi);

   for (t=0;t<L0;t++)
   {
      local_hsum(isx[t],&Q);
      qsl0[t+t0]=fact*Q;
   }

   if (NPROC>1)
   {
      MPI_Reduce(qsl0,qsl,N0,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(qsl,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      for (t=0;t<N0;t++)
         qsl[t]=qsl0[t];
   }

   Q=0.0;

   for (t=0;t<N0;t++)
      Q+=qsl[t];

   return Q;
}
