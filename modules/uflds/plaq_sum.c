
/*******************************************************************************
*
* File plaq_sum.c
*
* Copyright (C) 2005, 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Calculation of plaquette sums.
*
* The externally accessible functions are
*
*   double plaq_sum_dble(int icom)
*     Returns the sum of Re[tr{U(p)}] over all unoriented plaquettes p,
*     where U(p) is the product of the double-precision link variables
*     around p. If icom=1 the global sum of the local sums is returned
*     and otherwise just the local sum.
*
*   double plaq_wsum_dble(int icom)
*     Same as plaq_sum_dble(), but giving weight 1/2 to the contribution
*     of the space-like plaquettes at the boundaries of the lattice if
*     boundary conditions of type 0,1 or 2 are chosen.
*
*   double plaq_action_slices(double *asl)
*     Computes the time-slice sums asl[x0] of the tree-level O(a)-improved
*     plaquette action density of the double-precision gauge field. The
*     factor 1/g0^2 is omitted and the time x0 runs from 0 to NPROC0*L0-1.
*     The program returns the total action.
*
* Notes:
*
* The Wilson plaquette action density is defined so that it converges to the
* Yang-Mills action in the classical continuum limit with a rate proportional
* to a^2. In particular, at the boundaries of the lattice (if there are any),
* the space-like plaquettes are given the weight 1/2 and the contribution of
* a plaquette p in the bulk is 2*Re[tr{1-U(p)}].
*
* The time-slice sum asl[x0] computed by plaq_action_slices() includes the
* full contribution to the action of the space-like plaquettes at time x0 and
* 1/2 of the contribution of the time-like plaquettes at time x0 and x0-1.
*
* The programs in this module perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define PLAQ_SUM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "flags.h"
#include "su3fcts.h"
#include "lattice.h"
#include "uflds.h"
#include "global.h"

#define N0 (NPROC0*L0)

static int isA,isE[L0],isB[L0],init=0;
static double aslE[N0],aslB[N0];
static su3_dble *udb;
static su3_dble wd1 ALIGNED16;
static su3_dble wd2 ALIGNED16;


static double plaq_dble(int n,int ix)
{
   int ip[4];
   double sm;

   plaq_uidx(n,ix,ip);

   su3xsu3(udb+ip[0],udb+ip[1],&wd1);
   su3dagxsu3dag(udb+ip[3],udb+ip[2],&wd2);
   cm3x3_retr(&wd1,&wd2,&sm);

   return sm;
}


static void local_plaq_sum_dble(int iw)
{
   int bc,ix,t,n;
   double wp,pa;

   if (init<1)
   {
      isA=init_hsum(1);
      init=1;
   }

   bc=bc_type();

   if (iw==0)
      wp=1.0;
   else
      wp=0.5;

   udb=udfld();
   reset_hsum(isA);

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);
      pa=0.0;

      if ((t<(N0-1))||(bc!=0))
      {
         for (n=0;n<3;n++)
            pa+=plaq_dble(n,ix);
      }

      if (((t>0)||(bc==3))&&((t<(N0-1))||(bc!=0)))
      {
         for (n=3;n<6;n++)
            pa+=plaq_dble(n,ix);
      }
      else
      {
         for (n=3;n<6;n++)
            pa+=wp*plaq_dble(n,ix);
      }

      if ((t==(N0-1))&&((bc==1)||(bc==2)))
         pa+=9.0*wp;

      if (pa!=0.0)
         add_to_hsum(isA,&pa);
   }
}


double plaq_sum_dble(int icom)
{
   double p;

   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();

   local_plaq_sum_dble(0);

   if ((icom==1)&&(NPROC>1))
      global_hsum(isA,&p);
   else
      local_hsum(isA,&p);

   return p;
}


double plaq_wsum_dble(int icom)
{
   double p;

   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();

   local_plaq_sum_dble(1);

   if ((icom==1)&&(NPROC>1))
      global_hsum(isA,&p);
   else
      local_hsum(isA,&p);

   return p;
}


double plaq_action_slices(double *asl)
{
   int bc,ix,t,t0,n;
   double A,smE,smB;

   if (init<2)
   {
      if (init<1)
         isA=init_hsum(1);

      for (t=0;t<L0;t++)
      {
         isE[t]=init_hsum(1);
         isB[t]=init_hsum(1);
      }

      init=2;
   }

   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();

   bc=bc_type();
   t0=cpr[0]*L0;
   udb=udfld();

   for (t=0;t<L0;t++)
   {
      reset_hsum(isE[t]);
      reset_hsum(isB[t]);
   }

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);
      smE=0.0;
      smB=0.0;

      if ((t<(N0-1))||(bc!=0))
      {
         for (n=0;n<3;n++)
            smE+=(3.0-plaq_dble(n,ix));
      }

      if ((t>0)||(bc!=1))
      {
         for (n=3;n<6;n++)
            smB+=(3.0-plaq_dble(n,ix));
      }

      t-=t0;

      if (smE!=0.0)
         add_to_hsum(isE[t],&smE);
      if (smB!=0.0)
         add_to_hsum(isB[t],&smB);
   }

   for (t=0;t<N0;t++)
      asl[t]=0.0;

   for (t=0;t<L0;t++)
   {
      local_hsum(isE[t],&smE);
      asl[t+t0]=smE;
   }

   if (NPROC>1)
   {
      MPI_Reduce(asl,aslE,N0,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(aslE,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      for (t=0;t<N0;t++)
         aslE[t]=asl[t];
   }

   for (t=0;t<N0;t++)
      asl[t]=0.0;

   for (t=0;t<L0;t++)
   {
      local_hsum(isB[t],&smB);
      asl[t+t0]=smB;
   }

   if (NPROC>1)
   {
      MPI_Reduce(asl,aslB,N0,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(aslB,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      for (t=0;t<N0;t++)
         aslB[t]=asl[t];
   }

   if (bc!=3)
      asl[0]=aslE[0]+aslB[0];
   else
      asl[0]=aslE[0]+aslE[N0-1]+2.0*aslB[0];

   if (bc==0)
   {
      for (t=1;t<(N0-1);t++)
         asl[t]=aslE[t-1]+aslE[t]+2.0*aslB[t];

      asl[N0-1]=aslE[N0-2]+aslB[N0-1];
   }
   else
   {
      for (t=1;t<N0;t++)
         asl[t]=aslE[t-1]+aslE[t]+2.0*aslB[t];
   }

   reset_hsum(isA);

   if ((bc==1)||(bc==2))
      add_to_hsum(isA,aslE+N0-1);

   for (t=0;t<N0;t++)
      add_to_hsum(isA,asl+t);

   local_hsum(isA,&A);

   return A;
}
