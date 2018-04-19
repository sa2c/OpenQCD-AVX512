
/*******************************************************************************
*
* File sw_term.c
*
* Copyright (C) 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the SW term.
*
* The externally accessible functions are
*
*   int sw_term(ptset_t set)
*     Computes the SW term for the current double-precision gauge field
*     and assigns the matrix to the global double-precision SW field. The
*     matrices on the specified point set are then inverted and 0 or 1
*     is returned depending on whether all inversions were safe or not.
*
* Notes:
*
* The program sets the SW term to unity at global time
*
*  x0=0                (open, SF and open-SF boundary conditions),
*
*  x0=NPROC0*L0-1      (open boundary conditions).
*
* In all other cases, it is given by
*
*    c(x0)+csw*(i/4)*sigma_{mu nu}*Fhat_{mu nu}(x)
*
* where
*
*    c(x0) = 4+m0+cF[0]-1     if x0=1 (open, SF or open-SF bc),
*            4+m0+cF[1]-1     if x0=NPROCO*L0-2 (open bc),
*                             or x0=NPROC0*L0-1 (SF or open-SF bc),
*            4+m0             otherwise,
*
*    sigma_{mu nu}=(i/2)*[gamma_mu,gamma_nu],
*
* and Fhat_{mu nu} is the standard (clover) expression for the gauge field
* tensor as computed by the program ftensor() [tcharge/ftensor.c]. The upper
* and lower 6x6 blocks of the matrix are stored in the pauli_dble structures
* swd[2*ix] and swd[2*ix+1], where ix is the label of the point x.
*
* The quark mass m0 and the improvement coefficients csw and cF are obtained
* from the parameter data base by calling sw_parms() [flags/lat_parms.c]. Note
* that this program checks the flags data base and only computes those parts
* of the SW array that do not already have the correct values.
*
* This program performs global operations and must be called simultaneously
* on all processes.
*
*******************************************************************************/

#define SW_TERM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "tcharge.h"
#include "sw_term.h"
#include "global.h"

#define N0 (NPROC0*L0)

static double c1,c2,c3[2];
static u3_alg_dble X;
static const pauli_dble sw0={{1.0,1.0,1.0,1.0,1.0,1.0,
                              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                              0.0,0.0,0.0,0.0,0.0,0.0}};


static void u3_alg2pauli1(pauli_dble *m)
{
   (*m).u[10]=-X.c1;

   (*m).u[12]=-X.c5;
   (*m).u[13]= X.c4;
   (*m).u[14]=-X.c7;
   (*m).u[15]= X.c6;

   (*m).u[18]=-X.c5;
   (*m).u[19]=-X.c4;
   (*m).u[20]=-X.c2;

   (*m).u[22]=-X.c9;
   (*m).u[23]= X.c8;
   (*m).u[24]=-X.c7;
   (*m).u[25]=-X.c6;
   (*m).u[26]=-X.c9;
   (*m).u[27]=-X.c8;
   (*m).u[28]=-X.c3;
}


static void u3_alg2pauli2(pauli_dble *m)
{
   (*m).u[11] =X.c1;
   (*m).u[12]+=X.c4;
   (*m).u[13]+=X.c5;
   (*m).u[14]+=X.c6;
   (*m).u[15]+=X.c7;

   (*m).u[18]-=X.c4;
   (*m).u[19]+=X.c5;

   (*m).u[21] =X.c2;
   (*m).u[22]+=X.c8;
   (*m).u[23]+=X.c9;
   (*m).u[24]-=X.c6;
   (*m).u[25]+=X.c7;
   (*m).u[26]-=X.c8;
   (*m).u[27]+=X.c9;

   (*m).u[29] =X.c3;
}


static void u3_alg2pauli3(pauli_dble *m)
{
   (*m).u[ 0]=-X.c1;
   (*m).u[ 1]=-X.c2;
   (*m).u[ 2]=-X.c3;
   (*m).u[ 3]= X.c1;
   (*m).u[ 4]= X.c2;
   (*m).u[ 5]= X.c3;
   (*m).u[ 6]=-X.c5;
   (*m).u[ 7]= X.c4;
   (*m).u[ 8]=-X.c7;
   (*m).u[ 9]= X.c6;

   (*m).u[16]=-X.c9;
   (*m).u[17]= X.c8;

   (*m).u[30]= X.c5;
   (*m).u[31]=-X.c4;
   (*m).u[32]= X.c7;
   (*m).u[33]=-X.c6;
   (*m).u[34]= X.c9;
   (*m).u[35]=-X.c8;
}


static void set_swd(int vol,int ofs,u3_alg_dble **ft,pauli_dble *sw)
{
   int bc,ix,t;
   double c,*u;
   u3_alg_dble *ft0,*ft1,*ft2,*ft3,*ft4,*ft5;

   bc=bc_type();
   vol+=ofs;
   sw+=2*ofs;
   ft0=ft[0]+ofs;
   ft1=ft[1]+ofs;
   ft2=ft[2]+ofs;
   ft3=ft[3]+ofs;
   ft4=ft[4]+ofs;
   ft5=ft[5]+ofs;

   for (ix=ofs;ix<vol;ix++)
   {
      t=global_time(ix);

      if (((t==0)&&(bc!=3))||((t==(N0-1))&&(bc==0)))
      {
         sw[0]=sw0;
         sw[1]=sw0;
         sw+=2;
      }
      else
      {
         if ((t==1)&&(bc!=3))
            c=c3[0];
         else if (((t==(N0-2))&&(bc==0))||((t==(N0-1))&&((bc==1)||(bc==2))))
            c=c3[1];
         else
            c=c1;

         _u3_alg_mul_sub(X,c2,*ft3,*ft0);
         u3_alg2pauli1(sw);
         _u3_alg_mul_sub(X,c2,*ft4,*ft1);
         u3_alg2pauli2(sw);
         _u3_alg_mul_sub(X,c2,*ft5,*ft2);
         u3_alg2pauli3(sw);

         u=(*sw).u;
         u[0]+=c;
         u[1]+=c;
         u[2]+=c;
         u[3]+=c;
         u[4]+=c;
         u[5]+=c;
         sw+=1;

         _u3_alg_mul_add(X,c2,*ft3,*ft0);
         u3_alg2pauli1(sw);
         _u3_alg_mul_add(X,c2,*ft4,*ft1);
         u3_alg2pauli2(sw);
         _u3_alg_mul_add(X,c2,*ft5,*ft2);
         u3_alg2pauli3(sw);

         u=(*sw).u;
         u[0]+=c;
         u[1]+=c;
         u[2]+=c;
         u[3]+=c;
         u[4]+=c;
         u[5]+=c;
         sw+=1;
      }

      ft0+=1;
      ft1+=1;
      ft2+=1;
      ft3+=1;
      ft4+=1;
      ft5+=1;
   }
}


static int iswd(int vol,pauli_dble *sw)
{
   int ifail,n;
   pauli_dble *sm;

   ifail=0;
   sm=sw+vol;

   for (;sw<sm;sw++)
      ifail|=inv_pauli_dble(0.0,sw,sw);

   if (NPROC>1)
      MPI_Allreduce(&ifail,&n,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
   else
      n=ifail;

   return n;
}


int sw_term(ptset_t set)
{
   int iprms[1],ie,io,ifail;
   pauli_dble *sw;
   u3_alg_dble **ft;
   sw_parms_t swp;

   if (NPROC>1)
   {
      iprms[0]=(int)(set);
      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error(iprms[0]!=(int)(set),1,"sw_term [sw_term.c]",
            "Parameter is not global");
   }

   swp=sw_parms();
   c1=4.0+swp.m0;
   c2=-0.5*swp.csw;
   c3[0]=c1+swp.cF[0]-1.0;
   c3[1]=c1+swp.cF[1]-1.0;

   if (query_flags(SWD_UP2DATE)!=1)
   {
      ft=ftensor();
      sw=swdfld();
      set_swd(VOLUME,0,ft,sw);
      set_flags(COMPUTED_SWD);
   }

   ie=query_flags(SWD_E_INVERTED);
   io=query_flags(SWD_O_INVERTED);

   if ((ie==1)&&((set==NO_PTS)||(set==ODD_PTS)))
   {
      ft=ftensor();
      sw=swdfld();
      set_swd(VOLUME/2,0,ft,sw);
      ie=0;
   }

   if ((io==1)&&((set==NO_PTS)||(set==EVEN_PTS)))
   {
      ft=ftensor();
      sw=swdfld();
      set_swd(VOLUME/2,VOLUME/2,ft,sw);
      io=0;
   }

   ifail=0;

   if ((ie==0)&&((set==ALL_PTS)||(set==EVEN_PTS)))
   {
      sw=swdfld();
      ifail|=iswd(VOLUME,sw);
      ie=1;
   }

   if ((io==0)&&((set==ALL_PTS)||(set==ODD_PTS)))
   {
      sw=swdfld()+VOLUME;
      ifail|=iswd(VOLUME,sw);
      io=1;
   }

   set_flags(COMPUTED_SWD);

   if (ie==1)
      set_flags(INVERTED_SWD_E);

   if (io==1)
      set_flags(INVERTED_SWD_O);

   return ifail;
}
