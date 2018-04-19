
/*******************************************************************************
*
* File ltl_modes.c
*
* Copyright (C) 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the little modes.
*
* The externally accessible functions are
*
*   int set_ltl_modes(void)
*     Computes the little modes, the associated little-little Dirac
*     operator and its inverse. The program returns 0 if the inversion
*     was safe and 1 if not.
*
*   complex_dble *ltl_matrix(void)
*     Returns the pointer to an Ns x Ns matrix that represents the
*     *inverse* of the double-precision little-little Dirac operator.
*
* Notes:
*
* For a description of the little Dirac operator and the associated data
* structures see README.Aw. As usual, Ns denotes the number of deflation
* modes in each block of the DFL_BLOCKS grid.
*
* The inversion of a double-precision complex matrix is considered to be
* safe if and only if its Frobenius condition number is less than 10^6.
*
* All programs in this module may involve global communications and must
* be called simultaneously on all MPI processes.
*
*******************************************************************************/

#define LTL_MODES_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "vflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "block.h"
#include "dfl.h"
#include "little.h"
#include "global.h"

#define MAX_FROBENIUS 1.0e6

static int Ns=0,nv,nvh;
static complex **vs;
static complex_dble **vds,*Ads,*Bds,*Cds;


static void sum_vprod(int n,complex_dble *z,complex_dble *w)
{
   int k;

   if (NPROC>1)
   {
      MPI_Reduce(z,w,2*n,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(w,2*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      for (k=0;k<n;k++)
      {
         w[k].re=z[k].re;
         w[k].im=z[k].im;
      }
   }
}


static void alloc_matrices(void)
{
   int k,nmat;
   dfl_parms_t dfl;
   dfl_grid_t grd;

   dfl=dfl_parms();
   grd=dfl_geometry();

   Ns=dfl.Ns;
   nv=Ns*grd.nb;
   nvh=nv/2;
   vs=vflds();
   vds=vdflds();
   nmat=Ns*Ns;

   Ads=amalloc(3*nmat*sizeof(*Ads),ALIGN);

   error(Ads==NULL,1,"alloc_matrices [ltl_modes.c]",
         "Unable to allocate matrices");

   Bds=Ads+nmat;
   Cds=Bds+nmat;
   set_vd2zero(3*nmat,Ads);

   for (k=0;k<Ns;k++)
   {
      Ads[Ns*k+k].re=1.0;
      Bds[Ns*k+k].re=1.0;
   }
}


int set_ltl_modes(void)
{
   int k,l,nmat,ifail;
   double r,cn;
   complex **wv;
   complex_dble **wvd,z;

   if (Ns==0)
      alloc_matrices();

   wvd=reserve_wvd(2);
   wv=vs+Ns;

   for (k=0;k<Ns;k++)
   {
      assign_v2vd(nvh,wv[k],vds[k]);

      if (k>0)
      {
         for (l=0;l<k;l++)
            Cds[l]=vprod_dble(nvh,0,vds[l],vds[k]);

         sum_vprod(k,Cds,Cds+Ns);

         for (l=0;l<k;l++)
         {
            z.re=-Cds[Ns+l].re;
            z.im=-Cds[Ns+l].im;

            mulc_vadd_dble(nv,vds[k],vds[l],z);
         }
      }

      r=vnorm_square_dble(nvh,1,vds[k]);
      r=sqrt(r);

      if (r!=0.0)
      {
         vscale_dble(nvh,1.0/r,vds[k]);
         assign_vd2vd(nvh,vds[k],wvd[0]);
         Awhat_dble(wvd[0],wvd[1]);
         assign_vd2vd(nvh,wvd[1],vds[k]+nvh);
         assign_vd2v(nv,vds[k],vs[k]);
      }
      else
         error_root(1,1,"set_ltl_modes() [ltl_modes.c]",
                    "Degenerate little modes");
   }

   release_wvd();

   for (k=0;k<Ns;k++)
   {
      for (l=0;l<Ns;l++)
         Cds[Ns*k+l]=vprod_dble(nvh,0,vds[k],vds[l]+nvh);
   }

   nmat=Ns*Ns;
   sum_vprod(nmat,Cds,Ads);
   ifail=cmat_inv_dble(Ns,Ads,Bds,&cn);
   if (cn>MAX_FROBENIUS)
      ifail=1;

   return ifail;
}


complex_dble *ltl_matrix(void)
{
   if (Ns==0)
      alloc_matrices();

   return Bds;
}
