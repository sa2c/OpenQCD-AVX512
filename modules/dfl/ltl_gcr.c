
/*******************************************************************************
*
* File ltl_gcr.c
*
* Copyright (C) 2011-2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* GCR solver for the little Dirac equation Aw*psi=eta.
*
* The externally accessible functions are
*
*   double ltl_gcr(int nkv,int nmx,double res,double mu,
*                  complex_dble *eta,complex_dble *psi,int *status)
*     Obtains an approximate solution psi of the little Dirac equation for
*     given source eta using the even-odd preconditioned GCR algorithm. See
*     the notes for the explanation of the parameters of the program.
*
* Notes:
*
* The program is based on the flexible GCR algorithm for complex fields (see
* linsolv/fgcr4vd.c). The improvement coefficients, the quark mass in the SW
* term and the parameters of the deflation subspace are assumed to be set by
* set_lat_parms(), set_bc_parms(), set_sw_parms() and set_dfl_parms(). It is
* also taken for granted that the deflation subspace has been initialized by
* calling dfl_subspace().
*
* The parameters passed through the argument list are:
*
*   nkv     Maximal number of Krylov vectors generated before the GCR
*           algorithm is restarted.
*
*   nmx     Maximal total number of Krylov vectors that may be generated.
*
*   res     Desired maximal relative residue |eta-Aw*psi|/|eta| of the
*           calculated solution.
*
*   mu      Value of the twisted mass in the Dirac equation.
*
*   eta     Source field. eta is unchanged on exit if psi does not
*           coincide with eta (which is permissible).
*
*   psi     Calculated approximate solution of the little Dirac equation.
*
*   status  If the program is able to solve the little equation to the
*           desired accuracy, status reports the total number of Krylov
*           vectors that were required for the solution. Negative values
*           indicate that the program failed (-1: the algorithm did not
*           converge, -2: the inversion of the diagonal elements of the
*           little Dirac operator was not safe).
*
* When status>=-1, the program returns the norm of the residue of the
* calculated approximate solution of the even-odd preconditioned, globally
* deflated little Dirac equation. No action is performed if status=-2
* and the program returns 1.0.
*
* The even-odd preconditioned little Dirac operator is updated if it is
* not up-to-date. Evidently the solver is a global program that must be
* called on all processes simultaneously. The required workspaces are
*
*  complex             2*nkv+1
*  complex_dble        3
*
* (see utils/wspace.c).
*
*******************************************************************************/

#define LTL_GCR_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "flags.h"
#include "vflds.h"
#include "linalg.h"
#include "linsolv.h"
#include "little.h"
#include "dfl.h"
#include "global.h"

static int Ns=0,nv,nvh;
static double rvol;
static complex **vs;
static complex_dble **vds,*awd,*cs1,*cs2;


static void set_constants(void)
{
   dfl_parms_t dfl;
   dfl_grid_t grd;

   dfl=dfl_parms();
   grd=dfl_geometry();

   Ns=dfl.Ns;
   nv=Ns*grd.nb;
   nvh=nv/2;
   rvol=1.0/sqrt((double)(nv)*(double)(NPROC));

   vs=vflds();
   vds=vdflds();
   awd=ltl_matrix();

   cs1=amalloc(2*Ns*sizeof(*cs1),ALIGN);
   error(cs1==NULL,1,"set_constants [ltl_gcr.c]",
         "Unable to allocate auxiliary arrays");
   cs2=cs1+Ns;
}


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


static void Lvd(complex_dble *v,complex_dble *w)
{
   int i;
   complex_dble z;

   for (i=0;i<Ns;i++)
      cs1[i]=vprod_dble(nvh,0,vds[i],v);

   sum_vprod(Ns,cs1,cs2);
   cmat_vec_dble(Ns,awd,cs2,cs1);
   set_vd2zero(nvh,w);

   for (i=0;i<Ns;i++)
   {
      mulc_vadd_dble(nvh,w,vds[i],cs1[i]);
      z.re=-cs1[i].re;
      z.im=-cs1[i].im;
      mulc_vadd_dble(nvh,v,vds[i]+nvh,z);
   }
}


static void RLvd(complex_dble *v,complex_dble *w)
{
   int i;
   complex_dble z;

   for (i=0;i<Ns;i++)
      cs1[i]=vprod_dble(nvh,0,vds[i],w);

   sum_vprod(Ns,cs1,cs2);
   cmat_vec_dble(Ns,awd,cs2,cs1);

   for (i=0;i<Ns;i++)
   {
      z.re=-cs1[i].re;
      z.im=-cs1[i].im;
      mulc_vadd_dble(nvh,v,vds[i],z);
      mulc_vadd_dble(nvh,w,vds[i]+nvh,z);
   }
}


static void Lv(complex *v)
{
   int i;
   complex z;

   for (i=0;i<Ns;i++)
   {
      z=vprod(nvh,0,vs[i],v);
      cs1[i].re=(double)(z.re);
      cs1[i].im=(double)(z.im);
   }

   sum_vprod(Ns,cs1,cs2);
   cmat_vec_dble(Ns,awd,cs2,cs1);

   for (i=0;i<Ns;i++)
   {
      z.re=(float)(-cs1[i].re);
      z.im=(float)(-cs1[i].im);
      mulc_vadd(nvh,v,vs[i]+nvh,z);
   }
}


static void Mop(int k,complex *rho,complex *phi,complex *chi)
{
   assign_v2v(nvh,rho,phi);
   Awhat(phi,chi);
   Lv(chi);
}


static void Dop(complex_dble *v,complex_dble *w)
{
   Awhat_dble(v,w);
   RLvd(v,w);
}


double ltl_gcr(int nkv,int nmx,double res,double mu,
               complex_dble *eta,complex_dble *psi,int *status)
{
   int ifail;
   double rho,rho0,fact;
   complex **wv;
   complex_dble **wvd,z;

   if (Ns==0)
      set_constants();

   status[0]=0;
   rho0=sqrt(vnorm_square_dble(nv,1,eta));
   rho=rho0;
   ifail=set_Awhat(mu);

   if (ifail)
      status[0]=-2;
   else
   {
      wv=reserve_wv(2*nkv+1);
      wvd=reserve_wvd(3);

      assign_vd2vd(nvh,eta,wvd[0]);
      Awooinv_dble(eta,wvd[0]);
      Aweo_dble(wvd[0],wvd[0]);
      Aweeinv_dble(wvd[0],wvd[1]);
      fact=rvol*rho0;

      if (fact!=0.0)
      {
         vscale_dble(nvh,1.0/fact,wvd[1]);
         Lvd(wvd[1],wvd[2]);

         rho=fgcr4vd(nvh,1,Dop,Mop,wv,wvd,nkv,nmx,res,wvd[1],psi,status);

         z.re=1.0;
         z.im=0.0;
         mulc_vadd_dble(nvh,psi,wvd[2],z);
         vscale_dble(nvh,fact,psi);
         rho*=fact;
      }
      else
      {
         rho=0.0;
         set_vd2zero(nv,psi);
      }

      Awoe_dble(psi,wvd[0]);
      assign_vd2vd(nvh,eta+nvh,wvd[1]+nvh);
      z.re=-1.0;
      z.im=0.0;
      mulc_vadd_dble(nvh,wvd[1]+nvh,wvd[0]+nvh,z);
      Awooinv_dble(wvd[1],psi);

      release_wvd();
      release_wv();
   }

   if (status[0]<-1)
   {
      rho=rho0;
      set_vd2zero(nv,psi);
   }

   return rho;
}
