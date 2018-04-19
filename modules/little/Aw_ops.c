
/*******************************************************************************
*
* File Aw_ops.c
*
* Copyright (C) 2011, 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation and initialization of the little Dirac operator.
*
* The externally accessible functions are
*
*   Aw_t Awop(void)
*     Returns a structure containing the matrices that describe the
*     single-precision little Dirac operator.
*
*   Aw_t Awophat(void)
*     Returns a structure containing the matrices that describe the
*     single-precision even-odd preconditioned little Dirac operator.
*
*   Aw_dble_t Awop_dble(void)
*     Returns a structure containing the matrices that describe the
*     double-precision little Dirac operator.
*
*   Aw_dble_t Awophat_dble(void)
*     Returns a structure containing the matrices that describe the
*     double-precision even-odd preconditioned little Dirac operator.
*
*   void set_Aw(double mu)
*     Computes the single- and the double-precision little Dirac operator.
*     The SW term is updated if needed and the twisted mass is set to mu.
*     If the twisted-mass flag is set, the twisted-mass term is switched
*     on the odd sites of the lattice.
*
*   int set_Awhat(double mu)
*     Computes the single- and the double-precision even-odd preconditioned
*     little Dirac operator. The program calls set_Aw(mu) and thus updates
*     the operator w/o even-odd preconditioning too. The little modes are
*     updated as well (see ltl_modes.c). On exit the program returns 0 if
*     all matrix inversions were safe and 1 if not.
*
* Notes:
*
* For a description of the little Dirac operator and the associated data
* structures see README.Aw. The twisted-mass flag is retrieved from the
* parameter data base (see flags/lat_parms.c).
*
* The inversion of a double-precision complex matrix is considered to be
* safe if and only if its Frobenius condition number is less than 10^6.
*
* All programs in this module may involve global communications and must
* be called simultaneously on all processes.
*
*******************************************************************************/

#define AW_OPS_C

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
#define MAX_UPDATE 128

static int Ns=0,nb,nbh,nbbh;
static int *idx,(*inn)[8];
static int old_eo[2],nupd=0;
static double old_m0[2],old_mu[2];
static Aw_dble_t Awd={0},Awdhat={0};
static Aw_t Aws={0},Awshat={0};


static void set_constants(void)
{
   dfl_parms_t dfl;
   dfl_grid_t grd;

   dfl=dfl_parms();
   grd=dfl_geometry();

   Ns=dfl.Ns;
   nb=grd.nb;
   nbh=nb/2;
   nbbh=grd.nbb/2;
   idx=grd.idx;
   inn=grd.inn;
}


static void alloc_Awd(Aw_dble_t *Aw)
{
   int n,k,nmat,nbee,nboe;
   complex_dble **ww,*w;

   if (Ns==0)
      set_constants();

   nmat=Ns*Ns;
   nbee=0;
   nboe=0;
   if (Aw==(&Awd))
      nboe=nbbh;
   if (Aw==(&Awdhat))
      nbee=nbbh;
   n=18*nbh+nbee+2*nboe;
   ww=malloc(n*sizeof(*ww));
   w=amalloc(n*nmat*sizeof(*w),ALIGN);
   error((ww==NULL)||(w==NULL),1,"alloc_Awd [Aw_ops.c]",
         "Unable to allocate matrix arrays");

   for (k=0;k<n;k++)
      ww[k]=w+k*nmat;

   set_vd2zero(n*nmat,w);

   for (n=0;n<(nb+nbee);n++)
   {
      for (k=0;k<Ns;k++)
         ww[n][Ns*k+k].re=1.0;
   }

   (*Aw).Ns=Ns;
   (*Aw).nb=nb;
   (*Aw).Aee=ww;
   ww+=nbh+nbee;
   (*Aw).Aoo=ww;
   ww+=nbh;
   (*Aw).Aoe=ww;
   ww+=8*nbh+nboe;
   (*Aw).Aeo=ww;
}


static void alloc_Aws(Aw_t *Aw)
{
   int n,k,nmat;
   complex **ww,*w;

   if (Ns==0)
      set_constants();

   nmat=Ns*Ns;
   n=18*nbh;
   ww=malloc(n*sizeof(*ww));
   w=amalloc(n*nmat*sizeof(*w),ALIGN);
   error((ww==NULL)||(w==NULL),1,"alloc_Aws [Aw_ops.c]",
         "Unable to allocate matrix arrays");

   for (k=0;k<n;k++)
      ww[k]=w+k*nmat;

   set_v2zero(n*nmat,w);

   for (n=0;n<nb;n++)
   {
      for (k=0;k<Ns;k++)
         ww[n][Ns*k+k].re=1.0f;
   }

   (*Aw).Ns=Ns;
   (*Aw).nb=nb;
   (*Aw).Aee=ww;
   ww+=nbh;
   (*Aw).Aoo=ww;
   ww+=nbh;
   (*Aw).Aoe=ww;
   ww+=8*nbh;
   (*Aw).Aeo=ww;
}


Aw_dble_t Awop_dble(void)
{
   if (Awd.Ns==0)
      alloc_Awd(&Awd);

   return Awd;
}


Aw_dble_t Awophat_dble(void)
{
   if (Awdhat.Ns==0)
      alloc_Awd(&Awdhat);

   return Awdhat;
}


Aw_t Awop(void)
{
   if (Aws.Ns==0)
      alloc_Aws(&Aws);

   return Aws;
}


Aw_t Awophat(void)
{
   if (Awshat.Ns==0)
      alloc_Aws(&Awshat);

   return Awshat;
}


static void assign_Awd2Aw(Aw_dble_t *Bwd,Aw_t *Bw)
{
   int n;

   n=nbh*Ns*Ns;

   assign_vd2v(n,(*Bwd).Aee[0],(*Bw).Aee[0]);
   assign_vd2v(n,(*Bwd).Aoo[0],(*Bw).Aoo[0]);
   assign_vd2v(8*n,(*Bwd).Aoe[0],(*Bw).Aoe[0]);
   assign_vd2v(8*n,(*Bwd).Aeo[0],(*Bw).Aeo[0]);
}


static void update_Awdiag(double m0,double mu,int eo)
{
   int nbs,isw,vol,volh;
   int n,nsw,k,l;
   double dm0,dme,dmo;
   complex_dble w,*z;
   spinor_dble **sd;
   block_t *b;

   dm0=m0-old_m0[0];
   dme=mu-old_mu[0];

   if (eo==1)
   {
      if (old_eo[0]==1)
         dmo=0.0;
      else
         dmo=-old_mu[0];
   }
   else
   {
      if (old_eo[0]==1)
         dmo=mu;
      else
         dmo=dme;
   }

   b=blk_list(DFL_BLOCKS,&nbs,&isw);
   vol=(*b).vol;
   volh=vol/2;

   if ((nupd<MAX_UPDATE)&&(fabs(dm0)<1.0)&&(fabs(dme)<1.0)&&(fabs(dmo)<1.0))
   {
      for (n=0;n<nb;n++)
      {
         nsw=idx[n];

         if (nsw<nbh)
            z=Awd.Aee[nsw];
         else
            z=Awd.Aoo[nsw-nbh];

         sd=(*b).sd;

         for (k=0;k<Ns;k++)
         {
            z[Ns*k+k].re+=dm0;

            if (dme!=0.0)
            {
               for (l=k;l<Ns;l++)
               {
                  w=spinor_prod5_dble(volh,0,sd[k+1],sd[l+1]);

                  if (l!=k)
                  {
                     z[Ns*k+l].re-=dme*w.im;
                     z[Ns*k+l].im+=dme*w.re;
                     z[Ns*l+k].re+=dme*w.im;
                     z[Ns*l+k].im+=dme*w.re;
                  }
                  else
                     z[Ns*k+k].im+=dme*w.re;
               }
            }

            if (dmo!=0.0)
            {
               for (l=k;l<Ns;l++)
               {
                  w=spinor_prod5_dble(volh,0,sd[k+1]+volh,sd[l+1]+volh);

                  if (l!=k)
                  {
                     z[Ns*k+l].re-=dmo*w.im;
                     z[Ns*k+l].im+=dmo*w.re;
                     z[Ns*l+k].re+=dmo*w.im;
                     z[Ns*l+k].im+=dmo*w.re;
                  }
                  else
                     z[Ns*k+k].im+=dmo*w.re;
               }
            }
         }

         b+=1;
      }

      nupd+=1;
   }
   else
   {
      sw_term(NO_PTS);

      for (n=0;n<nb;n++)
      {
         assign_ud2udblk(DFL_BLOCKS,n);
         assign_swd2swdblk(DFL_BLOCKS,n,NO_PTS);
         nsw=idx[n];

         if (nsw<nbh)
            z=Awd.Aee[nsw];
         else
            z=Awd.Aoo[nsw-nbh];

         sd=(*b).sd;

         for (l=0;l<Ns;l++)
         {
            Dw_blk_dble(DFL_BLOCKS,n,mu,l+1,0);

            for (k=0;k<Ns;k++)
               z[Ns*k+l]=spinor_prod_dble(vol,0,sd[k+1],sd[0]);
         }

         b+=1;
      }

      nupd=0;
   }

   old_m0[0]=m0;
   old_mu[0]=mu;
   old_eo[0]=eo;
   assign_Awd2Aw(&Awd,&Aws);
}


void set_Aw(double mu)
{
   int n,nsw,msw,nbd,mbd,nu,ifc;
   int eo,nbs,isw,vol,ibn,k,l;
   double dprms[1],m0;
   complex_dble *z,*w,sp[2];
   spinor_dble **sd,**sde,**sdo;
   block_t *b;
   b2b_flds_t *b2b;
   sw_parms_t sw;
   tm_parms_t tm;

   if (NPROC>1)
   {
      dprms[0]=mu;
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      error(dprms[0]!=mu,1,
            "set_Aw [Aw_ops.c]","Parameters are not global");
   }

   if (Awd.Ns==0)
      alloc_Awd(&Awd);

   sw=sw_parms();
   m0=sw.m0;
   tm=tm_parms();
   eo=tm.eoflg;

   if (query_flags(AW_UP2DATE)==1)
   {
      if ((m0!=old_m0[0])||(mu!=old_mu[0])||(eo!=old_eo[0]))
         update_Awdiag(m0,mu,eo);
      return;
   }

   sw_term(NO_PTS);
   b=blk_list(DFL_BLOCKS,&nbs,&isw);

   for (n=0;n<nb;n++)
   {
      assign_ud2udblk(DFL_BLOCKS,n);
      assign_swd2swdblk(DFL_BLOCKS,n,NO_PTS);
      nsw=idx[n];

      if (nsw<nbh)
         z=Awd.Aee[nsw];
      else
         z=Awd.Aoo[nsw-nbh];

      vol=(*b).vol;
      sd=(*b).sd;

      for (l=0;l<Ns;l++)
      {
         Dw_blk_dble(DFL_BLOCKS,n,mu,l+1,0);

         for (k=0;k<Ns;k++)
            z[Ns*k+l]=spinor_prod_dble(vol,0,sd[k+1],sd[0]);
      }

      for (nu=0;nu<4;nu++)
      {
         b2b=b2b_flds(n,nu);
         msw=idx[(*b2b).n[1]];
         vol=(*b2b).vol;
         ibn=(*b2b).ibn;
         ifc=2*nu+1;

         sde=(*b2b).sde[0];
         sdo=(*b2b).sdo[1];

         if (msw>=nbh)
         {
            z=Awd.Aeo[8*(msw-nbh)+(ifc^0x1)];
            w=Awd.Aoe[8*(msw-nbh)+(ifc^0x1)];
         }
         else if (ibn)
         {
            mbd=inn[msw][ifc^0x1]-nb-nbbh;
            z=Awd.Aoe[8*nbh+mbd];
            w=Awd.Aeo[8*nbh+mbd];
         }
         else
         {
            z=Awd.Aoe[8*(nsw-nbh)+ifc];
            w=Awd.Aeo[8*(nsw-nbh)+ifc];
         }

         for (k=0;k<Ns;k++)
         {
            for (l=0;l<Ns;l++)
            {
               spinor_prod_gamma[nu](vol,sde[k],sdo[l],sp);

               z[k*Ns+l].re=(-0.5)*(sp[0].re-sp[1].re);
               z[k*Ns+l].im=(-0.5)*(sp[0].im-sp[1].im);

               w[l*Ns+k].re=(-0.5)*(sp[0].re+sp[1].re);
               w[l*Ns+k].im=( 0.5)*(sp[0].im+sp[1].im);
            }
         }

         sde=(*b2b).sde[1];
         sdo=(*b2b).sdo[0];

         if (nsw>=nbh)
         {
            z=Awd.Aoe[8*(nsw-nbh)+ifc];
            w=Awd.Aeo[8*(nsw-nbh)+ifc];
         }
         else if (ibn)
         {
            nbd=inn[nsw][ifc]-nb-nbbh;
            z=Awd.Aeo[8*nbh+nbd];
            w=Awd.Aoe[8*nbh+nbd];
         }
         else
         {
            z=Awd.Aeo[8*(msw-nbh)+(ifc^0x1)];
            w=Awd.Aoe[8*(msw-nbh)+(ifc^0x1)];
         }

         for (k=0;k<Ns;k++)
         {
            for (l=0;l<Ns;l++)
            {
               spinor_prod_gamma[nu](vol,sdo[k],sde[l],sp);

               if (ibn)
               {
                  z[k*Ns+l].re=(-0.5)*(sp[0].re-sp[1].re);
                  z[k*Ns+l].im=(-0.5)*(sp[0].im-sp[1].im);

                  w[l*Ns+k].re=(-0.5)*(sp[0].re+sp[1].re);
                  w[l*Ns+k].im=( 0.5)*(sp[0].im+sp[1].im);
               }
               else
               {
                  z[k*Ns+l].re+=(-0.5)*(sp[0].re-sp[1].re);
                  z[k*Ns+l].im+=(-0.5)*(sp[0].im-sp[1].im);

                  w[l*Ns+k].re+=(-0.5)*(sp[0].re+sp[1].re);
                  w[l*Ns+k].im+=( 0.5)*(sp[0].im+sp[1].im);
               }
            }
         }
      }

      b+=1;
   }

   cpAoe_ext_bnd();

   if (Aws.Ns==0)
      alloc_Aws(&Aws);
   assign_Awd2Aw(&Awd,&Aws);

   nupd=0;
   old_m0[0]=m0;
   old_mu[0]=mu;
   old_eo[0]=eo;
   set_flags(COMPUTED_AW);
}


int set_Awhat(double mu)
{
   int eo,n,m,ifc,ifail;
   double m0,cn;
   sw_parms_t sw;
   tm_parms_t tm;

   set_Aw(mu);
   sw=sw_parms();
   m0=sw.m0;
   tm=tm_parms();
   eo=tm.eoflg;

   if (query_flags(AWHAT_UP2DATE)==1)
   {
      if ((m0==old_m0[1])&&(mu==old_mu[1])&&(eo==old_eo[1]))
         return 0;
   }

   if (Awdhat.Ns==0)
      alloc_Awd(&Awdhat);

   ifail=0;

   for (n=0;n<nbh;n++)
   {
      ifail|=cmat_inv_dble(Ns,Awd.Aee[n],Awdhat.Aee[n],&cn);
      if (cn>MAX_FROBENIUS)
         ifail=1;
   }

   cpAee_int_bnd();

   for (n=0;n<nbh;n++)
   {
      for (ifc=0;ifc<8;ifc++)
      {
         m=inn[n+nbh][ifc];

         if (m>=nb)
            m-=nbh;

         cmat_mul_dble(Ns,Awdhat.Aee[m],Awd.Aeo[8*n+ifc],
                       Awdhat.Aeo[8*n+ifc]);
      }
   }

   for (n=0;n<nbh;n++)
   {
      ifail|=cmat_inv_dble(Ns,Awd.Aoo[n],Awdhat.Aoo[n],&cn);
      if (cn>MAX_FROBENIUS)
         ifail=1;

      for (ifc=0;ifc<8;ifc++)
         cmat_mul_dble(Ns,Awdhat.Aoo[n],Awd.Aoe[8*n+ifc],Awdhat.Aoe[8*n+ifc]);
   }

   if (Awshat.Ns==0)
      alloc_Aws(&Awshat);
   assign_Awd2Aw(&Awdhat,&Awshat);
   set_flags(COMPUTED_AWHAT);

   old_m0[1]=m0;
   old_mu[1]=mu;
   old_eo[1]=eo;
   ifail|=set_ltl_modes();
   MPI_Allreduce(&ifail,&n,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

   return n;
}
