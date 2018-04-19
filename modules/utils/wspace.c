
/*******************************************************************************
*
* File wspace.c
*
* Copyright (C) 2010, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Workspace allocation.
*
* The externally accessible functions are
*
*   void alloc_wud(int n)
*     Allocates a workspace of n double-precision gauge fields.
*
*   su3_dble **reserve_wud(int n)
*     Reserves a new workspace of n global double-precision gauge fields
*     and returns the array ud[0],..,ud[n-1] of the base addresses of the
*     fields in the workspace. No workspace is reserved and a NULL pointer
*     is returned if n<=0.
*
*   int release_wud(void)
*     Releases the workspace of global double-precision gauge fields that
*     was last reserved and returns the number of fields that are released.
*
*   int wud_size(void)
*     Returns the number of global double-precision gauge fields that
*     are currently reserved.
*
*   void alloc_wfd(int n)
*     Allocates a workspace of n double-precision force fields.
*
*   su3_alg_dble **reserve_wfd(int n)
*     Reserves a new workspace of n global double-precision force fields
*     and returns the array fd[0],..,fd[n-1] of the base addresses of the
*     fields in the workspace. No workspace is reserved and a NULL pointer
*     is returned if n<=0.
*
*   int release_wfd(void)
*     Releases the workspace of global double-precision force fields that
*     was last reserved and returns the number of fields that are released.
*
*   int wfd_size(void)
*     Returns the number of global double-precision force fields that
*     are currently reserved.
*
*   void alloc_ws(int n)
*     Allocates a workspace of n single-precision spinor fields.
*
*   spinor **reserve_ws(int n)
*     Reserves a new workspace of n global single-precision spinor fields
*     and returns the array s[0],..,s[n-1] of the base addresses of the
*     fields in the workspace. No workspace is reserved and a NULL pointer
*     is returned if n<=0.
*
*   int release_ws(void)
*     Releases the workspace of global single-precision spinor fields that
*     was last reserved and returns the number of fields that are released.
*
*   int ws_size(void)
*     Returns the number of global single-precision spinor fields that
*     are currently reserved.
*
*   void alloc_wsd(int n)
*     Allocates a workspace of n double-precision spinor fields.
*
*   spinor_dble **reserve_wsd(int n)
*     Reserves a new workspace of n global double-precision spinor fields
*     and returns the array sd[0],..,sd[n-1] of the base addresses of the
*     fields in the workspace. No workspace is reserved and a NULL pointer
*     is returned if n<=0.
*
*   int release_wsd(void)
*     Releases the workspace of global double-precision spinor fields that
*     was last reserved and returns the number of fields that are released.
*
*   int wsd_size(void)
*     Returns the number of global double-precision spinor fields that
*     are currently reserved.
*
*   void alloc_wv(int n)
*     Allocates a workspace of n single-precision vector fields.
*
*   complex **reserve_wv(int n)
*     Reserves a new workspace of n global single-precision vector fields
*     and returns the array v[0],..,v[n-1] of the base addresses of the
*     fields in the workspace. No workspace is reserved and a NULL pointer
*     is returned if n<=0.
*
*   int release_wv(void)
*     Releases the workspace of global single-precision vector fields that
*     was last reserved and returns the number of fields that are released.
*
*   int wv_size(void)
*     Returns the number of global single-precision vector fields that
*     are currently reserved.
*
*   void alloc_wvd(int n)
*     Allocates a workspace of n double-precision vector fields.
*
*   complex_dble **reserve_wvd(int n)
*     Reserves a new workspace of n global double-precision vector fields
*     and returns the array vd[0],..,vd[n-1] of the base addresses of the
*     fields in the workspace. No workspace is reserved and a NULL pointer
*     is returned if n<=0.
*
*   int release_wvd(void)
*     Releases the workspace of global double-precision vector fields that
*     was last reserved and returns the number of fields that are released.
*
*   int wvd_size(void)
*     Returns the number of global double-precision vector fields that
*     are currently reserved.
*
* Notes:
*
* By definition a workspace is a set of global fields which is reserved for
* the current program. This module allows to assign and release the required
* workspaces dynamically in such a way that the privacy of the workspaces is
* guaranteed.
*
* The total workspace is the set of all allocated fields of a given type. All
* programs must reserve the required workspaces and must release them before
* the program returns to the calling program. It is possible to reserve several
* workspaces of the same type in a given program. In this case, a corresponding
* number of calls to w*_release() is required to release all the workspaces
* that were reserved.
*
* The fields in the gauge and force workspaces have 4*VOLUME elements and can
* only store a copy of the field variables on the local lattice. In the case
* of the spinor and spinor_dble fields, the number of spinors allocated per
* field is NSPIN (see include/global.h). Vector fields are arrays of complex
* numbers on which the little Dirac operator acts. They have Ns*(nb+nbb/2)
* elements, where Ns is the number of local deflation modes, while nb and nbb
* are the numbers blocks in the DFL_BLOCKS grid and its exterior boundary.
*
* Workspaces can always be freed and reallocated by calling alloc_w*() if no
* fields in the workspace are in use. In particular, a workspace can be freed
* by calling alloc_w*(0). Except for the w*_size() programs, all programs in
* this module must be called simultaneously on all MPI processes.
*
*******************************************************************************/

#define WSPACE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "global.h"

static int nudt=0,iwud=0,nwudt=0,*nwud;
static su3_dble **wud0,**wud,ud0={{0.0}};

static int nfdt=0,iwfd=0,nwfdt=0,*nwfd;
static su3_alg_dble **wfd0,**wfd,fd0={0.0};

static int nst=0,iws=0,nwst=0,*nws;
static spinor **ws0,**ws,s0={{{0.0f}}};

static int nsdt=0,iwsd=0,nwsdt=0,*nwsd;
static spinor_dble **wsd0,**wsd,sd0={{{0.0}}};

static int nvt=0,iwv=0,nwvt=0,*nwv,nvec;
static complex **wv0,**wv,v0={0.0f};

static int nvdt=0,iwvd=0,nwvdt=0,*nwvd;
static complex_dble **wvd0,**wvd,vd0={0.0};


void alloc_wud(int n)
{
   int i;
   su3_dble *ud,*um;

   if (n==nudt)
      return;

   error_root(nwudt!=0,1,"alloc_wud [wspace.c]","Fields are in use");

   if (nudt>0)
   {
      free(nwud);
      afree(wud0[0]);
      free(wud0);
      nwud=NULL;
      wud0=NULL;
      wud=NULL;
   }

   nudt=n;
   iwud=0;
   nwudt=0;

   if (nudt>0)
   {
      nwud=malloc(nudt*sizeof(*nwud));
      wud0=malloc(2*nudt*sizeof(*wud0));
      wud=wud0+nudt;

      error((nwud==NULL)||(wud0==NULL),1,"alloc_wud [wspace.c]",
            "Unable to allocate index arrays");

      wud0[0]=amalloc(nudt*4*VOLUME*sizeof(**wud0),ALIGN);

      error(wud0[0]==NULL,1,"alloc_wud [wspace.c]",
            "Unable to allocate workspace");

      for (i=0;i<nudt;i++)
      {
         if (i>0)
            wud0[i]=wud0[i-1]+4*VOLUME;

         nwud[i]=0;
         wud[i]=NULL;
      }

      ud=wud0[0];
      um=ud+nudt*4*VOLUME;

      for (;ud<um;ud++)
         (*ud)=ud0;
   }
}


su3_dble **reserve_wud(int n)
{
   int iprms[1],i,ia;

   if (NPROC>1)
   {
      iprms[0]=n;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error(iprms[0]!=n,1,"reserve_wud [wspace.c]",
            "Parameter n is not global");
   }

   if (n>0)
   {
      error((nwudt+n)>nudt,1,"reserve_wud [wspace.c]",
            "Requested too many fields (tot=%d,use=%d,req=%d)",nudt,nwudt,n);

      ia=nwudt;
      nwud[iwud]=n;
      nwudt+=n;
      iwud+=1;

      for (i=ia;i<(ia+n);i++)
         wud[i]=wud0[i];

      return wud+ia;
   }
   else
      return NULL;
}


int release_wud(void)
{
   int n,i;

   if (nwudt==0)
      return 0;
   else
   {
      iwud-=1;
      n=nwud[iwud];
      nwudt-=n;
      nwud[iwud]=0;

      for (i=nwudt;i<(nwudt+n);i++)
         wud[i]=NULL;

      return n;
   }
}


int wud_size(void)
{
   return nwudt;
}


void alloc_wfd(int n)
{
   int i;
   su3_alg_dble *fd,*fm;

   if (n==nfdt)
      return;

   error_root(nwfdt!=0,1,"alloc_wfd [wspace.c]","Fields are in use");

   if (nfdt>0)
   {
      free(nwfd);
      afree(wfd0[0]);
      free(wfd0);
      nwfd=NULL;
      wfd0=NULL;
      wfd=NULL;
   }

   nfdt=n;
   iwfd=0;
   nwfdt=0;

   if (nfdt>0)
   {
      nwfd=malloc(nfdt*sizeof(*nwfd));
      wfd0=malloc(2*nfdt*sizeof(*wfd0));
      wfd=wfd0+nfdt;

      error((nwfd==NULL)||(wfd0==NULL),1,"alloc_wfd [wspace.c]",
            "Unable to allocate index arrays");

      wfd0[0]=amalloc(nfdt*4*VOLUME*sizeof(**wfd0),ALIGN);

      error(wfd0[0]==NULL,1,"alloc_wfd [wspace.c]",
            "Unable to allocate workspace");

      for (i=0;i<nfdt;i++)
      {
         if (i>0)
            wfd0[i]=wfd0[i-1]+4*VOLUME;

         nwfd[i]=0;
         wfd[i]=NULL;
      }

      fd=wfd0[0];
      fm=fd+nfdt*4*VOLUME;

      for (;fd<fm;fd++)
         (*fd)=fd0;
   }
}


su3_alg_dble **reserve_wfd(int n)
{
   int iprms[1],i,ia;

   if (NPROC>1)
   {
      iprms[0]=n;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error(iprms[0]!=n,1,"reserve_wfd [wspace.c]",
            "Parameter n is not global");
   }

   if (n>0)
   {
      error((nwfdt+n)>nfdt,1,"reserve_wfd [wspace.c]",
            "Requested too many fields (tot=%d,use=%d,req=%d)",nfdt,nwfdt,n);

      ia=nwfdt;
      nwfd[iwfd]=n;
      nwfdt+=n;
      iwfd+=1;

      for (i=ia;i<(ia+n);i++)
         wfd[i]=wfd0[i];

      return wfd+ia;
   }
   else
      return NULL;
}


int release_wfd(void)
{
   int n,i;

   if (nwfdt==0)
      return 0;
   else
   {
      iwfd-=1;
      n=nwfd[iwfd];
      nwfdt-=n;
      nwfd[iwfd]=0;

      for (i=nwfdt;i<(nwfdt+n);i++)
         wfd[i]=NULL;

      return n;
   }
}


int wfd_size(void)
{
   return nwfdt;
}


void alloc_ws(int n)
{
   int i;
   spinor *s,*sm;

   if (n==nst)
      return;

   error_root(nwst!=0,1,"alloc_ws [wspace.c]","Fields are in use");

   if (nst>0)
   {
      free(nws);
      afree(ws0[0]);
      free(ws0);
      nws=NULL;
      ws0=NULL;
      ws=NULL;
   }

   nst=n;
   iws=0;
   nwst=0;

   if (nst>0)
   {
      nws=malloc(nst*sizeof(*nws));
      ws0=malloc(2*nst*sizeof(*ws0));
      ws=ws0+nst;

      error((nws==NULL)||(ws0==NULL),1,"alloc_ws [wspace.c]",
            "Unable to allocate index arrays");

      ws0[0]=amalloc(nst*NSPIN*sizeof(**ws0),ALIGN);

      error(ws0[0]==NULL,1,"alloc_ws [wspace.c]",
            "Unable to allocate workspace");

      for (i=0;i<nst;i++)
      {
         if (i>0)
            ws0[i]=ws0[i-1]+NSPIN;

         nws[i]=0;
         ws[i]=NULL;
      }

      s=ws0[0];
      sm=s+nst*NSPIN;

      for (;s<sm;s++)
         (*s)=s0;
   }
}


spinor **reserve_ws(int n)
{
   int iprms[1],i,ia;

   if (NPROC>1)
   {
      iprms[0]=n;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error(iprms[0]!=n,1,"reserve_ws [wspace.c]",
            "Parameter n is not global");
   }

   if (n>0)
   {
      error((nwst+n)>nst,1,"reserve_ws [wspace.c]",
            "Requested too many fields (tot=%d,use=%d,req=%d)",nst,nwst,n);

      ia=nwst;
      nws[iws]=n;
      nwst+=n;
      iws+=1;

      for (i=ia;i<(ia+n);i++)
         ws[i]=ws0[i];

      return ws+ia;
   }
   else
      return NULL;
}


int release_ws(void)
{
   int n,i;

   if (nwst==0)
      return 0;
   else
   {
      iws-=1;
      n=nws[iws];
      nwst-=n;
      nws[iws]=0;

      for (i=nwst;i<(nwst+n);i++)
         ws[i]=NULL;

      return n;
   }
}


int ws_size(void)
{
   return nwst;
}


void alloc_wsd(int n)
{
   int i;
   spinor_dble *sd,*sm;

   if (n==nsdt)
      return;

   error_root(nwsdt!=0,1,"alloc_wsd [wspace.c]","Fields are in use");

   if (nsdt>0)
   {
      free(nwsd);
      afree(wsd0[0]);
      free(wsd0);
      nwsd=NULL;
      wsd0=NULL;
      wsd=NULL;
   }

   nsdt=n;
   iwsd=0;
   nwsdt=0;

   if (nsdt>0)
   {
      nwsd=malloc(nsdt*sizeof(*nwsd));
      wsd0=malloc(2*nsdt*sizeof(*wsd0));
      wsd=wsd0+nsdt;

      error((nwsd==NULL)||(wsd0==NULL),1,"alloc_wsd [wspace.c]",
            "Unable to allocate index arrays");

      wsd0[0]=amalloc(nsdt*NSPIN*sizeof(**wsd0),ALIGN);

      error(wsd0[0]==NULL,1,"alloc_wsd [wspace.c]",
            "Unable to allocate workspace");

      for (i=0;i<nsdt;i++)
      {
         if (i>0)
            wsd0[i]=wsd0[i-1]+NSPIN;

         nwsd[i]=0;
         wsd[i]=NULL;
      }

      sd=wsd0[0];
      sm=sd+nsdt*NSPIN;

      for (;sd<sm;sd++)
         (*sd)=sd0;
   }
}


spinor_dble **reserve_wsd(int n)
{
   int iprms[1],i,ia;

   if (NPROC>1)
   {
      iprms[0]=n;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error(iprms[0]!=n,1,"reserve_wsd [wspace.c]",
            "Parameter n is not global");
   }

   if (n>0)
   {
      error((nwsdt+n)>nsdt,1,"reserve_wsd [wspace.c]",
            "Requested too many fields (tot=%d,use=%d,req=%d)",nsdt,nwsdt,n);

      ia=nwsdt;
      nwsd[iwsd]=n;
      nwsdt+=n;
      iwsd+=1;

      for (i=ia;i<(ia+n);i++)
         wsd[i]=wsd0[i];

      return wsd+ia;
   }
   else
      return NULL;
}


int release_wsd(void)
{
   int n,i;

   if (nwsdt==0)
      return 0;
   else
   {
      iwsd-=1;
      n=nwsd[iwsd];
      nwsdt-=n;
      nwsd[iwsd]=0;

      for (i=nwsdt;i<(nwsdt+n);i++)
         wsd[i]=NULL;

      return n;
   }
}


int wsd_size(void)
{
   return nwsdt;
}


static void set_nvec(void)
{
   int *bs;
   dfl_parms_t dfl;

   dfl=dfl_parms();

   error_root(dfl.Ns==0,1,"set_nvec [wspace.c]",
              "Deflation subspace parameters are not set");

   bs=dfl.bs;
   nvec=VOLUME+FACE0*bs[0]+FACE1*bs[1]+FACE2*bs[2]+FACE3*bs[3];
   nvec/=(bs[0]*bs[1]*bs[2]*bs[3]);
   nvec*=dfl.Ns;
}


void alloc_wv(int n)
{
   int i;
   complex *v,*vm;

   if (n==nvt)
      return;

   error_root(nwvt!=0,1,"alloc_wv [wspace.c]","Fields are in use");

   if (nvt>0)
   {
      free(nwv);
      afree(wv0[0]);
      free(wv0);
      nwv=NULL;
      wv0=NULL;
      wv=NULL;
   }

   nvt=n;
   iwv=0;
   nwvt=0;

   if (nvt>0)
   {
      set_nvec();
      nwv=malloc(nvt*sizeof(*nwv));
      wv0=malloc(2*nvt*sizeof(*wv0));
      wv=wv0+nvt;

      error((nwv==NULL)||(wv0==NULL),1,"alloc_wv [wspace.c]",
            "Unable to allocate index arrays");

      wv0[0]=amalloc(nvt*nvec*sizeof(**wv0),ALIGN);

      error(wv0[0]==NULL,1,"alloc_wv [wspace.c]",
            "Unable to allocate workspace");

      for (i=0;i<nvt;i++)
      {
         if (i>0)
            wv0[i]=wv0[i-1]+nvec;

         nwv[i]=0;
         wv[i]=NULL;
      }

      v=wv0[0];
      vm=v+nvt*nvec;

      for (;v<vm;v++)
         (*v)=v0;
   }
}


complex **reserve_wv(int n)
{
   int iprms[1],i,ia;

   if (NPROC>1)
   {
      iprms[0]=n;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error(iprms[0]!=n,1,"reserve_wv [wspace.c]",
            "Parameter n is not global");
   }

   if (n>0)
   {
      error((nwvt+n)>nvt,1,"reserve_wv [wspace.c]",
            "Requested too many fields (tot=%d,use=%d,req=%d)",nvt,nwvt,n);

      ia=nwvt;
      nwv[iwv]=n;
      nwvt+=n;
      iwv+=1;

      for (i=ia;i<(ia+n);i++)
         wv[i]=wv0[i];

      return wv+ia;
   }
   else
      return NULL;
}


int release_wv(void)
{
   int n,i;

   if (nwvt==0)
      return 0;
   else
   {
      iwv-=1;
      n=nwv[iwv];
      nwvt-=n;
      nwv[iwv]=0;

      for (i=nwvt;i<(nwvt+n);i++)
         wv[i]=NULL;

      return n;
   }
}


int wv_size(void)
{
   return nwvt;
}


void alloc_wvd(int n)
{
   int i;
   complex_dble *vd,*vm;

   if (n==nvdt)
      return;

   error_root(nwvdt!=0,1,"alloc_wvd [wspace.c]","Fields are in use");

   if (nvdt>0)
   {
      free(nwvd);
      afree(wvd0[0]);
      free(wvd0);
      nwvd=NULL;
      wvd0=NULL;
      wvd=NULL;
   }

   nvdt=n;
   iwvd=0;
   nwvdt=0;

   if (nvdt>0)
   {
      set_nvec();
      nwvd=malloc(nvdt*sizeof(*nwvd));
      wvd0=malloc(2*nvdt*sizeof(*wvd0));
      wvd=wvd0+nvdt;

      error((nwvd==NULL)||(wvd0==NULL),1,"alloc_wvd [wspace.c]",
            "Unable to allocate index arrays");

      wvd0[0]=amalloc(nvdt*nvec*sizeof(**wvd0),ALIGN);

      error(wvd0[0]==NULL,1,"alloc_wvd [wspace.c]",
            "Unable to allocate workspace");

      for (i=0;i<nvdt;i++)
      {
         if (i>0)
            wvd0[i]=wvd0[i-1]+nvec;

         nwvd[i]=0;
         wvd[i]=NULL;
      }

      vd=wvd0[0];
      vm=vd+nvdt*nvec;

      for (;vd<vm;vd++)
         (*vd)=vd0;
   }
}


complex_dble **reserve_wvd(int n)
{
   int iprms[1],i,ia;

   if (NPROC>1)
   {
      iprms[0]=n;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error(iprms[0]!=n,1,"reserve_wvd [wspace.c]",
            "Parameter n is not global");
   }

   if (n>0)
   {
      error((nwvdt+n)>nvdt,1,"reserve_wvd [wspace.c]",
            "Requested too many fields (tot=%d,use=%d,req=%d)",nvdt,nwvdt,n);

      ia=nwvdt;
      nwvd[iwvd]=n;
      nwvdt+=n;
      iwvd+=1;

      for (i=ia;i<(ia+n);i++)
         wvd[i]=wvd0[i];

      return wvd+ia;
   }
   else
      return NULL;
}


int release_wvd(void)
{
   int n,i;

   if (nwvdt==0)
      return 0;
   else
   {
      iwvd-=1;
      n=nwvd[iwvd];
      nwvdt-=n;
      nwvd[iwvd]=0;

      for (i=nwvdt;i<(nwvdt+n);i++)
         wvd[i]=NULL;

      return n;
   }
}


int wvd_size(void)
{
   return nwvdt;
}
