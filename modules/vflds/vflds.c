
/*******************************************************************************
*
* File vflds.c
*
* Copyright (C) 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation and initialization of the global vector fields related to
* the deflation subspace.
*
* The externally accessible functions are
*
*   complex **vflds(void)
*     Returns the base address of the global single-precision vector fields
*     (see the notes). The fields are allocated and initialized to zero if
*     they are not already allocated.
*
*   complex_dble **vdflds(void)
*     Returns the base address of the global double-precision vector fields
*     (see the notes). The fields are allocated and initialized to zero if
*     they are not already allocated.
*
* Notes:
*
* The vector fields made available through the programs in this module
* are arrays complex numbers. Eventually they contain the global modes
* used to deflate the little Dirac operator that represents the action
* of the Wilson-Dirac operator in the deflation subspace.
*
* Each vector field has Ns*nb elements, where Ns is the number of local
* deflation modes and nb the number blocks in the DFL_BLOCKS grid. The
* elements of the fields are interpreted as the components of a spinor
* field along the deflation subspace.
*
* The programs vflds() and vdflds() return arrays of 2*Ns and Ns global
* vector fields, respectively. They may involve global communications and
* must therefore be called simultaneously on all processes.
*
*******************************************************************************/

#define VFLDS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "vflds.h"
#include "global.h"

static int Ns,nv=0;
static complex **vs=NULL,**v;
static complex_dble **vds=NULL,**vd;


static void vfld_size(void)
{
   int *bs;
   dfl_parms_t dfl;

   error_root(sizeof(complex)!=(2*sizeof(float)),1,
              "vfld_size [vflds.c]",
              "The complex structures are not properly packed");
   error_root(sizeof(complex_dble)!=(2*sizeof(double)),1,
              "vfld_size [vflds.c]",
              "The complex_dble structures are not properly packed");

   dfl=dfl_parms();
   bs=dfl.bs;
   Ns=dfl.Ns;

   error_root(dfl.Ns==0,1,"vfld_size [vflds.c]",
         "The deflation subspace parameters are not set");

   nv=VOLUME/(bs[0]*bs[1]*bs[2]*bs[3]);
   nv*=Ns;
}


static void alloc_vflds(void)
{
   int n;
   complex *w;

   if (nv==0)
      vfld_size();

   vs=malloc(4*Ns*sizeof(*vs));
   w=amalloc(2*Ns*nv*sizeof(*w),ALIGN);

   error((vs==NULL)||(w==NULL),1,"alloc_vflds [vflds.c]",
         "Unable to allocate vector fields");

   set_v2zero(2*Ns*nv,w);
   v=vs+2*Ns;

   for (n=0;n<(2*Ns);n++)
   {
      v[n]=NULL;
      vs[n]=w;
      w+=nv;
   }
}


static void alloc_vdflds(void)
{
   int n;
   complex_dble *wd;

   if (nv==0)
      vfld_size();

   vds=malloc(2*Ns*sizeof(*vds));
   wd=amalloc(Ns*nv*sizeof(*wd),ALIGN);

   error((vds==NULL)||(wd==NULL),1,"alloc_vdflds [vflds.c]",
         "Unable to allocate vector fields");

   set_vd2zero(Ns*nv,wd);
   vd=vds+Ns;

   for (n=0;n<Ns;n++)
   {
      vd[n]=NULL;
      vds[n]=wd;
      wd+=nv;
   }
}


complex **vflds(void)
{
   int n;

   if (vs==NULL)
      alloc_vflds();

   for (n=0;n<(2*Ns);n++)
      v[n]=vs[n];

   return v;
}


complex_dble **vdflds(void)
{
   int n;

   if (vds==NULL)
      alloc_vdflds();

   for (n=0;n<Ns;n++)
      vd[n]=vds[n];

   return vd;
}
