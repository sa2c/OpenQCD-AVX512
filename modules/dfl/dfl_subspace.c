
/*******************************************************************************
*
* File dfl_subspace.c
*
* Copyright (C) 2007, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Basic utility programs related to the deflation subspace.
*
* The externally accessible functions are
*
*   void dfl_sd2vd(spinor_dble *sd,complex_dble *vd)
*     Assigns the components of the global double-precision spinor field
*     sd along the deflation subspace to the double-precision vector
*     field vd.
*
*   void dfl_vd2sd(complex_dble *vd,spinor_dble *sd)
*     Assigns the element of the deflation subspace corresponding to the
*     double-precision vector field vd to the global double-precision spinor
*     field sd.
*
*   void dfl_sub_vd2sd(complex_dble *vd,spinor_dble *sd)
*     Subtracts the element of the deflation subspace corresponding to the
*     double-precision vector field vd from the global double-precision
*     spinor field sd.
*
*   void dfl_s2v(spinor *s,complex *v)
*     Assigns the components of the global single-precision spinor field
*     s along the deflation subspace to the single-precision vector
*     field v.
*
*   void dfl_v2s(complex *v,spinor *s)
*     Assigns the element of the deflation subspace corresponding to the
*     single-precision vector field v to the global single-precision spinor
*     field s.
*
*   void dfl_sub_v2s(complex *v,spinor *s)
*     Subtracts the element of the deflation subspace corresponding to the
*     double-precision vector field v from the global single-precision spinor
*     field s.
*
*   void dfl_subspace(spinor **mds)
*     Copies the global single-precision spinor fields mds[0],..,mds[Ns-1]
*     to the fields b.sd[1],..,b.sd[Ns] on the blocks b of the DFL_BLOCKS
*     grid. The block fields are then orthonormalized and are assigned to
*     the single-precision block fields b.s[1],..,b.s[Ns].
*      In this basis of fields, the modes mds[0],..,mds[Ns-1] are given by
*     fields vmds[0],..,vmds[Ns-1] of Ns*nb complex numbers, where nb is
*     the number of blocks in the block grid. These fields are assigned to
*     the last Ns single-precision vector fields of the array returned by
*     vflds() [vflds/vflds.c].
*
* Notes:
*
* The deflation subspace is spanned by the fields (*b).sd[1],..,(*b).sd[Ns]
* on the blocks b of the DFL_BLOCKS grid. The number Ns of fields is set by
* the program dfl_set_parms() [flags/dfl_parms.c].
*
* Any spinor field in the deflation subspace is a linear combination of the
* basis elements on the blocks. The associated complex coefficients form a
* vector field of the type described in vflds/vflds.c. Such fields are thus
* in one-to-one correspondence with the deflation modes. In particular, the
* deflation subspace contains the global spinor fields from which it was
* created by the program dfl_subspace().
*
* The program dfl_subspace() allocates the DFL_BLOCKS block grid if it is
* not already allocated. This program involves global operations and must be
* called simultaneously on all processes.
*
*******************************************************************************/

#define DFL_SUBSPACE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "linalg.h"
#include "sflds.h"
#include "block.h"
#include "vflds.h"
#include "dfl.h"
#include "global.h"


void dfl_sd2vd(spinor_dble *sd,complex_dble *vd)
{
   int Ns,nb,nbh,isw;
   int n,m,i,vol;
   block_t *b;
   spinor_dble **sdb;
   dfl_parms_t dfl;

   dfl=dfl_parms();
   Ns=dfl.Ns;
   b=blk_list(DFL_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;

   for (n=0;n<nb;n++)
   {
      if (n<nbh)
         m=n+isw*nbh;
      else
         m=n-isw*nbh;

      assign_sd2sdblk(DFL_BLOCKS,m,ALL_PTS,sd,0);
      sdb=b[m].sd;

      for (i=1;i<=Ns;i++)
      {
         (*vd)=spinor_prod_dble(vol,0,sdb[i],sdb[0]);
         vd+=1;
      }
   }
}


void dfl_vd2sd(complex_dble *vd,spinor_dble *sd)
{
   int Ns,nb,nbh,isw;
   int n,m,i,vol;
   block_t *b;
   spinor_dble **sdb;
   dfl_parms_t dfl;

   dfl=dfl_parms();
   Ns=dfl.Ns;
   b=blk_list(DFL_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;

   for (n=0;n<nb;n++)
   {
      if (n<nbh)
         m=n+isw*nbh;
      else
         m=n-isw*nbh;

      sdb=b[m].sd;
      set_sd2zero(vol,sdb[0]);

      for (i=1;i<=Ns;i++)
      {
         mulc_spinor_add_dble(vol,sdb[0],sdb[i],*vd);
         vd+=1;
      }

      assign_sdblk2sd(DFL_BLOCKS,m,ALL_PTS,0,sd);
   }
}


void dfl_sub_vd2sd(complex_dble *vd,spinor_dble *sd)
{
   int Ns,nb,nbh,isw;
   int n,m,i,vol;
   complex_dble z;
   block_t *b;
   spinor_dble **sdb;
   dfl_parms_t dfl;

   dfl=dfl_parms();
   Ns=dfl.Ns;
   b=blk_list(DFL_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;

   for (n=0;n<nb;n++)
   {
      if (n<nbh)
         m=n+isw*nbh;
      else
         m=n-isw*nbh;

      assign_sd2sdblk(DFL_BLOCKS,m,ALL_PTS,sd,0);
      sdb=b[m].sd;

      for (i=1;i<=Ns;i++)
      {
         z.re=-(*vd).re;
         z.im=-(*vd).im;
         mulc_spinor_add_dble(vol,sdb[0],sdb[i],z);
         vd+=1;
      }

      assign_sdblk2sd(DFL_BLOCKS,m,ALL_PTS,0,sd);
   }
}


void dfl_s2v(spinor *s,complex *v)
{
   int Ns,nb,nbh,isw;
   int n,m,i,vol;
   block_t *b;
   spinor **sb;
   dfl_parms_t dfl;

   dfl=dfl_parms();
   Ns=dfl.Ns;
   b=blk_list(DFL_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;

   for (n=0;n<nb;n++)
   {
      if (n<nbh)
         m=n+isw*nbh;
      else
         m=n-isw*nbh;

      assign_s2sblk(DFL_BLOCKS,m,ALL_PTS,s,0);
      sb=b[m].s;

      for (i=1;i<=Ns;i++)
      {
         (*v)=spinor_prod(vol,0,sb[i],sb[0]);
         v+=1;
      }
   }
}


void dfl_v2s(complex *v,spinor *s)
{
   int Ns,nb,nbh,isw;
   int n,m,i,vol;
   block_t *b;
   spinor **sb;
   dfl_parms_t dfl;

   dfl=dfl_parms();
   Ns=dfl.Ns;
   b=blk_list(DFL_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;

   for (n=0;n<nb;n++)
   {
      if (n<nbh)
         m=n+isw*nbh;
      else
         m=n-isw*nbh;

      sb=b[m].s;
      set_s2zero(vol,sb[0]);

      for (i=1;i<=Ns;i++)
      {
         mulc_spinor_add(vol,sb[0],sb[i],*v);
         v+=1;
      }

      assign_sblk2s(DFL_BLOCKS,m,ALL_PTS,0,s);
   }
}


void dfl_sub_v2s(complex *v,spinor *s)
{
   int Ns,nb,nbh,isw;
   int n,m,i,vol;
   complex z;
   block_t *b;
   spinor **sb;
   dfl_parms_t dfl;

   dfl=dfl_parms();
   Ns=dfl.Ns;
   b=blk_list(DFL_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;

   for (n=0;n<nb;n++)
   {
      if (n<nbh)
         m=n+isw*nbh;
      else
         m=n-isw*nbh;

      assign_s2sblk(DFL_BLOCKS,m,ALL_PTS,s,0);
      sb=b[m].s;

      for (i=1;i<=Ns;i++)
      {
         z.re=-(*v).re;
         z.im=-(*v).im;
         mulc_spinor_add(vol,sb[0],sb[i],z);
         v+=1;
      }

      assign_sblk2s(DFL_BLOCKS,m,ALL_PTS,0,s);
   }
}


void dfl_subspace(spinor **mds)
{
   int Ns,nb,nbh,isw;
   int n,m,i,j,vol;
   complex **vs,*v;
   complex_dble z;
   block_t *b;
   spinor **sb;
   spinor_dble **sdb;
   dfl_parms_t dfl;

   dfl=dfl_parms();
   Ns=dfl.Ns;

   error_root(Ns==0,1,"dfl_subspace [dfl_subspace.c]",
              "Deflation subspace parameters are not set");

   b=blk_list(DFL_BLOCKS,&nb,&isw);

   if (nb==0)
   {
      alloc_bgr(DFL_BLOCKS);
      b=blk_list(DFL_BLOCKS,&nb,&isw);
   }

   nbh=nb/2;
   vol=(*b).vol;
   vs=vflds()+Ns;

   for (n=0;n<nb;n++)
   {
      if (n<nbh)
         m=n+isw*nbh;
      else
         m=n-isw*nbh;

      sb=b[m].s;
      sdb=b[m].sd;

      for (i=1;i<=Ns;i++)
      {
         assign_s2sdblk(DFL_BLOCKS,m,ALL_PTS,mds[i-1],i);
         v=vs[i-1]+Ns*n;

         for (j=1;j<i;j++)
         {
            z=spinor_prod_dble(vol,0,sdb[j],sdb[i]);

            (*v).re=(float)(z.re);
            (*v).im=(float)(z.im);
            v+=1;

            z.re=-(z).re;
            z.im=-(z).im;
            mulc_spinor_add_dble(vol,sdb[i],sdb[j],z);
         }

         (*v).re=(float)(normalize_dble(vol,0,sdb[i]));
         (*v).im=0.0f;
         v+=1;

         for (j=(i+1);j<=Ns;j++)
         {
            (*v).re=0.0f;
            (*v).im=0.0f;
            v+=1;
         }

         assign_sd2s(vol,sdb[i],sb[i]);
      }
   }

   set_flags(ERASED_AW);
   set_flags(ERASED_AWHAT);
}
