
/*******************************************************************************
*
* File uidx.c
*
* Copyright (C) 2010, 2011, 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Labeling of the link variables on the faces of the local lattice.
*
* The externally accessible functions are
*
*   uidx_t *uidx(void)
*     Returns an array idx[4] of uidx_t structures containing the offsets
*     of the link variables at the faces of the local lattice.
*
*   void plaq_uidx(int n,int ix,int *ip)
*     Calculates the offsets ip[4] of the links in the (mu,nu)-plaquette at
*     the point on the local lattice with label ix. The indices (mu,nu) are
*     determined by the parameter n=0,..,5.
*
* Notes:
*
* The layout of the double-precision gauge field array and contents of the
* index structures returned by uidx() are described in the file README.uidx
* in this directory. The index arrays calculated by uidx() are determined
* by the local geometry of the lattice and are therefore independent of the
* boundary conditions.
*
* There are six planes
*
*  (mu,nu)={(0,1),(0,2),(0,3),(2,3),(3,1),(1,2)}
*
* labeled by an integer n running from 0 to 5 and the links in the
* (mu,nu)-plaquette at the point x are ordered such that
*
*   ip[0] -> U(x,mu)
*   ip[1] -> U(x+mu,nu)
*   ip[2] -> U(x,nu)
*   ip[3] -> U(x+nu,mu)
*
* In the program plaq_uidx() it is taken for granted that 0<=ix<VOLUME.
*
* If SF or open-SF boundary conditions are chosen, the offsets ip[4]
* returned by plaq_uidx() at global time NPROC0*L0-1 take into account the
* fact that the boundary values of the gauge field are stored at the end
* of the field array. On all MPI processes, and for all boundary conditions,
* the correct field variables are thus found at the calculated offsets.
*
*******************************************************************************/

#define UIDX_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "global.h"

#define N0 (NPROC0*L0)

static const int plns[6][2]={{0,1},{0,2},{0,3},{2,3},{3,1},{1,2}};
static int type,nfc[4],ofs[4],snu[4],init=0;
static uidx_t idx[4];


static void alloc_idx(void)
{
   int mu,nu0,nuk;
   int *iu0,*iuk;

   error(iup[0][0]==0,1,"alloc_idx [uidx.c]",
         "Geometry arrays are not set");

   type=bc_type();
   nfc[0]=FACE0/2;
   nfc[1]=FACE1/2;
   nfc[2]=FACE2/2;
   nfc[3]=FACE3/2;

   ofs[0]=VOLUME+(FACE0/2);
   ofs[1]=ofs[0]+(FACE0/2)+(FACE1/2);
   ofs[2]=ofs[1]+(FACE1/2)+(FACE2/2);
   ofs[3]=ofs[2]+(FACE2/2)+(FACE3/2);

   snu[0]=0;
   snu[1]=snu[0]+(FACE0/2);
   snu[2]=snu[1]+(FACE1/2);
   snu[3]=snu[2]+(FACE2/2);

   if (BNDRY>0)
   {
      iu0=malloc(7*(BNDRY/4)*sizeof(*iu0));
      error(iu0==NULL,1,"alloc_idx [uidx.c]",
            "Unable to allocate index array");
      iuk=iu0+(BNDRY/4);
   }
   else
   {
      iu0=NULL;
      iuk=NULL;
   }

   for (mu=0;mu<4;mu++)
   {
      nu0=nfc[mu];
      nuk=6*nfc[mu];

      idx[mu].nu0=nu0;
      idx[mu].nuk=nuk;

      if (nu0>0)
      {
         idx[mu].iu0=iu0;
         idx[mu].iuk=iuk;
         iu0+=nu0;
         iuk+=nuk;
      }
      else
      {
         idx[mu].iu0=NULL;
         idx[mu].iuk=NULL;
      }
   }
}


static int offset(int ix,int mu)
{
   int iy,ib;

   if (ix<(VOLUME/2))
   {
      iy=iup[ix][mu];

      if (iy<VOLUME)
         return 8*(iy-(VOLUME/2))+2*mu+1;
      else
      {
         ib=iy-ofs[mu]-(BNDRY/2);

         return 4*VOLUME+snu[mu]+ib;
      }
   }
   else
      return 8*(ix-(VOLUME/2))+2*mu;
}


static void set_idx(void)
{
   int mu,nu,k;
   int ib,iy,iz;
   int nu0,*iu0,*iuk;

   alloc_idx();

   for (mu=0;mu<4;mu++)
   {
      nu0=idx[mu].nu0;
      iu0=idx[mu].iu0;
      iuk=idx[mu].iuk;

      for (ib=0;ib<nu0;ib++)
      {
         iy=ib+ofs[mu]+(BNDRY/2);
         iz=map[iy-VOLUME];
         iu0[ib]=8*(iz-(VOLUME/2))+2*mu+1;
      }

      for (ib=0;ib<nu0;ib++)
      {
         iy=ib+ofs[mu];
         iz=map[iy-VOLUME];

         for (k=0;k<3;k++)
         {
            nu=k+(k>=mu);
            iuk[3*ib+k]=offset(iz,nu);
         }
      }

      for (ib=0;ib<nu0;ib++)
      {
         iy=ib+ofs[mu]+(BNDRY/2);
         iz=map[iy-VOLUME];

         for (k=0;k<3;k++)
         {
            nu=k+(k>=mu);
            iuk[3*(ib+nu0)+k]=offset(iz,nu);
         }
      }
   }

   init=1;
}


uidx_t *uidx(void)
{
   if (init==0)
      set_idx();

   return idx;
}


void plaq_uidx(int n,int ix,int *ip)
{
   int mu,nu;
   int iy,ic;

   if (init==0)
      set_idx();

   mu=plns[n][0];
   nu=plns[n][1];

   ip[0]=offset(ix,mu);

   if ((mu==0)&&(global_time(ix)==(N0-1))&&((type==1)||(type==2)))
   {
      ip[1]=4*VOLUME+7*(BNDRY/4)+nu-1;
   }
   else
   {
      iy=iup[ix][mu];

      if (iy<VOLUME)
         ip[1]=offset(iy,nu);
      else
      {
         if (iy<(VOLUME+(BNDRY/2)))
            ic=iy-VOLUME-nfc[mu];
         else
            ic=iy-VOLUME-(BNDRY/2);

         ip[1]=4*VOLUME+(BNDRY/4)+3*ic+nu-(nu>mu);
      }
   }

   ip[2]=offset(ix,nu);
   iy=iup[ix][nu];

   if (iy<VOLUME)
      ip[3]=offset(iy,mu);
   else
   {
      if (iy<(VOLUME+(BNDRY/2)))
         ic=iy-VOLUME-nfc[nu];
      else
         ic=iy-VOLUME-(BNDRY/2);

      ip[3]=4*VOLUME+(BNDRY/4)+3*ic+mu-(mu>nu);
   }
}
