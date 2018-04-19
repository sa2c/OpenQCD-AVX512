
/*******************************************************************************
*
* File force0.c
*
* Copyright (C) 2005, 2009-2014, 2016 Martin Luescher, John Bulava
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Action of the double-precision gauge field and associated force.
*
* The externally accessible functions are
*
*   void plaq_frc(void)
*     Computes the force deriving from the Wilson plaquette action,
*     omitting the prefactor 1/g0^2, and assigns the result to the MD
*     force field. In the case of open, SF or open-SF boundary conditions,
*     the boundary improvement coefficients are set to their tree-level
*     value independently of the values stored in the parameter data base.
*
*   void force0(double c)
*     Computes the force deriving from the gauge action, including the
*     prefactor 1/g0^2, multiplies the calculated force by c and assigns
*     the result to the MD force field. The coupling g0 and the other
*     parameters of the gauge action are retrieved from the parameter
*     data base.
*
*   double action0(int icom)
*     Computes the local part of the gauge action including the prefactor
*     1/g0^2. The coupling g0 and the other parameters of the action are
*     retrieved from the parameter data base. The program returns the sum
*     of the local parts of the action over all MPI processes if icom=1
*     and otherwise just the local part.
*
* Notes:
*
* See the notes doc/gauge_action.pdf for the definition of the gauge action
* and a description of the computation of the force deriving from it. The
* molecular-dynamics (MD) force field is the one returned by the program
* mdflds() (see mdflds/mdflds.c).
*
* On the links in the local lattice where the static link variables reside,
* the programs plaq_frc() and force0() set the force field to zero.
*
* The programs in this module perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define FORCE0_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "flags.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "mdflds.h"
#include "forces.h"
#include "global.h"

#define N0 (NPROC0*L0)

static const int plns[6][2]={{0,1},{0,2},{0,3},{2,3},{3,1},{1,2}};
static int nfc[8],ofs[8],hofs[8],ism,init=0;
static su3_dble *udb,*hdb;
static su3_dble wd[3] ALIGNED16;
static su3_dble vd[4] ALIGNED16;
static su3_alg_dble X ALIGNED16;


static void set_ofs(void)
{
   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   ofs[0]=VOLUME;
   ofs[1]=ofs[0]+(FACE0/2);
   ofs[2]=ofs[1]+(FACE0/2);
   ofs[3]=ofs[2]+(FACE1/2);
   ofs[4]=ofs[3]+(FACE1/2);
   ofs[5]=ofs[4]+(FACE2/2);
   ofs[6]=ofs[5]+(FACE2/2);
   ofs[7]=ofs[6]+(FACE3/2);

   hofs[0]=0;
   hofs[1]=hofs[0]+3*FACE0;
   hofs[2]=hofs[1]+3*FACE0;
   hofs[3]=hofs[2]+3*FACE1;
   hofs[4]=hofs[3]+3*FACE1;
   hofs[5]=hofs[4]+3*FACE2;
   hofs[6]=hofs[5]+3*FACE2;
   hofs[7]=hofs[6]+3*FACE3;

   init+=2;
}


static void set_staples(int n,int ix,int ia)
{
   int mu,nu,ifc;
   int iy,ib,ip[4];

   mu=plns[n][0];
   nu=plns[n][1];

   if (!ia)
   {
      iy=idn[ix][nu];

      if (iy<VOLUME)
      {
         plaq_uidx(n,iy,ip);

         su3xsu3(udb+ip[0],udb+ip[1],wd+2);
         su3dagxsu3(udb+ip[2],wd+2,vd);
      }
      else
      {
         ifc=2*nu;

         if (iy<(VOLUME+(BNDRY/2)))
            ib=iy-ofs[ifc];
         else
            ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

         vd[0]=hdb[hofs[ifc]+3*ib+mu-(mu>nu)];
      }
   }

   iy=iup[ix][mu];

   if (iy<VOLUME)
   {
      plaq_uidx(n,iy,ip);

      su3xsu3dag(udb+ip[1],udb+ip[3],wd+2);
      su3xsu3(udb+ip[0],wd+2,vd+1);
   }
   else
   {
      ifc=2*mu+1;

      if (iy<(VOLUME+(BNDRY/2)))
         ib=iy-ofs[ifc];
      else
         ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

      vd[1]=hdb[hofs[ifc]+3*ib+nu-(nu>mu)];
   }

   if (!ia)
   {
      iy=idn[ix][mu];

      if (iy<VOLUME)
      {
         plaq_uidx(n,iy,ip);

         su3xsu3(udb+ip[2],udb+ip[3],wd+2);
         su3dagxsu3(udb+ip[0],wd+2,vd+2);
      }
      else
      {
         ifc=2*mu;

         if (iy<(VOLUME+(BNDRY/2)))
            ib=iy-ofs[ifc];
         else
            ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

         vd[2]=hdb[hofs[ifc]+3*ib+nu-(nu>mu)];
      }
   }

   iy=iup[ix][nu];

   if (iy<VOLUME)
   {
      plaq_uidx(n,iy,ip);

      su3xsu3dag(udb+ip[3],udb+ip[1],wd+2);
      su3xsu3(udb+ip[2],wd+2,vd+3);
   }
   else
   {
      ifc=2*nu+1;

      if (iy<(VOLUME+(BNDRY/2)))
         ib=iy-ofs[ifc];
      else
         ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

      vd[3]=hdb[hofs[ifc]+3*ib+mu-(mu>nu)];
   }
}


void plaq_frc(void)
{
   int bc,n,ix,t,ip[4];
   double r;
   su3_alg_dble *fdb;
   mdflds_t *mdfs;

   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();

   bc=bc_type();
   udb=udfld();
   mdfs=mdflds();
   fdb=(*mdfs).frc;
   set_frc2zero();

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t<(N0-1))||(bc!=0))
      {
         for (n=0;n<3;n++)
         {
            plaq_uidx(n,ix,ip);

            su3xsu3dag(udb+ip[1],udb+ip[3],wd);
            su3dagxsu3(udb+ip[2],udb+ip[0],wd+1);

            if ((t<(N0-1))||(bc==3))
            {
               prod2su3alg(wd,wd+1,&X);
               _su3_alg_add_assign(*(fdb+ip[1]),X);
            }

            prod2su3alg(wd+1,wd,&X);
            _su3_alg_sub_assign(*(fdb+ip[3]),X);

            su3xsu3dag(wd,udb+ip[2],wd+1);
            prod2su3alg(udb+ip[0],wd+1,&X);
            _su3_alg_add_assign(*(fdb+ip[0]),X);

            if ((t>0)||(bc!=1))
            {
               _su3_alg_sub_assign(*(fdb+ip[2]),X);
            }
         }
      }

      if ((t>0)||(bc!=1))
      {
         r=1.0;

         if (((t==0)&&(bc!=3))||((t==(N0-1))&&(bc==0)))
            r=0.5;

         for (n=3;n<6;n++)
         {
            plaq_uidx(n,ix,ip);

            su3xsu3dag(udb+ip[1],udb+ip[3],wd);
            su3dagxsu3(udb+ip[2],udb+ip[0],wd+1);
            prod2su3alg(wd,wd+1,&X);
            _su3_alg_mul_add_assign(*(fdb+ip[1]),r,X);
            prod2su3alg(wd+1,wd,&X);
            _su3_alg_mul_sub_assign(*(fdb+ip[3]),r,X);

            su3xsu3dag(wd,udb+ip[2],wd+1);
            prod2su3alg(udb+ip[0],wd+1,&X);
            _su3_alg_mul_add_assign(*(fdb+ip[0]),r,X);
            _su3_alg_mul_sub_assign(*(fdb+ip[2]),r,X);
         }
      }
   }

   add_bnd_frc();
}


void force0(double c)
{
   int bc,n,ix,t,ip[4];
   double c0,c1,*cG;
   double r0,r1;
   su3_alg_dble *fdb;
   mdflds_t *mdfs;
   lat_parms_t lat;
   bc_parms_t bcp;

   lat=lat_parms();
   c*=(lat.beta/6.0);
   c0=lat.c0;
   c1=lat.c1;

   bcp=bc_parms();
   bc=bcp.type;
   cG=bcp.cG;

   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();

   udb=udfld();
   mdfs=mdflds();
   fdb=(*mdfs).frc;
   set_frc2zero();

   if (c0==1.0)
      hdb=NULL;
   else
   {
      if ((init&0x2)==0)
         set_ofs();

      if (query_flags(BSTAP_UP2DATE)!=1)
         set_bstap();
      hdb=bstap();
   }

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t<(N0-1))||(bc!=0))
      {
         r0=c*c0;
         r1=c*c1;

         if ((t==0)&&(bc==1))
            r0*=cG[0];
         else if ((t==(N0-1))&&(bc!=3))
            r0*=cG[1];

         for (n=0;n<3;n++)
         {
            plaq_uidx(n,ix,ip);

            su3xsu3dag(udb+ip[1],udb+ip[3],wd);
            su3dagxsu3(udb+ip[2],udb+ip[0],wd+1);

            if ((t<(N0-1))||(bc==3))
            {
               prod2su3alg(wd,wd+1,&X);
               _su3_alg_mul_add_assign(*(fdb+ip[1]),r0,X);
	    }

            prod2su3alg(wd+1,wd,&X);
            _su3_alg_mul_sub_assign(*(fdb+ip[3]),r0,X);

            su3xsu3dag(wd,udb+ip[2],wd+1);
            prod2su3alg(udb+ip[0],wd+1,&X);
            _su3_alg_mul_add_assign(*(fdb+ip[0]),r0,X);

            if ((t>0)||(bc!=1))
            {
               _su3_alg_mul_sub_assign(*(fdb+ip[2]),r0,X);
            }

            if (c0!=1.0)
            {
               set_staples(n,ix,0);

               if ((t==0)&&(bc==1))
               {
                  su3xsu3(wd+1,udb+ip[0],wd+2);
                  su3xsu3(udb+ip[0],wd+2,wd+2);

                  prod2su3alg(wd+1,wd+2,&X);
                  _su3_alg_mul_add_assign(*(fdb+ip[1]),r1,X);

                  prod2su3alg(wd+2,wd+1,&X);
                  _su3_alg_mul_add_assign(*(fdb+ip[0]),r1,X);

                  su3dagxsu3(udb+ip[2],wd+2,wd+2);

                  prod2su3alg(wd+2,wd,&X);
                  _su3_alg_mul_sub_assign(*(fdb+ip[3]),r1,X);
               }

               if ((t==(N0-1))&&(bc!=3))
               {
                  su3xsu3(wd+1,udb+ip[0],wd+2);
                  su3xsu3(udb+ip[0],wd+2,wd+2);

                  prod2su3alg(wd+2,wd+1,&X);
                  _su3_alg_mul_add_assign(*(fdb+ip[0]),r1,X);
                  _su3_alg_mul_sub_assign(*(fdb+ip[2]),r1,X);

                  su3dagxsu3(udb+ip[2],wd+2,wd+2);

                  prod2su3alg(wd+2,wd,&X);
                  _su3_alg_mul_sub_assign(*(fdb+ip[3]),r1,X);
               }

               if ((t<(N0-1))||(bc==3))
               {
                  prod2su3alg(wd+1,vd,&X);
                  _su3_alg_mul_add_assign(*(fdb+ip[1]),r1,X);
               }

               if ((t>0)||(bc!=1))
               {
                  prod2su3alg(vd,wd+1,&X);
                  _su3_alg_mul_sub_assign(*(fdb+ip[2]),r1,X);
               }

               su3dagxsu3(udb+ip[2],vd,wd+1);
               prod2su3alg(wd+1,wd,&X);
               _su3_alg_mul_sub_assign(*(fdb+ip[3]),r1,X);

               if ((t<(N0-2))||((t==(N0-2))&&(bc!=0))||(bc==3))
               {
                  su3xsu3dag(udb+ip[3],vd+1,wd+1);
                  su3xsu3dag(wd+1,udb+ip[0],wd+2);
                  prod2su3alg(udb+ip[2],wd+2,&X);
                  _su3_alg_mul_sub_assign(*(fdb+ip[0]),r1,X);

                  if ((t>0)||(bc!=1))
                  {
                     _su3_alg_mul_add_assign(*(fdb+ip[2]),r1,X);
                  }

                  prod2su3alg(wd+2,udb+ip[2],&X);
                  _su3_alg_mul_add_assign(*(fdb+ip[3]),r1,X);
               }

               if ((t>0)||(bc==3))
               {
                  su3xsu3dag(wd,vd+2,wd+1);
                  prod2su3alg(udb+ip[0],wd+1,&X);
                  _su3_alg_mul_add_assign(*(fdb+ip[0]),r1,X);

                  if ((t<(N0-1))||(bc==3))
                  {
                     prod2su3alg(wd+1,udb+ip[0],&X);
                     _su3_alg_mul_add_assign(*(fdb+ip[1]),r1,X);
                  }

                  su3dagxsu3(vd+2,udb+ip[0],wd+1);
                  prod2su3alg(wd+1,wd,&X);
                  _su3_alg_mul_sub_assign(*(fdb+ip[3]),r1,X);
               }

               su3xsu3dag(udb+ip[1],vd+3,wd);
               su3xsu3dag(wd,udb+ip[2],wd+1);
               prod2su3alg(udb+ip[0],wd+1,&X);
               _su3_alg_mul_add_assign(*(fdb+ip[0]),r1,X);

               if ((t>0)||(bc!=1))
               {
                  _su3_alg_mul_sub_assign(*(fdb+ip[2]),r1,X);
               }

               if ((t<(N0-1))||(bc==3))
               {
                  prod2su3alg(wd+1,udb+ip[0],&X);
                  _su3_alg_mul_add_assign(*(fdb+ip[1]),r1,X);
               }
            }
         }
      }

      if ((t>0)||(bc!=1))
      {
         r0=c*c0;
         r1=c*c1;

         if ((t==0)&&(bc!=3))
         {
            r0*=(0.5*cG[0]);
            r1*=(0.5*cG[0]);
         }
         else if ((t==(N0-1))&&(bc==0))
         {
            r0*=(0.5*cG[1]);
            r1*=(0.5*cG[1]);
         }

         for (n=3;n<6;n++)
         {
            plaq_uidx(n,ix,ip);

            su3xsu3dag(udb+ip[1],udb+ip[3],wd);
            su3dagxsu3(udb+ip[2],udb+ip[0],wd+1);
            prod2su3alg(wd,wd+1,&X);
            _su3_alg_mul_add_assign(*(fdb+ip[1]),r0,X);

            prod2su3alg(wd+1,wd,&X);
            _su3_alg_mul_sub_assign(*(fdb+ip[3]),r0,X);

            su3xsu3dag(wd,udb+ip[2],wd+1);
            prod2su3alg(udb+ip[0],wd+1,&X);
            _su3_alg_mul_add_assign(*(fdb+ip[0]),r0,X);
            _su3_alg_mul_sub_assign(*(fdb+ip[2]),r0,X);

            if (c0!=1.0)
            {
               set_staples(n,ix,0);

               prod2su3alg(wd+1,vd,&X);
               _su3_alg_mul_add_assign(*(fdb+ip[1]),r1,X);

               prod2su3alg(vd,wd+1,&X);
               _su3_alg_mul_sub_assign(*(fdb+ip[2]),r1,X);

               su3dagxsu3(udb+ip[2],vd,wd+1);
               prod2su3alg(wd+1,wd,&X);
               _su3_alg_mul_sub_assign(*(fdb+ip[3]),r1,X);

               su3xsu3dag(udb+ip[3],vd+1,wd+1);
               su3xsu3dag(wd+1,udb+ip[0],wd+2);
               prod2su3alg(udb+ip[2],wd+2,&X);
               _su3_alg_mul_sub_assign(*(fdb+ip[0]),r1,X);
               _su3_alg_mul_add_assign(*(fdb+ip[2]),r1,X);

               prod2su3alg(wd+2,udb+ip[2],&X);
               _su3_alg_mul_add_assign(*(fdb+ip[3]),r1,X);

               su3xsu3dag(wd,vd+2,wd+1);
               prod2su3alg(udb+ip[0],wd+1,&X);
               _su3_alg_mul_add_assign(*(fdb+ip[0]),r1,X);

               prod2su3alg(wd+1,udb+ip[0],&X);
               _su3_alg_mul_add_assign(*(fdb+ip[1]),r1,X);

               su3dagxsu3(vd+2,udb+ip[0],wd+1);
               prod2su3alg(wd+1,wd,&X);
               _su3_alg_mul_sub_assign(*(fdb+ip[3]),r1,X);

               su3xsu3dag(udb+ip[1],vd+3,wd);
               su3xsu3dag(wd,udb+ip[2],wd+1);
               prod2su3alg(udb+ip[0],wd+1,&X);
               _su3_alg_mul_add_assign(*(fdb+ip[0]),r1,X);
               _su3_alg_mul_sub_assign(*(fdb+ip[2]),r1,X);

               prod2su3alg(wd+1,udb+ip[0],&X);
               _su3_alg_mul_add_assign(*(fdb+ip[1]),r1,X);
            }
         }
      }
   }

   add_bnd_frc();
}


static void wloops(int n,int ix,int t,double c0,double *trU)
{
   int bc,ip[4];

   bc=bc_type();
   plaq_uidx(n,ix,ip);

   trU[0]=0.0;
   trU[1]=0.0;
   trU[2]=0.0;
   trU[3]=0.0;

   if ((n>=3)||(t<(N0-1))||(bc!=0))
   {
      su3dagxsu3(udb+ip[2],udb+ip[0],wd);
      su3xsu3dag(udb+ip[1],udb+ip[3],wd+1);
      cm3x3_retr(wd,wd+1,trU);
      trU[0]=3.0-trU[0];
   }

   if (c0!=1.0)
   {
      set_staples(n,ix,1);

      if ((n<3)&&(((t==0)&&(bc==1))||
                  ((t==(N0-1))&&((bc==1)||(bc==2)))))
      {
         su3xsu3(wd,wd+1,wd+1);
         cm3x3_retr(wd+1,wd+1,trU+3);
         trU[3]=3.0-trU[3];
      }

      if ((n>=3)||(t<(N0-1))||(bc!=0))
      {
         su3xsu3dag(udb+ip[1],vd+3,wd+1);
         cm3x3_retr(wd,wd+1,trU+1);
         trU[1]=3.0-trU[1];
      }

      if ((n>=3)||(t<(N0-2))||((t==(N0-2))&&(bc!=0))||(bc==3))
      {
         su3xsu3dag(vd+1,udb+ip[3],wd+1);
         cm3x3_retr(wd,wd+1,trU+2);
         trU[2]=3.0-trU[2];
      }
   }
}


double action0(int icom)
{
   int bc,ix,t,n;
   double c0,c1,*cG;
   double r0,r1,trU[4],act;
   lat_parms_t lat;
   bc_parms_t bcp;

   lat=lat_parms();
   c0=lat.c0;
   c1=lat.c1;

   bcp=bc_parms();
   bc=bcp.type;
   cG=bcp.cG;

   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();
   udb=udfld();

   if (c0==1.0)
      hdb=NULL;
   else
   {
      if ((init&0x2)==0)
         set_ofs();

      if (query_flags(BSTAP_UP2DATE)!=1)
         set_bstap();
      hdb=bstap();
   }

   if ((init&0x1)==0)
   {
      ism=init_hsum(1);
      init+=1;
   }

   reset_hsum(ism);

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);
      act=0.0;

      if ((t<(N0-1))||(bc!=0))
      {
         r0=c0;

         if ((t==0)&&(bc==1))
            r0*=cG[0];
         else if ((t==(N0-1))&&(bc!=3))
            r0*=cG[1];

         for (n=0;n<3;n++)
         {
            wloops(n,ix,t,c0,trU);
            act+=(r0*trU[0]+c1*(trU[1]+trU[2]+0.5*trU[3]));
         }
      }

      if ((t>0)||(bc!=1))
      {
         r0=c0;
         r1=c1;

         if ((t==0)&&(bc!=3))
         {
            r0*=(0.5*cG[0]);
            r1*=(0.5*cG[0]);
         }
         else if ((t==(N0-1))&&(bc==0))
         {
            r0*=(0.5*cG[1]);
            r1*=(0.5*cG[1]);
         }

         for (n=3;n<6;n++)
         {
            wloops(n,ix,t,c0,trU);
            act+=(r0*trU[0]+r1*(trU[1]+trU[2]));
         }
      }

      add_to_hsum(ism,&act);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(ism,&act);
   else
      local_hsum(ism,&act);

   return (lat.beta/3.0)*act;
}
