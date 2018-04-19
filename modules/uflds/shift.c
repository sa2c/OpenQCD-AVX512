
/*******************************************************************************
*
* File shift.c
*
* Copyright (C) 2006, 2009, 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Translation of the global double-precision gauge field.
*
* The externally accessible function is
*
*   int shift_ud(int *s)
*     Replaces the double-precision gauge field U(x,mu) by U(x-s,mu), where
*     s[4] is any given shift vector. The program returns the number of
*     elementary steps (translations by 1 lattice unit) that were performed.
*
* Notes:
*
* Shifts in the time direction are only permitted in the case of periodic
* boundary conditions. The required communication buffers are allocated
* automatically.
*
* The program shift_ud() acts globally and must be called simultaneously on
* all MPI processes with the same translation vector.
*
*******************************************************************************/

#define SHIFT_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "flags.h"
#include "lattice.h"
#include "uflds.h"
#include "global.h"

typedef struct
{
   int saddr,raddr;
   int iu,*idx;
} comlink_t;

typedef struct
{
   int saddr,raddr;
   int **idx;
} comstar_t;

static int init=0,bs[4],np[4],nlk[4],npt[4],ofs[4];
static su3_dble *sdbuf=NULL,*rdbuf;
static comlink_t comlink[4];
static comstar_t comstar[8];


static void set_const(void)
{
   int mu,ifc,iu;
   comlink_t *cl;
   comstar_t *cs;

   bs[0]=L0;
   bs[1]=L1;
   bs[2]=L2;
   bs[3]=L3;

   np[0]=NPROC0;
   np[1]=NPROC1;
   np[2]=NPROC2;
   np[3]=NPROC3;

   nlk[0]=FACE0/2;
   nlk[1]=FACE1/2;
   nlk[2]=FACE2/2;
   nlk[3]=FACE3/2;

   ofs[0]=FACE0/2;
   ofs[1]=ofs[0]+(FACE0+FACE1)/2;
   ofs[2]=ofs[1]+(FACE1+FACE2)/2;
   ofs[3]=ofs[2]+(FACE2+FACE3)/2;

   npt[0]=VOLUME/L0;
   npt[1]=VOLUME/L1;
   npt[2]=VOLUME/L2;
   npt[3]=VOLUME/L3;

   cl=comlink;
   iu=4*VOLUME;

   for (mu=0;mu<4;mu++)
   {
      (*cl).saddr=npr[2*mu];
      (*cl).raddr=npr[2*mu+1];
      (*cl).iu=iu;
      iu+=nlk[mu];
      cl+=1;
   }

   cs=comstar;

   for (ifc=0;ifc<8;ifc++)
   {
      (*cs).saddr=npr[ifc];
      (*cs).raddr=npr[ifc^0x1];
      cs+=1;
   }
}


static void alloc_idx(void)
{
   int mu,ifc,t,**id,*ib;
   comlink_t *cl;
   comstar_t *cs;

   cl=comlink;

   if (NPROC>1)
   {
      ib=malloc((BNDRY/4)*sizeof(*ib));
      error(ib==NULL,1,"alloc_idx [shift.c]",
            "Unable to allocate index arrays");
   }
   else
      ib=NULL;

   for (mu=0;mu<4;mu++)
   {
      if (np[mu]>1)
      {
         (*cl).idx=ib;
         ib+=nlk[mu];
      }
      else
         (*cl).idx=NULL;

      cl+=1;
   }

   cs=comstar;
   id=malloc(2*(L0+L1+L2+L3)*sizeof(*id));
   ib=malloc(16*VOLUME*sizeof(*ib));
   error((id==NULL)||(ib==NULL),1,"alloc_idx [shift.c]",
         "Unable to allocate index arrays");

   for (ifc=0;ifc<8;ifc++)
   {
      mu=ifc/2;
      (*cs).idx=id;
      id+=bs[mu];

      if ((ifc&0x1)==0)
      {
         for (t=0;t<bs[mu];t++)
         {
            (*cs).idx[t]=ib;
            ib+=4*npt[mu];
         }
      }
      else
      {
         for (t=0;t<bs[mu];t++)
            (*cs).idx[t]=(*(cs-1)).idx[bs[mu]-1-t];
      }

      cs+=1;
   }
}


static void fce_pts(void)
{
   int x0,x1,x2,x3,x[4],t;
   int mu,ix,iy,*idx[4];
   int *fce0,*fce1,*fcem;

   for (mu=0;mu<4;mu++)
      idx[mu]=(comstar[2*mu]).idx[1];

   ix=0;

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               x[0]=x0;
               x[1]=x1;
               x[2]=x2;
               x[3]=x3;
               iy=ipt[ix];
               ix+=1;

               for (mu=0;mu<4;mu++)
               {
                  if (x[mu]==0)
                     idx[mu][iy]=1;
                  else
                     idx[mu][iy]=0;
               }
            }
         }
      }
   }

   for (mu=0;mu<4;mu++)
   {
      fce0=(comstar[2*mu]).idx[0];

      for (ix=0;ix<VOLUME;ix++)
      {
         if (idx[mu][ix]==1)
         {
            *fce0=ix;
            fce0+=4;
         }
      }

      for (t=1;t<bs[mu];t++)
      {
         fce0=(comstar[2*mu]).idx[t-1];
         fce1=(comstar[2*mu]).idx[t];
         fcem=fce0+4*npt[mu];

         for (;fce0<fcem;fce0+=4,fce1+=4)
            *fce1=iup[*fce0][mu];
      }
   }
}


static int uofs(int ix,int mu)
{
   int iy,iz;

   if (ix>=(VOLUME/2))
      return 8*(ix-(VOLUME/2))+2*mu;

   iy=iup[ix][mu];

   if (iy<VOLUME)
      return 8*(iy-(VOLUME/2))+2*mu+1;

   iz=iy-(VOLUME+(BNDRY/2)+ofs[mu]);

   return comlink[mu].iu+iz;
}


static void set_idx(void)
{
   int mu,ix,iy,iz,t;
   int *idx,*idm;

   for (mu=0;mu<4;mu++)
   {
      if (np[mu]>1)
      {
         idx=(comlink[mu]).idx;

         for (ix=0;ix<nlk[mu];ix++)
         {
            iy=(BNDRY/2)+ofs[mu]+ix;
            iz=map[iy];
            (*idx)=8*(iz-(VOLUME/2))+2*mu+1;
            idx+=1;
         }
      }
   }

   fce_pts();

   for (mu=0;mu<4;mu++)
   {
      for (t=0;t<bs[mu];t++)
      {
         idx=(comstar[2*mu]).idx[t];
         idm=idx+4*npt[mu];

         for (;idx<idm;idx+=4)
         {
            ix=idx[0];
            idx[0]=uofs(ix,0);
            idx[1]=uofs(ix,1);
            idx[2]=uofs(ix,2);
            idx[3]=uofs(ix,3);
         }
      }
   }
}


static void alloc_udbufs(void)
{
   int mu,n;

   if (init==0)
   {
      set_const();
      alloc_idx();
      set_idx();
      init=1;
   }

   n=0;

   for (mu=0;mu<4;mu++)
   {
      if (npt[mu]>n)
         n=npt[mu];
   }

   sdbuf=amalloc(8*n*sizeof(su3_dble),ALIGN);
   error(sdbuf==NULL,1,"alloc_udbufs [shift.c]",
         "Unable to allocate communication buffers");

   rdbuf=sdbuf+4*n;
}


static void get_udlinks(void)
{
   int mu,*idx,*idm;
   int tag,ip,saddr,raddr,nbf;
   su3_dble *ub,*u,*sb,*rb;
   comlink_t *cl;
   MPI_Status stat;

   ub=udfld();
   ip=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   cl=comlink;

   for (mu=0;mu<4;mu++)
   {
      if (np[mu]>1)
      {
         u=sdbuf;
         idx=(*cl).idx;
         idm=idx+nlk[mu];

         for (;idx<idm;idx++,u++)
            *u=ub[*idx];

         tag=mpi_tag();
         nbf=18*nlk[mu];
         saddr=(*cl).saddr;
         raddr=(*cl).raddr;
         sb=sdbuf;
         rb=ub+(*cl).iu;

         if (ip==0)
         {
            MPI_Send(sb,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rb,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rb,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sb,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }
      }

      cl+=1;
   }
}


static void put_udlinks(void)
{
   int mu,*idx,*idm;
   int tag,ip,saddr,raddr,nbf;
   su3_dble *ub,*u,*sb,*rb;
   comlink_t *cl;
   MPI_Status stat;

   ub=udfld();
   ip=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   cl=comlink;

   for (mu=0;mu<4;mu++)
   {
      if (np[mu]>1)
      {
         tag=mpi_tag();
         nbf=18*nlk[mu];
         saddr=(*cl).raddr;
         raddr=(*cl).saddr;
         sb=ub+(*cl).iu;
         rb=sdbuf;

         if (ip==0)
         {
            MPI_Send(sb,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rb,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rb,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sb,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }

         u=sdbuf;
         idx=(*cl).idx;
         idm=idx+nlk[mu];

         for (;idx<idm;idx++,u++)
            ub[*idx]=*u;
      }

      cl+=1;
   }
}


static void get_udstars(int ifc)
{
   int mu,*idx,*idm;
   int tag,ip,saddr,raddr,nbf;
   su3_dble *ub,*u,*sb,*rb;
   comstar_t *cs;
   MPI_Status stat;

   ub=udfld();
   ip=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   cs=comstar+ifc;
   mu=ifc/2;

   if (np[mu]>1)
      u=sdbuf;
   else
      u=rdbuf;

   idx=(*cs).idx[0];
   idm=idx+4*npt[mu];

   for (;idx<idm;idx++,u++)
      *u=ub[*idx];

   if (np[mu]>1)
   {
      tag=mpi_tag();
      nbf=72*npt[mu];
      saddr=(*cs).saddr;
      raddr=(*cs).raddr;
      sb=sdbuf;
      rb=rdbuf;

      if (ip==0)
      {
         MPI_Send(sb,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         MPI_Recv(rb,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
      }
      else
      {
         MPI_Recv(rb,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         MPI_Send(sb,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
      }
   }
}


static void shift_udstars(int ifc)
{
   int mu,t;
   int *id0,*id1,*idm;
   su3_dble *ub,*u;
   comstar_t *cs;

   get_udstars(ifc);

   ub=udfld();
   cs=comstar+ifc;
   mu=ifc/2;

   for (t=0;t<(bs[mu]-1);t++)
   {
      id0=(*cs).idx[t];
      id1=(*cs).idx[t+1];
      idm=id0+4*npt[mu];

      for (;id0<idm;id0++,id1++)
         ub[*id0]=ub[*id1];
   }

   id0=(*cs).idx[bs[mu]-1];
   idm=id0+4*npt[mu];
   u=rdbuf;

   for (;id0<idm;id0++,u++)
      ub[*id0]=*u;
}


int shift_ud(int *s)
{
   int iprms[4],sr[4];
   int mu,ifc,t,n;

   if (NPROC>1)
   {
      iprms[0]=s[0];
      iprms[1]=s[1];
      iprms[2]=s[2];
      iprms[3]=s[3];

      MPI_Bcast(iprms,4,MPI_INT,0,MPI_COMM_WORLD);

      error((iprms[0]!=s[0])||(iprms[1]!=s[1])||
            (iprms[2]!=s[2])||(iprms[3]!=s[3]),1,
            "shift_ud [shift.c]","Shift vector is not global");
   }

   if (sdbuf==NULL)
      alloc_udbufs();

   for (mu=0;mu<4;mu++)
   {
      n=np[mu]*bs[mu];

      if (abs(s[mu])>(n/2))
      {
         sr[mu]=safe_mod(s[mu],n);

         if (sr[mu]>(n/2))
            sr[mu]-=n;
      }
      else
         sr[mu]=s[mu];
   }

   if ((sr[0]==0)&&(sr[1]==0)&&(sr[2]==0)&&(sr[3]==0))
      return 0;

   error_root((sr[0]!=0)&&(bc_type()!=3),1,"shift_ud [shift.c]",
              "Shifts in time are only permitted for periodic bc");

   error_root((sr[0]!=0)&&(query_flags(UD_PHASE_SET)!=0),1,"shift_ud [shift.c]",
              "Attempt to move sign-changed link variables in time");

   get_udlinks();
   n=0;

   for (mu=0;mu<4;mu++)
   {
      if (sr[mu]>=0)
         ifc=2*mu+1;
      else
         ifc=2*mu;

      for (t=0;t<abs(sr[mu]);t++)
      {
         shift_udstars(ifc);
         n+=1;
      }
   }

   put_udlinks();

   set_flags(SHIFTED_UD);

   return n;
}
