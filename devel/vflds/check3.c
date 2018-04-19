
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2007, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the communication programs cpv_int_bnd() and cpv_ext_bnd().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "linalg.h"
#include "dfl.h"
#include "vflds.h"
#include "global.h"

static int bs[4],Ns,nv,nvec;
static int nb,nbb,*nbbe,*nbbo,*obbe,*obbo;
static int (*inn)[8],*ipp;


static void set_field(complex *v)
{
   int n[4],no[4],c[4];
   int i0,i1,i2,i3,ibe,ibo;

   n[0]=L0/bs[0];
   n[1]=L1/bs[1];
   n[2]=L2/bs[2];
   n[3]=L3/bs[3];

   no[0]=cpr[0]*n[0];
   no[1]=cpr[1]*n[1];
   no[2]=cpr[2]*n[2];
   no[3]=cpr[3]*n[3];

   set_v2zero(nv,v);
   ibe=0;
   ibo=(n[0]*n[1]*n[2]*n[3])/2;

   for (i0=0;i0<n[0];i0++)
   {
      for (i1=0;i1<n[1];i1++)
      {
         for (i2=0;i2<n[2];i2++)
         {
            for (i3=0;i3<n[3];i3++)
            {
               c[0]=no[0]+i0;
               c[1]=no[1]+i1;
               c[2]=no[2]+i2;
               c[3]=no[3]+i3;

               if (((c[0]+c[1]+c[2]+c[3])&0x1)==0x0)
               {
                  v[ibe*Ns  ].re=(float)(c[0]);
                  v[ibe*Ns+1].re=(float)(c[1]);
                  v[ibe*Ns+2].re=(float)(c[2]);
                  v[ibe*Ns+3].re=(float)(c[3]);
                  ibe+=1;
               }
               else
               {
                  v[ibo*Ns  ].re=(float)(c[0]);
                  v[ibo*Ns+1].re=(float)(c[1]);
                  v[ibo*Ns+2].re=(float)(c[2]);
                  v[ibo*Ns+3].re=(float)(c[3]);
                  ibo+=1;
               }
            }
         }
      }
   }
}


static void random_iv(int n,complex *v)
{
   complex *vm;

   random_v(n,v,100.0f);
   vm=v+n;

   for (;v<vm;v++)
   {
      (*v).re=(float)(floor((double)((*v).re+0.5f)));
      (*v).im=(float)(floor((double)((*v).im+0.5f)));
   }
}


static int chk_ext_bnd(complex *v)
{
   int bc,np[4];
   int ifc,ib,ibb,mu,i,ie;
   float c[4],n[4];

   bc=bc_type();

   np[0]=NPROC0;
   np[1]=NPROC1;
   np[2]=NPROC2;
   np[3]=NPROC3;

   n[0]=(float)((NPROC0*L0)/bs[0]);
   n[1]=(float)((NPROC1*L1)/bs[1]);
   n[2]=(float)((NPROC2*L2)/bs[2]);
   n[3]=(float)((NPROC3*L3)/bs[3]);
   ie=0;

   for (ifc=0;ifc<8;ifc++)
   {
      if ((ifc>1)||
          ((ifc==0)&&(cpr[0]!=0))||
          ((ifc==1)&&(cpr[0]!=(NPROC0-1)))||
          (bc==3))
      {
         for (ibb=obbe[ifc];ibb<(obbe[ifc]+nbbe[ifc]);ibb++)
         {
            ib=ipp[ibb];

            for (mu=0;mu<4;mu++)
            {
               c[mu]=v[nv+ibb*Ns+mu].re-v[ib*Ns+mu].re;

               if (mu==(ifc/2))
               {
                  if ((ifc&0x1)==0x0)
                  {
                     c[mu]+=1.0f;

                     if (cpr[mu]==0)
                        c[mu]-=n[mu];
                  }
                  else
                  {
                     c[mu]-=1.0f;

                     if (cpr[mu]==(np[mu]-1))
                        c[mu]+=n[mu];
                  }
               }
            }

            if ((c[0]!=0.0f)||(c[1]!=0.0f)||(c[2]!=0.0f)||(c[3]!=0.0f))
               ie=1;
         }
      }
      else
      {
         for (ibb=obbe[ifc];ibb<(obbe[ifc]+nbbe[ifc]);ibb++)
         {
            for (i=0;i<Ns;i++)
            {
               if ((v[nv+Ns*ibb+i].re!=0.0f)||(v[nv+Ns*ibb+i].im!=0.0f))
                  ie=2;
            }
         }
      }
   }

   return ie;
}


static int chk_int_bnd(complex *v,complex *w)
{
   int bc,ifc,ib,ibb,ie;
   complex *vv,*ww,*vm;

   bc=bc_type();

   for (ifc=0;ifc<8;ifc++)
   {
      if ((ifc>1)||
          ((ifc==0)&&(cpr[0]!=0))||
          ((ifc==1)&&(cpr[0]!=(NPROC0-1)))||
          (bc==3))
      {
         for (ibb=obbo[ifc];ibb<(obbo[ifc]+nbbo[ifc]);ibb++)
         {
            ib=ipp[ibb];
            vv=v+ib*Ns;
            ww=w+ib*Ns;
            vm=vv+Ns;

            for (;vv<vm;vv++)
            {
               (*vv).re-=(*ww).re;
               (*vv).im-=(*ww).im;
               ww+=1;
            }
         }
      }
   }

   vv=v;
   ww=w;
   vm=vv+nv;
   ie=0;

   for (;vv<vm;vv++)
   {
      if (((*vv).re!=(*ww).re)||((*vv).im!=(*ww).im))
         ie=1;
      ww+=1;
   }

   return ie;
}


int main(int argc,char *argv[])
{
   int my_rank,bc,i,ie;
   float d;
   double phi[2],phi_prime[2],theta[3];
   complex **wv,z;
   dfl_grid_t dfl_grid;
   FILE *fin=NULL,*flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Check of the communication programs cpv_int_bnd() "
             "and cpv_ext_bnd()\n");
      printf("--------------------------------------------------"
             "-----------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n\n",bs[0],bs[1],bs[2],bs[3]);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>]");
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(0);

   start_ranlux(0,123456);
   geometry();
   Ns=6;
   set_dfl_parms(bs,Ns);

   dfl_grid=dfl_geometry();
   nb=dfl_grid.nb;
   nbb=dfl_grid.nbb;
   nbbe=dfl_grid.nbbe;
   nbbo=dfl_grid.nbbo;
   obbe=dfl_grid.obbe;
   obbo=dfl_grid.obbo;
   inn=dfl_grid.inn;
   ipp=dfl_grid.ipp;

   alloc_wv(4);
   wv=reserve_wv(4);

   nv=Ns*nb;
   nvec=Ns*(nb+nbb/2);
   z.re=-1.0f;
   z.im=0.0f;

   for (i=0;i<2;i++)
   {
      random_v(nvec,wv[i],1.0f);
      set_field(wv[i]);
      assign_v2v(nv,wv[i],wv[i+1]);
      cpv_int_bnd(wv[i]);
      mulc_vadd(nv,wv[i+1],wv[i],z);
      d=vnorm_square(nv,1,wv[i+1]);

      error_root(d!=0.0f,1,"main [check3.c]",
                 "cpv_int_bnd() modifies the input field on the local grid");

      ie=chk_ext_bnd(wv[i]);
      error(ie==1,1,"main [check3.c]",
            "Boundary values are incorrectly mapped by cpv_int_bnd()");
      error(ie==2,1,"main [check3.c]",
            "Boundary values are not set to zero where they should");

      random_iv(nvec,wv[i]);
      cpv_int_bnd(wv[i]);
      assign_v2v(nvec,wv[i],wv[i+1]);
      cpv_ext_bnd(wv[i]);
      mulc_vadd(nvec-nv,wv[i]+nv,wv[i+1]+nv,z);
      d=vnorm_square(nvec-nv,1,wv[i]+nv);

      error_root(d!=0.0f,1,"main [check3.c]",
                 "cpv_ext_bnd() modifies the input field on the boundary");

      ie=chk_int_bnd(wv[i],wv[i+1]);
      error(ie==1,1,"main [check3.c]",
            "Boundary values are incorrectly mapped by cpv_ext_bnd()");
   }

   if (my_rank==0)
   {
      printf("No errors detected\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
