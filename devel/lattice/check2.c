
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2010, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs set_bc() and check_bc().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "uflds.h"
#include "lattice.h"
#include "global.h"

#define N0 (NPROC0*L0)


static void new_fld(int ibnd)
{
   su3_dble *ud,*udm;

   ud=udfld();
   udm=ud+4*VOLUME;

   for (;ud<udm;ud++)
      random_su3_dble(ud);

   if (ibnd)
   {
      udm+=7*(BNDRY/4);

      if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
         udm+=3;

      for (;ud<udm;ud++)
         random_su3_dble(ud);
   }

   set_flags(UPDATED_UD);
   set_flags(UNSET_UD_PHASE);
}


static int cmp_ud(su3_dble *u,su3_dble *v)
{
   int i;
   double r[18];

   r[ 0]=(*u).c11.re-(*v).c11.re;
   r[ 1]=(*u).c11.im-(*v).c11.im;
   r[ 2]=(*u).c12.re-(*v).c12.re;
   r[ 3]=(*u).c12.im-(*v).c12.im;
   r[ 4]=(*u).c13.re-(*v).c13.re;
   r[ 5]=(*u).c13.im-(*v).c13.im;

   r[ 6]=(*u).c21.re-(*v).c21.re;
   r[ 7]=(*u).c21.im-(*v).c21.im;
   r[ 8]=(*u).c22.re-(*v).c22.re;
   r[ 9]=(*u).c22.im-(*v).c22.im;
   r[10]=(*u).c23.re-(*v).c23.re;
   r[11]=(*u).c23.im-(*v).c23.im;

   r[12]=(*u).c31.re-(*v).c31.re;
   r[13]=(*u).c31.im-(*v).c31.im;
   r[14]=(*u).c32.re-(*v).c32.re;
   r[15]=(*u).c32.im-(*v).c32.im;
   r[16]=(*u).c33.re-(*v).c33.re;
   r[17]=(*u).c33.im-(*v).c33.im;

   for (i=0;i<18;i++)
   {
      if (r[i]!=0.0)
         return 1;
   }

   return 0;
}


static int cmp_active(su3_dble *u,su3_dble *v)
{
   int ix,t,ifc,ie,bc;

   bc=bc_type();
   ie=0;

   for (ix=(VOLUME/2);ix<(VOLUME);ix++)
   {
      t=global_time(ix);

      for (ifc=0;ifc<8;ifc++)
      {
         if (((t>0)&&(t<(N0-1)))||
             ((t==0)&&((ifc==0)||((ifc==1)&&(bc!=0))||((ifc>=2)&&(bc!=1))))||
             ((t==(N0-1))&&(bc!=0)))
            ie|=cmp_ud(u,v);

         u+=1;
         v+=1;
      }
   }

   return ie;
}


static int check_diag(su3_dble *u)
{
   int i,ie;
   double r[18];
   complex_dble z;

   ie=0;

   r[ 0]=(*u).c11.re;
   r[ 1]=(*u).c11.im;
   r[ 2]=(*u).c12.re;
   r[ 3]=(*u).c12.im;
   r[ 4]=(*u).c13.re;
   r[ 5]=(*u).c13.im;

   r[ 6]=(*u).c21.re;
   r[ 7]=(*u).c21.im;
   r[ 8]=(*u).c22.re;
   r[ 9]=(*u).c22.im;
   r[10]=(*u).c23.re;
   r[11]=(*u).c23.im;

   r[12]=(*u).c31.re;
   r[13]=(*u).c31.im;
   r[14]=(*u).c32.re;
   r[15]=(*u).c32.im;
   r[16]=(*u).c33.re;
   r[17]=(*u).c33.im;

   ie|=(fabs(r[ 0]*r[ 0]+r[ 1]*r[ 1]-1.0)>(8.0*DBL_EPSILON));
   ie|=(fabs(r[ 8]*r[ 8]+r[ 9]*r[ 9]-1.0)>(8.0*DBL_EPSILON));
   ie|=(fabs(r[16]*r[16]+r[17]*r[17]-1.0)>(8.0*DBL_EPSILON));

   z.re=r[0]*r[8]-r[1]*r[9];
   z.im=r[0]*r[9]+r[1]*r[8];
   ie|=(fabs(z.re*r[16]-z.im*r[17]-1.0)>(16.0*DBL_EPSILON));
   ie|=(fabs(z.re*r[17]+z.im*r[16])>(16.0*DBL_EPSILON));

   for (i=0;i<18;i++)
   {
      if (((i>1)&&(i<8))||((i>9)&&(i<16)))
         ie|=(r[i]!=0.0);
   }

   return ie;
}


static int check_bval(su3_dble *u)
{
   int bc,ie,ifc;
   int ipt,npts,*pts;

   ie=0;
   bc=bc_type();

   if (bc==1)
   {
      pts=bnd_pts(&npts);

      if (npts>0)
      {
         pts+=(npts/2);

         ie|=check_diag(u+8*(pts[0]-(VOLUME/2))+2);
         ie|=check_diag(u+8*(pts[0]-(VOLUME/2))+4);
         ie|=check_diag(u+8*(pts[0]-(VOLUME/2))+6);

         for (ipt=0;ipt<(npts/2);ipt++)
         {
            for (ifc=2;ifc<8;ifc++)
               ie|=cmp_ud(u+8*(pts[0]-(VOLUME/2))+2*(ifc/2),
                          u+8*(pts[ipt]-(VOLUME/2))+ifc);
         }
      }
   }

   if (((bc==1)||(bc==2))&&(cpr[0]==(NPROC0-1)))
   {
      u+=4*VOLUME+7*(BNDRY/4);

      ie|=check_diag(u);
      ie|=check_diag(u+1);
      ie|=check_diag(u+2);
   }

   return ie;
}


int main(int argc,char *argv[])
{
   int my_rank,bc,ie;
   double phi[2],phi_prime[2],theta[3];
   double cG,cG_prime,cF,cF_prime;
   su3_dble *udb,**usv;
   bc_parms_t bcp;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      printf("\n");
      printf("Check of set_bc() and check_bc()\n");
      printf("--------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   cG=0.97;
   cG_prime=1.056;
   cF=0.82;
   cF_prime=1.12;
   theta[0]=0.38;
   theta[1]=-1.25;
   theta[2]=0.54;
   set_bc_parms(bc,cG,cG_prime,cF,cF_prime,phi,phi_prime,theta);
   print_bc_parms(3);

   start_ranlux(0,12345);
   geometry();
   alloc_wud(1);
   usv=reserve_wud(1);
   udb=udfld();

   ie=0;
   bcp=bc_parms();
   error(bcp.type!=bc,1,"main [check2.c]",
         "Type of boundary condition is not properly set");

   if (bc!=3)
   {
      ie|=(cG!=bcp.cG[0]);
      ie|=(cF!=bcp.cF[0]);
   }

   if (bc<=1)
   {
      ie|=(bcp.cG[0]!=bcp.cG[1]);
      ie|=(bcp.cF[0]!=bcp.cF[1]);
   }

   if (bc==2)
   {
      ie|=(cG_prime!=bcp.cG[1]);
      ie|=(cF_prime!=bcp.cF[1]);
   }

   if (bc==1)
   {
      ie|=(phi[0]!=bcp.phi[0][0]);
      ie|=(phi[1]!=bcp.phi[0][1]);
      ie|=(bcp.phi[0][2]!=-bcp.phi[0][0]-bcp.phi[0][1]);
   }

   if ((bc==1)||(bc==2))
   {
      ie|=(phi_prime[0]!=bcp.phi[1][0]);
      ie|=(phi_prime[1]!=bcp.phi[1][1]);
      ie|=(bcp.phi[1][2]!=-bcp.phi[1][0]-bcp.phi[1][1]);
   }

   ie|=(bcp.theta[0]!=theta[0]);
   ie|=(bcp.theta[1]!=theta[1]);
   ie|=(bcp.theta[2]!=theta[2]);

   error(ie,1,"main [check2.c]","Boundary parameters are not properly set");

   ie=check_bc(0.0);
   error(ie!=1,1,"main [check2.c]",
         "check_bc() gives the wrong answer");

   new_fld(0);
   ie=check_bc(0.0);
   error(((bc<2)&&(ie!=0))||((bc>=2)&&(ie!=1)),2,"main [check2.c]",
         "check_bc() gives the wrong answer");

   new_fld(1);
   ie=check_bc(0.0);
   error(((bc<3)&&(ie!=0))||((bc==3)&&(ie!=1)),2,"main [check2.c]",
         "check_bc() gives the wrong answer");

   cm3x3_assign(4*VOLUME,udb,usv[0]);
   set_bc();
   ie=check_bc(0.0);
   error(ie!=1,2,"main [check2.c]",
         "check_bc() gives the wrong answer");

   ie=cmp_active(udb,usv[0]);
   error(ie!=0,2,"main [check2.c]",
         "Active link variables are modified by set_bc()");

   ie=check_bval(udb);
   error(ie!=0,2,"main [check2.c]",
         "Boundary values are not properly set by set_bc()");

   if (my_rank==0)
   {
      printf("No errors detected --- all programs work correctly\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
