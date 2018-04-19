
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2010-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Basic checks on the implementation of the Wilson flow.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "flags.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "mdflds.h"
#include "linalg.h"
#include "forces.h"
#include "wflow.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static su3_alg_dble XX ALIGNED16;
static su3_dble mm ALIGNED16;
static su3_dble uu ALIGNED16;
static su3_dble vv ALIGNED16;


static double cmp_ud(su3_dble *u,su3_dble *v)
{
   int i;
   double r[18],dev,dmax;

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

   dmax=0.0;

   for (i=0;i<18;i+=2)
   {
      dev=r[i]*r[i]+r[i+1]*r[i+1];
      if (dev>dmax)
         dmax=dev;
   }

   return sqrt(dmax);
}


static double max_dev_ud(su3_dble *v)
{
   double d,dmax;
   su3_dble *u,*um;

   u=udfld();
   um=u+4*VOLUME;
   dmax=0.0;

   for (;u<um;u++)
   {
      d=cmp_ud(u,v);

      if (d>dmax)
         dmax=d;

      v+=1;
   }

   if (NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return dmax;
}


static double cmp_fd(su3_alg_dble *f,su3_alg_dble *g)
{
   int i;
   double r[8],dev,dmax;

   r[0]=(*f).c1-(*g).c1;
   r[1]=(*f).c2-(*g).c2;
   r[2]=(*f).c3-(*g).c3;
   r[3]=(*f).c4-(*g).c4;
   r[4]=(*f).c5-(*g).c5;
   r[5]=(*f).c6-(*g).c6;
   r[6]=(*f).c7-(*g).c7;
   r[7]=(*f).c8-(*g).c8;

   dmax=0.0;

   for (i=0;i<8;i++)
   {
      dev=fabs(r[i]);
      if (dev>dmax)
         dmax=dev;
   }

   return dmax;
}


static double max_dev_frc(su3_alg_dble *g)
{
   double d,dmax;
   su3_alg_dble *f,*fm;
   mdflds_t *mdfs;

   mdfs=mdflds();
   f=(*mdfs).frc;
   fm=f+4*VOLUME;
   dmax=0.0;

   for (;f<fm;f++)
   {
      d=cmp_fd(f,g);

      if (d>dmax)
         dmax=d;

      g+=1;
   }

   if (NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return dmax;
}


static int is_zero(su3_alg_dble *X)
{
   int ie;

   ie=((*X).c1==0.0);
   ie&=((*X).c2==0.0);
   ie&=((*X).c3==0.0);
   ie&=((*X).c4==0.0);
   ie&=((*X).c5==0.0);
   ie&=((*X).c6==0.0);
   ie&=((*X).c7==0.0);
   ie&=((*X).c8==0.0);

   return ie;
}


static int check_bnd_fld(su3_alg_dble *fld)
{
   int bc,npts,*pts,*ptm;
   int ix,t,ifc,ie;
   su3_alg_dble *f;

   bc=bc_type();
   pts=bnd_pts(&npts);
   ptm=pts+npts;
   pts+=(npts/2);
   ie=0;

   for (;pts<ptm;pts++)
   {
      ix=pts[0];
      t=global_time(ix);
      f=fld+8*(ix-(VOLUME/2));

      ie|=((t!=0)&&(t!=(N0-1)));
      ie|=((t==(N0-1))&&(bc!=0));

      if (bc==0)
      {
         if (t==0)
         {
            ie|=is_zero(f);
            ie|=(is_zero(f+1)^0x1);
         }
         else
         {
            ie|=(is_zero(f)^0x1);
            ie|=is_zero(f+1);
         }

         for (ifc=2;ifc<8;ifc++)
            ie|=is_zero(f+ifc);
      }
      else if (bc==1)
      {
         ie|=is_zero(f);
         ie|=is_zero(f+1);

         for (ifc=2;ifc<8;ifc++)
            ie|=(is_zero(f+ifc)^0x1);
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
            ie|=is_zero(f+ifc);
      }
   }

   return ie;
}


static int ofs(int ix,int mu)
{
   int iy;

   if (ix<(VOLUME/2))
   {
      iy=iup[ix][mu];

      return 8*(iy-(VOLUME/2))+2*mu+1;
   }
   else
      return 8*(ix-(VOLUME/2))+2*mu;
}


static double chkfrc(void)
{
   int x0,x1,x2,x3;
   int ix,iy,iz,iw,mu,nu;
   double d,dmax;
   su3_alg_dble *frc;
   su3_dble *udb;
   mdflds_t *mdfs;

   udb=udfld();
   mdfs=mdflds();
   dmax=0.0;

   for (x0=1;x0<(L0-2);x0++)
   {
      for (x1=1;x1<(L1-2);x1++)
      {
         for (x2=1;x2<(L2-2);x2++)
         {
            for (x3=1;x3<(L3-2);x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

               for (mu=0;mu<4;mu++)
               {
                  cm3x3_zero(1,&mm);
                  iy=iup[ix][mu];

                  for (nu=0;nu<4;nu++)
                  {
                     if (nu!=mu)
                     {
                        iz=iup[ix][nu];

                        su3xsu3dag(udb+ofs(iy,nu),udb+ofs(iz,mu),&uu);
                        su3xsu3dag(&uu,udb+ofs(ix,nu),&vv);
                        cm3x3_add(&vv,&mm);

                        iz=idn[ix][nu];
                        iw=idn[iy][nu];

                        su3dagxsu3(udb+ofs(iz,mu),udb+ofs(iz,nu),&uu);
                        su3dagxsu3(udb+ofs(iw,nu),&uu,&vv);
                        cm3x3_add(&vv,&mm);
                     }
                  }

                  prod2su3alg(udb+ofs(ix,mu),&mm,&XX);

                  if (ix<(VOLUME/2))
                     frc=(*mdfs).frc+8*(iy-(VOLUME/2))+2*mu+1;
                  else
                     frc=(*mdfs).frc+8*(ix-(VOLUME/2))+2*mu;

                  d=cmp_fd(&XX,frc);

                  if (d>dmax)
                     dmax=d;
               }
            }
         }
      }
   }

   if (NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return dmax;
}


static void scale_bnd_frc(su3_alg_dble *frc)
{
   int bc,ifc,npts,*pts,*ptm;
   su3_alg_dble *fr;

   bc=bc_type();

   if ((bc==0)||(bc==2))
   {
      pts=bnd_pts(&npts);
      ptm=pts+npts;
      pts+=(npts/2);

      for (;pts<ptm;pts++)
      {
         fr=frc+8*(pts[0]-(VOLUME/2));

         for (ifc=2;ifc<8;ifc++)
         {
            _su3_alg_mul_assign(fr[ifc],2.0);
         }
      }
   }
}


int main(int argc,char *argv[])
{
   int my_rank,bc,n,k,ie;
   double phi[2],phi_prime[2],theta[3];
   double eps,nplaq,act0,act1,dev0,dev1;
   su3_dble *udb,*u,*um,**usv;
   su3_alg_dble *frc,**fsv;
   mdflds_t *mdfs;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Basic checks on the implementation of the Wilson flow\n");
      printf("-----------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("n","%d\n",&n);
      read_line("eps","%lf",&eps);
      fclose(fin);

      printf("n = %d\n",n);
      printf("eps = %.3e\n\n",eps);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>]");
   }

   MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&eps,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);

   set_lat_parms(6.0,1.0,0,NULL,1.0);
   print_lat_parms();

   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(0);

   start_ranlux(0,1234);
   geometry();
   alloc_wud(1);
   alloc_wfd(2);
   mdfs=mdflds();
   usv=reserve_wud(1);
   fsv=reserve_wfd(1);
   udb=udfld();

   if (bc==0)
      nplaq=(double)(6*N0-6)*(double)(N1*N2*N3);
   else
      nplaq=(double)(6*N0)*(double)(N1*N2*N3);

   random_ud();
   act0=action0(1);
   act1=3.0*nplaq-plaq_wsum_dble(1);

   plaq_frc();
   ie=check_bnd_fld((*mdfs).frc);
   error(ie!=0,1,"main [check1.c]",
         "Force vanishes on an incorrect subset of links");
   assign_alg2alg(4*VOLUME,(*mdfs).frc,fsv[0]);
   force0(1.0);
   ie=check_bnd_fld((*mdfs).frc);
   error(ie!=0,1,"main [check1.c]",
         "Force vanishes on an incorrect subset of links");
   dev0=max_dev_frc(fsv[0]);

   if (my_rank==0)
   {
      printf("Random gauge field:\n");
      printf("Action (action0)   = %.15e\n",act0);
      printf("Action (plaq_wsum) = %.15e\n",2.0*act1);
      printf("Deviation of force = %.1e\n\n",dev0);
   }

   random_ud();
   cm3x3_assign(4*VOLUME,udb,usv[0]);
   plaq_frc();
   assign_alg2alg(4*VOLUME,(*mdfs).frc,fsv[0]);
   dev0=chkfrc();
   fwd_euler(1,eps);
   frc=fsv[0];
   scale_bnd_frc(frc);
   u=udb;
   um=u+4*VOLUME;

   for (;u<um;u++)
   {
      if (is_zero(frc)==0)
         expXsu3(eps,frc,u);
      frc+=1;
   }

   set_flags(UPDATED_UD);
   dev1=max_dev_ud(usv[0]);

   if (my_rank==0)
   {
      printf("Direct check of the generator: %.1e\n",dev0);
      printf("Check of the 1-step integration: %.1e\n\n",dev1);

      printf("Evolution of the Wilson action:\n\n");
   }

   random_ud();
   act0=3.0*nplaq-plaq_wsum_dble(1);

   if (my_rank==0)
      printf("k =  0: %.8e\n",2.0*act0/nplaq);

   for (k=1;k<=n;k++)
   {
      fwd_euler(1,eps);
      act1=3.0*nplaq-plaq_wsum_dble(1);

      error(((act1>act0)&&(eps>=0.0))||((act1<act0)&&(eps<=0.0)),1,
            "main [check1.c]","The Wilson action is not monotonic");

      act0=act1;

      if (my_rank==0)
         printf("k = %2d: %.8e\n",k,2.0*act0/nplaq);
   }

   ie=check_bc(0.0);
   error_root(ie!=1,1,"main [check3.c]",
              "Boundary values of the gauge field are not preserved");

   if (my_rank==0)
   {
      printf("\n");
      printf("Monotonicity check passed\n\n");
      fflush(stdout);
   }

   start_ranlux(0,1234);
   random_ud();
   fwd_euler(n,eps);
   ie=check_bc(0.0);
   error_root(ie!=1,1,"main [check3.c]",
              "Boundary values of the gauge field are not preserved");
   cm3x3_assign(4*VOLUME,udb,usv[0]);

   start_ranlux(0,1234);
   random_ud();
   fwd_rk2(n,eps);
   ie=check_bc(0.0);
   error_root(ie!=1,1,"main [check3.c]",
              "Boundary values of the gauge field are not preserved");
   dev0=max_dev_ud(usv[0]);
   cm3x3_assign(4*VOLUME,udb,usv[0]);

   if (my_rank==0)
      printf("Comparison of fwd_euler() and fwd_rk2(): |dU| = %.1e\n",
             dev0);

   start_ranlux(0,1234);
   random_ud();
   fwd_rk3(n,eps);
   ie=check_bc(0.0);
   error_root(ie!=1,1,"main [check3.c]",
              "Boundary values of the gauge field are not preserved");
   dev0=max_dev_ud(usv[0]);

   if (my_rank==0)
   {
      printf("Comparison of fwd_rk2() and fwd_rk3():   |dU| = %.1e\n\n",
             dev0);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
