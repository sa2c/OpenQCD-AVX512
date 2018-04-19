
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2005, 2008-2013, 2016 Martin Luescher, Filippo Palombi,
*                                     Stefan Schaefer
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of sw_frc() and hop_frc().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "mdflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "forces.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define MAX_LEVELS 8
#define BLK_LENGTH 8

static int cnt[MAX_LEVELS];
static double smx[MAX_LEVELS];


static int is_Xt_zero(u3_alg_dble *X)
{
   int ie;

   ie=1;
   ie&=((*X).c1==0.0);
   ie&=((*X).c2==0.0);
   ie&=((*X).c3==0.0);
   ie&=((*X).c4==0.0);
   ie&=((*X).c5==0.0);
   ie&=((*X).c6==0.0);
   ie&=((*X).c7==0.0);
   ie&=((*X).c8==0.0);
   ie&=((*X).c9==0.0);

   return ie;
}


static int is_Xv_zero(su3_dble *X)
{
   int ie;

   ie=1;
   ie&=((*X).c11.re==0.0);
   ie&=((*X).c11.im==0.0);
   ie&=((*X).c12.re==0.0);
   ie&=((*X).c12.im==0.0);
   ie&=((*X).c13.re==0.0);
   ie&=((*X).c13.im==0.0);

   ie&=((*X).c21.re==0.0);
   ie&=((*X).c21.im==0.0);
   ie&=((*X).c22.re==0.0);
   ie&=((*X).c22.im==0.0);
   ie&=((*X).c23.re==0.0);
   ie&=((*X).c23.im==0.0);

   ie&=((*X).c31.re==0.0);
   ie&=((*X).c31.im==0.0);
   ie&=((*X).c32.re==0.0);
   ie&=((*X).c32.im==0.0);
   ie&=((*X).c33.re==0.0);
   ie&=((*X).c33.im==0.0);

   return ie;
}


static int is_frc_zero(su3_alg_dble *f)
{
   int ie;

   ie=1;
   ie&=((*f).c1==0.0);
   ie&=((*f).c2==0.0);
   ie&=((*f).c3==0.0);
   ie&=((*f).c4==0.0);
   ie&=((*f).c5==0.0);
   ie&=((*f).c6==0.0);
   ie&=((*f).c7==0.0);
   ie&=((*f).c8==0.0);

   return ie;
}


static void check_Xtbnd(ptset_t set)
{
   int bc,ix,t,n,ie;
   int ia,ib;
   u3_alg_dble **xt;

   bc=bc_type();
   xt=xtensor();
   ie=0;
   ia=0;
   ib=VOLUME;

   if (set==EVEN_PTS)
      ib=(VOLUME/2);
   else if (set==ODD_PTS)
      ia=(VOLUME/2);
   else if (set==NO_PTS)
      ia=VOLUME;

   for (ix=0;ix<VOLUME;ix++)
   {
      if ((ix>=ia)&&(ix<ib))
      {
         t=global_time(ix);

         if (((t==0)&&(bc!=3))||((t==(N0-1))&&(bc==0)))
         {
            for (n=0;n<6;n++)
            {
               ie|=(is_Xt_zero(xt[n])^0x1);
               xt[n]+=1;
            }
         }
         else
         {
            for (n=0;n<6;n++)
            {
               ie|=is_Xt_zero(xt[n]);
               xt[n]+=1;
            }
         }
      }
      else
      {
         for (n=0;n<6;n++)
         {
            ie|=(is_Xt_zero(xt[n])^0x1);
            xt[n]+=1;
         }
      }
   }

   error(ie!=0,1,"check_Xtbnd [check4.c]",
         "X tensor field vanishes on an incorrect set of points");
}


static void check_Xvbnd(void)
{
   int bc,ix,t,ifc,ie;
   su3_dble *xv;

   bc=bc_type();
   xv=xvector();
   ie=0;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (((t==0)&&(bc!=3))||((t==(N0-1))&&(bc==0)))
      {
         for (ifc=0;ifc<8;ifc++)
         {
            ie|=(is_Xv_zero(xv)^0x1);
            xv+=1;
         }
      }
      else if ((t==1)&&(bc!=3))
      {
         ie|=is_Xv_zero(xv);
         xv+=1;

         ie|=(is_Xv_zero(xv)^0x1);
         xv+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            ie|=is_Xv_zero(xv);
            xv+=1;
         }
      }
      else if (((t==(N0-2))&&(bc==0))||((t==(N0-1))&&(bc!=3)))
      {
         ie|=(is_Xv_zero(xv)^0x1);
         xv+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            ie|=is_Xv_zero(xv);
            xv+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            ie|=is_Xv_zero(xv);
            xv+=1;
         }
      }
   }

   error(ie!=0,1,"check_Xvbnd [check4.c]",
         "X vector field vanishes on an incorrect set of links");
}


static void check_bnd_frc(void)
{
   int bc,ix,t,ifc,ie;
   su3_alg_dble *frc;
   mdflds_t *mdfs;

   bc=bc_type();
   mdfs=mdflds();
   frc=(*mdfs).frc;
   ie=0;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t==0)&&(bc==0))
      {
         ie|=is_frc_zero(frc);
         frc+=1;

         ie|=(is_frc_zero(frc)^0x1);
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            ie|=is_frc_zero(frc);
            frc+=1;
         }
      }
      else if ((t==0)&&(bc==1))
      {
         ie|=is_frc_zero(frc);
         frc+=1;

         ie|=is_frc_zero(frc);
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            ie|=(is_frc_zero(frc)^0x1);
            frc+=1;
         }
      }
      else if ((t==(N0-1))&&(bc==0))
      {
         ie|=(is_frc_zero(frc)^0x1);
         frc+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            ie|=is_frc_zero(frc);
            frc+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            ie|=is_frc_zero(frc);
            frc+=1;
         }
      }
   }

   error(ie!=0,1,"check_bnd_frc [check4.c]",
         "Force field vanishes on an incorrect set of links");
}


static void rot_ud(double eps)
{
   int bc,ix,t,ifc;
   su3_dble *u;
   su3_alg_dble *mom;
   mdflds_t *mdfs;

   bc=bc_type();
   mdfs=mdflds();
   mom=(*mdfs).mom;
   u=udfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         expXsu3(eps,mom,u);
         mom+=1;
         u+=1;

         if (bc!=0)
            expXsu3(eps,mom,u);
         mom+=1;
         u+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
               expXsu3(eps,mom,u);
            mom+=1;
            u+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            expXsu3(eps,mom,u);
         mom+=1;
         u+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            mom+=1;
            u+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            mom+=1;
            u+=1;
         }
      }
   }

   set_flags(UPDATED_UD);
}


static double action(int k,spinor_dble **phi)
{
   int l;
   spinor_dble **wsd;
   double act;

   wsd=reserve_wsd(2);
   sw_term(NO_PTS);
   assign_sd2sd(VOLUME,phi[0],wsd[0]);

   for (l=0;l<k;l++)
   {
      Dw_dble(0.0,wsd[0],wsd[1]);
      mulg5_dble(VOLUME,wsd[1]);
      scale_dble(VOLUME,0.125,wsd[1]);
      assign_sd2sd(VOLUME,wsd[1],wsd[0]);
   }

   act=spinor_prod_re_dble(VOLUME,0,phi[0],wsd[0]);
   release_wsd();

   return act;
}


static double dSdt(int k,spinor_dble **phi)
{
   int l;
   spinor_dble **wsd;
   mdflds_t *mdfs;

   wsd=reserve_wsd(k);
   sw_term(NO_PTS);
   assign_sd2sd(VOLUME,phi[0],wsd[0]);

   for (l=1;l<k;l++)
   {
      Dw_dble(0.0,wsd[l-1],wsd[l]);
      mulg5_dble(VOLUME,wsd[l]);
      scale_dble(VOLUME,0.125,wsd[l]);
   }

   set_frc2zero();
   set_xt2zero();
   set_xv2zero();

   for (l=0;l<k;l++)
   {
      add_prod2xt(-0.0625,wsd[l],wsd[k-l-1]);
      add_prod2xv(-0.0625,wsd[l],wsd[k-l-1]);
   }

   check_Xtbnd(ALL_PTS);
   check_Xvbnd();

   sw_frc(1.0);
   hop_frc(1.0);
   check_bnd_frc();
   release_wsd();

   mdfs=mdflds();

   return scalar_prod_alg(4*VOLUME,0,(*mdfs).mom,(*mdfs).frc);
}


static double action_det(ptset_t set)
{
   int bc,ie,io;
   int vol,ofs,ix,im,t,n;
   double c,p;
   complex_dble z;
   pauli_dble *m;
   sw_parms_t swp;

   if (set==NO_PTS)
      return 0.0;

   bc=bc_type();
   swp=sw_parms();

   if ((4.0+swp.m0)>1.0)
      c=pow(4.0+swp.m0,-6.0);
   else
      c=1.0;

   for (n=0;n<MAX_LEVELS;n++)
   {
      cnt[n]=0;
      smx[n]=0.0;
   }

   if (query_flags(SWD_UP2DATE)!=1)
      sw_term(NO_PTS);
   else
   {
      ie=query_flags(SWD_E_INVERTED);
      io=query_flags(SWD_O_INVERTED);

      if (((ie==1)&&((set==ALL_PTS)||(set==EVEN_PTS)))||
          ((io==1)&&((set==ALL_PTS)||(set==ODD_PTS))))
         sw_term(NO_PTS);
   }

   if (set==ODD_PTS)
      ofs=(VOLUME/2);
   else
      ofs=0;

   if (set==EVEN_PTS)
      vol=(VOLUME/2);
   else
      vol=VOLUME;

   ix=ofs;
   m=swdfld()+2*ofs;

   while (ix<vol)
   {
      im=ix+BLK_LENGTH;
      if (im>vol)
         im=vol;
      p=1.0;

      for (;ix<im;ix++)
      {
         t=global_time(ix);

         if (((t>0)||(bc==3))&&((t<(N0-1))||(bc!=0)))
         {
            z=det_pauli_dble(0.0,m);
            p*=(c*z.re);
            z=det_pauli_dble(0.0,m+1);
            p*=(c*z.re);
         }

         m+=2;
      }

      cnt[0]+=1;
      smx[0]-=log(fabs(p));

      for (n=1;(cnt[n-1]>=BLK_LENGTH)&&(n<MAX_LEVELS);n++)
      {
         cnt[n]+=1;
         smx[n]+=smx[n-1];

         cnt[n-1]=0;
         smx[n-1]=0.0;
      }
   }

   for (n=1;n<MAX_LEVELS;n++)
      smx[0]+=smx[n];

   return 2.0*smx[0];
}


static double dSdt_det(ptset_t set)
{
   int ifail;
   mdflds_t *mdfs;

   set_xt2zero();
   ifail=add_det2xt(2.0,set);
   error_root(ifail!=0,1,"dSdt_det [check4.c]",
              "Inversion of the SW term was not safe");
   check_Xtbnd(set);

   set_frc2zero();
   sw_frc(1.0);

   if (set==ALL_PTS)
      check_bnd_frc();

   mdfs=mdflds();

   return scalar_prod_alg(4*VOLUME,0,(*mdfs).mom,(*mdfs).frc);
}


int main(int argc,char *argv[])
{
   int my_rank,bc,k;
   double chi[2],chi_prime[2],theta[3];
   double eps,act0,act1,dsdt;
   double dev_frc,sig_loss,s[2],r[2];
   spinor_dble **phi;
   ptset_t set;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);

      printf("\n");
      printf("Check of sw_frc() and hop_frc()\n");
      printf("-------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>]");
   }

   set_lat_parms(5.5,1.0,0,NULL,1.782);
   print_lat_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   chi[0]=0.123;
   chi[1]=-0.534;
   chi_prime[0]=0.912;
   chi_prime[1]=0.078;
   theta[0]=0.38;
   theta[1]=-1.25;
   theta[2]=0.54;
   set_bc_parms(bc,1.0,1.0,0.953,1.203,chi,chi_prime,theta);
   print_bc_parms(1);

   start_ranlux(0,1245);
   geometry();

   set_sw_parms(-0.0123);
   alloc_wsd(6);
   phi=reserve_wsd(1);

   for (k=1;k<=4;k++)
   {
      random_ud();
      set_ud_phase();
      random_mom();
      random_sd(VOLUME,phi[0],1.0);
      bnd_sd2zero(ALL_PTS,phi[0]);
      dsdt=dSdt(k,phi);

      eps=5.0e-5;
      rot_ud(eps);
      act0=2.0*action(k,phi)/3.0;
      rot_ud(-eps);

      rot_ud(-eps);
      act1=2.0*action(k,phi)/3.0;
      rot_ud(eps);

      rot_ud(2.0*eps);
      act0-=action(k,phi)/12.0;
      rot_ud(-2.0*eps);

      rot_ud(-2.0*eps);
      act1-=action(k,phi)/12.0;
      rot_ud(2.0*eps);

      s[0]=dsdt-(act0-act1)/eps;
      s[1]=dsdt;

      if (NPROC>1)
      {
         MPI_Reduce(s,r,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
         MPI_Bcast(r,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }
      else
      {
         r[0]=s[0];
         r[1]=s[1];
      }

      dev_frc=fabs(r[0]/r[1]);
      sig_loss=-log10(fabs(1.0-act0/act1));

      if (my_rank==0)
      {
         printf("Calculation of the force for S=(phi,Q^%d*phi):\n",k);
         printf("Relative deviation of dS/dt = %.2e ",dev_frc);
         printf("[significance loss = %d digits]\n\n",(int)(sig_loss));
      }
   }

   if (my_rank==0)
      printf("Calculation of the force for S=-2*Tr{ln(SW term)}:\n");

   for (k=0;k<4;k++)
   {
      if (k==0)
         set=NO_PTS;
      else if (k==1)
         set=EVEN_PTS;
      else if (k==2)
         set=ODD_PTS;
      else
         set=ALL_PTS;

      random_ud();
      set_ud_phase();
      random_mom();
      dsdt=dSdt_det(set);

      eps=5.0e-4;
      rot_ud(eps);
      act0=2.0*action_det(set)/3.0;
      rot_ud(-eps);

      rot_ud(-eps);
      act1=2.0*action_det(set)/3.0;
      rot_ud(eps);

      rot_ud(2.0*eps);
      act0-=action_det(set)/12.0;
      rot_ud(-2.0*eps);

      rot_ud(-2.0*eps);
      act1-=action_det(set)/12.0;
      rot_ud(2.0*eps);

      s[0]=dsdt-(act0-act1)/eps;
      s[1]=dsdt;

      if (NPROC>1)
      {
         MPI_Reduce(s,r,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
         MPI_Bcast(r,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }
      else
      {
         r[0]=s[0];
         r[1]=s[1];
      }

      if (k>0)
      {
         dev_frc=fabs(r[0]/r[1]);
         sig_loss=-log10(fabs(1.0-act0/act1));
      }
      else
         dev_frc=fabs(r[0]);

      if (my_rank==0)
      {
         if (k==0)
            printf("set=NO_PTS:   ");
         else if (k==1)
            printf("set=EVEN_PTS: ");
         else if (k==2)
            printf("set=ODD_PTS:  ");
         else
            printf("set=ALL_PTS:  ");

         if (k>0)
         {
            printf("relative deviation of dS/dt = %.2e ",dev_frc);
            printf("[significance loss = %d digits]\n",(int)(sig_loss));
         }
         else
            printf("absolute deviation of dS/dt = %.2e\n",dev_frc);
      }
   }

   if (my_rank==0)
   {
      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
