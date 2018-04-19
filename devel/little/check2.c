
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the program b2b_flds().
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
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "dfl.h"
#include "little.h"
#include "global.h"

static int bs[4],Ns,bc;
static int l[4],np[4];
static const su3_dble ud0={{0.0}};


static void set_ud(void)
{
   su3_dble unity,*ud,*um;

   unity=ud0;
   unity.c11.re=1.0;
   unity.c22.re=1.0;
   unity.c33.re=1.0;
   ud=udfld();
   um=ud+4*VOLUME;

   for (;ud<um;ud++)
      (*ud)=unity;
}


static void set_sd(spinor_dble *sd)
{
   int x0,x1,x2,x3,ix;
   int y0,y1,y2,y3;

   set_sd2zero(VOLUME,sd);

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               y0=x0+cpr[0]*L0;
               y1=x1+cpr[1]*L1;
               y2=x2+cpr[2]*L2;
               y3=x3+cpr[3]*L3;

               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

               sd[ix].c1.c2.im=(double)(y0);
               sd[ix].c2.c2.im=(double)(y1);
               sd[ix].c3.c2.im=(double)(y2);
               sd[ix].c4.c2.im=(double)(y3);
            }
         }
      }
   }
}


static void chk_sde0(int mu,int vol,int ibn,int *bo,spinor_dble *sd)
{
   int a[4],b[4],y[4];
   int ix,nu,ie;

   for (nu=0;nu<4;nu++)
   {
      a[nu]=cpr[nu]*l[nu]+bo[nu];
      b[nu]=a[nu]+bs[nu];
   }

   a[mu]=cpr[mu]*l[mu]+bo[mu]+bs[mu]-1;
   if (ibn)
      a[mu]=safe_mod(a[mu]-l[mu],np[mu]*l[mu]);
   b[mu]=a[mu]+1;
   ie=0;

   for (ix=0;ix<vol;ix++)
   {
      y[0]=sd[ix].c1.c2.im;
      y[1]=sd[ix].c2.c2.im;
      y[2]=sd[ix].c3.c2.im;
      y[3]=sd[ix].c4.c2.im;

      if (((y[0]+y[1]+y[2]+y[3])&0x1)!=0)
         ie=1;

      for (nu=0;nu<4;nu++)
      {
         if ((y[nu]<a[nu])||(y[nu]>=b[nu]))
            ie=2;
      }
   }

   error(ie!=0,1,"chk_sde0 [check2.c]","Incorrect field components");
}


static void chk_sde1(int mu,int vol,int ibn,int *bo,spinor_dble *sd)
{
   int a[4],b[4],y[4];
   int ix,nu,ie;

   for (nu=0;nu<4;nu++)
   {
      a[nu]=cpr[nu]*l[nu]+bo[nu];
      b[nu]=a[nu]+bs[nu];
   }

   a[mu]=cpr[mu]*l[mu]+bo[mu];
   if (ibn)
      a[mu]=safe_mod(a[mu]+l[mu],np[mu]*l[mu]);
   b[mu]=a[mu]+1;
   ie=0;

   for (ix=0;ix<vol;ix++)
   {
      y[0]=sd[ix].c1.c2.im;
      y[1]=sd[ix].c2.c2.im;
      y[2]=sd[ix].c3.c2.im;
      y[3]=sd[ix].c4.c2.im;

      if (((y[0]+y[1]+y[2]+y[3])&0x1)!=0)
         ie=1;

      for (nu=0;nu<4;nu++)
      {
         if ((y[nu]<a[nu])||(y[nu]>=b[nu]))
            ie=2;
      }
   }

   error(ie!=0,1,"chk_sde1 [check2.c]","Incorrect field components");
}


static void chk_sdo0(int mu,int vol,int *bo,spinor_dble *sd)
{
   int a[4],b[4],y[4];
   int ix,nu,ie;

   for (nu=0;nu<4;nu++)
   {
      a[nu]=cpr[nu]*l[nu]+bo[nu];
      b[nu]=a[nu]+bs[nu];
   }

   a[mu]=cpr[mu]*l[mu]+bo[mu]+bs[mu]-1;
   b[mu]=a[mu]+1;
   ie=0;

   for (ix=0;ix<vol;ix++)
   {
      y[0]=sd[ix].c1.c2.im;
      y[1]=sd[ix].c2.c2.im;
      y[2]=sd[ix].c3.c2.im;
      y[3]=sd[ix].c4.c2.im;

      if (((y[0]+y[1]+y[2]+y[3])&0x1)!=1)
         ie=1;

      for (nu=0;nu<4;nu++)
      {
         if ((y[nu]<a[nu])||(y[nu]>=b[nu]))
            ie=2;
      }
   }

   error(ie!=0,1,"chk_sdo0 [check2.c]","Incorrect field components %d",ie);
}


static void chk_sdo1(int mu,int vol,int *bo,spinor_dble *sd)
{
   int a[4],b[4],y[4];
   int ix,nu,ie;

   for (nu=0;nu<4;nu++)
   {
      a[nu]=cpr[nu]*l[nu]+bo[nu];
      b[nu]=a[nu]+bs[nu];
   }

   a[mu]=cpr[mu]*l[mu]+bo[mu];
   b[mu]=a[mu]+1;
   ie=0;

   for (ix=0;ix<vol;ix++)
   {
      y[0]=sd[ix].c1.c2.im;
      y[1]=sd[ix].c2.c2.im;
      y[2]=sd[ix].c3.c2.im;
      y[3]=sd[ix].c4.c2.im;

      if (((y[0]+y[1]+y[2]+y[3])&0x1)!=1)
         ie=1;

      for (nu=0;nu<4;nu++)
      {
         if ((y[nu]<a[nu])||(y[nu]>=b[nu]))
            ie=2;
      }
   }

   error(ie!=0,1,"chk_sdo1 [check2.c]","Incorrect field components");
}


static void cmp_sde0_sdo1(int mu,int vol,int *bo,spinor_dble *sde,
                          spinor_dble *sdo)
{
   int ye[4],yo[4];
   int ix,nu,ie;

   ie=0;

   for (ix=0;ix<vol;ix++)
   {
      ye[0]=sde[ix].c1.c2.im;
      ye[1]=sde[ix].c2.c2.im;
      ye[2]=sde[ix].c3.c2.im;
      ye[3]=sde[ix].c4.c2.im;

      yo[0]=sdo[ix].c1.c2.im;
      yo[1]=sdo[ix].c2.c2.im;
      yo[2]=sdo[ix].c3.c2.im;
      yo[3]=sdo[ix].c4.c2.im;

      for (nu=0;nu<4;nu++)
      {
         if ((nu!=mu)&&(ye[nu]!=yo[nu]))
            ie=1;

         if ((nu==mu)&&(ye[nu]!=safe_mod(yo[nu]-1,np[mu]*l[mu])))
            ie=2;
      }
   }

   error(ie!=0,1,"cmp_sde0_sdo1 [check2.c]","Incorrect field components");
}


static void cmp_sde1_sdo0(int mu,int vol,int *bo,spinor_dble *sde,
                          spinor_dble *sdo)
{
   int ye[4],yo[4];
   int ix,nu,ie;

   ie=0;

   for (ix=0;ix<vol;ix++)
   {
      ye[0]=sde[ix].c1.c2.im;
      ye[1]=sde[ix].c2.c2.im;
      ye[2]=sde[ix].c3.c2.im;
      ye[3]=sde[ix].c4.c2.im;

      yo[0]=sdo[ix].c1.c2.im;
      yo[1]=sdo[ix].c2.c2.im;
      yo[2]=sdo[ix].c3.c2.im;
      yo[3]=sdo[ix].c4.c2.im;

      for (nu=0;nu<4;nu++)
      {
         if ((nu!=mu)&&(ye[nu]!=yo[nu]))
            ie=1;

         if ((nu==mu)&&(ye[nu]!=safe_mod(yo[nu]+1,np[mu]*l[mu])))
            ie=2;
      }
   }

   error(ie!=0,1,"cmp_sde0_sdo1 [check2.c]","Incorrect field components");
}


static void chk_b2b(int n,int mu,b2b_flds_t *b2b)
{
   int nb,isw;
   int *m,vol,ibn,ie;
   int *bo0,*bo1,nu,k;
   block_t *b,*b0,*b1;

   l[0]=L0;
   l[1]=L1;
   l[2]=L2;
   l[3]=L3;

   np[0]=NPROC0;
   np[1]=NPROC1;
   np[2]=NPROC2;
   np[3]=NPROC3;

   b=blk_list(DFL_BLOCKS,&nb,&isw);

   m=(*b2b).n;
   b0=b+m[0];
   b1=b+m[1];
   vol=(*b0).bb[2*mu+1].vol/2;
   ibn=(*b0).bb[2*mu+1].ibn;

   error((n!=m[0])||(vol!=(*b2b).vol)||(ibn!=(*b2b).ibn),1,
         "chk_b2b [check2.c]","Incorrect b2b.n, b2b.vol or b2b.ibn");

   bo0=(*b0).bo;
   bo1=(*b1).bo;
   ie=0;

   for (nu=0;nu<4;nu++)
   {
      if ((nu!=mu)&&(bo0[nu]!=bo1[nu]))
         ie=1;
   }

   if (bo1[mu]!=((bo0[mu]+bs[mu])%l[mu]))
      ie=2;

   error(ie!=0,1,"chk_b2b [check2.c]","Blocks are not neighbours");

   for (k=0;k<Ns;k++)
   {
      chk_sde0(mu,vol,ibn,bo0,(*b2b).sde[0][k]);
      chk_sde1(mu,vol,ibn,bo1,(*b2b).sde[1][k]);
      chk_sdo0(mu,vol,bo0,(*b2b).sdo[0][k]);
      chk_sdo1(mu,vol,bo1,(*b2b).sdo[1][k]);

      cmp_sde0_sdo1(mu,vol,bo1,(*b2b).sde[0][k],(*b2b).sdo[1][k]);
      cmp_sde1_sdo0(mu,vol,bo0,(*b2b).sde[1][k],(*b2b).sdo[0][k]);
   }
}


int main(int argc,char *argv[])
{
   int my_rank,nb,isw;
   int n,mu,k;
   double phi[2],phi_prime[2],theta[3];
   spinor_dble **wsd;
   b2b_flds_t *b2b;
   FILE *fin=NULL,*flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      fin=freopen("check2.in","r",stdin);

      printf("\n");
      printf("Check of the programs b2b_flds()\n");
      printf("--------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n\n",bs[0],bs[1],bs[2],bs[3]);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }

   set_lat_parms(5.5,1.0,0,NULL,1.978);
   print_lat_parms();

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,1.0,1.0,0.9012,1.2034,phi,phi_prime,theta);
   print_bc_parms(2);

   start_ranlux(0,123456);
   geometry();
   Ns=2;
   set_dfl_parms(bs,Ns);
   alloc_bgr(DFL_BLOCKS);
   blk_list(DFL_BLOCKS,&nb,&isw);

   alloc_wsd(Ns);
   wsd=reserve_wsd(Ns);
   set_ud();

   for (k=0;k<Ns;k++)
      set_sd(wsd[k]);

   for (n=0;n<nb;n++)
   {
      for (k=0;k<Ns;k++)
         assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,wsd[k],k+1);
   }

   for (n=0;n<nb;n++)
   {
      for (mu=0;mu<4;mu++)
      {
         b2b=b2b_flds(n,mu);
         chk_b2b(n,mu,b2b);
      }
   }

   if (my_rank==0)
   {
      printf("No errors detected\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
