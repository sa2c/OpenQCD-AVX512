
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005, 2007, 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the program that translates the double-precision gauge field.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "global.h"

static int my_rank,ipsnd,iprcv,*isnd;
static su3_dble *uold,*unew,*ubuf;


static void alloc_bufs(void)
{
   isnd=amalloc(NPROC*sizeof(*isnd),3);
   uold=amalloc(12*VOLUME*sizeof(*uold),ALIGN);
   error((isnd==NULL)||(uold==NULL),1,"alloc_bufs [check3.c]",
         "Unable to allocate auxiliary arrays");

   unew=uold+4*VOLUME;
   ubuf=unew+4*VOLUME;
}


static int range(int *dist,int *s,int *ra,int *rb)
{
   int io,l[4],nl[4];
   int mu,a,b;

   io=1;

   l[0]=L0;
   l[1]=L1;
   l[2]=L2;
   l[3]=L3;

   nl[0]=L0*NPROC0;
   nl[1]=L1*NPROC1;
   nl[2]=L2*NPROC2;
   nl[3]=L3*NPROC3;

   for (mu=0;mu<4;mu++)
   {
      a=dist[mu]+s[mu];
      b=a+l[mu];

      a=safe_mod(a,nl[mu]);
      b=safe_mod(b,nl[mu]);

      if (a==b)
      {
         ra[mu]=0;
         rb[mu]=l[mu];
      }
      else if (a<l[mu])
      {
         ra[mu]=a;
         rb[mu]=l[mu];
      }
      else if (b<l[mu])
      {
         ra[mu]=0;
         rb[mu]=b;
      }
      else
      {
         io=0;
         ra[mu]=0;
         rb[mu]=0;
      }
   }

   return io;
}


static void snd_pairs(int *dist)
{
   int c0,c1,c2,c3,ic;
   int bo0[4],bo1[4],ip0,ip1,dmy;

   for (ic=0;ic<NPROC;ic++)
      isnd[ic]=0;

   for (c0=0;c0<(NPROC0*L0);c0+=L0)
   {
      for (c1=0;c1<(NPROC1*L1);c1+=L1)
      {
         for (c2=0;c2<(NPROC2*L2);c2+=L2)
         {
            for (c3=0;c3<(NPROC3*L3);c3+=L3)
            {
               bo0[0]=c0;
               bo0[1]=c1;
               bo0[2]=c2;
               bo0[3]=c3;

               bo1[0]=c0+dist[0];
               bo1[1]=c1+dist[1];
               bo1[2]=c2+dist[2];
               bo1[3]=c3+dist[3];

               ipt_global(bo0,&ip0,&dmy);
               ipt_global(bo1,&ip1,&dmy);

               if (my_rank==ip0)
                  iprcv=ip1;
               if (my_rank==ip1)
                  ipsnd=ip0;

               if ((isnd[ip1]==0)&&(ip0!=ip1))
               {
                  isnd[ip1]=1;
                  isnd[ip0]=-1;
               }
            }
         }
      }
   }
}


static void snd_field(int *dist)
{
   int nbuf,tag,n;
   MPI_Status stat;

   snd_pairs(dist);
   nbuf=4*18*VOLUME;
   tag=mpi_tag();

   if (isnd[my_rank]==1)
   {
      MPI_Send(uold,nbuf,MPI_DOUBLE,ipsnd,tag,MPI_COMM_WORLD);
      MPI_Recv(ubuf,nbuf,MPI_DOUBLE,iprcv,tag,MPI_COMM_WORLD,&stat);
   }
   else if (isnd[my_rank]==-1)
   {
      MPI_Recv(ubuf,nbuf,MPI_DOUBLE,iprcv,tag,MPI_COMM_WORLD,&stat);
      MPI_Send(uold,nbuf,MPI_DOUBLE,ipsnd,tag,MPI_COMM_WORLD);
   }
   else
   {
      for (n=0;n<(4*VOLUME);n++)
         ubuf[n]=uold[n];
   }
}


static void save_field(su3_dble *u)
{
   int ix,iy,ip[4];
   su3_dble *ub;

   copy_bnd_ud();
   ub=udfld();

   for (ix=0;ix<VOLUME;ix++)
   {
      iy=ipt[ix];

      plaq_uidx(0,iy,ip);
      (*u)=ub[ip[0]];
      u+=1;
      (*u)=ub[ip[2]];
      u+=1;

      plaq_uidx(3,iy,ip);
      (*u)=ub[ip[0]];
      u+=1;
      (*u)=ub[ip[2]];
      u+=1;
   }
}


static int cmp_su3_dble(su3_dble *u,su3_dble *v)
{
   int k;
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

   for (k=0;k<18;k++)
   {
      if (r[k]!=0.0)
         return 1;
   }

   return 0;
}


static int cmp_field(int *s)
{
   int c0,c1,c2,c3,dist[4];
   int io,ra[4],rb[4];
   int x0,x1,x2,x3,ix;
   int y0,y1,y2,y3,iy;
   int mu,itest;

   itest=0;

   for (c0=0;c0<(NPROC0*L0);c0+=L0)
   {
      for (c1=0;c1<(NPROC1*L1);c1+=L1)
      {
         for (c2=0;c2<(NPROC2*L2);c2+=L2)
         {
            for (c3=0;c3<(NPROC3*L3);c3+=L3)
            {
               dist[0]=c0;
               dist[1]=c1;
               dist[2]=c2;
               dist[3]=c3;

               io=range(dist,s,ra,rb);

               if (io!=0)
               {
                  snd_field(dist);

                  for (x0=ra[0];x0<rb[0];x0++)
                  {
                     for (x1=ra[1];x1<rb[1];x1++)
                     {
                        for (x2=ra[2];x2<rb[2];x2++)
                        {
                           for (x3=ra[3];x3<rb[3];x3++)
                           {
                              y0=safe_mod(x0-s[0],L0);
                              y1=safe_mod(x1-s[1],L1);
                              y2=safe_mod(x2-s[2],L2);
                              y3=safe_mod(x3-s[3],L3);

                              ix=x3+L3*x2+L2*L3*x1+L1*L2*L3*x0;
                              iy=y3+L3*y2+L2*L3*y1+L1*L2*L3*y0;

                              for (mu=0;mu<4;mu++)
                                 itest|=cmp_su3_dble(unew+4*ix+mu,ubuf+4*iy+mu);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (NPROC>1)
   {
      io=itest;
      MPI_Reduce(&io,&itest,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&itest,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   return itest;
}


static void random_vec(int *svec)
{
   int mu,bs[4];
   double r[4];

   bs[0]=NPROC0*L0;
   bs[1]=NPROC1*L1;
   bs[2]=NPROC2*L2;
   bs[3]=NPROC3*L3;

   ranlxd(r,4);

   for (mu=0;mu<4;mu++)
   {
      svec[mu]=(int)((double)(bs[mu])*r[mu]);
      if (svec[mu]>(bs[mu]/2))
         svec[mu]-=bs[mu];
   }

   MPI_Bcast(svec,4,MPI_INT,0,MPI_COMM_WORLD);
}


int main(int argc,char *argv[])
{
   int bc,ie;
   int ifc,mu,s[4],n,itest;
   double phi[2],phi_prime[2],theta[3];
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);

      printf("\n");
      printf("Translation of the double-precision gauge field\n");
      printf("-----------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(0);

   geometry();
   alloc_bufs();

   if (my_rank==0)
      printf("Elementary shift vectors:\n\n");

   for (ifc=0;ifc<8;ifc++)
   {
      if ((ifc>1)||(bc==3))
      {
         random_ud();
         save_field(uold);

         s[0]=0;
         s[1]=0;
         s[2]=0;
         s[3]=0;
         mu=ifc/2;

         if ((ifc&0x1)==0)
            s[mu]=1;
         else
            s[mu]=-1;

         shift_ud(s);
         save_field(unew);
         itest=cmp_field(s);

         ie=check_bc(0.0);
         error_root(ie==0,1,"main [check3.c]","Boundary conditions changed");

         if (my_rank==0)
         {
            printf("Shift vector (% 3d,% 3d,% 3d,% 3d): ",
                   s[0],s[1],s[2],s[3]);

            if (itest==0)
               printf("ok\n");
            else
               printf("failed\n");
         }
      }
   }

   if (my_rank==0)
   {
      printf("\n");
      printf("Random shift vectors:\n\n");
   }

   for (n=0;n<8;n++)
   {
      random_ud();
      save_field(uold);

      random_vec(s);
      if (bc!=3)
         s[0]=0;
      shift_ud(s);
      save_field(unew);
      itest=cmp_field(s);

      ie=check_bc(0.0);
      error_root(ie==0,1,"main [check3.c]","Boundary conditions changed");

      if (my_rank==0)
      {
         printf("Shift vector (% 3d,% 3d,% 3d,% 3d): ",
                s[0],s[1],s[2],s[3]);

         if (itest==0)
            printf("ok\n");
         else
            printf("failed\n");
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
