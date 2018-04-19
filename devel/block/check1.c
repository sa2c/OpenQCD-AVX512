
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency checks on the geometry arrays in the known block grids.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "block.h"
#include "global.h"

static int ix_test[VOLUME+BNDRY];


static void test1(blk_grid_t grid,int *bs)
{
   int ix,iy,itest;
   int nb,isw,vol,*imb;
   block_t *b,*bm;

   for (ix=0;ix<VOLUME;ix++)
      ix_test[ix]=0;

   itest=0;
   b=blk_list(grid,&nb,&isw);
   bm=b+nb;

   error(((*b).bo[0]!=0)||((*b).bo[1]!=0)||((*b).bo[2]!=0)||((*b).bo[3]!=0),1,
         "test1 [check1.c]","Offset of the first block is incorrect");

   error((bs[0]!=((*b).bs[0]))||(bs[1]!=((*b).bs[1]))||
         (bs[2]!=((*b).bs[2]))||(bs[3]!=((*b).bs[3])),1,
         "test1 [check1.c]","Block sizes are not correctly assigned");

   for (;b<bm;b++)
   {
      vol=(*b).vol;
      imb=(*b).imb;

      if (vol!=(bs[0]*bs[1]*bs[2]*bs[3]))
         itest=1;

      for (ix=0;ix<vol;ix++)
      {
         iy=imb[ix];

         if ((iy<0)||(iy>=VOLUME))
            itest=2;
         else
         {
            ix_test[iy]+=1;
            if (ix_test[iy]>1)
               itest=3;
         }
      }
   }

   for (ix=0;ix<VOLUME;ix++)
      if (ix_test[ix]!=1)
         itest=4;

   error(itest==1,1,"test1 [check1.c]","b.vol is incorrect");
   error(itest==2,1,"test1 [check1.c]","b.imb is out of range");
   error(itest==3,1,"test1 [check1.c]","b.imb is not one-to-one");
   error(itest==4,1,"test1 [check1.c]","The blocks do not cover the lattice");
}


static void test2(blk_grid_t grid)
{
   int ix,iy,mu,itest;
   int x0,x1,x2,x3,x[4];
   int y0,y1,y2,y3,is;
   int nb,isw,vol,*bo,*bs,*imb;
   block_t *b,*bm;

   itest=0;
   b=blk_list(grid,&nb,&isw);
   bm=b+nb;

   for (;b<bm;b++)
   {
      vol=(*b).vol;
      bo=(*b).bo;
      bs=(*b).bs;
      imb=(*b).imb;

      for (x0=0;x0<bs[0];x0++)
      {
         for (x1=0;x1<bs[1];x1++)
         {
            for (x2=0;x2<bs[2];x2++)
            {
               for (x3=0;x3<bs[3];x3++)
               {
                  x[0]=x0;
                  x[1]=x1;
                  x[2]=x2;
                  x[3]=x3;

                  y0=bo[0]+x0;
                  y1=bo[1]+x1;
                  y2=bo[2]+x2;
                  y3=bo[3]+x3;

                  iy=ipt[y3+L3*y2+L2*L3*y1+L1*L2*L3*y0];
                  ix=ipt_blk(b,x);

                  if ((ix<0)||(ix>vol))
                     itest=1;
                  else
                  {
                     if (iy!=imb[ix])
                        itest=2;

                     is=(x0+x1+x2+x3+bo[0]+bo[1]+bo[2]+bo[3])%2;

                     if (((is==0)&&(ix>=(vol/2)))||((is!=0)&&(ix<(vol/2))))
                        itest=3;

                     for (mu=0;mu<4;mu++)
                     {
                        if ((x[mu]+1)<bs[mu])
                        {
                           if (imb[(*b).iup[ix][mu]]!=iup[iy][mu])
                              itest=4;
                        }
                        else
                        {
                           if ((*b).iup[ix][mu]!=vol)
                              itest=5;
                        }

                        if (x[mu]>0)
                        {
                           if (imb[(*b).idn[ix][mu]]!=idn[iy][mu])
                              itest=6;
                        }
                        else
                        {
                           if ((*b).idn[ix][mu]!=vol)
                              itest=7;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   error(itest==1,1,"test2 [check1.c]",
         "b.ipt is out of range");
   error(itest==2,1,"test2 [check1.c]",
         "The blocks are not properly embedded");
   error(itest==3,1,"test2 [check1.c]",
         "b.ipt does not respect the even-odd ordering");
   error(itest==4,1,"test2 [check1.c]",
         "b.iup is incorrect");
   error(itest==5,1,"test2 [check1.c]",
         "b.iup is incorrect at the block boundary");
   error(itest==6,1,"test2 [check1.c]",
         "b.idn is incorrect");
   error(itest==7,1,"test2 [check1.c]",
         "b.idn is incorrect at the block boundary");
}


static void test3(blk_grid_t grid)
{
   int bc,ix,iy,ie,itest;
   int nbp,nall,x[4];
   int nb,isw,vol,*bs,*imb;
   block_t *b,*bm;

   bc=bc_type();
   itest=0;
   nall=0;
   b=blk_list(grid,&nb,&isw);
   bm=b+nb;

   for (;b<bm;b++)
   {
      vol=(*b).vol;
      bs=(*b).bs;
      imb=(*b).imb;

      nbp=0;
      x[0]=0;
      x[1]=0;
      x[2]=0;
      x[3]=0;

      ix=ipt_blk(b,x);
      ix=imb[ix];
      if ((global_time(ix)==0)&&(bc!=3))
         nbp+=(bs[1]*bs[2]*bs[3]);

      x[0]=bs[0]-1;
      ix=ipt_blk(b,x);
      ix=imb[ix];
      if ((global_time(ix)==(NPROC0*L0-1))&&(bc==0))
         nbp+=(bs[1]*bs[2]*bs[3]);

      if ((*b).nbp!=nbp)
         itest=1;

      nall+=nbp;
      nbp=(*b).nbp;

      for (iy=0;iy<nbp;iy++)
      {
         ix=(*b).ibp[iy];

         if ((ix<0)||(ix>=vol))
            itest=2;

         if (iy>0)
         {
            if (ix<=(*b).ibp[iy-1])
               itest=3;
         }

         ix=imb[ix];
         ie=((global_time(ix)==0)&&(bc!=3));
         ie|=((global_time(ix)==(NPROC0*L0-1))&&(bc==0));

         if (ie==0)
            itest=4;
      }
   }

   if ((cpr[0]==0)&&(bc!=3))
      nall-=(L1*L2*L3);
   if ((cpr[0]==(NPROC0-1))&&(bc==0))
      nall-=(L1*L2*L3);

   error(itest==1,1,"test3 [check1.c]",
         "b.nbp is incorrect");
   error(itest==2,1,"test3 [check1.c]",
         "b.ibp is out of range");
   error(itest==3,1,"test3 [check1.c]",
         "b.ibp is not properly ordered");
   error(itest==4,1,"test3 [check1.c]",
         "The points b.ibp are not all on the boundary of the lattice");
   error(nall!=0,1,"test3 [check1.c]",
         "Incorrect total count of points at time 0 and NPROC0*L0-1");
}


static void test4(blk_grid_t grid)
{
   int ix,iy,ifc,mu,ib,itest;
   int nb,isw,vol,*bs,*imb;
   block_t *b,*bm;
   bndry_t *bb;

   for (ix=0;ix<(VOLUME+BNDRY);ix++)
      ix_test[ix]=0;

   itest=0;
   b=blk_list(grid,&nb,&isw);
   bm=b+nb;

   for (;b<bm;b++)
   {
      vol=(*b).vol;
      bs=(*b).bs;
      bb=(*b).bb;

      for (ifc=0;ifc<8;ifc++)
      {
         mu=ifc/2;

         if (((*bb).ifc!=ifc)||((*bb).vol!=(vol/bs[mu])))
            itest=1;

         for (ix=0;ix<(*bb).vol;ix++)
         {
            iy=(*bb).ipp[ix];

            if ((iy<0)||(iy>=vol))
               itest=2;

            iy=(*bb).imb[ix];

            if ((iy<0)||(iy>=(VOLUME+BNDRY)))
               itest=3;
            else
               ix_test[iy]+=1;
         }

         bb+=1;
      }
   }

   b=blk_list(grid,&nb,&isw);

   for (;b<bm;b++)
   {
      vol=(*b).vol;
      imb=(*b).imb;

      for (ix=0;ix<vol;ix++)
      {
         iy=imb[ix];
         ib=0;

         for (mu=0;mu<4;mu++)
         {
            if (((*b).iup[ix][mu]==vol)&&(iup[iy][mu]<VOLUME))
               ib+=1;
            if (((*b).idn[ix][mu]==vol)&&(idn[iy][mu]<VOLUME))
               ib+=1;
         }

         ix_test[iy]+=(1-ib);
      }
   }

   for (ix=0;ix<(VOLUME+BNDRY);ix++)
      if (ix_test[ix]!=1)
         itest=4;

   error(itest==1,1,"test4 [check1.c]",
         "bb.ifc and bb.vol are not correctly assigned");
   error(itest==2,1,"test4 [check1.c]",
         "bb.ipp is out of range");
   error(itest==3,1,"test4 [check1.c]",
         "bb.imb is out of range");
   error(itest==4,1,"test4 [check1.c]",
         "bb.imb is not one-to-one or not surjective");
}


static void test5(blk_grid_t grid)
{
   int ix,iy,iz,ifc,mu,is,itest;
   int nb,isw,vol,*bs,*imb;
   block_t *b,*bm;
   bndry_t *bb;

   for (ix=0;ix<(VOLUME+BNDRY);ix++)
      ix_test[ix]=0;

   itest=0;
   b=blk_list(grid,&nb,&isw);
   bm=b+nb;

   for (;b<bm;b++)
   {
      vol=(*b).vol;
      bs=(*b).bs;
      imb=(*b).imb;
      bb=(*b).bb;

      for (ifc=0;ifc<8;ifc++)
      {
         for (ix=0;ix<((*bb).vol/2);ix++)
         {
            iy=(*bb).ipp[ix];
            iz=(*bb).ipp[ix+((*bb).vol/2)];

            if ((iy<(vol/2))||(iz>=(vol/2)))
               itest=1;

            iy=(*bb).imb[ix];
            iz=(*bb).imb[ix+(*bb).vol/2];

            if ((iy>=(VOLUME+(BNDRY/2)))||((iy>=(VOLUME/2))&&(iy<VOLUME)))
               itest=2;

            if ((iz<(VOLUME/2))||((iz>=VOLUME)&&(iz<(VOLUME+(BNDRY/2)))))
               itest=2;
         }

         mu=ifc/2;

         for (ix=0;ix<(*bb).vol;ix++)
         {
            iy=(*bb).ipp[ix];

            if ((((ifc%2)==0)&&((*b).idn[iy][mu]!=vol))||
                (((ifc%2)==1)&&((*b).iup[iy][mu]!=vol)))
               itest=3;

            iz=(*bb).map[ix];

            if ((ifc%2)==0)
            {
               for (is=1;is<bs[mu];is++)
                  iz=(*b).idn[iz][mu];
            }
            else
            {
               for (is=1;is<bs[mu];is++)
                  iz=(*b).iup[iz][mu];
            }

            if (iy!=iz)
               itest=4;

            iz=imb[iy];
            iy=(*bb).imb[ix];

            if ((ifc%2)==0)
               iz=idn[iz][mu];
            else
               iz=iup[iz][mu];

            if (iy!=iz)
               itest=5;
         }

         bb+=1;
      }
   }

   error(itest==1,1,"test5 [check1.c]",
         "ipp does not respect the even-odd ordering");
   error(itest==2,1,"test5 [check1.c]",
         "imb does not respect the even-odd ordering");
   error(itest==3,1,"test5 [check1.c]",
         "Partner points are not on the interior boundary of the block");
   error(itest==4,1,"test5 [check1.c]",
         "The block map array is incorrect");
   error(itest==5,1,"test5 [check1.c]",
         "Embedding of the boundary is incorrect");
}


int main(int argc,char *argv[])
{
   int my_rank,bc,igr,bs[4];
   double phi[2],phi_prime[2],theta[3];
   blk_grid_t grid;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Checks on the geometry arrays in the known block grids\n");
      printf("------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n",bs[0],bs[1],bs[2],bs[3]);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>]");
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.38;
   theta[1]=1.57;
   theta[2]=-0.54;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(3);

   geometry();
   set_sap_parms(bs,0,1,1);
   set_dfl_parms(bs,2);
   grid=BLK_GRIDS;

   for (igr=0;igr<(int)(BLK_GRIDS);igr++)
   {
      if (igr==0)
         grid=SAP_BLOCKS;
      else if (igr==1)
         grid=DFL_BLOCKS;
      else
         error_root(1,1,"main [check1.c]","Unknown block grid");

      alloc_bgr(grid);

      test1(grid,bs);
      test2(grid);
      test3(grid);
      test4(grid);
      test5(grid);
   }

   if (my_rank==0)
   {
      printf("No errors detected\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
