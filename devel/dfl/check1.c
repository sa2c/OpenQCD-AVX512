
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2007, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the DFL_BLOCKS grid geometry arrays.
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
#include "dfl.h"
#include "global.h"


int main(int argc,char *argv[])
{
   int my_rank,bc,bs[4],nbs,isw;
   int nb,nbb,*nbbe,*nbbo,*obbe,*obbo;
   int (*inn)[8],*idx,*ipp,*map;
   int ix,iy,iz,ifc,ie;
   int *bo1,*bo2;
   int l[4],mu,is;
   double phi[2],phi_prime[2],theta[3];
   block_t *b;
   dfl_grid_t dfl_grid;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Check of the DFL_BLOCKS grid geometry arrays\n");
      printf("--------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n\n",bs[0],bs[1],bs[2],bs[3]);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>]");
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

   geometry();
   set_dfl_parms(bs,4);
   dfl_grid=dfl_geometry();
   nb=dfl_grid.nb;
   nbbe=dfl_grid.nbbe;
   nbbo=dfl_grid.nbbo;
   obbe=dfl_grid.obbe;
   obbo=dfl_grid.obbo;

   alloc_bgr(DFL_BLOCKS);
   b=blk_list(DFL_BLOCKS,&nbs,&isw);

   error((bs[0]!=(*b).bs[0])||(bs[1]!=(*b).bs[1])||(bs[2]!=(*b).bs[2])||
         (bs[3]!=(*b).bs[3])||(nb!=nbs),1,"main [check1.c]",
         "Block sizes bs are incorrectly set or incorrect block number");

   ie=0;
   nbb=(nbbe[0]+nbbo[0]);

   if (obbe[0]!=0)
      ie=1;
   if (obbo[0]!=(obbe[7]+nbbe[7]))
      ie=2;

   for (ifc=1;ifc<8;ifc++)
   {
      nbb+=(nbbe[ifc]+nbbo[ifc]);

      if (obbe[ifc]!=(obbe[ifc-1]+nbbe[ifc-1]))
         ie=1;
      if (obbo[ifc]!=(obbo[ifc-1]+nbbo[ifc-1]))
         ie=2;
   }

   error(nbb!=dfl_grid.nbb,1,"main [check1.c]","nbb is incorrect");
   error(ie==1,1,"main [check1.c]","Incorrect offsets obbe[ifc]");
   error(ie==2,1,"main [check1.c]","Incorrect offsets obbo[ifc]");

   inn=dfl_grid.inn;
   idx=dfl_grid.idx;
   ipp=dfl_grid.ipp;
   map=dfl_grid.map;
   iz=0;

   for (ifc=0;ifc<8;ifc++)
   {
      for (ix=obbe[ifc];ix<(obbe[ifc]+nbbe[ifc]);ix++)
      {
         iy=ipp[ix];

         if ((ix>obbe[ifc])&&(iy<=iz))
            ie=1;

         if (inn[iy][ifc]!=(nb+ix))
            ie=3;

         iz=iy;
      }

      for (ix=obbo[ifc];ix<(obbo[ifc]+nbbo[ifc]);ix++)
      {
         iy=ipp[ix];

         if ((ix>obbo[ifc])&&(iy<=iz))
            ie=2;

         if (inn[iy][ifc]!=(nb+ix))
            ie=3;

         iz=iy;
      }
   }

   error(ie==1,1,"main [check1.c]","Incorrect ipp at even boundary points");
   error(ie==2,1,"main [check1.c]","Incorrect ipp at odd boundary points");
   error(ie==3,1,"main [check1.c]","ipp and inn are inconsistent");

   for (ix=0;ix<nb;ix++)
   {
      if (idx[idx[ix]]!=ix)
         ie=1;

      if (((ix>0)&&(ix<(nb/2)))||(ix>(nb/2)))
      {
         if (idx[ix]!=idx[ix-1]+1)
            ie=2;
      }

      if (((ix==0)&&(isw==0))||((ix==(nb/2))&&(isw==1)))
      {
         bo1=b[idx[ix]].bo;

         for (mu=0;mu<4;mu++)
         {
            if (bo1[mu]!=0)
               ie=3;
         }
      }
   }

   error(ie==1,1,"main [check1.c]","Index array idx[ix] is not involutive");
   error(ie==2,1,"main [check1.c]","The ordering of idx[ix] is incorrect");
   error(ie==3,1,"main [check1.c]","Index of the first block is incorrect ");

   for (ix=0;ix<nb;ix++)
   {
      for (ifc=0;ifc<8;ifc++)
      {
         iy=inn[ix][ifc];

         if ((iy<0)||(iy>=(nb+nbb)))
            ie=1;
         else
         {
            if (iy<nb)
               iz=inn[iy][ifc^0x1];
            else
               iz=ipp[iy-nb];

            if (iz!=ix)
               ie=2;
         }
      }
   }

   error(ie==1,1,"main [check1.c]","Index inn[ix][ifc] is out of range");
   error(ie==2,1,"main [check1.c]","Neighbouring blocks are not paired");

   l[0]=L0;
   l[1]=L1;
   l[2]=L2;
   l[3]=L3;

   for (ix=0;ix<nb;ix++)
   {
      bo1=b[idx[ix]].bo;

      for (ifc=0;ifc<8;ifc++)
      {
         iy=inn[ix][ifc];

         if (iy<nb)
         {
            bo2=b[idx[iy]].bo;

            for (mu=0;mu<4;mu++)
            {
               if (mu!=(ifc/2))
               {
                  if (bo2[mu]!=bo1[mu])
                     ie=1;
               }
               else
               {
                  is=2*(ifc&0x1)-1;
                  is=(bo1[mu]+is*bs[mu]+l[mu])%l[mu];

                  if (bo2[mu]!=is)
                     ie=1;
               }
            }
         }
         else
         {
            mu=ifc/2;

            if ((((ifc&0x1)==1)&&((bo1[mu]+bs[mu])!=l[mu]))||
                (((ifc&0x1)==0)&&(bo1[mu]!=0)))
               ie=2;

            iy=map[iy-nb];
            bo2=b[idx[iy]].bo;

            for (mu=0;mu<4;mu++)
            {
               if (mu!=(ifc/2))
               {
                  if (bo2[mu]!=bo1[mu])
                     ie=3;
               }
               else
               {
                  if ((((ifc&0x1)==1)&&(bo2[mu]!=0))||
                      (((ifc&0x1)==0)&&(bo2[mu]!=(l[mu]-bs[mu]))))
                     ie=3;
               }
            }
         }
      }
   }

   error(ie==1,1,"main [check1.c]",
         "Index array inn[ix][ifc] is incorrect in the bulk");
   error(ie==2,1,"main [check1.c]",
         "Index array inn[ix][ifc] is incorrect at the boundary");
   error(ie==3,1,"main [check1.c]",
         "Index array map[ix] is incorrect");

   if (my_rank==0)
   {
      printf("No errors detected\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
