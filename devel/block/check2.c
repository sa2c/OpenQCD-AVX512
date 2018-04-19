
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Checks on the allocation and initialization of the gauge, Dirac and Weyl
* fields on the known block grids.
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
#include "block.h"
#include "global.h"

typedef union
{
   su3 u;
   float r[18];
} umat_t;

typedef union
{
   su3_dble u;
   double r[18];
} umat_dble_t;

typedef union
{
   spinor s;
   float r[24];
} spin_t;

typedef union
{
   spinor_dble s;
   double r[24];
} spin_dble_t;

typedef union
{
   weyl w;
   float r[12];
} wspin_t;

typedef union
{
   weyl_dble w;
   double r[12];
} wspin_dble_t;


static int check_u(int vol,su3 *u)
{
   int i;
   umat_t *m,*mm;

   m=(umat_t*)(u);
   mm=m+vol;

   for (;m<mm;m++)
   {
      for (i=0;i<18;i++)
      {
         if ((((i%8)==0)&&((*m).r[i]!=1.0f))||
             (((i%8)!=0)&&((*m).r[i]!=0.0f)))
            return 1;
      }
   }

   return 0;
}


static int check_ud(int vol,su3_dble *u)
{
   int i;
   umat_dble_t *m,*mm;

   m=(umat_dble_t*)(u);
   mm=m+vol;

   for (;m<mm;m++)
   {
      for (i=0;i<18;i++)
      {
         if ((((i%8)==0)&&((*m).r[i]!=1.0))||
             (((i%8)!=0)&&((*m).r[i]!=0.0)))
            return 1;
      }
   }

   return 0;
}


static int check_sw(int vol,pauli *sw)
{
   int i;
   pauli *sm;

   sm=sw+vol;

   for (;sw<sm;sw++)
   {
      for (i=0;i<36;i++)
      {
         if (((i<6)&&((*sw).u[i]!=1.0f))||
             (((i>=6)&&((*sw).u[i]!=0.0f))))
            return 1;
      }
   }

   return 0;
}


static int check_swd(int vol,pauli_dble *swd)
{
   int i;
   pauli_dble *sm;

   sm=swd+vol;

   for (;swd<sm;swd++)
   {
      for (i=0;i<36;i++)
      {
         if (((i<6)&&((*swd).u[i]!=1.0))||
             (((i>=6)&&((*swd).u[i]!=0.0))))
            return 1;
      }
   }

   return 0;
}


static int check_s(int ns,int vol,spinor **s)
{
   int k,i;
   spin_t *sp,*sm;

   for (k=0;k<ns;k++)
   {
      sp=(spin_t*)(s[k]);
      sm=sp+vol;

      for (;sp<sm;sp++)
      {
         for (i=0;i<24;i++)
         {
            if ((*sp).r[i]!=0.0f)
               return 1;
         }
      }
   }

   return 0;
}


static int check_sd(int nsd,int vol,spinor_dble **sd)
{
   int k,i;
   spin_dble_t *sp,*sm;

   for (k=0;k<nsd;k++)
   {
      sp=(spin_dble_t*)(sd[k]);
      sm=sp+vol;

      for (;sp<sm;sp++)
      {
         for (i=0;i<24;i++)
         {
            if ((*sp).r[i]!=0.0)
               return 1;
         }
      }
   }

   return 0;
}


static int check_w(int nw,int vol,weyl **w)
{
   int k,i;
   wspin_t *wp,*wm;

   for (k=0;k<nw;k++)
   {
      wp=(wspin_t*)(w[k]);
      wm=wp+vol;

      for (;wp<wm;wp++)
      {
         for (i=0;i<12;i++)
         {
            if ((*wp).r[i]!=0.0f)
               return 1;
         }
      }
   }

   return 0;
}


static int check_wd(int nwd,int vol,weyl_dble **wd)
{
   int k,i;
   wspin_dble_t *wp,*wm;

   for (k=0;k<nwd;k++)
   {
      wp=(wspin_dble_t*)(wd[k]);
      wm=wp+vol;

      for (;wp<wm;wp++)
      {
         for (i=0;i<12;i++)
         {
            if ((*wp).r[i]!=0.0)
               return 1;
         }
      }
   }

   return 0;
}


static int check_blk(block_t *b0,block_t *b,
                     int iu,int iud,int ns,int nsd,int shf)
{
   int vol;

   vol=(*b).vol;

   if ((*b).shf!=shf)
      return 1;

   if (shf&0x2)
   {
      if (((*b0).ipt!=(*b).ipt)||((*b0).iup!=(*b).iup)||
          ((*b0).idn!=(*b).idn))
         return 2;
   }

   if (iu==1)
   {
      if (((*b).u==NULL)||((shf&0x4)&&((*b0).u!=(*b).u)))
         return 3;
      if (check_u(4*vol,(*b).u))
         return 3;

      if (((*b).sw==NULL)||((shf&0x4)&&((*b0).sw!=(*b).sw)))
         return 3;
      if (check_sw(2*vol,(*b).sw))
         return 3;
   }
   else
   {
      if (((*b).u!=NULL)||((*b).sw!=NULL))
         return 2;
   }

   if (iud==1)
   {
      if (((*b).ud==NULL)||((shf&0x8)&&((*b0).ud!=(*b).ud)))
         return 4;
      if (check_ud(4*vol,(*b).ud))
         return 4;

      if (((*b).swd==NULL)||((shf&0x8)&&((*b0).swd!=(*b).swd)))
         return 4;
      if (check_swd(2*vol,(*b).swd))
         return 4;
   }
   else
   {
      if (((*b).ud!=NULL)||((*b).swd!=NULL))
         return 3;
   }

   if ((*b).ns!=ns)
      return 5;

   if (ns>0)
   {
      if (((*b).s==NULL)||((shf&0x10)&&((*b0).s!=(*b).s)))
         return 5;
      if (check_s(ns,vol,(*b).s))
         return 5;
   }
   else
   {
      if ((*b).s!=NULL)
         return 5;
   }

   if ((*b).nsd!=nsd)
      return 6;

   if (nsd>0)
   {
      if (((*b).sd==NULL)||((shf&0x20)&&((*b0).sd!=(*b).sd)))
         return 6;
      if (check_sd(nsd,vol,(*b).sd))
         return 6;
   }
   else
   {
      if ((*b).sd!=NULL)
         return 6;
   }

   return 0;
}


static int check_bnd(block_t *b0,block_t *b,
                     int iub,int iudb,int nw,int nwd,int shf)
{
   int vol,ifc;
   bndry_t *bb0,*bb;

   bb0=(*b).bb;
   bb=(*b).bb;

   for (ifc=0;ifc<8;ifc++)
   {
      vol=(*bb).vol;

      if (iub==1)
      {
         if (((*bb).u==NULL)||((shf&0x4)&&((*bb0).u!=(*bb).u)))
            return 7;
         if (check_u(vol,(*bb).u))
            return 7;
      }
      else
      {
         if ((*bb).u!=NULL)
            return 7;
      }

      if (iudb==1)
      {
         if (((*bb).ud==NULL)||((shf&0x8)&&((*bb0).ud!=(*bb).ud)))
            return 8;
         if (check_ud(vol,(*bb).ud))
            return 8;
      }
      else
      {
         if ((*bb).ud!=NULL)
            return 8;
      }

      if ((*bb).nw!=nw)
         return 9;

      if (nw>0)
      {
         if (((*bb).w==NULL)||((shf&0x40)&&((*bb0).w!=(*bb).w)))
            return 9;
         if (check_w(nw,vol,(*bb).w))
            return 9;
      }
      else
      {
         if ((*bb).w!=NULL)
            return 9;
      }

      if ((*bb).nwd!=nwd)
         return 10;

      if (nwd>0)
      {
         if (((*bb).wd==NULL)||((shf&0x80)&&((*bb0).wd!=(*bb).wd)))
            return 10;
         if (check_wd(nwd,vol,(*bb).wd))
            return 10;
      }
      else
      {
         if ((*bb).wd!=NULL)
            return 10;
      }

      bb0+=1;
      bb+=1;
   }

   return 0;
}

int main(int argc,char *argv[])
{
   int my_rank,n,n0,n1,n2,n3;
   int igr,bs[4],nb,isw,itest;
   int iu,iud,ns,nsd;
   int iub,iudb,nw,nwd;
   int shg,shu,shud,shs,shsd,shw,shwd,shf;
   block_t *b0,*b;
   blk_grid_t grid;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Checks on the allocation and initialization of the gauge,\n"
             "Dirac and Weyl fields on the known block grids.\n");
      printf("---------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n\n",bs[0],bs[1],bs[2],bs[3]);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   geometry();
   set_sap_parms(bs,0,1,1);
   set_dfl_parms(bs,2);
   grid=BLK_GRIDS;

   for (igr=0;igr<(int)(BLK_GRIDS);igr++)
   {
      iu=0;
      iud=0;
      ns=0;
      nsd=0;
      iub=0;
      iudb=0;
      nw=0;
      nwd=0;

      shg=1;
      shu=0;
      shud=0;
      shs=0;
      shsd=0;
      shw=0;
      shwd=0;

      if (igr==0)
      {
         grid=SAP_BLOCKS;

         iu=1;
         ns=3;
         nw=1;
         iub=1;
         shs=1;
      }
      else if (igr==1)
      {
         grid=DFL_BLOCKS;

         iud=1;
         ns=3;
         nsd=3;
         shud=1;
      }
      else
         error_root(1,1,"main [check2.c]","Unknown block grid");

      shf=0x1|(shg<<1)|(shu<<2)|(shud<<3)|(shs<<4)|(shsd<<5)|(shw<<6)|(shwd<<7);
      alloc_bgr(grid);
      print_grid_flags(grid);
      b0=blk_list(grid,&nb,&isw);

      n0=L0/bs[0];
      n1=L1/bs[1];
      n2=L2/bs[2];
      n3=L3/bs[3];
      n=n0*cpr[0]+n1*cpr[1]+n2*cpr[2]+n3*cpr[3];

      error((b0==NULL)||(nb!=(n0*n1*n2*n3))||(isw!=(n%2)),1,
            "main [check2.c]","Incorrect return values of blk_list");

      if (my_rank==0)
      {
         printf("Share flag on the blocks = %#x\n",(*b0).shf);
         printf("Should be                  %#x\n\n",shf);
      }

      itest=0;

      for (n=0;n<nb;n++)
      {
         b=b0+n;
         error(((*b).bs[0]!=bs[0])||((*b).bs[1]!=bs[1])||
               ((*b).bs[2]!=bs[2])||((*b).bs[3]!=bs[3]),1,
               "main [check2.c]","b.bs is incorrect");

         error(((*b).bo[0]<0)||(((*b).bo[0]+bs[0])>L0)||
               ((*b).bo[1]<0)||(((*b).bo[1]+bs[1])>L1)||
               ((*b).bo[2]<0)||(((*b).bo[2]+bs[2])>L2)||
               ((*b).bo[3]<0)||(((*b).bo[3]+bs[3])>L3),1,
               "main [check2.c]","b.bo is out of range");

         error((((*b).bo[0]%bs[0])!=0)||(((*b).bo[1]%bs[1])!=0)||
               (((*b).bo[2]%bs[2])!=0)||(((*b).bo[3]%bs[3])!=0),1,
               "main [check2.c]","b.bo is not an integer multiple of bs");

         n0=(*b).bo[0]/bs[0];
         n1=(*b).bo[1]/bs[1];
         n2=(*b).bo[2]/bs[2];
         n3=(*b).bo[3]/bs[3];

         isw=(n0+n1+n2+n3)%2;

         error(((isw==0)&&(n>=(nb/2)))||((isw==1)&&(n<(nb/2))),1,
               "main [check2.c]","Blocks are not locally even-odd ordered");

         itest=check_blk(b0,b,iu,iud,ns,nsd,shf);
         if (itest!=0)
            break;

         error((*b).bb==NULL,1,"main [check2.c]",
               "Block boundaries are not allocated");

         itest=check_bnd(b0,b,iub,iudb,nw,nwd,shf);
         if (itest!=0)
            break;
      }

      error(itest==1,1,"main [check2.c]","Unexpected share flag");
      error(itest==2,1,"main [check2.c]","Geometry arrays are not shared");
      error(itest==3,1,"main [check2.c]",
            "b.u or b.sw is not in the proper condition");
      error(itest==4,1,"main [check2.c]",
            "b.ud or b.swd is not in the proper condition");
      error(itest==5,1,"main [check2.c]",
            "b.s is not in the proper condition");
      error(itest==6,1,"main [check2.c]",
            "b.sd is not in the proper condition");
      error(itest==7,1,"main [check2.c]",
            "b.bb.u is not in the proper condition");
      error(itest==8,1,"main [check2.c]",
            "b.bb.ud is not in the proper condition");
      error(itest==9,1,"main [check2.c]",
            "b.bb.w is not in the proper condition");
      error(itest==10,1,"main [check2.c]",
            "b.bb.wd is not in the proper condition");
   }

   if (my_rank==0)
   {
      printf("No errors detected\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
