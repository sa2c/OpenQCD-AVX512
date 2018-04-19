
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2005, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of assign_s2sblk(),...,assign_sdblk2sd().
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
#include "sflds.h"
#include "sw_term.h"
#include "block.h"
#include "global.h"

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


static int cmp_s(spinor *r,spinor *s)
{
   int i;
   spin_t *rr,*rs;

   rr=(spin_t*)(r);
   rs=(spin_t*)(s);

   for (i=0;i<24;i++)
   {
      if ((*rr).r[i]!=(*rs).r[i])
         return 1;
   }

   return 0;
}


static int cmp_sd(spinor_dble *r,spinor_dble *s)
{
   int i;
   spin_dble_t *rr,*rs;

   rr=(spin_dble_t*)(r);
   rs=(spin_dble_t*)(s);

   for (i=0;i<24;i++)
   {
      if ((*rr).r[i]!=(*rs).r[i])
         return 1;
   }

   return 0;
}


static int check_sb(block_t *b,ptset_t set,int k,spinor *s)
{
   int x0,x1,x2,x3,x[4];
   int y0,y1,y2,y3;
   int ix,iy,is,n0,n1;

   for (x0=0;x0<(*b).bs[0];x0++)
   {
      for (x1=0;x1<(*b).bs[1];x1++)
      {
         for (x2=0;x2<(*b).bs[2];x2++)
         {
            for (x3=0;x3<(*b).bs[3];x3++)
            {
               x[0]=x0;
               x[1]=x1;
               x[2]=x2;
               x[3]=x3;

               y0=(*b).bo[0]+x0;
               y1=(*b).bo[1]+x1;
               y2=(*b).bo[2]+x2;
               y3=(*b).bo[3]+x3;

               ix=ipt_blk(b,x);
               iy=ipt[y3+L3*y2+L2*L3*y1+L1*L2*L3*y0];
               is=(y0+y1+y2+y3)%2;

               n0=((is==0)&&((set==ALL_PTS)||(set==EVEN_PTS)));
               n1=((is==1)&&((set==ALL_PTS)||(set==ODD_PTS)));

               if ((n0==1)||(n1==1))
               {
                  if  (cmp_s((*b).s[k]+ix,s+iy))
                     return 1;
               }
            }
         }
      }
   }

   return 0;
}


static int check_sdb(block_t *b,ptset_t set,int k,spinor_dble *sd)
{
   int x0,x1,x2,x3,x[4],n0,n1;
   int y0,y1,y2,y3;
   int ix,iy,is;

   for (x0=0;x0<(*b).bs[0];x0++)
   {
      for (x1=0;x1<(*b).bs[1];x1++)
      {
         for (x2=0;x2<(*b).bs[2];x2++)
         {
            for (x3=0;x3<(*b).bs[3];x3++)
            {
               x[0]=x0;
               x[1]=x1;
               x[2]=x2;
               x[3]=x3;

               y0=(*b).bo[0]+x0;
               y1=(*b).bo[1]+x1;
               y2=(*b).bo[2]+x2;
               y3=(*b).bo[3]+x3;

               ix=ipt_blk(b,x);
               iy=ipt[y3+L3*y2+L2*L3*y1+L1*L2*L3*y0];
               is=(y0+y1+y2+y3)%2;

               n0=((is==0)&&((set==ALL_PTS)||(set==EVEN_PTS)));
               n1=((is==1)&&((set==ALL_PTS)||(set==ODD_PTS)));

               if ((n0==1)||(n1==1))
               {
                  if (cmp_sd((*b).sd[k]+ix,sd+iy))
                     return 1;
               }
            }
         }
      }
   }

   return 0;
}


static int diff_s(int vol,spinor *s,spinor *r)
{
   spinor *sm;

   sm=s+vol;

   for (;s<sm;s++)
   {
      if (cmp_s(s,r))
         return 1;

      r+=1;
   }

   return 0;
}


static int diff_sd(int vol,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm;

   sm=s+vol;

   for (;s<sm;s++)
   {
      if (cmp_sd(s,r))
         return 1;

      r+=1;
   }

   return 0;
}


int main(int argc,char *argv[])
{
   int my_rank,nb,isw,iset;
   int bs[4],n,k,l,ns,nsd,vol;
   spinor **ps;
   spinor_dble **psd;
   ptset_t set;
   block_t *b;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check5.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Check of assign_s2sblk(),...,assign_sdblk2sd()\n");
      printf("----------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n",bs[0],bs[1],bs[2],bs[3]);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);

   start_ranlux(0,1234);
   geometry();
   alloc_ws(2);
   alloc_wsd(2);
   set_sap_parms(bs,0,1,1);
   set_dfl_parms(bs,2);
   alloc_bgr(SAP_BLOCKS);
   alloc_bgr(DFL_BLOCKS);

   ps=reserve_ws(2);
   psd=reserve_wsd(2);

   random_s(VOLUME,ps[0],1.0f);
   random_sd(VOLUME,psd[0],1.0);
   assign_s2s(VOLUME,ps[0],ps[1]);
   assign_sd2sd(VOLUME,psd[0],psd[1]);

   b=blk_list(SAP_BLOCKS,&nb,&isw);
   vol=(*b).vol;
   ns=(*b).ns;
   k=0;

   for (n=0;n<nb;n++)
   {
      for (iset=0;iset<=(int)(PT_SETS);iset++)
      {
         if (iset==0)
            set=ALL_PTS;
         else if (iset==1)
            set=EVEN_PTS;
         else if (iset==2)
            set=ODD_PTS;
         else
            set=NO_PTS;

         k=((k+1)%ns);
         l=((k+1)%ns);

         random_s(vol,(*b).s[k],0.1f);
         assign_s2s(vol,(*b).s[k],(*b).s[l]);
         assign_s2sblk(SAP_BLOCKS,n,set,ps[0],k);
         error(check_sb(b,set,k,ps[0]),1,
               "main [check5.c]","assign_s2sblk() is incorrect");

         if ((set==EVEN_PTS)||(set==NO_PTS))
            error(diff_s(vol/2,(*b).s[k]+vol/2,(*b).s[l]+vol/2),2,
                  "main [check5.c]",
                  "Unexpected change of the block field by assign_s2sblk()");
         if ((set==ODD_PTS)||(set==NO_PTS))
            error(diff_s(vol/2,(*b).s[k],(*b).s[l]),3,
                  "main [check5.c]",
                  "Unexpected change of the block field by assign_s2sblk()");

         error(diff_s(VOLUME,ps[0],ps[1]),1,
               "main [check5.c]","assign_s2sblk() changes the input field");

         random_s(vol,(*b).s[k],0.1f);
         assign_s2s(vol,(*b).s[k],(*b).s[l]);
         assign_sblk2s(SAP_BLOCKS,n,set,k,ps[0]);
         error(check_sb(b,set,k,ps[0]),1,
               "main [check5.c]","assign_sblk2s() is incorrect");

         error(diff_s(vol,(*b).s[k],(*b).s[l]),1,
               "main [check5.c]","assign_sblk2s() changes the input field");
         assign_s2sblk(SAP_BLOCKS,n,set,ps[1],k);
         assign_sblk2s(SAP_BLOCKS,n,set,k,ps[0]);

         error(diff_s(VOLUME,ps[0],ps[1]),1,
               "main [check5.c]",
               "Unexpected change of the global field by assign_s2sblk()");
      }

      b+=1;
   }

   b=blk_list(DFL_BLOCKS,&nb,&isw);
   vol=(*b).vol;
   nsd=(*b).nsd;
   k=0;

   for (n=0;n<nb;n++)
   {
      for (iset=0;iset<(int)(PT_SETS);iset++)
      {
         if (iset==0)
            set=ALL_PTS;
         else if (iset==1)
            set=EVEN_PTS;
         else if (iset==2)
            set=ODD_PTS;
         else
            set=NO_PTS;

         k=((k+1)%nsd);
         l=((k+1)%nsd);

         random_sd(vol,(*b).sd[k],0.1f);
         assign_sd2sd(vol,(*b).sd[k],(*b).sd[l]);
         assign_sd2sdblk(DFL_BLOCKS,n,set,psd[0],k);
         error(check_sdb(b,set,k,psd[0]),1,
               "main [check5.c]","assign_sd2sdblk() is incorrect");

         if ((set==EVEN_PTS)||(set==NO_PTS))
            error(diff_sd(vol/2,(*b).sd[k]+vol/2,(*b).sd[l]+vol/2),2,
                  "main [check5.c]",
                  "Unexpected change of the block field by assign_sd2sdblk()");
         if ((set==ODD_PTS)||(set==NO_PTS))
            error(diff_sd(vol/2,(*b).sd[k],(*b).sd[l]),3,
                  "main [check5.c]",
                  "Unexpected change of the block field by assign_sd2sdblk()");

         error(diff_sd(VOLUME,psd[0],psd[1]),1,
               "main [check5.c]","assign_sd2sdblk() changes the input field");

         random_sd(vol,(*b).sd[k],0.1f);
         assign_sd2sd(vol,(*b).sd[k],(*b).sd[l]);
         assign_sdblk2sd(DFL_BLOCKS,n,set,k,psd[0]);
         error(check_sdb(b,set,k,psd[0]),1,
               "main [check5.c]","assign_sdblk2sd() is incorrect");

         error(diff_sd(vol,(*b).sd[k],(*b).sd[l]),1,
               "main [check5.c]","assign_sdblk2sd() changes the input field");
         assign_sd2sdblk(DFL_BLOCKS,n,set,psd[1],k);
         assign_sdblk2sd(DFL_BLOCKS,n,set,k,psd[0]);

         error(diff_sd(VOLUME,psd[0],psd[1]),1,
               "main [check5.c]",
               "Unexpected change of the global field by assign_sd2sdblk()");
      }

      b+=1;
   }

   random_s(VOLUME,ps[0],1.0f);
   assign_s2s(VOLUME,ps[0],ps[1]);
   assign_s2sd(VOLUME,ps[0],psd[0]);

   b=blk_list(DFL_BLOCKS,&nb,&isw);
   vol=(*b).vol;
   nsd=(*b).nsd;
   k=0;

   for (n=0;n<nb;n++)
   {
      for (iset=0;iset<(int)(PT_SETS);iset++)
      {
         if (iset==0)
            set=ALL_PTS;
         else if (iset==1)
            set=EVEN_PTS;
         else if (iset==2)
            set=ODD_PTS;
         else
            set=NO_PTS;

         k=((k+1)%nsd);
         l=((k+1)%nsd);

         random_sd(vol,(*b).sd[k],0.1f);
         assign_sd2sd(vol,(*b).sd[k],(*b).sd[l]);
         assign_s2sdblk(DFL_BLOCKS,n,set,ps[0],k);
         error(check_sdb(b,set,k,psd[0]),1,
               "main [check5.c]","assign_s2sdblk() is incorrect");

         if ((set==EVEN_PTS)||(set==NO_PTS))
            error(diff_sd(vol/2,(*b).sd[k]+vol/2,(*b).sd[l]+vol/2),2,
                  "main [check5.c]",
                  "Unexpected change of the block field by assign_s2sdblk()");
         if ((set==ODD_PTS)||(set==NO_PTS))
            error(diff_sd(vol/2,(*b).sd[k],(*b).sd[l]),3,
                  "main [check5.c]",
                  "Unexpected change of the block field by assign_s2sdblk()");

         error(diff_s(VOLUME,ps[0],ps[1]),1,
               "main [check5.c]","assign_s2sdblk() changes the input field");
      }

      b+=1;
   }

   if (my_rank==0)
   {
      printf("\n");
      printf("No errors detected\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
