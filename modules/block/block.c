
/*******************************************************************************
*
* File block.c
*
* Copyright (C) 2005, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Basic allocation programs for blocks of lattice points.
*
* The externally accessible functions are
*
*   void alloc_blk(block_t *b,int *bo,int *bs,
*                  int iu,int iud,int ns,int nsd)
*     Sets the offset and side-lengths of the block b to bo[4] and bs[4],
*     respectively, and allocates the block fields depending on the values
*     of the other parameters. The single-precision gauge and SW fields are
*     allocated if iu=1, the double-precision gauge and SW fields if iud=1,
*     while ns and nsd are the numbers of single- and double-precision Dirac
*     fields that are allocated. All elements of the block are properly
*     initialized and the share flag b.shf is set to 0x0 (see the notes).
*
*   void alloc_bnd(block_t *b,int iu,int iud,int nw,int nwd)
*     Allocates the boundary structures b.bb in the block b and the fields
*     in there depending on the parameters iu,iud,nw and nwd. The single-
*     and double-precision gauge fields are allocated if iu=1 and iud=1,
*     respectively, while nw and nwd are the numbers of single- and double-
*     precision Weyl fields that are allocated. All elements of the block
*     are then properly initialized (see the notes).
*
*   void clone_blk(block_t *b,int shf,int *bo,block_t *c)
*     Sets the offset of the block c to bo[4] and its side lengths to
*     b.bs[4]. The fields in c are then allocated depending on the bits
*     b1,b2,..,b8 (counting from the lowest) of the share flag shf. The
*     relevant bits are:
*
*       b2=1: b.ipt,b.iup and b.idn are shared,
*       b3=1: b.u, b.bb.u and b.sw are shared,
*       b4=1: b.ud, b.bb.ud and b.swd are shared,
*       b5=1: b.s is shared,
*       b6=1: b.sd is shared.
*       b7=1: b.bb.w is shared,
*       b8=1: b.bb.wd is shared.
*
*     All fields that are not shared and are allocated on b are allocated
*     on c as well, while the pointers to the shared fields are set to those
*     of b. An error occurs if a field is shared according to the share flag
*     b.shf on b but not according to shf. Moreover, the offset differences
*     bo[mu]-b.bo[mu] must be integer multiples of b.bs[mu] for all mu. The
*     share flag c.shf is set to shf.
*
*   void free_blk(block_t *b)
*     Frees the arrays in the block b and in the boundaries b.bb that were
*     previously allocated by alloc_blk(), alloc_bnd() or clone_blk(). The
*     boundary structures are then freed too (if they were allocated) and
*     all entries in the block structure are set to 0 (or NULL).
*
*   int ipt_blk(block_t *b,int *x)
*     Returns the index of the lattice point in the block b with Cartesian
*     coordinates x[4] relative to the base point of b.
*
* Notes:
*
* The entries of the block and boundary structures are explained in the file
* README.block in this directory.
*
* It is currently not possible to allocate blocks that are not fully
* contained in the local lattice. Moreover, the block sizes must be even
* and not smaller than 4. The exterior boundaries of a block may, however,
* overlap with the lattices on the neighbouring processes. In all cases,
* the scalar elements of the structures and the geometry and field arrays
* are properly initialized (gauge and SW fields are set to 1, Dirac spinor
* and Weyl fields to 0).
*
* Block allocation is a global operation, i.e. alloc_blk(), alloc_bnd(),
* clone_blk() and free_blk() must be called on all processes simultaneously.
* The program ipt_blk() can be called locally.
*
* alloc_blk() and clone_blk() register the blocks as being allocated. In this
* way it is possible to exclude any misuses of the programs such as freeing
* an unallocated block (which could have unpredictable side-effects). An
* already allocated block is first freed and then reallocated by alloc_blk().
* Blocks b and their boundary structures b.bb cannot be freed or reallocated
* if the lowest bit of the share flag b.shf is equal to 1.
*
*******************************************************************************/

#define BLOCK_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "sflds.h"
#include "block.h"
#include "global.h"

static const su3 u0={{0.0f}};
static const su3_dble ud0={{0.0}};
static const pauli p0={{0.0f}};
static const pauli_dble pd0={{0.0}};
static const weyl w0={{{0.0f}}};
static const weyl_dble wd0={{{0.0}}};

struct ablk_t
{
   block_t *b;
   struct ablk_t *next;
};

static struct ablk_t *first=NULL;


static int ins_blk(block_t *b)
{
   struct ablk_t *p;

   p=malloc(sizeof(*p));

   if (p!=NULL)
   {
      (*p).b=b;
      (*p).next=first;
      first=p;

      return 0;
   }
   else
      return 1;
}


static void rmv_blk(block_t *b)
{
   struct ablk_t *p,*q;

   q=NULL;

   for (p=first;p!=NULL;p=(*p).next)
   {
      if ((*p).b==b)
      {
         if (q==NULL)
            first=(*p).next;
         else
            (*q).next=(*p).next;

         free(p);
         return;
      }

      q=p;
   }
}


static int fnd_blk(block_t *b)
{
   struct ablk_t *p;

   for (p=first;p!=NULL;p=(*p).next)
   {
      if ((*p).b==b)
         return 1;
   }

   return 0;
}


static void free_bnd(block_t *b)
{
   int shf;
   bndry_t *bb;

   shf=(*b).shf;
   bb=(*b).bb;

   if (bb==NULL)
      return;

   if (!(shf&0x2))
      free((*bb).ipp);

   free((*bb).imb);

   if ((!(shf&0x4))&&((*bb).u!=NULL))
      afree((*bb).u);

   if ((!(shf&0x8))&&((*bb).ud!=NULL))
      afree((*bb).ud);

   if ((!(shf&0x40))&&((*bb).nw>0))
   {
      afree((*bb).w[0]);
      free((*bb).w);
   }

   if ((!(shf&0x80))&&((*bb).nwd>0))
   {
      afree((*bb).wd[0]);
      free((*bb).wd);
   }

   free((*b).bb);
   (*b).bb=NULL;
}


void free_blk(block_t *b)
{
   int shf;

   if (fnd_blk(b)==0)
      return;

   shf=(*b).shf;
   error(shf&0x1,1,"free_blk [block.c]",
         "Protected block");

   free_bnd(b);

   free((*b).bo);
   (*b).bo=NULL;
   (*b).bs=NULL;

   if (!(shf&0x2))
   {
      free((*b).ipt);
      free((*b).iup);
   }

   free((*b).imb);

   (*b).vol=0;
   (*b).vbb=0;
   (*b).nbp=0;
   (*b).shf=0x0;
   (*b).ipt=NULL;
   (*b).imb=NULL;
   (*b).ibp=NULL;
   (*b).iup=NULL;
   (*b).idn=NULL;

   if ((!(shf&0x4))&&((*b).u!=NULL))
   {
      afree((*b).u);
      afree((*b).sw);
   }

   if ((!(shf&0x8))&&((*b).ud!=NULL))
   {
      afree((*b).ud);
      afree((*b).swd);
   }

   if ((!(shf&0x10))&&((*b).ns>0))
   {
      afree((*b).s[0]);
      free((*b).s);
   }

   if ((!(shf&0x20))&&((*b).nsd>0))
   {
      afree((*b).sd[0]);
      free((*b).sd);
   }

   (*b).ns=0;
   (*b).nsd=0;
   (*b).u=NULL;
   (*b).ud=NULL;
   (*b).sw=NULL;
   (*b).swd=NULL;
   (*b).s=NULL;
   (*b).sd=NULL;

   rmv_blk(b);
}


static void set_u2unity(int vol,su3 *u)
{
   su3 unity,*um;

   unity=u0;
   unity.c11.re=1.0f;
   unity.c22.re=1.0f;
   unity.c33.re=1.0f;

   um=u+vol;

   for (;u<um;u++)
      (*u)=unity;
}


static void set_ud2unity(int vol,su3_dble *ud)
{
   su3_dble unity,*um;

   unity=ud0;
   unity.c11.re=1.0;
   unity.c22.re=1.0;
   unity.c33.re=1.0;

   um=ud+vol;

   for (;ud<um;ud++)
      (*ud)=unity;
}


static void set_sw2unity(int vol,pauli *p)
{
   pauli unity,*pm;

   unity=p0;
   unity.u[0]=1.0f;
   unity.u[1]=1.0f;
   unity.u[2]=1.0f;
   unity.u[3]=1.0f;
   unity.u[4]=1.0f;
   unity.u[5]=1.0f;

   pm=p+vol;

   for (;p<pm;p++)
      (*p)=unity;
}


static void set_swd2unity(int vol,pauli_dble *pd)
{
   pauli_dble unity,*pm;

   unity=pd0;
   unity.u[0]=1.0;
   unity.u[1]=1.0;
   unity.u[2]=1.0;
   unity.u[3]=1.0;
   unity.u[4]=1.0;
   unity.u[5]=1.0;

   pm=pd+vol;

   for (;pd<pm;pd++)
      (*pd)=unity;
}


static void set_w2zero(int vol,weyl *w)
{
   weyl *wm;

   wm=w+vol;

   for (;w<wm;w++)
      (*w)=w0;
}


static void set_wd2zero(int vol,weyl_dble *wd)
{
   weyl_dble *wm;

   wm=wd+vol;

   for (;wd<wm;wd++)
      (*wd)=wd0;
}


static void new_blk(block_t *b,int *bo,int *bs,
                    int iu,int iud,int ns,int nsd,int shf)
{
   int n,mu;

   free_blk(b);

   error_root(
      (bo[0]<0)||((bo[0]+bs[0])>L0)||(bs[0]<4)||((bs[0]%2)!=0)||
      (bo[1]<0)||((bo[1]+bs[1])>L1)||(bs[1]<4)||((bs[1]%2)!=0)||
      (bo[2]<0)||((bo[2]+bs[2])>L2)||(bs[2]<4)||((bs[2]%2)!=0)||
      (bo[3]<0)||((bo[3]+bs[3])>L3)||(bs[3]<4)||((bs[3]%2)!=0),1,
      "new_blk [block.c]","Improper choice of block position or size");

   error_root((ns<0)||(nsd<0),1,"new_blk [block.c]",
              "Improper choice of the numbers of spinor fields");

   (*b).bo=malloc(8*sizeof(*(*b).bo));
   error((*b).bo==NULL,1,"new_blk [block.c]",
         "Unable to allocate size arrays");
   (*b).bs=(*b).bo+4;

   for (mu=0;mu<4;mu++)
   {
      (*b).bo[mu]=bo[mu];
      (*b).bs[mu]=bs[mu];
   }

   (*b).vol=bs[0]*bs[1]*bs[2]*bs[3];
   (*b).vbb=2*(bs[0]*bs[1]*bs[2]+bs[1]*bs[2]*bs[3]+
               bs[2]*bs[3]*bs[0]+bs[3]*bs[0]*bs[1]);
   (*b).nbp=0;

   if ((cpr[0]==0)&&(bo[0]==0)&&(bc_type()!=3))
      (*b).nbp+=bs[1]*bs[2]*bs[3];
   if ((cpr[0]==(NPROC0-1))&&((bo[0]+bs[0])==L0)&&(bc_type()==0))
      (*b).nbp+=bs[1]*bs[2]*bs[3];

   (*b).ns=ns;
   (*b).nsd=nsd;
   (*b).shf=shf;

   if (shf&0x2)
   {
      (*b).ipt=NULL;
      (*b).iup=NULL;
      (*b).idn=NULL;
   }
   else
   {
      (*b).ipt=malloc(((*b).vol+1)*sizeof(*(*b).ipt));
      (*b).iup=malloc(2*(*b).vol*sizeof(*(*b).iup));
      error(((*b).ipt==NULL)||((*b).iup==NULL),1,
            "new_blk [block.c]","Unable to allocate the geometry arrays");
      (*b).idn=(*b).iup+(*b).vol;
   }

   (*b).imb=malloc((((*b).vol+1)+(*b).nbp)*sizeof(*(*b).imb));
   (*b).ibp=(*b).imb+(*b).vol+1;

   if ((shf&0x4)||(iu!=1))
   {
      (*b).u=NULL;
      (*b).sw=NULL;
   }
   else
   {
      (*b).u=amalloc(4*(*b).vol*sizeof(*(*b).u),ALIGN);
      (*b).sw=amalloc(2*(*b).vol*sizeof(*(*b).sw),ALIGN);
      error(((*b).u==NULL)||((*b).sw==NULL),1,"new_blk [block.c]",
            "Unable to allocate the single-precision gauge field");
      set_u2unity(4*(*b).vol,(*b).u);
      set_sw2unity(2*(*b).vol,(*b).sw);
   }

   if ((shf&0x8)||(iud!=1))
   {
      (*b).ud=NULL;
      (*b).swd=NULL;
   }
   else
   {
      (*b).ud=amalloc(4*(*b).vol*sizeof(*(*b).ud),ALIGN);
      (*b).swd=amalloc(2*(*b).vol*sizeof(*(*b).swd),ALIGN);
      error(((*b).ud==NULL)||((*b).swd==NULL),1,"new_blk [block.c]",
            "Unable to allocate the double-precision gauge field");
      set_ud2unity(4*(*b).vol,(*b).ud);
      set_swd2unity(2*(*b).vol,(*b).swd);
   }

   if ((shf&0x10)||(ns==0))
      (*b).s=NULL;
   else
   {
      (*b).s=malloc(ns*sizeof(*(*b).s));
      error((*b).s==NULL,1,"new_blk [block.c]",
            "Unable to allocate the single-precision spinor fields");

      (*b).s[0]=amalloc(ns*((*b).vol+1)*sizeof(*((*b).s[0])),ALIGN);
      error((*b).s[0]==NULL,2,"new_blk [block.c]",
            "Unable to allocate the single-precision spinor fields");

      for (n=1;n<ns;n++)
         (*b).s[n]=(*b).s[n-1]+(*b).vol+1;

      set_s2zero(ns*((*b).vol+1),(*b).s[0]);
   }

   if ((shf&0x20)||(nsd==0))
      (*b).sd=NULL;
   else
   {
      (*b).sd=malloc(nsd*sizeof(*(*b).sd));
      error((*b).sd==NULL,1,"new_blk [block.c]",
            "Unable to allocate the double-precision spinor fields");

      (*b).sd[0]=amalloc(nsd*((*b).vol+1)*sizeof(*((*b).sd[0])),ALIGN);
      error((*b).sd[0]==NULL,2,"new_blk [block.c]",
            "Unable to allocate the pointer array (*b).sd");

      for (n=1;n<nsd;n++)
         (*b).sd[n]=(*b).sd[n-1]+(*b).vol+1;

      set_sd2zero(nsd*((*b).vol+1),(*b).sd[0]);
   }

   (*b).bb=NULL;

   error(ins_blk(b),1,"new_blk [block.c]",
         "Unable to register allocated block");
}


void alloc_blk(block_t *b,int *bo,int *bs,
               int iu,int iud,int ns,int nsd)
{
   int iprms[12],mu,ie;

   if (NPROC>1)
   {
      for (mu=0;mu<4;mu++)
      {
         iprms[mu]=bo[mu];
         iprms[4+mu]=bs[mu];
      }

      iprms[8]=iu;
      iprms[9]=iud;
      iprms[10]=ns;
      iprms[11]=nsd;

      MPI_Bcast(iprms,12,MPI_INT,0,MPI_COMM_WORLD);

      ie=0;

      for (mu=0;mu<4;mu++)
         if ((iprms[mu]!=bo[mu])||(iprms[4+mu]!=bs[mu]))
            ie=1;

      error((ie)||(iprms[8]!=iu)||(iprms[9]!=iud)||
            (iprms[10]!=ns)||(iprms[11]!=nsd),1,"alloc_blk [block.c]",
            "Parameters are not global");
   }

   error(iup[0][0]==0,1,"alloc_blk [block.c]",
         "The global geometry arrays are not set");

   new_blk(b,bo,bs,iu,iud,ns,nsd,0x0);
   blk_geometry(b);
   blk_imbed(b);
}


static void new_bnd(block_t *b,int iu,int iud,int nw,int nwd,int shf)
{
   int vol,ifc,n;
   int *bs,*ipp,*map,*imb;
   su3 *u;
   su3_dble *ud;
   weyl **w,*wb;
   weyl_dble **wd,*wdb;
   bndry_t *bb;

   error_root((nw<0)||(nwd<0),1,"new_bnd [block.c]",
              "Improper choice of the numbers of Weyl fields");

   free_bnd(b);
   bb=malloc(8*sizeof(*bb));
   error(bb==NULL,1,"new_bnd [block.c]",
         "Unable to allocate boundary structures");
   (*b).bb=bb;

   vol=(*b).vol;
   bs=(*b).bs;

   for (ifc=0;ifc<8;ifc++)
   {
      bb[ifc].ifc=ifc;
      bb[ifc].vol=vol/bs[ifc/2];
      bb[ifc].nw=nw;
      bb[ifc].nwd=nwd;
   }

   vol=(*b).vbb;

   if (shf&0x2)
   {
      for (ifc=0;ifc<8;ifc++)
      {
         bb[ifc].ipp=NULL;
         bb[ifc].map=NULL;
      }
   }
   else
   {
      ipp=malloc(2*(vol+8)*sizeof(*ipp));
      error(ipp==NULL,1,"new_bnd [block.c]",
            "Unable to allocate the geometry arrays");
      map=ipp+vol+8;

      for (ifc=0;ifc<8;ifc++)
      {
         bb[ifc].ipp=ipp;
         ipp+=(bb[ifc].vol+1);
         bb[ifc].map=map;
         map+=(bb[ifc].vol+1);
      }
   }

   imb=malloc((vol+8)*sizeof(*imb));
   error(imb==NULL,2,"new_bnd [block.c]",
         "Unable to allocate the geometry arrays");

   for (ifc=0;ifc<8;ifc++)
   {
      bb[ifc].imb=imb;
      imb+=(bb[ifc].vol+1);
   }

   if ((shf&0x4)||(iu!=1))
   {
      for (ifc=0;ifc<8;ifc++)
         bb[ifc].u=NULL;
   }
   else
   {
      u=amalloc(vol*sizeof(*u),ALIGN);
      error(u==NULL,1,"new_bnd [block.c]",
            "Unable to allocate the single-precision gauge field");
      set_u2unity(vol,u);

      for (ifc=0;ifc<8;ifc++)
      {
         bb[ifc].u=u;
         u+=bb[ifc].vol;
      }
   }

   if ((shf&0x8)||(iud!=1))
   {
      for (ifc=0;ifc<8;ifc++)
         bb[ifc].ud=NULL;
   }
   else
   {
      ud=amalloc(vol*sizeof(*ud),ALIGN);
      error(ud==NULL,1,"new_bnd [block.c]",
            "Unable to allocate the double-precision gauge field");
      set_ud2unity(vol,ud);

      for (ifc=0;ifc<8;ifc++)
      {
         bb[ifc].ud=ud;
         ud+=bb[ifc].vol;
      }
   }

   if ((shf&0x40)||(nw==0))
   {
      for (ifc=0;ifc<8;ifc++)
         bb[ifc].w=NULL;
   }
   else
   {
      w=malloc(8*nw*sizeof(*w));
      wb=amalloc(nw*vol*sizeof(*wb),ALIGN);
      error((w==NULL)||(wb==NULL),1,"new_bnd [block.c]",
            "Unable to allocate the single-precision Weyl fields");
      set_w2zero(nw*vol,wb);

      for (ifc=0;ifc<8;ifc++)
      {
         bb[ifc].w=w;

         for (n=0;n<nw;n++)
         {
            (*w)=wb+n*vol;
            w+=1;
         }

         wb+=bb[ifc].vol;
      }
   }

   if ((shf&0x80)||(nwd==0))
   {
      for (ifc=0;ifc<8;ifc++)
         bb[ifc].wd=NULL;
   }
   else
   {
      wd=malloc(8*nwd*sizeof(*wd));
      wdb=amalloc(nwd*vol*sizeof(*wdb),ALIGN);
      error((wd==NULL)||(wdb==NULL),1,"new_bnd [block.c]",
            "Unable to allocate the double-precision Weyl fields");
      set_wd2zero(nwd*vol,wdb);

      for (ifc=0;ifc<8;ifc++)
      {
         bb[ifc].wd=wd;

         for (n=0;n<nwd;n++)
         {
            (*wd)=wdb+n*vol;
            wd+=1;
         }

         wdb+=bb[ifc].vol;
      }
   }
}


void alloc_bnd(block_t *b,int iu,int iud,int nw,int nwd)
{
   int iprms[12],mu,ie;
   int *bo,*bs;

   error(fnd_blk(b)==0,1,"alloc_bnd [block.c]",
         "Block is not allocated");
   error((*b).shf&0x1,1,"alloc_bnd [block.c]",
         "Protected block");

   if (NPROC>1)
   {
      bo=(*b).bo;
      bs=(*b).bs;

      for (mu=0;mu<4;mu++)
      {
         iprms[mu]=bo[mu];
         iprms[4+mu]=bs[mu];
      }

      iprms[8]=iu;
      iprms[9]=iud;
      iprms[10]=nw;
      iprms[11]=nwd;

      MPI_Bcast(iprms,12,MPI_INT,0,MPI_COMM_WORLD);

      ie=0;

      for (mu=0;mu<4;mu++)
         if ((iprms[mu]!=bo[mu])||(iprms[4+mu]!=bs[mu]))
            ie=1;

      error((ie)||(iprms[8]!=iu)||(iprms[9]!=iud)||
            (iprms[10]!=nw)||(iprms[11]!=nwd),1,"alloc_bnd [block.c]",
            "Parameters are not global");
   }

   new_bnd(b,iu,iud,nw,nwd,0x0);
   bnd_geometry(b);
   bnd_imbed(b);
}


void clone_blk(block_t *b,int shf,int *bo,block_t *c)
{
   int iprms[23],mu,ie;
   int *bbo,*bs,bshf;
   int iu,iud,ns,nsd,iub,iudb,nw,nwd;
   int ib,ifc;

   error(fnd_blk(b)==0,1,"clone_blk [block.c]",
         "The block to be cloned is not allocated");

   bbo=(*b).bo;
   bs=(*b).bs;
   bshf=(*b).shf;
   iu=((*b).u!=NULL);
   iud=((*b).ud!=NULL);
   ns=(*b).ns;
   nsd=(*b).nsd;

   if ((*b).bb!=NULL)
   {
      iub=((*b).bb[0].u!=NULL);
      iudb=((*b).bb[0].ud!=NULL);
      nw=(*b).bb[0].nw;
      nwd=(*b).bb[0].nwd;
      ib=1;
   }
   else
   {
      iub=0;
      iudb=0;
      nw=0;
      nwd=0;
      ib=0;
   }

   if (NPROC>1)
   {
      for (mu=0;mu<4;mu++)
      {
         iprms[mu]=bbo[mu];
         iprms[4+mu]=bs[mu];
         iprms[8+mu]=bo[mu];
      }

      iprms[12]=bshf;
      iprms[13]=iu;
      iprms[14]=iud;
      iprms[15]=ns;
      iprms[16]=nsd;
      iprms[17]=iub;
      iprms[18]=iudb;
      iprms[19]=nw;
      iprms[20]=nwd;
      iprms[21]=ib;
      iprms[22]=shf;

      MPI_Bcast(iprms,23,MPI_INT,0,MPI_COMM_WORLD);

      ie=0;

      for (mu=0;mu<4;mu++)
      {
         if ((iprms[mu]!=bbo[mu])||
             (iprms[4+mu]!=bs[mu])||
             (iprms[8+mu]!=bo[mu]))
            ie=1;
      }

      error((ie)||(iprms[12]!=bshf)||(iprms[13]!=iu)||(iprms[14]!=iud)||
            (iprms[15]!=ns)||(iprms[16]!=nsd)||(iprms[17]!=iub)||
            (iprms[18]!=iudb)||(iprms[19]!=nw)||(iprms[20]!=nwd)||
            (iprms[21]!=ib)||(iprms[22]!=shf),1,"clone_blk [block.c]",
            "Parameters are not global");
   }

   error_root((bo[0]<0)||((bo[0]+bs[0])>L0)||((abs(bo[0]-bbo[0])%bs[0])!=0)||
              (bo[1]<0)||((bo[1]+bs[1])>L1)||((abs(bo[1]-bbo[1])%bs[1])!=0)||
              (bo[2]<0)||((bo[2]+bs[2])>L2)||((abs(bo[2]-bbo[2])%bs[2])!=0)||
              (bo[3]<0)||((bo[3]+bs[3])>L3)||((abs(bo[3]-bbo[3])%bs[3])!=0),1,
              "clone_blk [block.c]","Improper block offset");

   error_root(((bshf&0x2)&&(!(shf&0x2)))||
              ((bshf&0x4)&&(!(shf&0x4))&&(iu!=0))||
              ((bshf&0x8)&&(!(shf&0x8))&&(iud!=0))||
              ((bshf&0x10)&&(!(shf&0x10))&&(ns>0))||
              ((bshf&0x20)&&(!(shf&0x20))&&(nsd>0)),1,
              "clone_blk [block.c]","Share flag mismatch");

   new_blk(c,bo,bs,iu,iud,ns,nsd,shf);

   if (shf&0x2)
   {
      (*c).ipt=(*b).ipt;
      (*c).iup=(*b).iup;
      (*c).idn=(*b).idn;
   }

   if ((shf&0x4)&&(iu!=0))
   {
      (*c).u=(*b).u;
      (*c).sw=(*b).sw;
   }

   if ((shf&0x8)&&(iud!=0))
   {
      (*c).ud=(*b).ud;
      (*c).swd=(*b).swd;
   }

   if ((shf&0x10)&&(ns>0))
      (*c).s=(*b).s;

   if ((shf&0x20)&&(nsd>0))
      (*c).sd=(*b).sd;

   if (!(shf&0x2))
      blk_geometry(c);
   blk_imbed(c);

   if (ib)
   {
      error_root(((bshf&0x4)&&(!(shf&0x4))&&(iub!=0))||
                 ((bshf&0x8)&&(!(shf&0x8))&&(iudb!=0))||
                 ((bshf&0x40)&&(!(shf&0x40))&&(nw>0))||
                 ((bshf&0x80)&&(!(shf&0x80))&&(nwd>0)),2,
                 "clone_blk [block.c]","Share flag mismatch");

      new_bnd(c,iub,iudb,nw,nwd,shf);

      for (ifc=0;ifc<8;ifc++)
      {
         if (shf&0x2)
         {
            (*c).bb[ifc].ipp=(*b).bb[ifc].ipp;
            (*c).bb[ifc].map=(*b).bb[ifc].map;
         }

         if ((shf&0x4)&&(iub!=0))
            (*c).bb[ifc].u=(*b).bb[ifc].u;

         if ((shf&0x8)&&(iudb!=0))
            (*c).bb[ifc].ud=(*b).bb[ifc].ud;

         if ((shf&0x40)&&(nw>0))
            (*c).bb[ifc].w=(*b).bb[ifc].w;

         if ((shf&0x80)&&(nwd>0))
            (*c).bb[ifc].wd=(*b).bb[ifc].wd;
      }

      if (!(shf&0x2))
         bnd_geometry(c);
      bnd_imbed(c);
   }
}


int ipt_blk(block_t *b,int *x)
{
   int *bs,n,ix;

   bs=(*b).bs;

   n=((x[0]<0)||(x[0]>=bs[0]));
   ix=x[0];

   n|=((x[1]<0)||(x[1]>=bs[1]));
   ix=x[1]+bs[1]*ix;

   n|=((x[2]<0)||(x[2]>=bs[2]));
   ix=x[2]+bs[2]*ix;

   n|=((x[3]<0)||(x[3]>=bs[3]));
   ix=x[3]+bs[3]*ix;

   if (n==0)
      return (*b).ipt[ix];
   else
   {
      error_loc(1,1,"ipt_blk [block.c]","Point coordinates are out of range");
      return 0;
   }
}
