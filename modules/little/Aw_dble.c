
/*******************************************************************************
*
* File Aw_dble.c
*
* Copyright (C) 2007, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Application of the double-precision little Wilson-Dirac operator Aw
*
* The externally accessible functions are
*
*   void Aw_dble(complex_dble *v,complex_dble *w)
*     Applies the little Dirac operator to the field v and assigns the
*     result to the field w.
*
*   void Aweeinv_dble(complex_dble *v,complex_dble *w)
*     Applies the inverse of the even-even part of the little Dirac operator
*     to the field v and assigns the result to the field w on the even blocks.
*     On the odd blocks, w is unchanged.
*
*   void Awooinv_dble(complex_dble *v,complex_dble *w)
*     Applies the inverse of the odd-odd part of the little Dirac operator
*     to the field v and assigns the result to the field w on the odd blocks.
*     On the even blocks, w is unchanged.
*
*   void Awoe_dble(complex_dble *v,complex_dble *w)
*     Applies the odd-even part of the little Dirac operator to the field v
*     and assigns the result to the field w on the odd blocks. On the even
*     blocks, w is unchanged.
*
*   void Aweo_dble(complex_dble *v,complex_dble *w)
*     Applies the even-odd part of the little Dirac operator to the field v
*     and *subtracts* the result from the field w on the even blocks. On the 
*     odd blocks, w is unchanged.
*
*   void Awhat_dble(complex_dble *v,complex_dble *w)
*     Applies the even-odd preconditioned little Dirac operator to the field
*     v and assigns the result to the field w on the even blocks. On the odd
*     blocks, w is unchanged.
*
* Notes:
*
* The little Dirac operator and the associated data structures are described
* in the file README.Aw.
*
* The programs Aw_dble(), Awoe_dble() and Aweo_dble() take it for granted
* that the little Dirac operator is up-to-date, while the other programs,
* Aweeinv_dble(), Awooinv_dble() and Awhat_dble(), assume the even-odd
* preconditioned operator to be up-to-date (see Aw_ops.c).
*
* All programs in this module may perform global operations and should be
* called simultaneously on all processes.
*
*******************************************************************************/

#define AW_DBLE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "vflds.h"
#include "linalg.h"
#include "dfl.h"
#include "little.h"
#include "global.h"

static int Ns=0,nb,nbh;
static int nbbh,(*inn)[8];
static complex_dble *vs;


static void alloc_vs(void)
{
   dfl_parms_t dfl;
   dfl_grid_t grd;

   dfl=dfl_parms();
   grd=dfl_geometry();

   Ns=dfl.Ns;
   nb=grd.nb;
   nbh=nb/2;
   nbbh=grd.nbb/2;
   inn=grd.inn;
   
   vs=amalloc(Ns*sizeof(*vs),ALIGN);
   
   error(vs==NULL,1,"alloc_vs [Aw_dble.c]",
         "Unable to allocate auxiliary array");
}


static void apply_Aoe(int *nn,complex_dble **A,complex_dble *v)
{
   int ifc;
   
   cmat_vec_dble(Ns,*A,v+nn[0]*Ns,vs);
   A+=1;

   for (ifc=1;ifc<8;ifc++)
   {
      cmat_vec_assign_dble(Ns,*A,v+nn[ifc]*Ns,vs);
      A+=1;
   }
}


static void apply_Aeo(int *nn,complex_dble **A,complex_dble *v)
{
   int ifc;
   
   for (ifc=0;ifc<8;ifc++)
   {
      cmat_vec_assign_dble(Ns,*A,vs,v+nn[ifc]*Ns);
      A+=1;
   }
}


static void apply_Aee(complex_dble **A,complex_dble *v,complex_dble *w)
{
   complex_dble **Am;

   Am=A+nbh;
   
   for (;A<Am;A++)
   {
      cmat_vec_dble(Ns,*A,v,w);      
      v+=Ns;
      w+=Ns;
   }
}


static void apply_Aoo(complex_dble **A,complex_dble *v,complex_dble *w)
{
   complex_dble **Am;

   Am=A+nbh;
   v+=nbh*Ns;
   w+=nbh*Ns;

   for (;A<Am;A++)
   {
      cmat_vec_dble(Ns,*A,v,w);      
      v+=Ns;
      w+=Ns;
   }
}


void Aw_dble(complex_dble *v,complex_dble *w)
{
   int (*nn)[8],(*nm)[8];
   complex_dble *rv,*rw,*rs,*rm;
   complex_dble **Aeo,**Aoe;
   Aw_dble_t Aw;

   if (Ns==0)
      alloc_vs();
   
   Aw=Awop_dble();
   apply_Aee(Aw.Aee,v,w);
   apply_Aoo(Aw.Aoo,v,w);
   
   if (NPROC>1)
   {
      set_vd2zero(nbbh*Ns,w+nb*Ns);      
      cpvd_int_bnd(v);
   }
   
   Aoe=Aw.Aoe;
   Aeo=Aw.Aeo;
   rv=v+nbh*Ns;
   rw=w+nbh*Ns;
   
   nn=inn+nbh;
   nm=inn+nb;
   
   for (;nn<nm;nn++)
   {
      apply_Aoe(*nn,Aoe,v);

      rs=vs;
      rm=vs+Ns;

      for (;rs<rm;rs++)
      {
         (*rw).re+=(*rs).re;
         (*rw).im+=(*rs).im;

         (*rs).re=(*rv).re;
         (*rs).im=(*rv).im;
         
         rw+=1;
         rv+=1;
      }
      
      apply_Aeo(*nn,Aeo,w);      

      Aoe+=8;
      Aeo+=8;
   }

   if (NPROC>1)
      cpvd_ext_bnd(w);   
}


void Aweeinv_dble(complex_dble *v,complex_dble *w)
{
   Aw_dble_t Aw;

   if (Ns==0)
      alloc_vs();

   Aw=Awophat_dble();
   apply_Aee(Aw.Aee,v,w);
}

   
void Awooinv_dble(complex_dble *v,complex_dble *w)
{
   Aw_dble_t Aw;

   if (Ns==0)
      alloc_vs();
   
   Aw=Awophat_dble();   
   apply_Aoo(Aw.Aoo,v,w);
}


void Awoe_dble(complex_dble *v,complex_dble *w)
{
   int (*nn)[8],(*nm)[8];
   complex_dble *rw,*rs,*rm;
   complex_dble **Aoe;   
   Aw_dble_t Aw;

   if (Ns==0)
      alloc_vs();

   if (NPROC>1)
      cpvd_int_bnd(v);
      
   Aw=Awop_dble();
   Aoe=Aw.Aoe;
   rw=w+nbh*Ns;
   
   nn=inn+nbh;
   nm=inn+nb;
   
   for (;nn<nm;nn++)
   {
      apply_Aoe(*nn,Aoe,v);

      rs=vs;
      rm=vs+Ns;

      for (;rs<rm;rs++)
      {
         (*rw).re=(*rs).re;
         (*rw).im=(*rs).im;
         rw+=1;
      }
      
      Aoe+=8;
   }
}


void Aweo_dble(complex_dble *v,complex_dble *w)
{
   int (*nn)[8],(*nm)[8];
   complex_dble *rv,*rs,*rm;
   complex_dble **Aeo;  
   Aw_dble_t Aw;

   if (Ns==0)
      alloc_vs();

   if (NPROC>1)
      set_vd2zero(nbbh*Ns,w+nb*Ns); 

   Aw=Awop_dble();   
   Aeo=Aw.Aeo;
   rv=v+nbh*Ns;
   
   nn=inn+nbh;
   nm=inn+nb;
   
   for (;nn<nm;nn++)
   {
      rs=vs;
      rm=vs+Ns;

      for (;rs<rm;rs++)
      {
         (*rs).re=-(*rv).re;
         (*rs).im=-(*rv).im;
         rv+=1;
      }
      
      apply_Aeo(*nn,Aeo,w);      

      Aeo+=8;   
   }

   if (NPROC>1)
      cpvd_ext_bnd(w);   
}


void Awhat_dble(complex_dble *v,complex_dble *w)
{
   int (*nn)[8],(*nm)[8];
   complex_dble *rs,*rm;
   complex_dble **Aeo,**Aoe;
   Aw_dble_t Aw;

   if (Ns==0)
      alloc_vs();
   
   assign_vd2vd(nbh*Ns,v,w);
   
   if (NPROC>1)
   {
      set_vd2zero(nbbh*Ns,w+nb*Ns);
      cpvd_int_bnd(v);
   }

   Aw=Awophat_dble();
   Aoe=Aw.Aoe;
   Aeo=Aw.Aeo;
   
   nn=inn+nbh;
   nm=inn+nb;

   for (;nn<nm;nn++)
   {
      apply_Aoe(*nn,Aoe,v);

      rs=vs;
      rm=vs+Ns;

      for (;rs<rm;rs++)
      {
         (*rs).re=-(*rs).re;
         (*rs).im=-(*rs).im;
      }

      apply_Aeo(*nn,Aeo,w);      

      Aoe+=8;
      Aeo+=8;
   }

   if (NPROC>1)
      cpvd_ext_bnd(w);   
}
