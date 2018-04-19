
/*******************************************************************************
*
* File Pbnd.c
*
* Copyright (C) 2005, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic programs for the projector theta to the exterior boundary of a 
* block of lattice points (version for single-precision fields)
*
* The following are arrays of functions indexed by the face number
* ifc=0,..,7
*
*   void (*assign_s2w[8])(int *imb,int vol,spinor *s,weyl *r)
*     Applies the projector theta[ifc] to the spinor s[imb[ix]], 
*     ix=0,..,vol-1, and assigns the result to the weyl spinor r[ix].
*
*   void (*add_assign_w2s[8])(int *imb,int vol,weyl *s,spinor *r)
*     Expands the Weyl spinor s[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and adds psi to r[imb[ix]].
*
*   void (*sub_assign_w2s[8])(int *imb,int vol,weyl *s,spinor *r)
*     Expands the Weyl spinor s[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and subtracts psi from r[imb[ix]].
*
*   void (*mulg5_sub_assign_w2s[8])(int *imb,int vol,weyl *s,spinor *r)
*     Expands the Weyl spinor s[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and subtracts gamma5*psi from
*     r[imb[ix]].
*
* Notes:
*
* The projector theta was introduced in section 3.3 of 
*
*   Lattice QCD and the Schwarz alternating procedure, JHEP 0305 (2003) 052
*
* Block faces in the -0,+0,..,-3,+3 directions are labelled by an index
* ifc=0,..,7 and the projector on a given face is then given by
*
*   theta[ifc] = (1/2)*(1+gamma_mu) if ifc=2*mu,
*
*              = (1/2)*(1-gamma_mu) if ifc=2*mu+1
*  
* Dirac fields on the face that satisfy theta[ifc]*psi=psi are completely
* characterized by their first two Dirac components. In this way they are
* mapped to Weyl fields on the face in an invertible manner.
*
* The size and position of the faces is only implicitly defined through
* the parameter vol and the array *imb of the indices of the points on
* the face.
*
* None of these programs involves communications. They are general purpose
* routines that know nothing about the underiying geometry. In particular,
* they can be called locally. If SSE instructions are used, the fields must
* be aligned to a 16 byte boundary.
*
*******************************************************************************/

#define PBND_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "su3.h"
#include "sflds.h"

#if (defined x64)
#include "sse2.h"

#define _load_cst(c) \
__asm__ __volatile__ ("movss %0, %%xmm15 \n\t" \
                      "shufps $0x0, %%xmm15, %%xmm15" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm15")

#define _mul_cst() \
__asm__ __volatile__ ("mulps %%xmm15, %%xmm0 \n\t" \
                      "mulps %%xmm15, %%xmm1 \n\t" \
                      "mulps %%xmm15, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

static const float poh=0.5f;


static void assign_s2w0(int *imb,int vol,spinor *s,weyl *r)
{
   weyl *rm;
   spinor *si,*sin;

   _load_cst(poh);
   rm=r+vol;   
   si=s+(*imb);
   imb+=(r<(rm-1));
   sin=s+(*imb);

   for (;r<rm;r++)
   {
      _sse_pair_load((*si).c1,(*si).c2);
      _sse_pair_load_up((*si).c3,(*si).c4);               

      si=sin;
      imb+=(r<(rm-2));
      sin=s+(*imb);
      _prefetch_spinor(sin);      

      _sse_vector_sub();
      _mul_cst();
      _sse_pair_store((*r).c1,(*r).c2);
   }   
}


static void assign_s2w1(int *imb,int vol,spinor *s,weyl *r)
{
   weyl *rm;
   spinor *si,*sin;

   _load_cst(poh);   
   rm=r+vol;   
   si=s+(*imb);
   imb+=(r<(rm-1));
   sin=s+(*imb);

   for (;r<rm;r++)
   {
      _sse_pair_load((*si).c1,(*si).c2);
      _sse_pair_load_up((*si).c3,(*si).c4);      

      si=sin;
      imb+=(r<(rm-2));
      sin=s+(*imb);
      _prefetch_spinor(sin);      
      
      _sse_vector_add();
      _mul_cst();
      _sse_pair_store((*r).c1,(*r).c2);
   }   
}


static void assign_s2w2(int *imb,int vol,spinor *s,weyl *r)
{
   weyl *rm;
   spinor *si,*sin;

   _load_cst(poh);
   rm=r+vol;   
   si=s+(*imb);
   imb+=(r<(rm-1));
   sin=s+(*imb);

   for (;r<rm;r++)
   {
      _sse_pair_load((*si).c1,(*si).c2);
      _sse_pair_load_up((*si).c4,(*si).c3);      

      si=sin;
      imb+=(r<(rm-2));
      sin=s+(*imb);
      _prefetch_spinor(sin);      
      
      _sse_vector_i_sub();
      _mul_cst();
      _sse_pair_store((*r).c1,(*r).c2);
   }   
}


static void assign_s2w3(int *imb,int vol,spinor *s,weyl *r)
{
   weyl *rm;
   spinor *si,*sin;

   _load_cst(poh);   
   rm=r+vol;   
   si=s+(*imb);
   imb+=(r<(rm-1));
   sin=s+(*imb);

   for (;r<rm;r++)
   {
      _sse_pair_load((*si).c1,(*si).c2);
      _sse_pair_load_up((*si).c4,(*si).c3);      

      si=sin;
      imb+=(r<(rm-2));
      sin=s+(*imb);
      _prefetch_spinor(sin);      
      
      _sse_vector_i_add();
      _mul_cst();
      _sse_pair_store((*r).c1,(*r).c2);
   }   
}


static void assign_s2w4(int *imb,int vol,spinor *s,weyl *r)
{
   weyl *rm;
   spinor *si,*sin;

   _load_cst(poh);   
   rm=r+vol;   
   si=s+(*imb);
   imb+=(r<(rm-1));
   sin=s+(*imb);

   for (;r<rm;r++)
   {
      _sse_pair_load((*si).c1,(*si).c2);
      _sse_pair_load_up((*si).c4,(*si).c3);      

      si=sin;
      imb+=(r<(rm-2));
      sin=s+(*imb);
      _prefetch_spinor(sin);      
      
      _sse_vector_subadd();
      _mul_cst();
      _sse_pair_store((*r).c1,(*r).c2);      
   }   
}


static void assign_s2w5(int *imb,int vol,spinor *s,weyl *r)
{
   weyl *rm;
   spinor *si,*sin;

   _load_cst(poh);   
   rm=r+vol;   
   si=s+(*imb);
   imb+=(r<(rm-1));
   sin=s+(*imb);

   for (;r<rm;r++)
   {
      _sse_pair_load((*si).c1,(*si).c2);
      _sse_pair_load_up((*si).c4,(*si).c3);      

      si=sin;
      imb+=(r<(rm-2));
      sin=s+(*imb);
      _prefetch_spinor(sin);      
      
      _sse_vector_addsub();
      _mul_cst();
      _sse_pair_store((*r).c1,(*r).c2);
   }   
}


static void assign_s2w6(int *imb,int vol,spinor *s,weyl *r)
{
   weyl *rm;
   spinor *si,*sin;

   _load_cst(poh);   
   rm=r+vol;   
   si=s+(*imb);
   imb+=(r<(rm-1));
   sin=s+(*imb);

   for (;r<rm;r++)
   {
      _sse_pair_load((*si).c1,(*si).c2);
      _sse_pair_load_up((*si).c3,(*si).c4);      

      si=sin;
      imb+=(r<(rm-2));
      sin=s+(*imb);
      _prefetch_spinor(sin);      
      
      _sse_vector_i_subadd();
      _mul_cst();
      _sse_pair_store((*r).c1,(*r).c2);
   }   
}


static void assign_s2w7(int *imb,int vol,spinor *s,weyl *r)
{
   weyl *rm;
   spinor *si,*sin;

   _load_cst(poh);   
   rm=r+vol;   
   si=s+(*imb);
   imb+=(r<(rm-1));
   sin=s+(*imb);

   for (;r<rm;r++)
   {
      _sse_pair_load((*si).c1,(*si).c2);
      _sse_pair_load_up((*si).c3,(*si).c4);      

      si=sin;
      imb+=(r<(rm-2));
      sin=s+(*imb);
      _prefetch_spinor(sin);      
      
      _sse_vector_i_addsub();
      _mul_cst();
      _sse_pair_store((*r).c1,(*r).c2);
   }   
}


static void add_assign_w2s0(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_add();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void add_assign_w2s1(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_add();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_add();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void add_assign_w2s2(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_add();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_add();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void add_assign_w2s3(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_add();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_sub();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void add_assign_w2s4(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_add();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_subadd();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void add_assign_w2s5(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_add();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_addsub();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void add_assign_w2s6(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_add();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_addsub();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void add_assign_w2s7(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_add();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_subadd();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void sub_assign_w2s0(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      

      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;

      _sse_vector_add();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void sub_assign_w2s1(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      

      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;

      _sse_vector_sub();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void sub_assign_w2s2(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;

      _sse_vector_i_sub();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void sub_assign_w2s3(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_add();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void sub_assign_w2s4(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_addsub();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void sub_assign_w2s5(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_subadd();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void sub_assign_w2s6(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_subadd();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void sub_assign_w2s7(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_addsub();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void mulg5_sub_assign_w2s0(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {   
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      

      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void mulg5_sub_assign_w2s1(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {      
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_add();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void mulg5_sub_assign_w2s2(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {      
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_add();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void mulg5_sub_assign_w2s3(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {      
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_sub();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void mulg5_sub_assign_w2s4(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {      
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_subadd();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void mulg5_sub_assign_w2s5(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {      
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c4,(*ri).c3);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_addsub();
      _sse_pair_store((*ri).c4,(*ri).c3);
   }
}


static void mulg5_sub_assign_w2s6(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {      
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_addsub();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}


static void mulg5_sub_assign_w2s7(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri,*rin,*rim;

   sm=s+vol;
   rin=r+(*imb);
   imb+=(s<(sm-1));
   rim=r+(*imb);

   for (;s<sm;)
   {      
      _sse_pair_load_up((*s).c1,(*s).c2);
      _sse_pair_load((*rin).c1,(*rin).c2);      

      ri=rin;
      rin=rim;
      imb+=(s<(sm-2));
      rim=r+(*imb);
      _prefetch_spinor(rim);      
      
      _sse_vector_sub();
      _sse_pair_store((*ri).c1,(*ri).c2);
      _sse_pair_load((*ri).c3,(*ri).c4);      

      s+=4;
      _prefetch_weyl(s);
      s-=3;
      
      _sse_vector_i_subadd();
      _sse_pair_store((*ri).c3,(*ri).c4);
   }
}

#else

static void assign_s2w0(int *imb,int vol,spinor *s,weyl *r)
{
   float r1;
   weyl *rm;
   spinor *si;

   r1=0.5f;
   rm=r+vol;

   for (;r<rm;r++)
   {
      si=s+(*imb);
      imb+=1;
      _vector_sub((*r).c1,(*si).c1,(*si).c3);
      _vector_sub((*r).c2,(*si).c2,(*si).c4);       
      _vector_mul((*r).c1,r1,(*r).c1);
      _vector_mul((*r).c2,r1,(*r).c2);      
   }
}


static void assign_s2w1(int *imb,int vol,spinor *s,weyl *r)
{
   float r1;
   weyl *rm;
   spinor *si;

   r1=0.5f;
   rm=r+vol;

   for (;r<rm;r++)
   {
      si=s+(*imb);
      imb+=1;

      _vector_add((*r).c1,(*si).c1,(*si).c3);
      _vector_add((*r).c2,(*si).c2,(*si).c4);       
      _vector_mul((*r).c1,r1,(*r).c1);
      _vector_mul((*r).c2,r1,(*r).c2);       
   }
} 


static void assign_s2w2(int *imb,int vol,spinor *s,weyl *r)
{
   float r1;
   weyl *rm;
   spinor *si;

   r1=0.5f;
   rm=r+vol;

   for (;r<rm;r++)
   {
      si=s+(*imb);
      imb+=1;

      _vector_i_sub((*r).c1,(*si).c1,(*si).c4);
      _vector_i_sub((*r).c2,(*si).c2,(*si).c3);       
      _vector_mul((*r).c1,r1,(*r).c1);
      _vector_mul((*r).c2,r1,(*r).c2);      
   }
} 

static void assign_s2w3(int *imb,int vol,spinor *s,weyl *r)
{
   float r1;
   weyl *rm;
   spinor *si;

   r1=0.5f;
   rm=r+vol;

   for (;r<rm;r++)
   {
      si=s+(*imb);
      imb+=1;

      _vector_i_add((*r).c1,(*si).c1,(*si).c4);
      _vector_i_add((*r).c2,(*si).c2,(*si).c3);       
      _vector_mul((*r).c1,r1,(*r).c1);
      _vector_mul((*r).c2,r1,(*r).c2);         
   }
} 


static void assign_s2w4(int *imb,int vol,spinor *s,weyl *r)
{
   float r1;
   weyl *rm;
   spinor *si;

   r1=0.5f;
   rm=r+vol;

   for (;r<rm;r++)
   {
      si=s+(*imb);
      imb+=1;

      _vector_sub((*r).c1,(*si).c1,(*si).c4);
      _vector_add((*r).c2,(*si).c2,(*si).c3);       
      _vector_mul((*r).c1,r1,(*r).c1);
      _vector_mul((*r).c2,r1,(*r).c2);         
   }
} 


static void assign_s2w5(int *imb,int vol,spinor *s,weyl *r)
{
   float r1;
   weyl *rm;
   spinor *si;

   r1=0.5f;
   rm=r+vol;

   for (;r<rm;r++)
   {
      si=s+(*imb);
      imb+=1;

      _vector_add((*r).c1,(*si).c1,(*si).c4);
      _vector_sub((*r).c2,(*si).c2,(*si).c3);       
      _vector_mul((*r).c1,r1,(*r).c1);
      _vector_mul((*r).c2,r1,(*r).c2);      
   }
} 


static void assign_s2w6(int *imb,int vol,spinor *s,weyl *r)
{
   float r1;
   weyl *rm;
   spinor *si;

   r1=0.5f;
   rm=r+vol;

   for (;r<rm;r++)
   {
      si=s+(*imb);
      imb+=1;

      _vector_i_sub((*r).c1,(*si).c1,(*si).c3);
      _vector_i_add((*r).c2,(*si).c2,(*si).c4);       
      _vector_mul((*r).c1,r1,(*r).c1);
      _vector_mul((*r).c2,r1,(*r).c2);      
   }
} 


static void assign_s2w7(int *imb,int vol,spinor *s,weyl *r)
{
   float r1;
   weyl *rm;
   spinor *si;

   r1=0.5f;
   rm=r+vol;

   for (;r<rm;r++)
   {
      si=s+(*imb);
      imb+=1;

      _vector_i_add((*r).c1,(*si).c1,(*si).c3);
      _vector_i_sub((*r).c2,(*si).c2,(*si).c4);       
      _vector_mul((*r).c1,r1,(*r).c1);
      _vector_mul((*r).c2,r1,(*r).c2);        
   }
} 


static void add_assign_w2s0(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;
      _vector_add_assign((*ri).c1,(*s).c1);
      _vector_add_assign((*ri).c2,(*s).c2);
      _vector_sub_assign((*ri).c3,(*s).c1);
      _vector_sub_assign((*ri).c4,(*s).c2);
   }
}


static void add_assign_w2s1(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_add_assign((*ri).c1,(*s).c1);
      _vector_add_assign((*ri).c2,(*s).c2);
      _vector_add_assign((*ri).c3,(*s).c1);
      _vector_add_assign((*ri).c4,(*s).c2);         
   }
}


static void add_assign_w2s2(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_add_assign((*ri).c1,(*s).c1);
      _vector_add_assign((*ri).c2,(*s).c2);
      _vector_i_add_assign((*ri).c3,(*s).c2);
      _vector_i_add_assign((*ri).c4,(*s).c1);       
   }
}


static void add_assign_w2s3(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_add_assign((*ri).c1,(*s).c1);
      _vector_add_assign((*ri).c2,(*s).c2);
      _vector_i_sub_assign((*ri).c3,(*s).c2);
      _vector_i_sub_assign((*ri).c4,(*s).c1);      
   }
}


static void add_assign_w2s4(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_add_assign((*ri).c1,(*s).c1);
      _vector_add_assign((*ri).c2,(*s).c2);
      _vector_add_assign((*ri).c3,(*s).c2);
      _vector_sub_assign((*ri).c4,(*s).c1);      
   }
}


static void add_assign_w2s5(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_add_assign((*ri).c1,(*s).c1);
      _vector_add_assign((*ri).c2,(*s).c2);
      _vector_sub_assign((*ri).c3,(*s).c2);
      _vector_add_assign((*ri).c4,(*s).c1);      
   }
}


static void add_assign_w2s6(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_add_assign((*ri).c1,(*s).c1);
      _vector_add_assign((*ri).c2,(*s).c2);
      _vector_i_add_assign((*ri).c3,(*s).c1);
      _vector_i_sub_assign((*ri).c4,(*s).c2);      
   }
}


static void add_assign_w2s7(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_add_assign((*ri).c1,(*s).c1);
      _vector_add_assign((*ri).c2,(*s).c2);
      _vector_i_sub_assign((*ri).c3,(*s).c1);
      _vector_i_add_assign((*ri).c4,(*s).c2);      
   }
}


static void sub_assign_w2s0(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_add_assign((*ri).c3,(*s).c1);
      _vector_add_assign((*ri).c4,(*s).c2);
   }
}


static void sub_assign_w2s1(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_sub_assign((*ri).c3,(*s).c1);
      _vector_sub_assign((*ri).c4,(*s).c2);      
   }
}


static void sub_assign_w2s2(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_i_sub_assign((*ri).c3,(*s).c2);
      _vector_i_sub_assign((*ri).c4,(*s).c1);       
   }
}


static void sub_assign_w2s3(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_i_add_assign((*ri).c3,(*s).c2);
      _vector_i_add_assign((*ri).c4,(*s).c1);
   }
}


static void sub_assign_w2s4(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_sub_assign((*ri).c3,(*s).c2);
      _vector_add_assign((*ri).c4,(*s).c1);      
   }
}


static void sub_assign_w2s5(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_add_assign((*ri).c3,(*s).c2);
      _vector_sub_assign((*ri).c4,(*s).c1);      
   }
}


static void sub_assign_w2s6(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_i_sub_assign((*ri).c3,(*s).c1);
      _vector_i_add_assign((*ri).c4,(*s).c2);      
   }
}


static void sub_assign_w2s7(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_i_add_assign((*ri).c3,(*s).c1);
      _vector_i_sub_assign((*ri).c4,(*s).c2);      
   }
}


static void mulg5_sub_assign_w2s0(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_sub_assign((*ri).c3,(*s).c1);
      _vector_sub_assign((*ri).c4,(*s).c2);
   }
}


static void mulg5_sub_assign_w2s1(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_add_assign((*ri).c3,(*s).c1);
      _vector_add_assign((*ri).c4,(*s).c2);         
   }
}


static void mulg5_sub_assign_w2s2(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_i_add_assign((*ri).c3,(*s).c2);
      _vector_i_add_assign((*ri).c4,(*s).c1);       
   }
}


static void mulg5_sub_assign_w2s3(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_i_sub_assign((*ri).c3,(*s).c2);
      _vector_i_sub_assign((*ri).c4,(*s).c1);      
   }
}


static void mulg5_sub_assign_w2s4(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_add_assign((*ri).c3,(*s).c2);
      _vector_sub_assign((*ri).c4,(*s).c1);      
   }
}


static void mulg5_sub_assign_w2s5(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_sub_assign((*ri).c3,(*s).c2);
      _vector_add_assign((*ri).c4,(*s).c1);      
   }
}


static void mulg5_sub_assign_w2s6(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_i_add_assign((*ri).c3,(*s).c1);
      _vector_i_sub_assign((*ri).c4,(*s).c2);      
   }
}


static void mulg5_sub_assign_w2s7(int *imb,int vol,weyl *s,spinor *r)
{
   weyl *sm;
   spinor *ri;

   sm=s+vol;
   
   for (;s<sm;s++)
   {
      ri=r+(*imb);
      imb+=1;

      _vector_sub_assign((*ri).c1,(*s).c1);
      _vector_sub_assign((*ri).c2,(*s).c2);
      _vector_i_sub_assign((*ri).c3,(*s).c1);
      _vector_i_add_assign((*ri).c4,(*s).c2);      
   }
}

#endif

void (*assign_s2w[8])(int *imb,int vol,spinor *s,weyl *r) =
{assign_s2w0,assign_s2w1,assign_s2w2,assign_s2w3,
 assign_s2w4,assign_s2w5,assign_s2w6,assign_s2w7};

void (*add_assign_w2s[8])(int *imb,int vol,weyl *s,spinor *r) =
{add_assign_w2s0,add_assign_w2s1,add_assign_w2s2,add_assign_w2s3,
 add_assign_w2s4,add_assign_w2s5,add_assign_w2s6,add_assign_w2s7};

void (*sub_assign_w2s[8])(int *imb,int vol,weyl *s,spinor *r) =
{sub_assign_w2s0,sub_assign_w2s1,sub_assign_w2s2,sub_assign_w2s3,
 sub_assign_w2s4,sub_assign_w2s5,sub_assign_w2s6,sub_assign_w2s7};

void (*mulg5_sub_assign_w2s[8])(int *imb,int vol,weyl *s,spinor *r) =
{mulg5_sub_assign_w2s0,mulg5_sub_assign_w2s1,
 mulg5_sub_assign_w2s2,mulg5_sub_assign_w2s3,
 mulg5_sub_assign_w2s4,mulg5_sub_assign_w2s5,
 mulg5_sub_assign_w2s6,mulg5_sub_assign_w2s7};

