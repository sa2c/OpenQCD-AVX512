
/*******************************************************************************
*
* File Pbnd_dble.c
*
* Copyright (C) 2005, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic programs for the projector theta to the exterior boundary of a 
* block of lattice points (version for double-precision fields)
*
* The following are arrays of functions indexed by the face number
* ifc=0,..,7
*
*   void (*assign_sd2wd[8])(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
*     Applies the projector theta[ifc] to the spinor sd[imb[ix]], 
*     ix=0,..,vol-1, and assigns the result to the weyl spinor rd[ix]
*
*   void (*add_assign_wd2sd[8])(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
*     Expands the Weyl spinor sd[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and adds psi to rd[imb[ix]]
*
*   void (*sub_assign_wd2sd[8])(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
*     Expands the Weyl spinor sd[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and subtracts psi from rd[imb[ix]]
*
*   void (*mulg5_sub_assign_wd2sd[8])(int *imb,int vol,weyl_dble *sd,
*                                     spinor_dble *rd)
*     Expands the Weyl spinor sd[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and subtracts gamma5*psi from
*     rd[imb[ix]]
*
* Notes:
*
* The projector theta is described in the module sflds/Pbnd.c. The size and
* position of the faces is only implicitly defined through the parameter vol
* and the array *imb of the indices of the points on the face.
*
* None of these programs involves communications. They are general purpose
* routines that know nothing about the underlying geometry. In particular,
* they can be called locally. If SSE instructions are used, the fields must
* be aligned to a 16 byte boundary.
*
*******************************************************************************/

#define PBND_DBLE_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "su3.h"
#include "sflds.h"

#if (defined x64)
#include "sse2.h"

static const sse_double poh={0.5,0.5};


static void assign_sd2wd0(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   weyl_dble *rm;
   spinor_dble *si,*sin;

   rm=rd+vol;   
   si=sd+(*imb);
   imb+=(rd<(rm-1));
   sin=sd+(*imb);

   for (;rd<rm;rd++)
   {
      _sse_load_dble((*si).c1);
      _sse_load_up_dble((*si).c3);            
   
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c1);

      _sse_load_dble((*si).c2);
      _sse_load_up_dble((*si).c4);            

      si=sin;
      imb+=(rd<(rm-2));
      sin=sd+(*imb);
      _prefetch_spinor_dble(sin);      

      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c2);            
   }
}


static void assign_sd2wd1(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   weyl_dble *rm;
   spinor_dble *si,*sin;

   rm=rd+vol;   
   si=sd+(*imb);
   imb+=(rd<(rm-1));
   sin=sd+(*imb);

   for (;rd<rm;rd++)
   {   
      _sse_load_dble((*si).c1);
      _sse_load_up_dble((*si).c3);            

      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c1);

      _sse_load_dble((*si).c2);
      _sse_load_up_dble((*si).c4);            

      si=sin;
      imb+=(rd<(rm-2));
      sin=sd+(*imb);
      _prefetch_spinor_dble(sin);      
      
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c2);
   }
}


static void assign_sd2wd2(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   weyl_dble *rm;
   spinor_dble *si,*sin;

   rm=rd+vol;   
   si=sd+(*imb);
   imb+=(rd<(rm-1));
   sin=sd+(*imb);

   for (;rd<rm;rd++)
   {   
      _sse_load_dble((*si).c1);
      _sse_load_up_dble((*si).c4);

      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c1);

      _sse_load_dble((*si).c2);
      _sse_load_up_dble((*si).c3);

      si=sin;
      imb+=(rd<(rm-2));
      sin=sd+(*imb);
      _prefetch_spinor_dble(sin);      
      
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c2); 
   }
}


static void assign_sd2wd3(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   weyl_dble *rm;
   spinor_dble *si,*sin;

   rm=rd+vol;   
   si=sd+(*imb);
   imb+=(rd<(rm-1));
   sin=sd+(*imb);

   for (;rd<rm;rd++)
   {   
      _sse_load_dble((*si).c1);
      _sse_load_up_dble((*si).c4);

      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c1);

      _sse_load_dble((*si).c2);
      _sse_load_up_dble((*si).c3);

      si=sin;
      imb+=(rd<(rm-2));
      sin=sd+(*imb);
      _prefetch_spinor_dble(sin);      
      
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c2);
   }
}


static void assign_sd2wd4(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   weyl_dble *rm;
   spinor_dble *si,*sin;

   rm=rd+vol;   
   si=sd+(*imb);
   imb+=(rd<(rm-1));
   sin=sd+(*imb);

   for (;rd<rm;rd++)
   {   
      _sse_load_dble((*si).c1);
      _sse_load_up_dble((*si).c4);            

      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c1);

      _sse_load_dble((*si).c2);
      _sse_load_up_dble((*si).c3);            

      si=sin;
      imb+=(rd<(rm-2));
      sin=sd+(*imb);
      _prefetch_spinor_dble(sin);      
      
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c2);
   }
}


static void assign_sd2wd5(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   weyl_dble *rm;
   spinor_dble *si,*sin;

   rm=rd+vol;   
   si=sd+(*imb);
   imb+=(rd<(rm-1));
   sin=sd+(*imb);

   for (;rd<rm;rd++)
   {   
      _sse_load_dble((*si).c1);
      _sse_load_up_dble((*si).c4);            

      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c1);

      _sse_load_dble((*si).c2);
      _sse_load_up_dble((*si).c3);            

      si=sin;
      imb+=(rd<(rm-2));
      sin=sd+(*imb);
      _prefetch_spinor_dble(sin);      
      
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c2);
   }
}


static void assign_sd2wd6(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   weyl_dble *rm;
   spinor_dble *si,*sin;

   rm=rd+vol;   
   si=sd+(*imb);
   imb+=(rd<(rm-1));
   sin=sd+(*imb);

   for (;rd<rm;rd++)
   {   
      _sse_load_dble((*si).c1);
      _sse_load_up_dble((*si).c3);

      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c1);

      _sse_load_dble((*si).c2);
      _sse_load_up_dble((*si).c4);

      si=sin;
      imb+=(rd<(rm-2));
      sin=sd+(*imb);
      _prefetch_spinor_dble(sin);      
      
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c2);
   }
}


static void assign_sd2wd7(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   weyl_dble *rm;
   spinor_dble *si,*sin;

   rm=rd+vol;   
   si=sd+(*imb);
   imb+=(rd<(rm-1));
   sin=sd+(*imb);

   for (;rd<rm;rd++)
   {   
      _sse_load_dble((*si).c1);
      _sse_load_up_dble((*si).c3);

      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c1);

      _sse_load_dble((*si).c2);
      _sse_load_up_dble((*si).c4);

      si=sin;
      imb+=(rd<(rm-2));
      sin=sd+(*imb);
      _prefetch_spinor_dble(sin);      
      
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*rd).c2);
   }
}


static void add_assign_wd2sd0(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);            
   }
}


static void add_assign_wd2sd1(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);
   }
}


static void add_assign_wd2sd2(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {         
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);
   }
}


static void add_assign_wd2sd3(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3); 
   }
}


static void add_assign_wd2sd4(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);
   }
}


static void add_assign_wd2sd5(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {         
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3);
   }
}


static void add_assign_wd2sd6(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);
   }
}


static void add_assign_wd2sd7(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);
   }
}


static void sub_assign_wd2sd0(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {   
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      

      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);            
      
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;

      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);            
   }
}


static void sub_assign_wd2sd1(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);   
   }
}


static void sub_assign_wd2sd2(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3); 
   }
}


static void sub_assign_wd2sd3(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);
   }
}


static void sub_assign_wd2sd4(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3); 
   }
}


static void sub_assign_wd2sd5(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);
   }
}


static void sub_assign_wd2sd6(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);
   }
}


static void sub_assign_wd2sd7(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);
   }
}


static void mulg5_sub_assign_wd2sd0(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {      
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);            
   }
}


static void mulg5_sub_assign_wd2sd1(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {         
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4); 
   }
}


static void mulg5_sub_assign_wd2sd2(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {         
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3); 
   }
}


static void mulg5_sub_assign_wd2sd3(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {         
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3);
   }
}


static void mulg5_sub_assign_wd2sd4(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {         
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);  
   }
}


static void mulg5_sub_assign_wd2sd5(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {         
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3);
   }
}


static void mulg5_sub_assign_wd2sd6(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {         
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c4);
   }
}


static void mulg5_sub_assign_wd2sd7(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri,*rin,*rim;

   sm=sd+vol;
   rin=rd+(*imb);
   imb+=(sd<(sm-1));
   rim=rd+(*imb);

   for (;sd<sm;)
   {         
      _sse_load_up_dble((*sd).c1);
      _sse_load_dble((*rin).c1);            

      ri=rin;
      rin=rim;
      imb+=(sd<(sm-2));
      rim=rd+(*imb);
      _prefetch_spinor_dble(rim);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c1);
      _sse_load_dble((*ri).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c3);            
            
      _sse_load_up_dble((*sd).c2);
      _sse_load_dble((*ri).c2);            

      sd+=4;
      _prefetch_weyl_dble(sd);
      sd-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*ri).c2);
      _sse_load_dble((*ri).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*ri).c4);
   }
}

#else

static void assign_sd2wd0(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   double r1;
   weyl_dble *rm;
   spinor_dble *si;

   r1=0.5;
   rm=rd+vol;

   for (;rd<rm;rd++)
   {
      si=sd+(*imb);
      imb+=1;
      _vector_sub((*rd).c1,(*si).c1,(*si).c3);
      _vector_sub((*rd).c2,(*si).c2,(*si).c4);       
      _vector_mul((*rd).c1,r1,(*rd).c1);
      _vector_mul((*rd).c2,r1,(*rd).c2);      
   }
}


static void assign_sd2wd1(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   double r1;
   weyl_dble *rm;
   spinor_dble *si;

   r1=0.5;
   rm=rd+vol;

   for (;rd<rm;rd++)
   {
      si=sd+(*imb);
      imb+=1;
      _vector_add((*rd).c1,(*si).c1,(*si).c3);
      _vector_add((*rd).c2,(*si).c2,(*si).c4);       
      _vector_mul((*rd).c1,r1,(*rd).c1);
      _vector_mul((*rd).c2,r1,(*rd).c2); 
   }
}


static void assign_sd2wd2(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   double r1;
   weyl_dble *rm;
   spinor_dble *si;

   r1=0.5;
   rm=rd+vol;

   for (;rd<rm;rd++)
   {
      si=sd+(*imb);
      imb+=1;
      _vector_i_sub((*rd).c1,(*si).c1,(*si).c4);
      _vector_i_sub((*rd).c2,(*si).c2,(*si).c3);       
      _vector_mul((*rd).c1,r1,(*rd).c1);
      _vector_mul((*rd).c2,r1,(*rd).c2);     
   }
}


static void assign_sd2wd3(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   double r1;
   weyl_dble *rm;
   spinor_dble *si;

   r1=0.5;
   rm=rd+vol;

   for (;rd<rm;rd++)
   {
      si=sd+(*imb);
      imb+=1;
      _vector_i_add((*rd).c1,(*si).c1,(*si).c4);
      _vector_i_add((*rd).c2,(*si).c2,(*si).c3);       
      _vector_mul((*rd).c1,r1,(*rd).c1);
      _vector_mul((*rd).c2,r1,(*rd).c2);    
   }
}


static void assign_sd2wd4(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   double r1;
   weyl_dble *rm;
   spinor_dble *si;

   r1=0.5;
   rm=rd+vol;

   for (;rd<rm;rd++)
   {
      si=sd+(*imb);
      imb+=1;
      _vector_sub((*rd).c1,(*si).c1,(*si).c4);
      _vector_add((*rd).c2,(*si).c2,(*si).c3);       
      _vector_mul((*rd).c1,r1,(*rd).c1);
      _vector_mul((*rd).c2,r1,(*rd).c2);
   }
}


static void assign_sd2wd5(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   double r1;
   weyl_dble *rm;
   spinor_dble *si;

   r1=0.5;
   rm=rd+vol;

   for (;rd<rm;rd++)
   {
      si=sd+(*imb);
      imb+=1;
      _vector_add((*rd).c1,(*si).c1,(*si).c4);
      _vector_sub((*rd).c2,(*si).c2,(*si).c3);       
      _vector_mul((*rd).c1,r1,(*rd).c1);
      _vector_mul((*rd).c2,r1,(*rd).c2); 
   }
}


static void assign_sd2wd6(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   double r1;
   weyl_dble *rm;
   spinor_dble *si;

   r1=0.5;
   rm=rd+vol;

   for (;rd<rm;rd++)
   {
      si=sd+(*imb);
      imb+=1;
      _vector_i_sub((*rd).c1,(*si).c1,(*si).c3);
      _vector_i_add((*rd).c2,(*si).c2,(*si).c4);       
      _vector_mul((*rd).c1,r1,(*rd).c1);
      _vector_mul((*rd).c2,r1,(*rd).c2);  
   }
}


static void assign_sd2wd7(int *imb,int vol,spinor_dble *sd,weyl_dble *rd)
{
   double r1;
   weyl_dble *rm;
   spinor_dble *si;

   r1=0.5;
   rm=rd+vol;

   for (;rd<rm;rd++)
   {
      si=sd+(*imb);
      imb+=1;
      _vector_i_add((*rd).c1,(*si).c1,(*si).c3);
      _vector_i_sub((*rd).c2,(*si).c2,(*si).c4);       
      _vector_mul((*rd).c1,r1,(*rd).c1);
      _vector_mul((*rd).c2,r1,(*rd).c2);  
   }
}


static void add_assign_wd2sd0(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_add_assign((*ri).c1,(*sd).c1);
      _vector_add_assign((*ri).c2,(*sd).c2);
      _vector_sub_assign((*ri).c3,(*sd).c1);
      _vector_sub_assign((*ri).c4,(*sd).c2);
   }
}


static void add_assign_wd2sd1(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_add_assign((*ri).c1,(*sd).c1);
      _vector_add_assign((*ri).c2,(*sd).c2);
      _vector_add_assign((*ri).c3,(*sd).c1);
      _vector_add_assign((*ri).c4,(*sd).c2); 
   }
}

   
static void add_assign_wd2sd2(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_add_assign((*ri).c1,(*sd).c1);
      _vector_add_assign((*ri).c2,(*sd).c2);
      _vector_i_add_assign((*ri).c3,(*sd).c2);
      _vector_i_add_assign((*ri).c4,(*sd).c1);  
   }
}

   
static void add_assign_wd2sd3(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_add_assign((*ri).c1,(*sd).c1);
      _vector_add_assign((*ri).c2,(*sd).c2);
      _vector_i_sub_assign((*ri).c3,(*sd).c2);
      _vector_i_sub_assign((*ri).c4,(*sd).c1);
   }
}

   
static void add_assign_wd2sd4(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_add_assign((*ri).c1,(*sd).c1);
      _vector_add_assign((*ri).c2,(*sd).c2);
      _vector_add_assign((*ri).c3,(*sd).c2);
      _vector_sub_assign((*ri).c4,(*sd).c1);
   }
}

   
static void add_assign_wd2sd5(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_add_assign((*ri).c1,(*sd).c1);
      _vector_add_assign((*ri).c2,(*sd).c2);
      _vector_sub_assign((*ri).c3,(*sd).c2);
      _vector_add_assign((*ri).c4,(*sd).c1);
   }
}

   
static void add_assign_wd2sd6(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_add_assign((*ri).c1,(*sd).c1);
      _vector_add_assign((*ri).c2,(*sd).c2);
      _vector_i_add_assign((*ri).c3,(*sd).c1);
      _vector_i_sub_assign((*ri).c4,(*sd).c2);
   }
}

   
static void add_assign_wd2sd7(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_add_assign((*ri).c1,(*sd).c1);
      _vector_add_assign((*ri).c2,(*sd).c2);
      _vector_i_sub_assign((*ri).c3,(*sd).c1);
      _vector_i_add_assign((*ri).c4,(*sd).c2);
   }
}


static void sub_assign_wd2sd0(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_add_assign((*ri).c3,(*sd).c1);
      _vector_add_assign((*ri).c4,(*sd).c2);
   }
}


static void sub_assign_wd2sd1(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_sub_assign((*ri).c3,(*sd).c1);
      _vector_sub_assign((*ri).c4,(*sd).c2); 
   }
}


static void sub_assign_wd2sd2(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_i_sub_assign((*ri).c3,(*sd).c2);
      _vector_i_sub_assign((*ri).c4,(*sd).c1);  
   }
}


static void sub_assign_wd2sd3(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_i_add_assign((*ri).c3,(*sd).c2);
      _vector_i_add_assign((*ri).c4,(*sd).c1);
   }
}


static void sub_assign_wd2sd4(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_sub_assign((*ri).c3,(*sd).c2);
      _vector_add_assign((*ri).c4,(*sd).c1);
   }
}


static void sub_assign_wd2sd5(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_add_assign((*ri).c3,(*sd).c2);
      _vector_sub_assign((*ri).c4,(*sd).c1);
   }
}


static void sub_assign_wd2sd6(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_i_sub_assign((*ri).c3,(*sd).c1);
      _vector_i_add_assign((*ri).c4,(*sd).c2);
   }
}


static void sub_assign_wd2sd7(int *imb,int vol,weyl_dble *sd,spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_i_add_assign((*ri).c3,(*sd).c1);
      _vector_i_sub_assign((*ri).c4,(*sd).c2);
   }
}


static void mulg5_sub_assign_wd2sd0(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_sub_assign((*ri).c3,(*sd).c1);
      _vector_sub_assign((*ri).c4,(*sd).c2);
   }
}


static void mulg5_sub_assign_wd2sd1(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_add_assign((*ri).c3,(*sd).c1);
      _vector_add_assign((*ri).c4,(*sd).c2); 
   }
}


static void mulg5_sub_assign_wd2sd2(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_i_add_assign((*ri).c3,(*sd).c2);
      _vector_i_add_assign((*ri).c4,(*sd).c1); 
   }
}


static void mulg5_sub_assign_wd2sd3(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_i_sub_assign((*ri).c3,(*sd).c2);
      _vector_i_sub_assign((*ri).c4,(*sd).c1);
   }
}


static void mulg5_sub_assign_wd2sd4(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_add_assign((*ri).c3,(*sd).c2);
      _vector_sub_assign((*ri).c4,(*sd).c1);
   }
}


static void mulg5_sub_assign_wd2sd5(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_sub_assign((*ri).c3,(*sd).c2);
      _vector_add_assign((*ri).c4,(*sd).c1);
   }
}


static void mulg5_sub_assign_wd2sd6(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_i_add_assign((*ri).c3,(*sd).c1);
      _vector_i_sub_assign((*ri).c4,(*sd).c2);
   }
}


static void mulg5_sub_assign_wd2sd7(int *imb,int vol,weyl_dble *sd,
                                    spinor_dble *rd)
{
   weyl_dble *sm;
   spinor_dble *ri;

   sm=sd+vol;
   
   for (;sd<sm;sd++)
   {
      ri=rd+(*imb);
      imb+=1;
      _vector_sub_assign((*ri).c1,(*sd).c1);
      _vector_sub_assign((*ri).c2,(*sd).c2);
      _vector_i_sub_assign((*ri).c3,(*sd).c1);
      _vector_i_add_assign((*ri).c4,(*sd).c2);
   }
}

#endif

void (*assign_sd2wd[8])(int *imb,int vol,spinor_dble *sd,weyl_dble *rd) =
{assign_sd2wd0,assign_sd2wd1,assign_sd2wd2,assign_sd2wd3,
 assign_sd2wd4,assign_sd2wd5,assign_sd2wd6,assign_sd2wd7};

void (*add_assign_wd2sd[8])(int *imb,int vol,weyl_dble *sd,spinor_dble *rd) =
{add_assign_wd2sd0,add_assign_wd2sd1,add_assign_wd2sd2,add_assign_wd2sd3,
 add_assign_wd2sd4,add_assign_wd2sd5,add_assign_wd2sd6,add_assign_wd2sd7};

void (*sub_assign_wd2sd[8])(int *imb,int vol,weyl_dble *sd,spinor_dble *rd) =
{sub_assign_wd2sd0,sub_assign_wd2sd1,sub_assign_wd2sd2,sub_assign_wd2sd3,
 sub_assign_wd2sd4,sub_assign_wd2sd5,sub_assign_wd2sd6,sub_assign_wd2sd7};

void (*mulg5_sub_assign_wd2sd[8])(int *imb,int vol,weyl_dble *sd,
                                  spinor_dble *rd) =
{mulg5_sub_assign_wd2sd0,mulg5_sub_assign_wd2sd1,
 mulg5_sub_assign_wd2sd2,mulg5_sub_assign_wd2sd3,
 mulg5_sub_assign_wd2sd4,mulg5_sub_assign_wd2sd5,
 mulg5_sub_assign_wd2sd6,mulg5_sub_assign_wd2sd7};
