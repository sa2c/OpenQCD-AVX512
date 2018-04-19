
/*******************************************************************************
*
* File cgne.c
*
* Copyright (C) 2005, 2008, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic CG solver program for the lattice Dirac equation
*
* The externally accessible function is
*
*   double cgne(int vol,int icom,void (*Dop)(spinor *s,spinor *r),
*               void (*Dop_dble)(spinor_dble *s,spinor_dble *r),
*               spinor **ws,spinor_dble **wsd,int nmx,double res,
*               spinor_dble *eta,spinor_dble *psi,int *status)
*     Solution of the (normal) Dirac equation D^dag*D*psi=eta for given
*     source eta, using the CG algorithm. See the notes for the explanation
*     of the parameters of the program.

* Notes:
*
* This program uses single-precision arithmetic to reduce the execution
* time but obtains the solution with double-precision accuracy.
*
* The programs for the single- and double-precision implementations of the
* Dirac operator are assumed to have the following properties:
*
*   void Dop(spinor *s,spinor *r)
*     Application the operator D or its hermitian conjugate D^dag to the 
*     single-precision Dirac field s and assignement of the result to r.
*     D and D^dag are applied alternatingly, i.e. the first call of the
*     program applies D, the next call D^dag, then D again and so on. In
*     all cases the source field s is unchanged.
*
*   void Dop_dble(spinor *s,spinor *r)
*     Application the operator D or its hermitian conjugate D^dag to the 
*     double-precision Dirac field s and assignement of the result to r.
*     D and D^dag are applied alternatively, i.e. the first call of the
*     program applies D, the next call D^dag, then D again and so on. In
*     all cases the source field s is unchanged.
*
* The other parameters of the program cgne() are:
*
*   vol     Number of spinors in the Dirac fields.         
*
*   icom    Indicates whether the equation to be solved is a local
*           equation (icom=0) or a global one (icom=1). Scalar products
*           are summed over all MPI processes if icom=1, while no
*           communications are performed if icom=0.
*
*   nmx     Maximal total number of CG iterations that may be applied.
*
*   res     Desired maximal relative residue |eta-D^dag*D*psi|/|eta| of
*           the calculated solution.
*
*   ws      Array of at least 5 single-precision spinor fields (used
*           as work space).
*
*   wsd     Array of at least 2 double-precision spinor fields (used
*           as work space).
*
*   eta     Source field (unchanged on exit).
*
*   psi     Calculated approximate solution of the Dirac equation
*           D^dag*D*psi=eta.
*
*   status  On exit, this parameter reports the total number of CG
*           iterations that were required, or a negative value if the
*           program failed.
*
* Independently of whether the program succeeds in solving the Dirac equation
* to the desired accuracy, the program returns the norm of the residue of
* the field psi.
*
* Some debugging output is printed to stdout on process 0 if CGNE_DBG is
* defined at compilation time.
*
*******************************************************************************/

#define CGNE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "sflds.h"
#include "linalg.h"
#include "linsolv.h"
#include "global.h"

#define PRECISION_LIMIT ((double)(100.0f*FLT_EPSILON))

static float rsq,rsq_old,ai,bi;
static spinor *psx,*psr,*psp,*psap,*psw;
static spinor_dble *pdb,*pdx,*pdw,*pdv;

#if (defined x64)
#include "sse2.h"

static void update_g(int vol)
{
   float c;
   spinor *r,*s,*sm;

   c=-ai;
   
   __asm__ __volatile__ ("movss %0, %%xmm6 \n\t"
                         "shufps $0x0, %%xmm6, %%xmm6 \n\t"
                         "movaps %%xmm6, %%xmm7 \n\t"
                         "movaps %%xmm6, %%xmm8"
                         :
                         :
                         "m" (c)
                         :
                         "xmm6", "xmm7", "xmm8");

   r=psr;
   s=psap;
   sm=s+vol;

   for (;s<sm;)
   {
      _sse_spinor_load(*s);

      s+=4;
      _prefetch_spinor(s);
      s-=3;      
      
      __asm__ __volatile__ ("mulps %%xmm6, %%xmm0 \n\t"
                            "mulps %%xmm7, %%xmm1 \n\t"
                            "mulps %%xmm8, %%xmm2 \n\t"
                            "mulps %%xmm6, %%xmm3 \n\t"
                            "mulps %%xmm7, %%xmm4 \n\t"
                            "mulps %%xmm8, %%xmm5 \n\t"
                            "addps %0, %%xmm0 \n\t"
                            "addps %2, %%xmm1 \n\t"
                            "addps %4, %%xmm2"                            
                            :
                            :
                            "m" ((*r).c1.c1),
                            "m" ((*r).c1.c2),
                            "m" ((*r).c1.c3),
                            "m" ((*r).c2.c1),
                            "m" ((*r).c2.c2),
                            "m" ((*r).c2.c3)                            
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      r+=4;
      _prefetch_spinor(r);
      r-=4;
      
      __asm__ __volatile__ ("addps %0, %%xmm3 \n\t"
                            "addps %2, %%xmm4 \n\t"
                            "addps %4, %%xmm5"
                            :
                            :
                            "m" ((*r).c3.c1),
                            "m" ((*r).c3.c2),
                            "m" ((*r).c3.c3),
                            "m" ((*r).c4.c1),
                            "m" ((*r).c4.c2),
                            "m" ((*r).c4.c3)
                            :
                            "xmm3", "xmm4", "xmm5");
      
      _sse_spinor_store(*r);

      r+=1;
   }
}


static void update_xp(int vol)
{
   spinor *r,*s,*t,*tm;

   __asm__ __volatile__ ("movss %0, %%xmm6 \n\t"
                         "movss %1, %%xmm9 \n\t"
                         "shufps $0x0, %%xmm6, %%xmm6 \n\t"
                         "shufps $0x0, %%xmm9, %%xmm9 \n\t"
                         "movaps %%xmm6, %%xmm7 \n\t"
                         "movaps %%xmm9, %%xmm10 \n\t"
                         "movaps %%xmm6, %%xmm8 \n\t"                         
                         "movaps %%xmm9, %%xmm11"                         
                         :
                         :
                         "m" (ai),
                         "m" (bi)
                         :
                         "xmm6", "xmm7", "xmm8", "xmm9",
                         "xmm10", "xmm11");
   
   r=psr;
   s=psp;
   t=psx;
   tm=t+vol;

   for (;t<tm;t++)
   {
      _sse_spinor_load(*s);

      s+=4;
      _prefetch_spinor(s);
      s-=4;  
      
      __asm__ __volatile__ ("mulps %%xmm6, %%xmm0 \n\t"
                            "mulps %%xmm7, %%xmm1 \n\t"
                            "mulps %%xmm8, %%xmm2 \n\t"
                            "mulps %%xmm6, %%xmm3 \n\t"
                            "mulps %%xmm7, %%xmm4 \n\t"
                            "mulps %%xmm8, %%xmm5 \n\t"
                            "addps %0, %%xmm0 \n\t"
                            "addps %2, %%xmm1 \n\t"
                            "addps %4, %%xmm2"
                            :
                            :
                            "m" ((*t).c1.c1),
                            "m" ((*t).c1.c2),
                            "m" ((*t).c1.c3),
                            "m" ((*t).c2.c1),
                            "m" ((*t).c2.c2),
                            "m" ((*t).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      t+=4;
      _prefetch_spinor(t);
      t-=4;

      __asm__ __volatile__ ("addps %0, %%xmm3 \n\t"
                            "addps %2, %%xmm4 \n\t"
                            "addps %4, %%xmm5"
                            :
                            :
                            "m" ((*t).c3.c1),
                            "m" ((*t).c3.c2),
                            "m" ((*t).c3.c3),
                            "m" ((*t).c4.c1),
                            "m" ((*t).c4.c2),
                            "m" ((*t).c4.c3)
                            :
                            "xmm3", "xmm4", "xmm5");

      _sse_spinor_store(*t); 
      _sse_spinor_load(*s);

      r+=4;
      _prefetch_spinor(r);
      r-=4;
      
      __asm__ __volatile__ ("mulps %%xmm9, %%xmm0 \n\t"
                            "mulps %%xmm10, %%xmm1 \n\t"
                            "mulps %%xmm11, %%xmm2 \n\t"
                            "mulps %%xmm9, %%xmm3 \n\t"
                            "mulps %%xmm10, %%xmm4 \n\t"
                            "mulps %%xmm11, %%xmm5 \n\t"
                            "addps %0, %%xmm0 \n\t"
                            "addps %2, %%xmm1 \n\t"
                            "addps %4, %%xmm2"
                            :
                            :
                            "m" ((*r).c1.c1),
                            "m" ((*r).c1.c2),
                            "m" ((*r).c1.c3),
                            "m" ((*r).c2.c1),
                            "m" ((*r).c2.c2),
                            "m" ((*r).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");
      
      __asm__ __volatile__ ("addps %0, %%xmm3 \n\t"
                            "addps %2, %%xmm4 \n\t"
                            "addps %4, %%xmm5"
                            :
                            :
                            "m" ((*r).c3.c1),
                            "m" ((*r).c3.c2),
                            "m" ((*r).c3.c3),
                            "m" ((*r).c4.c1),
                            "m" ((*r).c4.c2),
                            "m" ((*r).c4.c3)
                            :
                            "xmm3", "xmm4", "xmm5");
      
      _sse_spinor_store(*s);

      r+=1;
      s+=1;
   }
}

#else

static void update_g(int vol)
{
   float c;
   spinor *r,*s,*sm;

   c=-ai;
   r=psr;
   s=psap;
   sm=s+vol;

   for (;s<sm;s++)
   {
      _vector_mulr_assign((*r).c1,c,(*s).c1);
      _vector_mulr_assign((*r).c2,c,(*s).c2);      
      _vector_mulr_assign((*r).c3,c,(*s).c3);
      _vector_mulr_assign((*r).c4,c,(*s).c4);      

      r+=1;
   }
}


static void update_xp(int vol)
{
   spinor *r,*s,*t,*tm;

   r=psr;
   s=psp;
   t=psx;
   tm=t+vol;

   for (;t<tm;t++)
   {
      _vector_mulr_assign((*t).c1,ai,(*s).c1);
      _vector_mulr_assign((*t).c2,ai,(*s).c2);
      _vector_mulr_assign((*t).c3,ai,(*s).c3);
      _vector_mulr_assign((*t).c4,ai,(*s).c4);

      _vector_mulr_add((*s).c1,bi,(*r).c1);
      _vector_mulr_add((*s).c2,bi,(*r).c2);
      _vector_mulr_add((*s).c3,bi,(*r).c3);
      _vector_mulr_add((*s).c4,bi,(*r).c4);      

      r+=1;
      s+=1;
   }
}

#endif

static void cg_init(int vol,int icom,spinor **ws,spinor_dble **wsd,
                    spinor_dble *eta,spinor_dble *psi)
{
   psx=ws[0];
   psr=ws[1];
   psp=ws[2];
   psap=ws[3];
   psw=ws[4];

   pdb=eta;
   pdx=psi;
   pdw=wsd[0];
   pdv=wsd[1];

   set_s2zero(vol,psx);
   assign_sd2s(vol,pdb,psr);
   assign_s2s(vol,psr,psp);
   set_sd2zero(vol,pdx);

   rsq=norm_square(vol,icom,psr);
}


static void cg_step(int vol,int icom,void (*Dop)(spinor *s,spinor *r))
{
   (*Dop)(psp,psw);
   (*Dop)(psw,psap);

   ai=rsq/norm_square(vol,icom,psw);
   update_g(vol);

   rsq_old=rsq;
   rsq=norm_square(vol,icom,psr);
   bi=rsq/rsq_old;
   update_xp(vol);
}


static void cg_reset(int vol,int icom,void (*Dop)(spinor *s,spinor *r),
                     void (*Dop_dble)(spinor_dble *s,spinor_dble *r))
{
   float r;
   complex z;

   (*Dop_dble)(pdx,pdw);
   (*Dop_dble)(pdw,pdv);

   diff_sd2s(vol,pdb,pdv,psr);
   rsq=norm_square(vol,icom,psr);

   assign_s2s(vol,psp,psw);
   assign_s2s(vol,psr,psp);

   z=spinor_prod(vol,icom,psr,psw);
   z.re=-z.re/rsq;
   z.im=-z.im/rsq;
   mulc_spinor_add(vol,psw,psr,z);

   (*Dop)(psw,psx);
   (*Dop)(psx,psap);

   r=norm_square(vol,icom,psx);
   z=spinor_prod(vol,icom,psap,psr);

   if ((z.re*z.re+z.im*z.im)<(2.0f*r*r))
   {
      z.re=-z.re/r;
      z.im=-z.im/r;
      mulc_spinor_add(vol,psp,psw,z);
   }

   set_s2zero(vol,psx);
}


double cgne(int vol,int icom,void (*Dop)(spinor *s,spinor *r),
            void (*Dop_dble)(spinor_dble *s,spinor_dble *r),
            spinor **ws,spinor_dble **wsd,int nmx,double res,
            spinor_dble *eta,spinor_dble *psi,int *status)
{
   int ncg,iprms[2];
   double xn,rn,tol,dprms[1];

   if ((icom==1)&&(NPROC>1))
   {
      iprms[0]=vol;
      iprms[1]=nmx;
      dprms[0]=res;

      MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      error((iprms[0]!=vol)||(iprms[1]!=nmx)||(dprms[0]!=res),1,
            "cgne [cgne.c]","Parameters are not global");

      error_root((vol<=0)||(nmx<1)||(res<=DBL_EPSILON),1,
                 "cgne [cgne.c]","Improper choice of vol,nmx or res");
   }
   else
   {
      if ((vol<=0)||(nmx<1)||(res<=DBL_EPSILON))
      {
         error_loc(1,1,"cgne [cgne.c]",
                   "Improper choice of vol,nmx or res");
         (*status)=0;
         return 1.0;
      }
   }

   cg_init(vol,icom,ws,wsd,eta,psi);
   rn=sqrt((double)(rsq));
   tol=res*rn;
   (*status)=0;

   xn=(double)(norm_square(vol,icom,psx));
   xn=sqrt(xn);

   while (rn>tol)
   {
#ifdef CGNE_DBG
      message("[cgne]: rn_old = %.2e\n",rn);
#endif
      ncg=0;

      for (;;)
      {
         cg_step(vol,icom,Dop);
         ncg+=1;
         (*status)+=1;

         xn=(double)(norm_square(vol,icom,psx));
         xn=sqrt(xn);
         rn=sqrt((double)(rsq));
#ifdef CGNE_DBG
         message("[cgne]: ncg = %d, xn = %.2e, rn = %.2e\n",(*status),xn,rn);
#endif         
         if ((rn<=tol)||(rn<=(PRECISION_LIMIT*xn))||(ncg>=100)||
             ((*status)>=nmx))
            break;
      }

      add_s2sd(vol,psx,pdx);
      xn=norm_square_dble(vol,icom,pdx);
      xn=sqrt(xn);
      cg_reset(vol,icom,Dop,Dop_dble);
      rn=sqrt((double)(rsq));

      if (((*status)>=nmx)&&(rn>tol))
      {
         (*status)=-1;
         break;
      }

      if ((100.0*DBL_EPSILON*xn)>tol)
      {
         (*status)=-2;
         break;
      }
   }

   return rn;
}
