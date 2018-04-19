
/*******************************************************************************
*
* File mdsteps.c
*
* Copyright (C) 2011, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Molecular-dynamics integrator
*
* The externally accessible functions are
* 
*   void set_mdsteps(void)
*     Constructs the integrator from the data available in the parameter
*     data base (see the notes). The integrator is stored internally in
*     the form of an array of elementary operations (force computations
*     and gauge-field update steps).
*
*   mdstep_t *mdsteps(int *nop,int *ntu,int *itu)
*     Returns the array of elementary operations that describe the current
*     integrator. On exit the program assigns the total number of operations
*     to nop and the index of the gauge-field update operation to itu.
*
*   void print_mdsteps(int ipr)
*     Prints some information on the current integrator to stdout on MPI
*     process 0. The program always prints the available information on
*     the different levels of the integrator. Whether further information
*     is printed depends on the 3 low bits of the print flat ipr:
*
*      if (ipr&0x1): Force descriptions
*
*      if (ipr&0x2): List of elementary operations
*
*      if (ipr&0x4): Integration time check
*
*     The full information is thus printed if ipr=0x7.
*
* Notes:
*
* The structure of the MD integrator is explained in the file README.mdint
* in this directory. It is assumed here that the parameters of the integrator
* have been entered to the parameter data base.
*
* An elementary update step is described by a structure of type mdstep_t
* with the following elements:
*
*  iop     Operation index (0<=iop<=itu+1). If iop<itu, the force number
*          iop is to be computed and to be assigned (gauge force) or added
*          (fermion forces) to the force field. If iop=itu, the momentum
*          and subsequently the gauge field are to be updated, using the
*          current force field. If iop=itu+1, the momentum field is to be
*          updated, using the current force, and the integration ends.
*
*  eps     Step sizes by which the forces (iop<itu) or the momentum field
*          in the update of the gauge field must be multiplied.
*
* The forces are described by the structures returned by force_parms(iop)
* if iop<itu (see flags/force_parms.c).
*
* In the operation list constructed by set_mdsteps(), the forces in each
* period from the last update of the gauge field to the next are ordered
* such that the gauge force comes first. The fermion forces are ordered
* according to their index.
*
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "update.h"
#include "global.h"

static int nsmx,nmds=0,iend=1;
static mdstep_t *mds=NULL,*mdw[3];


static void set_nsmx(int nlv)
{
   int ntu,ilv;
   int nfr,*ifr,i;
   mdint_parms_t mdp;

   iend=0;
   ntu=1;
   
   for (ilv=0;ilv<nlv;ilv++)
   {
      mdp=mdint_parms(ilv);

      if (mdp.integrator==LPFR)
         ntu*=mdp.nstep;
      else if (mdp.integrator==OMF2)
         ntu*=2*mdp.nstep;
      else if (mdp.integrator==OMF4)
         ntu*=5*mdp.nstep;
      else
         error_root(1,1,"set_nsmx [mdsteps.c]","Unknown integrator");

      nfr=mdp.nfr;
      ifr=mdp.ifr;
      
      for (i=0;i<nfr;i++)
      {
         if (ifr[i]>iend)
            iend=ifr[i];
      }
   }

   iend+=2;   
   nsmx=(ntu+1)*iend;
}


static void alloc_mds(void)
{
   int k;
   
   if (mds!=NULL)
      free(mds);

   mds=malloc(4*nsmx*sizeof(*mds));
   error(mds==NULL,1,"alloc_mds [mdsteps.c]",
         "Unable to allocate mdsteps array");

   for (k=0;k<3;k++)
      mdw[k]=mds+(k+1)*nsmx;
}


static void set_steps2zero(int n,mdstep_t *s)
{
   int i;

   for (i=0;i<n;i++)
   {
      s[i].iop=iend;
      s[i].eps=0.0;
   }
}


static void copy_steps(int n,double c,mdstep_t *s,mdstep_t *r)
{
   int i;

   for (i=0;i<n;i++)
   {
      r[i].iop=s[i].iop;
      r[i].eps=c*s[i].eps;     
   }
}


static void add_steps(int n,double c,mdstep_t *s,mdstep_t *r)
{
   int itu,i,j;

   itu=iend-1;
   
   for (i=0;i<n;i++)
   {
      for (j=0;r[j].iop<itu;j++)
      {
         if (r[j].iop==s[i].iop)
         {
            r[j].eps+=c*s[i].eps;
            break;
         }
      }

      if (r[j].iop>=itu)
      {
         r[j].iop=s[i].iop;
         r[j].eps=c*s[i].eps;
      }
   }
}


static int nfrc_steps(mdstep_t *s)
{
   int itu,n;

   itu=iend-1;
   n=0;
   
   while (s[n].iop<itu)
      n+=1;

   return n;
}


static int nall_steps(mdstep_t *s)
{
   int n;
   
   n=0;
   
   while (s[n].iop<iend)
      n+=1;

   return n;
}


static void expand_level(int ilv,double tau,mdstep_t *s,mdstep_t *ws)
{
   int nstep,nfr,*ifr;
   int itu,n,i,j;
   double r0,r1,r2,r3,r4,eps;
   mdint_parms_t mdp;

   mdp=mdint_parms(ilv);
   nstep=mdp.nstep;
   nfr=mdp.nfr;
   ifr=mdp.ifr;

   itu=iend-1;
   n=0;
   r0=mdp.lambda;
   r1=0.08398315262876693;
   r2=0.2539785108410595; 
   r3=0.6822365335719091; 
   r4=-0.03230286765269967;
   eps=tau/(double)(nstep);   

   set_steps2zero(nsmx,s);
   set_steps2zero(nsmx,ws);
   
   for (i=0;i<nfr;i++)
   {
      for (j=0;j<n;j++)
      {
         if (ifr[i]==ws[j].iop)
         {
            ws[j].eps+=eps;
            break;
         }
      }
      
      if (j==n)
      {
         ws[n].iop=ifr[i];
         ws[n].eps=eps;
         n+=1;
      }
   }

   if (mdp.integrator==LPFR)
   {
      copy_steps(n,0.5,ws,s);
      s+=n;

      for (i=1;i<=nstep;i++)
      {
         (*s).iop=itu;
         (*s).eps=eps;
         s+=1;
         if (i<nstep)
            copy_steps(n,1.0,ws,s);
         else
            copy_steps(n,0.5,ws,s);
         s+=n;
      }
   }
   else if (mdp.integrator==OMF2)
   {
      copy_steps(n,r0,ws,s);
      s+=n;

      for (i=1;i<=nstep;i++)
      {      
         (*s).iop=itu;
         (*s).eps=0.5*eps;
         s+=1;
         copy_steps(n,1.0-2.0*r0,ws,s);      
         s+=n;

         (*s).iop=itu;
         (*s).eps=0.5*eps;
         s+=1;
         if (i<nstep)
            copy_steps(n,2.0*r0,ws,s);
         else
            copy_steps(n,r0,ws,s);
         s+=n;
      }
   }
   else if (mdp.integrator==OMF4)
   {
      copy_steps(n,r1,ws,s);
      s+=n;

      for (i=1;i<=nstep;i++)
      {       
         (*s).iop=itu;
         (*s).eps=r2*eps;
         s+=1;
         copy_steps(n,r3,ws,s);
         s+=n;

         (*s).iop=itu;
         (*s).eps=r4*eps;
         s+=1;
         copy_steps(n,0.5-r1-r3,ws,s);
         s+=n;

         (*s).iop=itu;
         (*s).eps=(1.0-2.0*(r2+r4))*eps;
         s+=1;
         copy_steps(n,0.5-r1-r3,ws,s);
         s+=n;

         (*s).iop=itu;
         (*s).eps=r4*eps;
         s+=1;
         copy_steps(n,r3,ws,s);
         s+=n;

         (*s).iop=itu;
         (*s).eps=r2*eps;
         s+=1;
         if (i<nstep)         
            copy_steps(n,2.0*r1,ws,s);
         else
            copy_steps(n,r1,ws,s);
         s+=n;
      }
   }
}


static void insert_level(mdstep_t *s1,mdstep_t *s2,mdstep_t *r)
{
   int itu,n,nfrc,nall;
   double eps;

   set_steps2zero(nsmx,r);

   itu=iend-1;
   nfrc=nfrc_steps(s1);
   nall=nall_steps(s1+nfrc);

   n=nfrc_steps(s2);
   copy_steps(n,1.0,s2,r);
   s2+=n;   
                  
   while ((*s2).iop==itu)
   {
      eps=(*s2).eps;
      add_steps(nfrc,eps,s1,r);
      r+=nfrc_steps(r);
      copy_steps(nall,eps,s1+nfrc,r);
      r+=nall-nfrc;
      
      s2+=1;
      n=nfrc_steps(s2);
      add_steps(n,1.0,s2,r);
      s2+=n;
   }
}


static void swap_steps(mdstep_t *s,mdstep_t *r)
{
   int is;
   double rs;

   is=(*s).iop;
   (*s).iop=(*r).iop;
   (*r).iop=is;

   rs=(*s).eps;
   (*s).eps=(*r).eps;
   (*r).eps=rs;
}


static void sort_forces(void)
{
   int itu,n;
   int i,j,k,imn;
   mdstep_t *s;
   force_parms_t fp;

   itu=iend-1;
   s=mds;

   while ((*s).iop<iend)
   {
      n=nfrc_steps(s);
      k=0;

      for (i=0;i<n;i++)
      {
         fp=force_parms(s[i].iop);

         if (fp.force==FRG)
         {
            if (i>0)
               swap_steps(s,s+i);
            k+=1;
         }
      }

      error_root(k!=1,1,"sort_forces [mdsteps.c]",
                 "Incorrect gauge force count");

      for (i=1;i<n;i++)
      {
         imn=s[i].iop;
         k=i;
         
         for (j=(i+1);j<n;j++)
         {
            if (s[j].iop<imn)
            {
               imn=s[j].iop;
               k=j;
            }
         }

         if (k!=i)
            swap_steps(s+i,s+k);
      }
      
      s+=n;
      if ((*s).iop==itu)
         s+=1;
   }
}


void set_mdsteps(void)
{
   int nlv,ilv,n;
   double tau;
   hmc_parms_t hmc;

   hmc=hmc_parms();
   nlv=hmc.nlv;
   tau=hmc.tau;

   set_nsmx(nlv);
   alloc_mds();
   expand_level(nlv-1,tau,mds,mdw[0]);
      
   for (ilv=(nlv-2);ilv>=0;ilv--)
   {
      n=nall_steps(mds);
      copy_steps(n,1.0,mds,mdw[0]);      
      expand_level(ilv,1.0,mdw[1],mdw[2]);
      insert_level(mdw[1],mdw[0],mds);
   }

   sort_forces();
   nmds=nall_steps(mds)+1;   
}
   

mdstep_t *mdsteps(int *nop,int *itu)
{
   (*nop)=nmds;
   (*itu)=iend-1;

   return mds;
}


static void print_ops(void)
{
   int i,itu;
   
   printf("List of elementary operations:\n");

   itu=iend-1;
   
   for (i=0;i<nmds;i++)
   {
      if (mds[i].iop<itu)
         printf("TP: force %2d, eps = % .2e\n",mds[i].iop,mds[i].eps);
      else if (mds[i].iop==itu)
         printf("TU:           eps = % .2e\n",mds[i].eps);
      else if (mds[i].iop==iend)
         printf("END\n\n");
      else
         error_root(1,1,"print_ops [mdsteps.c]","Unkown operation");
   }
}


static void print_times(double tau)
{
   int i,j,it;
   double seps;
   
   printf("Total integration times:\n");

   for (i=0;i<iend;i++)
   {
      it=0;
      seps=0.0;

      for (j=0;j<nmds;j++)
      {
         if (mds[j].iop==i)
         {
            it=1;
            seps+=mds[j].eps;
         }
      }

      seps/=tau;
            
      if (i==(iend-1))
         printf("TU:       sum(eps)/tau = %.3e\n",seps);
      else if (it==1)
         printf("Force %2d: sum(eps)/tau = %.3e\n",i,seps);
   }

   printf("\n");
}


void print_mdsteps(int ipr)
{
   int my_rank;
   hmc_parms_t hmc;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   hmc=hmc_parms();

   if (my_rank==0)
   {
      printf("Molecular-dynamics integrator:\n\n");

      printf("Trajectory length = %.4e\n",hmc.tau);
      printf("Number of levels = %d\n\n",hmc.nlv);

      print_mdint_parms();

      if (ipr&0x1)
         print_force_parms();

      if (ipr&0x2)
         print_ops();

      if (ipr&0x4)
         print_times(hmc.tau);
   }
}
