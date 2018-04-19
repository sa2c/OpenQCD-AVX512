
/*******************************************************************************
*
* File counters.c
*
* Copyright (C) 2011, 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Solver iteration counters
*
* The externally accessible functions are
* 
*   void setup_counters(void)
*     Creates the counters required for counting the solver iteration
*     numbers along the molecular-dynamics trajectories. The solvers
*     for the Dirac equation used by the action and force programs in
*     the course of the HMC algorithm are inferred from the parameter
*     data base.
*
*   void clear_counters(void)
*     Sets all counters to zero.
*
*   void add2counter(char *type,int idx,int *status)
*     Adds the status numbers "status" to the counter characterized by 
*     "type" and index "idx" (see the notes).
*
*   int get_count(char *type,int idx,int *status)
*     Returns the number of times add2counter(type,idx,status) has been
*     called since the counter was last cleared. On exit the program
*     assigns the sum of the accumulated status values to the argument
*     status [the meaning of the arguments is otherwise the same is in 
*     the case of add2counter()].
*
*   void print_avgstat(char *type,int idx)
*     Prints the average status values of the counter specified by "type"
*     and "idx" [as in add2counter()] to stdout on MPI process 0.
*
*   void print_all_avgstat(void)
*     Prints the average status values of all known counters to stdout on
*     MPI process 0.
*
* Notes:
*
* In most cases, the computation of the fermion actions and forces requires
* the Dirac equation to be solved a number of times. Depending on the solver
* used, the number of status values returned by the solver program may vary.

* The counters administered by this module are set up for all fermion actions
* and forces that take part in the HMC algorithm according to the parameter 
* data base. In addition, the iteration numbers required for the solution of 
* the little Dirac equation in the course of the generation and updates of 
* the deflation subspace are monitored.
*
* The available counter types are "action", "field", "force" and "modes". In
* the first three cases, the index idx passed to add2counter() is the one of
* the action, pseudo-fermion field and force in the parameter data base (see
* flags/{action,force}_parms.c). If type="modes", there are counters for the
* solver iteration numbers required for the generation of the deflation sub-
* space (idx=0), the subspace updates (idx=1) and the subspace regeneration
* (idx=2; see dfl/dfl_sap_gcr.c). The status array passed to add2counter()
* is expected to contain all status values returned by the associated action,
* pseudo-fermion field generation, force and mode-generation program.
*
* When the HMC parameters or the specifications of the actions and forces
* are changed, it may be necessary to call setup_counters() again in order
* to ensure that all counters are properly set up. Except for the program
* setup_counters(), the programs in this module can be called locally.
*
*******************************************************************************/

#define COUNTERS_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "update.h"
#include "global.h"

typedef struct
{
   int n,ns;
   int *status;
} counter_t;

static int nac=0,nfd=0,nfr=0,nmd=0;
static counter_t *act=NULL,*fld=NULL,*frc=NULL,*mds=NULL;


static void free_cnt(int nc,counter_t *cnt)
{
   int i;

   for (i=0;i<nc;i++)
   {
      if (cnt[i].ns>0)
      {
         free(cnt[i].status);
         free(cnt);
         break;
      }
   }
}


static counter_t *alloc_cnt(int nc)
{
   int i;
   counter_t *cnt;
   
   if (nc>0)
   {
      cnt=malloc(nc*sizeof(*cnt));
      error(cnt==NULL,1,"alloc_cnt [counters.c]",
            "Unable to allocate counters");

      for (i=0;i<nc;i++)
      {
         cnt[i].n=0;
         cnt[i].ns=0;
         cnt[i].status=NULL;
      }

      return cnt;
   }
   else
      return NULL;
}


static void set_nc(void)
{
   int i,j,k;
   hmc_parms_t hmc;
   action_parms_t ap;
   mdint_parms_t mdp;
   force_parms_t fp;
   solver_parms_t sp;
   
   hmc=hmc_parms();
   nac=0;
   nfd=0;
   nfr=0;
   nmd=0;

   for (i=0;i<hmc.nact;i++)
   {
      j=hmc.iact[i];
      ap=action_parms(j);

      if ((ap.action==ACF_TM1)||
          (ap.action==ACF_TM1_EO)||
          (ap.action==ACF_TM1_EO_SDET)||
          (ap.action==ACF_TM2)||
          (ap.action==ACF_TM2_EO)||
          (ap.action==ACF_RAT)||
          (ap.action==ACF_RAT_SDET))
      {
         if (j>=nac)
            nac=j+1;

         sp=solver_parms(ap.isp[0]);
         if (sp.solver==DFL_SAP_GCR)
            nmd=3;

         if ((ap.action==ACF_TM2)||
             (ap.action==ACF_TM2_EO))
         {
            if (ap.ipf>=nfd)
               nfd=ap.ipf+1;

            sp=solver_parms(ap.isp[1]);
            if (sp.solver==DFL_SAP_GCR)
               nmd=3;
         }

         if ((ap.action==ACF_RAT)||
             (ap.action==ACF_RAT_SDET))
         {
            if (ap.ipf>=nfd)
               nfd=ap.ipf+1;
         }
      }
   }  

   for (i=0;i<hmc.nlv;i++)
   {
      mdp=mdint_parms(i);

      for (j=0;j<mdp.nfr;j++)
      {
         k=mdp.ifr[j];
         fp=force_parms(k);

         if ((fp.force==FRF_TM1)||
             (fp.force==FRF_TM1_EO)||
             (fp.force==FRF_TM1_EO_SDET)||
             (fp.force==FRF_TM2)||
             (fp.force==FRF_TM2_EO)||
             (fp.force==FRF_RAT)||
             (fp.force==FRF_RAT_SDET))
         {
            if (k>=nfr)
               nfr=k+1;

            sp=solver_parms(fp.isp[0]);
            if (sp.solver==DFL_SAP_GCR)
               nmd=3;
         }
      }
   }  
}


static void set_ns(void)
{
   int i,j,k;
   hmc_parms_t hmc;
   mdint_parms_t mdp;
   action_parms_t ap;
   force_parms_t fp;
   solver_parms_t sp;
   
   hmc=hmc_parms();

   for (i=0;i<hmc.nact;i++)
   {
      j=hmc.iact[i];
      ap=action_parms(j);

      if ((ap.action==ACF_TM1)||
          (ap.action==ACF_TM1_EO)||
          (ap.action==ACF_TM1_EO_SDET)||
          (ap.action==ACF_TM2)||
          (ap.action==ACF_TM2_EO)||
          (ap.action==ACF_RAT)||
          (ap.action==ACF_RAT_SDET))
      {
         sp=solver_parms(ap.isp[0]);

         if ((sp.solver==CGNE)||(sp.solver==MSCG)||(sp.solver==SAP_GCR))
            act[j].ns=1;
         else if (sp.solver==DFL_SAP_GCR)
            act[j].ns=2;
         else
            error_root(1,1,"set_ns [counters.c]","Unknown solver");
         
         if ((ap.action==ACF_TM2)||
             (ap.action==ACF_TM2_EO))
         {
            k=ap.ipf;
            sp=solver_parms(ap.isp[1]);

            if ((sp.solver==CGNE)||(sp.solver==MSCG)||(sp.solver==SAP_GCR))
               fld[k].ns=1;
            else if (sp.solver==DFL_SAP_GCR)
               fld[k].ns=2;
            else
               error_root(1,1,"set_ns [counters.c]","Unknown solver");
         }
         else if ((ap.action==ACF_RAT)||
                  (ap.action==ACF_RAT_SDET))
         {
            k=ap.ipf;
            sp=solver_parms(ap.isp[0]);

            if ((sp.solver==CGNE)||(sp.solver==MSCG)||(sp.solver==SAP_GCR))
               fld[k].ns=1;
            else if (sp.solver==DFL_SAP_GCR)
               fld[k].ns=2;
            else
               error_root(1,1,"set_ns [counters.c]","Unknown solver");
         }
      }
      else if (ap.action!=ACG)
         error_root(1,1,"set_ns [counters.c]","Unknown action");
   }  

   for (i=0;i<hmc.nlv;i++)
   {
      mdp=mdint_parms(i);

      for (j=0;j<mdp.nfr;j++)
      {
         k=mdp.ifr[j];
         fp=force_parms(k);

         if ((fp.force==FRF_TM1)||
             (fp.force==FRF_TM1_EO)||
             (fp.force==FRF_TM1_EO_SDET)||
             (fp.force==FRF_TM2)||
             (fp.force==FRF_TM2_EO)||
             (fp.force==FRF_RAT)||
             (fp.force==FRF_RAT_SDET))
         {
            sp=solver_parms(fp.isp[0]);

            if ((sp.solver==CGNE)||(sp.solver==MSCG))
               frc[k].ns=1;
            else if (sp.solver==SAP_GCR)
               frc[k].ns=2;
            else if (sp.solver==DFL_SAP_GCR)
               frc[k].ns=4;
            else
               error_root(1,1,"set_ns [counters.c]","Unknown solver");
         }
         else if (fp.force!=FRG)
            error_root(1,1,"set_ns [counters.c]","Unknown force");
      }
   }  

   if (nmd>0)
   {
      mds[0].ns=1;
      mds[1].ns=1;
      mds[2].ns=1;
   }
}


static void alloc_stat(int nc,counter_t *cnt)
{
   int i,ns,*stat;

   if (nc>0)
   {
      ns=0;

      for (i=0;i<nc;i++)
         ns+=cnt[i].ns;

      stat=malloc(ns*sizeof(*stat));
      error(stat==NULL,1,"alloc_stat [counters.c]",
            "Unable to allocate status arrays");

      for (i=0;i<nc;i++)
      {
         if (cnt[i].ns>0)
         {
            cnt[i].status=stat;
            stat+=cnt[i].ns;
         }
      }
   }
}


void setup_counters(void)
{
   free_cnt(nac,act);
   free_cnt(nfd,fld);
   free_cnt(nfr,frc);
   free_cnt(nmd,mds);

   set_nc();
   act=alloc_cnt(nac);
   fld=alloc_cnt(nfd);
   frc=alloc_cnt(nfr);
   mds=alloc_cnt(nmd);   

   set_ns();
   alloc_stat(nac,act);
   alloc_stat(nfd,fld);
   alloc_stat(nfr,frc);
   alloc_stat(nmd,mds);

   clear_counters();
}


static void set_cnt2zero(int nc,counter_t *cnt)
{
   int i,j,ns,*stat;
   
   for (i=0;i<nc;i++)
   {
      cnt[i].n=0;
      ns=cnt[i].ns;
      stat=cnt[i].status;
      
      for (j=0;j<ns;j++)
         stat[j]=0;
   }
}


void clear_counters(void)
{
   set_cnt2zero(nac,act);
   set_cnt2zero(nfd,fld);
   set_cnt2zero(nfr,frc);
   set_cnt2zero(nmd,mds);   
}


void add2counter(char *type,int idx,int *status)
{
   int i,nc,ns,*stat;
   counter_t *cnt;
   
   if (strcmp(type,"force")==0)
   {
      nc=nfr;
      cnt=frc;
   }
   else if (strcmp(type,"modes")==0)
   {
      nc=nmd;
      cnt=mds;
   }   
   else if (strcmp(type,"action")==0)
   {
      nc=nac;
      cnt=act;
   }
   else if (strcmp(type,"field")==0)
   {
      nc=nfd;
      cnt=fld;
   }   
   else
   {
      error_loc(1,1,"add2counter [counters.c]","Unknown counter type");
      return;
   }
   
   if ((idx>=0)&&(idx<nc))
   {
      cnt[idx].n+=1;
      ns=cnt[idx].ns;
      stat=cnt[idx].status;

      for (i=0;i<ns;i++)
         stat[i]+=status[i];
   }
   else
      error_loc(1,1,"add2counter [counters.c]","Counter index out of range");
}


int get_count(char *type,int idx,int *status)
{
   int i,nc,n,ns,*stat;
   counter_t *cnt;
   
   if (strcmp(type,"force")==0)
   {
      nc=nfr;
      cnt=frc;
   }
   else if (strcmp(type,"modes")==0)
   {
      nc=nmd;
      cnt=mds;
   }   
   else if (strcmp(type,"action")==0)
   {
      nc=nac;
      cnt=act;
   }
   else if (strcmp(type,"field")==0)
   {
      nc=nfd;
      cnt=fld;
   }   
   else
   {
      error_loc(1,1,"get_count [counters.c]","Unknown counter type");
      return 0;
   }
   
   if ((idx>=0)&&(idx<nc))
   {
      n=cnt[idx].n;
      ns=cnt[idx].ns;
      stat=cnt[idx].status;

      for (i=0;i<ns;i++)
         status[i]=stat[i];

      return n;
   }
   else
      error_loc(1,1,"get_count [counters.c]","Counter index out of range");

   return 0;
}


void print_avgstat(char *type,int idx)
{
   int my_rank;
   int i,n,ns,*stat;
   double r;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   
   if (my_rank==0)
   {
      n=0;
      ns=0;
      stat=NULL;
      
      if (strcmp(type,"force")==0)
      {
         if ((idx>=0)&&(idx<nfr))
         {
            n=frc[idx].n;
            ns=frc[idx].ns;
            stat=frc[idx].status;
            printf("Force  %2d: <status> = ",idx);
         }
      }
      else if (strcmp(type,"modes")==0)
      {
         if ((idx>=0)&&(idx<nmd))
         {         
            n=mds[idx].n;
            ns=mds[idx].ns;
            stat=mds[idx].status;
            printf("Modes  %2d: <status> = ",idx);
         }
      }
      else if (strcmp(type,"action")==0)
      {
         if ((idx>=0)&&(idx<nac))
         {          
            n=act[idx].n;
            ns=act[idx].ns;
            stat=act[idx].status;
            printf("Action %2d: <status> = ",idx);
         }
      }
      else if (strcmp(type,"field")==0)
      {
         if ((idx>=0)&&(idx<nac))
         {          
            n=fld[idx].n;
            ns=fld[idx].ns;
            stat=fld[idx].status;
            printf("Field  %2d: <status> = ",idx);
         }
      }
      else
      {
         error_loc(1,1,"print_avgstat [counters.c]","Unknown counter type");
         return;
      }

      if (ns>0)
      {
         if ((strcmp(type,"modes")==0)&&(idx==0))
         {
            n+=mds[2].n;

            if (n>0)
               r=1.0/(double)(n);
            else
               r=1.0;

            printf("%d",(int)((double)(stat[0]+mds[2].status[0])*r+0.5));

            if (mds[2].n>0)
               printf(" (no of regenerations = %d)",mds[2].n);
         }
         else
         {
            if (n>0)
               r=1.0/(double)(n);
            else
               r=1.0;
      
            printf("%d",(int)((double)(stat[0])*r+0.5));

            for (i=1;i<ns;i++)
               printf(",%d",(int)((double)(stat[i])*r+0.5));

            if ((strcmp(type,"modes")==0)&&(idx==1))
               printf(" (no of updates = %d)",n);
         }
         
         printf("\n");
      }
   }
}


void print_all_avgstat(void)
{
   int i;

   for (i=0;i<nac;i++)
   {
      if (act[i].ns>0)
         print_avgstat("action",i);
   }

   for (i=0;i<nfd;i++)
   {
      if (fld[i].ns>0)
         print_avgstat("field",i);
   }
   
   for (i=0;i<nfr;i++)
   {
      if (frc[i].ns>0)
         print_avgstat("force",i);
   }

   for (i=0;i<(nmd-1);i++)
   {
      if (mds[i].ns>0)
         print_avgstat("modes",i);
   }      
}
