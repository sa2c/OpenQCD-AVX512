
/*******************************************************************************
*
* File mdint_parms.c
*
* Copyright (C) 2011, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Molecular-dynamics integrator data base
*
* The externally accessible functions are
*
*   mdint_parms_t set_mdint_parms(int ilv,integrator_t integrator,double lambda,
*                                 int nstep,int nfr,int *ifr)
*     Sets the parameters of the molecular-dynamics integrator at level
*     ilv and returns a structure containing them (see the notes).
*
*   mdint_parms_t mdint_parms(int ilv)
*     Returns a structure containing the parameters of the integrator at
*     level ilv (see the notes).
*
*   void read_mdint_parms(int ilv)
*     On process 0, this program scans stdin for a line starting with the
*     string "[Level <int>]" (after any number of blanks), where <int> is
*     the integer value passed by the argument. An error occurs if no such
*     line or more than one is found. The lines 
*
*       integrator   <integrator_t>
*       lambda       <double>
*       nstep        <int>
*       forces       <int> [<int>]
*
*     are then read using read_line() [utils/mutils.c]. The line tagged
*     "lambda" is required only when the specified integrator is the 2nd
*     order OMF integrator. The line tagged "forces" must contain the
*     indices of the forces (separated by white space) that are to be
*     integrated at this level. On exit, the data are entered in the data
*     base by calling set_mdint_parms(ilv,...).
*
*   void print_mdint_parms(void)
*     Prints the parameters of the defined integrator levels to stdout
*     on MPI process 0.
*
*   void write_mdint_parms(FILE *fdat)
*     Writes the parameters of the defined integrator levels to the file
*     fdat on MPI process 0.
*
*   void check_mdint_parms(FILE *fdat)
*     Compares the parameters of the defined integrator levels with those
*     stored on the file fdat on MPI process 0, assuming the latter were
*     written to the file by the program write_mdint_parms().
*
* Notes:
*
* A structure of type mdint_parms_t contains the parameters of a hierarchical
* molecular-dynamics integrator at a specified level (see update/README.mdint).
* Its elements are
*
*   integrator   Elementary integrator used. This parameter is an enum
*                type with one of the following values:
*
*                  LPFR     Leapfrog integrator
*
*                  OMF2     2nd order Omelyan-Mryglod-Folk integrator
*
*                  OMF4     4th order Omelyan-Mryglod-Folk integrator
*
*   lambda       Parameter of the 2nd order OMF integrator
*
*   nstep        Number of times the elementary integrator is applied
*                at this level
*
*   nfr          Number of forces integrated at this level
*
*   ifr          Force indices ifr[i] (i=0,..,nfr-1)
*
* The parameter lambda is not used in the case of the leapfrog and the 4th
* order OMF integrator. Up to 32 integrator levels, labeled by an index
* ilv=0,1,..,31, can be specified.
*
* An example of valid section in an input file which can be read by calling
* read_mdint(3) is
*
*  [Level 3]
*  integrator OMF2
*  lambda     0.2
*  nstep      12
*  forces     2 4 5
*
* In this case, there are three forces with index 2, 4 and 5.
*
* The programs set_mdint_parms() and read_mdint_parms() perform global
* operations and must be called simultaneously on all MPI processes.
*
*******************************************************************************/

#define MDINT_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "global.h"

#define ILVMAX 32

static int init=0;
static mdint_parms_t mdp[ILVMAX+1]={{INTEGRATORS,0.0,0,0,NULL}};


static void init_mdp(void)
{
   int i;
   
   for (i=1;i<=ILVMAX;i++)
      mdp[i]=mdp[0];

   init=1;
}


static void alloc_ifr(int ilv,int nfr)
{
   int *ifr;

   if (mdp[ilv].nfr>0)
   {
      free(mdp[ilv].ifr);
      mdp[ilv].nfr=0;
      mdp[ilv].ifr=NULL;
   }

   if (nfr>0)
   {
      ifr=malloc(nfr*sizeof(*ifr));
      error(ifr==NULL,1,"alloc_ifr [mdint_parms.c]",
            "Unable to allocate index array");
      mdp[ilv].nfr=nfr;
      mdp[ilv].ifr=ifr;
   }
}


mdint_parms_t set_mdint_parms(int ilv,integrator_t integrator,double lambda,
                              int nstep,int nfr,int *ifr)
{
   int iprms[4],i,j,ie;
   double dprms[1];
   
   if (init==0)
      init_mdp();

   if (integrator!=OMF2)
      lambda=0.0;
   
   if (NPROC>1)
   {
      iprms[0]=ilv;
      iprms[1]=(int)(integrator);
      iprms[2]=nstep;
      iprms[3]=nfr;
      dprms[0]=lambda;

      MPI_Bcast(iprms,4,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      
      ie=0;
      ie|=(iprms[0]!=ilv);
      ie|=(iprms[1]!=(int)(integrator));
      ie|=(iprms[2]!=nstep);
      ie|=(iprms[3]!=nfr);
      ie|=(dprms[0]!=lambda);

      for (i=0;i<nfr;i++)
      {
         iprms[0]=ifr[i];
         
         MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);         

         ie|=(iprms[0]!=ifr[i]);
      }     
      
      error(ie!=0,1,"set_mdint_parms [mdint_parms.c]",
            "Parameters are not global");
   }

   ie=0;
   ie|=(ilv<0)||(ilv>=ILVMAX);
   ie|=(integrator==INTEGRATORS);
   ie|=(nstep<1);
   ie|=(nfr<1);
   
   for (i=0;i<nfr;i++)
      ie|=(ifr[i]<0);
   
   error_root(ie!=0,1,"set_mdint_parms [mdint_parms.c]",
              "Parameters are out of range");

   for (i=0;i<nfr;i++)
   {
      for (j=(i+1);j<nfr;j++)
         ie|=(ifr[i]==ifr[j]);
   }
   
   error_root(ie!=0,1,"set_mdint_parms [mdint_parms.c]",
              "Dublicate force indices");
   
   mdp[ilv].integrator=integrator;
   mdp[ilv].lambda=lambda;   
   mdp[ilv].nstep=nstep;
   alloc_ifr(ilv,nfr);
   
   for (i=0;i<nfr;i++)
      mdp[ilv].ifr[i]=ifr[i];
   
   return mdp[ilv];
}


mdint_parms_t mdint_parms(int ilv)
{
   if (init==0)
      init_mdp();      

   if ((ilv>=0)&&(ilv<ILVMAX))
      return mdp[ilv];
   else
   {
      error_loc(1,1,"mdint_parms [mdint_parms.c]",
                "Level index is out of range");
      return mdp[ILVMAX];
   }
}


void read_mdint_parms(int ilv)
{
   int my_rank,idi;
   int nstep,nfr,*ifr;
   double lambda;
   char line[NAME_SIZE];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      sprintf(line,"Level %d",ilv);
      find_section(line);

      read_line("integrator","%s",line);

      if (strcmp(line,"LPFR")==0)
      {
         idi=0;
         lambda=0.0;
      }
      else if (strcmp(line,"OMF2")==0)
      {
         idi=1;
         read_line("lambda","%lf",&lambda);
      }
      else if (strcmp(line,"OMF4")==0)
      {
         idi=2;
         lambda=0.0;
      }
      else
         error_root(1,1,"read_mdint_parms [mdint_parms.c]",
                    "Unknown integrator %s",line);
      
      read_line("nstep","%d",&nstep);
      error_root(nstep<1,1,"read_mdint [mdint_parms.c]",
                 "Parameter nstep out of range");
      
      nfr=count_tokens("forces");
      error_root(nfr==0,1,"read_mdint [mdint_parms.c]",
                 "No forces specified");
   }

   MPI_Bcast(&idi,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&lambda,1,MPI_DOUBLE,0,MPI_COMM_WORLD);      
   MPI_Bcast(&nstep,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nfr,1,MPI_INT,0,MPI_COMM_WORLD);
   
   ifr=malloc(nfr*sizeof(*ifr));
   error(ifr==NULL,1,"read_mdint [mdint_parms.c]",
         "Unable to allocate index array");
   if (my_rank==0)
      read_iprms("forces",nfr,ifr);

   MPI_Bcast(ifr,nfr,MPI_INT,0,MPI_COMM_WORLD);
   
   if (idi==0)
      set_mdint_parms(ilv,LPFR,lambda,nstep,nfr,ifr);
   else if (idi==1)
      set_mdint_parms(ilv,OMF2,lambda,nstep,nfr,ifr);
   else if (idi==2)
      set_mdint_parms(ilv,OMF4,lambda,nstep,nfr,ifr);
   
   free(ifr);
}


void print_mdint_parms(void)
{
   int my_rank,i,j,n;
   int nfr,*ifr;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   
   if ((my_rank==0)&&(init==1))
   {
      for (i=0;i<ILVMAX;i++)
      {
         if (mdp[i].integrator!=INTEGRATORS)
         {
            printf("Level %d:\n",i);

            if (mdp[i].integrator==LPFR)
               printf("Leapfrog integrator\n");
            else if (mdp[i].integrator==OMF2)
            {
               n=fdigits(mdp[i].lambda);
               printf("2nd order OMF integrator with lambda = %.*f\n",
                      IMAX(n,1),mdp[i].lambda);
            }
            else if (mdp[i].integrator==OMF4)
               printf("4th order OMF integrator\n");
            else
               printf("Unknown integrator\n");

            printf("Number of steps = %d\n",mdp[i].nstep);

            nfr=mdp[i].nfr;
            ifr=mdp[i].ifr;

            if (nfr>0)
            {
               printf("Forces =");
         
               for (j=0;j<nfr;j++)
                  printf(" %d",ifr[j]);

               printf("\n");
            }

            printf("\n");
         }
      }
   }
}


void write_mdint_parms(FILE *fdat)
{
   int my_rank,endian;
   int nfr,nmx,n,iw,i,j;
   stdint_t *istd;
   double dstd[1];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   
   if ((my_rank==0)&&(init==1))
   {
      nmx=0;
      
      for (i=0;i<ILVMAX;i++)
      {
         nfr=mdp[i].nfr;
         if (nfr>nmx)
            nmx=nfr;
      }
      
      istd=malloc((nmx+4)*sizeof(stdint_t));
      error_root(istd==NULL,1,"write_mdint_parms [mdint_parms.c]",
                 "Unable to allocate auxiliary array");
      
      for (i=0;i<ILVMAX;i++)
      {
         if (mdp[i].integrator!=INTEGRATORS)
         {
            nfr=mdp[i].nfr;
            
            istd[0]=(stdint_t)(i);            
            istd[1]=(stdint_t)(mdp[i].integrator);
            istd[2]=(stdint_t)(mdp[i].nstep);
            istd[3]=(stdint_t)(mdp[i].nfr);

            for (j=0;j<nfr;j++)
               istd[4+j]=(stdint_t)(mdp[i].ifr[j]);

            dstd[0]=mdp[i].lambda;
            n=4+nfr;            
            
            if (endian==BIG_ENDIAN)
            {
               bswap_int(n,istd);
               bswap_double(1,dstd);
            }

            iw=fwrite(istd,sizeof(stdint_t),n,fdat);         
            iw+=fwrite(dstd,sizeof(double),1,fdat);
            error_root(iw!=(n+1),1,"write_mdint_parms [mdint_parms.c]",
                       "Incorrect write count");
         }
      }

      free(istd);
   }
}


void check_mdint_parms(FILE *fdat)
{
   int my_rank,endian;
   int nfr,nmx,n,ir,ie,i,j;
   stdint_t *istd;
   double dstd[1];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   
   if ((my_rank==0)&&(init==1))
   {
      ie=0;
      nmx=0;
      
      for (i=0;i<ILVMAX;i++)
      {
         nfr=mdp[i].nfr;
         if (nfr>nmx)
            nmx=nfr;
      }

      istd=malloc((nmx+4)*sizeof(stdint_t));
      error_root(istd==NULL,1,"check_mdint_parms [mdint_parms.c]",
                 "Unable to allocate auxiliary array");
      
      for (i=0;i<ILVMAX;i++)
      {
         if (mdp[i].integrator!=INTEGRATORS)
         {
            nfr=mdp[i].nfr;
            n=4+nfr;
            
            ir=fread(istd,sizeof(stdint_t),n,fdat);         
            ir+=fread(dstd,sizeof(double),1,fdat);
            error_root(ir!=(n+1),1,"check_mdint_parms [mdint_parms.c]",
                       "Incorrect read count");

            if (endian==BIG_ENDIAN)
            {
               bswap_int(n,istd);
               bswap_double(1,dstd);
            }
            
            ie|=(istd[0]!=(stdint_t)(i));            
            ie|=(istd[1]!=(stdint_t)(mdp[i].integrator));
            ie|=(istd[2]!=(stdint_t)(mdp[i].nstep));
            ie|=(istd[3]!=(stdint_t)(mdp[i].nfr));

            for (j=0;j<nfr;j++)
               ie|=(istd[4+j]!=(stdint_t)(mdp[i].ifr[j]));

            ie|=(dstd[0]!=mdp[i].lambda);
         }
      }
         
      error_root(ie!=0,1,"check_mdint_parms [mdint_parms.c]",
                 "Parameters do not match");         
      free(istd);
   }
}
