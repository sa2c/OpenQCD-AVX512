
/*******************************************************************************
*
* File solver_parms.c
*
* Copyright (C) 2011, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Solver parameter data base
*
* The externally accessible functions are
*
*   solver_parms_t set_solver_parms(int isp,solver_t solver,
*                                   int nkv,int isolv,int nmr,int ncy,
*                                   int nmx,double res)
*     Sets the parameters in the solver parameter set number isp and returns
*     a structure containing them (see the notes).
*
*   solver_parms_t solver_parms(int isp)
*     Returns a structure containing the solver parameter set number
*     isp (see the notes).
*
*   void read_solver_parms(int isp)
*     On process 0, this program scans stdin for a line starting with the
*     string "[Solver <int>]" (after any number of blanks), where <int> is
*     the integer value passed by the argument. An error occurs if no such
*     line or more than one is found. The lines 
*
*       solver  <solver_t>
*       nkv     <int>
*       isolv   <int>
*       nmr     <int>
*       ncy     <int>
*       nmx     <int>
*       res     <double>
*
*     are then read one by one using read_line() [utils/mutils.c]. The 
*     lines with tags nkv,..,ncy may be absent in the case of the CGNE
*     and MSCG solvers (see the notes). The data are then added to the 
*     data base by calling set_solver_parms(isp,...).
*
*   void print_solver_parms(int *isap,int *idfl)
*     Prints the parameters of the defined solvers to stdout on MPI
*     process 0. On exit the flag isap is 1 or 0 depending on whether
*     one of the solvers makes use of the Schwarz Alternating Procedure
*     (SAP) or not. Similarly, the flag idfl is set 1 or 0 depending on 
*     whether deflation is used or not. On MPI processes other than 0,
*     the program does nothing and sets isap and idfl to zero.
*
*   void write_solver_parms(FILE *fdat)
*     Writes the parameters of the defined solvers to the file fdat on
*     MPI process 0.
*
*   void check_solver_parms(FILE *fdat)
*     Compares the parameters of the defined solvers with those stored
*     on the file fdat on MPI process 0, assuming the latter were written
*     to the file by the program write_solver_parms() (mismatches of the
*     maximal solver iteration number are not considered to be an error).
*
* Notes:
*
* The elements of a structure of type solver_parms_t are
*
*   solver  Solver program used. This parameter is an enum type with
*           one of the following values:
*
*            CGNE           Program tmcg() [forces/tmcg.c].
*
*            MSCG           Program tmcgm() [forces/tmcgm.c].
*
*            SAP_GCR        Program sap_gcr() [sap/sap_gcr.c].
*
*            DFL_SAP_GCR    Program dfl_sap_gcr() [dfl/dfl_sap_gcr.c].
*
*   nkv     Maximal number of Krylov vectors generated before the GCR
*           algorithm is restarted if solver=*_GCR.
*
*   isolv   Block solver to be used if solver=*SAP_GCR (0: plain MinRes,
*           1: eo-preconditioned MinRes).
*
*   nmr     Number of block solver iterations if solver=*SAP_GCR.
*
*   ncy     Number of SAP cycles to be applied if solver=*SAP_GCR.
*
*   nmx     Maximal number of CG iterations if solver={CGNE,MSCG} or 
*           maximal total number of Krylov vectors that may be generated 
*           if solver={SAP_GCR,DFL_SAP_GCR}.
*
*   res     Desired maximal relative residue of the calculated solution.
*
* Depending on the solver, some parameters are not used. These are set to
* zero by the program set_solver_parms() independently of the values of
* the arguments.
*
* Up to 32 solver parameter sets, labeled by an index isp=0,1,..,31, can
* be specified. Once a set is specified, it cannot be changed by calling
* set_solver_parms() again. Solver parameters must be globally the same.
*
* Except for solver_parms(), the programs in this module perform global
* operations and must be called simultaneously on all MPI processes.
*
*******************************************************************************/

#define SOLVER_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "global.h"

#define ISPMAX 32

static int init=0;
static solver_t solver[]={CGNE,MSCG,SAP_GCR,DFL_SAP_GCR};
static solver_parms_t sp[ISPMAX+1]={{SOLVERS,0,0,0,0,0,0.0}};


static void init_sp(void)
{
   int i;
   
   for (i=1;i<=ISPMAX;i++)
      sp[i]=sp[0];

   init=1;
}


solver_parms_t set_solver_parms(int isp,solver_t solver,
                                int nkv,int isolv,int nmr,int ncy,
                                int nmx,double res)
{
   int ie,iprms[7];
   double dprms[1];
   
   if (init==0)
      init_sp();

   if ((solver==CGNE)||(solver==MSCG))
   {
      nkv=0;
      isolv=0;
      nmr=0;
      ncy=0;
   }
   
   if (NPROC>1)
   {
      iprms[0]=isp;
      iprms[1]=(int)(solver);
      iprms[2]=nkv;
      iprms[3]=isolv;
      iprms[4]=nmr;      
      iprms[5]=ncy;
      iprms[6]=nmx;
      dprms[0]=res;

      MPI_Bcast(iprms,7,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      ie=0;
      ie|=(iprms[0]!=isp);
      ie|=(iprms[1]!=(int)(solver));
      ie|=(iprms[2]!=nkv);
      ie|=(iprms[3]!=isolv);
      ie|=(iprms[4]!=nmr);
      ie|=(iprms[5]!=ncy);
      ie|=(iprms[6]!=nmx);
      ie|=(dprms[0]!=res);
      
      error(ie!=0,1,"set_solver_parms [solver_parms.c]",
            "Parameters are not global");
   }

   ie=0;
   ie|=(isp<0)||(isp>=ISPMAX);
   ie|=(solver==SOLVERS);
   ie|=(nmx<1);

   if ((solver==SAP_GCR)||(solver==DFL_SAP_GCR))
   {
      ie|=(isolv<0)||(isolv>1);
      ie|=(nmr<1);
      ie|=(ncy<1);
   }

   error_root(ie!=0,1,"set_solver_parms [solver_parms.c]",
              "Parameters are out of range");
   
   error_root(sp[isp].solver!=SOLVERS,1,"set_solver_parms [solver_parms.c]",
              "Attempt to reset an already specified solver parameter set");

   sp[isp].solver=solver;
   sp[isp].nkv=nkv;
   sp[isp].isolv=isolv;
   sp[isp].nmr=nmr;
   sp[isp].ncy=ncy;
   sp[isp].nmx=nmx;
   sp[isp].res=res;

   return sp[isp];
}


solver_parms_t solver_parms(int isp)
{
   if (init==0)
      init_sp();

   if ((isp>=0)&&(isp<ISPMAX))
      return sp[isp];
   else
   {
      error_loc(1,1,"solver_parms [solver_parms.c]",
                "Solver index is out of range");
      return sp[ISPMAX];
   }
}


void read_solver_parms(int isp)
{
   int my_rank,ids;
   int nkv,isolv,nmr,ncy,nmx;
   double res;
   char line[NAME_SIZE];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   ids=0;
   nkv=0;
   isolv=0;
   nmr=0;
   ncy=0;
   nmx=0;
   res=1.0;
   
   if (my_rank==0)
   {
      sprintf(line,"Solver %d",isp);
      find_section(line);

      read_line("solver","%s",line);

      if (strcmp(line,"MSCG")==0)
         ids=1;
      else if (strcmp(line,"SAP_GCR")==0)
      {
         ids=2;
         read_line("nkv","%d",&nkv);
         read_line("isolv","%d",&isolv);
         read_line("nmr","%d",&nmr);
         read_line("ncy","%d",&ncy);
      }
      else if (strcmp(line,"DFL_SAP_GCR")==0)
      {
         ids=3;
         read_line("nkv","%d",&nkv);
         read_line("isolv","%d",&isolv);
         read_line("nmr","%d",&nmr);
         read_line("ncy","%d",&ncy);
      }
      else if (strcmp(line,"CGNE")!=0)
         error_root(1,1,"read_solver_parms [solver_parms.c]",
                    "Unknown solver %s",line);

      read_line("nmx","%d",&nmx);
      read_line("res","%lf",&res);      
   }

   if (NPROC>1)
   {     
      MPI_Bcast(&ids,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nkv,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&isolv,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nmr,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&ncy,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nmx,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);   
   }
   
   set_solver_parms(isp,solver[ids],nkv,isolv,nmr,ncy,nmx,res);
}


void print_solver_parms(int *isap,int *idfl)
{
   int my_rank,i;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   (*isap)=0;
   (*idfl)=0;
   
   if ((my_rank==0)&&(init==1))
   {
      for (i=0;i<ISPMAX;i++)
      {
         if (sp[i].solver!=SOLVERS)
         {
            printf("Solver %d:\n",i);
            
            if (sp[i].solver==CGNE)
            {
               printf("CGNE solver\n");
               printf("nmx = %d\n",sp[i].nmx);
               printf("res = %.1e\n\n",sp[i].res);
            }
            else if (sp[i].solver==MSCG)
            {
               printf("MSCG solver\n");
               printf("nmx = %d\n",sp[i].nmx);
               printf("res = %.1e\n\n",sp[i].res);
            }
            else if (sp[i].solver==SAP_GCR)
            {
               (*isap)=1;
               printf("SAP_GCR solver\n");
               printf("nkv = %d\n",sp[i].nkv);
               printf("isolv = %d\n",sp[i].isolv);
               printf("nmr = %d\n",sp[i].nmr);              
               printf("ncy = %d\n",sp[i].ncy);
               printf("nmx = %d\n",sp[i].nmx);
               printf("res = %.1e\n\n",sp[i].res);
            }
            else if (sp[i].solver==DFL_SAP_GCR)
            {
               (*isap)=1;
               (*idfl)=1;
               printf("DFL_SAP_GCR solver\n");
               printf("nkv = %d\n",sp[i].nkv);
               printf("isolv = %d\n",sp[i].isolv);
               printf("nmr = %d\n",sp[i].nmr);              
               printf("ncy = %d\n",sp[i].ncy);
               printf("nmx = %d\n",sp[i].nmx);
               printf("res = %.1e\n\n",sp[i].res);
            }
            else
               printf("UNKNOWN solver\n\n");
         }
      }
   }
}


void write_solver_parms(FILE *fdat)
{
   int my_rank,endian;
   int iw,i;
   stdint_t istd[7];
   double dstd[1];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   
   if ((my_rank==0)&&(init==1))
   {
      for (i=0;i<ISPMAX;i++)
      {
         if (sp[i].solver!=SOLVERS)
         {
            istd[0]=(stdint_t)(i);            
            istd[1]=(stdint_t)(sp[i].solver);
            istd[2]=(stdint_t)(sp[i].nmx);
            istd[3]=(stdint_t)(sp[i].nkv);
            istd[4]=(stdint_t)(sp[i].isolv);
            istd[5]=(stdint_t)(sp[i].nmr);
            istd[6]=(stdint_t)(sp[i].ncy);            
            dstd[0]=sp[i].res;

            if (endian==BIG_ENDIAN)
            {
               bswap_int(7,istd);
               bswap_double(1,dstd);
            }
            
            iw=fwrite(istd,sizeof(stdint_t),7,fdat);
            iw+=fwrite(dstd,sizeof(double),1,fdat);
            error_root(iw!=8,1,"write_solver_parms [solver_parms.c]",
                       "Incorrect write count");
         }
      }
   }
}


void check_solver_parms(FILE *fdat)
{
   int my_rank,endian;
   int ir,ie,i;
   stdint_t istd[7];
   double dstd[1];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   
   if ((my_rank==0)&&(init==1))
   {
      ie=0;
      
      for (i=0;i<ISPMAX;i++)
      {
         if (sp[i].solver!=SOLVERS)
         {
            ir=fread(istd,sizeof(stdint_t),7,fdat);
            ir+=fread(dstd,sizeof(double),1,fdat);
            error_root(ir!=8,1,"check_solver_parms [solver_parms.c]",
                       "Incorrect read count");

            if (endian==BIG_ENDIAN)
            {
               bswap_int(7,istd);
               bswap_double(1,dstd);
            }
            
            ie|=(istd[0]!=(stdint_t)(i));            
            ie|=(istd[1]!=(stdint_t)(sp[i].solver));
            ie|=(istd[3]!=(stdint_t)(sp[i].nkv));
            ie|=(istd[4]!=(stdint_t)(sp[i].isolv));
            ie|=(istd[5]!=(stdint_t)(sp[i].nmr));
            ie|=(istd[6]!=(stdint_t)(sp[i].ncy));
            ie|=(dstd[0]!=sp[i].res);
         }
      }
         
      error_root(ie!=0,1,"check_solver_parms [solver_parms.c]",
                 "Parameters do not match");         
   }
}
