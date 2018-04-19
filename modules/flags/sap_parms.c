
/*******************************************************************************
*
* File sap_parms.c
*
* Copyright (C) 2009, 2010, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* SAP parameters
*
* The externally accessible functions are
*
*   sap_parms_t set_sap_parms(int *bs,int isolv,int nmr,int ncy)
*     Sets the parameters of the SAP preconditioner. The parameters are
*
*       bs[4]         Sizes of the blocks in SAP_BLOCKS block grid.
*
*       isolv         Block solver to be used (0: plain MinRes,
*                     1: eo-preconditioned MinRes).
*
*       nmr           Number of block solver iterations.
*
*       ncy           Number of SAP cycles to be applied.
*
*     The return value is a structure that contains the parameters of the
*     SAP preconditioners. The block sizes bs[4] can only be set once, but
*     the values of the other parameters may be changed by calling the
*     program again.
*
*   sap_parms_t sap_parms(void)
*     Returns the parameters currently set for the SAP preconditioner.
*
*   void print_sap_parms(int ipr)
*     Prints the SAP parameters to stdout on MPI process 0. Depending
*     on whether ipr!=0 or 0, the full information is printed or only
*     the block size.
*
*   void write_sap_parms(FILE *fdat)
*     Writes the SAP parameters to the file fdat on MPI process 0.
*
*   void check_sap_parms(FILE *fdat)
*     Compares the SAP parameters with the values stored on the file fdat
*     on MPI process 0, assuming the latter were written to the file by
*     the program write_sap_parms().
*
* Notes:
*
* To ensure the consistency of the data base, the parameters must be set
* simultaneously on all processes. The type sap_parms_t is defined in the
* file flags.h.
*
*******************************************************************************/

#define SAP_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "global.h"

static sap_parms_t sap={{0,0,0,0},0,0,0};


static void check_block_size(int *bs)
{
   int n0,n1,n2,n3;
   
   error_root((bs[0]<4)||(bs[1]<4)||(bs[2]<4)||(bs[3]<4)||
              (bs[0]>L0)||(bs[1]>L1)||(bs[2]>L2)||(bs[3]>L3),1,
              "check_block_size [sap_parms.c]",
              "Block sizes are out of range");

   error_root((bs[0]%2)||(bs[1]%2)||(bs[2]%2)||(bs[3]%2),1,
              "check_block_size [sap_parms.c]",
              "Block sizes must be even");
   
   error_root((L0%bs[0])||(L1%bs[1])||(L2%bs[2])||(L3%bs[3]),1,
              "check_block_size [sap_parms.c]",
              "Blocks do not divide the local lattice");

   n0=L0/bs[0];
   n1=L1/bs[1];
   n2=L2/bs[2];
   n3=L3/bs[3];

   error_root(((NPROC0*n0)%2)||((NPROC1*n1)%2)||
              ((NPROC2*n2)%2)||((NPROC3*n3)%2),1,
              "check_block_size [sap_parms.c]",
              "There must be an even number of blocks in each direction");

   error_root((n0*n1*n2*n3)%2,1,
              "check_block_size [sap_parms.c]",
              "The number of blocks in the local lattice must be even");
}



sap_parms_t set_sap_parms(int *bs,int isolv,int nmr,int ncy)
{
   int iprms[7];

   if (NPROC>1)
   {
      iprms[0]=bs[0];
      iprms[1]=bs[1];
      iprms[2]=bs[2];
      iprms[3]=bs[3];
      iprms[4]=isolv;
      iprms[5]=nmr;
      iprms[6]=ncy;

      MPI_Bcast(iprms,7,MPI_INT,0,MPI_COMM_WORLD);

      error((iprms[0]!=bs[0])||(iprms[1]!=bs[1])||(iprms[2]!=bs[2])||
            (iprms[3]!=bs[3])||(iprms[4]!=isolv)||(iprms[5]!=nmr)||
            (iprms[6]!=ncy),1,
            "set_sap_parms [sap_parms.c]","Parameters are not global");
   }

   if (sap.ncy>0)
   {
      error_root((bs[0]!=sap.bs[0])||(bs[1]!=sap.bs[1])||
                 (bs[2]!=sap.bs[2])||(bs[3]!=sap.bs[3]),1,
                 "set_sap_parms [sap_parms.c]","bs[4] may be set only once");
   }
   else
   {      
      check_block_size(bs);
      sap.bs[0]=bs[0];
      sap.bs[1]=bs[1];
      sap.bs[2]=bs[2];
      sap.bs[3]=bs[3];
   }

   error_root((isolv<0)||(isolv>1)||(nmr<1)||(ncy<1),1,
              "set_sap_parms [sap_parms.c]",
              "Improper value of isolv, nmr or ncy");
      
   sap.isolv=isolv;
   sap.nmr=nmr;
   sap.ncy=ncy;
   
   return sap;
}


sap_parms_t sap_parms(void)
{
   return sap;
}


void print_sap_parms(int ipr)
{
   int my_rank;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      if (ipr)
      {
         printf("SAP parameters:\n");
         printf("bs = %d %d %d %d\n",
                sap.bs[0],sap.bs[1],sap.bs[2],sap.bs[3]);
         printf("isolv = %d\n",sap.isolv);
         printf("nmr = %d\n",sap.nmr);
         printf("ncy = %d\n\n",sap.ncy);         
      }
      else
      {
         printf("SAP block size:\n");
         printf("bs = %d %d %d %d\n\n",
                sap.bs[0],sap.bs[1],sap.bs[2],sap.bs[3]);
      }
   }
}


void write_sap_parms(FILE *fdat)
{
   int my_rank,endian;
   int i,iw;
   stdint_t istd[7];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   
   if (my_rank==0)
   {
      for (i=0;i<4;i++)
         istd[i]=(stdint_t)(sap.bs[i]);
      
      istd[4]=(stdint_t)(sap.isolv);
      istd[5]=(stdint_t)(sap.nmr);
      istd[6]=(stdint_t)(sap.ncy);
      
      if (endian==BIG_ENDIAN)
         bswap_int(7,istd);

      iw=fwrite(istd,sizeof(stdint_t),7,fdat);         
      error_root(iw!=7,1,"write_sap_parms [sap_parms.c]",
                 "Incorrect write count");
   }
}


void check_sap_parms(FILE *fdat)
{
   int my_rank,endian;
   int i,ir,ie;
   stdint_t istd[7];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   
   if (my_rank==0)
   {
      ir=fread(istd,sizeof(stdint_t),7,fdat);         
      error_root(ir!=7,1,"check_sap_parms [sap_parms.c]",
                 "Incorrect read count");         

      if (endian==BIG_ENDIAN)
         bswap_int(7,istd);

      ie=0;

      for (i=0;i<4;i++)
         ie|=(istd[i]!=(stdint_t)(sap.bs[i]));

      ie|=(istd[4]!=(stdint_t)(sap.isolv));
      ie|=(istd[5]!=(stdint_t)(sap.nmr));
      ie|=(istd[6]!=(stdint_t)(sap.ncy));
         
      error_root(ie!=0,1,"check_sap_parms [sap_parms.c]",
                 "Parameters do not match");
   }
}
