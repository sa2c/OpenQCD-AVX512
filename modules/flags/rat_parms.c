
/*******************************************************************************
*
* File rat_parms.c
*
* Copyright (C) 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Rational function parameter data base
*
* The externally accessible functions are
*
*   rat_parms_t set_rat_parms(int irp,int degree,double *range)
*     Sets the parameters in the rational function parameter set number
*     irp and returns a structure containing them (see the notes).
*
*   rat_parms_t rat_parms(int irp)
*     Returns a structure containing the rational function parameter set
*     number irp (see the notes).
*
*   void read_rat_parms(int irp)
*     On process 0, this program scans stdin for a line starting with the
*     string "[Rational <int>]" (after any number of blanks), where <int> is
*     the integer value passed by the argument. An error occurs if no such
*     line or more than one is found. The lines 
*
*       degree  <int>
*       range   <double> <double>
*
*     are then read using read_line() [utils/mutils.c] and the data are 
*     entered into the data base by calling set_rat_parms().
*
*   void print_rat_parms(void)
*     Prints the defined rational function parameter sets to stdout on MPI
*     process 0.
*
*   void write_rat_parms(FILE *fdat)
*     Writes the defined rational function parameter sets to the file fdat 
*     on MPI process 0.
*
*   void check_rat_parms(FILE *fdat)
*     Compares the defined rational function parameter sets with those 
*     on the file fdat on MPI process 0, assuming the latter were written
*     to the file by the program write_rat_parms().
*
* Notes:
*
* Currently only Zolotorev rational functions are supported (see the modules
* ratfcts/zolotarev.c and ratfcts/ratfcts.c). The elements of a structure of
* type rat_parms_t are
*
*   degree       Degree of the rational function
*
*   range[2]     Lower and upper end of the approximation range (see
*                ratfcts/ratfcts.c)
*
* Up to 32 parameter sets, labeled by an index irp=0,1,..,31, can be
* specified. Once a set is defined, it cannot be changed by calling
* set_rat_parms() again. Rational function parameters must be globally
* the same.
*
* Except for rat_parms(), the programs in this module perform global
* operations and must be called simultaneously on all MPI processes.
*
*******************************************************************************/

#define RAT_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "global.h"

#define IRPMAX 32

static int init=0;
static rat_parms_t rp[IRPMAX+1]={{0,{0.0,0.0}}};


static void init_rp(void)
{
   int irp;
   
   for (irp=1;irp<=IRPMAX;irp++)
      rp[irp]=rp[0];

   init=1;
}


rat_parms_t set_rat_parms(int irp,int degree,double *range)
{
   int ie,iprms[2];
   double dprms[2];
   
   if (init==0)
      init_rp();

   if (NPROC>1)
   {
      iprms[0]=irp;
      iprms[1]=degree;
      dprms[0]=range[0];
      dprms[1]=range[1];
      
      MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

      ie=0;
      ie|=(iprms[0]!=irp);
      ie|=(iprms[1]!=degree);
      ie|=(dprms[0]!=range[0]);
      ie|=(dprms[1]!=range[1]);      
      
      error(ie!=0,1,"set_rat_parms [rat_parms.c]",
            "Parameters are not global");
   }

   ie=0;
   ie|=((irp<0)||(irp>=IRPMAX));
   ie|=(degree<1);
   ie|=(range[0]>=range[1]);
   ie|=(range[0]<=0.0);

   error_root(ie!=0,1,"set_rat_parms [rat_parms.c]",
              "Parameters are out of range");
   
   error_root(rp[irp].degree!=0,1,"set_rat_parms [rat_parms.c]",
              "Attempt to reset an already specified parameter set");

   rp[irp].degree=degree;
   rp[irp].range[0]=range[0];
   rp[irp].range[1]=range[1];

   return rp[irp];
}


rat_parms_t rat_parms(int irp)
{
   if (init==0)
      init_rp();

   if ((irp>=0)&&(irp<IRPMAX))
      return rp[irp];
   else
   {
      error_loc(1,1,"rat_parms [rat_parms.c]",
                "Rational function index is out of range");
      return rp[IRPMAX];
   }
}


void read_rat_parms(int irp)
{
   int my_rank;
   int degree;
   double range[2];
   char line[NAME_SIZE];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   
   if (my_rank==0)
   {
      sprintf(line,"Rational %d",irp);
      find_section(line);

      read_line("degree","%d",&degree);
      read_line("range","%lf %lf",range,range+1);      
   }

   if (NPROC>1)
   {
      MPI_Bcast(&degree,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(range,2,MPI_DOUBLE,0,MPI_COMM_WORLD);   
   }
   
   set_rat_parms(irp,degree,range);
}


void print_rat_parms(void)
{
   int my_rank,irp,n[2];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   
   if ((my_rank==0)&&(init==1))
   {
      for (irp=0;irp<IRPMAX;irp++)
      {
         if (rp[irp].degree!=0)
         {
            printf("Rational %d:\n",irp);
            printf("degree = %d\n",rp[irp].degree);
            n[0]=fdigits(rp[irp].range[0]);
            n[1]=fdigits(rp[irp].range[1]);            
            printf("range = [%.*f,%.*f]\n\n",IMAX(n[0],1),rp[irp].range[0],
                   IMAX(n[1],1),rp[irp].range[1]);
         }
      }
   }
}


void write_rat_parms(FILE *fdat)
{
   int my_rank,endian;
   int iw,irp;
   stdint_t istd[2];
   double dstd[2];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   
   if ((my_rank==0)&&(init==1))
   {
      for (irp=0;irp<IRPMAX;irp++)
      {
         if (rp[irp].degree!=0)
         {
            istd[0]=(stdint_t)(irp);            
            istd[1]=(stdint_t)(rp[irp].degree);
            dstd[0]=rp[irp].range[0];
            dstd[1]=rp[irp].range[1];
            
            if (endian==BIG_ENDIAN)
            {
               bswap_int(2,istd);
               bswap_double(2,dstd);
            }
            
            iw=fwrite(istd,sizeof(stdint_t),2,fdat);
            iw+=fwrite(dstd,sizeof(double),2,fdat);

            error_root(iw!=4,1,"write_rat_parms [rat_parms.c]",
                       "Incorrect write count");
         }
      }
   }
}


void check_rat_parms(FILE *fdat)
{
   int my_rank,endian;
   int ir,irp,ie;
   stdint_t istd[2];
   double dstd[2];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
      
   if ((my_rank==0)&&(init==1))
   {
      ie=0;
      
      for (irp=0;irp<IRPMAX;irp++)
      {
         if (rp[irp].degree!=0)
         {
            ir=fread(istd,sizeof(stdint_t),2,fdat);
            ir+=fread(dstd,sizeof(double),2,fdat);

            error_root(ir!=4,1,"check_rat_parms [rat_parms.c]",
                       "Incorrect read count");

            if (endian==BIG_ENDIAN)
            {
               bswap_int(2,istd);
               bswap_double(2,dstd);
            }
            
            ie|=(istd[0]!=(stdint_t)(irp));            
            ie|=(istd[1]!=(stdint_t)(rp[irp].degree));
            ie|=(dstd[0]!=rp[irp].range[0]);
            ie|=(dstd[1]!=rp[irp].range[1]);
         }
      }
         
      error_root(ie!=0,1,"check_rat_parms [rat_parms.c]",
                 "Parameters do not match");         
   }
}
