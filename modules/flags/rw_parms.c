
/*******************************************************************************
*
* File rw_parms.c
*
* Copyright (C) 2012-2014 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Reweighting factor parameter data base.
*
* The externally accessible functions are
*
*   rw_parms_t set_rw_parms(int irw,rwfact_t rwfact,int im0,int nsrc,
*                           int irp,int nfct,double *mu,int *np,int *isp)
*     Sets the parameters in the reweighting factor parameter set number
*     irw and returns a structure containing them (see the notes).
*
*   rw_parms_t rw_parms(int irw)
*     Returns a structure containing the reweighting factor parameter set
*     number irw (see the notes).
*
*   void read_rw_parms(int irw)
*     On process 0, this program scans stdin for a line starting with the
*     string "[Reweighting factor <int>]" (after any number of blanks), where
*     <int> is the integer value passed through the argument. An error occurs
*     if no such line or more than one is found. The lines
*
*       rwfact   <rwfact_t>
*       im0      <int>
*       nsrc     <int>
*       irp      <int>
*       mu       <double> [<double>]
*       np       <int> [<int>]
*       isp      <int> [<int>]
*
*     are then read using read_line() [utils/mutils.c] and the data are
*     added to the data base by calling set_rw_parms(irw,...). Depending
*     on the value of "rwfact", some lines are not read and can be omitted
*     in the input file. The number of items on the lines with tag "mu",
*     "np" and "isp" depends on the reweighting factor too (see the notes).
*
*   void print_rw_parms(void)
*     Prints the defined reweighting factor parameter sets to stdout on
*     MPI process 0.
*
*   void write_rw_parms(FILE *fdat)
*     Writes the defined reweighting factor parameter sets to the file fdat
*     on MPI process 0.
*
*   void check_rw_parms(FILE *fdat)
*     Compares the defined reweighting factor parameter sets with those
*     on the file fdat on MPI process 0, assuming the latter were written
*     to the file by the program write_rw_parms().
*
* Notes:
*
* The elements of a structure of type rw_parms_t are:
*
*   rwfact  Reweighting factor program used. This parameter is an enum
*           type with one of the following values:
*
*            RWTM1       (program rwtm1() [update/rwtm.c]),
*
*            RWTM1_EO    (program rwtm1eo() [update/rwtmeo.c]),
*
*            RWTM2       (program rwtm2() [update/rwtm.c]),
*
*            RWTM2_EO    (program rwtm2eo() [update/rwtmeo.c]),
*
*            RWRAT       (program rwrat() [update/rwrat.c]).
*
*   im0     Index of the bare sea quark mass in the parameter data base
*           (see flags/lat_parms.c).
*
*   nsrc    Number N of random source fields to be used for the stochastic
*           estimation of the reweighting factor. If the latter is split
*           into a product factors, N random fields are used for each of
*           them.
*
*   irp     Rational function parameter set index. Only relevant if
*           rwfact=RWRAT.
*
*   nfct    If rwfact=RWTM*: Number of Hasenbusch factors into which the
*           reweighting factor is decomposed;
*           If rwfact=RWRAT: Number of rational factors into which the
*           rational function is decomposed.
*
*   mu      Array of twisted masses that define the Hasenbusch factors
*           (nfct elements; 0<mu[0]<mu[1]<..<mu[nfct-1]). Only relevant
*           if rwfact=RWTM* and otherwise set to NULL.
*
*   np      Array of the numbers of poles of the rational factors into
*           which the rational function is decomposed. Only relevant if
*           rwfact=RWRAT and otherwise set to NULL.
*
*   isp     Array of solver parameter set indices describing the solvers
*           for the Dirac equation to be used. The array must have nfct
*           elements that correspond to the factors of the reweighting
*           factor if rwfact=RWTM* or, if rwfact=RWRAT, to the rational
*           functions into which the rational function is decomposed.
*
* Valid examples of parameter sections that can be read by read_rw_parms()
* are
*
*   [Reweighting factor 1]
*    rwfact   RWTM1           # Or RWTM1_EO, RWTM2, RWTM2_EO
*    im0      0
*    nsrc     12
*    mu       0.001 0.003     # Implies a decomposition in 2 factors
*    isp      3               # Solver no 3 will be used for all factors
*
*   [Reweighting factor 4]
*    rwfact   RWRAT
*    im0      1
*    nsrc     4
*    irp      0
*    np       6 2 2           # Implies a decomposition in 3 rational parts
*    isp      3 5 6           # For each part, a separate solver is used
*
* The number of solver parameter sets specified on the line with tag "isp"
* need not match the number of factors or rational parts specified on the
* lines with tag "mu" or "np". If fewer are given, the list is padded with
* the last entry on the line. Superfluos entries are ignored.
*
* Up to 32 reweighting parameter sets, labeled by an index irw=0,1,..,31,
* can be specified. Once a set is defined, it cannot be changed by calling
* set_rw_parms() again. All parameters must be globally the same.
*
* Except for rw_parms(), the programs in this module perform global operations
* and must be called simultaneously on all MPI processes.
*
*******************************************************************************/

#define RW_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "global.h"

#define IRWMAX 32

static int init=0;
static rwfact_t rwfact[]={RWTM1,RWTM1_EO,RWTM2,RWTM2_EO,RWRAT};
static rw_parms_t rw[IRWMAX+1]={{RWFACTS,0,0,0,0,NULL,NULL,NULL}};


static void init_rw(void)
{
   int irw;

   for (irw=1;irw<=IRWMAX;irw++)
      rw[irw]=rw[0];

   init=1;
}


rw_parms_t set_rw_parms(int irw,rwfact_t rwfact,int im0,int nsrc,
                        int irp,int nfct,double *mu,int *np,int *isp)
{
   int iprms[6],i,ie;
   double dprms[1];

   if (init==0)
      init_rw();

   error_root((rwfact!=RWTM1)&&(rwfact!=RWTM1_EO)&&
              (rwfact!=RWTM2)&&(rwfact!=RWTM2_EO)&&(rwfact!=RWRAT),1,
              "set_rw_parms [rw_parms.c]","Unknown type of reweighting factor");

   if (rwfact!=RWRAT)
      irp=0;

   if (NPROC>1)
   {
      iprms[0]=irw;
      iprms[1]=(int)(rwfact);
      iprms[2]=im0;
      iprms[3]=nsrc;
      iprms[4]=irp;
      iprms[5]=nfct;

      MPI_Bcast(iprms,6,MPI_INT,0,MPI_COMM_WORLD);

      ie=0;
      ie|=(iprms[0]!=irw);
      ie|=(iprms[1]!=(int)(rwfact));
      ie|=(iprms[2]!=im0);
      ie|=(iprms[3]!=nsrc);
      ie|=(iprms[4]!=irp);
      ie|=(iprms[5]!=nfct);

      error(ie!=0,1,"set_rw_parms [rw_parms.c]",
            "Parameters are not global");
   }

   ie=0;
   ie|=((irw<0)||(irw>=IRWMAX));
   ie|=(im0<0);
   ie|=(nsrc<1);
   ie|=(irp<0);
   ie|=(nfct<1);

   error_root(ie!=0,1,"set_rw_parms [rw_parms.c]",
              "Parameters are out of range");

   if (NPROC>1)
   {
      if (rwfact!=RWRAT)
      {
         for (i=0;i<nfct;i++)
         {
            dprms[0]=mu[i];
            iprms[0]=isp[i];

            MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

            ie|=(dprms[0]!=mu[i]);
            ie|=(iprms[0]!=isp[i]);
         }

         error(ie!=0,1,"set_rw_parms [rw_parms.c]",
               "Parameters mu or isp are not global");
      }
      else
      {
         for (i=0;i<nfct;i++)
         {
            iprms[0]=np[i];
            iprms[1]=isp[i];

            MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

            ie|=(iprms[0]!=np[i]);
            ie|=(iprms[1]!=isp[i]);
         }

         error(ie!=0,1,"set_rw_parms [rw_parms.c]",
               "Parameters np or isp are not global");
      }
   }

   error_root(rw[irw].rwfact!=RWFACTS,1,"set_rw_parms [rw_parms.c]",
              "Attempt to reset an already specified parameter set");

   rw[irw].rwfact=rwfact;
   rw[irw].im0=im0;
   rw[irw].nsrc=nsrc;
   rw[irw].irp=irp;
   rw[irw].nfct=nfct;

   if (rwfact!=RWRAT)
   {
      rw[irw].mu=malloc(nfct*sizeof(*mu));
      rw[irw].np=NULL;
      rw[irw].isp=malloc(nfct*sizeof(*isp));

      error_root((rw[irw].mu==NULL)||(rw[irw].isp==NULL),1,
                 "set_rw_parms [rw_parms.c]",
                 "Unable to allocate parameter arrays");

      for (i=0;i<nfct;i++)
      {
         rw[irw].mu[i]=mu[i];
         rw[irw].isp[i]=isp[i];

         if (i==0)
            ie|=(mu[i]<=0.0);
         else
            ie|=(mu[i]<=mu[i-1]);
      }

      error_root(ie!=0,1,"set_rw_parms [rw_parms.c]",
                 "The twisted masses must be in ascending order");
   }
   else
   {
      rw[irw].np=malloc(2*nfct*sizeof(*np));
      rw[irw].isp=rw[irw].np+nfct;
      rw[irw].mu=NULL;

      error_root(rw[irw].np==NULL,1,"set_rw_parms [rw_parms.c]",
                 "Unable to allocate parameter arrays");

      for (i=0;i<nfct;i++)
      {
         rw[irw].np[i]=np[i];
         rw[irw].isp[i]=isp[i];
      }
   }

   return rw[irw];
}


rw_parms_t rw_parms(int irw)
{
   if (init==0)
      init_rw();

   if ((irw>=0)&&(irw<IRWMAX))
      return rw[irw];
   else
   {
      error_loc(1,1,"rw_parms [rw_parms.c]",
                "Reweighting factor index is out of range");
      return rw[IRWMAX];
   }
}


void read_rw_parms(int irw)
{
   int my_rank,n,i;
   int idr,im0,nsrc,irp,nfct;
   int *np,*isp;
   double *mu;
   char line[NAME_SIZE];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      sprintf(line,"Reweighting factor %d",irw);
      find_section(line);

      read_line("rwfact","%s",line);

      if (strcmp(line,"RWTM1")==0)
         idr=0;
      else if (strcmp(line,"RWTM1_EO")==0)
         idr=1;
      else if (strcmp(line,"RWTM2")==0)
         idr=2;
      else if (strcmp(line,"RWTM2_EO")==0)
         idr=3;
      else if (strcmp(line,"RWRAT")==0)
         idr=4;
      else
      {
         idr=5;
         error_root(1,1,"read_rw_parms [rw_parms.c]",
                    "Unknown reweighting factor %s",line);
      }

      read_line("im0","%d",&im0);
      read_line("nsrc","%d",&nsrc);

      if (idr<4)
      {
         irp=0;
         nfct=count_tokens("mu");
         error_root(nfct<1,1,"read_rw_parms [rw_parms.c]",
                    "No data on line with tag mu");
      }
      else
      {
         read_line("irp","%d",&irp);
         nfct=count_tokens("np");
         error_root(nfct<1,1,"read_rw_parms [rw_parms.c]",
                    "No data on line with tag np");
      }
   }
   else
   {
      idr=0;
      im0=0;
      nsrc=0;
      irp=0;
      nfct=0;
   }

   if (NPROC>1)
   {
      MPI_Bcast(&idr,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&im0,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nsrc,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&irp,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nfct,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   if (idr<4)
   {
      mu=malloc(nfct*sizeof(*mu));
      np=NULL;
      isp=malloc(nfct*sizeof(*isp));
      error((mu==NULL)||(isp==NULL),1,"read_rw_parms [rw_parms.c]",
            "Unable to allocated data arrays");
   }
   else
   {
      mu=NULL;
      np=malloc(2*nfct*sizeof(*np));
      isp=np+nfct;
      error(np==NULL,1,"read_rw_parms [rw_parms.c]",
            "Unable to allocated data arrays");
   }

   if (my_rank==0)
   {
      if (idr<4)
         read_dprms("mu",nfct,mu);
      else
         read_iprms("np",nfct,np);

      n=count_tokens("isp");
      error_root(n<1,1,"read_rw_parms [rw_parms.c]",
            "No data on the line with tag isp");

      if (n>nfct)
         n=nfct;
      read_iprms("isp",n,isp);

      for (i=n;i<nfct;i++)
         isp[i]=isp[n-1];
   }

   if (NPROC>1)
   {
      if (idr<4)
         MPI_Bcast(mu,nfct,MPI_DOUBLE,0,MPI_COMM_WORLD);
      else
         MPI_Bcast(np,nfct,MPI_INT,0,MPI_COMM_WORLD);

      MPI_Bcast(isp,nfct,MPI_INT,0,MPI_COMM_WORLD);
   }

   set_rw_parms(irw,rwfact[idr],im0,nsrc,irp,nfct,mu,np,isp);

   if (idr<4)
   {
      free(mu);
      free(isp);
   }
   else
      free(np);
}


void print_rw_parms(void)
{
   int my_rank,irw,idr,nfct,n,i;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if ((my_rank==0)&&(init==1))
   {
      for (irw=0;irw<IRWMAX;irw++)
      {
         if (rw[irw].rwfact!=RWFACTS)
         {
            printf("Reweighting factor %d:\n",irw);
            idr=0;

            if (rw[irw].rwfact==RWTM1)
               printf("RWTM1 factor\n");
            else if (rw[irw].rwfact==RWTM1_EO)
               printf("RWTM1_EO factor\n");
            else if (rw[irw].rwfact==RWTM2)
               printf("RWTM2 factor\n");
            else if (rw[irw].rwfact==RWTM2_EO)
               printf("RWTM2_EO factor\n");
            else if (rw[irw].rwfact==RWRAT)
            {
               idr=1;
               printf("RWRAT factor\n");
            }

            printf("im0 = %d\n",rw[irw].im0);
            printf("nsrc = %d\n",rw[irw].nsrc);
            nfct=rw[irw].nfct;

            if (idr==0)
            {
               printf("mu =");

               for (i=0;i<nfct;i++)
               {
                  n=fdigits(rw[irw].mu[i]);
                  printf(" %.*f",IMAX(n,1),rw[irw].mu[i]);
               }

               printf("\n");
            }
            else
            {
               printf("irp = %d\n",rw[irw].irp);
               printf("np =");

               for (i=0;i<nfct;i++)
                  printf(" %d",rw[irw].np[i]);

               printf("\n");
            }

            printf("isp =");

            for (i=0;i<nfct;i++)
               printf(" %d",rw[irw].isp[i]);

            printf("\n\n");
         }
      }
   }
}


void write_rw_parms(FILE *fdat)
{
   int my_rank,endian;
   int iw,irw,nfct,i;
   stdint_t istd[6];
   double dstd[1];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();

   if ((my_rank==0)&&(init==1))
   {
      for (irw=0;irw<IRWMAX;irw++)
      {
         if (rw[irw].rwfact!=RWFACTS)
         {
            istd[0]=(stdint_t)(irw);
            istd[1]=(stdint_t)(rw[irw].rwfact);
            istd[2]=(stdint_t)(rw[irw].im0);
            istd[3]=(stdint_t)(rw[irw].nsrc);
            istd[4]=(stdint_t)(rw[irw].irp);
            istd[5]=(stdint_t)(rw[irw].nfct);

            if (endian==BIG_ENDIAN)
               bswap_int(6,istd);

            iw=fwrite(istd,sizeof(stdint_t),6,fdat);
            nfct=rw[irw].nfct;

            if (rw[irw].rwfact==RWRAT)
            {
               for (i=0;i<nfct;i++)
               {
                  istd[0]=(stdint_t)(rw[irw].np[i]);

                  if (endian==BIG_ENDIAN)
                     bswap_int(1,istd);

                  iw+=fwrite(istd,sizeof(stdint_t),1,fdat);
               }
            }
            else
            {
               for (i=0;i<nfct;i++)
               {
                  dstd[0]=rw[irw].mu[i];

                  if (endian==BIG_ENDIAN)
                     bswap_double(1,dstd);

                  iw+=fwrite(dstd,sizeof(double),1,fdat);
               }
            }

            for (i=0;i<nfct;i++)
            {
               istd[0]=(stdint_t)(rw[irw].isp[i]);

               if (endian==BIG_ENDIAN)
                  bswap_int(1,istd);

               iw+=fwrite(istd,sizeof(stdint_t),1,fdat);
            }

            error_root(iw!=(6+2*nfct),1,"write_rw_parms [rw_parms.c]",
                       "Incorrect write count");
         }
      }
   }
}


void check_rw_parms(FILE *fdat)
{
   int my_rank,endian;
   int ir,irw,nfct,i,ie;
   stdint_t istd[6];
   double dstd[1];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();

   if ((my_rank==0)&&(init==1))
   {
      ie=0;

      for (irw=0;irw<IRWMAX;irw++)
      {
         if (rw[irw].rwfact!=RWFACTS)
         {
            ir=fread(istd,sizeof(stdint_t),6,fdat);

            error_root(ir!=6,1,"check_rw_parms [rw_parms.c]",
                       "Incorrect read count");

            if (endian==BIG_ENDIAN)
               bswap_int(6,istd);

            ie|=(istd[0]!=(stdint_t)(irw));
            ie|=(istd[1]!=(stdint_t)(rw[irw].rwfact));
            ie|=(istd[2]!=(stdint_t)(rw[irw].im0));
            ie|=(istd[3]!=(stdint_t)(rw[irw].nsrc));
            ie|=(istd[4]!=(stdint_t)(rw[irw].irp));
            ie|=(istd[5]!=(stdint_t)(rw[irw].nfct));

            error_root(ie!=0,1,"check_rw_parms [rw_parms.c]",
                       "Parameters do not match");

            nfct=rw[irw].nfct;

            if (rw[irw].rwfact==RWRAT)
            {
               for (i=0;i<nfct;i++)
               {
                  ir+=fread(istd,sizeof(stdint_t),1,fdat);

                  if (endian==BIG_ENDIAN)
                     bswap_int(1,istd);

                  ie|=(istd[0]!=(stdint_t)(rw[irw].np[i]));
               }
            }
            else
            {
               for (i=0;i<nfct;i++)
               {
                  ir+=fread(dstd,sizeof(double),1,fdat);

                  if (endian==BIG_ENDIAN)
                     bswap_double(1,dstd);

                  ie|=(dstd[0]!=rw[irw].mu[i]);
               }
            }

            for (i=0;i<nfct;i++)
            {
               ir+=fread(istd,sizeof(stdint_t),1,fdat);

               if (endian==BIG_ENDIAN)
                  bswap_int(1,istd);

               ie|=(istd[0]!=(stdint_t)(rw[irw].isp[i]));
            }

            error_root(ir!=(6+2*nfct),1,"check_rw_parms [rw_parms.c]",
                       "Incorrect read count");
            error_root(ie!=0,1,"check_rw_parms [rw_parms.c]",
                       "Parameters do not match");
         }
      }
   }
}
