
/*******************************************************************************
*
* File hsum.c
*
* Copyright (C) 2015, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Hierarchical summation programs.
*
*   int init_hsum(int n)
*     Creates a new instance of a hierarchical sum and returns the id of
*     the instance. The parameter n>=1 specifies the number of streams
*     of double-precision floating point numbers that are to be summed
*     in parallel.
*
*   void reset_hsum(int id)
*     Resets the hiearchical sum with the given id to zero.
*
*   void add_to_hsum(int id,double *x)
*     Adds the numbers x[0],..,x[n-1] to the hiearchical sum with the
*     given id. More precisely, the number x[k] is added to the k'th
*     sum accumulated in parallel. The number n of parallel sums is
*     the one set by init_hsum().
*
*   void local_hsum(int id,double *sx)
*     Evaluates the hiearchical sum with the given id and assigns the
*     k'th sum computed in parallel to sx[k] (k=0,..,n-1). The number
*     n of parallel sums is the one set by init_hsum().
*
*   void global_hsum(int id,double *sx)
*     Same as local_hsum(), but additionally sums the results over all
*     MPI processes. The calculated sums are guaranteed to be exactly
*     the same on all processes.
*
* Except for global_hsum(), all programs in this module can be locally
* called.
*
* While the argument id of the program global_hsum() need not be globally
* the same, the number of parallel sums processed by the hiearchical sums
* with the specified id's may not vary from process to process.
*
* The maximal number of hsum instances is equal to INT_MAX. Any attempt to
* initialize more than INT_MAX hsums generates an error.
*
*******************************************************************************/

#define HSUM_C

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "mpi.h"
#include "utils.h"
#include "global.h"

#define BLK_LENGTH 8
#define MAX_LEVELS 12

typedef struct
{
   int n;
   int cnt[MAX_LEVELS];
   double *sm[MAX_LEVELS];
} hsum_t;

static int nhs=0,nhsmx=0,nmx=0;
static double *sloc=NULL;
static hsum_t **hs=NULL;


static hsum_t *alloc_hsum(int n)
{
   int i;
   double *s;
   hsum_t *h;

   h=malloc(1*sizeof(*h));
   s=malloc(MAX_LEVELS*n*sizeof(*s));

   if (n>nmx)
   {
      if (sloc!=NULL)
         free(sloc);
      sloc=malloc(n*sizeof(*sloc));
      nmx=n;
   }

   error_loc((h==NULL)||(s==NULL)||(sloc==NULL),1,"alloc_hsum [hsum.c]",
             "Unable to allocate hsum_t structure");

   (*h).n=n;

   for (i=0;i<MAX_LEVELS;i++)
   {
      (*h).sm[i]=s;
      s+=n;
   }

   return h;
}


int init_hsum(int n)
{
   int i,mx;
   hsum_t **hsn;

   error_loc(n<1,1,"init_hsum [hsum.c]","Parameter n is out of range");

   if (nhs==nhsmx)
   {
      error_loc(nhsmx==INT_MAX,1,"init_hsum [hsum.c]",
                "Maximal number of hsum instances reached");

      if (nhsmx<(INT_MAX-16))
         mx=nhsmx+16;
      else
         mx=INT_MAX;

      hsn=malloc(mx*sizeof(*hsn));
      error_loc(hsn==NULL,1,"init_hsum [hsum.c]",
                "Unable to allocate hsum_t structures");

      for (i=0;i<nhsmx;i++)
         hsn[i]=hs[i];

      for (;i<mx;i++)
         hsn[i]=NULL;

      if (hs!=NULL)
         free(hs);
      hs=hsn;
      nhsmx=mx;
   }

   hs[nhs]=alloc_hsum(n);
   nhs+=1;

   return nhs-1;
}


void reset_hsum(int id)
{
   int n,i,j,*cnt;
   double **sm;

   if ((id>=0)&&(id<nhs))
   {
      n=hs[id][0].n;
      cnt=hs[id][0].cnt;
      sm=hs[id][0].sm;

      for (i=0;i<MAX_LEVELS;i++)
      {
         cnt[i]=0;

         for (j=0;j<n;j++)
            sm[i][j]=0.0;
      }
   }
   else
      error_loc(1,1,"reset_hsum [hsum.c]","Parameter id is out of range");
}


void add_to_hsum(int id,double *x)
{
   int n,i,j,*cnt;
   double **sm;

   if ((id>=0)&&(id<nhs))
   {
      n=hs[id][0].n;
      cnt=hs[id][0].cnt;
      sm=hs[id][0].sm;

      for (j=0;j<n;j++)
         sm[0][j]+=x[j];

      cnt[0]+=1;

      for (i=1;(cnt[i-1]>=BLK_LENGTH)&&(i<MAX_LEVELS);i++)
      {
         for (j=0;j<n;j++)
         {
            sm[i][j]+=sm[i-1][j];
            sm[i-1][j]=0.0;
         }

         cnt[i]+=1;
         cnt[i-1]=0;
      }
   }
   else
      error_loc(1,1,"add_to_hsum [hsum.c]","Parameter id is out of range");
}


void local_hsum(int id,double *sx)
{
   int n,i,j,*cnt;
   double **sm;

   if ((id>=0)&&(id<nhs))
   {
      n=hs[id][0].n;
      cnt=hs[id][0].cnt;
      sm=hs[id][0].sm;

      for (j=0;j<n;j++)
         sx[j]=sm[0][j];

      for (i=1;i<MAX_LEVELS;i++)
      {
         if (cnt[i]>0)
         {
            for (j=0;j<n;j++)
               sx[j]+=sm[i][j];
         }
      }
   }
   else
      error_loc(1,1,"local_hsum [hsum.c]","Parameter id is out of range");
}


void global_hsum(int id,double *sx)
{
   int n,iprms[1];

   if ((id>=0)&&(id<nhs))
      n=hs[id][0].n;
   else
      n=0;

   error_root(n==0,1,"global_hsum [hsum.c]","Parameter id is out of range");

   if (NPROC>1)
   {
      iprms[0]=n;
      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error(iprms[0]!=n,1,"global_hsum [hsum.c]",
            "Number of parallel sums is not global");

      local_hsum(id,sloc);
      MPI_Reduce(sloc,sx,n,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(sx,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
      local_hsum(id,sx);
}
