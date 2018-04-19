
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Writes the state of ranlxs to a file together with the next 147 random
* numbers. Then reads the data back in and checks the correct reinitialization
* of the generator
*
*******************************************************************************/

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "utils.h"

#define N 147


static void write_state(void)
{
   FILE *fp;
   int k,ns,*state;
   float base,r[N];

   ns=rlxs_size();
   state=malloc(ns*sizeof(int));
   base=(float)(ldexp(1.0,24));
   rlxs_init(1,1234567);

   for (k=0;k<10;k++) 
      ranlxs(r,N);

   rlxs_get(state);
   ranlxs(r,N);

   fp=fopen(".tmp","w");

   for (k=0;k<ns;k++)
      fprintf(fp,"%d\n",state[k]);

   for (k=0;k<N;k++)
      fprintf(fp,"%12.1f\n",base*r[k]);

   fclose(fp);
}


static void read_state(void)
{
   int ir,k,ns,*state;
   float base,r[N],r_old[N];
   FILE *fp;
   
   ns=rlxs_size();
   state=malloc(ns*sizeof(int));   
   base=(float)(ldexp(1.0,24));

   fp=fopen(".tmp","r");
   error(fp==NULL,1,"read_state [check2.c]","Unable to open file");
   ir=0;

   for (k=0;k<ns;k++)
      ir+=fscanf(fp,"%d",&state[k]);

   for (k=0;k<N;k++)
      ir+=fscanf(fp,"%f",&r_old[k]);

   error(ir!=(ns+N),1,"read_state [check2.c]","Read error");
   fclose(fp);
   remove(".tmp");
   
   rlxs_reset(state);
   ranlxs(r,N);
   
   for (k=0;k<N;k++)
   {
      if (r_old[k]!=(base*r[k]))
      {
         printf("\n");
         printf("Error: state of ranlxs has not been properly reset\n");
         printf("\n");
         exit(1);
      }
   }
}


int main(void)
{
   write_state();
   read_state();
   
   printf("\n");
   printf("State properly reset, generated random numbers are correct\n\n");
   exit(0);   
}
