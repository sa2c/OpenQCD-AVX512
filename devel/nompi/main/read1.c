
/*******************************************************************************
*
* File read1.c
*
* Copyright (C) 2010-2014 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Reads and evaluates data from the *.dat files created by the programs qcd1
* and ym1. The file to be read has to be specified on the command line.
*
* This program writes the history of the MD energy deficit dH, the acceptance
* flag iac and the average plaquette to the file <run name>.run1.dat in the
* plots directory. In addition, some information about the distribution of dH
* and the integrated autocorrelation time of the plaquette are printed to
* stdout.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "extras.h"

typedef struct
{
   int nt,iac;
   double dH,avpl;
} dat_t;

static int nms,nfirst,nlast,neff;
static dat_t *adat;


static int read_dat(int n,dat_t *ndat,FILE *fin)
{
   int i,ir,ic,endian;
   stdint_t istd[2];
   double dstd[2];

   endian=endianness();
   ic=0;

   for (i=0;i<n;i++)
   {
      ir=fread(istd,sizeof(stdint_t),2,fin);
      ir+=fread(dstd,sizeof(double),2,fin);

      if (ir!=4)
         return ic;

      if (endian==BIG_ENDIAN)
      {
         bswap_int(2,istd);
         bswap_double(2,dstd);
      }

      (*ndat).nt=(int)(istd[0]);
      (*ndat).iac=(int)(istd[1]);

      (*ndat).dH=dstd[0];
      (*ndat).avpl=dstd[1];

      ic+=1;
      ndat+=1;
   }

   return ic;
}


static void read_file(char *fin)
{
   long ipos;
   dat_t ndat;
   FILE *fdat;

   fdat=fopen(fin,"rb");
   error(fdat==NULL,1,"read_file [read1.c]","Unable to open data file");

   printf("Read data from file %s\n\n",fin);
   ipos=ftell(fdat);
   nms=0;

   while (read_dat(1,&ndat,fdat)==1)
      nms+=1;

   error(nms==0,1,"read_file [read1.c]",
         "Empty data file");

   adat=amalloc(nms*sizeof(*adat),3);
   error(adat==NULL,1,"read_file [read1.c]",
         "Unable to allocate data array");

   fseek(fdat,ipos,SEEK_SET);
   error(read_dat(nms,adat,fdat)!=nms,1,"read_file [read1.c]",
         "Error while reading data file");
   fclose(fdat);
}


static void select_range(void)
{
   int n,no,nf,nl;
   int np,dn,ie;

   printf("There are %d measurements (trajectories no %d - %d).\n",
          nms,adat[0].nt,adat[nms-1].nt);
   printf("Range [nfirst,nlast] of trajectories to analyse: ");
   scanf("%d %d",&nfirst,&nlast);

   nf=0;
   nl=0;

   for (n=0;n<nms;n++)
   {
      no=adat[n].nt;

      if (no<nfirst)
         nf+=1;

      if (no<=nlast)
         nl+=1;
   }

   nfirst=nf;
   nlast=nl;
   neff=nlast-nfirst;

   printf("Keep %d measurements (trajectories no %d - %d).\n\n",
          neff,adat[nfirst].nt,adat[nlast-1].nt);

   error(neff<2,1,"select_range [read1.c]",
         "Selected range contains less than 2 measurements");

   np=adat[nfirst].nt;
   dn=adat[nfirst+1].nt-adat[nfirst].nt;

   if (dn<=0)
      ie=1;
   else
   {
      ie=0;

      for (n=(nfirst+1);n<nlast;n++)
      {
         no=adat[n].nt;

         if ((no-np)!=dn)
         {
            ie=2;
            break;
         }

         np=no;
      }
   }

   error(ie!=0,1,"select_range [read1.c]",
         "Varying trajectory number separation in selected range");
}


static double tail(int n,double *a,double amx)
{
   int ia,ic;

   ic=0;

   for (ia=0;ia<n;ia++)
      ic+=(fabs(a[ia])>amx);

   return (double)(ic)/(double)(n);
}


static double f(int nx,double x[])
{
   return x[0];
}


static void print_plot(char *fin)
{
   int n,ims;
   char base[NAME_SIZE],plt_file[NAME_SIZE],*p;
   dat_t *ndat;
   FILE *fout;

   p=strstr(fin,".dat");
   error(p==NULL,1,"print_plot [read1.c]","Unexpected data file name");
   n=p-fin;

   p=strrchr(fin,'/');
   if (p==NULL)
      p=fin;
   else
      p+=1;
   n-=(p-fin);

   error(n>=NAME_SIZE,1,"print_plot [read1.c]","File name is too long");
   strncpy(base,p,n);
   base[n]='\0';

   error(name_size("plots/%s.run1.dat",base)>=NAME_SIZE,1,
         "print_plot [read1.c]","File name is too long");
   sprintf(plt_file,"plots/%s.run1.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_plot [read1.c]",
         "Unable to open output file");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ym1 or qcd1\n");
   fprintf(fout,"# ---------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# Number of measurements = %d\n",nms);
   fprintf(fout,"#\n");
   fprintf(fout,"# nt:   trajectory number\n");
   fprintf(fout,"# dH:   MD energy deficit\n");
   fprintf(fout,"# iac:  acceptance flag\n");
   fprintf(fout,"#\n");
   fprintf(fout,"#  nt         dH      iac    <tr{U(p)}>\n");
   fprintf(fout,"#\n");

   ndat=adat;

   for (ims=0;ims<nms;ims++)
   {
      fprintf(fout," %5d  ",(*ndat).nt);
      fprintf(fout," % .4e  ",(*ndat).dH);
      fprintf(fout," %1d  ",(*ndat).iac);
      fprintf(fout," %.8e",(*ndat).avpl);
      fprintf(fout,"\n");

      ndat+=1;
   }

   fclose(fout);

   printf("Data printed to file %s\n\n",plt_file);
}


int main(int argc,char *argv[])
{
   int n;
   double *a,abar;

   error(argc!=2,1,"main [read1.c]","Syntax: read1 <filename>");

   printf("\n");
   printf("HMC simulation of QCD\n");
   printf("---------------------\n\n");

   read_file(argv[1]);
   select_range();

   a=malloc(neff*sizeof(double));
   error(a==NULL,1,"main [read1.c]",
         "Unable to allocate data array");

   for (n=0;n<neff;n++)
      a[n]=fabs(adat[nfirst+n].dH);

   printf("Fraction of trajectories with |dH| larger than\n\n");
   printf("    1.0: %.4f\n",tail(neff,a,1.0));
   printf("    2.0: %.4f\n",tail(neff,a,2.0));
   printf("   10.0: %.4f\n",tail(neff,a,10.0));
   printf("  100.0: %.4f\n",tail(neff,a,100.0));
   printf(" 1000.0: %.4f\n",tail(neff,a,1000.0));
   printf("\n");

   for (n=0;n<neff;n++)
      a[n]=exp(-adat[nfirst+n].dH);

   printf("<exp(-dH)> = %.3f (%.3f)\n",
          average(neff,a),sigma0(neff,a));

   for (n=0;n<neff;n++)
     {
       if (adat[nfirst+n].dH>0.0)
	 a[n]=exp(-adat[nfirst+n].dH);
       else
	 a[n]=1.0;
     }

   printf("<min{1,exp(-dH)}> = %.3f (%.3f)\n",
          average(neff,a),sigma0(neff,a));

   for (n=0;n<neff;n++)
      a[n]=(double)(adat[nfirst+n].iac);

   printf("<iac> = %.3f (%.3f)\n\n",
          average(neff,a),sigma0(neff,a));

   for (n=0;n<neff;n++)
      a[n]=adat[nfirst+n].avpl;

   printf("The integrated autocorrelation time and the associated"
          "\nstatistical error sigma of the plaquette is estimated ");

   if (neff>100)
      printf("using the\nnumerically determined "
             "autocorrelation function.\n\n");
   else
      printf("by binning the\ndata and by calculating "
             "the jackknife errors of the binned series.\n\n");

   printf("The autocorrelation times are given in numbers of measurements\n"
          "separated by %d trajectories.\n\n",
          adat[nfirst+1].nt-adat[nfirst].nt);

   if (neff>=100)
      abar=print_auto(neff,a);
   else
      abar=print_jack(1,neff,&a,f);

   printf(" <tr{U(p)}> = %1.6f\n\n",abar);

   print_plot(argv[1]);
   exit(0);
}
