/*******************************************************************************
*
* File read2.c
*
* Copyright (C) 2012-2014 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Reads and evaluates data from the data files created by the program ms1.
* The file to be read has to be specified on the command line.
*
* This program writes the history of the measured normalized reweighting
* factors to the file <run name>.run2.dat in the plots directory. The
* associated integrated  autocorrelation times are estimated and printed
* to stdout.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "extras.h"

static struct
{
   int nrw;
   int *nfct,*nsrc;
} file_head;

static struct
{
   int nc;
   double ***sqn,***lnr;
} data;

static int endian;
static int first,last,step,nms;
static double ***avrw,***lnrw,*avtot,*lntot;


static void read_file_head(FILE *fdat)
{
   int nrw,*nfct,*nsrc;
   int ir,ie,irw;
   stdint_t istd[1];

   ir=fread(istd,sizeof(stdint_t),1,fdat);
   error(ir!=1,1,"read_file_head [read2.c]",
         "Incorrect read count");

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   nrw=(int)(istd[0]);
   error(nrw<1,1,"read_file_head [read2.c]",
         "nrw is out of range");

   nfct=malloc(2*nrw*sizeof(*nfct));
   error(nfct==NULL,1,"read_file_head [read2.c]",
         "Unable to allocate data arrays");
   nsrc=nfct+nrw;
   ie=0;

   for (irw=0;irw<nrw;irw++)
   {
      ir+=fread(istd,sizeof(stdint_t),1,fdat);

      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);

      nfct[irw]=(int)(istd[0]);
      ie|=(nfct[irw]<1);
   }

   for (irw=0;irw<nrw;irw++)
   {
      ir+=fread(istd,sizeof(stdint_t),1,fdat);

      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);

      nsrc[irw]=(int)(istd[0]);
      ie|=(nsrc[irw]<1);
   }

   error(ir!=(1+2*nrw),1,"read_file_head [read2.c]",
         "Incorrect read count");
   error(ie!=0,1,"read_file_head [read2.c]",
         "Unexpected values of nfct or nsrc");

   file_head.nrw=nrw;
   file_head.nfct=nfct;
   file_head.nsrc=nsrc;
}


static void alloc_data(void)
{
   int nrw,*nfct,*nsrc;
   int i,irw,ifct,n1,n2,n3;
   double ***ppp,**pp,*p;

   nrw=file_head.nrw;
   nfct=file_head.nfct;
   nsrc=file_head.nsrc;
   n1=nrw;
   n2=0;
   n3=0;

   for (irw=0;irw<nrw;irw++)
   {
      n2+=nfct[irw];
      n3+=(nfct[irw]*nsrc[irw]);
   }

   ppp=malloc(2*n1*sizeof(*ppp));
   pp=malloc(2*n2*sizeof(*pp));
   p=malloc(2*n3*sizeof(*p));
   error((ppp==NULL)||(pp==NULL)||(p==NULL),1,"alloc_data [read2.c]",
         "Unable to allocate data arrays");

   data.sqn=ppp;
   data.lnr=ppp+nrw;

   for (i=0;i<2;i++)
   {
      for (irw=0;irw<nrw;irw++)
      {
         (*ppp)=pp;
         ppp+=1;

         for (ifct=0;ifct<nfct[irw];ifct++)
         {
            (*pp)=p;
            pp+=1;
            p+=nsrc[irw];
         }
      }
   }
}


static int read_data(FILE *fdat)
{
   int ir,n;
   int nrw,*nfct,*nsrc,irw,ifct,isrc;
   stdint_t istd[1];
   double dstd[1];

   ir=fread(istd,sizeof(stdint_t),1,fdat);

   if (ir!=1)
      return 0;

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   data.nc=(int)(istd[0]);

   nrw=file_head.nrw;
   nfct=file_head.nfct;
   nsrc=file_head.nsrc;
   n=0;

   for (irw=0;irw<nrw;irw++)
   {
      for (ifct=0;ifct<nfct[irw];ifct++)
      {
         for (isrc=0;isrc<nsrc[irw];isrc++)
         {
            ir+=fread(dstd,sizeof(double),1,fdat);

            if (endian==BIG_ENDIAN)
               bswap_double(1,dstd);

            data.sqn[irw][ifct][isrc]=dstd[0];
         }

         for (isrc=0;isrc<nsrc[irw];isrc++)
         {
            ir+=fread(dstd,sizeof(double),1,fdat);

            if (endian==BIG_ENDIAN)
               bswap_double(1,dstd);

            data.lnr[irw][ifct][isrc]=dstd[0];
         }

         n+=nsrc[irw];
      }
   }

   error(ir!=(1+2*n),1,"read_data [read2.c]",
         "Read error or incomplete data record");

   return 1;
}


static void cnfg_range(FILE *fdat,int *fst,int *lst,int *stp)
{
   int nc,ie;

   (*fst)=0;
   (*lst)=0;
   (*stp)=1;
   nc=0;
   ie=0;

   while (read_data(fdat))
   {
      nc+=1;

      if (nc==1)
         (*fst)=data.nc;
      else if (nc==2)
         (*stp)=data.nc-(*fst);
      else
         ie|=((data.nc-(*lst))!=(*stp));

      (*lst)=data.nc;
   }

   error(nc==0,1,"cnfg_range [read2.c]","No data records on data file");
   error(ie!=0,1,"cnfg_range [read2.c]","Non-contiguous configuration numbers");
}


static void select_cnfg_range(FILE *fdat)
{
   int fst,lst,stp;

   cnfg_range(fdat,&fst,&lst,&stp);

   printf("Available configuration range: %d - %d by %d\n",
          fst,lst,stp);
   printf("Select first,last,step: ");
   scanf("%d",&first);
   scanf(",");
   scanf("%d",&last);
   scanf(",");
   scanf("%d",&step);
   printf("\n");

   error((step%stp)!=0,1,"select_cnfg_range [read2.c]",
         "Step must be a multiple of the configuration separation");

   if (first<fst)
   {
      first=first+((fst-first)/step)*step;
      if (first<fst)
         first+=step;
   }

   if (last>lst)
   {
      last=last-((last-lst)/step)*step;
      if (last>lst)
         last-=step;
   }

   error((last<first)||(((last-first)%step)!=0)||(((first-fst)%stp)!=0),1,
         "select_cnfg_range [read2.c]","Improper configuration range");

   printf("Selected configuration range: %d - %d by %d\n\n",
          first,last,step);

   nms=(last-first)/step+1;
}


static void alloc_avrw(void)
{
   int nrw,*nfct,n1,n2,n3;
   int i,irw,ifct;
   double ***ppp,**pp,*p;

   nrw=file_head.nrw;
   nfct=file_head.nfct;
   n1=nrw;
   n2=0;

   for (irw=0;irw<nrw;irw++)
      n2+=nfct[irw];

   n3=n2*nms+nms;

   ppp=malloc(2*n1*sizeof(*ppp));
   pp=malloc(2*n2*sizeof(*pp));
   p=malloc(2*n3*sizeof(*p));
   error((ppp==NULL)||(pp==NULL)||(p==NULL),1,"alloc_avrw [read2.c]",
         "Unable to allocate data arrays");
   avrw=ppp;
   lnrw=ppp+n1;

   for (i=0;i<2;i++)
   {
      for (irw=0;irw<nrw;irw++)
      {
         (*ppp)=pp;
         ppp+=1;

         for (ifct=0;ifct<nfct[irw];ifct++)
         {
            (*pp)=p;
            pp+=1;
            p+=nms;
         }
      }
   }

   avtot=p;
   lntot=p+nms;
}


static void data2avrw(int ims)
{
   int nrw,*nfct,*nsrc;
   int irw,ifct,isrc;
   double lnm,rw,*lnr;

   nrw=file_head.nrw;
   nfct=file_head.nfct;
   nsrc=file_head.nsrc;

   avtot[ims]=1.0;
   lntot[ims]=0.0;

   for (irw=0;irw<nrw;irw++)
   {
      for (ifct=0;ifct<nfct[irw];ifct++)
      {
         lnr=data.lnr[irw][ifct];
         lnm=lnr[0];

         for (isrc=1;isrc<nsrc[irw];isrc++)
         {
            if (lnr[isrc]<lnm)
               lnm=lnr[isrc];
         }

         lnrw[irw][ifct][ims]=lnm;
         lntot[ims]+=lnm;

         rw=0.0;

         for (isrc=0;isrc<nsrc[ifct];isrc++)
            rw+=exp(lnm-lnr[isrc]);

         rw/=(double)(nsrc[ifct]);

         avrw[irw][ifct][ims]=rw;
         avtot[ims]*=rw;
      }
   }
}


static void normalize_avrw(void)
{
   int nrw,*nfct;
   int irw,ifct,ims;
   double lnma,rwa,*lnm,*rw;

   nrw=file_head.nrw;
   nfct=file_head.nfct;

   for (irw=0;irw<nrw;irw++)
   {
      for (ifct=1;ifct<nfct[irw];ifct++)
      {
         for (ims=0;ims<nms;ims++)
         {
            lnrw[irw][0][ims]+=lnrw[irw][ifct][ims];
            avrw[irw][0][ims]*=avrw[irw][ifct][ims];
         }
      }

      lnm=lnrw[irw][0];
      lnma=lnm[0];

      for (ims=1;ims<nms;ims++)
      {
         if (lnm[ims]<lnma)
            lnma=lnm[ims];
      }

      rw=avrw[irw][0];

      for (ims=0;ims<nms;ims++)
         rw[ims]*=exp(lnma-lnm[ims]);

      rwa=0.0;

      for (ims=0;ims<nms;ims++)
         rwa+=rw[ims];

      rwa/=(double)(nms);

      for (ims=0;ims<nms;ims++)
         rw[ims]/=rwa;
   }

   lnma=lntot[0];

   for (ims=1;ims<nms;ims++)
   {
      if (lntot[ims]<lnma)
         lnma=lntot[ims];
   }

   for (ims=0;ims<nms;ims++)
      avtot[ims]*=exp(lnma-lntot[ims]);

   rwa=0.0;

   for (ims=0;ims<nms;ims++)
      rwa+=avtot[ims];

   rwa/=(double)(nms);

   for (ims=0;ims<nms;ims++)
      avtot[ims]/=rwa;
}


static void read_file(char *fin)
{
   int nc,ims;
   long ipos;
   FILE *fdat;

   fdat=fopen(fin,"rb");
   error(fdat==NULL,1,"read_file [read2.c]","Unable to open data file");
   printf("Read data from file %s\n\n",fin);

   endian=endianness();
   read_file_head(fdat);
   alloc_data();

   ipos=ftell(fdat);
   select_cnfg_range(fdat);
   fseek(fdat,ipos,SEEK_SET);
   alloc_avrw();
   ims=0;

   while ((ims<nms)&&(read_data(fdat)))
   {
      nc=data.nc;

      if ((nc>=first)&&(nc<=last)&&(((nc-first)%step)==0))
      {
         data2avrw(ims);
         ims+=1;
      }
   }

   fclose(fdat);
   error((ims!=nms)||(data.nc!=last),1,"read_file [read2.c]",
         "Incorrect read count");

   normalize_avrw();
}


static double f(int nx,double x[])
{
   return x[0];
}


static void print_plot(char *fin)
{
   int n,nrw,irw,ims;
   char base[NAME_SIZE],plt_file[NAME_SIZE],*p;
   FILE *fout;

   p=strstr(fin,".ms1.dat");
   error(p==NULL,1,"print_plot [read2.c]","Unexpected data file name");
   n=p-fin;

   p=strrchr(fin,'/');
   if (p==NULL)
      p=fin;
   else
      p+=1;
   n-=(p-fin);

   error(n>=NAME_SIZE,1,"print_plot [read2.c]","File name is too long");
   strncpy(base,p,n);
   base[n]='\0';

   error(name_size("plots/%s.run2.dat",base)>=NAME_SIZE,1,
         "print_plot [read2.c]","File name is too long");
   sprintf(plt_file,"plots/%s.run2.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_plot [read2.c]",
         "Unable to open output file");

   nrw=file_head.nrw;

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms1\n");
   fprintf(fout,"# -------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# Number of measurements = %d\n",nms);
   fprintf(fout,"#\n");
   fprintf(fout,"# nc:   Configuration number\n");
   fprintf(fout,"# W:    Normalized reweighting factors\n");
   fprintf(fout,"#\n");
   fprintf(fout,"#  nc");

   for (irw=0;irw<nrw;irw++)
      fprintf(fout,"       W[%d] ",irw);

   if (nrw==1)
      fprintf(fout,"\n");
   else
      fprintf(fout,"       W[all]\n");

   fprintf(fout,"#\n");

   for (ims=0;ims<nms;ims++)
   {
      fprintf(fout," %5d  ",first+ims*step);

      for (irw=0;irw<nrw;irw++)
         fprintf(fout,"  %.4e",avrw[irw][0][ims]);

      if (nrw==1)
         fprintf(fout,"\n");
      else
         fprintf(fout,"  %.4e\n",avtot[ims]);
   }

   fclose(fout);

   printf("Data printed to file %s\n\n",plt_file);
}


int main(int argc,char *argv[])
{
   int nrw,irw,*nfct,*nsrc;

   error(argc!=2,1,"main [read2.c]","Syntax: read2 <filename>");

   printf("\n");
   printf("History of reweighting factors\n");
   printf("------------------------------\n\n");

   read_file(argv[1]);
   nrw=file_head.nrw;
   nfct=file_head.nfct;
   nsrc=file_head.nsrc;

   printf("The total number of measurements is %d.\n",nms);
   printf("Integrated autocorrelation times and associated errors are ");
   printf("estimated\n");

   if (nms>100)
      printf("using the numerically determined autocorrelation function.\n");
   else
      printf("by binning and calculating jackknife errors.\n");

   printf("Autocorrelation times are given in numbers of measurements.\n\n");

   for (irw=0;irw<nrw;irw++)
   {
      printf("Reweighting factor no %d:\n",irw);
      if (nfct[irw]>1)
         printf("Factorized into %d factors.\n",nfct[irw]);
      if (nsrc[irw]>1)
         printf("Using %d random source fields.\n\n",nsrc[irw]);
      else
         printf("Using 1 random source field.\n\n");

      if (nms>=100)
         print_auto(nms,avrw[irw][0]);
      else
         print_jack(1,nms,avrw[irw],f);

      printf("\n");
   }

   if (nrw!=1)
   {
      printf("Product of all reweighting factors:\n\n");

      if (nms>=100)
         print_auto(nms,avtot);
      else
         print_jack(1,nms,&avtot,f);

      printf("\n");
   }

   print_plot(argv[1]);
   exit(0);
}
