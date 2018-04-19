
/*******************************************************************************
*
* File sarchive.c
*
* Copyright (C) 2007, 2008, 2011, 2013, 2014 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs to read and write global double-precision spinor fields.
*
* The externally accessible functions are
*
*   void write_sfld(char *out,spinor_dble *sd)
*     Writes the lattice sizes, the process grid sizes, the coordinates
*     of the calling process, the square of the norm of the spinor field
*     sd and the local part of the latter to the file "out".
*
*   void read_sfld(char *in,spinor_dble *sd)
*     Reads the local part of the spinor field sd from the file "in",
*     assuming the field was written to the file by write_sfld().
*
*   void export_sfld(char *out,spinor_dble *sd)
*     Writes the lattice sizes and the spinor field sd to the file "out"
*     from process 0 in the universal format specified below (see the
*     notes).
*
*   void import_sfld(char *in,spinor_dble *sd)
*     Reads the spinor field sd from the file "in" on process 0, assuming
*     the field was written to the file in the universal format (see the
*     notes).
*
* Notes:
*
* The spinor fields are assumed to be global quark fields as described in
* main/README.global. Only their physical components (i.e. the spinors on
* the local lattices) are written and read.
*
* The program export_sfld() first writes the global lattice sizes and the
* square-norm of the spinor field. Then follow the spinors at the first
* lattice point, the second point, and so on, in the order given by the
* index
*
*   ix=x3+N3*x2+N2*N3*x1+N1*N2*N3*x0,
*
* where N0,N1,N2,N3 are the (global) lattice sizes and (x0,x1,x2,x3) the
* Cartesian coordinates of the points (0<=x0<N0,...,0<=x3<N3).
*
* Independently of the machine, the export function writes the data to the
* output file in little-endian byte order. Integers and double-precision
* numbers on the output file occupy 4 and 8 bytes, respectively, the latter
* being formatted according to the IEEE-754 standard. The import function
* assumes the data on the input file to be little endian and converts them
* to big-endian order if the machine is big endian. Exported fields can
* thus be safely exchanged between different machines.
*
* All programs in this module involve communications and must be called
* simultaneously on all MPI processes.
*
*******************************************************************************/

#define SARCHIVE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "lattice.h"
#include "linalg.h"
#include "archive.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int endian;
static spinor_dble *sbuf=NULL;


void write_sfld(char *out,spinor_dble *sd)
{
   int ldat[16],iw;
   double norm;
   FILE *fout=NULL;

   error(sd==NULL,1,"write_sfld [sarchive.c]",
         "Attempt to access unallocated memory space");
   error(iup[0][0]==0,1,"write_sfld [sarchive.c]",
         "Geometry arrays are not set");

   fout=fopen(out,"wb");
   error_loc(fout==NULL,1,"write_sfld [sarchive.c]",
             "Unable to open output file");

   ldat[0]=NPROC0;
   ldat[1]=NPROC1;
   ldat[2]=NPROC2;
   ldat[3]=NPROC3;

   ldat[4]=L0;
   ldat[5]=L1;
   ldat[6]=L2;
   ldat[7]=L3;

   ldat[8]=NPROC0_BLK;
   ldat[9]=NPROC1_BLK;
   ldat[10]=NPROC2_BLK;
   ldat[11]=NPROC3_BLK;

   ldat[12]=cpr[0];
   ldat[13]=cpr[1];
   ldat[14]=cpr[2];
   ldat[15]=cpr[3];

   iw=fwrite(ldat,sizeof(int),16,fout);
   norm=norm_square_dble(VOLUME,0,sd);
   iw+=fwrite(&norm,sizeof(double),1,fout);
   iw+=fwrite(sd,sizeof(spinor_dble),VOLUME,fout);

   error_loc(iw!=(17+VOLUME),1,"write_sfld [sarchive.c]",
             "Incorrect write count");
   fclose(fout);
}


void read_sfld(char *in,spinor_dble *sd)
{
   int ldat[16],ir,ie;
   double norm0,norm1,eps;
   FILE *fin=NULL;

   error(sd==NULL,1,"read_sfld [sarchive.c]",
         "Attempt to access unallocated memory space");
   error(iup[0][0]==0,1,"read_sfld [sarchive.c]",
         "Geometry arrays are not set");

   fin=fopen(in,"rb");
   error_loc(fin==NULL,1,"read_sfld [sarchive.c]",
             "Unable to open input file");

   ir=fread(ldat,sizeof(int),16,fin);

   ie=0;
   ie|=((ldat[0]!=NPROC0)||(ldat[1]!=NPROC1)||
        (ldat[2]!=NPROC2)||(ldat[3]!=NPROC3));
   ie|=((ldat[4]!=L0)||(ldat[5]!=L1)||
        (ldat[6]!=L2)||(ldat[7]!=L3));
   ie|=((ldat[8]!=NPROC0_BLK)||(ldat[9]!=NPROC1_BLK)||
        (ldat[10]!=NPROC2_BLK)||(ldat[11]!=NPROC3_BLK));
   ie|=((ldat[12]!=cpr[0])||(ldat[13]!=cpr[1])||
        (ldat[14]!=cpr[2])||(ldat[15]!=cpr[3]));
   error(ie!=0,1,"read_sfld [sarchive.c]","Unexpected lattice data");

   ir+=fread(&norm0,sizeof(double),1,fin);
   ir+=fread(sd,sizeof(spinor_dble),VOLUME,fin);

   error_loc(ir!=(17+VOLUME),1,"read_sfld [sarchive.c]",
             "Incorrect read count");
   fclose(fin);

   norm1=norm_square_dble(VOLUME,0,sd);
   eps=sqrt(64.0*(double)(VOLUME))*DBL_EPSILON;
   error_loc(fabs(norm1-norm0)>(eps*norm0),1,"read_sfld [sarchive.c]",
             "Incorrect square norm");
}


static void check_machine(void)
{
   error_root(sizeof(stdint_t)!=4,1,"check_machine [sarchive.c]",
              "Size of a stdint_t integer is not 4");
   error_root(sizeof(double)!=8,1,"check_machine [sarchive.c]",
              "Size of a double is not 8");
   error_root(sizeof(spinor_dble)!=192,1,"check_machine [sarchive.c]",
              "The spinor_dble structures are not properly packed");

   endian=endianness();
   error_root(endian==UNKNOWN_ENDIAN,1,"check_machine [sarchive.c]",
              "Unkown endianness");
}


static void alloc_sbuf(void)
{
   error(iup[0][0]==0,1,"alloc_sbuf [sarchive.c]",
         "Geometry arrays are not set");
   sbuf=amalloc(L3*sizeof(spinor_dble),ALIGN);
   error(sbuf==NULL,1,"alloc_sbuf [sarchive.c]",
         "Unable to allocate auxiliary array");
}


static void get_spinors(int iy,spinor_dble *sd)
{
   int y3,iz;
   spinor_dble *sb;

   sb=sbuf;
   iy*=L3;

   for (y3=0;y3<L3;y3++)
   {
      iz=ipt[iy+y3];
      (*sb)=sd[iz];
      sb+=1;
   }
}


static void set_spinors(int iy,spinor_dble *sd)
{
   int y3,iz;
   spinor_dble *sb;

   sb=sbuf;
   iy*=L3;

   for (y3=0;y3<L3;y3++)
   {
      iz=ipt[iy+y3];
      sd[iz]=(*sb);
      sb+=1;
   }
}


void export_sfld(char *out,spinor_dble *sd)
{
   int my_rank,np[4],n,iw;
   int iwa,dmy,tag0,tag1;
   int x0,x1,x2,x3,y0,y1,y2,ix,iy;
   stdint_t lsize[4];
   double norm;
   MPI_Status stat;
   FILE *fout=NULL;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (sbuf==NULL)
   {
      check_machine();
      alloc_sbuf();
   }

   error(sd==NULL,1,"export_sfld [sarchive.c]",
         "Attempt to access unallocated memory space");

   dmy=1;
   tag0=mpi_tag();
   tag1=mpi_tag();
   norm=norm_square_dble(VOLUME,1,sd);

   if (my_rank==0)
   {
      fout=fopen(out,"wb");
      error_root(fout==NULL,1,"export_sfld [sarchive.c]",
                 "Unable to open output file");

      lsize[0]=N0;
      lsize[1]=N1;
      lsize[2]=N2;
      lsize[3]=N3;

      if (endian==BIG_ENDIAN)
      {
         bswap_int(4,lsize);
         bswap_double(1,&norm);
      }

      iw=fwrite(lsize,sizeof(stdint_t),4,fout);
      iw+=fwrite(&norm,sizeof(double),1,fout);
      error_root(iw!=5,1,"export_sfld [sarchive.c]","Incorrect write count");
   }

   iwa=0;

   for (ix=0;ix<(N0*N1*N2);ix++)
   {
      x0=ix/(N1*N2);
      x1=(ix/N2)%N1;
      x2=ix%N2;

      y0=x0%L0;
      y1=x1%L1;
      y2=x2%L2;
      iy=y2+L2*y1+L1*L2*y0;

      np[0]=x0/L0;
      np[1]=x1/L1;
      np[2]=x2/L2;

      for (x3=0;x3<N3;x3+=L3)
      {
         np[3]=x3/L3;
         n=ipr_global(np);
         if (my_rank==n)
            get_spinors(iy,sd);

         if (n>0)
         {
            if (my_rank==0)
            {
               MPI_Send(&dmy,1,MPI_INT,n,tag0,MPI_COMM_WORLD);
               MPI_Recv(sbuf,L3*24,MPI_DOUBLE,n,tag1,MPI_COMM_WORLD,&stat);
            }
            else if (my_rank==n)
            {
               MPI_Recv(&dmy,1,MPI_INT,0,tag0,MPI_COMM_WORLD,&stat);
               MPI_Send(sbuf,L3*24,MPI_DOUBLE,0,tag1,MPI_COMM_WORLD);
            }
         }

         if (my_rank==0)
         {
            if (endian==BIG_ENDIAN)
               bswap_double(L3*24,(double*)(sbuf));
            iw=fwrite(sbuf,sizeof(spinor_dble),L3,fout);
            iwa|=(iw!=L3);
         }
      }
   }

   if (my_rank==0)
   {
      error_root(iwa!=0,1,"export_sfld [sarchive.c]","Incorrect write count");
      fclose(fout);
   }
}


void import_sfld(char *in,spinor_dble *sd)
{
   int my_rank,np[4],n,ir;
   int ira,dmy,tag0,tag1;
   int x0,x1,x2,x3,y0,y1,y2,ix,iy;
   stdint_t lsize[4];
   double norm0,norm1,eps;
   MPI_Status stat;
   FILE *fin=NULL;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (sbuf==NULL)
   {
      check_machine();
      alloc_sbuf();
   }

   error(sd==NULL,1,"import_sfld [sarchive.c]",
         "Attempt to access unallocated memory space");

   dmy=1;
   tag0=mpi_tag();
   tag1=mpi_tag();

   if (my_rank==0)
   {
      fin=fopen(in,"rb");
      error_root(fin==NULL,1,"import_sfld [sarchive.c]",
                 "Unable to open input file");

      ir=fread(lsize,sizeof(stdint_t),4,fin);
      ir+=fread(&norm0,sizeof(double),1,fin);
      error_root(ir!=5,1,"import_sfld [sarchive.c]","Incorrect read count");

      if (endian==BIG_ENDIAN)
      {
         bswap_int(4,lsize);
         bswap_double(1,&norm0);
      }

      error_root((lsize[0]!=N0)||(lsize[1]!=N1)||(lsize[2]!=N2)||
                 (lsize[3]!=N3),1,"import_sfld [sarchive.c]",
                 "Lattice sizes do not match");
   }
   else
      norm0=0.0;

   ira=0;

   for (ix=0;ix<(N0*N1*N2);ix++)
   {
      x0=ix/(N1*N2);
      x1=(ix/N2)%N1;
      x2=ix%N2;

      y0=x0%L0;
      y1=x1%L1;
      y2=x2%L2;
      iy=y2+L2*y1+L1*L2*y0;

      np[0]=x0/L0;
      np[1]=x1/L1;
      np[2]=x2/L2;

      for (x3=0;x3<N3;x3+=L3)
      {
         np[3]=x3/L3;
         n=ipr_global(np);

         if (my_rank==0)
         {
            ir=fread(sbuf,sizeof(spinor_dble),L3,fin);
            ira|=(ir!=L3);

            if (endian==BIG_ENDIAN)
               bswap_double(L3*24,(double*)(sbuf));
         }

         if (n>0)
         {
            if (my_rank==0)
            {
               MPI_Send(sbuf,L3*24,MPI_DOUBLE,n,tag1,MPI_COMM_WORLD);
               MPI_Recv(&dmy,1,MPI_INT,n,tag0,MPI_COMM_WORLD,&stat);
            }
            else if (my_rank==n)
            {
               MPI_Recv(sbuf,L3*24,MPI_DOUBLE,0,tag1,MPI_COMM_WORLD,&stat);
               MPI_Send(&dmy,1,MPI_INT,0,tag0,MPI_COMM_WORLD);
            }
         }

         if (my_rank==n)
            set_spinors(iy,sd);
      }
   }

   if (my_rank==0)
   {
      error_root(ira!=0,1,"import_sfld [sarchive.c]","Incorrect read count");
      fclose(fin);
   }

   norm1=norm_square_dble(VOLUME,1,sd);
   eps=sqrt(64.0*(double)(N0*N1)*(double)(N2*N3))*DBL_EPSILON;
   error_root(fabs(norm1-norm0)>(eps*norm0),1,"import_sfld [sarchive.c]",
              "Incorrect square norm");
}
