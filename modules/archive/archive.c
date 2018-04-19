
/*******************************************************************************
*
* File archive.c
*
* Copyright (C) 2005, 2007, 2009-2014 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs to read and write gauge-field configurations.
*
* The externally accessible functions are
*
*   void write_cnfg(char *out)
*     Writes the lattice sizes, the process grid sizes, the coordinates
*     of the calling process, the state of the random number generators,
*     the local plaquette sum and the local double-precision gauge field
*     to the file "out".
*
*   void read_cnfg(char *in)
*     Reads the local double-precision gauge field from the file "in",
*     assuming it was written to the file by the program write_cnfg().
*     The program then resets the random number generator and checks
*     that the restored field is compatible with the chosen boundary
*     conditions.
*
*   void export_cnfg(char *out)
*     Writes the lattice sizes and the global double-precision gauge
*     field to the file "out" from process 0 in the universal format
*     specified below (see the notes).
*
*   void import_cnfg(char *in)
*     Reads the global double-precision gauge field from the file "in"
*     on process 0, assuming the field was written to the file in the
*     universal format. The field is periodically extended if needed
*     and the program then checks that the configuration is compatible
*     with the chosen boundary conditions (see the notes).
*
* Notes:
*
* The program export_cnfg() first writes the lattice sizes and the average
* of the plaquette Re(tr{U(p)}) to the output file. Then follow the 8 link
* variables in the directions +0,-0,...,+3,-3 at the first odd point, the
* second odd point, and so on. The order of the point (x0,x1,x2,x3) with
* Cartesian coordinates in the range 0<=x0<N0,...,0<=x3<N3 is determined by
* the index
*
*   ix=x3+N3*x2+N2*N3*x1+N1*N2*N3*x0,
*
* where N0,N1,N2,N3 are the global lattice sizes (N0=NPROC0*L0, etc.). The
* average plaquette is calculated by summing the plaquette values over all
* plaquettes in the lattice, including the space-like ones at time N0 if
* SF or open-SF boundary conditions are chosen, and dividing the sum by
* 6*N0*N1*N2*N3.
*
* Independently of the machine, the export function writes the data to the
* output file in little-endian byte order. Integers and double-precision
* numbers on the output file occupy 4 and 8 bytes, respectively, the latter
* being formatted according to the IEEE-754 standard. The import function
* assumes the data on the input file to be little endian and converts them
* to big-endian order if the machine is big endian. Exported configurations
* can thus be safely exchanged between different machines.
*
* If the current lattice sizes N0,..,N3 are larger than the lattice sizes
* n0,..,n3 read from the configuration file, and if N0,..,N3 are integer
* multiples of n0,..,n3, the program import_cnfg() periodically extends the
* imported field. An extension in the time direction is only possible with
* periodic boundary conditions. Note that the boundary values of the link
* variables (rather than the angles characterizing them) must match in the
* case of SF and open-SF boundary conditions (see doc/gauge_action.pdf).
*
* Compatibility of a configuration with the chosen boundary conditions is
* established by calling check_bc() [lattice/bcnds.c], with a tolerance on
* the boundary link variables of 64.0*DBL_EPSILON, and by checking that the
* average plaquette coincides with the value read from the configuration
* file. On exit both read_cnfg() and import_cnfg() set the boundary values
* of the field (if any) to the ones stored in the parameter data base so
* as to guarantee that they are bit-identical to the latter.
*
* All programs in this module may involve global communications and must be
* called simultaneously on all processes.
*
*******************************************************************************/

#define ARCHIVE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "archive.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int endian,ns,nd,*state=NULL;
static su3_dble *ubuf=NULL,*vbuf,*udb;


static void alloc_state(void)
{
   int n;

   ns=rlxs_size();
   nd=rlxd_size();

   if (ns<nd)
      n=nd;
   else
      n=ns;

   state=malloc(n*sizeof(int));
   error(state==NULL,1,"alloc_state [archive.c]",
         "Unable to allocate auxiliary array");
}


void write_cnfg(char *out)
{
   int ldat[16],iw;
   double plaq;
   FILE *fout;

   error_root(query_flags(UD_PHASE_SET)==1,1,"write_cnfg [archive.c]",
             "Attempt to write phase-modified gauge field");
   fout=fopen(out,"wb");
   error_loc(fout==NULL,1,"write_cnfg [archive.c]",
             "Unable to open output file");

   if (state==NULL)
      alloc_state();

   udb=udfld();
   plaq=plaq_sum_dble(0);

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
   rlxs_get(state);
   iw+=fwrite(state,sizeof(int),ns,fout);
   rlxd_get(state);
   iw+=fwrite(state,sizeof(int),nd,fout);
   iw+=fwrite(&plaq,sizeof(double),1,fout);
   iw+=fwrite(udb,sizeof(su3_dble),4*VOLUME,fout);

   error_loc(iw!=(17+ns+nd+4*VOLUME),1,"write_cnfg [archive.c]",
             "Incorrect write count");
   fclose(fout);
}


void read_cnfg(char *in)
{
   int ldat[16],ir,ie;
   double nplaq,plaq0,plaq1,eps;
   FILE *fin;

   if (state==NULL)
      alloc_state();

   udb=udfld();
   unset_ud_phase();

   fin=fopen(in,"rb");
   error_loc(fin==NULL,1,"read_cnfg [archive.c]",
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
   error(ie!=0,1,"read_cnfg [archive.c]","Unexpected lattice data");

   ir+=fread(state,sizeof(int),ns,fin);
   rlxs_reset(state);
   ir+=fread(state,sizeof(int),nd,fin);
   rlxd_reset(state);
   ir+=fread(&plaq0,sizeof(double),1,fin);
   ir+=fread(udb,sizeof(su3_dble),4*VOLUME,fin);

   error_loc(ir!=(17+ns+nd+4*VOLUME),1,"read_cnfg [archive.c]",
             "Incorrect read count");
   fclose(fin);

   set_flags(UPDATED_UD);
   ie=check_bc(64.0*DBL_EPSILON);
   error_root(ie!=1,1,"read_cnfg [archive.c]",
              "Incompatible boundary conditions");

   ie=0;
   nplaq=(double)(6*L0*L1)*(double)(L2*L3);
   eps=sqrt(nplaq)*DBL_EPSILON;
   plaq0/=nplaq;
   plaq1=plaq_sum_dble(0)/nplaq;
   ie|=(fabs(plaq1-plaq0)>eps);
   set_bc();
   plaq1=plaq_sum_dble(0)/nplaq;
   ie|=(fabs(plaq1-plaq0)>eps);
   error_loc(ie!=0,1,"read_cnfg [archive.c]",
             "Incorrect average plaquette");
}


static void check_machine(void)
{
   error_root(sizeof(stdint_t)!=4,1,"check_machine [archive.c]",
              "Size of a stdint_t integer is not 4");
   error_root(sizeof(double)!=8,1,"check_machine [archive.c]",
              "Size of a double is not 8");

   endian=endianness();
   error_root(endian==UNKNOWN_ENDIAN,1,"check_machine [archive.c]",
              "Unkown endianness");
}


static void alloc_ubuf(int my_rank)
{
   if (my_rank==0)
   {
      ubuf=amalloc(4*(L3+N3)*sizeof(su3_dble),ALIGN);
      vbuf=ubuf+4*L3;
   }
   else
   {
      ubuf=amalloc(4*L3*sizeof(su3_dble),ALIGN);
      vbuf=NULL;
   }

   error(ubuf==NULL,1,"alloc_ubuf [archive.c]",
         "Unable to allocate auxiliary array");
}


static void get_links(int iy)
{
   int y3,ifc;
   su3_dble *u,*v;

   v=ubuf;
   iy*=L3;

   if (ipt[iy]<(VOLUME/2))
      iy+=1;

   for (y3=0;y3<L3;y3+=2)
   {
      u=udb+8*(ipt[iy+y3]-(VOLUME/2));

      for (ifc=0;ifc<8;ifc++)
      {
         v[0]=u[0];
         v+=1;
         u+=1;
      }
   }
}


static void set_links(int iy)
{
   int y3,ifc;
   su3_dble *u,*v;

   v=ubuf;
   iy*=L3;

   if (ipt[iy]<(VOLUME/2))
      iy+=1;

   for (y3=0;y3<L3;y3+=2)
   {
      u=udb+8*(ipt[iy+y3]-(VOLUME/2));

      for (ifc=0;ifc<8;ifc++)
      {
         u[0]=v[0];
         v+=1;
         u+=1;
      }
   }
}


void export_cnfg(char *out)
{
   int my_rank,np[4],n,iw;
   int iwa,dmy,tag0,tag1;
   int x0,x1,x2,x3,y0,y1,y2,ix,iy;
   stdint_t lsize[4];
   double nplaq,plaq;
   MPI_Status stat;
   FILE *fout=NULL;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   error_root(query_flags(UD_PHASE_SET)==1,1,"export_cnfg [archive.c]",
             "Attempt to export phase-modified gauge field");

   if (ubuf==NULL)
   {
      check_machine();
      alloc_ubuf(my_rank);
   }

   dmy=1;
   tag0=mpi_tag();
   tag1=mpi_tag();
   nplaq=(double)(6*N0*N1)*(double)(N2*N3);
   plaq=plaq_sum_dble(1)/nplaq;

   if (my_rank==0)
   {
      fout=fopen(out,"wb");
      error_root(fout==NULL,1,"export_cnfg [archive.c]",
                 "Unable to open output file");

      lsize[0]=(stdint_t)(N0);
      lsize[1]=(stdint_t)(N1);
      lsize[2]=(stdint_t)(N2);
      lsize[3]=(stdint_t)(N3);

      if (endian==BIG_ENDIAN)
      {
         bswap_int(4,lsize);
         bswap_double(1,&plaq);
      }

      iw=fwrite(lsize,sizeof(stdint_t),4,fout);
      iw+=fwrite(&plaq,sizeof(double),1,fout);

      error_root(iw!=5,1,"export_cnfg [archive.c]","Incorrect write count");
   }

   iwa=0;
   udb=udfld();

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
            get_links(iy);

         if (n>0)
         {
            if (my_rank==0)
            {
               MPI_Send(&dmy,1,MPI_INT,n,tag0,MPI_COMM_WORLD);
               MPI_Recv(ubuf,4*L3*18,MPI_DOUBLE,n,tag1,MPI_COMM_WORLD,&stat);
            }
            else if (my_rank==n)
            {
               MPI_Recv(&dmy,1,MPI_INT,0,tag0,MPI_COMM_WORLD,&stat);
               MPI_Send(ubuf,4*L3*18,MPI_DOUBLE,0,tag1,MPI_COMM_WORLD);
            }
         }

         if (my_rank==0)
         {
            if (endian==BIG_ENDIAN)
               bswap_double(4*L3*18,ubuf);
            iw=fwrite(ubuf,sizeof(su3_dble),4*L3,fout);
            iwa|=(iw!=(4*L3));
         }
      }
   }

   if (my_rank==0)
   {
      error_root(iwa!=0,1,"export_cnfg [archive.c]",
                 "Incorrect write count");
      fclose(fout);
   }
}


void import_cnfg(char *in)
{
   int my_rank,np[4],ir,ie;
   int ira,dmy,tag0,tag1;
   int n0,n1,n2,n3,nc0,nc1,nc2,nc3;
   int x0,x1,x2,y0,y1,y2,y3,c0,c1,c2,ix,iy,ic;
   int n,k,l;
   stdint_t lsize[4];
   double nplaq,plaq0,plaq1,eps;
   MPI_Status stat;
   FILE *fin=NULL;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (ubuf==NULL)
   {
      check_machine();
      alloc_ubuf(my_rank);
   }

   dmy=1;
   tag0=mpi_tag();
   tag1=mpi_tag();
   udb=udfld();
   unset_ud_phase();

   if (my_rank==0)
   {
      fin=fopen(in,"rb");
      error_root(fin==NULL,1,"import_cnfg [archive.c]",
                 "Unable to open input file");

      ir=fread(lsize,sizeof(stdint_t),4,fin);
      ir+=fread(&plaq0,sizeof(double),1,fin);
      error_root(ir!=5,1,"import_cnfg [archive.c]","Incorrect read count");

      if (endian==BIG_ENDIAN)
      {
         bswap_int(4,lsize);
         bswap_double(1,&plaq0);
      }

      np[0]=(int)(lsize[0]);
      np[1]=(int)(lsize[1]);
      np[2]=(int)(lsize[2]);
      np[3]=(int)(lsize[3]);

      error_root((np[0]<1)||((N0%np[0])!=0)||
                 (np[1]<1)||((N1%np[1])!=0)||
                 (np[2]<1)||((N2%np[2])!=0)||
                 (np[3]<1)||((N3%np[3])!=0),1,"import_cnfg [archive.c]",
                 "Unexpected or incompatible lattice sizes");

      error_root((np[0]!=N0)&&(bc_type()!=3),1,"import_cnfg [archive.c]",
                 "Periodic extension in time is only possible when\n"
                 "periodic boundary conditions are chosen");
   }
   else
   {
      np[0]=0;
      np[1]=0;
      np[2]=0;
      np[3]=0;
      plaq0=0.0;
   }

   MPI_Bcast(np,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&plaq0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   n0=np[0];
   n1=np[1];
   n2=np[2];
   n3=np[3];

   nc0=N0/n0;
   nc1=N1/n1;
   nc2=N2/n2;
   nc3=N3/n3;
   ira=0;

   for (ix=0;ix<(n0*n1*n2);ix++)
   {
      x0=ix/(n1*n2);
      x1=(ix/n2)%n1;
      x2=ix%n2;

      if (my_rank==0)
      {
         n=4*n3;
         ir=fread(vbuf,sizeof(su3_dble),n,fin);
         ira|=(ir!=n);

         if (endian==BIG_ENDIAN)
            bswap_double(n*18,vbuf);

         for (k=1;k<nc3;k++)
         {
            for (l=0;l<n;l++)
               vbuf[k*n+l]=vbuf[l];
         }
      }

      for (ic=0;ic<(nc0*nc1*nc2);ic++)
      {
         c0=ic/(nc1*nc2);
         c1=(ic/nc2)%nc1;
         c2=ic%nc2;

         y0=x0+c0*n0;
         y1=x1+c1*n1;
         y2=x2+c2*n2;
         iy=(y2%L2)+L2*(y1%L1)+L1*L2*(y0%L0);

         np[0]=y0/L0;
         np[1]=y1/L1;
         np[2]=y2/L2;

         for (y3=0;y3<N3;y3+=L3)
         {
            np[3]=y3/L3;
            n=ipr_global(np);

            if (n>0)
            {
               if (my_rank==0)
               {
                  MPI_Send(vbuf+4*y3,4*L3*18,MPI_DOUBLE,n,tag1,MPI_COMM_WORLD);
                  MPI_Recv(&dmy,1,MPI_INT,n,tag0,MPI_COMM_WORLD,&stat);
               }
               else if (my_rank==n)
               {
                  MPI_Recv(ubuf,4*L3*18,MPI_DOUBLE,0,tag1,MPI_COMM_WORLD,&stat);
                  MPI_Send(&dmy,1,MPI_INT,0,tag0,MPI_COMM_WORLD);
               }
            }
            else if (my_rank==0)
               for (l=0;l<(4*L3);l++)
                  ubuf[l]=vbuf[4*y3+l];

            if (my_rank==n)
               set_links(iy);
         }
      }
   }

   if (my_rank==0)
   {
      error_root(ira!=0,1,"import_cnfg [archive.c]","Incorrect read count");
      fclose(fin);
   }

   set_flags(UPDATED_UD);
   ie=check_bc(64.0*DBL_EPSILON);
   error_root(ie!=1,1,"import_cnfg [archive.c]",
              "Incompatible boundary conditions");

   ie=0;
   nplaq=(double)(6*N0*N1)*(double)(N2*N3);
   eps=sqrt(nplaq)*DBL_EPSILON;
   plaq1=plaq_sum_dble(1)/nplaq;
   ie|=(fabs(plaq1-plaq0)>eps);
   set_bc();
   plaq1=plaq_sum_dble(1)/nplaq;
   ie|=(fabs(plaq1-plaq0)>eps);
   error_root(ie!=0,1,"import_cnfg [archive.c]",
              "Incorrect average plaquette");
}
