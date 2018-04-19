
/*******************************************************************************
*
* File ym1.c
*
* Copyright (C) 2010-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* HMC simulation program for the SU(3) gauge theory.
*
* Syntax: ym1 -i <filename> [-noloc] [-noexp] [-rmold] [-noms]
*                           [-c <filename> [-a [-norng]]]
*
* For usage instructions see the file README.ym1.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "flags.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "archive.h"
#include "forces.h"
#include "update.h"
#include "wflow.h"
#include "tcharge.h"
#include "version.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

typedef struct
{
   int nt,iac;
   double dH,avpl;
} dat_t;

static struct
{
   int dn,nn,tmax;
   double eps;
} file_head;

static struct
{
   int nt;
   double **Wsl,**Ysl,**Qsl;
} data;

static int my_rank,noloc,noexp,rmold,noms,norng;
static int scnfg,append,endian;
static int level,seed;
static int nth,ntr,dtr_log,dtr_ms,dtr_cnfg;
static int ipgrd[2],flint;
static double *Wact,*Yact,*Qtop;

static char line[NAME_SIZE];
static char log_dir[NAME_SIZE],dat_dir[NAME_SIZE];
static char loc_dir[NAME_SIZE],cnfg_dir[NAME_SIZE];
static char log_file[NAME_SIZE],log_save[NAME_SIZE];
static char par_file[NAME_SIZE],par_save[NAME_SIZE];
static char dat_file[NAME_SIZE],dat_save[NAME_SIZE];
static char msdat_file[NAME_SIZE],msdat_save[NAME_SIZE];
static char rng_file[NAME_SIZE],rng_save[NAME_SIZE];
static char cnfg_file[NAME_SIZE],end_file[NAME_SIZE];
static char nbase[NAME_SIZE],cnfg[NAME_SIZE];
static FILE *fin=NULL,*flog=NULL,*fdat=NULL,*fend=NULL;

static lat_parms_t lat;
static bc_parms_t bcp;
static hmc_parms_t hmc;


static int write_dat(int n,dat_t *ndat)
{
   int i,iw,ic;
   stdint_t istd[2];
   double dstd[2];

   ic=0;

   for (i=0;i<n;i++)
   {
      istd[0]=(stdint_t)((*ndat).nt);
      istd[1]=(stdint_t)((*ndat).iac);

      dstd[0]=(*ndat).dH;
      dstd[1]=(*ndat).avpl;

      if (endian==BIG_ENDIAN)
      {
         bswap_int(2,istd);
         bswap_double(2,dstd);
      }

      iw=fwrite(istd,sizeof(stdint_t),2,fdat);
      iw+=fwrite(dstd,sizeof(double),2,fdat);

      if (iw!=4)
         return ic;

      ic+=1;
      ndat+=1;
   }

   return ic;
}


static int read_dat(int n,dat_t *ndat)
{
   int i,ir,ic;
   stdint_t istd[2];
   double dstd[2];

   ic=0;

   for (i=0;i<n;i++)
   {
      ir=fread(istd,sizeof(stdint_t),2,fdat);
      ir+=fread(dstd,sizeof(double),2,fdat);

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


static void alloc_data(void)
{
   int nn,tmax;
   int in;
   double **pp,*p;

   nn=file_head.nn;
   tmax=file_head.tmax;

   pp=amalloc(3*(nn+1)*sizeof(*pp),3);
   p=amalloc(3*(nn+1)*(tmax+1)*sizeof(*p),4);

   error((pp==NULL)||(p==NULL),1,"alloc_data [ym1.c]",
         "Unable to allocate data arrays");

   data.Wsl=pp;
   data.Ysl=pp+nn+1;
   data.Qsl=pp+2*(nn+1);

   for (in=0;in<(3*(nn+1));in++)
   {
      *pp=p;
      pp+=1;
      p+=tmax;
   }

   Wact=p;
   p+=nn+1;
   Yact=p;
   p+=nn+1;
   Qtop=p;
}


static void write_file_head(void)
{
   int iw;
   stdint_t istd[3];
   double dstd[1];

   istd[0]=(stdint_t)(file_head.dn);
   istd[1]=(stdint_t)(file_head.nn);
   istd[2]=(stdint_t)(file_head.tmax);
   dstd[0]=file_head.eps;

   if (endian==BIG_ENDIAN)
   {
      bswap_int(3,istd);
      bswap_double(1,dstd);
   }

   iw=fwrite(istd,sizeof(stdint_t),3,fdat);
   iw+=fwrite(dstd,sizeof(double),1,fdat);

   error_root(iw!=4,1,"write_file_head [ym1.c]",
              "Incorrect write count");
}


static void check_file_head(void)
{
   int ir;
   stdint_t istd[3];
   double dstd[1];

   ir=fread(istd,sizeof(stdint_t),3,fdat);
   ir+=fread(dstd,sizeof(double),1,fdat);

   error_root(ir!=4,1,"check_file_head [ym1.c]",
              "Incorrect read count");

   if (endian==BIG_ENDIAN)
   {
      bswap_int(3,istd);
      bswap_double(1,dstd);
   }

   error_root(((int)(istd[0])!=file_head.dn)||
              ((int)(istd[1])!=file_head.nn)||
              ((int)(istd[2])!=file_head.tmax)||
              (dstd[0]!=file_head.eps),1,"check_file_head [ym1.c]",
              "Unexpected value of dn,nn,tmax or eps");
}


static void write_data(void)
{
   int iw,nn,tmax;
   int in,t;
   stdint_t istd[1];
   double dstd[1];

   istd[0]=(stdint_t)(data.nt);

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   iw=fwrite(istd,sizeof(stdint_t),1,fdat);

   nn=file_head.nn;
   tmax=file_head.tmax;

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         dstd[0]=data.Wsl[in][t];

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         iw+=fwrite(dstd,sizeof(double),1,fdat);
      }
   }

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         dstd[0]=data.Ysl[in][t];

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         iw+=fwrite(dstd,sizeof(double),1,fdat);
      }
   }

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         dstd[0]=data.Qsl[in][t];

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         iw+=fwrite(dstd,sizeof(double),1,fdat);
      }
   }

   error_root(iw!=(1+3*(nn+1)*tmax),1,"write_data [ym1.c]",
              "Incorrect write count");
}


static int read_data(void)
{
   int ir,nn,tmax;
   int in,t;
   stdint_t istd[1];
   double dstd[1];

   ir=fread(istd,sizeof(stdint_t),1,fdat);

   if (ir!=1)
      return 0;

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   data.nt=(int)(istd[0]);

   nn=file_head.nn;
   tmax=file_head.tmax;

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         ir+=fread(dstd,sizeof(double),1,fdat);

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         data.Wsl[in][t]=dstd[0];
      }
   }

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         ir+=fread(dstd,sizeof(double),1,fdat);

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         data.Ysl[in][t]=dstd[0];
      }
   }

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         ir+=fread(dstd,sizeof(double),1,fdat);

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         data.Qsl[in][t]=dstd[0];
      }
   }

   error_root(ir!=(1+3*(nn+1)*tmax),1,"read_data [ym1.c]",
              "Read error or incomplete data record");

   return 1;
}


static void read_dirs(void)
{
   if (my_rank==0)
   {
      find_section("Run name");
      read_line("name","%s",nbase);

      find_section("Directories");
      read_line("log_dir","%s",log_dir);
      read_line("dat_dir","%s",dat_dir);
      if (noloc==0)
         read_line("loc_dir","%s",loc_dir);
      else
         loc_dir[0]='\0';
      if ((noexp==0)||((scnfg)&&(strchr(cnfg,'*')==NULL)))
         read_line("cnfg_dir","%s",cnfg_dir);
      else
         cnfg_dir[0]='\0';
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(log_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(dat_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(loc_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
}


static void setup_files(void)
{
   check_dir_root(log_dir);
   check_dir_root(dat_dir);

   if (noloc==0)
      check_dir(loc_dir);

   if (noexp==0)
      check_dir_root(cnfg_dir);

   error_root(name_size("%s/%s.log~",log_dir,nbase)>=NAME_SIZE,1,
              "setup_files [ym1.c]","log_dir name is too long");
   error_root(name_size("%s/%s.ms.dat~",dat_dir,nbase)>=NAME_SIZE,1,
              "setup_files [ym1.c]","dat_dir name is too long");

   sprintf(log_file,"%s/%s.log",log_dir,nbase);
   sprintf(par_file,"%s/%s.par",dat_dir,nbase);
   sprintf(dat_file,"%s/%s.dat",dat_dir,nbase);
   sprintf(msdat_file,"%s/%s.ms.dat",dat_dir,nbase);
   sprintf(rng_file,"%s/%s.rng",dat_dir,nbase);
   sprintf(end_file,"%s/%s.end",log_dir,nbase);
   sprintf(log_save,"%s~",log_file);
   sprintf(par_save,"%s~",par_file);
   sprintf(dat_save,"%s~",dat_file);
   sprintf(msdat_save,"%s~",msdat_file);
   sprintf(rng_save,"%s~",rng_file);
}


static void read_lat_parms(void)
{
   double beta,c0;

   if (my_rank==0)
   {
      find_section("Lattice parameters");
      read_line("beta","%lf",&beta);
      read_line("c0","%lf",&c0);
   }

   MPI_Bcast(&beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&c0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   lat=set_lat_parms(beta,c0,0,NULL,1.0);

   if (append)
      check_lat_parms(fdat);
   else
      write_lat_parms(fdat);
}


static void read_bc_parms(void)
{
   int bc;
   double cG,cG_prime;
   double phi[2],phi_prime[2],theta[3];

   if (my_rank==0)
   {
      find_section("Boundary conditions");
      read_line("type","%d",&bc);

      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;
      cG=1.0;
      cG_prime=1.0;

      if (bc==1)
         read_dprms("phi",2,phi);

      if ((bc==1)||(bc==2))
         read_dprms("phi'",2,phi_prime);

      if (bc!=3)
         read_line("cG","%lf",&cG);

      if (bc==2)
         read_line("cG'","%lf",&cG_prime);
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(phi,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(phi_prime,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cG,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cG_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;

   bcp=set_bc_parms(bc,cG,cG_prime,1.0,1.0,phi,phi_prime,theta);

   if (append)
      check_bc_parms(fdat);
   else
      write_bc_parms(fdat);
}


static void read_hmc_parms(void)
{
   int iact[1];
   double tau;

   if (my_rank==0)
   {
      find_section("Trajectory length");
      read_line("tau","%lf",&tau);
   }

   MPI_Bcast(&tau,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   iact[0]=0;
   hmc=set_hmc_parms(1,iact,0,0,NULL,1,tau);

   if (append)
      check_hmc_parms(fdat);
   else
      write_hmc_parms(fdat);
}


static void read_integrator(void)
{
   int nstep,imd,ifr[1];
   double lambda;

   if (my_rank==0)
   {
      find_section("MD integrator");
      read_line("integrator","%s",line);
      lambda=0.0;

      if (strcmp(line,"LPFR")==0)
         imd=(int)(LPFR);
      else if (strcmp(line,"OMF2")==0)
      {
         imd=(int)(OMF2);
         read_line("lambda","%lf",&lambda);
      }
      else if (strcmp(line,"OMF4")==0)
         imd=(int)(OMF4);
      else
         error_root(1,1,"read_integrator [ym1.c]","Unknown integrator");

      read_line("nstep","%d",&nstep);
   }

   MPI_Bcast(&imd,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&lambda,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&nstep,1,MPI_INT,0,MPI_COMM_WORLD);

   ifr[0]=0;

   if (imd==(int)(LPFR))
      set_mdint_parms(0,LPFR,lambda,nstep,1,ifr);
   else if (imd==(int)(OMF2))
      set_mdint_parms(0,OMF2,lambda,nstep,1,ifr);
   else if (imd==(int)(OMF4))
      set_mdint_parms(0,OMF4,lambda,nstep,1,ifr);

   set_action_parms(0,ACG,0,0,NULL,NULL,NULL);
   set_force_parms(0,FRG,0,0,NULL,NULL,NULL,NULL);

   if (append)
   {
      check_mdint_parms(fdat);
      check_action_parms(fdat);
      check_force_parms(fdat);
   }
   else
   {
      write_mdint_parms(fdat);
      write_action_parms(fdat);
      write_force_parms(fdat);
   }
}


static void read_schedule(void)
{
   int ie,ir,iw;
   stdint_t istd[3];

   if (my_rank==0)
   {
      find_section("MD trajectories");
      read_line("nth","%d",&nth);
      read_line("ntr","%d",&ntr);
      read_line("dtr_log","%d",&dtr_log);
      if (noms==0)
         read_line("dtr_ms","%d",&dtr_ms);
      else
         dtr_ms=0;
      read_line("dtr_cnfg","%d",&dtr_cnfg);

      error_root((append!=0)&&(nth!=0),1,"read_schedule [ym1.c]",
                 "Continuation run: nth must be equal to zero");

      ie=0;
      ie|=(nth<0);
      ie|=(ntr<1);
      ie|=(dtr_log<1);
      ie|=(dtr_log>dtr_cnfg);
      ie|=((dtr_cnfg%dtr_log)!=0);
      ie|=((nth%dtr_cnfg)!=0);
      ie|=((ntr%dtr_cnfg)!=0);

      if (noms==0)
      {
         ie|=(dtr_ms<dtr_log);
         ie|=(dtr_ms>dtr_cnfg);
         ie|=((dtr_ms%dtr_log)!=0);
         ie|=((dtr_cnfg%dtr_ms)!=0);
      }

      error_root(ie!=0,1,"read_schedule [ym1.c]",
                 "Improper value of nth,ntr,dtr_log,dtr_ms or dtr_cnfg");
   }

   MPI_Bcast(&nth,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ntr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dtr_log,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dtr_ms,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dtr_cnfg,1,MPI_INT,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      if (append)
      {
         ir=fread(istd,sizeof(stdint_t),3,fdat);
         error_root(ir!=3,1,"read_schedule [ym1.c]",
                    "Incorrect read count");

         if (endian==BIG_ENDIAN)
            bswap_int(3,istd);

         ie=0;
         ie|=(istd[0]!=(stdint_t)(dtr_log));
         ie|=(istd[1]!=(stdint_t)(dtr_ms));
         ie|=(istd[2]!=(stdint_t)(dtr_cnfg));

         error_root(ie!=0,1,"read_schedule [ym1.c]",
                    "Parameters do not match previous run");
      }
      else
      {
         istd[0]=(stdint_t)(dtr_log);
         istd[1]=(stdint_t)(dtr_ms);
         istd[2]=(stdint_t)(dtr_cnfg);

         if (endian==BIG_ENDIAN)
            bswap_int(3,istd);

         iw=fwrite(istd,sizeof(stdint_t),3,fdat);
         error_root(iw!=3,1,"read_schedule [ym1.c]",
                    "Incorrect write count");
      }
   }
}


static void read_wflow_parms(void)
{
   int nstep,dnms,ie,ir,iw;
   stdint_t istd[3];
   double eps,dstd[1];

   if (my_rank==0)
   {
      if (append)
      {
         ir=fread(istd,sizeof(stdint_t),1,fdat);
         error_root(ir!=1,1,"read_wflow_parms [ym1.c]",
                    "Incorrect read count");

         if (endian==BIG_ENDIAN)
            bswap_int(1,istd);

         error_root(istd[0]!=(stdint_t)(noms==0),1,"read_wflow_parms [ym1.c]",
                    "Attempt to mix measurement with other runs");
      }
      else
      {
         istd[0]=(stdint_t)(noms==0);

         if (endian==BIG_ENDIAN)
            bswap_int(1,istd);

         iw=fwrite(istd,sizeof(stdint_t),1,fdat);
         error_root(iw!=1,1,"read_wflow_parms [ym1.c]",
                    "Incorrect write count");
      }

      if (noms==0)
      {
         find_section("Wilson flow");
         read_line("integrator","%s",line);
         read_line("eps","%lf",&eps);
         read_line("nstep","%d",&nstep);
         read_line("dnms","%d",&dnms);

         if (strcmp(line,"EULER")==0)
            flint=0;
         else if (strcmp(line,"RK2")==0)
            flint=1;
         else if (strcmp(line,"RK3")==0)
            flint=2;
         else
            error_root(1,1,"read_wflow_parms [ym1.c]","Unkown integrator");

         error_root((dnms<1)||(nstep<dnms)||((nstep%dnms)!=0),1,
                    "read_wflow_parms [ym1.c]",
                    "nstep must be a multiple of dnms");
      }
      else
      {
         flint=0;
         eps=0.0;
         nstep=1;
         dnms=1;
      }
   }

   MPI_Bcast(&flint,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&eps,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&nstep,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dnms,1,MPI_INT,0,MPI_COMM_WORLD);

   file_head.dn=dnms;
   file_head.nn=nstep/dnms;
   file_head.tmax=N0;
   file_head.eps=eps;

   if ((my_rank==0)&&(noms==0))
   {
      if (append)
      {
         ir=fread(istd,sizeof(stdint_t),3,fdat);
         ir+=fread(dstd,sizeof(double),1,fdat);
         error_root(ir!=4,1,"read_wflow_parms [ym1.c]",
                    "Incorrect read count");

         if (endian==BIG_ENDIAN)
         {
            bswap_int(3,istd);
            bswap_double(1,dstd);
         }

         ie=0;
         ie|=(istd[0]!=(stdint_t)(flint));
         ie|=(istd[1]!=(stdint_t)(nstep));
         ie|=(istd[2]!=(stdint_t)(dnms));
         ie|=(dstd[0]!=eps);

         error_root(ie!=0,1,"read_wflow_parms [ym1.c]",
                    "Parameters do not match previous run");
      }
      else
      {
         istd[0]=(stdint_t)(flint);
         istd[1]=(stdint_t)(nstep);
         istd[2]=(stdint_t)(dnms);
         dstd[0]=eps;

         if (endian==BIG_ENDIAN)
         {
            bswap_int(3,istd);
            bswap_double(1,dstd);
         }

         iw=fwrite(istd,sizeof(stdint_t),3,fdat);
         iw+=fwrite(dstd,sizeof(double),1,fdat);
         error_root(iw!=4,1,"read_wflow_parms [ym1.c]",
                    "Incorrect write count");
      }
   }
}


static void read_infile(int argc,char *argv[])
{
   int ifile;

   if (my_rank==0)
   {
      flog=freopen("STARTUP_ERROR","w",stdout);

      ifile=find_opt(argc,argv,"-i");
      noloc=find_opt(argc,argv,"-noloc");
      noexp=find_opt(argc,argv,"-noexp");
      rmold=find_opt(argc,argv,"-rmold");
      noms=find_opt(argc,argv,"-noms");
      scnfg=find_opt(argc,argv,"-c");
      append=find_opt(argc,argv,"-a");
      norng=find_opt(argc,argv,"-norng");
      endian=endianness();

      error_root((ifile==0)||(ifile==(argc-1))||(scnfg==(argc-1))||
                 ((append!=0)&&(scnfg==0)),1,"read_infile [ym1.c]",
                 "Syntax: ym1 -i <filename> [-noloc] [-noexp] "
                 "[-rmold] [-noms] [-c <filename> [-a [-norng]]]");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [ym1.c]",
                 "Machine has unknown endianness");

      error_root((noexp)&&(noloc),1,"read_infile [ym1.c]",
		 "The concurrent use of -noloc and -noexp is not permitted");

      if (scnfg)
      {
         strncpy(cnfg,argv[scnfg+1],NAME_SIZE-1);
         cnfg[NAME_SIZE-1]='\0';
      }
      else
         cnfg[0]='\0';

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [ym1.c]",
                 "Unable to open input file");
   }

   MPI_Bcast(&noloc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noexp,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&rmold,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noms,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&scnfg,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&append,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&norng,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&endian,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      find_section("Random number generator");
      read_line("level","%d",&level);
      read_line("seed","%d",&seed);
   }

   MPI_Bcast(&level,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);

   read_dirs();
   setup_files();

   if (my_rank==0)
   {
      if (append)
         fdat=fopen(par_file,"rb");
      else
         fdat=fopen(par_file,"wb");

      error_root(fdat==NULL,1,"read_infile [ym1.c]",
                 "Unable to open parameter file");
   }

   read_lat_parms();
   read_bc_parms();
   read_hmc_parms();
   read_schedule();
   read_integrator();
   read_wflow_parms();

   if (my_rank==0)
   {
      fclose(fin);
      fclose(fdat);

      if (append==0)
         copy_file(par_file,par_save);
   }
}


static void check_old_log(int ic,int *nl,int *icnfg)
{
   int ir,isv;
   int np[4],bp[4];

   fend=fopen(log_file,"r");
   error_root(fend==NULL,1,"check_old_log [ym1.c]",
              "Unable to open log file");
   (*nl)=0;
   (*icnfg)=0;
   ir=1;
   isv=0;

   while (fgets(line,NAME_SIZE,fend)!=NULL)
   {
      if (strstr(line,"process grid")!=NULL)
      {
         ir&=(sscanf(line,"%dx%dx%dx%d process grid, %dx%dx%dx%d",
                     np,np+1,np+2,np+3,bp,bp+1,bp+2,bp+3)==8);

         ipgrd[0]=((np[0]!=NPROC0)||(np[1]!=NPROC1)||
                   (np[2]!=NPROC2)||(np[3]!=NPROC3));
         ipgrd[1]=((bp[0]!=NPROC0_BLK)||(bp[1]!=NPROC1_BLK)||
                   (bp[2]!=NPROC2_BLK)||(bp[3]!=NPROC3_BLK));
      }
      else if (strstr(line,"Trajectory no")!=NULL)
      {
         ir&=(sscanf(line,"Trajectory no %d",nl)==1);
         isv=0;
      }
      else if (strstr(line,"Configuration no")!=NULL)
      {
         ir&=(sscanf(line,"Configuration no %d",icnfg)==1);
         isv=1;
      }
   }

   fclose(fend);

   error_root(ir!=1,1,"check_old_log [ym1.c]","Incorrect read count");

   error_root(ic!=(*icnfg),1,"check_old_log [ym1.c]",
              "Continuation run:\n"
              "Initial configuration is not the last one of the previous run");

   error_root(isv==0,1,"check_old_log [ym1.c]",
              "Continuation run:\n"
              "The log file extends beyond the last configuration save");
}


static void check_old_dat(int nl)
{
   int nt;
   dat_t ndat;

   fdat=fopen(dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_dat [ym1.c]",
              "Unable to open data file");
   nt=0;

   while (read_dat(1,&ndat)==1)
      nt=ndat.nt;

   fclose(fdat);

   error_root(nt!=nl,1,"check_old_dat [ym1.c]",
              "Continuation run: Incomplete or too many data records");
}


static void check_old_msdat(int nl)
{
   int ic,ir,nt,pnt,dnt;

   fdat=fopen(msdat_file,"rb");
   error_root(fdat==NULL,1,"check_old_msdat [ym1.c]",
              "Unable to open data file");

   check_file_head();

   nt=0;
   dnt=0;
   pnt=0;

   for (ic=0;;ic++)
   {
      ir=read_data();

      if (ir==0)
      {
         error_root(ic==0,1,"check_old_msdat [ym1.c]",
                    "No data records found");
         break;
      }

      nt=data.nt;

      if (ic==1)
      {
         dnt=nt-pnt;
         error_root(dnt<1,1,"check_old_msdat [ym1.c]",
                    "Incorrect trajectory separation");
      }
      else if (ic>1)
         error_root(nt!=(pnt+dnt),1,"check_old_msdat [ym1.c]",
                    "Trajectory sequence is not equally spaced");

      pnt=nt;
   }

   fclose(fdat);

   error_root((nt!=nl)||((ic>1)&&(dnt!=dtr_ms)),1,
              "check_old_msdat [ym1.c]","Last trajectory numbers "
              "or the trajectory separations do not match");
}


static void check_files(int *nl,int *icnfg)
{
   int icmax,ic;

   ipgrd[0]=0;
   ipgrd[1]=0;

   if (my_rank==0)
   {
      if (noloc)
         error_root(cnfg[strlen(cnfg)-1]=='*',1,
                    "check_files [ym1.c]","Attempt to read an "
                    "imported configuration when -noloc is set");

      if (append)
      {
         error_root(strstr(cnfg,nbase)!=cnfg,1,"check_files [ym1.c]",
                    "Continuation run:\n"
                    "Run name does not match the previous one");
         error_root(sscanf(cnfg+strlen(nbase),"n%d",&ic)!=1,1,
                    "check_files [ym1.c]","Continuation run:\n"
                    "Unable to read configuration number from file name");

         check_old_log(ic,nl,icnfg);
         check_old_dat(*nl);
         if (noms==0)
            check_old_msdat(*nl);

         (*icnfg)+=1;
      }
      else
      {
         fin=fopen(log_file,"r");
         fdat=fopen(dat_file,"rb");

         if (noms==0)
            fend=fopen(msdat_file,"rb");
         else
            fend=NULL;

         error_root((fin!=NULL)||(fdat!=NULL)||(fend!=NULL),1,
                    "check_files [ym1.c]",
                    "Attempt to overwrite old *.log or *.dat file");

         if (noms==0)
         {
            fdat=fopen(msdat_file,"wb");
            error_root(fdat==NULL,1,"check_files [ym1.c]",
                       "Unable to open measurement data file");
            write_file_head();
            fclose(fdat);
         }

         (*nl)=0;
         (*icnfg)=1;
      }

      icmax=(*icnfg)+(ntr-nth)/dtr_cnfg;

      if (noloc==0)
         error_root(name_size("%s/%sn%d_%d",loc_dir,nbase,icmax,NPROC-1)>=
                    NAME_SIZE,1,"check_files [ym1.c]",
                    "loc_dir name is too long");

      if (noexp==0)
         error_root(name_size("%s/%sn%d",cnfg_dir,nbase,icmax)>=NAME_SIZE,1,
                    "check_files [ym1.c]","cnfg_dir name is too long");

      if (scnfg)
      {
         if (cnfg[strlen(cnfg)-1]=='*')
            error_root(name_size("%s/%s%d",loc_dir,cnfg,NPROC-1)>=NAME_SIZE,1,
                       "check_files [ym1.c]","loc_dir name is too long");
         else
            error_root(name_size("%s/%s",cnfg_dir,cnfg)>=NAME_SIZE,1,
                       "check_files [ym1.c]","cnfg_dir name is too long");
      }
   }

   MPI_Bcast(nl,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(icnfg,1,MPI_INT,0,MPI_COMM_WORLD);
}


static void init_rng(int icnfg)
{
   int ic;

   if (append)
   {
      if (cnfg[strlen(cnfg)-1]!='*')
      {
         if (norng)
            start_ranlux(level,seed^(icnfg-1));
         else
         {
            ic=import_ranlux(rng_file);
            error_root(ic!=(icnfg-1),1,"init_rng [ym1.c]",
                       "Configuration number mismatch (*.rng file)");
         }
      }
   }
   else
      start_ranlux(level,seed);
}


static void init_ud(void)
{
   char *p;

   if (scnfg)
   {
      if (cnfg[strlen(cnfg)-1]!='*')
      {
         sprintf(cnfg_file,"%s/%s",cnfg_dir,cnfg);
         import_cnfg(cnfg_file);
      }
      else
      {
         sprintf(line,"%s/%s",loc_dir,cnfg);
         p=line+strlen(line)-1;
         p[0]='\0';
         sprintf(cnfg_file,"%s_%d",line,my_rank);
         read_cnfg(cnfg_file);
      }
   }
   else
      random_ud();
}


static void store_ud(su3_dble *usv)
{
   su3_dble *udb;

   udb=udfld();
   cm3x3_assign(4*VOLUME,udb,usv);
}


static void recall_ud(su3_dble *usv)
{
   su3_dble *udb;

   udb=udfld();
   cm3x3_assign(4*VOLUME,usv,udb);
   set_flags(UPDATED_UD);
}


static void set_data(int nt)
{
   int in,dn,nn;
   double eps;

   data.nt=nt;
   dn=file_head.dn;
   nn=file_head.nn;
   eps=file_head.eps;

   for (in=0;in<nn;in++)
   {
      Wact[in]=plaq_action_slices(data.Wsl[in]);
      Yact[in]=ym_action_slices(data.Ysl[in]);
      Qtop[in]=tcharge_slices(data.Qsl[in]);

      if (flint==0)
         fwd_euler(dn,eps);
      else if (flint==1)
         fwd_rk2(dn,eps);
      else
         fwd_rk3(dn,eps);
   }

   Wact[in]=plaq_action_slices(data.Wsl[in]);
   Yact[in]=ym_action_slices(data.Ysl[in]);
   Qtop[in]=tcharge_slices(data.Qsl[in]);
}


static void print_info(int icnfg)
{
   int n;
   long ip;
   mdint_parms_t mdp;

   if (my_rank==0)
   {
      ip=ftell(flog);
      fclose(flog);

      if (ip==0L)
         remove("STARTUP_ERROR");

      if (append)
         flog=freopen(log_file,"a",stdout);
      else
         flog=freopen(log_file,"w",stdout);

      error_root(flog==NULL,1,"print_info [ym1.c]","Unable to open log file");

      if (append)
         printf("Continuation run, start from configuration %s\n\n",cnfg);
      else
      {
         printf("\n");
         printf("Simulation of the SU(3) gauge theory\n");
         printf("------------------------------------\n\n");

         if (scnfg)
            printf("New run, start from configuration %s\n\n",cnfg);
         else
            printf("New run, start from random configuration\n\n");

         printf("Using the HMC algorithm\n");
         printf("Program version %s\n",openQCD_RELEASE);
      }

      if (endian==LITTLE_ENDIAN)
         printf("The machine is little endian\n");
      else
         printf("The machine is big endian\n");
      if (noloc)
         printf("The local disks are not used\n");
      if (noexp)
         printf("The generated configurations are not exported\n");
      if (rmold)
         printf("Old configurations are deleted\n");
      printf("\n");

      if ((ipgrd[0]!=0)&&(ipgrd[1]!=0))
         printf("Process grid and process block size changed:\n");
      else if (ipgrd[0]!=0)
         printf("Process grid changed:\n");
      else if (ipgrd[1]!=0)
         printf("Process block size changed:\n");

      if ((append==0)||(ipgrd[0]!=0)||(ipgrd[1]!=0))
      {
         printf("%dx%dx%dx%d lattice, ",N0,N1,N2,N3);
         printf("%dx%dx%dx%d local lattice\n",L0,L1,L2,L3);
         printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
         printf("%dx%dx%dx%d process block size\n\n",
                NPROC0_BLK,NPROC1_BLK,NPROC2_BLK,NPROC3_BLK);
      }

      if (append==0)
      {
         n=fdigits(lat.beta);
         printf("beta = %.*f\n",IMAX(n,1),lat.beta);
         n=fdigits(lat.c0);
         printf("c0 = %.*f, ",IMAX(n,1),lat.c0);
         n=fdigits(lat.c1);
         printf("c1 = %.*f\n\n",IMAX(n,1),lat.c1);

         print_bc_parms(1);
      }

      printf("Random number generator:\n");

      if (append)
      {
         if (cnfg[strlen(cnfg)-1]!='*')
         {
            if (norng)
               printf("level = %d, seed = %d, effective seed = %d\n\n",
                      level,seed,seed^(icnfg-1));
            else
            {
               printf("State of ranlxs and ranlxd reset to the\n");
               printf("last exported state\n\n");
            }
         }
         else
         {
            printf("State of ranlxs and ranlxd read from\n");
            printf("initial field-configuration file\n\n");
         }
      }
      else
         printf("level = %d, seed = %d\n\n",level,seed);

      if (append)
      {
         printf("Trajectories:\n");
         printf("ntr = %d\n\n",ntr);
      }
      else
      {
         printf("Trajectories:\n");
         n=fdigits(hmc.tau);
         printf("tau = %.*f\n",IMAX(n,1),hmc.tau);

         mdp=mdint_parms(0);

         if (mdp.integrator==LPFR)
            printf("Leapfrog integrator\n");
         else if (mdp.integrator==OMF2)
            printf("2nd order OMF integrator with lambda = %.4e\n",mdp.lambda);
         else if (mdp.integrator==OMF4)
            printf("4th order OMF integrator\n");

         printf("Number of steps = %d\n\n",mdp.nstep);

         printf("nth = %d, ntr = %d\n",nth,ntr);

         if (noms)
         {
            printf("dtr_log = %d, dtr_cnfg = %d\n\n",
                   dtr_log,dtr_cnfg);
            printf("Wilson flow observables are not measured\n\n");
         }
         else
         {
            printf("dtr_log = %d, dtr_ms = %d, dtr_cnfg = %d\n\n",
                   dtr_log,dtr_ms,dtr_cnfg);
            printf("Online measurement of Wilson flow observables\n\n");

            printf("Wilson flow:\n");
            if (flint==0)
               printf("Euler integrator\n");
            else if (flint==1)
               printf("2nd order RK integrator\n");
            else
               printf("3rd order RK integrator\n");
            n=fdigits(file_head.eps);
            printf("eps = %.*f\n",IMAX(n,1),file_head.eps);
            printf("nstep = %d\n",file_head.dn*file_head.nn);
            printf("dnms = %d\n\n",file_head.dn);
         }
      }

      fflush(flog);
   }
}


static void print_log(dat_t *ndat)
{
   if (my_rank==0)
   {
      printf("Trajectory no %d\n",(*ndat).nt);
      printf("dH = %+.1e, ",(*ndat).dH);
      printf("iac = %d\n",(*ndat).iac);
      printf("Average plaquette = %.6f\n",(*ndat).avpl);
   }
}


static void save_dat(int n,double siac,double wtcyc,double wtall,dat_t *ndat)
{
   int iw;

   if (my_rank==0)
   {
      fdat=fopen(dat_file,"ab");
      error_root(fdat==NULL,1,"save_dat [ym1.c]",
                 "Unable to open data file");

      iw=write_dat(1,ndat);
      error_root(iw!=1,1,"save_dat [ym1.c]",
                 "Incorrect write count");
      fclose(fdat);

      printf("Acceptance rate = %1.2f\n",siac/(double)(n+1));
      printf("Time per trajectory = %.2e sec (average = %.2e sec)\n\n",
             wtcyc/(double)(dtr_log),wtall/(double)(n+1));
      fflush(flog);
   }
}


static void save_msdat(int n,double wtms,double wtmsall)
{
   int nms,in,dn,nn,din;
   double eps;

   if (my_rank==0)
   {
      fdat=fopen(msdat_file,"ab");
      error_root(fdat==NULL,1,"save_msdat [ym1.c]",
                 "Unable to open data file");
      write_data();
      fclose(fdat);

      nms=(n+1-nth)/dtr_ms+(nth>0);
      dn=file_head.dn;
      nn=file_head.nn;
      eps=file_head.eps;

      din=nn/10;
      if (din<1)
         din=1;

      printf("Measurement run:\n\n");

      for (in=0;in<=nn;in+=din)
         printf("n = %3d, t = %.2e, Wact = %.6e, Yact = %.6e, Q = % .2e\n",
                in*dn,eps*(double)(in*dn),Wact[in],Yact[in],Qtop[in]);

      printf("\n");
      printf("Configuration fully processed in %.2e sec ",wtms);
      printf("(average = %.2e sec)\n",wtmsall/(double)(nms));
      printf("Measured data saved\n\n");
      fflush(flog);
   }
}


static void save_cnfg(int icnfg)
{
   int ie;

   ie=query_flags(UD_PHASE_SET);
   if (ie==0)
      ie=check_bc(0.0)^0x1;
   error_root(ie!=0,1,"save_cnfg [ym1.c]",
              "Phase-modified field or unexpected boundary values");

   if (noloc==0)
   {
      sprintf(cnfg_file,"%s/%sn%d_%d",loc_dir,nbase,icnfg,my_rank);
      write_cnfg(cnfg_file);
   }

   if (noexp==0)
   {
      sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
      export_cnfg(cnfg_file);
   }

   if (my_rank==0)
   {
      if ((noloc==0)&&(noexp==0))
         printf("Configuration no %d saved on the local disks "
                "and exported\n\n",icnfg);
      else if (noloc==0)
         printf("Configuration no %d saved on the local disks\n\n",icnfg);
      else if (noexp==0)
         printf("Configuration no %d exported\n\n",icnfg);
   }
}


static void check_endflag(int *iend)
{
   if (my_rank==0)
   {
      fend=fopen(end_file,"r");

      if (fend!=NULL)
      {
         fclose(fend);
         remove(end_file);
         (*iend)=1;
         printf("End flag set, run stopped\n\n");
      }
      else
         (*iend)=0;
   }

   MPI_Bcast(iend,1,MPI_INT,0,MPI_COMM_WORLD);
}


static void remove_cnfg(int icnfg)
{
   if ((rmold)&&(icnfg>=1))
   {
      if (noloc==0)
      {
         sprintf(cnfg_file,"%s/%sn%d_%d",loc_dir,nbase,icnfg,my_rank);
         remove(cnfg_file);
      }

      if ((noexp==0)&&(my_rank==0))
      {
         sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
         remove(cnfg_file);
      }
   }
}


int main(int argc,char *argv[])
{
   int nl,icnfg;
   int n,iend,iac;
   double act0[2],act1[2],w0[2],w1[2],npl,siac;
   double wt1,wt2,wtcyc,wtall,wtms,wtmsall;
   su3_dble **usv;
   dat_t ndat;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile(argc,argv);
   if (noms==0)
      alloc_data();
   check_files(&nl,&icnfg);
   geometry();
   alloc_wud(1);
   if ((noms==0)&&(flint))
         alloc_wfd(1);

   print_info(icnfg);
   set_mdsteps();
   init_rng(icnfg);
   init_ud();

   if (bc_type()==0)
      npl=(double)(6*(N0-1)*N1)*(double)(N2*N3);
   else
      npl=(double)(6*N0*N1)*(double)(N2*N3);

   iend=0;
   siac=0.0;
   wtcyc=0.0;
   wtall=0.0;
   wtms=0.0;
   wtmsall=0.0;

   for (n=0;(iend==0)&&(n<ntr);n++)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      iac=run_hmc(act0,act1);

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      siac+=(double)(iac);
      wtcyc+=(wt2-wt1);

      if (((ntr-n-1)%dtr_log)==0)
      {
         w0[0]=(act1[0]-act0[0])+(act1[1]-act0[1]);
         w0[1]=plaq_wsum_dble(0)/npl;

         MPI_Reduce(w0,w1,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
         MPI_Bcast(w1,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

         ndat.nt=nl+n+1;
         ndat.iac=iac;
         ndat.dH=w1[0];
         ndat.avpl=w1[1];

         print_log(&ndat);
         wtall+=wtcyc;
         save_dat(n,siac,wtcyc,wtall,&ndat);
         wtcyc=0.0;

         if ((noms==0)&&((n+1)>=nth)&&(((ntr-n-1)%dtr_ms)==0))
         {
            MPI_Barrier(MPI_COMM_WORLD);
            wt1=MPI_Wtime();

            usv=reserve_wud(1);
            store_ud(usv[0]);
            set_data(nl+n+1);
            recall_ud(usv[0]);
            release_wud();

            MPI_Barrier(MPI_COMM_WORLD);
            wt2=MPI_Wtime();

            wtms=wt2-wt1;
            wtmsall+=wtms;
            save_msdat(n,wtms,wtmsall);
         }
      }

      if (((n+1)>=nth)&&(((ntr-n-1)%dtr_cnfg)==0))
      {
         save_cnfg(icnfg);
         export_ranlux(icnfg,rng_file);
         check_endflag(&iend);

         if (my_rank==0)
         {
            fflush(flog);
            copy_file(log_file,log_save);
            copy_file(dat_file,dat_save);
            if (noms==0)
               copy_file(msdat_file,msdat_save);
            copy_file(rng_file,rng_save);
         }

         remove_cnfg(icnfg-1);
         icnfg+=1;
      }
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
