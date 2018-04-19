
/*******************************************************************************
*
* File ms4.c
*
* Copyright (C) 2012, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of quark propagators.
*
* Syntax: ms4 -i <input file> [-noexp]
*
* For usage instructions see the file README.ms4.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "archive.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "dirac.h"
#include "sap.h"
#include "dfl.h"
#include "forces.h"
#include "version.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

#define MAX(n,m) \
   if ((n)<(m)) \
      (n)=(m)

static int my_rank,noexp,endian;
static int first,last,step;
static int level,seed,x0,nsrc;
static int *rlxs_state=NULL,*rlxd_state=NULL;
static double mus;

static char log_dir[NAME_SIZE],loc_dir[NAME_SIZE];
static char cnfg_dir[NAME_SIZE],sfld_dir[NAME_SIZE];
static char log_file[NAME_SIZE],log_save[NAME_SIZE],end_file[NAME_SIZE];
static char cnfg_file[NAME_SIZE],sfld_file[NAME_SIZE],nbase[NAME_SIZE];
static FILE *fin=NULL,*flog=NULL,*fend=NULL;

static lat_parms_t lat;
static bc_parms_t bcp;


static void read_dirs(void)
{
   if (my_rank==0)
   {
      find_section("Run name");
      read_line("name","%s",nbase);

      find_section("Directories");
      read_line("log_dir","%s",log_dir);

      if (noexp)
      {
         read_line("loc_dir","%s",loc_dir);
         cnfg_dir[0]='\0';
      }
      else
      {
         read_line("cnfg_dir","%s",cnfg_dir);
         loc_dir[0]='\0';
      }

      read_line("sfld_dir","%s",sfld_dir);

      find_section("Configurations");
      read_line("first","%d",&first);
      read_line("last","%d",&last);
      read_line("step","%d",&step);

      find_section("Random number generator");
      read_line("level","%d",&level);
      read_line("seed","%d",&seed);

      error_root((last<first)||(step<1)||(((last-first)%step)!=0),1,
                 "read_dirs [ms4.c]","Improper configuration range");
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(log_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(loc_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(sfld_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(&first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&step,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(&level,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);
}


static void setup_files(void)
{
   if (noexp)
      error_root(name_size("%s/%sn%d_%d",loc_dir,nbase,last,NPROC-1)>=NAME_SIZE,
                 1,"setup_files [ms4.c]","loc_dir name is too long");
   else
      error_root(name_size("%s/%sn%d",cnfg_dir,nbase,last)>=NAME_SIZE,
                 1,"setup_files [ms4.c]","cnfg_dir name is too long");

   check_dir_root(sfld_dir);
   error_root(name_size("%s/%sn%d.s%d",sfld_dir,nbase,last,nsrc-1)>=NAME_SIZE,
              1,"setup_files [ms4.c]","sfld_dir name is too long");

   check_dir_root(log_dir);
   error_root(name_size("%s/%s.ms4.log~",log_dir,nbase)>=NAME_SIZE,
              1,"setup_files [ms4.c]","log_dir name is too long");

   sprintf(log_file,"%s/%s.ms4.log",log_dir,nbase);
   sprintf(end_file,"%s/%s.ms4.end",log_dir,nbase);
   sprintf(log_save,"%s~",log_file);
}


static void read_lat_parms(void)
{
   double kappa,csw;

   if (my_rank==0)
   {
      find_section("Dirac operator");
      read_line("kappa","%lf",&kappa);
      read_line("mu","%lf",&mus);
      read_line("csw","%lf",&csw);

      find_section("Source fields");
      read_line("x0","%d",&x0);
      read_line("nsrc","%d",&nsrc);

      error_root((x0<0)||(x0>=N0),1,"read_lat_parms [ms4.c]",
                 "Specified time x0 is out of range");
      error_root(nsrc<1,1,"read_lat_parms [ms4.c]",
                 "The number of source fields must be at least 1");
   }

   MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&mus,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   MPI_Bcast(&x0,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nsrc,1,MPI_INT,0,MPI_COMM_WORLD);

   lat=set_lat_parms(0.0,1.0,1,&kappa,csw);
   set_sw_parms(sea_quark_mass(0));
}


static void read_bc_parms(void)
{
   int bc;
   double cF,cF_prime;
   double phi[2],phi_prime[2],theta[3];

   if (my_rank==0)
   {
      find_section("Boundary conditions");
      read_line("type","%d",&bc);

      error_root(((x0==0)&&(bc!=3))||((x0==(N0-1))&&(bc==0)),1,
                 "read_bc_parms [ms4.c]","Incompatible choice of boundary "
                 "conditions and source time");

      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;
      cF=1.0;
      cF_prime=1.0;

      if (bc==1)
         read_dprms("phi",2,phi);

      if ((bc==1)||(bc==2))
         read_dprms("phi'",2,phi_prime);

      if (bc!=3)
         read_line("cF","%lf",&cF);

      if (bc==2)
         read_line("cF'","%lf",&cF_prime);

      read_dprms("theta",3,theta);
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(phi,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(phi_prime,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(theta,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

   bcp=set_bc_parms(bc,1.0,1.0,cF,cF_prime,phi,phi_prime,theta);
}


static void read_sap_parms(void)
{
   int bs[4];

   if (my_rank==0)
   {
      find_section("SAP");
      read_line("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   set_sap_parms(bs,1,4,5);
}


static void read_dfl_parms(void)
{
   int bs[4],Ns;
   int ninv,nmr,ncy,nkv,nmx;
   double kappa,mu,res;

   if (my_rank==0)
   {
      find_section("Deflation subspace");
      read_line("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3);
      read_line("Ns","%d",&Ns);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_parms(bs,Ns);

   if (my_rank==0)
   {
      find_section("Deflation subspace generation");
      read_line("kappa","%lf",&kappa);
      read_line("mu","%lf",&mu);
      read_line("ninv","%d",&ninv);
      read_line("nmr","%d",&nmr);
      read_line("ncy","%d",&ncy);
   }

   MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&mu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&ninv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncy,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_gen_parms(kappa,mu,ninv,nmr,ncy);

   if (my_rank==0)
   {
      find_section("Deflation projection");
      read_line("nkv","%d",&nkv);
      read_line("nmx","%d",&nmx);
      read_line("res","%lf",&res);
   }

   MPI_Bcast(&nkv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   set_dfl_pro_parms(nkv,nmx,res);
}


static void read_solver(void)
{
   solver_parms_t sp;

   read_solver_parms(0);
   sp=solver_parms(0);

   if ((sp.solver==SAP_GCR)||(sp.solver==DFL_SAP_GCR))
      read_sap_parms();

   if (sp.solver==DFL_SAP_GCR)
      read_dfl_parms();
}


static void read_infile(int argc,char *argv[])
{
   int ifile;

   if (my_rank==0)
   {
      flog=freopen("STARTUP_ERROR","w",stdout);

      ifile=find_opt(argc,argv,"-i");
      endian=endianness();

      error_root((ifile==0)||(ifile==(argc-1)),1,"read_infile [ms4.c]",
                 "Syntax: ms4 -i <input file> [-noexp]");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [ms4.c]",
                 "Machine has unknown endianness");

      noexp=find_opt(argc,argv,"-noexp");

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [ms4.c]",
                 "Unable to open input file");
   }

   MPI_Bcast(&endian,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noexp,1,MPI_INT,0,MPI_COMM_WORLD);

   read_dirs();
   setup_files();
   read_lat_parms();
   read_bc_parms();
   read_solver();

   if (my_rank==0)
      fclose(fin);
}


static void check_files(void)
{
   if (my_rank==0)
   {
      fin=fopen(log_file,"r");
      error_root(fin!=NULL,1,"check_files [ms4.c]",
                 "Attempt to overwrite old *.log file");
   }
}


static void print_info(void)
{
   int isap,idfl,n;
   long ip;

   if (my_rank==0)
   {
      ip=ftell(flog);
      fclose(flog);

      if (ip==0L)
         remove("STARTUP_ERROR");

      flog=freopen(log_file,"w",stdout);
      error_root(flog==NULL,1,"print_info [ms4.c]","Unable to open log file");
      printf("\n");

      printf("Computation of quark propagators\n");
      printf("--------------------------------\n\n");

      printf("Program version %s\n",openQCD_RELEASE);

      if (endian==LITTLE_ENDIAN)
         printf("The machine is little endian\n");
      else
         printf("The machine is big endian\n");
      if (noexp)
         printf("Configurations are read in imported file format\n\n");
      else
         printf("Configurations are read in exported file format\n\n");

      printf("%dx%dx%dx%d lattice, ",N0,N1,N2,N3);
      printf("%dx%dx%dx%d local lattice\n",L0,L1,L2,L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d process block size\n",
             NPROC0_BLK,NPROC1_BLK,NPROC2_BLK,NPROC3_BLK);
      printf("SF boundary conditions on the quark fields\n\n");

      printf("Random number generator:\n");
      printf("level = %d, seed = %d\n\n",level,seed);

      printf("Dirac operator:\n");
      n=fdigits(lat.kappa[0]);
      printf("kappa = %.*f\n",IMAX(n,6),lat.kappa[0]);
      n=fdigits(mus);
      printf("mu = %.*f\n",IMAX(n,1),mus);
      n=fdigits(lat.csw);
      printf("csw = %.*f\n\n",IMAX(n,1),lat.csw);
      print_bc_parms(2);

      printf("Source fields:\n");
      printf("x0 = %d\n",x0);
      printf("nsrc = %d\n\n",nsrc);

      print_solver_parms(&isap,&idfl);

      if (isap)
         print_sap_parms(0);

      if (idfl)
         print_dfl_parms(0);

      printf("Configurations no %d -> %d in steps of %d\n\n",
             first,last,step);
      fflush(flog);
   }
}


static void dfl_wsize(int *nws,int *nwv,int *nwvd)
{
   dfl_parms_t dp;
   dfl_pro_parms_t dpp;

   dp=dfl_parms();
   dpp=dfl_pro_parms();

   MAX(*nws,dp.Ns+2);
   MAX(*nwv,2*dpp.nkv+2);
   MAX(*nwvd,4);
}


static void wsize(int *nws,int *nwsd,int *nwv,int *nwvd)
{
   int nsd;
   solver_parms_t sp;

   (*nws)=0;
   (*nwsd)=0;
   (*nwv)=0;
   (*nwvd)=0;

   sp=solver_parms(0);
   nsd=2;

   if (sp.solver==CGNE)
   {
      MAX(*nws,5);
      MAX(*nwsd,nsd+3);
   }
   else if (sp.solver==SAP_GCR)
   {
      MAX(*nws,2*sp.nkv+1);
      MAX(*nwsd,nsd+2);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      MAX(*nws,2*sp.nkv+2);
      MAX(*nwsd,nsd+4);
      dfl_wsize(nws,nwv,nwvd);
   }
   else
      error_root(1,1,"wsize [ms4.c]",
                 "Unknown or unsupported solver");
}


static void random_source(spinor_dble *eta)
{
   int y0,iy,ix;

   set_sd2zero(VOLUME,eta);
   y0=x0-cpr[0]*L0;

   if ((y0>=0)&&(y0<L0))
   {
      for (iy=0;iy<(L1*L2*L3);iy++)
      {
         ix=ipt[iy+y0*L1*L2*L3];
         random_sd(1,eta+ix,1.0);
      }
   }
}


static void solve_dirac(spinor_dble *eta,spinor_dble *psi,int *status)
{
   solver_parms_t sp;
   sap_parms_t sap;

   sp=solver_parms(0);

   if (sp.solver==CGNE)
   {
      mulg5_dble(VOLUME,eta);

      tmcg(sp.nmx,sp.res,mus,eta,eta,status);

      error_root(status[0]<0,1,"solve_dirac [ms4.c]",
                 "CGNE solver failed (status = %d)",status[0]);

      Dw_dble(-mus,eta,psi);
      mulg5_dble(VOLUME,psi);
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      sap_gcr(sp.nkv,sp.nmx,sp.res,mus,eta,psi,status);

      error_root(status[0]<0,1,"solve_dirac [ms4.c]",
                 "SAP_GCR solver failed (status = %d)",status[0]);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mus,eta,psi,status);

      error_root((status[0]<0)||(status[1]<0),1,
                 "solve_dirac [ms4.c]","DFL_SAP_GCR solver failed "
                 "(status = %d,%d,%d)",status[0],status[1],status[2]);
   }
   else
      error_root(1,1,"solve_dirac [ms4.c]",
                 "Unknown or unsupported solver");
}


static void propagator(int nc,int *status,double *time)
{
   int isrc,l,stat[3];
   double wt1,wt2,wtsum;
   spinor_dble *eta,*psi,**wsd;

   wsd=reserve_wsd(2);
   eta=wsd[0];
   psi=wsd[1];
   wtsum=0.0;

   for (l=0;l<3;l++)
   {
      status[l]=0;
      stat[l]=0;
   }

   for (isrc=0;isrc<nsrc;isrc++)
   {
      random_source(eta);

      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      solve_dirac(eta,psi,stat);

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      wtsum+=(wt2-wt1);

      for (l=0;l<2;l++)
         status[l]+=stat[l];

      status[2]+=(stat[2]!=0);

      sprintf(sfld_file,"%s/%sn%d.s%d",sfld_dir,nbase,nc,isrc);
      export_sfld(sfld_file,psi);
   }

   for (l=0;l<2;l++)
      status[l]=(status[l]+(nsrc/2))/nsrc;

   (*time)=wtsum/(double)(nsrc);

   release_wsd();
}


static void save_ranlux(void)
{
   int nlxs,nlxd;

   if (rlxs_state==NULL)
   {
      nlxs=rlxs_size();
      nlxd=rlxd_size();

      rlxs_state=malloc((nlxs+nlxd)*sizeof(int));
      rlxd_state=rlxs_state+nlxs;

      error(rlxs_state==NULL,1,"save_ranlux [ms4.c]",
            "Unable to allocate state arrays");
   }

   rlxs_get(rlxs_state);
   rlxd_get(rlxd_state);
}


static void restore_ranlux(void)
{
   rlxs_reset(rlxs_state);
   rlxd_reset(rlxd_state);
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


int main(int argc,char *argv[])
{
   int nc,iend,status[3];
   int nws,nwsd,nwv,nwvd,n;
   double wt1,wt2,wtavg;
   double wts,wtsavg;
   dfl_parms_t dfl;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile(argc,argv);
   check_files();
   print_info();
   dfl=dfl_parms();

   start_ranlux(level,seed);
   geometry();

   wsize(&nws,&nwsd,&nwv,&nwvd);
   alloc_ws(nws);
   alloc_wsd(nwsd);
   alloc_wv(nwv);
   alloc_wvd(nwvd);

   iend=0;
   wtavg=0.0;
   wtsavg=0.0;

   for (nc=first;(iend==0)&&(nc<=last);nc+=step)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      if (my_rank==0)
         printf("Configuration no %d\n",nc);

      if (noexp)
      {
         save_ranlux();
         sprintf(cnfg_file,"%s/%sn%d_%d",loc_dir,nbase,nc,my_rank);
         read_cnfg(cnfg_file);
         restore_ranlux();
      }
      else
      {
         sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,nc);
         import_cnfg(cnfg_file);
      }

      set_ud_phase();

      if (dfl.Ns)
      {
         dfl_modes(status);
         error_root(status[0]<0,1,"main [ms4.c]",
                    "Deflation subspace generation failed (status = %d)",
                    status[0]);

         if (my_rank==0)
            printf("Deflation subspace generation: status = %d\n",status[0]);
      }

      propagator(nc,status,&wts);
      wtsavg+=wts;

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      wtavg+=(wt2-wt1);

      if (my_rank==0)
      {
         printf("Computation of propagator completed\n");

         if (dfl.Ns)
         {
            printf("status = %d,%d",status[0],status[1]);

            if (status[2])
               printf(" (no of subspace regenerations = %d)\n",status[2]);
            else
               printf("\n");
         }
         else
            printf("status = %d\n",status[0]);

         n=(nc-first)/step+1;

         printf("Dirac equation: ");
         printf("%.2e sec per solution (average %.2e sec)\n",
                wts,wtsavg/(double)(n));
         printf("Configuration no %d fully processed in %.2e sec ",
                nc,wt2-wt1);
         printf("(average = %.2e sec)\n\n",wtavg/(double)(n));

         fflush(flog);
         copy_file(log_file,log_save);
      }

      check_endflag(&iend);
   }

   if (my_rank==0)
   {
      fflush(flog);
      copy_file(log_file,log_save);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
