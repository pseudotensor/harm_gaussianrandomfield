// icc -O3 -o inits inits.c

/* 
 * This code generates numbers that can be used as initial
 * conditions for velocity in an isothermal MHD turbulence run.
 *
 * The output, sent to stdout, consists of three columns:
 * vx,vy,vz.  Columns are in C-order, i.e. [0][0][0], [0][0][1],...  
 *
 * The velocity field is specified by the following parameters:
 *
 * nx,ny,nz: size of arrays
 * amp: rms velocity amplitude.
 * nodiv: (1 or 0) (don't) include the compressive part
 *    of the velocity field.  If 1, the field is guaranteed
 *    to be divergence-free to truncation error, *not* machine
 *    precision.
 * kmin: 
 * kmax: P(k) is nonzero only for kmin < k < kmax.
 *
 * In addition, a function P(k) gives the power spectrum.  
 * P(k) \sim k^{-11/3} is Kolmogorov.
 *       
 * ffts are done with the fftw package.
 * assumes dx = dy = dz
 *
 * cfg 10-1-01
 *
 *
 * modified to correctly deal with random number generation ordering when changing resolution so that can do convergence testing
 *
 * modified to study randomized field
 *
 *
 */




//#include "global.h"

#define TESTPOT 0
// which test potential to try 0=none


#define DOLOOP 1
// whether to do loop over realizations


#define OUTPUTTYPE 1
// 0: output files as below
// 1: no data files outputted, just output of statistics for each realization


#define WHICHVARIABLE 1
// 0: velocity
// 1: vector potential


#define REALTYPE 0 // output type
// 0: float
// 1: double
#define READIN 1
// 0: create dump file (no header)
// 1: create readin-velocity file for chandran() (ignore REALTYPE, just output double)

// per CPU size, NOT total size


#include <time.h>
#include <stdio.h>
#include <math.h>
#include <rfftw.h>

//#define MYRAND ((double)rand()/(double)RAND_MAX)
#define MYRAND ranc(0)
int nx,ny,nz ;
double amp ;
int nodiv ;
double kmin,kmax ;
long rannum;
double INDEX,EQK;

int main(int argc, char *argv[])
{
	int ijk,i,j,k ;
	int di,dj,dk;
	int l;
	int myid;
	int ncpux1,ncpux2,ncpux3,numprocs,mycpupos[3+1];
	double ph,a,a_re,a_im,ranc(int seed),K,return_K(int i, int j, int k),
	  rand_P(double K),return_kx(int i, int j, int k),
	  return_ky(int i, int j, int k),return_kz(int i, int j, int k),
	  kx,ky,kz,kdotv_re,kdotv_im,rms,rmstot ;
	int condition_k(double k);
	double rand_ph(double k);
	static fftw_complex *kvx,*kvy,*kvz ;
	static fftw_real *vx,*vy,*vz ;
	static fftw_real *vx2,*vy2,*vz2 ;

	static fftw_real *Bx,*By,*Bz ;
	static fftw_real *rho,*mofr,*press;
	double *rvsr,*rhovsr,*mofrvsr,*pressvsr;

	static rfftwnd_plan p, pinv ;
	rfftwnd_plan rfftw3d_create_plan(int nx, int ny, int nz, 
		fftw_direction dir, int flags) ;

	char filename[200];
	double dumlf;
	float duml;
	int dumi;
	double ke,ke1,ke2,ke3,rho0,Lx,Ly,Lz,bx0,by0,bz0;
	double va,vxa,vya,vza;
	double myrms=0,rmsvx=0;
	int seed;
	int argnum;
	char basefilename[200];
	int n1,n2,n3;
	FILE **files;

	double bsq,betaact,beta,norm,pmax,bsqmax;
	double pavg,petot,bsqavg;
	double helicity,betot,alphahel;

	FILE*averydata;
	double x,y,z,radius;
	int betanormtype;
	double stellarradius;
	double INDEXI,INDEXF;

	FILE *outputfile;


	// realizations
	int par1,par2,par3,par4,par5,par6,par7;
	int par1i,par2i,par3i,par4i,par5i,par6i,par7i;
	int par1f,par2f,par3f,par4f,par5f,par6f,par7f;
	int par1d,par2d,par3d,par4d,par5d,par6d,par7d;

	double dx,dy,dz,dV,Vol;

	int STENCILTYPE,SUBSTENCILTYPE;

	double Axavg,Ayavg,Azavg;

	double divB1,divB2,divB3,divB,divBnorm,divBavg,divBmax,divBavgnorm,divBmaxnorm,divBsum,divBnormsum,divBavgnorm2;
	int divBimax,divBjmax,divBkmax,divBimaxnorm,divBjmaxnorm,divBkmaxnorm;
	double divBnormalization;

	double divA1,divA2,divA3,divA,divAnorm,divAavg,divAmax,divAavgnorm,divAmaxnorm,divAsum,divAnormsum,divAavgnorm2;
	int divAimax,divAjmax,divAkmax,divAimaxnorm,divAjmaxnorm,divAkmaxnorm;
	double divAnormalization;

	rannum=0;

	if(argc!=7+1){
	  fprintf(stderr,"wrong number of args: got=%d expected=%d\n",argc-1,7);
	  fprintf(stderr,"./inits ncpux1 ncpux2 ncpux3 n1 n2 n3 basefilename\n");
	  exit(1);
	}
       

	argnum=1;

	ncpux1=atoi(argv[argnum++]);
	ncpux2=atoi(argv[argnum++]);
	ncpux3=atoi(argv[argnum++]);
	numprocs=ncpux1*ncpux2*ncpux3;

	n1 = atoi(argv[argnum++]);
	n2 = atoi(argv[argnum++]);
	n3 = atoi(argv[argnum++]);
	
	nx = n1*ncpux1 ;
	ny = n2*ncpux2 ;
	nz = n3*ncpux3 ;


	// allocate files for multiple-CPU files
	
	files=(FILE**)malloc(sizeof(FILE*)*ncpux1*ncpux2*ncpux3);
	if(files==NULL){
	  fprintf(stderr,"cannot allocate file name list\n");
	  exit(1);
	}
	strcpy(basefilename,argv[argnum++]);



	/* make some space for the transforms */
	kvx = (fftw_complex *) 	calloc(nx*ny*(nz/2+2) , sizeof(fftw_complex)) ;
	kvy = (fftw_complex *) 	calloc(nx*ny*(nz/2+2) , sizeof(fftw_complex)) ;
	kvz = (fftw_complex *) 	calloc(nx*ny*(nz/2+2) , sizeof(fftw_complex)) ;

	vx = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	vy = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	vz = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;

	vx2 = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	vy2 = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	vz2 = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;

	Bx = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	By = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	Bz = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;

	rho = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	mofr = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
        press = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;

	if(
	   (kvx==NULL)||(kvy==NULL)||(kvz==NULL)||
	   (vx==NULL)||(vy==NULL)||(vz==NULL)||
	   (vx2==NULL)||(vy2==NULL)||(vz2==NULL)||
	   (Bx==NULL)||(By==NULL)||(Bz==NULL)||
	   (rho==NULL)||(mofr==NULL)||(press==NULL)
	   )
	  {
	    fprintf(stderr,"Couldn't allocate memory: %d %d %d \n %d %d %d \n %d %d %d \n %d %d %d \n %d %d %d\n",kvx,kvy,kvz,vx,vy,vz,vx2,vy2,vz2,Bx,By,Bz,rho,mofr,press);
	    exit(1);
	  }
	/* make an fftw plan; don't measure performance */
	// fast start but can be slower to run many similar plans
	//	pinv = rfftw3d_create_plan(nx,ny,nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE) ;
	// slow to start but faster if doing many similar plans
	pinv = rfftw3d_create_plan(nx,ny,nz, FFTW_COMPLEX_TO_REAL, FFTW_MEASURE) ;


#if(WHICHVARIABLE==1)
	// need size of box before reading in data
	// make slightly bigger than radius of 1 so reach outside star
	Lx=2.2;
	Ly=2.2;
	Lz=2.2;
	Vol=Lx*Ly*Lz;

	stellarradius=1.0; // should be consistent with data read in


	// read in rho,mofr, and press as functions of radius
#define NAVERY 23850

	rvsr=(double*)malloc(NAVERY*sizeof(double));
	rhovsr=(double*)malloc(NAVERY*sizeof(double));
	mofrvsr=(double*)malloc(NAVERY*sizeof(double));
	pressvsr=(double*)malloc(NAVERY*sizeof(double));
	if((rvsr==NULL)||(rhovsr==NULL)||(mofrvsr==NULL)||(pressvsr==NULL)){
	  fprintf(stderr,"couldn't allocate memory for rhovsr mofrvsr or pressvsr\n");
	  exit(1);
	}
	fprintf(stderr,"Open and read avery's file...\n");
	// open and read avery's data file
	averydata=fopen("A0III_1.4","rt");
	if(averydata==NULL){
	  fprintf(stderr,"couldn't open avery's file\n");
	  exit(1);
	}
#define LINESTOSKIP 9

	for(i=1;i<=LINESTOSKIP;i++) while(fgetc(averydata)!='\n'); // skip lines
	for(i=0;i<NAVERY;i++){
	  fscanf(averydata,"%lf %lf %lf %lf\n",&rvsr[i],&mofrvsr[i],&rhovsr[i],&pressvsr[i]);
	  while(fgetc(averydata)!='\n'); // skip rest of line
	}
	fclose(averydata);
	fprintf(stderr,"Closing avery's file.\n");

	fprintf(stderr,"Interpolate avery's data...\n");
	// interpolate avery data and assign to 3D grid
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;

	  // zone centered where -Lx/2,-Ly/2,-Lz/2 are edges
	  dx=Lx/(double)nx;
	  dy=Ly/(double)ny;
	  dz=Lz/(double)nz;
	  dV=dx*dy*dz;

	  // centered x,y,z
	  x=-Lx/2+(i+0.5)*dx;
	  y=-Ly/2+(j+0.5)*dy;
	  z=-Lz/2+(k+0.5)*dz;

	  radius=sqrt(x*x+y*y+z*z);

	  // lookup table
	  l=0;
	  while(rvsr[l]<radius){
	    l++;
	    if(l>=NAVERY) break;
	  }
	  // assign interpolated value
	  if(l>=NAVERY-1){
	    rho[ijk]=0.0;
	    mofr[ijk]=0.0;
	    press[ijk]=0.0;
	  }
	  else{
	    // linear interpolation: places values at center of zones
	    rho[ijk]=rhovsr[l]+(rhovsr[l+1]-rhovsr[l])/(rvsr[l+1]-rvsr[l])*(radius-rvsr[l]);
	    mofr[ijk]=mofrvsr[l]+(mofrvsr[l+1]-mofrvsr[l])/(rvsr[l+1]-rvsr[l])*(radius-rvsr[l]);
	    press[ijk]=pressvsr[l]+(pressvsr[l+1]-pressvsr[l])/(rvsr[l+1]-rvsr[l])*(radius-rvsr[l]);
	  }
	}
	fprintf(stderr,"Done interplating avery's data.\n");

#if(0) // check data
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  fprintf(stderr,"%d %d %d : %21.15g %21.15g %21.15g\n",i,j,k,mofr[ijk],rho[ijk],press[ijk]);
	}

#endif


#endif

	//////////////////////////////////////////////////
	//
	// now we can loop over possible realizations
	//
	//////////////////////////////////////////////////

	outputfile=fopen(basefilename,"wt");
	if(outputfile==NULL){
	  fprintf(stderr,"couldn't open outputfile\n");
	  exit(1);
	}
	// output header
	fprintf(outputfile,"#N-INDEX,INDEX,seed,EQK,kmin,kmax,nodiv,betanormtype");
	fprintf(outputfile,",betaact,helicity,betot,alphahel");
	fprintf(outputfile,",divBavg,divBavgnorm,divBmax,divBmaxnorm,divBavgnorm2");
	fprintf(outputfile,",divAavg,divAavgnorm,divAmax,divAmaxnorm,divAavgnorm2");
	fprintf(outputfile,"\n");
	fflush(outputfile);


	

/* this returns the power spectrum */
#if(WHICHVARIABLE==0)
	rho0=1.0; // just for ke
	bx0=3.0; // Watson bx=10,3 (new best)
	by0=bz0=0.0;
	Lx=1.0; // just for ke
	Ly=1.0;  // just for ke
	Lz=1.0;  // just for ke
	Vol=Lx*Ly*Lz;
	dx=Lx/(double)nx;
	dy=Ly/(double)ny;
	dz=Lz/(double)nz;
	dV=dx*dy*dz;
	

	// Watson
	// sub-alfvenic
	// subsonic (at best mach=1)
	// many waves within box, not just 1
	// rms velocity along field larger than rms velocity perp field (want other way)

	EQK=2.0;
	amp = 4.0; // Watson amp=4
	nodiv = 1 ;	
	kmin = 2.0 ; //  Watson kmin=2	
	kmax = 32.0 ; // Watson kmax=31
	INDEX=-(11./3.);
	INDEXI=INDEXF=INDEX;
	seed=70;

	par1i=par1f=0; par1d=1;
	par2i=par2f=seed; par2d=1;
	par3i=par3f=(int)EQK; par3d=1;
	par4i=par4f=(int)kmin; par4d=1;
	par5i=par5f=(int)kmax; par5d=1;
	par6i=par6f=nodiv; par6d=1;
	par7i=par7f=0; par7d=1; // not used

	STENCILTYPE=0;
	SUBSTENCILTYPE=0;

#endif
#if(WHICHVARIABLE==1)

	// STENCILTYPE==1 and SUBSTENCIL==0 preserved divB=0 to machine precision
	STENCILTYPE=1;
	SUBSTENCILTYPE=0;

#if(0)
	nodiv = 0 ;	
	// kmin*L actually, where L is Lx Ly or Lz for each direction.
	kmin = 0.0 ;
	// kmax*L actually
	//kmax = 32.0 ;
	kmax = 16.0 ;
	INDEX=-10.0;
	seed=3;
	EQK=0.0;
	betanormtype=0; // 0=use peaks 1=use averages
#else
	nodiv = 1 ;	
	// kmin*L actually, where L is Lx Ly or Lz for each direction.
	kmin = 1.0 ;
	// kmax*L actually
	//kmax = 32.0 ;
	kmax = 4.0 ;
	INDEX=-8.18181818181818;
	seed=8;
	EQK=0.0;
	betanormtype=1; // 0=use peaks 1=use averages
#endif

	beta = 100.0; 

#if(DOLOOP)
	INDEXI=-10.0; INDEXF=10.0;
	par1i=0; par1f=11; par1d=1;
	par2i=1; par2f=20; par2d=1;
	par3i=0; par3f=3; par3d=1;
	par4i=0; par4f=3; par4d=1;
	par5i=4; par5f=255; par5d=4;
	par6i=0; par6f=1; par6d=1;
	par7i=0; par7f=1; par7d=1;
#else
	INDEXI=INDEXF=INDEX;
	par1i=par1f=0; par1d=1; // set by INDEX
	par2i=par2f=seed; par2d=1;
	par3i=par3f=(int)EQK; par3d=1;
	par4i=par4f=(int)kmin; par4d=1;
	par5i=par5f=(int)kmax; par5d=2;
	par6i=par6f=nodiv; par6d=1;
	par7i=par7f=betanormtype; par7d=1;
#endif




#endif

	

	for(par1=par1i;par1<=par1f;par1+=par1d){
	for(par2=par2i;par2<=par2f;par2+=par2d){
	for(par3=par3i;par3<=par3f;par3+=par3d){
	for(par4=par4i;par4<=par4f;par4+=par4d){
        for(par5=par5i;par5<=par5f;par5*=par5d){ // *=
	for(par6=par6i;par6<=par6f;par6+=par6d){
	for(par7=par7i;par7<=par7f;par7+=par7d){

	  // power index is a float
	  if(par1f==par1i) INDEX=INDEXI;
	  else INDEX=INDEXI+(INDEXF-INDEXI)/(par1f-par1i)*(par1-par1i);
	  seed=par2;
	  EQK=(double)par3;
	  kmin=(double)par4;
	  kmax=(double)par5;
	  nodiv=par6;
	  betanormtype=par7;


	  fprintf(stderr,"Generating power spectrum\n");
	  fprintf(stderr,"index=%g seed=%d eqk=%g kmin=%g kmax=%g nodiv=%d betanormtype=%d\n",INDEX,seed,EQK,kmin,kmax,nodiv,betanormtype); fflush(stderr);


	// check kmax for resolution
	if((nx/2<=kmax)||(ny/2<=kmax)||(nz/2<=kmax)){
	  fprintf(stderr,"kmax chosen poorly, need N>%d=2*kmax\n",(int)(2*kmax));
	  if(DOLOOP==0) exit(1);
	  else continue; // go on to next one
	}




	//	myrms=2.31029310017193;
	//rmsvx=0;
	//seed=-1;
	//	while(fabs(myrms-rmsvx)>1E-6){
	//	if(1)
	//	{
	  //seed++;
	  //	  seed=1;
	//ranc(1); // new Watson alternative seed
	//	ranc(121); // Original Watson seed
	//ranc(119); // Original Watson seed
	//	ranc(119);
	//ranc(7);
	//srand(7);

	// set first random number
	ranc(seed);

	/* fill in k-space arrays */
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz/2+1;k++) {
		K = return_K(i,j,k) ;
		//if(K==kmax) fprintf(stderr,"i=%d j=%d k=%d K=%g\n",i,j,k,K);

		ph=rand_ph(K);/* random phase */
		a = rand_P(K) ;

		a_re = a*cos(ph) ;
		a_im = a*sin(ph) ;

		ijk = k + (nz/2 + 1)*(j + ny*i) ;
		kvx[ijk].re = a_re ;
		kvx[ijk].im = a_im ;
	}
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz/2+1;k++) {
		K = return_K(i,j,k) ;

		ph=rand_ph(K);/* random phase */
		a = rand_P(K) ;

		a_re = a*cos(ph) ;
		a_im = a*sin(ph) ;

		ijk = k + (nz/2 + 1)*(j + ny*i) ;
		kvy[ijk].re = a_re ;
		kvy[ijk].im = a_im ;
	}
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz/2+1;k++) {
		K = return_K(i,j,k) ;

		ph=rand_ph(K);/* random phase */
		a = rand_P(K) ;

		a_re = a*cos(ph) ;
		a_im = a*sin(ph) ;

		ijk = k + (nz/2 + 1)*(j + ny*i) ;
		kvz[ijk].re = a_re ;
		kvz[ijk].im = a_im ;
	}

	/* for testing */
#if 0
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz/2+1;k++) {
		ijk = k + (nz/2 + 1)*(j + ny*i) ;
		ky = return_ky(i,j,k) ;
		kvx[ijk].re = exp(-0.5*ky*ky/(2.*2.)) ;
	}
	i = 0 ; j = 15 ; k = 0 ;
	ijk = k + (nz/2 + 1)*(j + ny*i) ;
	kvx[ijk].im = 1. ;
#endif

	/* project out divvy part of velocity field if desired */
	if(nodiv) {
		for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
		for(k=0;k<nz/2+1;k++) {
			kx = return_kx(i,j,k) ;
			ky = return_ky(i,j,k) ;
			kz = return_kz(i,j,k) ;
			K = return_K(i,j,k) ;

			ijk = k + (nz/2 + 1)*(j + ny*i) ;
			kdotv_re = kx*kvx[ijk].re +
				ky*kvy[ijk].re +
				kz*kvz[ijk].re ;
			kdotv_im = kx*kvx[ijk].im +
				ky*kvy[ijk].im +
				kz*kvz[ijk].im ;

			if(K != 0.) {
				kvx[ijk].re -= kx*kdotv_re/(K*K) ;
				kvx[ijk].im -= kx*kdotv_im/(K*K) ;

				kvy[ijk].re -= ky*kdotv_re/(K*K) ;
				kvy[ijk].im -= ky*kdotv_im/(K*K) ;

				kvz[ijk].re -= kz*kdotv_re/(K*K) ;
				kvz[ijk].im -= kz*kdotv_im/(K*K) ;
			}
		}
	}

	/* inverse transform */
	rfftwnd_one_complex_to_real(pinv, kvx, vx) ;
	rfftwnd_one_complex_to_real(pinv, kvy, vy) ;
	rfftwnd_one_complex_to_real(pinv, kvz, vz) ;


	fprintf(stderr,"Calculate and Output interesting things\n");
	


#if(WHICHVARIABLE==0)
	// normalize velocity by rms
	rms = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
		ijk = k + nz*(j + ny*i) ;
		rms += vx[ijk]*vx[ijk] + 
			vy[ijk]*vy[ijk] + 
			vz[ijk]*vz[ijk] ;
	}
	rms = sqrt(rms/(nx*ny*nz)) ;

	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
		ijk = k + nz*(j + ny*i) ;
		vx[ijk]*=amp/rms;
		vy[ijk]*=amp/rms;
		vz[ijk]*=amp/rms;
	}


#endif



#if(WHICHVARIABLE==0)
	//////////////////////////
	// diagnostics/statistics of data
	//
	va=vxa=vya=vza=0.0;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
		ijk = k + nz*(j + ny*i) ;
		vxa += vx[ijk];
		vya += vy[ijk];
		vza += vz[ijk];
	}
	vxa = vxa/(nx*ny*nz) ;
	vya = vya/(nx*ny*nz) ;
	vza = vza/(nx*ny*nz) ;


	rms = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
		ijk = k + nz*(j + ny*i) ;
		rms += (vx[ijk]-vxa)*(vx[ijk]-vxa) + 
		  (vy[ijk]-vya)*(vy[ijk]-vya) + 
		  (vz[ijk]-vza)*(vz[ijk]-vza) ;
	}
	rmstot = rms = sqrt(rms/(nx*ny*nz)) ;

	fprintf(stdout,"finalrms: %21.15g avg: %21.15g \n",rms,va); fflush(stdout);
	rms = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
		ijk = k + nz*(j + ny*i) ;
		rms += (vx[ijk]-vxa)*(vx[ijk]-vxa);
	}
	rms = sqrt(rms/(nx*ny*nz)) ;

	fprintf(stdout,"finalrms(vx): %21.15g avg: %21.15g \n",rms,vxa); fflush(stdout);
	rmsvx=rms;




	rms = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
		ijk = k + nz*(j + ny*i) ;
		rms += (vy[ijk]-vya)*(vy[ijk]-vya);
			
	}
	rms = sqrt(rms/(nx*ny*nz)) ;

	fprintf(stdout,"finalrms(vy): %21.15g avg: %21.15g \n",rms,vya); fflush(stdout);
	rms = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
		ijk = k + nz*(j + ny*i) ;
		rms += (vz[ijk]-vza)*(vz[ijk]-vza) ;
	}
	rms = sqrt(rms/(nx*ny*nz)) ;

	fprintf(stdout,"finalrms(vz): %21.15g avg: %21.15g \n",rms,vza); fflush(stdout);

	ke = 0. ;
	ke1=ke2=ke3=0.;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  ke1+=0.5*rho0*(vx[ijk]*vx[ijk])*dV;
	  ke2+=0.5*rho0*(vy[ijk]*vy[ijk])*dV;
	  ke3+=0.5*rho0*(vz[ijk]*vz[ijk])*dV;
	}
	ke=ke1+ke2+ke3;

	// output interesting things
	fprintf(stdout,"rms: %21.15g amp: %21.15g ke: %21.15g ke1: %21.15g ke2: %21.15g ke3: %21.15g\n",rmstot,amp,ke,ke1,ke2,ke3); fflush(stdout);


#endif




#if(WHICHVARIABLE==1)

	///////////////
	//
	// compute max and average of pressure
        pmax = 0. ;
        pavg = 0. ;
	petot = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  if(press[ijk]>pmax) pmax=press[ijk];
	  petot+=press[ijk]*dV;
	}
	pavg=petot/Vol; // average pressure

	fprintf(stderr,"pavg=%21.15g petot=%21.15g pmax=%21.15g\n",pavg,petot,pmax);

	//////////////////
	//
	// first shape the vector potential to vanish according to pressure of star

	///////////////////
	//
	// Need boundary zones -- so assume this includes 1 boundary zone
	//
	// Field only exists inside 1...n-1
	//
	//
	//
	//////////////////////////////

	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  // norm doesn't really matter since normalized later by beta, but divide by pmax for fun.
	  vx[ijk]*=press[ijk]/pmax;
	  vy[ijk]*=press[ijk]/pmax;
	  vz[ijk]*=press[ijk]/pmax;
	}

#if(TESTPOT==1) // test vector potentials
	//////////////////
	//

	// Bx = Az,y - Ay,z
	// By = Ax,z - Az,x
	// Bz = -Ax,y + Ay,x

	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;

	  // centered x,y,z
	  x=-Lx/2+(i+0.5)*dx;
	  y=-Ly/2+(j+0.5)*dy;
	  z=-Lz/2+(k+0.5)*dz;

	  // norm doesn't really matter since normalized later by beta, but divide by pmax for fun.
	  vx[ijk]=0;
	  vy[ijk]=0;
	  vz[ijk]=x;
	}
#endif

	////////////////////////
	//
	// compute alternative vector potential

	if(STENCILTYPE==0){
	  // already centered appropriately
	}
	else if(STENCILTYPE==1){
	  // move vector potential from center to corner

	  for(i=1;i<nx;i++) for(j=1;j<ny;j++) for(k=1;k<nz;k++) {
	    ijk = k + nz*(j + ny*i) ;
	    di=ny*nz;
	    dj=nz;
	    dk=1;
	    // centered A to corner A
	    vx2[ijk]=0.125*(  (vx[ijk-dk]+vx[ijk-dk -dj]+vx[ijk-dk -di]+vx[ijk-dk -dj-di]) + (vx[ijk]+vx[ijk -dj]+vx[ijk -di]+vx[ijk -dj-di]) );
	    vy2[ijk]=0.125*(  (vy[ijk-dk]+vy[ijk-dk -dj]+vy[ijk-dk -di]+vy[ijk-dk -dj-di]) + (vy[ijk]+vy[ijk -dj]+vy[ijk -di]+vy[ijk -dj-di]) );
	    vz2[ijk]=0.125*(  (vz[ijk-di]+vz[ijk-di -dj]+vz[ijk-di -dk]+vz[ijk-di -dj-dk]) + (vz[ijk]+vz[ijk -dj]+vz[ijk -dk]+vz[ijk -dj-dk]) );
	  }
	}
	else if(STENCILTYPE==2){ // ZEUS super-staggered placement of vector potential

	}




	//////////////////////////
	//
	// compute field (B)

	// initialize all zones to 0
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  Bx[ijk]=0.0;
	  By[ijk]=0.0;
	  Bz[ijk]=0.0;
	}
	// set rest using vector potential
	// B^i = \epsilon^{ijk} \partial_j A_k

	// only works if field really 0 on boundaries
	if(STENCILTYPE==0){
	  // for centered A, basic centered difference to get B at center so A.B easy, divB=0 not on any stencil?
	  for(i=1;i<nx-1;i++) for(j=1;j<ny-1;j++) for(k=1;k<nz-1;k++) {
	    ijk = k + nz*(j + ny*i) ;
	    Bx[ijk]=(vz[ijk+nz]-vz[ijk-nz])/(2.0*dy)-(vy[ijk+1]-vy[ijk-1])/(2.0*dz); // Bx = Az,y - Ay,z
	    By[ijk]=(vx[ijk+1]-vx[ijk-1])/(2.0*dz)-(vz[ijk+nz*ny]-vz[ijk-nz*ny])/(2.0*dx); // By = Ax,z - Az,x
	    Bz[ijk]=(-vx[ijk+nz]-vz[ijk-nz])/(2.0*dy)+(vy[ijk+nz*ny]-vy[ijk-nz*ny])/(2.0*dx); // Bz = -Ax,y + Ay,x
	  }
	}
	else if(STENCILTYPE==1){
	  // corners for A, B at center,  this keeps a definition of divB=0 (HARM stencil) divB is at corner

	  for(i=0;i<nx-1;i++) for(j=0;j<ny-1;j++) for(k=0;k<nz-1;k++) {
	    ijk = k + nz*(j + ny*i) ;
	    di=ny*nz;
	    dj=nz;
	    dk=1;

	    if(SUBSTENCILTYPE==0){
	      // Bx = Az,y - Ay,z
	      Bx[ijk]+=(  (vz2[ijk+dj]+vz2[ijk+dj +dk]+vz2[ijk+dj +di]+vz2[ijk+dj +dk+di]) - (vz2[ijk]+vz2[ijk +dk]+vz2[ijk +di]+vz2[ijk +dk+di]) )/(4.0*dy);
	      Bx[ijk]-=(  (vy2[ijk+dk]+vy2[ijk+dk +dj]+vy2[ijk+dk +di]+vy2[ijk+dk +dj+di]) - (vy2[ijk]+vy2[ijk +dj]+vy2[ijk +di]+vy2[ijk +dj+di]) )/(4.0*dz);
	      
	      // By = Ax,z - Az,x
	      By[ijk]+=(  (vx2[ijk+dk]+vx2[ijk+dk +dj]+vx2[ijk+dk +di]+vx2[ijk+dk +dj+di]) - (vx2[ijk]+vx2[ijk +dj]+vx2[ijk +di]+vx2[ijk +dj+di]) )/(4.0*dz);
	      By[ijk]-=(  (vz2[ijk+di]+vz2[ijk+di +dj]+vz2[ijk+di +dk]+vz2[ijk+di +dj+dk]) - (vz2[ijk]+vz2[ijk +dj]+vz2[ijk +dk]+vz2[ijk +dj+dk]) )/(4.0*dx);

	      // Bz = -Ax,y + Ay,x
	      Bz[ijk]-=(  (vx2[ijk+dj]+vx2[ijk+dj +dk]+vx2[ijk+dj +di]+vx2[ijk+dj +dk+di]) - (vx2[ijk]+vx2[ijk +dk]+vx2[ijk +di]+vx2[ijk +dk+di]) )/(4.0*dy);
	      Bz[ijk]+=(  (vy2[ijk+di]+vy2[ijk+di +dj]+vy2[ijk+di +dk]+vy2[ijk+di +dj+dk]) - (vy2[ijk]+vy2[ijk +dj]+vy2[ijk +dk]+vy2[ijk +dj+dk]) )/(4.0*dx);
	    }
	    else if(SUBSTENCILTYPE==1){
	      // use centered vx,vy,vz directly (then only requires 2 values instead of 4)
	      // same as STENCILTYPE==0
	      Bx[ijk]=(vz[ijk+nz]-vz[ijk-nz])/(2.0*dy)-(vy[ijk+1]-vy[ijk-1])/(2.0*dz); // Bx = Az,y - Ay,z
	      By[ijk]=(vx[ijk+1]-vx[ijk-1])/(2.0*dz)-(vz[ijk+nz*ny]-vz[ijk-nz*ny])/(2.0*dx); // By = Ax,z - Az,x
	      Bz[ijk]=(-vx[ijk+nz]-vz[ijk-nz])/(2.0*dy)+(vy[ijk+nz*ny]-vy[ijk-nz*ny])/(2.0*dx); // Bz = -Ax,y + Ay,x
	    }
	  }
	}
	else if(STENCILTYPE==2){
	  // ZEUS version, A in special corners, B at spatial edges, and divB at center
	  for(i=1;i<nx-1;i++) for(j=1;j<ny-1;j++) for(k=1;k<nz-1;k++) {
	    ijk = k + nz*(j + ny*i) ;
	    // not defined yet
	  }
	}

#if(0) // look at B
	for(i=1;i<nx-1;i++) for(j=1;j<ny-1;j++) for(k=1;k<nz-1;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  fprintf(stderr,"ijk=%d %d %d : B = %21.15g %21.15g %21.15g\n",i,j,k,Bx[ijk],By[ijk],Bz[ijk]);
	  //	  if((Bx[ijk]>0.0)||(Bz[ijk]>0.0)){ fprintf(stderr,"GOD\n"); fflush(stderr); }
	}
#endif

	////////////////////////////
	//
	// compute max and average of bsq
        bsqmax = 0. ;
	bsqavg = 0. ;
	betot=0.0;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  bsq=Bx[ijk]*Bx[ijk]+By[ijk]*By[ijk]+Bz[ijk]*Bz[ijk];	  
	  if(bsq>bsqmax) bsqmax=bsq;
	  betot+=bsq*dV;
	}
	bsqavg=betot/Vol; // average bsq
	betot*=0.5;  // total energy

	if(betanormtype==0){
	  betaact=pmax/(0.5*bsqmax);
	}
	else{
	  betaact=pavg/(0.5*bsqavg);
	}
	norm=sqrt(betaact/beta);
	fprintf(stderr,"betaact=%21.15g norm=%21.15g\n",betaact,norm);

	///////////////////////
	//
	// normalize field and vector potential
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  Bx[ijk]*=norm;
	  By[ijk]*=norm;
	  Bz[ijk]*=norm;
	  // normalize vector potential so consistent with field
	  vx[ijk]*=norm;
	  vy[ijk]*=norm;
	  vz[ijk]*=norm;
	  // normalize special vector potential so consistent with field
	  vx2[ijk]*=norm;
	  vy2[ijk]*=norm;
	  vz2[ijk]*=norm;
	}


	////////////////////////////
	//
	// compute divA
	if((STENCILTYPE==0)||(STENCILTYPE==1)){ // B is at center, so just place divB at corner.
	  // stenciltype=0 same as stenciltype=1 (which is like HARM)
	  divAavg=0.0;
	  divAmax=0.0;
	  divAavgnorm=0.0;
	  divAmaxnorm=0.0;
	  divAimax=-1;
	  divAjmax=-1;
	  divAkmax=-1;
	  divAimaxnorm=-1;
	  divAjmaxnorm=-1;
	  divAkmaxnorm=-1;
	  divAsum=0.0;
	  divAnormsum=0.0;


	  for(i=2;i<nx-1;i++) for(j=2;j<ny-1;j++) for(k=2;k<nz-1;k++) {
	    ijk = k + nz*(j + ny*i) ;
	    di=ny*nz;
	    dj=nz;
	    dk=1;

	    divA1=((vx[ijk]+vx[ijk-dj]+vx[ijk-dk]+vx[ijk-dj-dk])-(vx[ijk-di]+vx[ijk-di-dj]+vx[ijk-di-dk]+vx[ijk-di-dj-dk]))/(4.0*dx);
	    divA2=((vy[ijk]+vy[ijk-di]+vy[ijk-dk]+vy[ijk-di-dk])-(vy[ijk-dj]+vy[ijk-dj-di]+vy[ijk-dj-dk]+vy[ijk-dj-di-dk]))/(4.0*dy);
	    divA3=((vz[ijk]+vz[ijk-dj]+vz[ijk-di]+vz[ijk-dj-di])-(vz[ijk-dk]+vz[ijk-dk-dj]+vz[ijk-dk-di]+vz[ijk-dk-dj-di]))/(4.0*dz);
	    divA=fabs(divA1+divA2+divA3);
	    // fprintf(stderr,"%d %d %d divA: %21.15g : %21.15g %21.15g %21.15g : %21.15g %21.15g\n",i,j,k,divA,divA1,divA2,divA3,dy,vy[ijk]);

	    if(divA==0.0){ divAnorm=0.0; divAnormalization=0.0; }
	    else{
	      divAnormalization=fabs( (((vx[ijk]+vx[ijk-dj]+vx[ijk-dk]+vx[ijk-dj-dk])+(vx[ijk-di]+vx[ijk-di-dj]+vx[ijk-di-dk]+vx[ijk-di-dj-dk]))/(4.0*dx))+( ((vy[ijk]+vy[ijk-di]+vy[ijk-dk]+vy[ijk-di-dk])-(vy[ijk-dj]+vy[ijk-dj-di]+vy[ijk-dj-dk]+vy[ijk-dj-di-dk]))/(4.0*dy))+( ((vz[ijk]+vz[ijk-dj]+vz[ijk-di]+vz[ijk-dj-di])-(vz[ijk-dk]+vz[ijk-dk-dj]+vz[ijk-dk-di]+vz[ijk-dk-dj-di]))/(4.0*dz)));
	      if(divAnormalization==0.0) divAnorm=0.0;
	      else divAnorm=divA/divAnormalization;
	    }

	    divAsum+=divA;
	    divAnormsum+=divAnormalization;

	    divAavg+=divA;
	    divAavgnorm+=divAnorm;

	    if(divA>divAmax){
	      divAmax=divA;
	      divAimax=i;
	      divAjmax=j;
	      divAkmax=k;
	    }
	    if(divAnorm>divAmaxnorm){
	      divAmaxnorm=divAnorm;
	      //	      fprintf(stderr,"got here: %d %d %d : %21.15g %21.15g\n",i,j,k,divAmaxnorm,divAnorm);
	      divAimaxnorm=i;
	      divAjmaxnorm=j;
	      divAkmaxnorm=k;
	    }
	  }
	  divAavg/=(nx*ny*nz);
	  divAavgnorm/=(nx*ny*nz);
	  divAavgnorm2=divAsum/divAnormsum;
	}
	else if(STENCILTYPE==2){ // ZEUS form of divA
	  // NOT DONE YET
	}
	//	fprintf(stderr,"2 divAavg=%g divAavgnorm=%g divAmax=%g divAmaxnorm=%g\n",divAavg,divAavgnorm,divAmax,divAmaxnorm);


	////////////////////////////
	//
	// compute divB
	if((STENCILTYPE==0)||(STENCILTYPE==1)){ // B is at center, so just place divB at corner.
	  // stenciltype=0 same as stenciltype=1 (which is like HARM)
	  divBavg=0.0;
	  divBmax=0.0;
	  divBavgnorm=0.0;
	  divBmaxnorm=0.0;
	  divBimax=-1;
	  divBjmax=-1;
	  divBkmax=-1;
	  divBimaxnorm=-1;
	  divBjmaxnorm=-1;
	  divBkmaxnorm=-1;
	  divBsum=0.0;
	  divBnormsum=0.0;

	  for(i=2;i<nx-1;i++) for(j=2;j<ny-1;j++) for(k=2;k<nz-1;k++) {
	    ijk = k + nz*(j + ny*i) ;
	    di=ny*nz;
	    dj=nz;
	    dk=1;

	    divB1=((Bx[ijk]+Bx[ijk-dj]+Bx[ijk-dk]+Bx[ijk-dj-dk])-(Bx[ijk-di]+Bx[ijk-di-dj]+Bx[ijk-di-dk]+Bx[ijk-di-dj-dk]))/(4.0*dx);
	    divB2=((By[ijk]+By[ijk-di]+By[ijk-dk]+By[ijk-di-dk])-(By[ijk-dj]+By[ijk-dj-di]+By[ijk-dj-dk]+By[ijk-dj-di-dk]))/(4.0*dy);
	    divB3=((Bz[ijk]+Bz[ijk-dj]+Bz[ijk-di]+Bz[ijk-dj-di])-(Bz[ijk-dk]+Bz[ijk-dk-dj]+Bz[ijk-dk-di]+Bz[ijk-dk-dj-di]))/(4.0*dz);
	    divB=fabs(divB1+divB2+divB3);
	    // fprintf(stderr,"%d %d %d divB: %21.15g : %21.15g %21.15g %21.15g : %21.15g %21.15g\n",i,j,k,divB,divB1,divB2,divB3,dy,By[ijk]);

	    if(divB==0.0){ divBnorm=0.0; divBnormalization=0.0; }
	    else{
	      divBnormalization=fabs( (((Bx[ijk]+Bx[ijk-dj]+Bx[ijk-dk]+Bx[ijk-dj-dk])+(Bx[ijk-di]+Bx[ijk-di-dj]+Bx[ijk-di-dk]+Bx[ijk-di-dj-dk]))/(4.0*dx))+( ((By[ijk]+By[ijk-di]+By[ijk-dk]+By[ijk-di-dk])-(By[ijk-dj]+By[ijk-dj-di]+By[ijk-dj-dk]+By[ijk-dj-di-dk]))/(4.0*dy))+( ((Bz[ijk]+Bz[ijk-dj]+Bz[ijk-di]+Bz[ijk-dj-di])-(Bz[ijk-dk]+Bz[ijk-dk-dj]+Bz[ijk-dk-di]+Bz[ijk-dk-dj-di]))/(4.0*dz)));
	      if(divBnormalization==0.0) divBnorm=0.0;
	      else divBnorm=divB/divBnormalization;
	    }

	    divBsum+=divB;
	    divBnormsum+=divBnormalization;

	    divBavg+=divB;
	    divBavgnorm+=divBnorm;

	    if(divB>divBmax){
	      divBmax=divB;
	      divBimax=i;
	      divBjmax=j;
	      divBkmax=k;
	    }
	    if(divBnorm>divBmaxnorm){
	      divBmaxnorm=divBnorm;
	      //	      fprintf(stderr,"got here: %d %d %d : %21.15g %21.15g\n",i,j,k,divBmaxnorm,divBnorm);
	      divBimaxnorm=i;
	      divBjmaxnorm=j;
	      divBkmaxnorm=k;
	    }
	  }
	  divBavg/=(nx*ny*nz);
	  divBavgnorm/=(nx*ny*nz);
	  divBavgnorm2=divBsum/divBnormsum;
	}
	else if(STENCILTYPE==2){ // ZEUS form of divB
	  // NOT DONE YET
	}
	//	fprintf(stderr,"2 divBavg=%g divBavgnorm=%g divBmax=%g divBmaxnorm=%g\n",divBavg,divBavgnorm,divBmax,divBmaxnorm);


	////////////////////
	//
	// compute max and average of bsq AGAIN
        bsqmax = 0. ;
	bsqavg = 0. ;
	betot=0.0;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  bsq=Bx[ijk]*Bx[ijk]+By[ijk]*By[ijk]+Bz[ijk]*Bz[ijk];	  
	  if(bsq>bsqmax) bsqmax=bsq;
	  betot+=bsq*dV;
	}
	bsqavg=betot/Vol; // average bsq
	betot*=0.5;  // total energy

	if(betanormtype==0){
	  betaact=pmax/(0.5*bsqmax);
	}
	else{
	  betaact=pavg/(0.5*bsqavg);
	}
	norm=sqrt(betaact/beta);
	fprintf(stderr,"new betaact=%21.15g new norm=%21.15g\n",betaact,norm);


	/////////////////////////////////////////
	//
	// compute helicity
	helicity = 0. ;

	if(STENCILTYPE==0){
	  for(i=0;i<nx;i++) for(j=0;j<ny;j++) for(k=0;k<nz;k++) {
	    ijk = k + nz*(j + ny*i) ;
	    helicity+=Bx[ijk]*vx[ijk]+By[ijk]*vy[ijk]+Bz[ijk]*vz[ijk];	  
	  }
	}
	else if(STENCILTYPE==1){ // HARM form of divB
	  // corners for A, B at center,  this keeps a definition of divB=0 (HARM stencil) divB is at corner.

	  // B only exists for 1 ... n-1
	  // Aavg only exists for 0 ... n-2

	  // so calculation is for 1 ... n-2

	  // Thus need to average A to location of B to take A.B
	  for(i=1;i<nx-1;i++) for(j=1;j<ny-1;j++) for(k=1;k<nz-1;k++) {
	    ijk = k + nz*(j + ny*i) ;
	    di=ny*nz;
	    dj=nz;
	    dk=1;
	    
	    if(SUBSTENCILTYPE==0){
	      // corner A averaged to location of centered B
	      Axavg=0.125*(  (vx2[ijk+dk]+vx2[ijk+dk +dj]+vx2[ijk+dk +di]+vx2[ijk+dk +dj+di]) + (vx2[ijk]+vx2[ijk +dj]+vx2[ijk +di]+vx2[ijk +dj+di]) );
	      Ayavg=0.125*(  (vy2[ijk+dk]+vy2[ijk+dk +dj]+vy2[ijk+dk +di]+vy2[ijk+dk +dj+di]) + (vy2[ijk]+vy2[ijk +dj]+vy2[ijk +di]+vy2[ijk +dj+di]) );
	      Azavg=0.125*(  (vz2[ijk+di]+vz2[ijk+di +dj]+vz2[ijk+di +dk]+vz2[ijk+di +dj+dk]) + (vz2[ijk]+vz2[ijk +dj]+vz2[ijk +dk]+vz2[ijk +dj+dk]) );
	    }
	    else if(SUBSTENCILTYPE==1){ // directly use originally centered vector potential
	      Axavg=vx[ijk];
	      Ayavg=vy[ijk];
	      Azavg=vz[ijk];
	    }

	    helicity+=(Bx[ijk]*Axavg+By[ijk]*Ayavg+Bz[ijk]*Azavg)*dV;

	  }
	}
	else if(STENCILTYPE==2){ // ZEUS form of divB
	  // NOT DONE YET
	}

	alphahel=stellarradius*betot/helicity; // R\alpha



	//////////////////////
	//
	// output interesting things
	fprintf(stdout,"seed %d betaact %21.15g helicity %21.15g betot %21.15g alphahel=%21.15g\n",seed,betaact,helicity,betot,alphahel);


#endif

	// generally interesting things
	fprintf(stderr,"seed=%d\n",seed); fflush(stdout);
	fprintf(stderr,"nx=%d ny=%d nz=%d\n",nx,ny,nz);
	fprintf(stderr,"amp=%g nodiv=%d kmin=%g kmax=%g bx0=%g EQK=%g\n",amp,nodiv,kmin,kmax,bx0,EQK);
	fprintf(stderr,"divBavg=%g divBavgnorm=%g divBmax=%g divBmaxnorm=%g divBavgnorm2=%g\n",divBavg,divBavgnorm,divBmax,divBmaxnorm,divBavgnorm2);
	fprintf(stderr,"divBimax=%d divBjmax=%d divBkmax=%d  : divBimaxnorm=%d divBjmaxnorm=%d divBkmaxnorm=%d\n",divBimax,divBjmax,divBkmax,divBimaxnorm,divBjmaxnorm,divBkmaxnorm);
	fprintf(stderr,"divAavg=%g divAavgnorm=%g divAmax=%g divAmaxnorm=%g divAavgnorm2=%g\n",divAavg,divAavgnorm,divAmax,divAmaxnorm,divAavgnorm2);
	fprintf(stderr,"divAimax=%d divAjmax=%d divAkmax=%d  : divAimaxnorm=%d divAjmaxnorm=%d divAkmaxnorm=%d\n",divAimax,divAjmax,divAkmax,divAimaxnorm,divAjmaxnorm,divAkmaxnorm);
	fflush(stderr);

	// output to file some things
	
	// par1 used to track INDEX number
	fprintf(outputfile,"%d %21.15g %d %21.15g %21.15g %21.15g %d %d ",par1,INDEX,seed,EQK,kmin,kmax,nodiv,betanormtype);
	fprintf(outputfile,"%21.15g %21.15g %21.15g %21.15g ",betaact,helicity,betot,alphahel);
	fprintf(outputfile,"%21.15g %21.15g %21.15g %21.15g %21.15g ",divBavg,divBavgnorm,divBmax,divBmaxnorm,divBavgnorm2);
	fprintf(outputfile,"%21.15g %21.15g %21.15g %21.15g %21.15g ",divAavg,divAavgnorm,divAmax,divAmaxnorm,divAavgnorm2);
	fprintf(outputfile,"\n");
	fflush(outputfile);



	}}}}}}} // assume no files outputted when going over realizations





	//////////////////////////////////////////////////
	//
	// output files (if any)
	//
	//////////////////////////////////////////////////



#if(OUTPUTTYPE==0)

	// open file(s)
	for(myid=0;myid<numprocs;myid++){
	  if(numprocs>1){
	    sprintf(filename,"%s%s.%04d",basefilename,".in",myid);
	  }
	  else{
	    sprintf(filename,"%s%s",basefilename,".in");
	  }
	  if((files[myid]=fopen(filename,"w"))==NULL){
	    fprintf(stdout,"Can't open %s\n",filename);
	  }
	}



	/* output scaled velocities */
	for(myid=0;myid<numprocs;myid++){ // assumes each CPU same size
	  mycpupos[1]=myid%ncpux1;
	  mycpupos[2]=(int)((myid%(ncpux1*ncpux2))/ncpux1);
	  mycpupos[3]=(int)(myid/(ncpux1*ncpux2));
	  for(i=mycpupos[1]*n1;i<(mycpupos[1]+1)*n1;i++) for(j=mycpupos[2]*n2;j<(mycpupos[2]+1)*n2;j++) for(k=mycpupos[3]*n3;k<(mycpupos[3]+1)*n3;k++) {
	    if(READIN){
	      // gammie order
	      ijk = k + nz*(j + ny*i) ;
	      dumi=i-mycpupos[1]*n1; // relative offset
	      fwrite(&dumi,sizeof(int),1,files[myid]);
	      dumi=j-mycpupos[2]*n2; // relative offset
	      fwrite(&dumi,sizeof(int),1,files[myid]);
	      dumi=k-mycpupos[3]*n3; // relative offset
	      fwrite(&dumi,sizeof(int),1,files[myid]);
	      dumlf=vx[ijk];
	      fwrite(&dumlf,sizeof(double),1,files[myid]);
	      dumlf=vy[ijk];
	      fwrite(&dumlf,sizeof(double),1,files[myid]);
	      dumlf=vz[ijk];
	      fwrite(&dumlf,sizeof(double),1,files[myid]);
	    }
	    else {
	      // c-order
	      ijk = i + nx*(j + ny*k) ;
	      if(REALTYPE==1){
		dumlf=vx[ijk];
		fwrite(&dumlf,sizeof(double),1,files[myid]);
		dumlf=vy[ijk];
		fwrite(&dumlf,sizeof(double),1,files[myid]);
		dumlf=vz[ijk];
		fwrite(&dumlf,sizeof(double),1,files[myid]);
		dumlf=bx0;
		fwrite(&dumlf,sizeof(double),1,files[myid]);
		dumlf=by0;
		fwrite(&dumlf,sizeof(double),1,files[myid]);
		dumlf=bz0;
		fwrite(&dumlf,sizeof(double),1,files[myid]);
		dumlf=rho0;
		fwrite(&dumlf,sizeof(double),1,files[myid]);
	      }
	      else{
		duml=vx[ijk];
		fwrite(&duml,sizeof(float),1,files[myid]);
		duml=vy[ijk];
		fwrite(&duml,sizeof(float),1,files[myid]);
		duml=vz[ijk];
		fwrite(&duml,sizeof(float),1,files[myid]);
		duml=bx0;
		fwrite(&duml,sizeof(float),1,files[myid]);
		duml=by0;
		fwrite(&duml,sizeof(float),1,files[myid]);
		duml=bz0;
		fwrite(&duml,sizeof(float),1,files[myid]);
		duml=rho0;
		fwrite(&duml,sizeof(float),1,files[myid]);
	      }
	      
	    }
	    /*
	    fprintf(files[myid],"%d %d %d %21.15g %21.15g %21.15g\n",
		    i,j,k,
		    vx[ijk]*amp/rms,
		    vy[ijk]*amp/rms,
		    vz[ijk]*amp/rms) ;
	    */
	  }
	  fclose(files[myid]);
	}

#endif


	// should really clear memory

	/* done! */
}

// wavenumber in Frenchy form (i.e. really an integer and not 2\pi times an integer)

/* these functions returns the wavenumber corresponding to a
 * point on the rfftw grid */
double return_kx(int i, int j, int k)
{
	return((double)((i + nx/2)%nx - nx/2)) ;
}
double return_ky(int i, int j, int k)
{
	return((double)((j + ny/2)%ny - ny/2)) ;
}
double return_kz(int i, int j, int k)
{
	return((double)((k + nz/2)%nz - nz/2)) ;
}
double return_K(int i, int j, int k)
{
	double kx,ky,kz,return_kx(int i, int j, int k),
		return_ky(int i, int j, int k),return_kz(int i, int j, int k) ;
	kx = return_kx(i,j,k) ;
	ky = return_ky(i,j,k) ;
	kz = return_kz(i,j,k) ;

	return(sqrt(kx*kx + ky*ky + kz*kz)) ;
}

double rand_ph(double k)
{
  int cond;
  int condition_k(double k);
  double ranc(int seed);
  
  cond=condition_k(k);

#if(0)
  if(cond==1){ fprintf(stdout,"ph: %ld %g\n",rannum++,k);}
#endif

  if(cond==1) return(2.*M_PI*MYRAND) ;	
  else if(cond==2) return(2.*M_PI*0.5);	
  else return(0.);
}

/* this returns a gaussian-distributed variable
 * with expected amplitude = power spectrum,
 * following NR's gaussian generator */
double rand_P(double k)
{
  double ranc(int seed);
  int condition_k(double k);
  double P(double k) ;
  double p,x1,x2;
  double x10,x20;
  int cond;

	x10=0.5;
	x20=0.5;

	cond=condition_k(k);

#if(0)
	if(cond==1){ fprintf(stdout,"P: %ld %g\n",rannum++,k);}
	if(cond==1){ fprintf(stdout,"P: %ld %g\n",rannum++,k);}
#endif
	
	if(cond==1){
	  p = P(k) ;
	  x1 = MYRAND ;
	  x2 = MYRAND ;
	  p *= sqrt(-2.*log(x1))*cos(2.*M_PI*x2) ;
	}
	else if(cond==2){
	  p = P(k) ;
	  x1 = x10;
	  x2 = x20;
	  p *= sqrt(-2.*log(x1))*cos(2.*M_PI*x2) ;
	}
	else {
	  p = 0.0;
	}

	return(p) ;
}


double P(double k)
{
  return(pow(k,INDEX)) ;
}

// 0: no power/phase
// 1: power phase randomly
// 2: power phase non-randomly
int condition_k(double k)
{
  if(k>EQK){
    if( k >= kmin && k <= kmax && k > 0)
      return(1) ;
    else
      return(0) ;
  }
  else{
    if( k >= kmin && k <= kmax && k > 0)
      return(2) ;
    else
      return(0) ;
  }
}

