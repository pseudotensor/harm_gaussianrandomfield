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
 */




//#include "global.h"



#define OUTPUTTYPE 1
// 0: output files as below
// 1: no data files outputted, just output of statistics for each realization

#define NUMREALIZATIONS 3 // number of realizations to loop over

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

// ignore it = 0.0
// 2.0 : used for good Watson run
#define EQK (2.0)

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

int main(int argc, char *argv[])
{
	int ijk,i,j,k ;
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
	double helicity;

	FILE*averydata;
	double x,y,z,radius;




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

	
	files=(FILE**)malloc(sizeof(FILE*)*ncpux1*ncpux2*ncpux3);
	if(files==NULL){
	  fprintf(stderr,"cannot allocate file name list\n");
	  exit(1);
	}


	strcpy(basefilename,argv[argnum++]);

#if(WHICHVARIABLE==0)
	rho0=1.0; // just for ke
	bx0=3.0; // Watson bx=10,3 (new best)
	by0=bz0=0.0;
	Lx=1.0; // just for ke
	Ly=1.0;  // just for ke
	Lz=1.0;  // just for ke


	// Watson
	// sub-alfvenic
	// subsonic (at best mach=1)
	// many waves within box, not just 1
	// rms velocity along field larger than rms velocity perp field (want other way)

	amp = 4.0; // Watson amp=4
	nodiv = 1 ;	
	kmin = 2.0 ; //  Watson kmin=2	
	kmax = 32.0 ; // Watson kmax=31
#endif
#if(WHICHVARIABLE==1)
	Lx=1.0; // just for ke
	Ly=1.0;  // just for ke
	Lz=1.0;  // just for ke


	beta = 100.0; 
	nodiv = 0 ;	
	kmin = 2.0 ; //  Watson kmin=2	
	kmax = 32.0 ; // Watson kmax=31

#endif
  

	fprintf(stdout,"nx=%d ny=%d nz=%d\n",nx,ny,nz);
	fprintf(stdout,"amp=%g nodiv=%d kmin=%g kmax=%g bx0=%g EQK=%g\n",amp,nodiv,kmin,kmax,bx0,EQK);
	fflush(stdout);

	if((nx/2<=kmax)||(ny/2<=kmax)||(nz/2<=kmax)){
	  fprintf(stderr,"kmax chosen poorly, need N>%d=2*kmax\n",(int)(2*kmax));
	  exit(1);
	}

	/* make some space for the transforms */
	kvx = (fftw_complex *) 	calloc(nx*ny*(nz/2+2) , sizeof(fftw_complex)) ;
	kvy = (fftw_complex *) 	calloc(nx*ny*(nz/2+2) , sizeof(fftw_complex)) ;
	kvz = (fftw_complex *) 	calloc(nx*ny*(nz/2+2) , sizeof(fftw_complex)) ;

	vx = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	vy = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	vz = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;

	Bx = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	By = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	Bz = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;

	rho = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	mofr = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
        press = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;

	/* make an fftw plan; don't measure performance */
	pinv = rfftw3d_create_plan(nx,ny,nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE) ;


#if(WHICHVARIABLE==1)
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

	// open and read avery's data file
	averydata=fopen("A0III_1.4","rt");
	if(averydata==NULL){
	  fprintf(stderr,"couldn't open avery's file\n");
	  exit(1);
	}
#define LINESTOSKIP 9

	for(i=1;i<=LINESTOSKIP;i++) while(fgetc(averydata)!='\n'); // skip lines
	for(i=0;i<NAVERY;i++){
	  fscanf(averydata,"%lf %lf %lf %lf\n",&rvsr[i],&rhovsr[i],&mofrvsr[i],&pressvsr[i]);
	  while(fgetc(averydata)!='\n'); // skip rest of line
	}
	fclose(averydata);


	// interpolate avery data and assign to 3D grid
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;

	  x=-Lx/2+i*Lx/nx;
	  y=-Ly/2+j*Ly/ny;
	  z=-Lz/2+k*Lz/nz;

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
	    // linear interpolation
	    rho[ijk]=rhovsr[l]+(rhovsr[l+1]-rhovsr[l])/(rvsr[l+1]-rvsr[l])*(radius-rvsr[l]);
	    mofr[ijk]=mofrvsr[l]+(mofrvsr[l+1]-mofrvsr[l])/(rvsr[l+1]-rvsr[l])*(radius-rvsr[l]);
	    press[ijk]=pressvsr[l]+(pressvsr[l+1]-pressvsr[l])/(rvsr[l+1]-rvsr[l])*(radius-rvsr[l]);
	  }
	}

#endif

	
	//	myrms=2.31029310017193;
	//rmsvx=0;
	//seed=-1;
	//	while(fabs(myrms-rmsvx)>1E-6){
	if(1)
	{
	  //seed++;
	  seed=117;
	  //	  seed=1;
	//ranc(1); // new Watson alternative seed
	//	ranc(121); // Original Watson seed
	//ranc(119); // Original Watson seed
	//	ranc(119);
	//ranc(7);
	//srand(7);
	  ranc(seed);
	  fprintf(stdout,"seed=%d\n",seed); fflush(stdout);

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


	


#if(WHICHVARIABLE==0)
	/* scale */
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
	}
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
	  ke1+=0.5*rho0*(vx[ijk]*vx[ijk]);
	  ke2+=0.5*rho0*(vy[ijk]*vy[ijk]);
	  ke3+=0.5*rho0*(vz[ijk]*vz[ijk]);
	}
	ke1*=(Lx/((double)nx))*(Ly/((double)ny))*(Lz/((double)nz));
	ke2*=(Lx/((double)nx))*(Ly/((double)ny))*(Lz/((double)nz));
	ke3*=(Lx/((double)nx))*(Ly/((double)ny))*(Lz/((double)nz));
	ke=ke1+ke2+ke3;

#endif


#if(WHICHVARIABLE==1)



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
	for(i=1;i<nx-1;i++)
	for(j=1;j<ny-1;j++)
	for(k=1;k<nz-1;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  Bx[ijk]=0.5*((vz[ijk+nz]-vz[ijk-nz])-(vy[ijk+1]-vy[ijk-1]));
	  By[ijk]=0.5*((vx[ijk+1]-vx[ijk-1])-(vz[ijk+nz*ny]-vz[ijk-nz*ny]));
	  Bz[ijk]=0.5*((vx[ijk+nz]-vz[ijk-nz])-(vy[ijk+nz*ny]-vy[ijk-nz*ny]));
	}

	// compute max bsq
        bsqmax = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  bsq=Bx[ijk]*Bx[ijk]+By[ijk]*By[ijk]+Bz[ijk]*Bz[ijk];	  
	  if(bsq>bsqmax) bsqmax=bsq;
	}
	// compute max p_g
        pmax = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  if(press[ijk]>pmax) pmax=press[ijk];
	}
	betaact=pmax/(0.5*bsqmax);
	norm=sqrt(betaact/beta);
	fprintf(stderr,"betaact=%21.15g norm=%21.15g\n",betaact,norm);

	// normalize
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  Bx[ijk]*=norm;
	  By[ijk]*=norm;
	  Bz[ijk]*=norm;
	}

	// compute max bsq
        bsqmax = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  bsq=Bx[ijk]*Bx[ijk]+By[ijk]*By[ijk]+Bz[ijk]*Bz[ijk];	  
	  if(bsq>bsqmax) bsqmax=bsq;
	}
	betaact=pmax/(0.5*bsqmax);
	norm=sqrt(betaact/beta);
	fprintf(stderr,"new betaact=%21.15g new norm=%21.15g\n",betaact,norm);

	// compute helicity
	helicity = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  helicity+=Bx[ijk]*vx[ijk]+By[ijk]*vy[ijk]+Bz[ijk]*vz[ijk];	  
	}
	helicity*=(Lx/((double)nx))*(Ly/((double)ny))*(Lz/((double)nz));
	fprintf(stdout,"seed %d betaact %21.15g helicity %21.15g\n",seed,betaact,helicity);


#endif


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

#if(WHICHVARIABLE==0)

	fprintf(stdout,"rms: %21.15g amp: %21.15g ke: %21.15g ke1: %21.15g ke2: %21.15g ke3: %21.15g\n",rmstot,amp,ke,ke1,ke2,ke3); fflush(stdout);
#endif


	/* done! */
}

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

  if(cond==1){ fprintf(stdout,"ph: %ld %g\n",rannum++,k);}

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

	if(cond==1){ fprintf(stdout,"P: %ld %g\n",rannum++,k);}
	if(cond==1){ fprintf(stdout,"P: %ld %g\n",rannum++,k);}
	
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


/* this returns the power spectrum */
#define INDEX	(11./3.)
double P(double k)
{
  return(pow(k,-INDEX)) ;
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

