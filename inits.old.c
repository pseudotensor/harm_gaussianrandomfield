

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

// per CPU size, NOT total size
#define N1 (32)
#define N2 (32)
#define N3 (32)

#define NCPUX1 (1)
#define NCPUX2 (1)
#define NCPUX3 (1)

#include <stdio.h>
#include <math.h>
#include <rfftw.h>

int nx,ny,nz ;
double amp ;
int nodiv ;
double kmin,kmax ;

int chandran_init(int argc, char *argv[])
{
	int ijk,i,j,k ;
	int l;
	int myid;
	int ncpux1,ncpux2,ncpux3,numprocs,mycpupos[3+1];
	double ph,a,a_re,a_im,ranc(int seed),K,return_K(int i, int j, int k),
		rand_P(double K),return_kx(int i, int j, int k),
		return_ky(int i, int j, int k),return_kz(int i, int j, int k),
		kx,ky,kz,kdotv_re,kdotv_im,rms ;
	static fftw_complex *kvx,*kvy,*kvz ;
	static fftw_real *vx,*vy,*vz ;
	static rfftwnd_plan p, pinv ;
	rfftwnd_plan rfftw3d_create_plan(int nx, int ny, int nz, 
		fftw_direction dir, int flags) ;

	FILE * files[NCPUX1*NCPUX2*NCPUX3];
	char filename[200];
	double dumlf;
	int dumi;
	double ke,rho,Lx,Ly,Lz;

	/* set parameters */
	ranc(7);

	ncpux1=NCPUX1;
	ncpux2=NCPUX2;
	ncpux3=NCPUX3;
	numprocs=ncpux1*ncpux2*ncpux3;

	nx = N1*ncpux1 ;
	ny = N2*ncpux2 ;
	nz = N3*ncpux3 ;

	rho=1.0; // just for ke
	Lx=1.0; // just for ke
	Ly=1.0;  // just for ke
	Lz=1.0;  // just for ke

	amp = 0.2 ;
	//amp=0.0;
	//amp = 1.0;
	nodiv = 0 ;
	kmin = 0 ;
	kmax = 1 ;

	/* make some space for the transforms */
	kvx = (fftw_complex *) 	calloc(nx*ny*(nz/2+2) , sizeof(fftw_complex)) ;
	kvy = (fftw_complex *) 	calloc(nx*ny*(nz/2+2) , sizeof(fftw_complex)) ;
	kvz = (fftw_complex *) 	calloc(nx*ny*(nz/2+2) , sizeof(fftw_complex)) ;

	vx = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	vy = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;
	vz = (fftw_real *) calloc(nx*ny*(2*(nz/2+1)) , sizeof(fftw_real)) ;

	/* make an fftw plan; don't measure performance */
	pinv = rfftw3d_create_plan(nx,ny,nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE) ;

	/* fill in k-space arrays */
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz/2+1;k++) {
		K = return_K(i,j,k) ;

		ph = 2.*M_PI*ranc(0) ;	/* random phase */
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

		ph = 2.*M_PI*ranc(0) ;	/* random phase */
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

		ph = 2.*M_PI*ranc(0) ;	/* random phase */
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


	ke = 0. ;
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++) {
	  ijk = k + nz*(j + ny*i) ;
	  ke+=0.5*rho*(vx[ijk]*vx[ijk]+vy[ijk]*vy[ijk]+vz[ijk]*vz[ijk]);
	}
	ke*=(Lx/((double)nx));
	ke*=(Ly/((double)ny));
	ke*=(Lz/((double)nz));


	// open file(s)
	for(myid=0;myid<numprocs;myid++){
	  if(numprocs>1){
	    sprintf(filename,"%s%s.%04d","velocities",".in",myid);
	  }
	  else{
	    sprintf(filename,"%s%s","velocities",".in");
	  }
	  if((files[myid]=fopen(filename,"w"))==NULL){
	    fprintf(stderr,"Can't open %s\n",filename);
	  }
	}



	/* output scaled velocities */
	for(myid=0;myid<numprocs;myid++){ // assumes each CPU same size
	  mycpupos[1]=myid%ncpux1;
	  mycpupos[2]=(int)((myid%(ncpux1*ncpux2))/ncpux1);
	  mycpupos[3]=(int)(myid/(ncpux1*ncpux2));
	  for(i=mycpupos[1]*N1;i<(mycpupos[1]+1)*N1;i++) for(j=mycpupos[2]*N2;j<(mycpupos[2]+1)*N2;j++) for(k=mycpupos[3]*N3;k<(mycpupos[3]+1)*N3;k++) {
	    ijk = k + nz*(j + ny*i) ;

	    dumi=i-mycpupos[1]*N1; // relative offset
	    fwrite(&dumi,sizeof(int),1,files[myid]);
	    dumi=j-mycpupos[2]*N2; // relative offset
	    fwrite(&dumi,sizeof(int),1,files[myid]);
	    dumi=k-mycpupos[3]*N3; // relative offset
	    fwrite(&dumi,sizeof(int),1,files[myid]);
	    dumlf=vx[ijk];
	    fwrite(&dumlf,sizeof(double),1,files[myid]);
	    dumlf=vy[ijk];
	    fwrite(&dumlf,sizeof(double),1,files[myid]);
	    dumlf=vz[ijk];
	    fwrite(&dumlf,sizeof(double),1,files[myid]);
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
	fprintf(stderr,"rms: %21.15g amp: %21.15g ke: %21.15g\n",rms,amp,ke); fflush(stderr);


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

/* this returns a gaussian-distributed variable
 * with expected amplitude = power spectrum,
 * following NR's gaussian generator */
double rand_P(double k)
{
	double P(double k) ;
	double p,x1,x2,ranc(int seed) ;

	p = P(k) ;
	x1 = ranc(0) ;
	x2 = ranc(0) ;
	p *= sqrt(-2.*log(x1))*cos(2.*M_PI*x2) ;

	return(p) ;
}

/* this returns the power spectrum */
#define INDEX	(11./3.)
double P(double k)
{
	if(k >= kmin && k <= kmax && k > 0)
		return(pow(k,-INDEX)) ;
	else
		return(0.) ;
}
