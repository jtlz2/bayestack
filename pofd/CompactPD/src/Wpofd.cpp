//============================================================================
// Name        : Cpofd.cpp
// Author      : Song Chen
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Wpofd.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <fftw3.h>
#include <math.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>



#define PI 3.1415926
#define float_error 1e-6

#define amp_I 1.0
#define amp_II 1.0

int findbin(double flux, double Delta, int shift){
	return int(flux/Delta)+shift;
}

double interp( double xi, int Ndata, double* data ){

	double x[2][Ndata];
	memcpy(x, data, sizeof(double)*2*Ndata);



	if(xi>=x[0][0]&&xi<=x[0][Ndata-1]){

	gsl_interp_accel *acc
	      = gsl_interp_accel_alloc ();
	    gsl_spline *spline
	      = gsl_spline_alloc (gsl_interp_cspline, Ndata);

	    gsl_spline_init (spline, x[0], x[1], Ndata);


	       double yi = gsl_spline_eval (spline, xi, acc);
	        //printf ("%g %g\n", xi, yi);

	    gsl_spline_free (spline);
	    gsl_interp_accel_free (acc);
	    return yi;
	}else{
		printf("error: The x=%e should not beyond data[%e~%e] ! \n",xi,x[0][0],x[0][Ndata-1]);
		exit(0);
	}

}

double interp_adv( double xi, int Ndata, double* data, int kmin, int kmax ){

	double x[2][Ndata];
	memcpy(x, data, sizeof(double)*2*Ndata);





	if(xi>=x[0][kmin]&&xi<=x[0][kmax]){

	gsl_interp_accel *acc
	      = gsl_interp_accel_alloc ();
	    gsl_spline *spline
	      = gsl_spline_alloc (gsl_interp_cspline, kmax-kmin+1);

	    gsl_spline_init (spline, &x[0][kmin], &x[1][kmin], kmax-kmin+1);


	       double yi = gsl_spline_eval (spline, xi, acc);
	        //printf ("%g %g\n", xi, yi);

	    gsl_spline_free (spline);
	    gsl_interp_accel_free (acc);
	    return yi;
	}else{
		printf("error: The x=%e should not beyond data[%e~%e] ! \n",xi,x[0][kmin],x[0][kmax]);
		exit(0);
	}

}


double np_roll(double * in, double *out, int Nmax, int shift){

for(int i=0; i<Nmax; i++){
	if(i<Nmax-1-shift){
		out[i]=in[i+shift];
	}else{
		out[i]=in[i+shift-(Nmax-1)];
	}

}

	return 0.0;
}


double input_fw(int Mod,int Nmax, double dmin, double Delta, double * in)
{
	if (Mod==FFTW_FORWARD){
	   for(int j=0;j<Nmax;j++){

	            	   if(j==0){
	            		   in[j]=0.0+dmin;
	            	   }else {
	            		   in[j]=dmin+j*Delta;
	            	   }


	  }
	}else if(Mod==FFTW_BACKWARD){
		for(int j=0;j<Nmax;j++){


				    	   if(j==0){
				    		   in[j]=0.0+dmin;

				    	   }else if(j<0.5*Nmax+1) {
				    		   in[j]=dmin+j*1.0/(Nmax*Delta);

				    	   }else{
				    		   in[j]=-(Nmax-j)*1.0/(Nmax*Delta)-dmin;
				    	   }


				    }

	}else{
		printf("error: wrong input !\n");
		return -1.0;

	}


	return 0.0;
}




double chi2_res( int Nbins, double * xbin_min, double * xbin_max, double *data, double Delta, int N, double * out, double sigma_noise){



	double * flux=new double[N];

	double * outshift;

	outshift =(double *) malloc(sizeof(double)*N);

	double sigma=10*sigma_noise;

	if(input_fw(FFTW_FORWARD,N,-sigma, Delta,flux)==-1.0) exit(0);

	int shift=N-1-int((sigma/Delta));

	//printf("shift %d\n",shift);
	np_roll(out,outshift,N,shift);
	//printf("save mode: FORWARD !!\n");

int xshift=int((sigma/Delta));
	int Ncut= findbin(xbin_max[Nbins-1]*1e-6,Delta,xshift);
printf("Ncut %d\n",Ncut);

	double marray[2][Ncut+2];
	//memcpy(array[0],flux, sizeof(double)*N);
	//memcpy(array[1], out, sizeof(double)*N);


	for(int i=0;i<Ncut+2;i++){
		marray[0][i]=flux[i];
		marray[1][i]=outshift[i]/(N*Delta)/amp_II;
	}


    int bin_split=10;
    double bin_integration_element=0.0;

    double bin_delta;
    double chi2=0.0;
    int count=0;
    int totalpixel=0;

    for(int i=0;i<Nbins;i++){
    	totalpixel+=data[i];
    }
    printf("total pixels %d \n",totalpixel);

    int kmin,kmax;

    for(int i=0;i<Nbins;i++){
    	bin_integration_element=0.0;
    	bin_delta=(xbin_max[i]-xbin_min[i])/(bin_split-1);

    	kmin=findbin(xbin_min[i]*1e-6,Delta,xshift);
    	kmax=findbin(xbin_max[i]*1e-6,Delta,xshift);

    	if((kmin-10)>0){
    		kmin=kmin-10;
    	}else{
    		kmin=0;
    	}

    	if((kmax+10)<Ncut+1){
    		kmax=kmax+10;
    	}else{
    		kmax=Ncut+1;
    	}

    	for(int j=0;j<bin_split;j++){
    		double inputflux=(xbin_min[i]+j*bin_delta)*1e-6;
    		bin_integration_element+=interp_adv( inputflux,Ncut+2,(double *)marray, kmin, kmax)*bin_delta*1e-6;

    	}
    	// covariance mertix ??
    	/*if(i==30){
    		printf("data[i] %lf\n",data[i]);
    		printf("theory %lf\n",bin_integration_element*totalpixel);
    		exit(0);
    	}*/

    	if(data[i]>0.0){
    		count++;
    	chi2+=pow((log10(data[i])-log10(bin_integration_element*totalpixel))/(0.434*(sqrt(data[i])/data[i])),2);
    	}
//    	printf("#%d chi2 =%e \n",i,chi2);
    	if(isnan(chi2))printf("data %e\n",data[i]);
		}
    printf("count=%d\n",count);

	free(flux);
	free(outshift);
	return chi2;
}


double Set_beam(void * ParamsArray){

	struct PD_params *p
				    = (struct PD_params *) ParamsArray;

	double FHWM=p->PSFresultionFWHM;
	double pixelsize_asec=p->pixelsize;

	double  b_fwhmpx = FHWM/pixelsize_asec;   //BEAM FWHM IN PIXELS

	int size=(int)b_fwhmpx*3.;
	double * bx = new double[size]; //CREATE X ARRAY TO CALC GAUSSIAN, PROBABLY 2 OR 3 TIMES THE SIZE
	for(int i=0;i<size;i++){
		bx[i] = i-(size-1)*1.0/2;
	}
	double * beam = new double[size*size];
	p->m_beam_size=size;
	p->m_beam=(double *) malloc(sizeof(double)*p->m_beam_size*p->m_beam_size);

	printf("Beam FWHM in Pixels is %f\n",b_fwhmpx);
	printf("Beam size is: %d x %d\n",p->m_beam_size,p->m_beam_size);

	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++){
	    beam[j+i*size] = exp(-(pow(bx[i],2)+pow(bx[j],2))*1.0/(2.*pow(b_fwhmpx/2.3548,2)));
		}
	}

	memcpy(p->m_beam, beam, sizeof(double)*size*size);

	free(beam);
	return 0.0;
}

double Get_dNdS_per_px(double flux, double px_size, double log10x,double log10y){

	double dnds=0.0;

	double arrary[2][8]={{log10x,-6.7,-6.3,-5.533,-4.766,-4.0,-3.25,-1.9},{log10y,15.0620,14.4342,13.4563,12.1586,10.2764,8.55,6.314}};

	if(flux>=0.0)dnds=pow(10,interp(log10(flux),8,(double *)arrary));

	double pxsz_sr = pow(px_size,2)/4.25452e10;
	return dnds*pxsz_sr;

}


double Get_Rx(double x, void * params){

	struct PD_params *p
			    = (struct PD_params *) params;



	int Ntsize=p->m_beam_size*p->m_beam_size;
	double * flux=new double[Ntsize];
    double * newy=new double[Ntsize];


    int check_count=0;
	for(int i=0;i<Ntsize;i++){

		  flux[i]=x/p->m_beam[i];

		  if((flux[i]>=p->source_min)&&(flux[i]<=p->source_max)){
			   newy[i]=Get_dNdS_per_px(flux[i],p->pixelsize,p->last_interplot_log10x,p->last_interplot_log10y);
		  		}else{
		  			newy[i]=-100.0;
		  			check_count++;
		  		}
	}

	if(check_count==Ntsize){
	//	printf("error: Get_Rx( %f )  flux > source_max !  %f \n",log10(x),Get_dNdS_per_px(x,1.5667));
    //   exit(0);
		free(flux);
		free(newy);
		return 0.0;
	}


	double rx=0.0;
	for(int i=0;i<Ntsize;i++){

	if((newy[i]>0.0)&&(p->m_beam[i]>1.0e-6)) rx=rx+newy[i]/p->m_beam[i];

	}


free(flux);
free(newy);

	return rx;
}


double CompactPD_LH(int Nbins, double * DataArray, void * ParamsArray ) {

     //int N=pow(2,18);
	int N=pow(2,18);

     double x[3][Nbins];
     	memcpy(x, DataArray, sizeof(double)*3*Nbins);


	fftw_complex *inbetween;
	double *f, *w, *out,*in;

	fftw_plan p1,p2;

	w =(double *) malloc(sizeof(double)*N);
	f =(double *) malloc(sizeof(double)*N);



	inbetween =(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N);

	in =(double *) malloc(sizeof(double)*N);
	out =(double *) malloc(sizeof(double)*N);

	//make a plan "estimate(guess) or measure(find an opptimal way)"
    p1 = fftw_plan_dft_r2c_1d(N, in, inbetween, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_c2r_1d(N, inbetween, out, FFTW_ESTIMATE);

	//FFTW_MEASURE

	struct PD_params *p
		    = (struct PD_params *) ParamsArray;



	double Delta=(p->d_max)/(N-1);

    Set_beam(ParamsArray);
    //it changes because of 1.0 show or not show in beam

    if(input_fw(FFTW_FORWARD,N,0.0,Delta,f)==-1.0) exit(0);


    double dx=Delta;


    for (int i = 0; i < N; i++)
    		{

    	       if(f[i]< p->d_min ){
    				in[i]=0.0;

    			}else{

    		       in[i]=Get_Rx(f[i],ParamsArray)*dx*amp_I;
    		 			}
    		}


	fftw_execute(p1); // execute the plan

	double nbar=inbetween[0][0];

	double noise=p->sigma_noise; //unit Jy

	double temp_real;

    if(input_fw(FFTW_BACKWARD,N,0.0,Delta,w)==-1.0) exit(0);

		  for (int i = 0; i < N; i++){
			  temp_real=inbetween[i][0];
		inbetween[i][0]=cos(inbetween[i][1]/amp_I)*exp(temp_real/amp_I-nbar/amp_I-0.5*pow(noise*w[i]*2*PI,2))*amp_II;
		inbetween[i][1]=sin(inbetween[i][1]/amp_I)*exp(temp_real/amp_I-nbar/amp_I-0.5*pow(noise*w[i]*2*PI,2))*amp_II;
		}


	fftw_execute(p2); // execute the plan


	double chi2=0.0;
	chi2=chi2_res(Nbins, x[0],x[1],x[2],Delta,N,out,p->sigma_noise);



	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);




	free(in);
	free(f);
	free(out);
	free(w);
	fftw_free(inbetween);

	free(p->m_beam);


	return chi2;

}
