/*
 * Main.cpp
 *
 *  Created on: Jun 9, 2016
 *      Author: ubuntu
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "Wpofd.h"



int main(void){

	int Nbins=249;
	double DataArray[3][Nbins];
	double result[Nbins];

	FILE *fpr;
	fpr=fopen("./out/histo.txt","r");

	     double temp1,temp2,temp3;
	     int Ncount=0;

		  while ((Ncount<=Nbins)&&(fscanf(fpr,"%lf  %lf  %lf\n", &temp1, &temp2, &temp3)!=EOF)){

			  if(Ncount<10)printf(" Read: %e %e %e\n",temp1, temp2, temp3);


			  DataArray[0][Ncount]=temp1;
			  DataArray[1][Ncount]=temp2;
			  DataArray[2][Ncount]=temp3;
		      Ncount++;


		  }
		  fclose(fpr);


		  printf("Ncount =%d\n", Ncount);



		 struct PD_params ParamsArray;
		 ParamsArray.d_max=1e-3; //in unit  Jy
		 ParamsArray.d_min=pow(10,-7.32);
		 ParamsArray.source_max=8.5e-5;
		 ParamsArray.source_min=ParamsArray.d_min;

		 ParamsArray.interplot_length=8;

		 double inf = -INFINITY;
		 printf("-inf is %f\n",-inf);
		 // interplot the dN/dS with {{log10x} {log10y}}
		 //double arrary[2][8]={{-7.32,-6.7,-6.3,-5.533,-4.766,-4.0,-3.25,-1.9},{16.1737,15.0620,14.4342,13.4563,12.1586,10.2764,8.55,6.314}};
		 double arrary[2][8]={{-7.69897,     -4.91506643,  -4.61439373,  -4.43842163,  -4.31354249,
				   -4.21666824,  -4.13751083,  -4.07058107},{-inf,         -inf,         11.69216698,  11.25442674,  10.94378302,10.70280301,  -inf,  -inf}};
		 ParamsArray.interplot_pointer=(double *) malloc(sizeof(double)*ParamsArray.interplot_length*2);
		 memcpy(ParamsArray.interplot_pointer, arrary, sizeof(double)*ParamsArray.interplot_length*2);

		 ParamsArray.PSFresultionFWHM=6.0;
		 ParamsArray.pixelsize=1.0;
		 ParamsArray.sigma_noise=17e-6;




	CompactPD_LH(Ncount, (double *) DataArray, result, &ParamsArray );
	CompactPD_LH(Ncount, (double *) DataArray, result, &ParamsArray );
/*	for(int i=0;i<Ncount; i++){
	printf("x-mean %lf  log10_data %lf log10_model %lf\n", (DataArray[0][i]+DataArray[1][i])*0.5,log10(DataArray[2][i]),result[i]);
	}*/

	free(ParamsArray.interplot_pointer);
	return 0;
}


