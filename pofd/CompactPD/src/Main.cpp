/*
 * Main.cpp
 *
 *  Created on: Jun 9, 2016
 *      Author: ubuntu
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "Wpofd.h"



int main(void){

	int Nbins=249;
	double DataArray[3][Nbins];

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
		 ParamsArray.last_interplot_log10x=-7.32;
		 ParamsArray.last_interplot_log10y=16.1737;
		 ParamsArray.PSFresultionFWHM=6.0;
		 ParamsArray.pixelsize=1.0;
		 ParamsArray.sigma_noise=17e-6;




	double chi2= CompactPD_LH(Ncount, (double *) DataArray, &ParamsArray );
	printf("the calculated chi2 value is  %lf\n", chi2);


	return 0;
}


