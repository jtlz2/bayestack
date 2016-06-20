/*
 * Wpofd.h
 *
 *  Created on: Jun 9, 2016
 *      Author: ubuntu
 */

#ifndef WPOFD_H_
#define WPOFD_H_

struct PD_params{
  double d_max;
  double d_min;
  double source_max;
  double source_min;
  double faint_slop;
  double last_interplot_log10x;
  double last_interplot_log10y;
  double PSFresultionFWHM;
  double pixelsize;
  double sigma_noise;

  double * m_beam;
  int m_beam_size;


  };



double CompactPD_LH(int Nbins, double * DataArray, void * ParamsArray );


#endif /* WPOFD_H_ */
