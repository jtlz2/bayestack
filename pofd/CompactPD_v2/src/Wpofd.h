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
  double PSFresultionFWHM;
  double pixelsize;
  double sigma_noise;

  double * m_beam;
  int m_beam_size;

  int interplot_length;
  double * interplot_pointer;

  };



double CompactPD_LH(int Nbins, double * DataArray, double * result, void * ParamsArray );


#endif /* WPOFD_H_ */
