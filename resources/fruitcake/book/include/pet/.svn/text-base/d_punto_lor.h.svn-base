#ifndef __D_PUNTO_LOR_H
#define __D_PUNTO_LOR_H

#include <util/utiles_vectores.h>

struct lor_params
{
  float n;
  float seno;  
  float coseno;
  float phi;
  float z_med;
  float mod_vlor;
};

float d_punto_lor(float xpix,float ypix,float zpix,float n,float seno,float coseno,float phi,float z_med,float mod_vlor);
float calcul_pes(float tbin,float dplor,float factor);
float d_punto_lor_ss(struct coord3d pto,float n,float seno,float coseno,float tan_phi,float z_med,float mod_vlor);
void calcula_fact_sum_inc(int Nbins,float n0,float inc_n,float ddet,float *fact_sum_inc);

#endif //__D_PUNTO_LOR_H
