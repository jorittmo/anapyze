
#ifndef __MP_I4_H
#define __MP_I4_H

#include <string.h>

typedef struct
{
  /* componentes sparse */
  float *ar;
  int *ja;
  int *ia;
  
  /* elementos sparse, no nulos pues */
  int ne;
  
  /* dimensiones de la imagen a reconstruir y tamanhos */
  int n_fil,n_col,Ntallsima;
  float tpix,lvox;
  /* numero detectores en z (num. sinogramas directos) y tamanho cada uno*/
  int ndet_z;
  float tdet_z;
  
  /* dimensiones de cada sinograma y tbin */ 
  int Nbins,Nang;
  float tbin;
  
  /* distancia entre detectores planos */
  float ddet;
  
  /* numero de subsets */
  int Nsubsets;
  
  /* diferencia máxima entre anillos */
  int axial_diff;
		  
  /* especifica si está en formato inverso, referido a pixeles de la imagen es 1 sino 0 */
  int inv;
  
  /* especifica si está aplicada la simetría sym (1:SI, 0:NO) */
  int sym;

}sparse_i4;


/* se han introducido nuevos campos - paguiar */

void projeccio(float *lc,float *q,int p,sparse_i4 *mp,int nbins);

void projeccio_indexs(float *lcr,float *q,int i0,int i1,sparse_i4 *mp);

void projeccio_seq(float *seqima,float *seqlc,sparse_i4 *mp,int MM,int nima,int nbins);

void per_mt(float *lc,float *q2,sparse_i4 *mp,int npixels);

void per_mt_index(float *lc,float *q2,int i0,int i1,sparse_i4 *mp,int npixels);

void per_mt_index_append(float *lc,float *q2,int i0,int i1,sparse_i4 *mp);

void pesos(sparse_i4 *mp,char *nom_fitxer,int *P,int *N);

void pesos_xOS(sparse_i4 *mp,char *nom_fitxer);

void escriu_mpes(char *nom_fitxer,int *PN,int *MM,sparse_i4 *mp);

void free_mpes(sparse_i4 *f);

void vector_pesos(float *vp,char *tprec,sparse_i4 *mp,int nbins);

void vector_pesosxOS(float *vp,sparse_i4 *mp,int nvox);

void vector_pesos_os_ord(float **vp_os,int OS,int itepr,sparse_i4 *mp,int npixels);

void vector_pesos_os_reord(float **vp_os,int OS,sparse_i4 *mp,int npixels,int nproj,int *ind,int P);

void ima_sensibilidad_spect(float *sens,sparse_i4 *mp);

void error_mp(int nerr,char *text);

#endif //__MP_I4_H


