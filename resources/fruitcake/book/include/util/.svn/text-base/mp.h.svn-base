/* mp.h es la version anterior a mp_i4.h en la que se define mp en lugar del nuevo mp_i4, las diferencias son importantes porque 
los elementos *ja de mp_i4 son int y aqui son unsigned short int */ 


#ifndef __MP_I4_H
#define __MP_I4_H

#ifndef __MP_H
#define __MP_H

typedef struct
{
  float *ar;
  unsigned short int *ja;
  int *ia;
  size_t n_fil,n_col,ne;
}sparse;


void projeccio_old(float *lc,float *q,int p,sparse *mp,int nbins);

void projeccio_indexs_old(float *lcr,float *q,int i0,int i1,sparse *mp);

void projeccio_seq_old(float *seqima,float *seqlc,sparse *mp,int MM,int nima,int nbins);

void projeccio_seq2_old(float *seqima,float *seqlc,sparse *mp,int MM,int nima,int nbins);

void projeccio_seq3_old(float *seqima,float *seqlc,sparse *mp,int MM,int nima,int nbins);

void per_mt_old(float *lc,float *q2,sparse *mp,int npixels);

void per_mt_index_old(float *lc,float *q2,int i0,int i1,sparse *mp,int npixels);

void per_mt_index_append_old(float *lc,float *q2,int i0,int i1,sparse *mp);

void pesos_old(sparse *mp,char *nom_fitxer,int *P,int *N);

void free_mpes_old(sparse *f);

void vector_pesos_old(float *vp,char *tprec,sparse *mp,int nbins);

void vector_pesos_os_ord_old(float **vp_os,int OS,int itepr,sparse *mp,int npixels);

void vector_pesos_os_reord_old(float **vp_os,int OS,sparse *mp,int npixels,int nproj,int *ind,int P);

void error_mp_old(int nerr,char *text);
  
#endif //__MP_H
#endif //__MP_I4_H
