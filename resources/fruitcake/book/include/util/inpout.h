#ifndef INPOUT__H
#define INPOUT__H

#include <stdio.h>
/*******************************************************************************
*                                                                              *
*          AQUESTES RUTINES LLEGEIXEN UNA IMATGE O COLECCIO D'IMATGES,         *
*          PROJECCIONS O COLECCIO DE PROJECCIONS EN QUALSEVOL FORMAT           *
*          I HO RETORNA A UN PUNTER FLOAT. TAMBE GRAVA IMATGES FLOAT           *
*          EN FITXERS DE QUALSEVOL FORMAT.                                     *
*									       *
*******************************************************************************/

struct imagen
{
  short int nfil;
  short int ncol;
  short int nima;
  char tipo[2];
  float tpix;
  float tcorte;
  float offset;  
  float *datos;

 /* Campos para imagen tipo sinograma, obligatorios para formato interfile */
  short int span; 
  short int nsegments;
  short nrings;
  int *axial_coord;
  int *min_rd; 
  int *max_rd; 
  float dring;

};


void token_dades_fitxer_ima(char *fitxerdades,char *nomfitxer,int *nfil, int *ncol, int *nima,char *tipoin);

void input_fitxer(int nfil,int ncol,int nima,float *le,char *nom_fitxer,char tipoin[2]);

void input_fitxer_offset(int nfil,int ncol,int nima,float *le,char *nom_fitxer,char tipoin[2],int offset);

void input_fitxer_uchar(int nfil,int ncol,int nima,unsigned char *ue,char *nom_fitxer,char tipoin[2]);

void input_llesca(int nfil,int ncol,int nima,float *le,int p,int sumar,char *nom_fitxer,char tipoin[2]);

void input_selecc_ima_de_seq(char *fitxer,char tinp[3],float *seqima,int M,int N,int primera,int nima,int salt);

void subseq(float *seqima,float *seqima2,int MM,int primera,int nima,int salt);

void output_ima(int n_fil,int n_col,float *q,FILE *file,char tipout[3]);

void output_ima3(int n_fil,int n_col,int SL,float *q,FILE *file,char *tipout);  

void output_ima_nf(int n_fil,int n_col,float *q,char *nom_fitxer,char tipout[3]);

void output_seq_nf(int n_fil,int n_col,int nima,float *q,char *nom_fitxer,char tipout[3]);

void normalitzacio_imatge_vm(float *q,float *qout,int nitems,float vm);

void normalitzacio_imatge(float *q,float *qout,int nitems);

void numeroc(float *l,float *res,int nitems);

void numerocp(float *l,float *res,int nitems);

void dispersio(float *l,int vm,float *res,int nitems);

void omple(float *ima,int nitems,float valor);

void omple_cilindre(float *ima,int nfil,int ntalls,float valor);

void omple_cilindre_mes_petit(float *ima,int nfil,int ntalls,float valor, int rest);

char *itoa(int n,char *s);

void lee_imagen_hdr(const char *nom_fitxer_hdr, struct imagen *ima);

void lee_template_hdr(const char *nom_fitxer_hdr, struct imagen *ima);

void lee_imagen_hv(const char *nom_fitxer_hv, struct imagen *ima);

void guarda_imagen_hdr(const char *nom_fitxer_hdr, struct imagen *ima);

void print_hdr(const char *nom_fitxer_hdr);

void guarda_imagen_hs(char *nombre_archivo_hs,struct imagen *ima);

void guarda_imagen_hv(char *nombre_archivo_hv,struct imagen *ima);

void error_io(int nerr,char *text);

#endif /*__INPOUT_H*/
