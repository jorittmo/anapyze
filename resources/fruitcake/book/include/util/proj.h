#ifndef _PROJ_H_
#define _PROJ_H_

void separa_tall(float *p,int n,float *le,int P,int N,int NT);

void log_projeccio(float *le1,float *le,int PN);

void scatering(float *le,int conv,int PN,int P,int N);

void filtrar_vect(char *fitxer_filtre,float *projeccio,int M,int N);

float sigma(float distancia,int colimador,float t_pixel);

void fft(float *data,int nn,int isign);

float ran1(int *idum);

float gammln(float xx);

float poidev(float xm,int *idum);

void soroll_projeccio(int ini,float *projeccio,int tamano);

void nomf(char *fitxer,char *f_res,int p);

void reduccio_pixels(float *lsp,float *le,int Psp,int P,int N);

void reduccio_angles(float *lsp,int Nsp,int Psp,int N);

void normalitzacio_nc(float *q,int nitems,float nc);

void no_negativitat(float *q,int nitems);

#endif
