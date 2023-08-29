#ifndef NRUTIL__H
#define NRUTIL__H

/**********************************************************************\
*                                                                      *
*                                                                      *
*                        N R U T I L . C                               *
*                        ===============                               *
*                                                                      *
*                                                                      *
*             Funciones para definir vectores  y matrices de dos       *
*       sub¡ndices con extremos arbitrarios para los sub¡ndices.       *
*                                                                      *
*       Se debe incluir <alloc.h> y definir el vector como puntero     *
*       (p. ej. float *a)  y las  matrices  como  punteros  dobles     *
*       (p. ej. float **a)                                             *
*                                                                      *
*       Llamadas para vectores char, int, float y double:              *
*       a = cvector (nl, nh);                                          *
*       a = ivector (nl, nh);                                          *
*       a = fvector (nl, nh);                                          *
*       a = dvector (nl, nh);                                          *
*       Siendo nl el sub¡ndice inferior y nh el superior               *
*                                                                      *
*       Llamadas para matrices char, int, float y double:              *
*       a = cmatrix (nrl, nrh, ncl, nch);                              *
*       a = imatrix (nrl, nrh, ncl, nch);                              *
*       a = fmatrix (nrl, nrh, ncl, nch);                              *
*       a = dmatrix (nrl, nrh, ncl, nch);                              *
*       Siendo nrl y nrh los sub¡ndices inferior y superior de filas   *
*       y ncl nch los de columnas.                                     *
*                                                                      *
*       Estas funciones no precisan del modelo "huge". A cambio de     *
*       esto, las matrices que generan no forman un bloque seguido.    *
*       Las  siguientes  funciones s¡ lo forman, pero  precisan del    *
*       modelo "huge":                                                 *
*                                                                      *
*       a = hcmatrix (nrl, nrh, ncl, nch);                             *
*       a = himatrix (nrl, nrh, ncl, nch);                             *
*       a = hfmatrix (nrl, nrh, ncl, nch);                             *
*       a = hdmatrix (nrl, nrh, ncl, nch);                             *
*                                                                      *
*       a debe estar declarado como    huge **a                        *
*                                                                      *
\**********************************************************************/

void nrerror (char error_text[]);

char *cvector (int nl,int nh);

unsigned char *ucvector (int nl,int nh);

unsigned short int *usvector (int nl,int nh);

short int *sivector (int nl,int nh);

int *ivector (int nl,int nh);

unsigned int *uivector (int nl,int nh);

float *vector (int nl,int nh);

double *dvector (int nl,int nh);

char **cmatrix (int nrl,int nrh,int ncl,int nch);

short int **simatrix (int nrl,int nrh,int ncl,int nch);

int **imatrix (int nrl,int nrh,int ncl,int nch);

float **matrix (int nrl,int nrh,int ncl,int nch);

double **dmatrix (int nrl,int nrh,int ncl,int nch);

float **submatrix (float **a,int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl);

float ***f3tensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh);

void free_cvector (char *v,int nl,int nh);

void free_ucvector (unsigned char *v,int nl,int nh);

void free_usvector (unsigned short int *v,int nl,int nh);

void free_sivector (short int *v,unsigned int nl,unsigned int nh);

void free_ivector (int *v,int nl,int nh);

void free_uivector (unsigned int *v,unsigned int nl,unsigned int nh);

void free_vector (float *v,int nl,int nh);

void free_dvector (double *v,int nl,int nh);

void free_cmatrix (char **m,int nrl,int nrh,int ncl,int nch);

void free_imatrix (int **m,int nrl,int nrh,int ncl,int nch);

void free_matrix (float **m,int nrl,int nrh,int ncl,int nch);

void free_dmatrix (double **m,int nrl,int nrh,int ncl,int nch);

void free_submatrix (float **b,int nrl,int nrh,int ncl,int nch);

float **convert_matrix (float *a,int nrl,int nrh,int ncl,int nch);

void free_convert_matrix (float **b,int nrl,int nrh,int ncl,int nch);

void free_f3tensor(float ***t,int nrl,int nrh,int ncl,int nch,int ndl,int ndh);

void piksrt(int n, float arr[]);

void ordena_vector(int n, float *vector);

void ordena_cuatro_vectores(int n, float *vector1,float *vector2,float *vector3,float *vector4);

int signo(float a,float limite);

void calcula_max_ima(float *seq,int nitems,float *max);

void calcula_min_ima(float *seq,int nitems,float *min);

#endif /*NRUTIL__H*/
