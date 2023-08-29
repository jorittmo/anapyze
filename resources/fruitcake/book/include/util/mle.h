#ifndef MLE__H
#define MLE__H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <util/mp_i4.h>
#define EOS '\0'
#define maxim(a,b) ((a)>=(b)?(a):(b))

extern int MM,PN,M,N,P;
extern sparse_i4 mp;

void mle(char *fitxer_resultats,float *le,float *vp,char tipoin[3],char tipout[3],float *n_contes,int NITER,int gravar,float W,char alg[4],char *CV,int inisor,int norma);

void rec_mle(float *qa,float *qb,float *le,float *lc,float *vp,float W,float *n_contes,char alg[4],int norma,int iteracio,FILE *file5);

void alg_mle(float W,float *q,float *q2,float *vp,char alg[4]);

void mleos(char *fitxer_resultats,float *le,char tipoin[3],char tipout[3],float *n_contes,int NITER,int gravar,float W,char alg[4],char *CV,int inisor,int OS,int norma);

void rec_mle_os(float *qa,float *qb,float *le,float *lcp,float *lcr,float **vp_os,float W,float *n_contes,char alg[4],int NdOS,int *indexs,int OS,FILE *file,char *tipout,int norma);

void alg_mle_os(float W,float *q,float *q2,float *vp,char alg[4]);

void calcul_indexs(int *indexs,int OS,int NdOS);

void calcul_dif(int *dif,int OS);

void error_rec(int nerr,char *text);

#endif //__MLE_H
