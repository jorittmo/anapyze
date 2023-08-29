#ifndef POCS__H
#define POCS__H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <util/mp_i4.h>
#define EOS '\0'
#define maxim(a,b) ((a)>=(b)?(a):(b))

extern int MM,PN,M,N;
extern float *q;
extern sparse_i4 mp;

void pocs(char *fitxer_resultats,float *le,float *vp,char tipout[3],int NITER,int primera,int gravar,float W,float nc);

void rec_pocs(float *le,float *lc,float *vp,float W,float nc,int *index,char *masc);

void omple_mascara(char *masc,float *le,float *lc);

#endif //__POCS_H
