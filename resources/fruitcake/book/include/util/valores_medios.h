#ifndef VALORES__MEDIOS__H
#define VALORES__MEDIOS__H

#include <stdio.h>

void vm_ds_vector(float *vector,int nitems,float *vm,float *ds);
int maxims_vector(float *vector,int nitems,float *maxim,int primera);
int minims_vector(float *vector,int nitems,float *minim,int primera);
int canvi_signe(float *vector,int nitems,int primera);
void escritura_vm(float **vm,float **desv_st, int NPARAM, int NIMA,char *fitxer3, int nd,float *iter);
void error_vm(int nerr,char *text);

#endif /*VALORES_MEDIOS_H*/
