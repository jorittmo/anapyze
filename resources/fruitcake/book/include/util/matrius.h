#ifndef MATRIUS__H
#define MATRIUS__H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void trasp(double **inp,int nf,int nc,double **out);
/*ESTA ES LA EXPRESION ORIGINAL EN LA QUE LAS VARIABLES SON DE TIPO DOUBLE*/
void mulmat(double **inp1,int nf1,int nc1,double **inp2,int nc2,double **out);

void mulmat_float(float **inp1,int nf1,int nc1,float **inp2,int nc2,float **out);

#endif /*MATRIUS__H*/
