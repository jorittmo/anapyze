#ifndef FACTOR__H
#define FACTOR__H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <util/inpout.h>
#include <util/nrutil.h>


void fc_max_no_masc(int tipo,struct imagen *bas,struct imagen *ict,double *factor);

void fc_max_masc(int tipo,struct imagen *bas,struct imagen *ict,struct imagen *masc, double *factor);

void determinant(double a1,double a2,double a3,double a4,double a5,double a6,double a7,double a8,double a9,double *det);

#endif /*FACTOR__H*/
