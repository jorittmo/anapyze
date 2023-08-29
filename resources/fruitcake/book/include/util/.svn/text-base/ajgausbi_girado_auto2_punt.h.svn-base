#ifndef __AJGAUSBI_GIRADO_AUTO2_PUNT_H
#define __AJGAUSBI_GIRADO_AUTO2_PUNT_H


// Este programa hace un ajuste bidimensional utiliza una c cte que
//es el valor mínimo de la parte central de la imagen y no se adapta. 
#include <stdio.h>

/*================================CALCULO DEL CENTRO DE MASA=====================================*/
void centroide(int ntt,double tp, char nom[], double *x1cm, double *y1cm, double *sg1x, double *sg1y, int n3, int *punts, int ntalls, int tall, double xx, double yy, double *r1);

/* ===============================  CALCULA CENTRO DE GAUSSIANA  ================================= */
void calcentro3(int ntt, double vectx, double vecty, double *xcm, double *ycm, double *sgx, double *sgy, float *v, int n, double *rr);


/*===================== CALCULA amplada, mitjana i sigma guassiana amb  SIMPLEX ===================*/	
void ajgausbi(double *a,double *mrad,double *mtg,double *sigmax,double *sigmay, double *c,double *r);

/* ======================================= SUMA RESIDUS========================================= */
double suma_mgbi(double *x);

#endif /*__AJGAUSBI_GIRADO_AUTO2_PUNT_H */
