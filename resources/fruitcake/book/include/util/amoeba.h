#ifndef __AMOEBA_H
#define __AMOEBA_H

#include <math.h>

void amoeba(double **p,double *y,int ndim,double ftol,double funk(),int *nfunk);


double amotry(double **p,double *y,double *psum,int ndim,double funk(),int ihi,int *nfunk,double fac);


#endif /*__AMOEBA_H */
