#ifndef __CONVOLUCION_H
#define __CONVOLUCION_H

#include <util/inpout.h>

void tam_kernel(int *tam_kernel, float *sigma);
void gaussian_kernel(float sigma, struct imagen *ima);
void convolucion(struct imagen *ima,struct imagen *ima2,struct imagen *ima_kernel);

#endif /*__CONVOLUCION_H */




