/******************************************************************************************************

Esta librería guarda funciones relacionadas con STIR y que son utiles para generar headers, cambiar formato
y otras utilidades.
 
paguiar 2005

*******************************************************************************************************/
#ifndef __STIR_H
#define __STIR_H

#include <util/inpout.h>

int num_imagenes(int slices,int segment);
void rellena_campos_if(struct imagen *ima,int span,int nrings);
void compresion_axial(struct imagen *ima,struct imagen *ima2,int span);
void compresion_axial_simple(struct imagen *ima,struct imagen *ima2,int axial_f);
void compresion_angular(struct imagen *ima,struct imagen *ima2,int mashing);
void compresion_transax(struct imagen *ima,struct imagen *ima2,int transax_f);
void skip_ring(struct imagen *ima,struct imagen *ima2, int ring_number_to_skip);
void change_axial_diff(struct imagen *ima,struct imagen *ima2, int axial_diff);
#endif //__STIR_H










