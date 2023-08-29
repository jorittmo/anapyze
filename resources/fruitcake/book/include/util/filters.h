/*
		filters.h contiene funciones que aplican filtros a imagenes

*/
#ifndef FILTERS__H
#define FILTERS__H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <util/inpout.h>

void aplica_filtro_mediana(struct imagen *ima_e, struct imagen *ima_s, int f);
void aplica_filtro_cubo_cruz(struct imagen *ima_e, struct imagen *ima_s, int f);

#endif /*FILTERS__H*/
