#ifndef GIRA_ELEMENTOS_MATRIZ__H
#define GIRA_ELEMENTOS_MATRIZ__H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <util/inpout.h>

/* ============= GIRO Y/O DESPLAZAMIENTO LEYENDO DIRECTAMENTE LOS ELEMNTOS DE LA MATRIZ DE ROTACION ============ */

void giroydesplazo(struct imagen *ima,float xd,float yd,float zd,float xc,float yc,float zc,float  *prom_a, struct imagen *ima2);

#endif /*GIRA_ELEMENTOSMATRIZ__H*/
