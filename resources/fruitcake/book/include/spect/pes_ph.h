#ifndef __PES_PH_H
#define __PES_PH_H

/*se calculan las contribuciones del pixel poyectado al bin de la proyeccion 
considerando q un cuadrado al proyectarse es un cuadrado
las distancias se calculan en valores absolutos*/

float calcul_pes_ph(float dz,float dy,float d1,float costheta,float tpp,float tbin,float tbind2,float tppd2);

#endif //__PES_PH_H
