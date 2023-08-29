#include <spect/pes_ph.h>

/*se calculan las contribuciones del pixel poyectado al bin de la proyeccion 
considerando q un cuadrado al proyectarse es un cuadrado
las distancias se calculan en valores absolutos*/

float calcul_pes_ph(float dz,float dy,float d1,float costheta,float tpp,float tbin,float tbind2,float tppd2)
{

float lz,ly,pes;

if (dz>tppd2+tbind2 || dy>tppd2+tbind2) return(0.);

if (dz<=tbind2-tppd2)  lz=tpp;
	
else if (dz<=tppd2-tbind2)  lz=tbin;
	
else lz=tbind2+tppd2-dz;	

if (dy<=tbind2-tppd2) ly=tpp;

else if (dy<=tppd2-tbind2) ly=tbin; 
	
else ly=tbind2+tppd2-dy;


pes=lz*ly/(tpp*tpp*d1*d1)*costheta; //normalitzat a 1 mm del pinhole

return(pes);

}

