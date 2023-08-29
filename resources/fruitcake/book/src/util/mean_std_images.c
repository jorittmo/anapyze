
/*******************************************************************************
*                                                                                                           
*  Aquesta rutina calcula les imatges suma i suma de quadrats                                                
*  de n imatges   20.1.05
*  ******************************************************************************/

/* ULL! vectors suma i suma2 s'han d'inicialitzar a l'altre programa, Pablo*/
#include <math.h>
#include <stdio.h>
#include <util/mean_std_images.h>

/* Esta funcion devuelve la suma de inp y el vector suma (inicializado previamente) y la suma de cuadrados */
void suma_ima(int Nitems, float *inp, float *suma, float *suma2)
{
int i;
for (i=0;i<Nitems;i++){
	if(inp[i]>30000) inp[i]=0;
	suma[i]+= inp[i];
	suma2[i]+= inp[i]*inp[i];
	}

}

/* Esta funcion devuelve mean y std a partir de outputs de la funcion anterior */
void mean_std(int Nitems, float *suma, float *suma2, int Nima)
{
int i;
for (i=0;i<Nitems;i++){
	suma2[i]-= suma[i]*suma[i]/(float)Nima;
	suma2[i]= sqrt(suma2[i]/(float)(Nima-1));
	suma[i]/= (float)Nima;
	}
}

