#ifndef MEAN_STD_IMAGES__H
#define MEAN_STD_IMAGES__H
/*******************************************************************************
*                                                                                                           
*  Aquesta rutina calcula les imatges suma i suma de quadrats                                                
*  de n imatges   20.1.05
*  ******************************************************************************/

/* ULL! vectors suma i suma2 s'han d'inicialitzar a l'altre programa, Pablo*/


/* Esta funcion devuelve la suma de dos imagenes y la suma de cuadrados */
void suma_ima(int Nitems, float *inp, float *suma, float *suma2);

/* Esta funcion devuelve mean y std a partir de outputs de la funcion anterior */
void mean_std(int Nitems, float *suma, float *suma2, int Nima);


#endif /*__MEAN_STD_IMAGES__H*/
