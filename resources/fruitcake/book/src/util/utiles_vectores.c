/*
Libreria de funciones utiles para vectores:
paguiar, Noviembre 2003*/

#include <stdlib.h>
#include <math.h>
#include <util/utiles_vectores.h>

/*============================= PRODUCTO VECTORIAL =================================================
Calcula el producto vectorial de dos vectores v1 x v2, el resultado es otro vector*/
void producto_vectorial(float u1,float u2,float u3,float v1,float v2,float v3, float *R1,float *R2,float *R3)
{
/*Calculo el determinante de dos vectores v1 y v2*/
*R1=u2*v3-u3*v2;
*R2=u3*v1-u1*v3;
*R3=u1*v2-u2*v1;
}

/*==================================== MODULO_VECTOR ==========================================================
Funcion que calcula el modulo de un vector */

float modulo_vector(float v1,float v2,float v3)
{
float modulo;
modulo=sqrt((v1*v1)+(v2*v2)+(v3*v3));
return modulo;
}

/*============================= PRODUCTO VECTORIAL =================================================
Calcula el producto vectorial de dos vectores v1 x v2, el resultado es otro vector modificada para struct paguiar feb2006*/
void producto_vectorial_ss(struct coord3d u,struct coord3d v, struct coord3d *r)
{
/*Calculo el determinante de dos vectores v1 y v2*/
(*r).x=u.y*v.z-u.z*v.y;
(*r).y=u.z*v.x-u.x*v.z;
(*r).z=u.x*v.y-u.y*v.x;
}

/*==================================== MODULO_VECTOR ==========================================================
Funcion que calcula el modulo de un vector modificada para struct paguiar feb2006*/

float modulo_vector_ss(struct coord3d v)
{
float modulo;
modulo=sqrt((v.x*v.x)+(v.y*v.y)+(v.z*v.z));
return modulo;
}
