/*calcula el centro de masas de una imagen*/
/*adaptada para utilizar desde la libreria*/
/*cris febrero 2006*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <util/nrutil.h>
#include <util/inpout.h>

//*********************************calcula centro de masas********************************************//
void cdms(struct imagen *ima,float *x,float *y,float *z)
{

int Mf,Mc,N;
int i,j,k,l;
double a,t=0;



Mf=ima->nfil;
Mc=ima->ncol;
N=ima->nima;


 *x=0;*y=0;*z=0;

   for(k=0;k<N;k++) {
   for(i=0;i<Mf;i++) {
   for(j=0;j<Mc;j++) {
	             l=j+Mc*i+Mf*Mc*k;
		     a=ima->datos[l];
	             /*a=le[l];*/
	             *x+=a*i;
        	     *y+=a*j;
 		     *z+=a*k;
 		     t+=a;
 		    }}} 	 
   *x/=t; 
   *y/=t;
   *z/=t;
    *x+=.5; 
   *y+=.5;
   *z+=.5;


}
