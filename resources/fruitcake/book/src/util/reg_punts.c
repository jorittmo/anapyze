#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <util/matrius.h>
#include <util/nrutil.h>
#include <util/inpout_double.h>
#include <util/inpout.h>
#include <util/reg_punts.h>

/* ====================================== GIR I DESPLACAMENT ================ */
/* lee angulos de Euler y traslaciones, calcula matriz de rotaci�n, el vector de traslaciones y aplica primero la rotaci�n y luego la traslaci�n a la 
imagen de entrada */
/*
cfalcon ?
ccrespo, paguiar: 
- modificacion para aplicar a imagenes no cuadradas. 
- inclusi�n de lectura/escritura de header Analyze
Junio-Julio 2005
*/
void girdesp3d(struct imagen *ima,double xd,double yd,double zd,double xc,double yc,double zc,double fi,double theta,double tzeta,struct imagen *ima2)
{
int i,j,k,l,i0,j0,k0,l1,l2,l3,l4,l5,l6,l7,l8,U,V,W,MM,MMM;
double x,y,z,dx,dy,dz;
double **punt, **punt_g,**mat_gir;

FILE *file;

V=ima->ncol;
U=ima->nfil;
W=ima->nima;
	
MM=ima->nfil*ima->ncol;
MMM=ima->nfil*ima->ncol*ima->nima;

	if (fi==0. && theta==0. && tzeta==0. && xd==0. && yd==0. && zd==0.)
	{
	for(i=0;i<MMM;i++) ima2->datos[i]=ima->datos[i];
	return;
	}

punt = dmatrix(1,3,1,1);		/*   reserva memoria	*/
punt_g = dmatrix(1,3,1,1);	/*   reserva memoria	*/
mat_gir = dmatrix(1,3,1,3);	/*   reserva memoria	*/

for(i=0;i<MMM;i++) ima2->datos[i]=0.;

/*la transformacion de las coordenadas del s�lido a un sistema de ejes del espacio*/

mat_gir[1][1]=cos(tzeta)*cos(fi)-cos(theta)*sin(fi)*sin(tzeta);
mat_gir[1][2]=-sin(tzeta)*cos(fi)-cos(theta)*sin(fi)*cos(tzeta);
mat_gir[1][3]=sin(theta)*sin(fi);
mat_gir[2][1]=cos(tzeta)*sin(fi)+cos(theta)*cos(fi)*sin(tzeta);
mat_gir[2][2]=-sin(tzeta)*sin(fi)+cos(theta)*cos(fi)*cos(tzeta);
mat_gir[2][3]=-sin(theta)*cos(fi);
mat_gir[3][1]=sin(theta)*sin(tzeta);
mat_gir[3][2]=sin(theta)*cos(tzeta);
mat_gir[3][3]=cos(theta);

/*el primer elemento marca las filas y el segundo las columnas*/

file=fopen("matriz_de_giro","w");

fprintf(file,"%f %f %f\n %f %f %f\n %f %f %f\n",mat_gir[1][1],mat_gir[1][2],mat_gir[1][3],mat_gir[2][1],mat_gir[2][2],mat_gir[2][3],mat_gir[3][1],mat_gir[3][2],mat_gir[3][3]);
fclose(file);

	for(k=0;k<W;k++) {       /*bucle imagen*/
	for(i=0;i<U;i++) {       /*bucle para las filas de cada corte*/
	for(j=0;j<V;j++) {       /*bucle para las columnas de cada imagen*/
   		     
   		punt[1][1]=(double)i-xc+.5-xd;
   		punt[2][1]=(double)j-yc+.5-yd;
   		punt[3][1]=(double)k-zc+.5-zd;
   		     
   		mulmat(mat_gir,3,3,punt,1,punt_g);
		 
   		x=punt_g[1][1]+xc;
   		y=punt_g[2][1]+yc;
   		z=punt_g[3][1]+zc;
   		     
   		i0=(int) floor(x-.5);
   		j0=(int) floor(y-.5);
	        k0=(int) floor(z-.5);
                
                
	        
		
	        dx=x-i0-.5;
   		dy=y-j0-.5;
	        dz=z-k0-.5;

		if(i0>-1 && i0<U && j0>-1 && j0<V && k0>-1 && k0<W){  		             
	          l=j+V*i+MM*k;
		  l1=j0+V*i0+MM*k0;
		  l2=(j0+1)+V*i0+MM*k0;
		  l3=j0+V*(i0+1)+MM*k0; 
	          l4=(j0+1)+V*(i0+1)+MM*k0;
	          l5=j0+V*i0+MM*(k0+1);
		  l6=(j0+1)+V*i0+MM*(k0+1);
		  l7=j0+V*(i0+1)+MM*(k0+1); 
	          l8=(j0+1)+V*(i0+1)+MM*(k0+1);

	        	ima2->datos[l]=ima->datos[l1]*(1.-dx)*(1.-dy)*(1.-dz);
	         	ima2->datos[l]+=ima->datos[l2]*(1.-dx)*dy*(1.-dz);
	        	ima2->datos[l]+=ima->datos[l3]*dx*(1.-dy)*(1.-dz);
	        	ima2->datos[l]+=ima->datos[l4]*dx*dy*(1.-dz);
	      	        ima2->datos[l]+=ima->datos[l5]*(1.-dx)*(1.-dy)*dz;
	        	ima2->datos[l]+=ima->datos[l6]*(1.-dx)*dy*dz;
	         	ima2->datos[l]+=ima->datos[l7]*dx*(1.-dy)*dz; 
	                ima2->datos[l]+=ima->datos[l8]*dx*dy*dz;

				 //if(ima2->datos[l] > 1.) printf("\n i0=%d j0=%d k0=%d i=%d j=%d k=%d\n", i0,j0,k0,i,j,k);	             
		}
 
 	}
	}
	} 

free_dmatrix(punt,1,3,1,1);		/*   allibera memoria	*/
free_dmatrix(punt_g,1,3,1,1);		/*   allibera memoria	*/
free_dmatrix(mat_gir,1,3,1,3);	/*   allibera memoria	*/

}


/* ====================================== GIR I DESPLACAMENT INVERS================ */
void girdesp3d_invers(struct imagen *ima,double xd,double yd,double zd,double xc,double yc,double zc,double fi,double theta,double tzeta,struct imagen *ima2)
{
int i,j,k,l,i0,j0,k0,l1,l2,l3,l4,l5,l6,l7,l8,U,V,W,MM,MMM;
double x,y,z,dx,dy,dz;
double **punt, **punt_g,**mat_gir;
FILE *file;

V=ima->nfil;
U=ima->ncol;
W=ima->nima;
MM=ima->nfil*ima->ncol;
MMM=ima->nfil*ima->ncol*ima->nima;

	if (fi==0. && theta==0. && tzeta==0. && xd==0. && yd==0. && zd==0.)
	{
	for(i=0;i<MMM;i++) ima2->datos[i]=ima->datos[i];
	return;
	}

punt = dmatrix(1,3,1,1);			/*   reserva memoria	*/
punt_g = dmatrix(1,3,1,1);		/*   reserva memoria	*/
mat_gir = dmatrix(1,3,1,3);		/*   reserva memoria	*/

for(i=0;i<MMM;i++) ima2->datos[i]=0.;

mat_gir[1][1]=cos(tzeta)*cos(fi)-cos(theta)*sin(fi)*sin(tzeta);
mat_gir[2][1]=-sin(tzeta)*cos(fi)-cos(theta)*sin(fi)*cos(tzeta);
mat_gir[3][1]=sin(theta)*sin(fi);
mat_gir[1][2]=cos(tzeta)*sin(fi)+cos(theta)*cos(fi)*sin(tzeta);
mat_gir[2][2]=-sin(tzeta)*sin(fi)+cos(theta)*cos(fi)*cos(tzeta);
mat_gir[3][2]=-sin(theta)*cos(fi);
mat_gir[1][3]=sin(theta)*sin(tzeta);
mat_gir[2][3]=sin(theta)*cos(tzeta);
mat_gir[3][3]=cos(theta);
file=fopen("matriz_de_giro_inversa","w");

fprintf(file,"%f %f %f\n %f %f %f\n %f %f %f\n",mat_gir[1][1],mat_gir[1][2],mat_gir[1][3],mat_gir[2][1],mat_gir[2][2],mat_gir[2][3],mat_gir[3][1],mat_gir[3][2],mat_gir[3][3]);
fclose(file);


	for
(k=0;k<W;k++) {
	for(i=0;i<U;i++) {
	for(j=0;j<V;j++) {
   		     
   		punt[1][1]=(double)i-xc+.5;
   		punt[2][1]=(double)j-yc+.5;
   		punt[3][1]=(double)k-zc+.5;
   		     
   		mulmat(mat_gir,3,3,punt,1,punt_g);
		 
   		x=punt_g[1][1]+xc+xd;
   		y=punt_g[2][1]+yc+yd;
   		z=punt_g[3][1]+zc+zd;
   		     
   		i0=(int)floor(x-.5);
   		j0=(int)floor(y-.5);
	        k0=(int)floor(z-.5);

	        dx=x-i0-.5;
   		dy=y-j0-.5;
	        dz=z-k0-.5;

		if(i0>-1 && i0<U-1 && j0>-1 && j0<V-1 && k0>-1 && k0<W-1){  		             
		        
		   	l=j+U*i+MM*k;
			l1=j0+U*i0+MM*k0;
		   	l2=(j0+1)+U*i0+MM*k0;
		  	l3=j0+U*(i0+1)+MM*k0; 
	                l4=(j0+1)+U*(i0+1)+MM*k0;
	                l5=j0+U*i0+MM*(k0+1);
		    	l6=(j0+1)+U*i0+MM*(k0+1);
		        l7=j0+U*(i0+1)+MM*(k0+1); 
	                l8=(j0+1)+U*(i0+1)+MM*(k0+1);

	        	ima2->datos[l]=ima->datos[l1]*(1.-dx)*(1.-dy)*(1.-dz);
	         	ima2->datos[l]+=ima->datos[l2]*(1.-dx)*dy*(1.-dz);
	        	ima2->datos[l]+=ima->datos[l3]*dx*(1.-dy)*(1.-dz);
	        	ima2->datos[l]+=ima->datos[l4]*dx*dy*(1.-dz);
	      	        ima2->datos[l]+=ima->datos[l5]*(1.-dx)*(1.-dy)*dz;
	        	ima2->datos[l]+=ima->datos[l6]*(1.-dx)*dy*dz;
	         	ima2->datos[l]+=ima->datos[l7]*dx*(1.-dy)*dz; 
	                ima2->datos[l]+=ima->datos[l8]*dx*dy*dz;	             
		}
 	}
	}
	} 

free_dmatrix(punt,1,3,1,1);		/*   allibera memoria	*/
free_dmatrix(punt_g,1,3,1,1);		/*   allibera memoria	*/
free_dmatrix(mat_gir,1,3,1,3);	/*   allibera memoria	*/

}

/*****************************************************************************************************************************************/

/* girdesp3d_2 hace lo mismo que girdesp3d pero lee prom_a (elementos de matriz de rotaci�n) en lugar de los �ngulos de Euler 
ccrespo Octubre 2006
*/
void girdesp3d_2(struct imagen *ima,float xd,float yd,float zd,float xc,float yc,float zc,float *prom_a,struct imagen *ima2)
{
int i,j,k,l,i0,j0,k0,l1,l2,l3,l4,l5,l6,l7,l8,U,V,W,MM,MMM;
float x,y,z,dx,dy,dz;
float **punt, **punt_g,**mat_gir;

FILE *file;

V=ima->ncol;
U=ima->nfil;
W=ima->nima;

MM=ima->nfil*ima->ncol;
MMM=ima->nfil*ima->ncol*ima->nima;

/*printf("\n\nvalor traslacion: %f\n",xd);
printf("\n\nvalor traslacion: %f\n",yd);
printf("\n\nvalor traslacion: %f\n",zd);*/

	
punt = matrix(1,3,1,1);		/*   reserva memoria	*/
punt_g = matrix(1,3,1,1); 	       /*   reserva memoria	*/
mat_gir = matrix(1,3,1,3);	      /*   reserva memoria	*/

for(i=0;i<MMM;i++) ima2->datos[i]=0.;

/*la transformacion de las coordenadas del s�lido a un sistema de ejes del espacio*/
mat_gir[1][1]=prom_a[0];
mat_gir[1][2]=prom_a[1];
mat_gir[1][3]=prom_a[2];
mat_gir[2][1]=prom_a[3];
mat_gir[2][2]=prom_a[4];
mat_gir[2][3]=prom_a[5];
mat_gir[3][1]=prom_a[6];
mat_gir[3][2]=prom_a[7];
mat_gir[3][3]=prom_a[8];

/*el primer elemento marca las filas y el segundo las columnas*/

file=fopen("matriz_de_giro","w");

fprintf(file,"%f %f %f\n %f %f %f\n %f %f %f\n",mat_gir[1][1],mat_gir[1][2],mat_gir[1][3],mat_gir[2][1],mat_gir[2][2],mat_gir[2][3],mat_gir[3][1],mat_gir[3][2],mat_gir[3][3]);
fclose(file);

	for(k=0;k<W;k++) {       /*bucle imagen*/
	for(i=0;i<U;i++) {       /*bucle para las filas de cada corte*/
	for(j=0;j<V;j++) {       /*bucle para las columnas de cada imagen*/
   		     
   		punt[1][1]=(float)i-xc+.5-xd;
   		punt[2][1]=(float)j-yc+.5-yd;
   		punt[3][1]=(float)k-zc+.5-zd;
   		     
   		mulmat_float(mat_gir,3,3,punt,1,punt_g);

   		x=punt_g[1][1]+xc;
   		y=punt_g[2][1]+yc;
   		z=punt_g[3][1]+zc;
 
   		i0=(int) floorf(x-.5);
   		j0=(int) floorf(y-.5);
	        k0=(int) floorf(z-.5); 
	        
		
                dx=x-i0-.5;
   		dy=y-j0-.5;
	        dz=z-k0-.5;

		if(i0>-1 && i0<U-1 && j0>-1 && j0<V-1 && k0>-1 && k0<W-1){  		             
	          l=j+U*i+MM*k;
		  l1=j0+U*i0+MM*k0;
		  l2=(j0+1)+U*i0+MM*k0;
		  l3=j0+U*(i0+1)+MM*k0; 
	          l4=(j0+1)+U*(i0+1)+MM*k0;
	          l5=j0+U*i0+MM*(k0+1);
		  l6=(j0+1)+U*i0+MM*(k0+1);
		  l7=j0+U*(i0+1)+MM*(k0+1); 
	          l8=(j0+1)+U*(i0+1)+MM*(k0+1);

	        	ima2->datos[l]=ima->datos[l1]*(1.-dx)*(1.-dy)*(1.-dz);
	         	ima2->datos[l]+=ima->datos[l2]*(1.-dx)*dy*(1.-dz);
	        	ima2->datos[l]+=ima->datos[l3]*dx*(1.-dy)*(1.-dz);
	        	ima2->datos[l]+=ima->datos[l4]*dx*dy*(1.-dz);
	      	        ima2->datos[l]+=ima->datos[l5]*(1.-dx)*(1.-dy)*dz;
	        	ima2->datos[l]+=ima->datos[l6]*(1.-dx)*dy*dz;
	         	ima2->datos[l]+=ima->datos[l7]*dx*(1.-dy)*dz; 
	                ima2->datos[l]+=ima->datos[l8]*dx*dy*dz;
		}
 
 	}
	}
	} 

free_matrix(punt,1,3,1,1);		/*   allibera memoria	*/
free_matrix(punt_g,1,3,1,1);		/*   allibera memoria	*/
free_matrix(mat_gir,1,3,1,3);	/*   allibera memoria	*/

}


