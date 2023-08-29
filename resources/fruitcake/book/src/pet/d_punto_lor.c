#include <stdio.h>
#include <math.h>
#include <util/utiles_vectores.h> 

/*===================================== D_PUNTO_LOR =======================================
Calcula la distancia (x,y,z) de un punto a una lor(a una recta)
xpix, y pix, zpix son las coordenadas del punto y los demas parametros hacen referencia
a la recta respecto la cual voy a calcular la distancia.
paguiar, Noviembre 2003*/
float d_punto_lor(float xpix,float ypix,float zpix,float n,float seno,float coseno,float phi,float z_med,float mod_vlor)
{
float u1,u2,u3,v1,v2,v3,X1,X2,X3,dplor,mod_X;

/* Vector desde punto del objeto a un punto de la lor, -nsin(teta),ncos(teta),z0 */
// LOR point
//x=n**seno
//y=n*coseno
//z=z_med

u1=xpix-n*seno;
u2=ypix-n*coseno;
//u1=xpix-n*cos(phi)*seno;
//u2=ypix-n*cos(phi)*coseno;
u3=zpix-z_med;


/* Vector director de la lor */
v1=-coseno;
v2=seno;
v3=tan(phi);



/* Producto vectorial del segmento que une el punto objeto con un punto cualquiera de la LOR y modulo */
producto_vectorial(u1,u2,u3,v1,v2,v3,&X1,&X2,&X3);
mod_X=modulo_vector(X1,X2,X3);
/* Formula de la distancia de un punto P a una recta en el espacio, |AP x v|/|v|, A pto de recta y v vector director */
dplor=mod_X;


/* OJO!!! importante, aunque interviene mod_vlor, en realidad la dplor es independiente del modulo de la lor
pues también interviene en el cálculo de mod_X */
return(dplor);
}


/*===================================== D_PUNTO_LOR =======================================
Calcula la distancia (x,y,z) de un punto a una lor(a una recta)
igual a la anterior pero trabajando con estructuras paguiar febrero2006*/
float d_punto_lor_ss(struct coord3d pto,float n,float seno,float coseno,float tan_phi,float z_med,float mod_vlor)
{
struct coord3d dist;
struct coord3d uni;
struct coord3d pvect;
float dplor,mod_X;

/* Vector desde punto del objeto a un punto de la lor, -nsin(teta),ncos(teta),z0 */
dist.x=pto.x+n*seno;
dist.y=pto.y-n*coseno;
dist.z=pto.z-z_med;

uni.x=coseno;
uni.y=seno;
uni.z=tan_phi;

/* Producto vectorial del segmento que une el punto objeto con un punto cualquiera de la LOR y modulo */
producto_vectorial_ss(dist,uni,&pvect);
mod_X=modulo_vector_ss(pvect);

/* Formula de la distancia de un punto P a una recta en el espacio, |AP x v|/|v|, A pto de recta y v vector director */
dplor=mod_X/mod_vlor;

return(dplor);
}

/*========================================= CALCUL_PES ==============================
Calcula el peso del voxel según la distancia a la LOR es mayor/menor que el tamaño 
del bin. P.Aguiar, Noviembre 2003 */
float calcul_pes(float tbin,float dplor,float factor)
{
float pes;
//printf("dplor:%f pes:%f factor:%f",dplor,pes,factor);
pes=1.-dplor/(2*tbin);
pes*=factor;
return(pes);
}


/*==================================== CALCULA FACT SUM INC ====================================*/
// calcula factores debidos a la incidencia transaxial en PET con detectores planos enfrentados
void calcula_fact_sum_inc(int Nbins,float n0,float inc_n,float ddet,float *fact_sum_inc)
{
int in;
float n,n_i,ang_i;

	for(in=0,n=n0;in<Nbins;in++,n+=inc_n)
	{
		fact_sum_inc[in]=0.;
		
		// Adding equivalent angular positions
		n_i=n;
		while(n_i<=0 && n_i>=n0) 
		{
			ang_i=atan((n_i-n)/0.5*ddet);
			fact_sum_inc[in]+=cos(ang_i*M_PI/180.);
			n_i-=inc_n;
		}
		
		n_i=n;	
		while(n_i>0 && n_i<=-n0) 
		{
			ang_i=atan((n_i-n)/0.5*ddet);
			fact_sum_inc[in]+=cos(ang_i*M_PI/180.);
			n_i+=inc_n;
		}
	}
}


