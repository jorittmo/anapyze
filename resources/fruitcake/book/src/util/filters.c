#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <util/filters.h>
#include <util/inpout.h>
#include <util/nrutil.h>


void cambia_channels_sinogram_hdr(struct imagen *ima_e, struct imagen *ima_s, int ncol, char *interp)
{

int i,j,k,ie,je,ke,ieref,jeref,keref,iieref,jjeref,kkeref,le,ls;
int nitemss;
float dimxe,dimze;
float x,y,z,xref,yref,zref;
float pesx,pesy,pesz,sum,sum_pes;

/* Tamaño de la matriz de entrada (debe ser igual a la de salida)*/
dimxe=ima_e->tpix*ima_e->ncol;
//dimye=ima_e->tpix*ima_e->nfil;
dimze=ima_e->tcorte*ima_e->nima;
//nitemse=ima_e->nfil*ima_e->ncol*ima_e->nima;

/* Definicion de la imagen de salida */
ima_s->ncol=ncol; /* Nuevo numero de columnas */
ima_s->nfil=ima_e->nfil;
ima_s->nima=ima_e->nima;
ima_s->Dfo=ima_e->Dfo;
ima_s->Dfd=ima_e->Dfd;
strcpy(ima_s->tipo,ima_e->tipo);

/* Nuevos tamaños de pixel y de corte debidos al cambio de elementos de matriz */
ima_s->tpix=dimxe/ima_s->ncol; /* Ojo como 1º version solo pixeles cuadrados */
ima_s->tcorte=dimze/ima_s->nima;
ima_s->offset=0;
nitemss=ima_s->nfil*ima_s->ncol*ima_s->nima;
ima_s->datos=(float *)calloc(nitemss,sizeof(float));


/* Vamos a escribir la imagen de salida con lo que le corresponde */
for(i=1;i<ima_s->nfil-1;i++)
{
for(j=1;j<ima_s->ncol-1;j++)
{
for(k=0;k<ima_s->nima;k++)
{
	
	/* pto i,j,k de salida en coordenadas espaciales */
	y=j*ima_s->tpix+0.5*ima_s->tpix;
	
	/* conversion de coordenadas espaciales x,y,z a indices de imagen de entrada por vecino mas proximo */	
	ie=i;;
	je=floor(y/ima_e->tpix);
	ke=k;
    
    //printf("i:%d j:%d k:%d  --- ie:%d je:%d ke:%d\n",i,j,k,ie,je,ke);

    
	/* conversion a indice absoluto para imagen de salida */
	ls=(k*ima_s->nfil*ima_s->ncol)+(i*ima_s->ncol)+j; 

    // Interp vecino
	if(strncmp(interp,"vecino",6)==0) 
	{
		le=(ke*ima_e->nfil*ima_e->ncol)+(ie*ima_e->ncol)+je; /* conversion a indice absoluto */
 		ima_s->datos[ls]=ima_e->datos[le];
	}

	
    // Interp novecino
	if(strncmp(interp,"novecino",8)==0) 
	{
		xref=x-0.5*ima_e->tpix;
		yref=y-0.5*ima_e->tpix;
		zref=z-0.5*ima_e->tcorte;
		
		ieref=floor(xref/ima_e->tpix);
		jeref=floor(yref/ima_e->tpix);
		keref=floor(zref/ima_e->tcorte);
		
		sum=0;
		sum_pes=0;
		for(iieref=ieref;iieref<=ieref+1;iieref++)
		{
			for(jjeref=jeref;jjeref<=jeref+1;jjeref++)
			{		
				for(kkeref=keref;kkeref<=keref+1;kkeref++)
				{
					pesx=1.-fabs(xref-(iieref*ima_e->tpix))/ima_e->tpix;
					pesy=1.-fabs(yref-(jjeref*ima_e->tpix))/ima_e->tpix;
					pesz=1.-fabs(zref-(kkeref*ima_e->tcorte))/ima_e->tcorte;			
					le=(kkeref*ima_e->nfil*ima_e->ncol)+(iieref*ima_e->ncol)+jjeref; /* conversion a indice absoluto */
					sum+=pesx*pesy*pesz*ima_e->datos[le];
					sum_pes+=pesx*pesy*pesz;
				}
			}
		}
		
		ima_s->datos[ls]=sum/sum_pes;
	}	
}
}
}
for(k=0;k<ima_s->nima;k++) 
{
for(i=0;i<ima_s->nfil;i++)
{
    ima_s->datos[(k*ima_s->nfil*ima_s->ncol)+(i*ima_s->ncol)]=1.;
    ima_s->datos[(k*ima_s->nfil*ima_s->ncol)+(i*ima_s->ncol)+ima_s->ncol-1]=1.;
}

for(j=0;j<ima_s->ncol;j++)
{
    ima_s->datos[(k*ima_s->nfil*ima_s->ncol)+j]=1.;
    ima_s->datos[(k*ima_s->nfil*ima_s->ncol)+((ima_s->nfil-1)*ima_s->ncol)+j]=1.;
}

}




}


//////////////////////////
// Actualizado para cualquier dimension de matriz y para que funcione VECINO y NOVECINO
// Pablo&Jesus 2 Mayo 2017


void cambia_matriz_imagen_hdr(char *fseq1,char *fres,int ncol,int nfil,int nima,char *interp)
{

int i,j,k,ie,je,ke,ieref,jeref,keref,iieref,jjeref,kkeref,le,ls;
int nitemss;
float dimxe,dimye,dimze;
float x,y,z,xref,yref,zref;
float pesx,pesy,pesz,sum,sum_pes;
struct imagen ima_e;
struct imagen ima_s;

lee_imagen_hdr(fseq1,&ima_e);

if(strncmp(ima_e.tipo,"1b",2)==0) strncpy(ima_e.tipo,"1B",2);	

/* Tamaño de la matriz de entrada (debe ser igual a la de salida)*/
dimxe=ima_e.tpix*ima_e.ncol;
dimye=ima_e.tpiy*ima_e.nfil;
dimze=ima_e.tcorte*ima_e.nima;
//nitemse=ima_e.nfil*ima_e.ncol*ima_e.nima;

/* Definicion de la imagen de salida */
ima_s.ncol=ncol; /* Nuevo numero de columnas */
ima_s.nfil=nfil; /* Nuevo numero de filas */
ima_s.nima=nima; /* Nuevo numero de cortes */
strcpy(ima_s.tipo,ima_e.tipo);

/* Nuevos tamaños de pixel y de corte debidos al cambio de elementos de matriz */
ima_s.tpix=dimxe/ima_s.ncol; /* Ojo como 1º version solo pixeles cuadrados */
ima_s.tpiy=dimye/ima_s.nfil; /* Ojo como 1º version solo pixeles cuadrados */
ima_s.tcorte=dimze/ima_s.nima;
ima_s.offset=0;
nitemss=ima_s.nfil*ima_s.ncol*ima_s.nima;
ima_s.datos=(float *)calloc(nitemss,sizeof(float));


/* Vamos a escribir la imagen de salida con lo que le corresponde */
for(k=1;k<ima_s.nima-1;k++)
{
for(i=1;i<ima_s.nfil-1;i++)
{
for(j=1;j<ima_s.ncol-1;j++)
{

	
	/* pto i,j,k de salida en coordenadas espaciales */
	x=j*ima_s.tpix+0.5*ima_s.tpix;
	y=i*ima_s.tpiy+0.5*ima_s.tpiy;
	z=k*ima_s.tcorte+0.5*ima_s.tcorte;
	
	/* conversion de coordenadas espaciales x,y,z a indices de imagen de entrada por vecino mas proximo */	
	je=floor(x/ima_e.tpix);
	ie=floor(y/ima_e.tpiy);
	ke=floor(z/ima_e.tcorte);
    
   	//printf("i:%d j:%d k:%d  --- ie:%d je:%d ke:%d\n",i,j,k,ie,je,ke);

    
	/* conversion a indice absoluto para imagen de salida */
	ls=(k*ima_s.nfil*ima_s.ncol)+(i*ima_s.ncol)+j; 

    // Interp vecino
	if(strncmp(interp,"vecino",6)==0) 
	{
		le=(ke*ima_e.nfil*ima_e.ncol)+(ie*ima_e.ncol)+je; /* conversion a indice absoluto */
 		ima_s.datos[ls]=ima_e.datos[le];
	}

	
    // Interp novecino
	if(strncmp(interp,"novecino",8)==0) 
	{
		xref=x-0.5*ima_e.tpix;
		yref=y-0.5*ima_e.tpiy;
		zref=z-0.5*ima_e.tcorte;
		
		jeref=floor(xref/ima_e.tpix);
		ieref=floor(yref/ima_e.tpiy);
		keref=floor(zref/ima_e.tcorte);
		
		sum=0;
		sum_pes=0;
		for(jjeref=jeref;jjeref<=jeref+1;jjeref++)
		{
			for(iieref=ieref;iieref<=ieref+1;iieref++)
			{		
				for(kkeref=keref;kkeref<=keref+1;kkeref++)
				{
					pesx=1.-fabs(xref-(jjeref*ima_e.tpix))/ima_e.tpix;
					pesy=1.-fabs(yref-(iieref*ima_e.tpiy))/ima_e.tpiy;
					pesz=1.-fabs(zref-(kkeref*ima_e.tcorte))/ima_e.tcorte;			
					le=(kkeref*ima_e.nfil*ima_e.ncol)+(iieref*ima_e.ncol)+jjeref; /* conversion a indice absoluto */
					sum+=pesx*pesy*pesz*ima_e.datos[le];
					sum_pes+=pesx*pesy*pesz;
				}
			}
		}
		
		ima_s.datos[ls]=sum/sum_pes;
	}	
}
}
}

guarda_imagen_hdr(fres,&ima_s);
}




/////////////////////////////////////////7


void aplica_filtro_mediana(struct imagen *ima_e,struct imagen *ima_s,int f)
{
int i,j,k,m,t;
float *v;
if(f==0) v=(float*)calloc(7,sizeof(float));
if(f==1) v=(float*)calloc(27,sizeof(float));


for(k=1;k<ima_e->nima-1;k++){
if(ima_e->nima==1) k=0;	
	for(j=1;j<ima_e->nfil-1;j++){
		for(i=1;i<ima_e->ncol-1;i++){

		v[0]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i]; /* punto central */
		v[1]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i];
		v[2]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i];
		v[3]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i-1];
		v[4]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i+1];
		v[5]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i];
		v[6]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i];

		if(f==0){
			t=7; /* numero de voxeles a comparar (filtro de 7 elementos)*/
			m=2; /* voxel mediano */
			}
		
		else if (f==1) {
		
		/* voxeles sobre corte actual */				
		v[7]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i-1];
		v[8]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i+1];
		v[9]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i-1];
		v[10]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i+1];

		/* voxeles sobre corte anterior */
		v[11]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i];
		v[12]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i-1];
		v[13]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i+1];
		v[14]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i];
		v[15]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i-1];
		v[16]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i+1];
		v[17]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i-1];
		v[18]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i+1];	
		
		/* voxeles sobre corte posterior */
		v[19]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i];
		v[20]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i-1];
		v[21]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i+1];
		v[22]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i];
		v[23]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i-1];
		v[24]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i+1];
		v[25]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i-1];
		v[26]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i+1];
		
		t=27; /* numero de voxeles a comparar (filtro de 27 elementos)*/
		m=13; /* voxel mediano */
		}
					
		ordena_vector(t, v);
		ima_s->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i]=v[m];
		}
	}
}
}

void aplica_filtro_cubo_cruz(struct imagen *ima_e,struct imagen *ima_s,int f)
{
	int le;
	int i,j,k,t,m;
	float *v;
	if(f==0) v=(float*)calloc(7,sizeof(float));
	if(f==1) v=(float*)calloc(27,sizeof(float));

	for(i=0;i<ima_e->nfil;i++){
		for(j=0;j<ima_e->ncol;j++){
			for(k=0;k<ima_e->nima;k++){
				
				le=(k*ima_e->nfil*ima_e->ncol)+(i*ima_e->ncol)+j; /* conversion a indice absoluto */
				//if(ima_e->datos[le]!=0) printf("%f - ",ima_e->datos[le]);
				/*Cubo de 27 elementos, mayo2006 paguiar*/
				if(i==0 || i==ima_e->nfil-1 || j==0 || j==ima_e->ncol-1 || k==0 || k==ima_e->nima-1)
				{
					ima_s->datos[le]=0;
				}
				else
				{
					
					v[0]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i]; /* punto central */
					v[1]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i];
					v[2]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i];
					v[3]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i-1];
					v[4]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i+1];
					v[5]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i];
					v[6]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i];
					t=7;
		
					if (f==1) {
		
						/* voxeles sobre corte actual */				
						v[7]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i-1];
						v[8]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i+1];
						v[9]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i-1];
						v[10]=ima_e->datos[(k*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i+1];

						/* voxeles sobre corte anterior */
						v[11]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i];
						v[12]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i-1];
						v[13]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i+1];
						v[14]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i];
						v[15]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i-1];
						v[16]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i+1];
						v[17]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i-1];
						v[18]=ima_e->datos[((k-1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i+1];	
		
						/* voxeles sobre corte posterior */
						v[19]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i];
						v[20]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i-1];
						v[21]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j-1)*ima_e->ncol)+i+1];
						v[22]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i];
						v[23]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i-1];
						v[24]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+((j+1)*ima_e->ncol)+i+1];
						v[25]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i-1];
						v[26]=ima_e->datos[((k+1)*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i+1];
		
						t=27; /* numero de voxeles a comparar (filtro de 27 elementos)*/
					}

					ima_s->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i]=0;
 					for(m=0;m<t;m++)
 					{
 					ima_s->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i]+=v[m];
 					}
					ima_s->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i]=ima_s->datos[(k*ima_e->nfil*ima_e->ncol)+(j*ima_e->ncol)+i]/t;
				}
			
			}
		}
	}
}
