#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <util/enc_centros.h>

void encuentra(int ntt, int *puntos, int *vectx, int *vecty, float *v, float *minimo, float *ampli, int rad_px, int fov)
{
int j,i,d,posx,posy,num=0,primero=0,selec=1,nos=0,kk,yini,yfin,xini,xfin;
double sv,sxv,syv,min=-1.,max,distan;
float *nv;
float *les;
char f_res[]="nuevov2.img";
FILE *file2;

minimo[0]=10000.;
sv=sxv=syv=0.;
d=ntt*ntt;
nv=(float*)calloc(sizeof(float),d);
for(i=0;i<d;i++){ nv[i]=v[i]; }
//v1=(double*)calloc(8,d);


//Calculamos el pixel con mayor número de cuentas y lo dividimos entre 10.
//Usamos este valor para restarlo a la matriz de datos y así eliminar el fondo.

//Calculamos la posición x e y del c.m.

for(kk=0;kk<ntt;kk++){
	max=-1.; posx=0; posy=0;
//	printf(" \n Vamos por el intento numero %i \n ",kk);
	if (fov!=1) {
		for(i=ntt/fov;i<(fov-1)*ntt/fov;i++){	
			for(j=ntt/fov;j<(fov-1)*ntt/fov;j++){
			
				if (i>ntt/fov && i<(fov-1)*ntt/fov) {
					if (j>ntt/fov && j<(fov-1)*ntt/fov) {
						if (nv[ntt*i+j]<minimo[0] && nv[ntt*i+j]>-1.) {
							minimo[0]=nv[ntt*i+j];
						}
					}
				}
	
				if (nv[ntt*i+j]>max) {
					max=nv[ntt*i+j];
					posx=j;
					posy=i;
				}		   

			}
		}
	}else {
		for(i=4;i<ntt-4;i++){	
			for(j=4;j<ntt-4;j++){
			
				if (i>4 && i<(ntt-4)) {
					if (j>4 && j<(ntt-4)) {
						if (nv[ntt*i+j]<minimo[0] && nv[ntt*i+j]>-1.) {
							minimo[0]=nv[ntt*i+j];
						}
					}
				}
	
				if (nv[ntt*i+j]>max) {
					max=nv[ntt*i+j];
					posx=j;
					posy=i;
				}		   

			}
		}
	}

	if (nv[ntt*posy+posx+1]==0. || nv[ntt*posy+posx-1]==0. || nv[ntt*(posy+1)+posx]==0. || nv[ntt*(posy-1)+posx]==0. || nv[ntt*(posy+1)+posx+1]==0. || nv[ntt*(posy+1)+posx-1]==0. || nv[ntt*(posy-1)+posx+1]==0. || nv[ntt*(posy-1)+posx-1]==0.){
		nv[ntt*posy+posx]=-1.;

		continue;
	}



	if (primero==0) {
		min=max/3.;
		primero=1;
			printf(" \n El valor del min es: %f \n ",min);
	}

	if (min>=max)	{
		printf(" \n Min mayor que max \n ");
		file2=fopen(f_res,"wb");
		les=(float *)calloc(sizeof(float),d);

		for (i=0;i<d;i++){les[i] = (float)nv[i];}
		
		fwrite((void*)les,sizeof(float),d,file2);
		free(les);
		fclose(file2);

		return;
	}
	yini=posy-rad_px; yfin=posy+rad_px; xini=posx-rad_px; xfin=posx+rad_px;
/*	if (yini<0) yini=0;
	if (xini<0) xini=0;
	if (yfin>=ntt) yfin=ntt-1;
	if (xfin>=ntt) xfin=ntt-1;
*/

		for(i=yini;i<=yfin;i++){
			for(j=xini;j<=xfin;j++){
				distan=(double) ((j-posx)*(j-posx)+(i-posy)*(i-posy));
				distan=sqrt(distan);
				if (distan>=rad_px){
					if (nv[ntt*i+j]>max/2. || nv[ntt*i+j]==-1.){
						if (selec!=0){
							selec=0;
							nos++;
						}
					}
				}
				nv[ntt*i+j]=-1.;

			}
		}

		if (selec==1){
			vectx[num]=posx;
			vecty[num]=posy;
			ampli[num]=(float) max;
			num++;
			printf("Hola \n ");
			printf("Tenemos %i %i %f \n ", posx+1, posy+1, max); 

		} else {
			selec=1;
		}

		max=-1.;
		posx=0;
		posy=0;
		puntos[0]=num;


	if (nos>35)	{
		printf(" \n Demasiados nos \n ");
		file2=fopen(f_res,"wb");
		les=(float *)calloc(sizeof(float),d);

		for (i=0;i<d;i++){les[i] = (float)nv[i];}
		
		fwrite((void*)les,sizeof(float),d,file2);
		free(les);
		fclose(file2);
		return;
	}

}

printf(" \n Final del bucle \n ");
return;
}
