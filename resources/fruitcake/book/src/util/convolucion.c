/*
Libreria de convolucion:
paguiar, Julio 2006*/

#include <stdlib.h>
#include <math.h>
#include <util/convolucion.h>


void tam_kernel(int *tam_kernel, float *sigma)
{
	int i,q;
	q=0;
	for(i=0;!q && i<10000;i++)
	{			
		if(exp(-((i*i))/(*sigma))<0.25)
		{
			*tam_kernel=i+5;
			q=1;		
		}	
	}	
}


// FUNCIONES GAUSSIAN_KERNEL Y CONVOLUCION PARA STRUCT IMAGEN (más abajo para float *q)
void gaussian_kernel(float sigma, struct imagen *ima)
{
	int i,j,k,ind;	
	float x,y,z; /*coordenandas de la gaussiana o del kernel */
	
	/* Control para que sea matriz cuadrada */
	if(ima->nfil != ima->ncol)
	{
		printf("El tamaño filas/columnas del kernel debe ser el mismo (cuadrado)\n");
		exit(0);
	}
	
	/* Control para que sea matriz cuadrada */
	if(ima->nima != ima->nfil && ima->nima != 1)
	{
		printf("El tamaño en z debe ser el mismo que en x-y para el kernel o bien 1 (convolucion 2D)\n");
		exit(0);
	}
	
	for(k=0;k<ima->nima;k++)
	{
		for(j=0;j<ima->nfil;j++)
		{
			for(i=0;i<ima->ncol;i++)
			{
				ind=(k*ima->nfil*ima->ncol)+(j*ima->ncol)+i;
				x=i-ima->ncol/2;
				y=j-ima->nfil/2;
				z=k-ima->nima/2;
				ima->datos[ind]=exp(-((x*x)+(y*y)+(z*z))/sigma);
				if(ima->datos[ind]<0.001) ima->datos[ind]=0;
			}
		}
	}
}

void convolucion(struct imagen *ima, struct imagen *ima2,struct imagen *ima_kernel)
{	
	int i,j,k,ik,jk,kk,ind,ind_i,ind1,ind2,ind3,ind_kernel,n;
	
	/*recorre imagen*/
	for(k=0;k<ima->nima;k++)
	{
		for(j=0;j<ima->nfil;j++)
		{
			for(i=0;i<ima->ncol;i++)
			{
				ind_i=(k*ima->nfil*ima->ncol)+(j*ima->ncol)+i;
				ima2->datos[ind_i]=0;	
				n=0;
				
				/*recorre kernel*/
				for(kk=-0.5-ima_kernel->nima/2;kk<0.5+ima_kernel->nima/2;kk++)
				{
					for(jk=-0.5-ima_kernel->nfil/2;jk<0.5+ima_kernel->nfil/2;jk++)
					{
						for(ik=-0.5-ima_kernel->ncol/2;ik<0.5+ima_kernel->ncol/2;ik++)
						{
							/* Indice sobre el kernel */
							ind_kernel=((kk+ima_kernel->nima/2)*ima_kernel->nfil*ima_kernel->ncol)+((jk+ima_kernel->nfil/2)*ima_kernel->ncol)+ik+ima_kernel->ncol/2;
							
							/* Indice sobre la imagen */
							ind1=k+kk;
							ind2=j+jk;
							ind3=i+ik;
							ind=(ind1*ima->nfil*ima->ncol)+(ind2*ima->ncol)+ind3; /*indice para buscar valor que entra en el promedio del pto en ind_i */
							
							if(ind1<0 || ind1>=ima->nima || ind2<0 || ind2>=ima->nfil || ind3<0 || ind3>=ima->ncol)
							{
								/*no hace nada*/
							}
							else
							{
								ima2->datos[ind_i]+=ima_kernel->datos[ind_kernel]*ima->datos[ind];
								n++;
							}		
						}
					}
				}			
			}
		}
	}

}



// FUNCIONES GAUSSIAN_KERNEL Y CONVOLUCION PARA FLOAT *Q (más arriba para struct imagen)
void gaussian_kernel_q(float sigma, float *kernel, int knfil, int kncol, int knima)
{
	int i,j,k,ind;	
	float x,y,z; /*coordenandas de la gaussiana o del kernel */
	
	/* Control para que sea matriz cuadrada */
	if(knfil != kncol)
	{
		printf("El tamaño filas/columnas del kernel debe ser el mismo (cuadrado)\n");
		exit(0);
	}
	
	/* Control para que sea matriz cuadrada */
	if(knima != knfil && knima != 1)
	{
		printf("El tamaño en z debe ser el mismo que en x-y para el kernel o bien 1 (convolucion 2D)\n");
		exit(0);
	}
	
	for(k=0;k<knima;k++)
	{
		for(j=0;j<knfil;j++)
		{
			for(i=0;i<kncol;i++)
			{
				ind=(k*knfil*kncol)+(j*kncol)+i;
				x=i-kncol/2;
				y=j-knfil/2;
				z=k-knima/2;
				kernel[ind]=exp(-((x*x)+(y*y)+(z*z))/sigma);
				if(kernel[ind]<0.001) kernel[ind]=0;
			}
		}
	}
}

void convolucion_q(float *q, float *cq,float *kernel, int nfil, int ncol, int nima, int knfil, int kncol, int knima)
{	
	int i,j,k,ik,jk,kk,ind,ind_i,ind1,ind2,ind3,ind_kernel,n;
	
	/*recorre imagen*/
	for(k=0;k<nima;k++)
	{
		for(j=0;j<nfil;j++)
		{
			for(i=0;i<ncol;i++)
			{
				ind_i=(k*nfil*ncol)+(j*ncol)+i;
				cq[ind_i]=0;	
				n=0;
				
				/*recorre kernel*/
				for(kk=-0.5-knima/2;kk<0.5+knima/2;kk++)
				{
					for(jk=-0.5-knfil/2;jk<0.5+knfil/2;jk++)
					{
						for(ik=-0.5-kncol/2;ik<0.5+kncol/2;ik++)
						{
							/* Indice sobre el kernel */
							ind_kernel=((kk+knima/2)*knfil*kncol)+((jk+knfil/2)*kncol)+ik+kncol/2;
							
							/* Indice sobre la imagen */
							ind1=k+kk;
							ind2=j+jk;
							ind3=i+ik;
							ind=(ind1*nfil*ncol)+(ind2*ncol)+ind3; /*indice para buscar valor que entra en el promedio del pto en ind_i */
							
							if(ind1<0 || ind1>=nima || ind2<0 || ind2>=nfil || ind3<0 || ind3>=ncol)
							{
								/*no hace nada*/
							}
							else
							{
								cq[ind_i]+=kernel[ind_kernel]*q[ind];
								n++;
							}		
						}
					}
				}			
			}
		}
	}

}
