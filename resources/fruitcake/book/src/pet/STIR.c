/******************************************************************************************************

Esta librería guarda funciones relacionadas con STIR y que son utiles para generar headers, cambiar formato
y otras utilidades.
 
paguiar 2005

*******************************************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <pet/STIR.h>

/*Calcula numero maximo de imagenes en función del -segment- para span 1*/
int num_imagenes(int slices,int segment)
{
int i=1;
int num_ima=0;
int num_ima_corr=0;
int num_ima_eff=0;
num_ima=slices*((2*segment)+1);
while(i<=segment)
{ // por cada segmento que se amplia se aumentan 2*slices pero se pierden 2 de cada vez
num_ima_corr=num_ima_corr+i;
i=i+1;
}
num_ima_eff=num_ima-(2*num_ima_corr);
return(num_ima_eff);
}







/*  Rellena los campos exclusivos de sinogramas interfile a partir del valor del SPAN

---> ojo los campos propios de la imagen no los modifica paguiar 2006 

  short int nfil;
  short int ncol;
  short int nima;
  char tipo[2];
  float tpix;
  float tcorte;
  float offset;  
  float *datos;

 Campos para imagen tipo sinograma, obligatorios para formato interfile
  short int span; 
  short int nsegments;
  short nrings;
  int *axial_coord;
  int *min_rd; 
  int *max_rd; 
  float dring;
  
*/
void rellena_campos_if(struct imagen *ima,int span,int nrings)
{
	int segment_max,segment_min,nsegments,max_rd,min_rd;
	int N,ctr,n,i_ini,j_ini,i,j,i2,j2;
	ima->span=span;
	ima->nrings=nrings;
	segment_max=floor((ima->nrings-((span+1)/2))/span);
	segment_min=(-1)*segment_max;
	nsegments=2*floor((nrings-((span+1)/2))/span)+1;
	printf("segment_max:%d segment_min:%d nsegments:%d",segment_max,segment_min,nsegments);
	ima->nsegments=nsegments;
	ima->dring=82.5; /* por defecto */

	N=0;
	for(n=segment_min;n<=segment_max;n++)
	{
		ctr=0;
		max_rd=((span+1)/2)+(span*n)-1;
		min_rd=max_rd-span+1;
		
		/* relleno campos de ring difference */
		ima->max_rd[n+segment_max]=max_rd;
		ima->min_rd[n+segment_max]=min_rd;
		
		if(n<=0) 
		{
			j_ini=0;
			i_ini=max_rd*(-1);
		}
		if(n>0) 
		{
			i_ini=0;
			j_ini=min_rd;
		}
		if(n==0) 
		{
			i_ini=0;
			j_ini=0;
		}
		
		i=i_ini;j=j_ini;
		while(i<nrings && j<nrings)
		{
			i2=i;
			j2=j;
			if(n<0)
			{
				while(i2>=0 && i2<nrings && j2>=0 && j2<nrings && j2-i2>=min_rd && j2-i2<=max_rd)
				{
					//printf("(%d,%d) + ",i2,j2);
					i2++;
					j2--;
				}
				if(span==1)
				{
					i++;j++;
				}
				else
				{
					if(ctr==0) 
					{
						i++;
						ctr=1;
					}
					else 
					{
						j++;
						ctr=0;
					}
				}
			}
			
			if(n>0)
			{
				while(i2>=0 && i2<nrings && j2>=0 && j2<nrings && j2-i2>=min_rd && j2-i2<=max_rd)
				{
					//printf("(%d,%d) + ",i2,j2);
					i2--;
					j2++;
				}
				if(span==1)
				{
					i++;j++;
				}
				else
				{
					if(ctr==0) 
					{
						j++;
						ctr=1;
					}
					else 
					{
						i++;
						ctr=0;
					}
				}
			}

			if(n==0)
			{
				while(i2>=0 && i2<nrings && j2>=0 && j2<nrings && j2-i2>=min_rd && j2-i2<=max_rd)
				{
					//printf("(%d,%d) + ",i2,j2);
					i2--;
					j2++;
				}
				if(span==1)
				{
					i++;j++;
				}
				else
				{
					if(i<max_rd) 
					{
						i++;
					}
					else if(j>=nrings-max_rd)
					{
						j++;
					}
					else
					{
						if(ctr==0) 
						{
							j++;
							ctr=1;
						}
						else 
						{
							i++;
							ctr=0;
						}
					}
				}
				
			}	
			ima->axial_coord[n+segment_max]=ima->axial_coord[n+segment_max]+1;
		}
	}
}

/*Aplica transformacion axial para un span dado*/
void compresion_axial(struct imagen *ima,struct imagen *ima2,int span)
{
	int nrings,segment_max,segment_min,nsegments,nitems,k,N;
	int ctr,n,min_rd,max_rd,i_ini,j_ini,i,j,i2,j2;
	nrings=sqrt(ima->nima);
	segment_max=floor((nrings-((span+1)/2))/span);
	segment_min=(-1)*segment_max;
	nsegments=2*segment_max+1;
	nitems=ima->nfil*ima->ncol;
	N=0;
	for(n=segment_min;n<=segment_max;n++)
	{

		ctr=0;
		max_rd=((span+1)/2)+(span*n)-1;
		min_rd=max_rd-span+1;	
		
		if(n<=0) 
		{
			j_ini=0;
			i_ini=max_rd*(-1);
		}
		if(n>0) 
		{
			i_ini=0;
			j_ini=min_rd;
		}
		if(n==0) 
		{
			i_ini=0;
			j_ini=0;
		}
		
		i=i_ini;j=j_ini;
		while(i<nrings && j<nrings)
		{
			i2=i;
			j2=j;
			if(n<0)
			{
				while(i2>=0 && i2<nrings && j2>=0 && j2<nrings && j2-i2>=min_rd && j2-i2<=max_rd)
				{
					printf("(%d,%d) + ",i2,j2);
					for(k=0;k<nitems;k++) ima2->datos[N*nitems+k]+=ima->datos[(j*nrings+i)*nitems+k];
					i2++;
					j2--;
				}
				N++;
				if(span==1)
				{
					i++;j++;
				}
				else
				{
					if(ctr==0) 
					{
						i++;
						ctr=1;
					}
					else 
					{
						j++;
						ctr=0;
					}
				}
			}
			
			if(n>0)
			{
				while(i2>=0 && i2<nrings && j2>=0 && j2<nrings && j2-i2>=min_rd && j2-i2<=max_rd)
				{
					printf("(%d,%d) + ",i2,j2);
					for(k=0;k<nitems;k++) ima2->datos[N*nitems+k]+=ima->datos[(j*nrings+i)*nitems+k];
					i2--;
					j2++;
				}
				N++;
				if(span==1)
				{
					i++;j++;
				}
				else
				{
					if(ctr==0) 
					{
						j++;
						ctr=1;
					}
					else 
					{
						i++;
						ctr=0;
					}
				}
			}

			if(n==0)
			{
				while(i2>=0 && i2<nrings && j2>=0 && j2<nrings && j2-i2>=min_rd && j2-i2<=max_rd)
				{
					printf("(%d,%d) + ",i2,j2);
					for(k=0;k<nitems;k++) ima2->datos[N*nitems+k]+=ima->datos[(j*nrings+i)*nitems+k];
					i2--;
					j2++;
				}
				N++;
				if(span==1)
				{
					i++;j++;
				}
				else
				{
					if(i<max_rd) 
					{
						i++;
					}
					else if(j>=nrings-max_rd)
					{
						j++;
					}
					else
					{
						if(ctr==0) 
						{
							j++;
							ctr=1;
						}
						else 
						{
							i++;
							ctr=0;
						}
					}
				}
				
			}	
			printf("\n");
		}
	}
}


/*Aplica compresion axial simple, reduce numero axial detectors, no aplica SPAN */
void compresion_axial_simple(struct imagen *ima,struct imagen *ima2,int axial_f)
{
	int n,i,j,k1,k2,kk1,kk2,l,l2,sqrt_nima,flag_k1,flag_k2;
	ima2->nima=ima->nima/(axial_f*axial_f);
	sqrt_nima=sqrt(ima->nima);
	
	flag_k1=flag_k2=0;
	
	for(n=0;n<ima2->nfil*ima2->ncol*ima2->nima;n++)
	{
		ima2->datos[n]=0;
	}
	
	
	for(k1=0,kk1=0;k1<sqrt_nima;k1++)
	{
		for(k2=0,kk2=0;k2<sqrt_nima;k2++)
		{
			printf("adding %d (%d,%d) in %d (%d,%d)\n",(k1*sqrt_nima)+k2,k1,k2,(kk1*sqrt_nima/axial_f)+kk2,kk1,kk2);
			for(j=0;j<ima->nfil;j++)
			{
				for(i=0;i<ima->ncol;i++)
				{
					l=(((k1*sqrt_nima)+k2)*ima->ncol*ima->nfil)+(j*ima->ncol)+i;
					l2=(((kk1*sqrt_nima/axial_f)+kk2)*ima2->ncol*ima2->nfil)+(j*ima2->ncol)+i;
					ima2->datos[l2]+=ima->datos[l]/axial_f;
				}
			}
			
			flag_k2++;
			if(flag_k2==axial_f) 
			{
				kk2++;
				flag_k2=0;
			}
		}
			
		flag_k1++;
		if(flag_k1==axial_f) 
		{
			kk1++;
			flag_k1=0;
		}
		
	}
}




/*Aplica compresion angular para un mashing dado dado*/
void compresion_angular(struct imagen *ima,struct imagen *ima2,int mashing)
{
	int n,i,j,jj,k,l,l2,flag;

	ima2->nfil=ima->nfil/mashing;
	flag=0;
	
	for(n=0;n<ima2->nfil*ima2->ncol*ima2->nima;n++)
	{
		ima2->datos[n]=0.;
	}
	
	
	for(k=0;k<ima->nima;k++)
	{
		
		for(j=0,jj=0;j<ima->nfil;j++)
		{	
			for(i=0;i<ima->ncol;i++)
			{
				l=(k*ima->ncol*ima->nfil)+(j*ima->ncol)+i;
				l2=(k*ima2->ncol*ima2->nfil)+(jj*ima2->ncol)+i;
				ima2->datos[l2]+=ima->datos[l]/mashing;
				
			}
			flag++;
			if(flag==mashing) 
			{
				jj++;
				flag=0;
			}
		}
	}
}




/*Aplica compresion transaxial */
void compresion_transax(struct imagen *ima,struct imagen *ima2,int transax_f)
{
	int n,i,ii,j,k,l,l2,flag;

	ima2->nfil=ima->nfil/transax_f;
	
	for(n=0;n<ima2->nfil*ima2->ncol*ima2->nima;n++)
	{
		ima2->datos[n]=0.;
	}
	
	
	for(k=0;k<ima->nima;k++)
	{
		flag=-1;
		for(j=0;j<ima->nfil;j++)
		{	
			for(i=0,ii=0;i<ima->ncol;i++)
			{
				l=(k*ima->ncol*ima->nfil)+(j*ima->ncol)+i;
				l2=(k*ima2->ncol*ima2->nfil)+(j*ima2->ncol)+ii;
				ima2->datos[l2]+=ima->datos[l]/transax_f;
				flag++;
				if(flag==transax_f) 
				{
					ii++;
					flag=0;
				}
			}
			
		}
	}
}

// elimina el anillo del conjunto de sinogramas
void skip_ring(struct imagen *ima,struct imagen *ima2, int ring_number_to_skip)
{
int i,j,k1,k2,l,ll;
int nrings;

ll=0;
nrings=sqrt(ima->nima);

for(k1=0;k1<nrings;k1++)
{
	for(k2=0;k2<nrings;k2++)
	{	
		if(k1!=ring_number_to_skip && k2!=ring_number_to_skip )
		{	
			for(j=0;j<ima->nfil;j++)
			{	
				for(i=0;i<ima->ncol;i++)
				{
	
					l=((k2+(k1*nrings))*ima->ncol*ima->nfil)+(j*ima->ncol)+i;
					ima2->datos[ll]=ima->datos[l];
					ll++;
				}
			}
		}
	}
}

}


// elimina sinogramas con axial diff mayor que un valor dado
void change_axial_diff(struct imagen *ima,struct imagen *ima2, int axial_diff)
{
int i,j,k1,k2,l,ll;
int nrings,n;

ll=0;
nrings=sqrt(ima->nima);
n=0;
for(k1=0;k1<nrings;k1++)
{
	for(k2=0;k2<nrings;k2++)
	{	
		if(abs(k1-k2)<=axial_diff)
		{	
			//printf("\nProcesando %d-%d (%d to %d)",k1,k2,(k1*nrings)+k2,n);
			n++;
			for(j=0;j<ima->nfil;j++)
			{	
				for(i=0;i<ima->ncol;i++)
				{
	
					l=((k2+(k1*nrings))*ima->ncol*ima->nfil)+(j*ima->ncol)+i;
					ima2->datos[ll]=ima->datos[l];
					ll++;
				}
			}
		}
	}
}
//printf("\n\n");
}


