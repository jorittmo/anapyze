#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <util/nrutil.h>

void nrerror (char error_text[])
{
        void exit ();

        printf ("\nNumerical Recipes run-time error...\n");
        printf ("%s\n",error_text);
        printf ("...now exiting to system...\n");
        exit (1);
}



char *cvector (int nl,int nh)
{
        char *v;

        v = (char *) malloc ((unsigned)   (nh-nl+1) *sizeof (char));
        if  (!v)  nrerror ("allocation failure in cvector () ");
        return v-nl;
}

unsigned char *ucvector (int nl,int nh)
{
        unsigned char *v;

        v = (unsigned char *) malloc ((unsigned)   (nh-nl+1) *sizeof (char));
        if  (!v)  nrerror ("allocation failure in ucvector () ");
        return v-nl;
}

unsigned short int *usvector (int nl,int nh)
{
        unsigned short int *v;

        v = (unsigned short int *) malloc ((unsigned)   (nh-nl+1) *sizeof (short int));
        if  (!v)  nrerror ("allocation failure in usvector () ");
        return v-nl;
}

short int *sivector (int nl,int nh)
{
	short int *v;

	v = (short int *) malloc ((unsigned)   (nh-nl+1) *sizeof (short int));
	if  (!v)  nrerror ("allocation failure in sivector () ");
	return v-nl;
}

int *ivector (int nl,int nh)
{
        int *v;

        v = (int *) malloc ((unsigned)   (nh-nl+1) *sizeof (int));
        if  (!v)  nrerror ("allocation failure in ivector () ");
        return v-nl;
}

unsigned int *uivector (int nl,int nh)
{
        unsigned int *v;

        v = (unsigned int *) malloc ((unsigned)   (nh-nl+1) *sizeof (int));
        if  (!v)  nrerror ("allocation failure in uivector () ");
        return v-nl;
}

float *vector (int nl,int nh)
{
        float *v;

        v = (float *) malloc ((unsigned)   (nh-nl+1) *sizeof (float));
        if  (!v)  nrerror ("allocation failure in vector () ");
        return v-nl;
}

double *dvector (int nl,int nh)
{
        double *v;

        v = (double *) malloc ((unsigned)   (nh-nl+1) *sizeof (double));
        if  (!v)  nrerror ("allocation failure in dvector () ");
        return v-nl;
}


char **cmatrix (int nrl,int nrh,int ncl,int nch)
{
        char **m;
        int    i;

        m = (char **) malloc ((unsigned)   (nrh-nrl+1) *sizeof (char*));
        if  (!m)  nrerror ("allocation failure 1 in imatrix () ");
        m -= nrl;

        for (i = nrl;i <= nrh;i++)  {
                m[i] = (char *) malloc ((unsigned)   (nch-ncl+1) *sizeof (char));
                if  (!m[i])  nrerror ("allocation failure 2 in imatrix () ");
                m[i] -= ncl;
        }
        return m;
}

short int **simatrix (int nrl,int nrh,int ncl,int nch)
{
        int i;
        short int **m;

        m = (short int **) malloc ((unsigned) (nrh-nrl+1) *sizeof (short int*));
        if  (!m)  nrerror ("allocation failure 1 in simatrix () ");
        m -= nrl;

        for (i = nrl;i <= nrh;i++)  {
                m[i] = (short int *) malloc((unsigned) (nch-ncl+1) *sizeof (short int));
                if  (!m[i])  nrerror ("allocation failure 2 in simatrix () ");
                m[i] -= ncl;
        }
        return m;
}


int **imatrix (int nrl,int nrh,int ncl,int nch)
{
        int i,**m;

        m = (int **) malloc ((unsigned)   (nrh-nrl+1) *sizeof (int*));
        if  (!m)  nrerror ("allocation failure 1 in imatrix () ");
        m -= nrl;

        for (i = nrl;i <= nrh;i++)  {
                m[i] = (int *) malloc ((unsigned)   (nch-ncl+1) *sizeof (int));
                if  (!m[i])  nrerror ("allocation failure 2 in imatrix () ");
                m[i] -= ncl;
        }
        return m;
}



float **matrix (int nrl,int nrh,int ncl,int nch)
{
        int i;
        float **m;

        m = (float **)  malloc ((unsigned)   (nrh-nrl+1) *sizeof (float*));
        if  (!m)  nrerror ("allocation failure 1 in matrix () ");
        m -= nrl;

        for (i = nrl;i <= nrh;i++)  {
                m[i] = (float *)  malloc ((unsigned)   (nch-ncl+1) *sizeof (float));
                if  (!m[i])  nrerror ("allocation failure 2 in matrix () ");
                m[i] -= ncl;
        }
        return m;
}

double **dmatrix (int nrl,int nrh,int ncl,int nch)
{
        int i;
        double **m;

        m = (double **)  malloc ((unsigned)   (nrh-nrl+1) *sizeof (double*));
        if  (!m)  nrerror ("allocation failure 1 in dmatrix () ");
        m -= nrl;

        for (i = nrl;i <= nrh;i++)  {
                m[i] = (double *)  malloc ((unsigned)   (nch-ncl+1) *sizeof (double));
                if  (!m[i])  nrerror ("allocation failure 2 in dmatrix () ");
                m[i] -= ncl;
        }
        return m;
}


float **submatrix (float **a,int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl)
{
        int i,j;
        float **m;

        m = (float **)  malloc ((unsigned)   (oldrh-oldrl+1) *sizeof (float*));
        if  (!m)  nrerror ("allocation failure in submatrix () ");
        m -= newrl;

        for (i = oldrl,j = newrl;i <= oldrh;i++,j++)  m[j] = a[i]+oldcl-newcl;

        return m;
}

float ***f3tensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{ 
   int i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
   float ***t;
   
   t=(float ***)malloc((size_t)((nrow+1)*sizeof(float**)));
   if (!t) nrerror("allocation faillure 1 in f3tensor()");
   t++;
   t-=nrl;
   
   t[nrl]=(float **)malloc((size_t)((nrow*ncol+1)*sizeof(float*)));
   if (!t[nrl]) nrerror("allocation faillure 2 in f3tensor()");
   t[nrl]+=1;
   t[nrl]-=ncl;
   
   t[nrl][ncl]=(float *)malloc((size_t)((nrow*ncol*ndep+1)*sizeof(float)));
   if (!t[nrl][ncl]) nrerror("allocation faillure 3 in f3tensor()");
   t[nrl][ncl]+=1;
   t[nrl][ncl]-=ndl;
   
   for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
   for(i=nrl+1;i<=nrh;i++){
      t[i]=t[i-1]+ncol;
      t[i][ncl]=t[i-1][ncl]+ncol*ndep;
      for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
      }
    
    return t;
}      

void free_cvector (char *v,int nl,int nh)
{
        free ((char*)   (v+nl));
}

void free_ucvector (unsigned char *v,int nl,int nh)
{
        free ((char*)   (v+nl));
}

void free_usvector (unsigned short int *v,int nl,int nh)
{
        free ((char*)   (v+nl));
}

void free_sivector (short int *v,unsigned int nl,unsigned int nh)
{
	free ((char*)   (v+nl));
}

void free_ivector (int *v,int nl,int nh)
{
        free ((char*)   (v+nl));
}

void free_uivector (unsigned int *v,unsigned int nl,unsigned int nh)
{
        free ((char*)   (v+nl));
}

void free_vector (float *v,int nl,int nh)
{
        free ((char*)   (v+nl));
}

void free_dvector (double *v,int nl,int nh)
{
        free ((char*)   (v+nl));
}


void free_cmatrix (char **m,int nrl,int nrh,int ncl,int nch)
{
        int i;

        for (i = nrh;i >= nrl;i--)  free ((char*)   (m[i]+ncl));
        free ((char*)   (m+nrl));
}


void free_imatrix (int **m,int nrl,int nrh,int ncl,int nch)
{
        int i;

        for (i = nrh;i >= nrl;i--)  free ((char*)   (m[i]+ncl));
        free ((char*)   (m+nrl));
}


void free_matrix (float **m,int nrl,int nrh,int ncl,int nch)
{
        int i;

        for (i = nrh;i >= nrl;i--)  free ((char*)   (m[i]+ncl));
        free ((char*)   (m+nrl));
}


void free_dmatrix (double **m,int nrl,int nrh,int ncl,int nch)
{
        int i;

        for (i = nrh;i >= nrl;i--)  free ((char*)   (m[i]+ncl));
        free ((char*)   (m+nrl));
}


void free_submatrix (float **b,int nrl,int nrh,int ncl,int nch)
{
        free ((char*)   (b+nrl));
}


float **convert_matrix (float *a,int nrl,int nrh,int ncl,int nch)
{
        int i,j,nrow,ncol;
        float **m;

        nrow = nrh-nrl+1;
        ncol = nch-ncl+1;
        m  =  (float **)  malloc ((unsigned)   (nrow) *sizeof (float*));
        if  (!m)  nrerror ("allocation failure in convert_matrix () ");
        m -= nrl;
        for (i = 0,j = nrl;i <= nrow-1;i++,j++)  m[j] = a+ncol*i-ncl;
        return m;
}


void free_convert_matrix (float **b,int nrl,int nrh,int ncl,int nch)
{
        free ((char*)   (b+nrl));
}


void free_f3tensor(float ***t,int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{
   
    free((char*)(t[nrl][ncl]+ndl-1));
    free((char*)(t[nrl]+ncl-1));
    free((char*)(t+nrl-1));
    
}


// Numerical recipes: Sin validar (ver???)  11-Marzo-2005
void piksrt(int n, float arr[])
{
int i,j;
float a;

for(j=1;j<n-1;j++){
a=arr[j];
i=j-1;
while(i>=0 && arr[i]>a){
arr[i+1]=arr[i];
i--;
}
arr[i+1]=a;
}
}

/* se ordenan componentes de un vector , n- numero de componentes */
void ordena_vector(int n, float *vector)
{
int i,j;
float temp;

    for(i=0; i<(n-1); i++) {
        for (j=i+1; j<n; j++) {
            if(vector[j]<vector[i]) {
                temp=vector[j];
                vector[j]=vector[i];
                vector[i]=temp;
            }
        }
    }
}

// Se ordenan vector 2-3-4 segun valores crecientes del vector 1
void ordena_cuatro_vectores(int n, float *vector1,float *vector2,float *vector3,float *vector4)
{
int i,j;
float temp1,temp2,temp3,temp4;

    for(i=0; i<(n-1); i++) {
        for (j=i+1; j<n; j++) {
            if(vector1[j]<vector1[i]) {
                temp1=vector1[j];
			 temp2=vector2[j];
			 temp3=vector3[j];
			 temp4=vector4[j];
			 vector1[j]=vector1[i];
                	 vector2[j]=vector2[i];
                	 vector3[j]=vector3[i];
			 vector4[j]=vector4[i];
			 vector1[i]=temp1;
			 vector2[i]=temp2;
			 vector3[i]=temp3;
			 vector4[i]=temp4;
            }
        }
    }
}
//Otra forma: #define SIGNE(a) (a<-EPSILON?-1:(a>EPSILON?1:0))
// se traza el intervalo --return -1----(-limite)---return 0----(limite)----return 1-----
int signo(float a,float limite)
{
int e;
if(a<-limite) e=-1;
if(a>limite) e=1;
if(a>=-limite && a<=limite) e=0;
return(e);
}

// calcula maximo de una imagen, forma muy lenta, MEJORAR con numerical recipes paguiar abril-2005
void calcula_max_ima(float *seq,int nitems,float *max)
{
int i;
*max=seq[0];
for(i=1;i<=nitems-1;i++){
if(seq[i]>=*max) *max=seq[i];
}
}

// calcula maximo de una imagen, forma muy lenta, MEJORAR con numerical recipes paguiar abril-2005
void calcula_min_ima(float *seq,int nitems,float *min)
{
int i;
*min=seq[0];
for(i=1;i<=nitems-1;i++){
if(seq[i]<=*min) *min=seq[i];
}
}

// calcula el histograma 1D de una imagen 2D (solo valores positivos)
int *histogram (struct imagen *ima,float interv)
{
	int i,k,ctrl,nitems,*hist;
	float tam_interv,limit;
	float min[1],max[1];

	nitems=ima->nfil*ima->ncol*ima->nima;
	hist=(int *)calloc(interv,sizeof(int)); /*calloc inicializa los valores a cero*/

	calcula_max_ima(ima->datos,nitems,max);
	calcula_min_ima(ima->datos,nitems,min);
	if(min[0]<0 && min[0]!=-999)
	{
		printf("Error: negative values in histogram min:%f\n",min[0]);
		exit(1);
	}

	if(max[0]==0 && min[0]==0)
	{
		printf("Error: only zero values\n\n");
		exit(1);
	}
	

	// Creating histogram for positive values (-999 values are skipped because are out of a predefined ROI)
	if(min[0]>=0 || min[0]==-999)
	{
		if(min[0]==-999)
		{
			if(max[0]>4000) printf("Warning some values >4000 (max: %f)\n",max[0]);
			min[0]=0;
			max[0]=4000;
		}

		tam_interv=(max[0]-min[0])/interv;
		k=0;
		ctrl=0;
		limit=tam_interv;

		for(i=0;i<nitems;i++)
		{	
			// Parece redundante el control min/max pero algunas imagenes de SPM tienen -nan y es necesario no procesarlos (son ceros en realidad) - jul2012
			if(ima->datos[i]!=-999 && ima->datos[i]>=min[0] && ima->datos[i]<=max[0])
			{
				while(ctrl<=0)
				{	
					if(ima->datos[i]<=limit+(0.00001*limit)) // To be validated (DEC-2012) ****
					{
						ctrl=1;
						hist[k]=hist[k]+1;
					}
					k++;
					limit=limit+tam_interv;
					
				}
				ctrl=0;
				k=0;
				limit=tam_interv;
			}
		}
	}
	return(hist);
}


