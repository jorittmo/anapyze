#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void trasp(double **inp,int nf,int nc,double **out)
{
	int i,j;
    	for(i=1;i<nf+1;i++){for(j=1;j<nc+1;j++)  out[j][i]=inp[i][j];}    
}

/*ESTA ES LA EXPRESION ORIGINAL EN LA QUE LAS VARIABLES DE LA FUNCION SON DE TIPO DOUBLE*/
void mulmat(double **inp1,int nf1,int nc1,double **inp2,int nc2,double **out)
{
int i,j,k;  
for(i=1;i<nf1+1;i++){
    	for(j=1;j<nc2+1;j++){	
			     out[i][j]=0;
    	 		     for(k=1;k<nc1+1;k++)
    	 		     {
    	 		       out[i][j]+=inp1[i][k]*inp2[k][j];
    	 		     }
    	}
}     	
}

void mulmat_float(float **inp1,int nf1,int nc1,float **inp2,int nc2,float **out)
{
int i,j,k;  
for(i=1;i<nf1+1;i++){
    	for(j=1;j<nc2+1;j++){	
			     out[i][j]=0;
    	 		     for(k=1;k<nc1+1;k++)
    	 		     {
    	 		       out[i][j]+=inp1[i][k]*inp2[k][j];
    	 		     }
    	}
}     	
}
