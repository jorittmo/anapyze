#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <util/valores_medios.h>

/*==================================== VMIG i DST ============================*/

void vm_ds_vector(float *vector,int nitems,float *vm,float *ds)
{
int j;
/*int index;*/

for(j=0,*vm=*ds=0;j<nitems;j++){
    *vm+=vector[j];
    *ds+=vector[j]*vector[j];
    }
*vm/=(float)nitems;
*ds/=nitems;
*ds-=*vm* *vm;
if(*ds<=0)*ds=0;
else *ds=sqrt(*ds);    
}
    
/*==================================== MAXIMS ================================*/

int maxims_vector(float *vector,int nitems,float *maxim,int primera)
{
int index,j;

for(j=primera+1,index=primera,*maxim=vector[primera];j<nitems;j++)
    if(vector[j]>*maxim){
      *maxim=vector[j];
      /*index==j;*/
	index=j;
      }
return index;
}

/*==================================== MINIMS ================================*/

int minims_vector(float *vector,int nitems,float *minim,int primera)
{
int index,j;

for(j=primera+1,index=primera,*minim=vector[primera];j<nitems;j++)
    if(vector[j]<*minim){
      *minim=vector[j];
      /*index==j;*/
	index=j;
	
      }
return index;
}

/*==================================== CANVI_SIGNRE ==========================*/

int canvi_signe(float *vector,int nitems,int primera)
{
int j,signe;
/*int index;*/
signe=vector[primera]/(fabs(vector[primera]));
for(j=primera,signe=vector[primera]/fabs(vector[primera]);j<nitems;j++)
   if(vector[j]==0 || (vector[j]/fabs(vector[j]))!=signe) return(j);

return primera;
}

/*=================================== ESCRITURA_PROMIGS ======================*/

void escritura_vm(float **vm,float **desv_st, int NPARAM, int NIMA,char *fitxer3, int nd,float *iter)

{
int i,j;
FILE *file;
   
if((file=fopen(fitxer3,"w"))==NULL) error_vm(30,fitxer3);

for(i=0;i<NIMA;i++){
    fprintf(file,"%4.0f",iter[i]);
    for(j=0;j<NPARAM;j++) fprintf (file,"\t%f",vm[j][i]);
    fprintf(file,"\n");
    }
if(nd!=0){
  for(i=0;i<NIMA;i++){
     fprintf(file,"%4.0f",iter[i]);
     for(j=0;j<NPARAM;j++) fprintf (file,"\t%f",vm[j][i]+nd*desv_st[j][i]); 
     fprintf(file,"\n");
     }
  for(i=0;i<NIMA;i++){
     fprintf(file,"%4.0f",iter[i]);
     for(j=0;j<NPARAM;j++) fprintf (file,"\t%f",vm[j][i]-nd*desv_st[j][i]); 
     fprintf(file,"\n");
     }
  }
fclose(file);
}

/*=========================================== ERRORS =========================*/

void error_vm(int nerr,char *text)
{
switch(nerr){
   case 10: printf("\nVMIG: No hi ha mes que 1 fitxer, no es poden fer promitjos\n");
            break;
   case 20: printf("\nVMIG: ERROR OBRIR F_DADES\n");
            printf("El fitxer %s no existeix o no es pot obrir",text); break;
   case 30: printf("\nVMIG: ERROR AL FITXER D'ESCRITURA\n"); break;
            printf("El fitxer %s no es pot obrir",text); break;
   default: printf("\nVMIG: Error en el numero d'error\n");
   }
exit(0);
}
