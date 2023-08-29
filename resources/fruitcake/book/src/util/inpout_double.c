#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <util/inpout_double.h>
#include <util/nrutil.h>
#define EOS '\0'
/* Esta librería es temporal a la espera de integrar estas cuatro funciones en la libreria
inpout.h substituyendo a las existentes, pues son mas generales */


/*=================================== INPUT FITXER ==========================*/

void input_fitxer_db(int nfil,int ncol,int nima,double *dq,char *nom_fitxer,char tipoin[3])
{
FILE *file;
if((file=fopen(nom_fitxer,"r"))==NULL) error_io_double(102);
llegir_fitxer(nfil,ncol,nima,dq,file,tipoin);
fclose(file);
}

/*=================================== INPUT FITXER OFFSET ==========================*/
// Esta funcion copia de archivo a un puntero saltando un offset  -paguiar Abril 2005-

void input_fitxer_db_offset(int nfil,int ncol,int nima,double *dq,char *nom_fitxer,char tipoin[3],int offset)
{
FILE *file;
if((file=fopen(nom_fitxer,"r"))==NULL) error_io_double(102);
llegir_fitxer_offset(nfil,ncol,nima,dq,file,tipoin,offset);  //paguiar, Abril2005
fclose(file);
}
      

/*------------------------------------- LLEGIR FITXER ------------------------*/

void llegir_fitxer(int nfil,int ncol,int nima,double *dq,FILE *file,char tipoin[3])
{
int nitems,j;
unsigned short int *uq;
unsigned char *cq;
int *iq;
float *fq;

nitems=nfil*ncol*nima;
if(strncmp(tipoin,"1b",2)==0){
  cq=ucvector(0,nitems-1);
  fread(cq,1,nitems,file);
  for(j=0;j<nitems;j++) *(dq+j)=(double) *(cq+j);
  free_ucvector(cq,0,nitems-1);
  }
else{
  if(strncmp(tipoin,"2b",2)==0){  
    uq=usvector(0,nitems-1);
    fread(uq,2,nitems,file);
    for(j=0;j<nitems;j++) {*(dq+j)=(double) *(uq+j);
    }
    for(j=0;j<nitems;){
    printf("%f - %d\n",*(dq+j),*(uq+j));
    j=j+(nitems/5);
    }
    free_usvector(uq,0,nitems-1);
    }  
           	else{
    	if(strncmp(tipoin,"i4",2)==0){  
     	  iq=ivector(0,nitems-1);
      	  fread(iq,4,nitems,file);
      	  for(j=0;j<nitems;j++) *(dq+j)=(double) *(iq+j);
      	  free_ivector(iq,0,nitems-1);
      	}
  	  else{
      	  if(strncmp(tipoin,"fl",2)==0) {  
      	      fq=vector(0,nitems-1);
              fread(fq,4,nitems,file);
      	      for(j=0;j<nitems;j++) *(dq+j)=(double) *(fq+j);
              free_vector(fq,0,nitems-1);
      	     }
      
      else{
        if(strncmp(tipoin,"db",2)==0){  
          fread(dq,8,nitems,file);
          } 
        else error_io_double(103);
  } } } }
}            

/*------------------------------------- LLEGIR FITXER OFFSET------------------------*/
// Se aumenta un parametro offset para poder leer imagenes con header -paguiar Abril 2005-


void llegir_fitxer_offset(int nfil,int ncol,int nima,double *dq,FILE *file,char tipoin[3],int offset)
{
int nitems,j;
unsigned short int *uq;
unsigned char *cq;
int *iq;
float *fq;


// skip header  (paguiar Abril 2005
if (fseek(file, offset, SEEK_SET) != 0) printf("\nError skipping header (llegir_fitxer_offset())\n"); 

nitems=nfil*ncol*nima;
if(strncmp(tipoin,"1b",2)==0){
  cq=ucvector(0,nitems-1);
  fread(cq,1,nitems,file);
  for(j=0;j<nitems;j++) *(dq+j)=(double) *(cq+j);
  free_ucvector(cq,0,nitems-1);
  }
else{
  if(strncmp(tipoin,"2b",2)==0){  
    uq=usvector(0,nitems-1);
    fread(uq,2,nitems,file);
    for(j=0;j<nitems;j++) {*(dq+j)=(double) *(uq+j);
    }
    for(j=0;j<nitems;){
    printf("%f - %d\n",*(dq+j),*(uq+j));
    j=j+(nitems/5);
    }
    free_usvector(uq,0,nitems-1);
    }  
           	else{
    	if(strncmp(tipoin,"i4",2)==0){  
     	  iq=ivector(0,nitems-1);
      	  fread(iq,4,nitems,file);
      	  for(j=0;j<nitems;j++) *(dq+j)=(double) *(iq+j);
      	  free_ivector(iq,0,nitems-1);
      	}
  	  else{
      	  if(strncmp(tipoin,"fl",2)==0) {  
      	      fq=vector(0,nitems-1);
              fread(fq,4,nitems,file);
      	      for(j=0;j<nitems;j++) *(dq+j)=(double) *(fq+j);
              free_vector(fq,0,nitems-1);
      	     }
      
      else{
        if(strncmp(tipoin,"db",2)==0){  
          fread(dq,8,nitems,file);
          } 
        else error_io_double(103);
  } } } }
}            

    
/*====================================== OUTPUT IMATGE =======================*/

void output_ima_db(int nfil,int ncol,double *q,char *nom_fitxer,char tipout[3])
{
FILE *file;
if((file=fopen(nom_fitxer,"w"))==NULL) error_io_double(106);
gravar_ima(nfil,ncol,q,file,tipout);
fclose(file);
}
  
/*------------------------------------ GRAVAR IMATGE --------------------------*/

void gravar_ima(int nfil,int ncol,double *q,FILE *file,char tipout[3])
{
int i,n_items;
unsigned short int *sq;
unsigned char *cq;
int *iq;
float *fq;

n_items=nfil*ncol;

if(strncmp(tipout,"1b",2)==0){
   cq=(unsigned char *)malloc(n_items);
   for(i=0;i<n_items;i++) *(cq+i)=(unsigned char) *(q+i);
   fwrite(cq,1,n_items,file);
   free(cq);
   }
else{
   if(strncmp(tipout,"2b",2)==0){
      sq=(unsigned short int *)malloc(n_items*2);
      for(i=0;i<n_items;i++) *(sq+i)=(unsigned short int) *(q+i);
      fwrite(sq,2,n_items,file);
      free(sq);
      }
   else{
      if(strncmp(tipout,"i4",2)==0){
         iq=(int *)malloc(n_items*4);
         for(i=0;i<n_items;i++) *(iq+i)=(int) *(q+i);
         fwrite(iq,4,n_items,file);
         free(iq);
         }
      else{
         if(strncmp(tipout,"fl",2)==0) {
         fq=(float *)malloc(n_items*4);
         for(i=0;i<n_items;i++) *(fq+i)=(float) *(q+i);
         fwrite(fq,4,n_items,file);
         free(fq);
         }
         
      else{
            if(strncmp(tipout,"db",2)==0) {
            fwrite(q,sizeof(double),n_items,file);
            }
            else error_io_double(101); 
   }  }  }  }     
}

/*======================================= NUMERO DE CONTES ===================*/

void ncuentas(double *l,double *res,int nitems)
{
int i;
 *res=0;
 for(i=0;i<nitems;i++) *res+=*(l+i);
}

/*=========================================== ERROR =========================*/

void error_io_double(int nerr)
{
static char *formats[]={
     "\n\nEls formats s'han d'entrar segons la seguent clau:",
	"\nUnsigned char: 1b",
	"\nShort int: 2b",
	"\nInteger: i4",
	"\nfloat: fl",
	"\nDouble: db"
	};
int i;

switch(nerr){
    case 101: printf("\n\nERROR en el format d'escriptura");
              for(i=0;i<6;i++) printf(formats[i]); 
              break;
    case 102: printf("\n\nError al obrir el fitxer on llegir"); break;
    case 103: printf("\n\nError en el format del fitxer de lectura");
              for(i=0;i<6;i++) printf(formats[i]);
              break;
    case 104: printf("\n\nError en el format del fitxer de la sequencia");
              for(i=0;i<6;i++) printf(formats[i]);
              break;
    case 105: printf("\n\nError al obrir el fitxer de la sequÉncia"); break;
    case 106: printf("\n\nError al obrir el fitxer on escriure"); break;
    default: printf("\n\nError del numero d'error en la funciÆ error_io_double()"); 
    }
    
exit(0);
} 
