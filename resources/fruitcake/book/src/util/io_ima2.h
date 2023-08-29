/*******************************************************************************
*                                                                              * 
*          AQUESTA RUTINA LLEGEIX UNA IMATGE O COLECCIO D'IMATGES              *
*           EN QUALSEVOL FORMAT I HO RETORNA A UN PUNTER DOUBLE.               *
*                                                                              * 
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <util/io_ima2.h>

/*=========================================== ERROR =========================*/
void error_io(int nerr)
{
static char *formats[]={
        "\n\nEls formats s'han d'entrar segons la seguent clau:",
	"\nUnsigned char: 1b",
	"\nShort int: 2b",
	"\nInteger: i4",
	"\nFloat: fl",
	"\nDouble: db"
	};
int i;

switch(nerr){
    case 101: printf("\n\nERROR en el format d'escritura");
              for(i=0;i<6;i++) printf(formats[i]); 
              break;
    case 102: printf("\n\nError al obrir el fitxer on llegir"); break;
    case 103: printf("\n\nError en el format del fitxer de lectura");
              for(i=0;i<6;i++) printf(formats[i]);
              break;
    case 104: printf("\n\nError en el format del fitxer de la sequencia");
              for(i=0;i<6;i++) printf(formats[i]);
              break;
    case 105: printf("\n\nError al obrir el fitxer de la sequencia"); break;
    case 106: printf("\n\nError al obrir el fitxer on escriure"); break;
    default: printf("\n\nError del numero d'error en la funcio error_io()"); 
    }
    
exit(0);
}    



void input_imatge(int Mfil,int Mcol,int N,double *dq,int tipo,char *nom_fitxer)
{

int i,MMN,MMN1;
float *q;
unsigned int *iq;
unsigned short *sq;
unsigned char *cq;
FILE *file;
        	
MMN=Mfil*Mcol*N;
MMN1=MMN-1;

file=fopen(nom_fitxer,"r");
if(file==NULL) error_io(102);

switch(tipo){
       case 1:cq=ucvector(0,MMN1);
                fread(cq,1,MMN,file);
                for(i=0;i<MMN;i++) dq[i]=(double)cq[i];
                free_ucvector(cq,0,MMN1);
                break;
       case 2:sq=usvector(0,MMN1);
                fread(sq,2,MMN,file);
                for(i=0;i<MMN;i++) dq[i]=(double)sq[i];
                free_usvector(sq,0,MMN1);  
                break;
       case 3:iq=uivector(0,MMN1);
                fread(iq,4,MMN,file);
                for(i=0;i<MMN;i++) dq[i]=(double)iq[i];
                free_uivector(iq,0,MMN1);
                break; 
       case 4:q=vector(0,MMN1);
                fread(q,4,MMN,file);
                for(i=0;i<MMN;i++) dq[i]=(double)q[i];
                free_vector(q,0,MMN1);
                break;
       case 5:fread(dq,8,MMN,file);
                break;          
                }
fclose(file);
}                 




/*============================================= O U T P U T _ I M A T G E=====*/
void output_imatge(int Mfil,int Mcol,int N,double *dq,char tipo,char *nom_fitxer,double offset,double factor_escala)
{
FILE *file;

file=fopen(nom_fitxer,"w");
if(file==NULL) error_io(106);

grabar_imatge(Mfil,Mcol,N,dq,tipo,file,offset,factor_escala);

fclose(file);

}



/*============================================= A P P E N D _ I M A T G E=====*/
void append_imatge(int Mfil,int Mcol,int N,double *dq,char tipo,char *nom_fitxer,double offset,double factor_escala)
{
FILE *file;

file=fopen(nom_fitxer,"a");
if(file==NULL) error_io(106);

grabar_imatge(Mfil,Mcol,N,dq,tipo,file,offset,factor_escala);

fclose(file);

}




/*============================================ GRABAR_IMATGE =================*/
void grabar_imatge(int Mfil,int Mcol,int N,double *dq,char tipo,FILE *file,double offset,double factor_escala)
{

int i,MMN,MMN1;
float *q;
unsigned int *iq;
unsigned short *sq;
unsigned char *cq;

MMN=Mfil*Mcol*N;
MMN1=MMN-1;
for(i=0;i<MMN;i++) dq[i]=factor_escala*(dq[i]-offset);

switch(tipo){
       case 1:cq=ucvector(0,MMN1);
                for(i=0;i<MMN;i++) cq[i]=(unsigned char)dq[i];
                fwrite(cq,1,MMN,file);
                free_ucvector(cq,0,MMN1);
                break;
       case 2:sq=usvector(0,MMN1);
                for(i=0;i<MMN;i++) sq[i]=(unsigned short)dq[i];
                fwrite(sq,2,MMN,file);
                free_usvector(sq,0,MMN1);
                break;
       case 3:iq=uivector(0,MMN1);
                for(i=0;i<MMN;i++) iq[i]=(unsigned int)dq[i];
                fwrite(iq,4,MMN,file);
                free_uivector(iq,0,MMN1);
                break; 
       case 4:q=vector(0,MMN1);
                for(i=0;i<MMN;i++) q[i]=(float)dq[i];
                fwrite(q,4,MMN,file);
                free_vector(q,0,MMN1);
                break;
       case 5:fwrite(dq,8,MMN,file);
                break;          
                }
}

