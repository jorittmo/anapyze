/*******************************************************************************
*                                                                              *
*          AQUESTA LLIBRERIA CONTE UNA SERIE DE RUTINES PER CONVERTIR          *
*          SEQUENCIES D'IMATGES EN UNA SOLA IMATGE, PER SEPARAR UNA IMATGE     *
*          D'UNA SEQUENCIA O UNA LLESCA O TALL D'UNA SEQUENCIA D'IMATGES       *
*          DE FITXERS DE QUALSEVOL FORMAT. FORMATS 1b, 2b, i4, fl, db.         *
*									       *
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <util/sequen.h>
#define EOS '\0'

/*========================================SEPARA_1b===========================*/

void separa_imatge_1b(unsigned char *seq,unsigned char *q,int NFIL,int NCOL,int NIMA,int n)
{
int i,nitems,nnitems;

nitems=NFIL*NCOL;
nnitems=n*nitems;

for(i=0;i<nitems;i++) *(q+i)=*(seq+nnitems+i);
}

/*---------------------------------------separa_2b----------------------------*/

void separa_imatge_2b(unsigned short int *seq,unsigned short int *q,int NFIL,int NCOL,int NIMA,int n)
{
int i,nitems,nnitems;

nitems=NFIL*NCOL;
nnitems=n*nitems;

for(i=0;i<nitems;i++) *(q+i)=*(seq+nnitems+i);
}

/*---------------------------------------separa_i4----------------------------*/

void separa_imatge_i4(int *seq,int *q,int NFIL,int NCOL,int NIMA,int n)
{
int i,nitems,nnitems;

nitems=NFIL*NCOL;
nnitems=n*nitems;

for(i=0;i<nitems;i++) *(q+i)=*(seq+nnitems+i);
}

/*---------------------------------------separa_fl----------------------------*/

void separa_imatge_fl(float *seq,float *q,int NFIL,int NCOL,int NIMA,int n)
{
int i,nitems,nnitems;

nitems=NFIL*NCOL;
nnitems=n*nitems;

for(i=0;i<nitems;i++) *(q+i)=*(seq+nnitems+i);
}

/*---------------------------------------separa_db----------------------------*/

void separa_imatge_db(double *seq,double *q,int NFIL,int NCOL,int NIMA,int n)
{
int i,nitems,nnitems;

nitems=NFIL*NCOL;
nnitems=n*nitems;

for(i=0;i<nitems;i++) *(q+i)=*(seq+nnitems+i);
}

   
/*=======================================LLESCA_1b============================*/

void separa_llesca_1b(unsigned char *seq,unsigned char *q,int NFIL,int NCOL,int NIMA,int n)
{
int i,j,nitems,nncol;

nitems=NFIL*NCOL;
nncol=n*NCOL;

for(i=0;i<NIMA;i++) for(j=0;j<NCOL;j++) *(q+i*NCOL+j)=*(seq+i*nitems+nncol+j);
}
   
/*---------------------------------------llesca_2b----------------------------*/
   
void separa_llesca_2b(short int *seq,short int *q,int NFIL,int NCOL,int NIMA,int n)
{
int i,j,nitems,nncol;

nitems=NFIL*NCOL;
nncol=n*NCOL;

for(i=0;i<NIMA;i++) for(j=0;j<NCOL;j++) *(q+i*NCOL+j)=*(seq+i*nitems+nncol+j);
}

/*---------------------------------------llesca_i4----------------------------*/

void separa_llesca_i4(int *seq,int *q,int NFIL,int NCOL,int NIMA,int n)
{
int i,j,nitems,nncol;

nitems=NFIL*NCOL;
nncol=n*NCOL;

for(i=0;i<NIMA;i++) for(j=0;j<NCOL;j++) *(q+i*NCOL+j)=*(seq+i*nitems+nncol+j);
}
  
/*---------------------------------------llesca_fl----------------------------*/

void separa_llesca_fl(float *seq,float *q,int NFIL,int NCOL,int NIMA,int n)
{
int i,j,nitems,nncol;

nitems=NFIL*NCOL;
nncol=n*NCOL;

for(i=0;i<NIMA;i++) for(j=0;j<NCOL;j++) *(q+i*NCOL+j)=*(seq+i*nitems+nncol+j);
}

/*---------------------------------------llesca_db----------------------------*/

void separa_llesca_db(double *seq,double *q,int NFIL,int NCOL,int NIMA,int n)
{
int i,j,nitems,nncol;

nitems=NFIL*NCOL;
nncol=n*NCOL;

for(i=0;i<NIMA;i++) for(j=0;j<NCOL;j++) *(q+i*NCOL+j)=*(seq+i*nitems+nncol+j);
}
   
/*========================================SEQUEN_1b===========================*/


void sequencia_a_ima_1b(unsigned char *seq,unsigned char *q,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay)
{
int i,j,k,l,nitems,nixnit,ncnix,nitcdn,nitprim,index,index2,nitemax;

nitems=NFIL*NCOL;
ncnix=NCOL*nimax;
nixnit=nimax*nitems;
nitcdn=nitems*cada_n;
nitprim=nitems*primera;
nitemax=nitems*NIMA;

for(k=0;k<nimay;k++) 
  for(l=0;l<NFIL;l++)
    for(j=0;j<nimax;j++)
      for(i=0;i<NCOL;i++){
        index=j*NCOL+k*nixnit+l*ncnix+i;
        if((index2=nitprim+nitcdn*(j+nimax*k)+l*NCOL+i)>=nitemax) *(q+index)=0;
        else *(q+index)= *(seq+index2);
        }
}

/*--------------------------------------sequen_2b-----------------------------*/

void sequencia_a_ima_2b(short int *seq,short int *q,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay)
{
int i,j,k,l,nitems,nixnit,ncnix,nitcdn,nitprim,index,index2,nitemax;

nitems=NFIL*NCOL;
ncnix=NCOL*nimax;
nixnit=nimax*nitems;
nitcdn=nitems*cada_n;
nitprim=nitems*primera;
nitemax=nitems*NIMA;

for(k=0;k<nimay;k++) 
  for(l=0;l<NFIL;l++)
    for(j=0;j<nimax;j++)
      for(i=0;i<NCOL;i++){
        index=j*NCOL+k*nixnit+l*ncnix+i;
        if((index2=nitprim+nitcdn*(j+nimax*k)+l*NCOL+i)>=nitemax) *(q+index)=0;
        else *(q+index)= *(seq+index2);
        }
}

/*--------------------------------------sequen_i4-----------------------------*/

void sequencia_a_ima_i4(int *seq,int *q,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay)
{
int i,j,k,l,nitems,nixnit,ncnix,nitcdn,nitprim,index,index2,nitemax;

nitems=NFIL*NCOL;
ncnix=NCOL*nimax;
nixnit=nimax*nitems;
nitcdn=nitems*cada_n;
nitprim=nitems*primera;
nitemax=nitems*NIMA;

for(k=0;k<nimay;k++) 
  for(l=0;l<NFIL;l++)
    for(j=0;j<nimax;j++)
      for(i=0;i<NCOL;i++){
        index=j*NCOL+k*nixnit+l*ncnix+i;
        if((index2=nitprim+nitcdn*(j+nimax*k)+l*NCOL+i)>=nitemax) *(q+index)=0;
        else *(q+index)= *(seq+index2);
        }
}

/*--------------------------------------sequen_fl-----------------------------*/

void sequencia_a_ima_fl(float *seq,float *q,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay)
{
int i,j,k,l,nitems,nixnit,ncnix,nitcdn,nitprim,index,index2,nitemax;

nitems=NFIL*NCOL;
ncnix=NCOL*nimax;
nixnit=nimax*nitems;
nitcdn=nitems*cada_n;
nitprim=nitems*primera;
nitemax=nitems*NIMA;

for(k=0;k<nimay;k++) 
  for(l=0;l<NFIL;l++)
    for(j=0;j<nimax;j++)
      for(i=0;i<NCOL;i++){
        index=j*NCOL+k*nixnit+l*ncnix+i;
        if((index2=nitprim+nitcdn*(j+nimax*k)+l*NCOL+i)>=nitemax) *(q+index)=0;
        else *(q+index)= *(seq+index2);
        }
}        

/*--------------------------------------sequen_db-----------------------------*/

void sequencia_a_ima_db(double *seq,double *q,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay)
{
int i,j,k,l,nitems,nixnit,ncnix,nitcdn,nitprim,index,index2,nitemax;

nitems=NFIL*NCOL;
ncnix=NCOL*nimax;
nixnit=nimax*nitems;
nitcdn=nitems*cada_n;
nitprim=nitems*primera;
nitemax=nitems*NIMA;

for(k=0;k<nimay;k++) 
  for(l=0;l<NFIL;l++)
    for(j=0;j<nimax;j++)
      for(i=0;i<NCOL;i++){
        index=j*NCOL+k*nixnit+l*ncnix+i;
        if((index2=nitprim+nitcdn*(j+nimax*k)+l*NCOL+i)>=nitemax) *(q+index)=0;
        else *(q+index)= *(seq+index2);
        }
}

/*----------------------------------ima_a_seq_fl------------------------------*/

void ima_a_seq_fl(float *q,float *seq,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay)
{
int i,j,k,l,nitems,nixnit,ncnix,nitcdn,nitprim,index,index2,nitemax;

nitems=NFIL*NCOL;
ncnix=NCOL*nimax;
nixnit=nimax*nitems;
nitcdn=nitems*cada_n;
nitprim=nitems*primera;
nitemax=nitems*NIMA;

for(k=0;k<nimay;k++) 
  for(l=0;l<NFIL;l++)
    for(j=0;j<nimax;j++)
      for(i=0;i<NCOL;i++){
        index=j*NCOL+k*nixnit+l*ncnix+i;
        index2=nitprim+nitcdn*(j+nimax*k)+l*NCOL+i;
        *(seq+index2)= *(q+index);
        }
} 


/*--------------------------------transposar matriu-----------------------------*/

void trasposa_matriu_fl(float *q1,float *q2,int Nfil,int Ncol)
{
int indexq1,indexq2,i,j;

for(i=0,indexq1=0;i<Nfil;i++){
     indexq2=i;
     for(j=0;j<Ncol;j++,indexq1++,indexq2+=Nfil){
          q2[indexq2]=q1[indexq1];
          }
     }
}  

/*------------------------------copia_vectors------------------------------*/

void copia_vectors_fl(float *seq,float *seq2,int Nitems)
{
register int i;

for(i=0;i<Nitems;i++) seq[i]=seq2[i];

}
