#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <pet/pocs3d_pet.h>
#include <util/mp_i4.h>
#include <util/inpout.h>
#include <pet/mle3dxOS_pet.h>
#define EOS '\0'
#define maxim(a,b) ((a)>=(b)?(a):(b))

/*============================= MASCARA ====================================*/

void omple_mascara(int *masc,float *le,float *lc)
{
register int i,j;

for(i=0;i<Nvox;i++) *(masc+i)=1;
for(i=0;i<Nbt;i++) if(*(le+i)==0) for(j=mp.ia[i];j<mp.ia[i+1];j++) *(masc+mp.ja[j])=0;
}

/*==================================== POCS ==================================*/

void pocs(char *fitxer_resultats,float *le,float *vp,char tipout[3],int NITER,float W,float nc)
{  
int i,j,k,iteracio,*index;
FILE *file;
float *lcp;
int *masc;

if((lcp=(float *)calloc(Nbt,sizeof(float)))==NULL) error_rec(24,"");
if((masc=(int *)calloc(Nvox,sizeof(int)))==NULL) error_rec(24,"");
if((index=(int *)calloc(Nbt,sizeof(float)))==NULL) error_rec(24,"");
if((file=fopen(fitxer_resultats,"w"))==NULL) error_rec(31,fitxer_resultats);
iteracio=0;

/*calcul_indexs(index,Nbt,1);*/

for(i=0,j=0;i<Nbp;i++) for(k=0;k<Nang;k++,j++) *(index+j)=k*Nbp+i;

omple_mascara(masc,le,lcp);

while(iteracio<NITER) {
      iteracio++; 
      printf("\n*POCS iter: %d *",iteracio);
      rec_pocs(q,le,lcp,vp,W,index,masc,nc);
      output_ima(Nvox,1,q,file,tipout);
      }
fclose(file);
}
   
/*============================= REC_POCS ====================================*/

void rec_pocs(float *q2,float *le,float *lc,float *vp,float W,int *index,int *masc,int nc)
{
int i,j,k;
float factor,nc2;

for(k=0,i=index[k];k<Nbt;k++,i=index[k]){
   projeccio(lc,q2,i,&mp,Nbt);
   for(j= *(mp.ia+i); j< *(mp.ia+i+1); j++){
     if(masc[*(mp.ja+j)]){
       *(q2+ *(mp.ja+j))+=W*( *(le+i)- *(lc+i) )/ *(vp+i)* *(mp.ar+j);
       if(*(q2+*(mp.ja+j))<0) *(q2+*(mp.ja+j))=0;
       }
     }  
   }

projeccio(lc,q2,-1,&mp,Nbt);
numeroc(lc,&nc2,Nbt);
factor=nc/ nc2;
for(i=0;i<Nvox;i++) *(q2+i)*=factor;
}

