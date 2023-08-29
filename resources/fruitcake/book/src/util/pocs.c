#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <util/pocs.h>
#include <util/mle.h>
#include <util/inpout.h>
#define EOS '\0'
#define maxim(a,b) ((a)>=(b)?(a):(b))

/*==================================== POCS ==================================*/

void pocs(char *fitxer_resultats,float *le,float *vp,char tipout[3],int NITER,int primera,int gravar,float W,float nc)
{

int iteracio,*index;
FILE *file;
float *lcp;
char *masc;

lcp=(float *)calloc(PN,4);
masc=(char *)calloc(MM,1);

/*index=(int *)calloc(PN,4);*/
if((file=fopen(fitxer_resultats,"w"))==NULL) error_rec(31,fitxer_resultats);
iteracio=0;
/*calcul_indexs(index,PN,1);*/ /*Per actualització no sequencial sino respecte bins separats*/
omple_mascara(masc,le,lcp); /* per no calcular on les projeccions son zero */

while(iteracio<NITER) {
      printf("\n*POCS iter: %d *",iteracio);
      rec_pocs(le,lcp,vp,W,nc,index,masc);
      if( ((iteracio%gravar)==0) && (iteracio>=primera) ) output_ima(M,M,q,file,tipout);
      iteracio++;
      }

fclose(file);
}

/*============================= REC_POCS ====================================*/

void rec_pocs(float *le,float *lc,float *vp,float W,float nc,int *index,char *masc)
{
int i,j,iq;
float factor,nc2;

/*for(k=0,i=index[k];k<PN;k++,i=index[k]){*/ /*canviar els comentaris d'aquesta fila per la seguent si es vol fer a bins alternats*/
for(i=0;i<PN;i++){
   projeccio(lc,q,i,&mp,PN);
   for(j= *(mp.ia+i); j< *(mp.ia+i+1); j++){
     iq= *(mp.ja+j);
     if(masc[iq]){
       q[iq]+=W*( *(le+i)- *(lc+i) )/ *(vp+i)* *(mp.ar+j);
       if ( q[iq] <0) q[iq]=0.;
       }
     }
   }

projeccio(lc,q,-1,&mp,PN);
numeroc(lc,&nc2,PN);
factor=nc/ nc2;
for(i=0;i<MM;i++) q[i]*=factor;
}

/*============================= MASCARA ====================================*/

void omple_mascara(char *masc,float *le,float *lc)
{
register int i,j;

for(i=0;i<MM;i++) masc[i]=1;
for(i=0;i<PN;i++) if(*(le+i)==0) for(j=mp.ia[i];j<mp.ia[i+1];j++) *(masc+mp.ja[j])=0;
}
