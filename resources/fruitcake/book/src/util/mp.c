#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <util/mp.h>

/*============================= PROJECCIO ====================================*/

void projeccio_old(float *lc,float *q,int p,sparse *mp,int nbins)
{
register int i,j;

if(p==-1)i=0;
    else i=p;

do{
  *(lc+i)=0;
  for(j=*(mp->ia+i);j<*(mp->ia+i+1);j++) *(lc+i)+= *(q+*(mp->ja+j))* *(mp->ar+j);
  i++;
  }
while(p==-1 && i<nbins);
}

/*============================= PROJECCIO_INDEX====================================*/

void projeccio_indexs_old(float *lcr,float *q,int i0,int i1,sparse *mp)
{
register int i,j,k;

for(i=i0,k=0;i<i1;i++,k++){
  *(lcr+k)=0;
  for(j=*(mp->ia+i);j<*(mp->ia+i+1);j++) *(lcr+k)+= *(q+*(mp->ja+j))* *(mp->ar+j);
  }
}

/*============================= PROJECCIO_SEQ ================================*/

void projeccio_seq_old(float *seqima,float *seqlc,sparse *mp,int MM,int nima,int nbins)
{
register int i,j,k,kPN,kMM;

for(k=0,kMM=0,kPN=0;k<nima;k++,kMM+=MM,kPN+=nbins){
  i=0;
  do{
    *(seqlc+kPN+i)=0;
    for(j=*(mp->ia+i);j<*(mp->ia+i+1);j++) *(seqlc+kPN+i)+= *(seqima+kMM+*(mp->ja+j))* *(mp->ar+j);
    i++;
    }
  while(i<nbins);
  }
}
/*============================= PROJECCIO_SEQ ================================*/

void projeccio_seq2_old(float *seqima,float *seqlc,sparse *mp,int MM,int nima,int nbins)
{
register int i,j, j1,k;

for(k=0;k<nima;k++,seqima+=MM){
  i=0;
  do{
    *seqlc=0;
    j1=*(mp->ia+i+1);
    for(j=*(mp->ia+i);j<j1;j++) *seqlc += *(seqima+ *(mp->ja+j))* *(mp->ar+j);
    ++seqlc;
    i++;
    }
  while(i<nbins);
  }
}

/*============================= PROJECCIO_SEQ ================================*/

void projeccio_seq3_old(float *seqima,float *seqlc,sparse *mp,int MM,int nima,int nbins)
{
register int i,j,j1,k;
float *punt_ar;
unsigned short int *punt_ja;
int *punt_ia;

for(k=0;k<nima;k++,seqima+=MM){
  i=0;
  punt_ar=mp->ar;
  punt_ja=mp->ja;
  punt_ia=mp->ia;
  do{
    *seqlc=0;
    j1=*(punt_ia+1);
    for(j=*(punt_ia);j<j1;j++,punt_ja++,punt_ar++) *seqlc += *(seqima+ *punt_ja)* *punt_ar;
    ++seqlc;
    ++i;
    ++punt_ia;
    } while(i<nbins);
  }
}

/*=========================== MULTIPLICAR PER TRANSPOSTA =====================*/

void per_mt_old(float *lc,float *q2,sparse *mp,int npixels)
{
register int i,fila;

for(i=0;i<npixels;i++) *(q2+i)=0;
fila=0;

for(i=0;i<mp->ne;i++){
   if(i>= *(mp->ia+fila+1)) fila++;
   *(q2+ *(mp->ja+i))+= *(mp->ar+i)* *(lc+fila);
   }
}

/*=========================== MULTIPLICAR PER TRANSPOSTA INDEX ===============*/

void per_mt_index_old(float *lc,float *q2,int i0,int i1,sparse *mp,int npixels)
{
register int i,l,index0,index1;

for(i=0;i<npixels;i++) *(q2+i)=0;
l=0;
index0=(int)(*(mp->ia+i0));
index1=(int)(*(mp->ia+i1));
for(i=index0;i<index1;i++){
   if(i>= *(mp->ia+i0+l+1)) l++; 
   *(q2+ *(mp->ja+i))+= *(mp->ar+i)* *(lc+l);
   }
}

/*=========================== MULTIPLICAR PER TRANSPOSTA INDEX APPEND ========*/

void per_mt_index_append_old(float *lc,float *q2,int i0,int i1,sparse *mp)
{
register int i,l,index0,index1;

l=0;
index0=(int)(*(mp->ia+i0));
index1=(int)(*(mp->ia+i1));
for(i=index0;i<index1;i++){
   if(i>= *(mp->ia+i0+l+1)) l++; 
   *(q2+ *(mp->ja+i))+= *(mp->ar+i)* *(lc+l);
   }
}

/*===========================  P E S O S  ====================================*/

void pesos_old(sparse *mp,char *nom_fitxer,int *P,int *N)
{
FILE *mat;
register int i,nitems;
 
if((mat = fopen(nom_fitxer,"r"))==NULL) error_mp_old(11,nom_fitxer);

printf("Llegint la matriu de pesos ...");

fread (&mp->n_fil,sizeof(size_t),1,mat);
fread (&mp->n_col,sizeof(size_t),1,mat);
fread (&mp->ne,sizeof(size_t),1,mat);
nitems=mp->ne;
if(!(mp->ar=(float *)malloc(nitems*sizeof(float))))error_mp_old(12,"");
if(!(mp->ja=(unsigned short int*)malloc(nitems*sizeof(short int))))error_mp_old(12,"");
if(!(mp->ia=(int*)malloc((mp->n_fil+1)*sizeof(int))))error_mp_old(12,"");

fread((float*)mp->ar,sizeof(float),nitems,mat);
if(!*P && !*N) error_mp_old(13,"");
if(!*N) *N=mp->n_fil/ *P;
fread ((short int*)mp->ja,sizeof(short int),nitems,mat);
fread((int *)mp->ia,sizeof(int),(size_t)(mp->n_fil+1),mat);
printf("   fet\n");

for(i=0;i<nitems;i++) *(mp->ar+i)/=*N;
fclose (mat);  
}

/*=============================== FREE_MPES ==================================*/

void free_mpes_old(sparse *f)
{
   free((float*)f->ia);
   free((short int*)f->ja);
   free((int*)f->ar);
}

/*====================================== VECTOR PESOS ========================*/

void vector_pesos_old(float *vp,char *tprec,sparse *mp,int nbins)
{
register int i,j;

if( strncmp(tprec,"pocs",4)==0) for(i=0;i<nbins;i++) 
  for(j=*(mp->ia+i);j<*(mp->ia+i+1);j++) *(vp+i)+=*(mp->ar+j)* *(mp->ar+j);

if(strncmp(tprec,"mle",3)==0) for(i=0;i<mp->ne;i++)
                                           *(vp+ *(mp->ja+i)) += *(mp->ar+i);

}

/*====================================== VECTOR PESOS OS ORD =================*/

void vector_pesos_os_ord_old(float **vp_os,int OS,int itepr,sparse *mp,int npixels)
{
register int i,j,k,fila;

for(i=0;i<OS;i++) for(j=0;j<npixels;j++) vp_os[i][j]=0;

for(i=0,k=0,fila=0;i<mp->ne;i++){
   if(i>= *(mp->ia+fila+1)){
      fila++;
      k=(fila/itepr);
      }
   vp_os[k][*(mp->ja+i)] += *(mp->ar+i);
   }
}

/*====================================== VECTOR PESOS OS REORD ===============*/

void vector_pesos_os_reord_old(float **vp_os,int OS,sparse *mp,int npixels,int nproj,int *ind,int P)
{
register int i,j,k,fila;

for(i=0;i<OS;i++) for(j=0;j<npixels;j++) vp_os[i][j]=0;

for(i=0,fila=0,k=0;i<mp->ne;i++){
   if(i>= *(mp->ia+fila+1)){
      fila++;
      k=fila/P;
      }
   for(j=0;j<nproj;j++) if(k==*(ind+j)){
      vp_os[j*OS/nproj][*(mp->ja+i)] += *(mp->ar+i);
      break;
      }
   }
}

/*=========================================== ERROR =========================*/

void error_mp_old(int nerr,char *text)
{
switch(nerr){
    case 11: printf("\n\nError al obrir fitxer de la matriu de pesos");
             printf("\nLa matriu de pesos  %s  no existeix\n",text); break;
    case 12: printf("\n\nError al reservar memoria per a mp"); break;
    case 13: printf("\n\nError al llegir mpes, P i N son zero"); break;
    default: printf("\n\nError del numero d'error en la funcio error_mp()"); 
    }
    
exit(0);
}    

