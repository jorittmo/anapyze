#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <util/mp_i4.h>

/*============================= PROJECCIO ====================================*/

void projeccio(float *lc,float *q,int p,sparse_i4 *mp,int nbins)
{
register int i,j;

if(p==-1)i=0;
    else i=p;

do{
  lc[i]=0;
  i++;
  }
while(p==-1 && i<nbins);
    

if(p==-1)i=0;
    else i=p;    
    
do{
  //for(j=*(mp->ia+i);j<*(mp->ia+i+1);j++) *(lc+i)+= *(q+*(mp->ja+j))* *(mp->ar+j);
  for(j=mp->ia[i];j<mp->ia[i+1];j++) 
  {
/*	  k=mp->ja[j]/Npix;
	  l=i/Nbins;
	  k2=l/Nrings;
	  k1=l%Nrings;
	  ind_imagen2D=mp->ja[j]%Npix;
	  		  */
	  	  
	  lc[i]+= q[mp->ja[j]]*mp->ar[j];
/*	  i'=
	  j'=ind_imagen2D+(Ncortes-k-1)*Npix;
	  lc[i']+=q[j']*mp->ar[j];
  */}
  i++;
  }
while(p==-1 && i<nbins);
}

/*============================= PROJECCIO_INDEX============================== */

void projeccio_indexs(lcr,q,i0,i1,mp)
float *lcr,*q;
int i0,i1;
sparse_i4 *mp;
{
register int i,j,k;

for(i=i0,k=0;i<i1;i++,k++){
  *(lcr+k)=0;
  for(j=*(mp->ia+i);j<*(mp->ia+i+1);j++) *(lcr+k)+= *(q+*(mp->ja+j))* *(mp->ar+j);
  }
}

/*=============================== PROJECCIO_SEQ ================================*/

void projeccio_seq(float *seqima,float *seqlc,sparse_i4 *mp,int MM,int nima,int nbins)
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

/*=========================== MULTIPLICAR PER TRANSPOSTA =====================*/

void per_mt(float *lc,float *q2,sparse_i4 *mp,int npixels)
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

void per_mt_index(float *lc,float *q2,int i0,int i1,sparse_i4 *mp,int npixels)
{
register int i,l,index0,index1;

for(i=0;i<npixels;i++) *(q2+i)=0;
l=0;
index0=(*(mp->ia+i0));
index1=(*(mp->ia+i1));
for(i=index0;i<index1;i++){
   if(i>= *(mp->ia+i0+l+1)) l++; 
   *(q2+ *(mp->ja+i))+= *(mp->ar+i)* *(lc+l);
   }
}

/*=========================== MULTIPLICAR PER TRANSPOSTA INDEX APPEND ========*/

void per_mt_index_append(float *lc,float *q2,int i0,int i1,sparse_i4 *mp)
{
register int i,l,index0,index1;

l=0;
index0=(*(mp->ia+i0));
index1=(*(mp->ia+i1));
for(i=index0;i<index1;i++){
   if(i>= *(mp->ia+i0+l+1)) l++; 
   *(q2+ *(mp->ja+i))+= *(mp->ar+i)* *(lc+l);
   }
}

/*===========================  P E S O S  ====================================*/

void pesos(sparse_i4 *mp,char *nom_fitxer,int *P,int *N)
{
FILE *mat;
register int nitems;
 int i;
 
if((mat = fopen(nom_fitxer,"rb"))==NULL) error_mp(11,nom_fitxer);

printf("Llegint %s  ...",nom_fitxer);

fread (&mp->n_fil,sizeof(size_t),1,mat);
fread (&mp->n_col,sizeof(size_t),1,mat);
fread (&mp->ne,sizeof(size_t),1,mat);
nitems=mp->ne;
if(!(mp->ar=(float *)malloc(nitems*sizeof(float))))error_mp(12,"");
if(!(mp->ja=(int*)malloc(nitems*sizeof(int))))error_mp(12,"");
if(!(mp->ia=(int*)malloc((mp->n_fil+1)*sizeof(int))))error_mp(12,"");

fread((float*)mp->ar,sizeof(float),nitems,mat);
if(!*P && !*N) error_mp(13,"");
if(!*N) *N=mp->n_fil/ *P;
fread ((int*)mp->ja,sizeof(int),nitems,mat);
fread((int *)mp->ia,sizeof(int),(size_t)(mp->n_fil+1),mat);
printf("   fet\n");

for(i=0;i<nitems;i++) *(mp->ar+i)/=*N;
fclose (mat);  
}



/*===========================  PESOSxOS ====================================*/

void pesos_xOS(sparse_i4 *mp,char *nom_fitxer)
{
FILE *mat;
register int nitems;
 
if((mat = fopen(nom_fitxer,"r"))==NULL) error_mp(11,nom_fitxer);

printf("Llegint  %s ...",nom_fitxer);

fread (&mp->n_fil,sizeof(size_t),1,mat);
fread (&mp->n_col,sizeof(size_t),1,mat);
fread (&mp->ne,sizeof(size_t),1,mat);
nitems=mp->ne;

fread((float*)mp->ar,sizeof(float),nitems,mat);
fread ((int*)mp->ja,sizeof(int),nitems,mat);
fread((int *)mp->ia,sizeof(int),(size_t)(mp->n_fil+1),mat);
printf("   fet\n");

/*for(i=0;i<nitems;i++) *(mp->ar+i)/=*N;*/
fclose (mat);  
}


/*======================  ESCRIU MATRIU ======================= */


void escriu_mpes(char *nom_fitxer,int *PN,int *MM,sparse_i4 *mp)
{
FILE *mat;

if((mat=fopen(nom_fitxer,"w"))==NULL) error_mp(30,nom_fitxer);
fwrite (PN,sizeof(size_t),1,mat);
fwrite (MM,sizeof(size_t),1,mat);
fwrite (&(mp->ne),sizeof(size_t),1,mat);
fwrite (mp->ar,sizeof(float),(size_t)mp->ne,mat);
fwrite (mp->ja,sizeof(int),(size_t)mp->ne,mat);
fwrite (mp->ia,sizeof(int),(size_t)(*PN+1),mat);
fclose (mat);
}

/*=============================== FREE_MPES ==================================*/

void free_mpes(sparse_i4 *f)
{
   free(f->ar);
   free(f->ja);
   free(f->ia);
}

/*====================================== VECTOR PESOS ========================*/

void vector_pesos(float *vp,char *tprec,sparse_i4 *mp,int nbins)
{
register int i,j;

if( strncmp(tprec,"pocs",4)==0) for(i=0;i<nbins;i++) 
  for(j=*(mp->ia+i);j<*(mp->ia+i+1);j++) *(vp+i)+=*(mp->ar+j)* *(mp->ar+j);

if(strncmp(tprec,"mle",3)==0) for(i=0;i<mp->ne;i++)
                                           *(vp+ *(mp->ja+i)) += *(mp->ar+i);

}

/*====================================== VECTOR PESOSxOS =====================*/

void vector_pesosxOS(float *vp,sparse_i4 *mp,int nvox)
{
register int i;

// inicializa vector de pesos 
for(i=0;i<nvox; i++) *(vp+i)=0;

// para cada elemento de la imagen suma la contribución de todos los elementos de 
// las proyecciones, se obtiene la llamada matriz de sensibilidad (dimensiones de imagen)
for(i=0;i<mp->ne;i++) *(vp+ *(mp->ja+i)) += *(mp->ar+i);
}

/*====================================== VECTOR PESOS OS ORD =================*/

void vector_pesos_os_ord(float **vp_os,int OS,int itepr,sparse_i4 *mp,int npixels)
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

void vector_pesos_os_reord(float **vp_os,int OS,sparse_i4 *mp,int npixels,int nproj,int *ind,int P)
{
register int i,j,k,fila;

for(i=0;i<OS;i++) for(j=0;j<npixels;j++) vp_os[i][j]=0;

for(i=0,fila=0,k=0;i<mp->ne;i++){
   if(i>= *(mp->ia+fila+1)){
      fila++;
      k=fila/P;
      }
   for(j=0;j<nproj;j++) if(k==*(ind+j)){
      vp_os[(j*OS)/nproj][*(mp->ja+i)] += *(mp->ar+i);
      break;
      }
   }
}


/*====================================== Obtiene imagen sensibilidad =====================*/

void ima_sensibilidad_spect(float *sens,sparse_i4 *mp)
{
int i,nvox;


nvox=mp->n_fil*mp->n_col*mp->Ntallsima;

// inicializa vector  
for(i=0;i<nvox;i++) sens[i]=0.;

// se recorre el vector mp.ja que contiene posiciones de la imagen con contribución no nula y se va sumando
// su contribución (contenida en la posición correspondiente de mp.ar)
for(i=0;i<mp->ne;i++) 
{
sens[mp->ja[i]]+=mp->ar[i];
//if(mp->ja[i]>=143360) printf("\n%d",mp->ja[i]);
}
printf("\n********************************************************************************\n");
}

/*=========================================== ERROR =========================*/

void error_mp(int nerr,char *text)
{
switch(nerr){
    case 11: printf("\n\nError al obrir fitxer de la matriu de pesos");
             printf("\nLa matriu de pesos  %s  no existeix\n",text); break;
    case 12: printf("\n\nError al reservar memoria per a mp"); break;
    case 13: printf("\n\nError al llegir mpes, P i N son zero"); break;
    case 30: printf("\n\nError al obrir el fitxer per escriure mpes\n%s  no es pot obreir\n",text); break;
    default: printf("\n\nError del numero d'error en la funcio error_mp()"); 
    }
    
exit(0);
}    
