#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <spect/mle3dxOS.h>
#include <spect/indexOS.h>
#include <util/inpout.h>

#define EPSILON 1E-2

/*============================= MLExOS inicializando al cilindro inscrito =========================================*/

void mlemxOS(q,Nfil,Ntalls,Nvox,le,lc,Nang,Nbp,NbOS,n_contes,OS,fmatriu_pes,NITER,norma,fitxer_ima,tipout,le_scat)

char *fmatriu_pes,*fitxer_ima,tipout[3];
int Nfil,Ntalls,Nvox,Nang,Nbp,NbOS,NITER,norma;
float *q,*le,*lc,*le_scat,*n_contes;

{
int iteracio;		/* index indica la iteracio actual */
int NangOS;		/* nombre d'angles de projeccio a cada subset */  
int Nitems_max;		/* nombre maxims d'elements a les matrius de pesos. Per alloc */                
float *q2;		/* imatge auxiliar */
float *vp;		/* vector pesos (suma per columnes de la matriu mp) */
sparse_i4 mp;		/* variable de la matriu de pesos */

FILE *file;

NangOS=Nang/OS;

/* .... check dimensions maximas de les matrius de pesos ... */
check_nitems_max_mp(fmatriu_pes,OS,&Nitems_max);

/* .... reserva de memoria per a la matriu de pesos ... */
if((mp.ar=(float *)calloc(Nitems_max,sizeof(float)))==NULL) error_mle3dxOS(23,"");  
if((mp.ja=(int *)calloc(Nitems_max,sizeof(int)))==NULL) error_mle3dxOS(23,"");
if((mp.ia=(int *)calloc(NbOS+1,sizeof(int)))==NULL) error_mle3dxOS(23,"");    

/* .... reserva de memoria de les variables auxiliars ... */
if((q2=(float *)calloc(Nvox,sizeof(float)))==NULL) error_mle3dxOS(24,"");  
if((vp=(float *)calloc(Nvox,sizeof(float)))==NULL) error_mle3dxOS(24,"");  

/* .... inicialitza imatge a promig de n_comptes (ULL: nom≈s a cilindre circumscrit) ... */
omple_cilindre(q,Nfil,Ntalls,*n_contes/Nvox);

/* .... obrir fitxer de resultats ... */
if((file=fopen(fitxer_ima,"w"))==NULL) error_mle3dxOS(31,fitxer_ima);

/* .... iteracions de la reconstruccio ... */
for(iteracio=0; iteracio<NITER; iteracio ++){
    
    printf("\n*MLE iteracio: %d *\n",iteracio);
     
    /* .... calcul de la reconstruccio ... */   
    rec_mlexOS(q,q2,le,lc,vp,n_contes,OS,Nvox,NbOS,NangOS,norma,fmatriu_pes,&mp,le_scat);
      
    /* .... escritura de la imatge obtinguda a la iteracio actual ... */
    output_ima(Nvox,1,q,file,tipout);
    }

/* .... tancar fitxer de resultats ... */
fclose(file);

free_mpes(&mp);
free(q2);  free(vp);
}


/*============================= MLExOS con mascara =========================================*/

void mlemxOS_mascara(q,Nfil,Ntalls,Nvox,le,lc,Nang,Nbp,NbOS,n_contes,OS,fmatriu_pes,NITER,norma,fitxer_ima,tipout,le_scat,num_pix_masc)

char *fmatriu_pes,*fitxer_ima,tipout[3];
int Nfil,Ntalls,Nvox,Nang,Nbp,NbOS,NITER,norma;
float *q,*le,*lc,*le_scat,*n_contes;
float num_pix_masc;

{
int iteracio;		/* index indica la iteracio actual */
int NangOS;		/* nombre d'angles de projeccio a cada subset */  
int Nitems_max;		/* nombre maxims d'elements a les matrius de pesos. Per alloc */                
int i;                  /*indice para recorrer la mascara*/
float *q2;		/* imatge auxiliar */
float *vp;		/* vector pesos (suma per columnes de la matriu mp) */
sparse_i4 mp;		/* variable de la matriu de pesos */

FILE *file;

NangOS=Nang/OS;

/* .... check dimensions maximas de les matrius de pesos ... */
check_nitems_max_mp(fmatriu_pes,OS,&Nitems_max);

/* .... reserva de memoria per a la matriu de pesos ... */
if((mp.ar=(float *)calloc(Nitems_max,sizeof(float)))==NULL) error_mle3dxOS(23,"");  
if((mp.ja=(int *)calloc(Nitems_max,sizeof(int)))==NULL) error_mle3dxOS(23,"");
if((mp.ia=(int *)calloc(NbOS+1,sizeof(int)))==NULL) error_mle3dxOS(23,"");    

/* .... reserva de memoria de les variables auxiliars ... */
if((q2=(float *)calloc(Nvox,sizeof(float)))==NULL) error_mle3dxOS(24,"");  
if((vp=(float *)calloc(Nvox,sizeof(float)))==NULL) error_mle3dxOS(24,"");  

/*imagen sobre la que se inicializa el proceso iterativo de la reconstruccion*/
for(i=0;i<Nvox;i++){
if(q[i]>EPSILON) q[i]=*n_contes/num_pix_masc;
}

/* .... obrir fitxer de resultats ... */
if((file=fopen(fitxer_ima,"w"))==NULL) error_mle3dxOS(31,fitxer_ima);

/* .... iteracions de la reconstruccio ... */
for(iteracio=0; iteracio<NITER; iteracio ++){
    
    printf("\n*MLE iteracio: %d *\n",iteracio);
     
    /* .... calcul de la reconstruccio ... */   
    rec_mlexOS(q,q2,le,lc,vp,n_contes,OS,Nvox,NbOS,NangOS,norma,fmatriu_pes,&mp,le_scat);
      
    /* .... escritura de la imatge obtinguda a la iteracio actual ... */
    output_ima(Nvox,1,q,file,tipout);
    }

/* .... tancar fitxer de resultats ... */
fclose(file);

free_mpes(&mp);
free(q2);  free(vp);
}

/*============================= REC_MLExOS ===================================*/

void rec_mlexOS(q,q2,le,lc,vp,n_contes,OS,Nvox,NbOS,NangOS,norma,fmatriu_pes,mp,le_scat)
float *q,*q2,*le,*lc,*vp,*n_contes,*le_scat;
int OS,NangOS,NbOS,Nvox,norma;
char *fmatriu_pes;
sparse_i4 *mp;
{ 
register int i,k;
int index_le;			/* index sobre les projeccions experimentals */
float n_contes2[1],factor; 		/* nom de les comptes calculades i factor de correcio */
char fmatriu_pesxOS[240];		/* nom de la matriu de pesos. S'actualitza a cada subset */

for(k=0;k<OS;k++){
   printf("\nOS::%d\t",k);
     
   /* ..... generacio del nom de la matriu de pesos .... */
   generar_nom_mpes(fmatriu_pesxOS,fmatriu_pes,k,OS);
     
   /* .... lectura de la matriu de pesos corresponent al subset ... */
   pesos_xOS(mp,fmatriu_pesxOS);
      
   /* .... calcul del vector de pesos de la matriu (suma de columnes) ... */
   vector_pesosxOS(vp,mp,Nvox);
      
   /* .... inicialitzacio de la matriu auxiliar ... */ 
   for(i=0;i<Nvox;i++) q2[i]=0;  
  
   /* .... calcul de la projeccio de la imatge actual ... */
   projeccio(lc,q,-1,mp,NbOS);

   /* .... quocient de projeccions ...*/     
   for(i=0,index_le=k*NbOS ; i<NbOS ; i++,index_le++){
       /* .... donar un valor minim a zeros per evitar dividir per 0 ... */
       if(lc[i]==0) lc[i]=1.; 
       lc[i]= le[index_le]/(lc[i]+le_scat[index_le]);   
       }

   /* .... producte per la matriu trasposta ... */
   per_mt(lc,q2,mp,Nvox);

   /* .... actualitzacio de la imatge ... */   
   for (i=0;i<Nvox;i++) if (vp[i]!= 0)  q[i] *= q2[i]/vp[i];
      
   /* .... normalitzar per fer que projeccio calculada i experimental tinguin el mateix n_comptes ... */ 
   if(norma){ 
      
      projeccio(lc,q,-1,mp,NbOS);     
      numeroc(lc,n_contes2,NbOS);
      factor=n_contes[k]/ n_contes2[0];
      for(i=0;i<Nvox;i++) q[i]*=factor;
      printf("*N* %f\n",factor);
      }
   }
}

/* ============================= NUMERO DE COMPTES PER SUBSET ==============*/

void numerocpxOS(le,n_comptes,OS,NbOS)
float *le,*n_comptes;
int OS,NbOS;
{
int i,j,k;
float nct=0.;

for(k=0;k<OS;k++) n_comptes[k]=0;

printf("\n\nComptes a les projeccions:\n");

for(k=0,j=0;k<OS;k++){
    for(i=0;i<NbOS;i++,j++) n_comptes[k]+=le[j];
    printf("OS: %d\t%8.2f\n",k, n_comptes[k]);
    nct+= n_comptes[k];
    }
printf("total comptes a les projeccions: %f\n",nct);
}

/* ============================= CHECK NITEMS MAX MPES ==============*/

void check_nitems_max_mp(fmatriu_pes,OS,Nitems_max)
char *fmatriu_pes;
int OS,*Nitems_max;
{
char fmatriu_pesxOS[240];
int j,aux;
FILE *mat;

*Nitems_max=0;

for(j=0;j<OS;j++){
   generar_nom_mpes(fmatriu_pesxOS,fmatriu_pes,j,OS);
   
   if((mat = fopen(fmatriu_pesxOS,"r"))==NULL) error_mp(11,fmatriu_pesxOS);
   
   fread (&aux,sizeof(size_t),1,mat);
   fread (&aux,sizeof(size_t),1,mat);
   fread (&aux,sizeof(size_t),1,mat);
   
   if(aux > *Nitems_max) *Nitems_max=aux;
   
   fclose(mat);
   }
}   

void error_mle3dxOS(int nerr,const char *text)
{
printf("\n\nERROR EN EL PROGRAMA rec3dxOS");

switch(nerr){
    case 11: printf("\n\nError al obrir el fitxer de dades");
             printf("\nEl fitxer  %s  no es pot obrir\n",text); break;
    case 12: printf("\n\nError al token: falta el parametre: %s\n",text); break; 
    case 23: printf("\n\nProblemes en malloc a matriu de pesos"); break;
    case 24: printf("\n\nProblemes en malloc a mlexOS"); break;
    case 31: printf("\n\nError al obrir el fitxer de resultats");
             printf("\nEl fitxer  %s  no es pot obrir\n",text); break;
    default: printf("\n\nError del numero d'error en la funcio error_mle3dxOS()"); 
    }
    
exit(0);
}    
