#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include <util/mle.h>
#include <util/inpout.h>
#include <util/crv.h>
#include <util/nrutil.h>
#define EOS '\0'
#define maxim(a,b) ((a)>=(b)?(a):(b))


/*============================= MLE =========================================*/

void mle(char *fitxer_resultats,float *le,float *vp,char tipoin[3],char tipout[3],float *n_contes,int NITER,int gravar,float W,char alg[4],char *CV,int inisor,int norma)
{  
int iteracio;
FILE *file,*file2,*file3,*file4,*file5;
float *q2,*qs,*qp,*lep,*les,*lcs,*lcp;
double acpl,adpl,acsl,adsl;
char fres_a[80],fres_b[80];
extern float *q;

q2=(float *)calloc(MM,4);
lcp=(float *)calloc(PN,4);
iteracio=0;

file5=fopen("factors.norm","w");

if(strncmp(CV,"no",2)==0){
   omple(q,MM,*n_contes/MM); 
   if((file=fopen(fitxer_resultats,"w"))==NULL) error_rec(31,fitxer_resultats);
   while(iteracio<NITER){
      iteracio++;
      printf("\n*MLE iter: %d *",iteracio);
      rec_mle(q,q2,le,lcp,vp,W,n_contes,alg,norma,iteracio,file5);
      if((iteracio%gravar)==0) output_ima(M,M,q,file,tipout);
      }
    /*  for(i=0;i<MM;i++) printf("\n %d  %f",i,q[i]);*/
      
   fclose(file);     
   }
else{
   qp=(float *)calloc(MM,4);
   qs=(float *)calloc(MM,4);
   lep=(float *)calloc(PN,4);
   les=(float *)calloc(PN,4);
   lcs=(float *)calloc(PN,4);
   acpl=adpl=acsl=adsl=0.;
   strcpy(fres_a,fitxer_resultats);
   strncat(fres_a,"a",1);
   strcpy(fres_b,fitxer_resultats);
   strncat(fres_b,"b",1);
   if((file2=fopen(fres_a,"w"))==NULL) error_rec(31,fres_a);
   if((file3=fopen(fres_b,"w"))==NULL) error_rec(31,fres_b);
   if((file4=fopen("cross_v.res","w"))==NULL) error_rec(31,"cross_v.res");
   proj_segones(le,lep,les,CV,inisor,tipoin,n_contes);  
   omple(qp,MM,*n_contes/MM);
   omple(qs,MM,*n_contes/MM);
   while(iteracio<NITER) {
      iteracio++;
      printf("\n*MLE-CV iter: %d *",iteracio);
      rec_mle(qp,q2,lep,lcp,vp,W,n_contes,alg,norma,iteracio,file5);
      rec_mle(qs,q2,les,lcs,vp,W,n_contes,alg,norma,iteracio,file5);
      /*for(i=0;i<MM;i++) *(q+i)=*(qs+i)+ *(qp+i);*/
      projeccio(lcp,qp,-1,&mp,PN);
      projeccio(lcs,qs,-1,&mp,PN);
      cros_val(lcp,lcs,lep,les,&acpl,&adpl,&acsl,&adsl,file4,iteracio);
      if((iteracio%gravar)==0){
        /* output_ima(M,M,q,file,tipout);*/
         output_ima(M,M,qp,file2,tipout);
         output_ima(M,M,qs,file3,tipout);
         }
      }
   fclose(file4); fclose(file3); fclose(file2); 
   free(qs); free(qp); free(lcs); free(lcp); free(les); free(lep);
   }   
fclose(file5);   
free(q2);

}

/*============================= REC_MLE =====================================*/

void rec_mle(float *qa,float *qb,float *le,float *lc,float *vp,float W,float *n_contes,char alg[4],int norma,int iteracio,FILE *file5)
{
extern sparse_i4 mp;
register int i,j;
float n_contes2,factor;


projeccio(lc,qa,-1,&mp,PN); 

for(j=0;j<PN;j++){
      if(*(lc+j)==0) {
  	/*printf("\n j=%d lc=%f le=%f ",j,lc[j],le[j]);*/
  
   	*(lc+j)=1.;/* 14.2.01 canviat, hi havia un 1e-10*/
   	}
   	
  /* if(le[j]<5.) le[j]=0.;
   printf("\n j=%d lc=%f le=%f ",j,lc[j],le[j]);*/	
   *(lc+j)= *(le+j)/ *(lc+j); 
 /*  printf("\n j=%d lc=%f le=%f ",j,lc[j],le[j]);*/
   }
per_mt(lc,qb,&mp,MM);
alg_mle(W,qa,qb,vp,alg);
if(norma!=0){
  printf("\t*N*");
  projeccio(lc,qa,-1,&mp,PN);
  numeroc(lc,&n_contes2,PN);
  factor=*n_contes/ n_contes2;
  fprintf(file5,"%d\t%f\t%f\n",iteracio,factor,pow(factor,(-1/W))-1);
  for(i=0;i<MM;i++){
       *(qa+i)*=factor;
       if(*(qa+i)<0) *(qa+i)=0;
       }   
  }
}  

/*============================= ALGORITME MLE ================================*/

void alg_mle(float W,float *q,float *q2,float *vp,char alg[4])
{

register int i;

if((strncmp(alg,"adi",3))==0) 
                       for (i=0;i<MM;i++) *(q+i) *=(1+W*( *(q2+i) / *(vp+i) -1));
else{
  if((strncmp(alg,"mult",4))==0){
      if(W==1.){
         for (i=0;i<MM;i++) if (*(vp+i) != 0)  *(q+i) *=( *(q2+i) / *(vp+i) );
          }
         
         else if(W==2.){
          for (i=0;i<MM;i++) if (*(vp+i) != 0)  *(q+i) *=( *(q2+i) / *(vp+i) ) * ( *(q2+i) / *(vp+i) );
          }
      else for (i=0;i<MM;i++) if(*(vp+i) != 0) *(q+i) *=pow(( *(q2+i) / *(vp+i) ),W);
      }                    
  else error_rec(33,alg);
  }
} 

/*============================= MLE-OS =========================================*/
void mleos(char *fitxer_resultats,float *le,char tipoin[3],char tipout[3],float *n_contes,int NITER,int gravar,float W,char alg[4],char *CV,int inisor,int OS,int norma)
{  
int i,iteracio,NdOS,*indexs;
FILE *file,*file2,*file3,*file4;
/*float *q2,*qs,*qp,*lep,*les,*lcs,*lcp,csl,cpl,dsl,dpl,as,ap,**vp_os,itepr,**lcr,segons; */
float *q2,*qs,*qp,*lep,*les,*lcs,*lcp,**vp_os,itepr,*lcr;  /* 28/02/01 *lcr pointer vector no matriu */
double acpl,adpl,acsl,adsl;
char fres_a[80],fres_b[80];
extern float *q;

if((PN%OS)!=0) error_rec(41," ");
NdOS=N/OS;
itepr=PN/OS;
q2=(float *)malloc(MM*4);
lcp=(float *)calloc(PN,4);
lcr=(float *)calloc(P,4);
vp_os=matrix(0,OS-1,0,MM-1);
indexs=(int*)calloc(N,4);
calcul_indexs(indexs,OS,NdOS); 
vector_pesos_os_reord(vp_os,OS,&mp,MM,N,indexs,P);
for(i=0;i<N;i++) indexs[i]*=P;

if((file=fopen(fitxer_resultats,"w"))==NULL) error_rec(31,fitxer_resultats);
iteracio=0;

printf("\n\nOS::%d\n",OS);

if(strncmp(CV,"no",2)==0){   
   omple(q,MM,*n_contes/MM);
   while(iteracio<NITER){
      iteracio++;
      printf("\n*MLE_OS iter: %d *",iteracio);
   /*   rec_mle_os(q,q2,le,lcp,lcr,vp_os,W,n_contes,alg,NdOS,indexs,OS);     */
      rec_mle_os(q,q2,le,lcp,lcr,vp_os,W,n_contes,alg,NdOS,indexs,OS,file,tipout,norma); /*canviat 1/3/01 */
      if((iteracio%gravar)==0) output_ima(M,M,q,file,tipout);
      }
   }
else{
   qp=(float *)calloc(MM,4);
   qs=(float *)calloc(MM,4);
   lep=(float *)calloc(PN,4);
   les=(float *)calloc(PN,4);
   lcs=(float *)calloc(PN,4);
   acpl=adpl=acsl=adsl=0.;
   strcpy(fres_a,fitxer_resultats);
   strncat(fres_a,"a",1);
   strcpy(fres_b,fitxer_resultats);
   strncat(fres_b,"b",1);
   if((file2=fopen(fres_a,"w"))==NULL) error_rec(31,fres_a);
   if((file3=fopen(fres_b,"w"))==NULL) error_rec(31,fres_b);
   if((file4=fopen("cross_v.res","w"))==NULL) error_rec(31,"cross_v.res");
   proj_segones(le,lep,les,CV,inisor,tipoin,n_contes);  
   omple(qp,MM,*n_contes/MM);
   omple(qs,MM,*n_contes/MM);
   while(iteracio<NITER) {
      iteracio++;
      printf("\n*MLE_OS-CV iter: %d *",iteracio);
      rec_mle_os(qp,q2,lep,lcp,lcr,vp_os,W,n_contes,alg,NdOS,indexs,OS,file2,tipout,norma);
      rec_mle_os(qs,q2,les,lcs,lcr,vp_os,W,n_contes,alg,NdOS,indexs,OS,file3,tipout,norma);
      /*for(i=0;i<MM;i++) *(q+i)=*(qs+i)+ *(qp+i);*/
      projeccio(lcp,qp,-1,&mp,PN);
      projeccio(lcs,qs,-1,&mp,PN);
      cros_val(lcp,lcs,lep,les,&acpl,&adpl,&acsl,&adsl,file4,iteracio);
      if((iteracio%gravar)==0){
         /*output_ima(M,M,q,file,tipout);*/
         output_ima(M,M,qp,file2,tipout);
         output_ima(M,M,qs,file3,tipout);
         }
      }
   fclose(file4); fclose(file3); fclose(file2); 
   free(qs); free(qp); free(les); free(lep); free(lcs);
   }   
   
fclose(file);
free(q2); free(lcp); free(lcr); /*  free(lcp); 1/3/01/ aquest alliberament de memòria és la causa de l'error de segmentation fault  */

}

/*============================= REC_MLE_OS ===================================*/
void rec_mle_os(float *qa,float *qb,float *le,float *lcp,float *lcr,float **vp_os,float W,float *n_contes,char alg[4],int NdOS,int *indexs,int OS,FILE *file,char *tipout,int norma)
{
extern sparse_i4 mp;
int i,j,k,kNdOS;
float n_contes2[1],factor;

for(k=0,kNdOS=0;k<OS;k++,kNdOS+=NdOS){
   for(i=0;i<MM;i++) *(qb+i)=0;  
   for(i=0;i<NdOS;i++){
      projeccio_indexs(lcr,qa,indexs[i+kNdOS],indexs[i+kNdOS]+P,&mp);
      for(j=0;j<P;j++){
         if(*(lcr+j)==0) *(lcr+j)=1.; /* 14.2.01 canviat, hi havia un 1e-10*/
         *(lcr+j)= *(le+indexs[i+kNdOS]+j)/ *(lcr+j); 
         }
      per_mt_index_append(lcr,qb,indexs[i+kNdOS],indexs[i+kNdOS]+P,&mp);
      }
   alg_mle_os(W,qa,qb,vp_os[k],alg);
   
   if(norma!=0){
     printf("\t*N*");
     projeccio(lcp,qa,-1,&mp,PN);
     numeroc(lcp,n_contes2,PN);
     factor=*n_contes/ n_contes2[0];
     for(i=0;i<MM;i++){
        *(qa+i)*=factor;
        if(*(qa+i)<0) *(qa+i)=0;
        } 
     }  
  }
}

/*============================= ALGORITME MLE OS =============================*/

void alg_mle_os(float W,float *q,float *q2,float *vp,char alg[4])
{

register int i;

if((strncmp(alg,"adi",3))==0) 
           for (i=0;i<MM;i++) {if (*(vp+i) != 0) *(q+i) *=(1+W*( *(q2+i) / *(vp+i) -1));}
else{
  if((strncmp(alg,"mult",4))==0){
     if(W==1.){
         for (i=0;i<MM;i++) if(*(vp+i) != 0) *(q+i) *=( *(q2+i) / *(vp+i) );
         }
     else if(W==2.){
       for (i=0;i<MM;i++) if(*(vp+i) != 0) *(q+i) *=( *(q2+i) / *(vp+i) ) * ( *(q2+i) / *(vp+i) );
       }
     else for (i=0;i<MM;i++) if(*(vp+i) != 0) *(q+i) *=pow(( *(q2+i) / *(vp+i) ),W);
     }                                         
  else error_rec(33,alg);
  }
} 

/*================================== INDEXS ==================================*/

#define maxim(a,b) ((a)>=(b)?(a):(b))
#define minim(a,b) ((a)<=(b)?(a):(b))
#define abs(a) ((a)>=0?(a):(-a))

void calcul_indexs(int *indexs,int OS,int NdOS)
{
int i,j,k,im,m,n,*ple,*iOS,*a,*sa,*dif;
    
iOS=(int *)calloc(OS,sizeof(int));
ple=(int *)calloc(OS,sizeof(int));
a=(int *)calloc(OS,sizeof(int));
sa=(int *)calloc(OS,sizeof(int));
dif=(int *)calloc(OS+1,sizeof(int));

for(i=0;i<OS;i++) iOS[i]=ple[i]=0;
ple[0]=1;

calcul_dif(dif,OS);

im=0;
for(k=1;k<OS;k++){
  for(i=1;i<OS;i++){
    if(!ple[i]){
      j=i-im;
      a[i]=dif[abs(j)];
      }
    }  
  for(i=0;i<OS;i++){
      a[i]*=(1-ple[i]);
      sa[i]*=(1-ple[i]);
      sa[i]+=a[i];
      }
  m=n=0;
  for(i=0;i<OS;i++) m=maxim(m,sa[i]);
  for(i=1;i<OS;i++)if(!ple[i]) m=minim(m,sa[i]);
  for(i=1;i<OS;i++) if(sa[i]==m) n=maxim(n,a[i]);
  for(i=OS-1;i>0;i--) if(sa[i]==m) if(a[i]<=n){
       n=a[i];
       im=i;
       }
  iOS[k]=im;
  ple[im]=1;  
  }

for(i=0;i<OS;i++) for(j=0;j<NdOS;j++) indexs[j+i*NdOS]=iOS[i]+OS*j;

free(a);free(sa);free(dif);free(ple);
}

/*==================== CALCUL_DIF ==============================*/

void calcul_dif(int *dif,int OS)
{
int i,OS2;
OS2=OS*OS;

for(i=1;i<=OS;i++) dif[i]+=2*(i*(i-OS))+OS2;
}

/*=========================================== ERROR =========================*/

    
 void error_rec(int nerr,char *text)
{
printf("\n\nERROR EN EL PROGRAMA rec");
switch(nerr){
    case 11: printf("\n\nError al obrir el fitxer de dades");
             printf("\nEl fitxer  %s  no es pot obrir\n",text); break;
    case 12: printf("\n\nError al token: falta el parametre: %s\n",text); break; 
    case 21: printf("\n\nError al seleccionar la reconstruccio"); 
             printf("\nLa reconstruccio  %s  no existeix o no esta implementada\n",text); break;
    case 23: printf("\n\nProblemes en malloc a main"); break;
    case 31: printf("\n\nError al obrir el fitxer de resultats");
             printf("\nEl fitxer  %s  no es pot obrir\n",text); break;
    case 33: printf("\n\nError al seleccionar tipus algorisme mle");
             printf("\nL'algorisme %s no existeix o no esta implementat (nomes adi/mult)\n",text);break;
    case 41: printf("\n\nError: el numero de subsets ha de ser congruent amb 60");
    default: printf("\n\nError del numero d'error en la funcio error_rec()"); 
    }
    
exit(0);
}  
