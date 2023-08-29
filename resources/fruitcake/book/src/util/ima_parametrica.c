#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <util/inpout.h>
#include <util/ima_parametrica.h>


/*calculo de las diferencias de dos imagenes realineadas*/
/*adaptacion del imafun del regbrain; en ester caso no se aplica ninguna mascara*/

/************************calculo de imagen diferencia sin empleo de una mascara***********************************************/

int diferen(struct imagen *im1,struct imagen *im2,struct imagen *it,int tipo)
{
	
	int N,M,P,MM,MMM;
	int i,j,k,l,i1,j1,k1,ll,imin,imax;
	float s1,s2;

	
	if (im1->nfil != im2->nfil || im1->ncol != im2->ncol || im1->nima != im2->nima)
	{
	printf("entra");
		return -1;
	}

N=im1->nfil;
M=im1->ncol;
P=im1->nima;
MM=N*M;
MMM=N*M*P;

imin=(1-tipo)/2;
imax=1-imin;

/*tipo hace referencia al número de pixeles que se utilizan para hacer la resta*/
for (l=0;l<MMM;l++) it->datos[l]=0.;

for (i=tipo/2;i<N-tipo/2;i++) {
    for (j=tipo/2;j<M-tipo/2;j++) {
        for (k=tipo/2;k<P-tipo/2;k++) {
	    s1=0.;s2=0.;
	    l=i+j*M+k*MM;
	
	    if (im2->datos[l]*im1->datos[l]>0.) {
	 
	       for (i1=imin;i1<imax;i1++) {
		   for (j1=imin;j1<imax;j1++) {
		       for (k1=imin;k1<imax;k1++) {
			   ll=i+i1+(j+j1)*M+(k+k1)*MM;
			
			   if (im2->datos[ll]*im1->datos[ll]>0.) {
			
			      s1+=im1->datos[ll];s2+=im2->datos[ll];
			   }
		       }
		   }
	       }
	        /*imagen de diferencias positivas y negativas*/
		  it->datos[l]=100.*(s2-s1)/s1+100.;
	          if (it->datos[l]>200.) it->datos[l]=200.;
	          if (it->datos[l]<0.) it->datos[l]=0.;
	    }
	                else it->datos[l]=0.;
}
}

}
return (0);
}

/*====================================== màscara imatge basal per aplicar a la imafun========================*/
int mask_imafun(struct imagen *dif,struct imagen *bas,struct imagen *ict)
{

int P,N,M,MM,MMM;
int i,index=0;
double maxb=0.,maxa=0.,thresh=.4,vmigb=0.,vmiga=0.; /*4.7.01 nivell per on tallar el màxim per fer màscara binària*imafun*/
 
N=bas->nfil;
M=bas->ncol;
P=bas->nima;

MM=N*M;
MMM=N*M*P;

for(i=0;i<MMM;i++) {
	if(bas->datos[i]>0.1) {
		vmigb+=bas->datos[i];
		vmiga+=ict->datos[i];
		index++;
	}
}

vmigb/=index;vmiga/=index;
printf("\n vmigb=%f vmiga=%f\n",vmigb,vmiga);
for(i=0;i<MMM;i++) {
	if(bas->datos[i]>maxb) maxb=bas->datos[i];		/*buscar màxim estudi basal*/
	if(ict->datos[i]>maxa) maxa=ict->datos[i];
}

printf("\n maxb=%f  thres=%f maxa=%f\n",maxb,maxb*thresh,maxa);
for(i=0;i<MMM;i++) {
	if(bas->datos[i]<thresh*maxb) dif->datos[i]=0.;
}
return(0);
}

/*======================================separamos diferencias positivas y negativas==============================================*/
int separa_pos_neg(struct imagen *it,struct imagen *it1,struct imagen *it2)
{

int N,M,P,MMM;
int l;
	
N=it->nfil;
M=it->ncol;
P=it->nima;

MMM=N*M*P;

for (l=0;l<MMM;l++) it1->datos[l]=0.;
for (l=0;l<MMM;l++) it2->datos[l]=0.;

for (l=0;l<MMM;l++){
  if (it->datos[l]!=0.){
       if(it->datos[l]>0. && it->datos[l]<100.) it1->datos[l]=it->datos[l];
       if(it->datos[l]>=100. && it->datos[l]<=200.) it2->datos[l]=it->datos[l];
                       }
 else if (it->datos[l]==0.){
      it1->datos[l]=it->datos[l];
      it2->datos[l]=it->datos[l];
                     }
}
return(0);
}
