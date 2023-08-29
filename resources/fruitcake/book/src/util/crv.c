#include <string.h>
#include <stdlib.h>
#include <util/crv.h>

/*============================== SEGONES_PROJ ================================*/

void proj_segones(float *le,float *lep,float *les,char CV[80],int inisor,char tipoin[3],float nc[1])
{
int i;
nc[0]/=2.;
if((strncmp(CV,"pp",2))==0){ 
   rifa1(le,lep,inisor);
   rifa1(le,les,inisor+100);
   }
else if((strncmp(CV,"ps",2))==0) rifa(le,lep,les,inisor);
else if((strncmp(CV,"ud",2))==0) rifa_nou(le,lep,les,inisor);
else{
   for(i=0;i<PN;i++) lep[i]=le[i];
   input_fitxer(P,N,1,les,CV,tipoin);
   nc[0]*=2.;
   }   
}

/*============================= CROSS-VAL ====================================*/

void cros_val(float *lcp,float *lcs,float *lep,float *les,double *acpl,double *adpl,double *acsl,double *adsl,FILE *file2,int iteracio)
{
double dpl,cpl,dsl,csl,ap,as,loglepj,loglesj,leploglcpj,lesloglcpj,leploglcsj,lesloglcsj,lcpj,lcsj,logfact();
int j;

dpl=cpl=dsl=csl=0.;
loglepj=loglesj=leploglcpj=lesloglcpj=leploglcsj=lesloglcsj=lcpj=lcsj=0;

for(j=0;j<PN;j++){
      if(*(lcp+j)==0) *(lcp+j)=1e-20;
      if(*(lcs+j)==0) *(lcs+j)=1e-20;
      ap=log (lcp[j]);
      as=log (lcs[j]);
      loglepj+=logfact(lep[j]);
      loglesj+=logfact(les[j]);
      leploglcpj+=ap*(double)lep[j];
      lesloglcpj+=ap*(double)les[j];
      leploglcsj+=as*(double)lep[j];
      lesloglcsj+=as*(double)les[j];
      lcpj+=(double)(lcp[j]);
      lcsj+=(double)(lcs[j]);
      }
dpl+=(double)(-lcpj+leploglcpj-loglepj );
cpl+=(double)(-lcpj+lesloglcpj-loglesj );
dsl+=(double)(-lcsj+lesloglcsj-loglesj );
csl+=(double)(-lcsj+leploglcsj-loglepj );
   
ap=(cpl- *acpl)/(dpl- *adpl);
as=(csl- *acsl)/(dsl- *adsl);
printf("\tCVR1:%7.4f    CVR2:%7.4f",ap,as);

/*fprintf(file2,"%3.0f\t%f\t%f\n",(float)iteracio,ap,as);*/

fprintf(file2,"%3.0f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",(float)iteracio,ap,as,cpl,dpl,csl,dsl,lcpj,lcsj);

/*fprintf(file,"%3.0f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n",(float)iteracio,cpl,csl,dpl,dsl);

 if(as<=0 || ap<=0) {
      output_ima(M,M,q,file,tipout);
      exit(0);
      } */
   
*acpl=cpl;
*adpl=dpl;
*acsl=csl;
*adsl=dsl;
}
/*=================================== RIFA_NOU ==================================*/

void rifa_nou(float *le,float *lp,float *ls,int i0)
{
int vpix,i,vmax;
long idum;
float r;
float ran2();

idum=-i0;
r=ran2(&idum);
for(i=0;i<PN;i++){
   vmax=(int)floor(le[i]);
   for(vpix=0;vpix<vmax;vpix++){
   	r=ran2(&idum);
   	if( r < 0.5 ) lp[i]++;
   	else ls[i]++;
        }
   }

}
/*=================================== RIFA ==================================*/

void rifa(float *le,float *lp,float *ls,int i0)
{
float ran1();
int vpix,i,vmax,a;

srand(i0);
for(i=0;i<PN;i++){
   vmax=(int)floor(le[i]);
   for(vpix=0;vpix<vmax;vpix++){
   	a=rand();
   	if( (a%2 ) ==0 ) lp[i]++;
   	else ls[i]++;
        }
   }
}

/*=================================== RIFA1 =================================*/

void rifa1(float *le,float *lp,int i0)
{

float ran1();
int vpix,i,vmax,a;

for(i=0;i<PN;i++){
   vmax=(int)floor(le[i]);
   srand(vmax+i+i0);
   for(vpix=0;vpix<vmax;vpix++){
   	a=rand();
   	if( (a%2 ) ==0 ) lp[i]++;
        }
   }
}
  
/*==================== LOGARITME DEL FACTORIAL (aprox. Stirling) ===========*/

double logfact(double num)
{

double res,a;
int i,inum;

res=1.;
inum=(int) num;
if(inum==0) res=0.;
else{
    if(inum<=10 && inum>0){
	 for(i=1;i<=inum;i++) res *= (double)i;
	 res=log(res);
         }
    else {
    		a=log(num);
		res=0.5*log(2*M_PI)+(0.5*a)+num*a-num;
		}
    }
    
return(res);
}

/*============================================= RAN2 =========================*/

#include <math.h>

#define MO 714025
#define IA 1366
#define IC 150889

float ran2(long *idum)
{
	static long iy,ir[98];
	static int iff=0;
	int j;
	void nrerror();

	if (*idum < 0 || iff == 0) {
		iff=1;
		if ((*idum=(IC-(*idum)) % MO) < 0) *idum = -(*idum);
		for (j=1;j<=97;j++) {
			*idum=(IA*(*idum)+IC) % MO;
			ir[j]=(*idum);
		}
		*idum=(IA*(*idum)+IC) % MO;
		iy=(*idum);
	}
	j=1 + 97.0*iy/MO;
	if (j > 97 || j < 1) nrerror("RAN2: This cannot happen.");
	iy=ir[j];
	*idum=(IA*(*idum)+IC) % MO;
	ir[j]=(*idum);
	return (float) iy/MO;
}

#undef MO
#undef IA
#undef IC

