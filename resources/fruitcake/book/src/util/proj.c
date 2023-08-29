/* Carles Falcon */

#include <util/proj.h>
#include <util/nrutil.h>
#include <util/inpout.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*=========================================SEPARA_TALL========================*/

void separa_tall(float *p,int n,float *le,int P,int N,int NT)
{
register int i,k,i1,i2;

for(k=0,i1=0,i2=n*P;k<N;k++,i2+=P*(NT-1)) 
	{
		
		for(i=0;i<P;i++,i1++,i2++) 
		{
			le[i1]=p[i2];
		}
	}
}

/*=========================================LOG DE LA PROJ========================*/

void log_projeccio(float *le1,float *le,int PN)
{
int i;
for(i=0;i<PN;i++) 	le[i]=-log(le1[i]);
}

/*======================================= PROJ-SCATTERING ====================*/

void scatering(float *le,int conv,int PN,int P,int N)
{

float *sinc,*vec;
int i,k,kP,P4,P2;
P4=P*4;
P2=P*2;

sinc=vector(0,P4-1);
vec=vector(0,P4-1);

for(i=0;i<P4;i++) sinc[i]=0;

for(i=2;i<P ;i+=2){
  sinc[i]=0.035*exp(-0.2*(float)(i/2));
  sinc[P4-i]=sinc[i];
  }

sinc[0]=0.035;
fft(sinc-1,P2,1);

for(i=2;i<P4;i+=2) sinc[i]=sqrt(sinc[i]*sinc[i]+sinc[i+1]*sinc[i+1]);
    
for(k=0;k<N;k++){
  kP=k*P;
  for(i=0;i<P4;i++) vec[i]=0;
  for(i=0;i<P;i++) vec[2*i]=le[kP+i];
  
  fft(vec-1,P2,1);
  
      if(conv==1){
        for(i=0;i<P4;i+=2){
            vec[i]*=sinc[i]/(P2);
            vec[i+1]*=sinc[i]/(P2);
            }
       fft(vec-1,P2,-1);
       for(i=0;i<P;i++) le[kP+i]+=vec[2*i];
       }
       
      else{
        for(i=0;i<P4;i+=2){
            vec[i]/=(sinc[i]+1)*(P2);
            vec[i+1]/=(sinc[i]+1)*(P2);
            }
        fft(vec-1,P2,-1);
        for(i=0;i<P;i++) le[kP+i]=vec[2*i];
        }
  }
  
free(sinc);
free(vec);
}

/*====================================== FILTRAR_VECT ========================*/

void filtrar_vect(char *fitxer_filtre,float *projeccio,int M,int N)
{
float *sinc,*vector,*h,s,freq_tall,expo,MM8c2,c,distancia,t_pixel;
int i,k,kM,M4,M2,colimador;
FILE *file;
char filtre[10];

M4=M*4;
M2=M*2;

sinc=(float *)malloc(M4*4);
h=(float *)malloc(M2*4);
vector=(float *)malloc(M4*4);

if((file=fopen(fitxer_filtre,"r"))==NULL) 
 printf("\n error en el fitxer de filtre\n\n"), exit(0);

fscanf(file,"%s",filtre); printf("\nFILTER: %s\n",filtre);
h[0]=1;

if(strcmp(filtre,"rampa")==0) 
{
	for(i=0;i<M2;i++) h[i]=1;
}
else
{
  if(strcmp(filtre,"shepp")==0){
            for(i=1;i<M+1;i++){
               h[i]=sin(M_PI*i/(float)M)/M_PI/(float)i*M;
               h[M2-i]=h[i];
               }
            }   
  else{
    if(strcmp(filtre,"butt")==0){ 
            fscanf(file,"%f",&freq_tall);
            fscanf(file,"%f",&expo);           
            for(i=1;i<M+1;i++){
               h[i]=1.*sqrt(1+pow((float)i/(2*(M-1)*freq_tall),2*expo));
               h[M2-i]=h[i];
               }
            }
    else{
      if(strcmp(filtre,"hann")==0){
            fscanf(file,"%f",&freq_tall);
            fscanf(file,"%f",&expo);
            for(i=1;i<M+1;i++){
               h[i]=pow(.5+.5*cos(M_PI*pow((float)i/(M-1),freq_tall)),expo);
               h[M2-i]=h[i];
               }
            }
       else{
         if(strcmp(filtre,"metz")==0){
            fscanf(file,"%f",&expo);    
            fscanf(file,"%f",&distancia); 
            fscanf(file,"%d",&colimador); 
            fscanf(file,"%f",&t_pixel);   
            c=sigma(distancia,colimador,t_pixel);
            MM8c2=M_PI*M_PI*c*c/M2/M;
            for(i=1;i<M+1;i++){
               s= exp(-(float)i*(float)i*MM8c2);
               h[i]=1/s*( 1-pow(1-s*s,expo) );
               h[M2-i]=h[i];
               }
            }
         else  printf("\n error filter type\n\n"), exit(0);;
         }     
       }       
    }
}
   
fclose(file);
     
for(i=0;i<M4;i++) sinc[i]=0;

for(i=2;i<M;i+=4){
  sinc[i]=-1./(i/2)/(i/2)/M_PI/M_PI;
  sinc[M4-i]=sinc[i];
  }

sinc[0]=.25;
fft(sinc-1,M2,1);

for(k=0;k<N;k++){
  kM=k*M;
  for(i=0;i<M4;i++) vector[i]=0;
  for(i=0;i<M;i++) vector[2*i]=projeccio[kM+i];
  fft(vector-1,M2,1);
  for(i=0;i<M4;i+=2){
      vector[i]*=sinc[i]/(M2)*h[i/2];
      vector[i+1]*=sinc[i]/(M2)*h[i/2];
      }
  fft(vector-1,M2,-1);
  for(i=0;i<M;i++) projeccio[kM+i]=vector[2*i];
  }
  
free(sinc); free(vector); free(h);
}









/////////////////////////7


void filtrar_ramp(float *projeccio,int M,int N)
{
float *sinc,*vector,*h;
int i,k,kM,M4,M2;

M4=M*4;
M2=M*2;

sinc=(float *)malloc(M4*4);
h=(float *)malloc(M2*4);
vector=(float *)malloc(M4*4);

h[0]=1;
for(i=0;i<M2;i++) h[i]=1;

for(i=0;i<M4;i++) sinc[i]=0;

for(i=2;i<M;i+=4){
  sinc[i]=-1./(i/2)/(i/2)/M_PI/M_PI;
  sinc[M4-i]=sinc[i];
  }

sinc[0]=.25;
fft(sinc-1,M2,1);

for(k=0;k<N;k++){
  kM=k*M;
  for(i=0;i<M4;i++) vector[i]=0;
  for(i=0;i<M;i++) vector[2*i]=projeccio[kM+i];
  fft(vector-1,M2,1);
  for(i=0;i<M4;i+=2){
      vector[i]*=sinc[i]/(M2)*h[i/2];
      vector[i+1]*=sinc[i]/(M2)*h[i/2];
      }
  fft(vector-1,M2,-1);
  for(i=0;i<M;i++) projeccio[kM+i]=vector[2*i];
  }
  
free(sinc); free(vector); free(h);
}















/*---------------------------------- sigma ----------------------------------*/

float sigma(float distancia,int colimador,float t_pixel)
{

float A,B,s;

switch(colimador){
     case 3:A=0.0275;
            B=0.2/t_pixel;
            break;
     case 4:A=0.0172; 
            B=0.2/t_pixel;
            break;
     default:printf("colimador no implementat\n"), exit(0);
     }
     
s=A*distancia/t_pixel+B;

return(s);
}
            


/*=================================== F F T ===============================*/


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fft(float *data,int nn,int isign)
{
int n,mmax,m,j,istep,i;
float wtemp,wr,wpr,wpi,wi,theta;
float tempr,tempi;

n=nn << 1;
j=1;
for(i=1;i<n;i+=2){
   if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
   }
   m=n >> 1;
   while (m>=2 && j>m) {
       j -= m;
       m >>= 1;
   }
 j += m;
}

mmax=2;

while (n > mmax) {
    istep=2*mmax;
    theta=6.28318530717959/(isign*mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m=1;m<mmax;m+=2) {
        for (i=m;i<=n;i+=istep) {
            j=i+mmax;
            tempr=wr*data[j]-wi*data[j+1];
            tempi=wr*data[j+1]+wi*data[j];
            data[j]=data[i]-tempr;
            data[j+1]=data[i+1]-tempi;
            data[i] += tempr;
            data[i+1] += tempi;
		  //if(j>2*nn) printf("j:%d\n",j);
        }
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
}
}

/*==================================================== SOROLL ================*/

/* Rutines de n�meros aleatoris
   Numerical Recipes in C.
   generador de numeros amb distribucio normal */

#define M1    259200
#define IA1   7141
#define IC1   54773
#define RM1   (1.0/M1)
#define M2    134456
#define IA2   8121
#define IC2   28411
#define RM2   (1.0/M2)
#define M3    243000
#define IA3   4561
#define IC3   51349

float ran1(int *idum)
{
   static long ix1,ix2,ix3;
   static float r[98];
   float temp;
   static int iff=0;
   int j;

   if (*idum < 0 || iff==0)
   {
      iff=1;
      ix1=(IC1-(*idum)) % M1;
      ix1=(IA1*ix1+IC1) % M1;
      ix2=ix1 % M2;
      ix1=(IA1*ix1+IC1) % M1;
      ix3=ix1 % M3;
      for (j=1;j<=97;j++)
      {
         ix1=(IA1*ix1+IC1) % M1;
         ix2=(IA2*ix2+IC2) % M2;
         r[j]=(ix1+ix2*RM2)*RM1;
      }
      *idum=1;
   }
   ix1=(IA1*ix1+IC1) % M1;
   ix2=(IA2*ix2+IC2) % M2;
   ix3=(IA3*ix3+IC3) % M3;
   j=1+((97*ix3)/M3);
   if (j>97 || j<1)
   {
      printf("RAN1: error\n");
      exit(1);
   }
   temp=r[j];
   r[j]=(ix1+ix2*RM2)*RM1;
   return temp;
}


/* retorna el valor de ln(gamma(xx)) per xx>0 */

float gammln(float xx)
{
   float x,tmp,ser;
   static float cof[6]=
   {
      76.18009173, -86.50532033, 24.01409822, -1.231739516,
      0.120858003e-2, -0.536382e-5
   };
   int j;

   x=xx-1.0;
   tmp=x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser=1.0;
   for (j=0;j<=5;j++)
   {
      x+=1.0;
      ser+=cof[j]/x;
   }
   return (-tmp+log(2.50662827465*ser));
}

   
/* 
   generador de numeros aleatoris amb distribucio de Poisson 
   de parametre xm
*/
      
float poidev(float xm,int *idum)
{
   static float sq,alxm,g,oldm=(-1.0);
   float em,t,y;
   
   if (xm < 12.0)
   {
      if (xm!=oldm)
      {
         oldm=xm;
         g=exp(-xm);
      }
      em = -1;
      t=1.0;
      do
      {
         em += 1.0;
         t *= ran1(idum);
      } while (t>g);
   }
   else
   {
      if (xm!=oldm)
      {
         oldm=xm;
         sq=sqrt(2.0*xm);
         alxm=log(xm);
         g=xm*alxm-gammln(xm+1.0);
      }
      do
      {
         do
         {
            y=tan(M_PI*ran1(idum));
            em=sq*y+xm;
         } while (em<0.0);
         em=floor(em);
         t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      }  while (ran1(idum) > t);
   }   
   return em;
}    

/*###########################################################################*/
/* modificada por paguiar 2004 porque sobraban par�metros de entrada */

void soroll_projeccio(int ini,float *projeccio,int tamano)
{ 
  int i,idum;
  if (ini==0) return;
  if (ini > 0) idum = -ini;
  for (i=0;i<tamano;i++) projeccio[i]=poidev(projeccio[i],&idum);
}
/*===================================== NOM FITXER RESULTATS ================*/

void nomf(char *fitxer,char *f_res,int p)
{

char h[2];

strcpy(fitxer,f_res);
strcat(fitxer,".");
strcat(fitxer,itoa(p,h));

}   

/*=================================== REDUCCIO_PIXELS ========================*/

void reduccio_pixels(float *lsp,float *le,int Psp,int P,int N)
{
int i,j,k,kP,kPsp,ppp;

if((Psp%P)!=0){ printf("\nErroR: El numero de pixels d'una a de ser divisible pel el numero de pixels de l'altra"); exit(0);}
ppp=Psp/P;

for(k=0,kP=0,kPsp=0;k<N;k++,kP+=P,kPsp+=Psp){
   for(j=0;j<P;j++){
      le[kP+j]=0;
      for(i=0;i<ppp;i++) le[kP+j]+=lsp[kPsp+j*ppp+i];
      }
   }   
}   

/*=================================== REDUCCIO_ANGLES ========================*/

void reduccio_angles(float *lsp,int Nsp,int Psp,int N)
{
int j,k,kPsp,kPsppas,pas;

if((Nsp%N)!=0){ printf("\nError: El numero d'angles ha de ser divisible pel de l'altra"); exit(0);}
pas=Nsp/N;

for(k=0,kPsp=0,kPsppas=0;k<N;k++,kPsp+=Psp,kPsppas+=Psp*pas){
   for(j=0;j<Psp;j++) lsp[kPsp+j]=lsp[kPsppas+j];
   }   
}

/*================================== NORMALITZACIO N_CONTES===================*/

void normalitzacio_nc(float *q,int nitems,float nc)
{
int i;
float suma,factor;

suma=0;
for(i=0;i<nitems;i++) suma+=*(q+i);
factor=nc/suma;
for(i=0;i<nitems;i++) *(q+i)*=factor;

}

/*================================== NO_NEGATIVITAT ==========================*/

void no_negativitat(float *q,int nitems)
{
int i;

for(i=0;i<nitems;i++) if(*(q+i)<0) *(q+i)=0;

}
