// Este programa hace un ajuste bidimensional utiliza una c cte que
//es el valor mínimo de la parte central de la imagen y no se adapta. 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <util/amoeba.h>
#include <util/ajgausbi_girado_auto2_punt.h>
#include <util/enc_centros.h>
#include <util/nrutil.h>
#define pi 4.*atan(1.)

double tp,x1cm,y1cm,sgx1;
double ang=0.;
double *sf;
double *t;
int n = 0,ndim = 3,ndim2,npar=14, z=0;

/*================================CALCULO DEL CENTRO DE MASA=====================================*/
void centroide(int ntt,double tp, char nom[], double *x1cm, double *y1cm, double *sg1x, double *sg1y, int n3, int *punts, int ntalls, int tall, double xx, double yy, double *r1)
{
FILE *ima;

//unsigned int *iq;
double xcm=0.,ycm=0.,sgx=0.,sgy=0.,rr=0.;
float *v,*vtot, *ampli;
int i=0,d=0,puntos=0,*vectx,*vecty,kk,ini=0,fin=0;



//Abrimos la imagen que nos dan desde el programa principal.		
ima=fopen(nom,"rb");

d=ntt*ntt;
//y1=ntt/3;
//y2=2*ntt/3;

if(ima==NULL){   printf("\n Problemes amb els fitxers de dades!"); printf("\n");exit(-1);   }

//Reservamos memoria.
//iq=(unsigned int*)calloc(sizeof(unsigned int),d);
vtot=(float*)calloc(4,d*ntalls);
v=(float*)calloc(4,d);

n=n3;

//Leemos el fichero y lo guardamos en vtot
// Copiamos los datos del corte a utilizar sobre v
fread(vtot,4,d*ntalls,ima);
ini = d*(tall-1);
fin = d*tall;
for(kk=ini;kk<fin;kk++) v[kk-ini]=vtot[kk];
free(vtot);


//for(i=0;i<d;i++) v[i]=(double)iq[i];
//free(iq);

vectx=(int*)calloc(sizeof(int),ntt);
vecty=(int*)calloc(sizeof(int),ntt);
ampli=(float*)calloc(sizeof(float),ntt);

for (i=0;i<ntt;i++) 
{
	vectx[i]=0;
	vecty[i]=0;
	ampli[i]=0.;
}

//Llamamos a calcentro, y le pasamos el nº total de filas, la fila
//inicial y final y el vector de datos. Mientras que nos devolvera 
//x.c.m y la y.c.m. 
//encuentra(ntt,&puntos,vectx,vecty,v,&minimo,ampli,rad_pix,fov);
//printf(" \n Se encontraron %i puntos \n ",puntos);
//printf(" \n Nº de puntos es  %i \n ", puntos);

//for (i=0;i<puntos;i++) {

x1cm[i]=0.;
y1cm[i]=0.;
sg1x[i]=0.;
sg1y[i]=0.;

//for (i=0;i<puntos;i++) {

printf(" \n Vamos por el pto %i \n ", i);
calcentro3(ntt,xx,yy,&xcm,&ycm,&sgx,&sgy,v,n3,&rr);
printf(" \n Vamos por el pto\n ");
x1cm[0]=xcm;
y1cm[0]=ycm;
//	x1cm[i]=(x1cm[i]-149.5)*tp;
//	y1cm[i]=(y1cm[i]-149.5)*tp;
//	y1cm[i]=-y1cm[i];
sg1x[0]=sgx;
sg1y[0]=sgy;
r1[0]=rr;
sg1x[0]=sg1x[0]*tp;
//sg1y[i]=sgy;
sg1y[0]=sg1y[0]*tp;

//}

//Cerramos el archivo y liberamos espacio en memoria.
punts[0]=puntos;
/*
printf(" \n La componente x es: \n ");
for (i=0;i<puntos;i++) {printf(" %f ",x1cm[i]);}
printf(" \n La componente y es: \n ");
for (i=0;i<puntos;i++) {printf(" %f ",y1cm[i]);}
*/
fclose(ima);
free(v);

}







/* ===============================  CALCULA CENTRO DE GAUSSIANA  ================================= */
void calcentro3(int ntt, double vectx, double vecty, double *xcm, double *ycm, double *sgx, double *sgy, float *v, int n, double *rr)
{

extern double ang;

int j,i,x1,x2,y1,y2,d,n2;
double sv,sxv,syv;
double a,c,mx,my,sigmax,sigmay,r,maxc,maxf;
float *les, posx, posy, mod;
double mrad, mtg;
FILE *perfil, *punto;

sv=sxv=syv=0.;
maxc=0.; maxf=0.;

d=ntt*ntt;
n2=n*n;

//printf(" \n%i %i \n",vectx,vecty);
//Utilizamos esta posición como referencia para definir un cuadrado
//alrededor del punto y ajustar una gaussiana.
x1=vectx-n/2; x2=vectx+n/2; y1=vecty-n/2; y2=vecty+n/2;
if (x1<0) {x1=0; x2=n;}
if (x2>ntt) {x1=ntt-n; x2=ntt;}
if (y1<0) {y1=0; y2=n;}
if (y2>ntt) {y1=ntt-n; y2=ntt;}
sf=(double *)calloc(sizeof(double),n2);

//sc=(double *)calloc(sizeof(double),n);
posx=(float)vectx-(float)ntt/2.+1.;
posy=(float)ntt/2.-(float)vecty-1.;

if (posx == 0.) {ang=1.5708;} else{ ang=atan(-posy/posx); }

ang=ang*180.;

ang=pi/ang;
ang=1./ang;


//ang=ang*180./pi;
mod=sqrt(posx*posx+posy*posy);

printf(" \n El angulo, mod, posx y posy son: %f %f %f %f",ang,mod,posx,posy);
//if (mod<(float)ntt/12.) { ang=0.; }
ang=ang*pi/180.;
//printf(" \n El angulo es: %f",ang);
///////////
//ang=-0.7853981633974483;
//////////
									/*sumar filas i columnas*/
for(i=0;i<n2;i++){	sf[i]=0.; 	}		/*inicialitzacion*/
					
for(i=y1;i<y2;i++){	
	for(j=x1;j<x2;j++){  			
		          sf[j-x1+n*(i-y1)]=(double)v[ntt*i+j];
		          }  
		}


//Muestra el punto
punto=fopen("punt.img","wb");
les=(float *)calloc(sizeof(float),n2);
for (i=0;i<n2;i++){les[i] = (float)sf[i];}		
fwrite((void*)les,sizeof(float),n2,punto);
free(les);
fclose(punto);
//////******//////

perfil=fopen("perfil_fwhm.xls","wb");

fprintf(perfil,"\n");
for (i=0;i<n2;i++) {fprintf(perfil,"%f \t",sf[i]);}
fclose(perfil);
sigmax=2.; sigmay=2.; mx=vectx-x1-1.; my=vecty-y1-1.; 
/////////
c=0.;
a=2.;

/////////
//printf ("La c antes es: %f",c);
//my=-my;
mrad=mx*cos(ang)+my*sin(ang);
mtg=my*cos(ang)-mx*sin(ang);
//printf(" \n%f %f %f\n",a,mx,my);
//valors_inicials(0,sf,&a,&m,&sigma,&biaix,&curtosi,tp,n);

ajgausbi(&a,&mrad,&mtg,&sigmax,&sigmay,&c,&r);



//sigmax=sigma;rx=r;mx=m;ax=a;
/*for(i=0;i<n;i++) sf[i]=sc[i];
valors_inicials(0,sf,&a,&m,&sigma,&biaix,&curtosi,tp,n);
ajgaus(&a,&m,&sigma,&r);
*/
mx=mrad*cos(ang)-mtg*sin(ang);
my=mrad*sin(ang)+mtg*cos(ang);

//printf ("La c es: %f",c);
xcm[0]=mx+ (double) x1+1.;
ycm[0]=my+ (double) y1+1.;
sgx[0]=sigmax;
sgy[0]=sigmay;
rr[0]=r;
//sgy[0]=sigma;


return;
}









/*===================== CALCULA amplada, mitjana i sigma guassiana amb  SIMPLEX ===================*/	

void ajgausbi(double *a,double *mrad,double *mtg,double *sigmax,double *sigmay, double *c,double *r)
{
extern double suma_mgbi(double *), *sf;
extern int n;
extern double ang;
	int i,j,npar=6,iter,n2,cx,cy;
    	double ftol=1.e-7,vmax;
    	double *y,**p,*lambda,*x,sx,sx2,sy,sy2,sxy,sigma0x,sigma0y,c0,vg,distx,disty;
    	double crad,ctg,sdif;

	

y =dvector(1,npar+1);				/*   reserva memoria			*/
lambda=dvector(1,npar);				/*   reserva memoria			*/
p=dmatrix(1,npar+1,1,npar);			/*   reserva memoria			*/
x =dvector(1,npar);					/*   reserva memoria			*/

lambda[1]=fabs(*a)/2.;lambda[2]=fabs(*mrad)/2.;lambda[3]=fabs(*mtg)/2.;lambda[4]=fabs(*sigmax)/2.;lambda[5]=fabs(*sigmay)/2.;
//lambda[6]=fabs(*c)/2.;
lambda[6]=0.;
for(i=1;i<=npar+1;i++){ p[i][1]=*a;p[i][2]=*mrad;p[i][3]=*mtg;p[i][4]=*sigmax;p[i][5]=*sigmay;p[i][6]=*c;}
for(i=2;i<=npar+1;i++) p[i][i-1]+=lambda[i-1];

for(i=1;i<=npar+1;i++)
	{
	for(j=1;j<=npar;j++)  x[j]=p[i][j]; 
	y[i]=(*suma_mgbi)(x);
	}

amoeba(p,y,npar,ftol,suma_mgbi,&iter);

j=1;vmax=10000.;

for (i=1;i<=npar+1;i++) 
	{
	if (y[i]<vmax) 
		{
		vmax=y[i];
		j=i;
		}
	}

*a=p[j][1];
*mrad=p[j][2];
*mtg=p[j][3];
*sigmax=p[j][4];
*sigmay=p[j][5];
*c=p[j][6];
n2=n*n;
sx=sy=sx2=sy2=sxy=sdif=0.;
sigma0x=*sigmax;
sigma0y=*sigmay;
c0=*c;
for(i=0;i<n2;i++)
	{
		cx= i%n;
		cy= i/n;

		crad=(double)cx*cos(ang)+(double)cy*sin(ang);
		ctg=(double)cy*cos(ang)-(double)cx*sin(ang);
		
			distx=(crad-mrad[0])*(crad-mrad[0]);
			disty=(ctg-mtg[0])*(ctg-mtg[0]);
			
	sx+=sf[i];
	sx2+=sf[i]*sf[i];
	vg=*a*exp(-distx/(2.*sigma0x*sigma0x))*exp(-disty/(2.*sigma0y*sigma0y))+c0;
	sy+=vg;
	sy2+=vg*vg;
	sxy+=sf[i]*vg;
/*	printf(" \n sf: %f\n ",sf[i]);
	printf(" \n vg: %f\n ",vg);
	printf(" \n sy2: %f\n ",sy2);
	printf(" \n sxy: %f\n ",sxy);
*/
	}
	printf("\n ddd  22");
	printf(" \n sx: %f\n ",sx);
	printf(" \n sx2: %f\n ",sx2);
	printf(" \n sy: %f\n ",sy);
	printf(" \n sy2: %f\n ",sy2);
	printf(" \n sxy: %f\n ",sxy);
	printf(" \n n: %i\n ",n);

//*r=(sxy/(double)n-sx*sy)/(sqrt((sx2/(double)n-sx*sx)*(sy2/(double)n-sy*sy)));
//*r=(sxy-sx*sy/(float)n)/(sqrt((sx2-sx*sx/(float)n)*(sy2-sy*sy/(float)n)));
r[0]=(sxy-sx*sy/(float)n2)/(sqrt(sx2-sx*sx/(float)n2)*sqrt(sy2-sy*sy/(float)n2));
//*r=((double)n2*sxy-sx*sy)/(sqrt((double)n2*sx2-sx*sx)*sqrt((double)n2*sy2-sy*sy));
r[0]=r[0]*r[0];	

	
free_dvector(y,1,npar+1);			/*   allibera memoria			*/
free_dvector(lambda,1,npar);		/*   allibera memoria			*/
free_dmatrix (p,1,npar+1,1,npar);	/*   allibera memoria			*/
free_dvector(x,1,npar);				/*   allibera memoria			*/

}



/* ======================================= SUMA RESIDUS========================================= */
double suma_mgbi(double *x)
{
	extern double ang;
	extern int n;	
        extern double *sf;
        
        double a,c,sigmax,sigmay,s,v,mrad,mtg,distx,disty,vx,vy;
        int k,n2,cx,cy;
	double crad,ctg;
	//printf(" \n ang: %f \n ",ang);
	//printf(" \n n: %i\n ",n);	
		n2=n*n;
        a=x[1];mrad=x[2];mtg=x[3];sigmax=x[4];sigmay=x[5];
		c=x[6];

        s=0.;
        for (k=0;k<n2;k++){
			cx= k%n;
			cy= k/n;
			crad=(double)cx*cos(ang)+(double)cy*sin(ang);
			ctg=(double)cy*cos(ang)-(double)cx*sin(ang);
		
			distx=(crad-mrad)*(crad-mrad);
			disty=(ctg-mtg)*(ctg-mtg);
	
			vx=a*exp(-distx/(2.*sigmax*sigmax));
			vy=exp(-disty/(2.*sigmay*sigmay));
            v=vx*vy-sf[k]+c;
            s+=v*v;

        }
			
	
        return (s);
}



