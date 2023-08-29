/*Versió: 3/10/00
  fcmax.h q calcula el factor de comptes entre l'activada i la basal 
  ajustant una paràbola al perfil de quocients
  Optimització de cerca del màxim (no falsos pics a l'extrem de l'histograma)*/
/*adaptación a la nueva libreria*/
/*estamos eliminando la mascara*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <util/nrutil.h>
#include <util/factor.h>

/****************calcula el factor entre dos imagenes sin aplicar ninguna mascara *****************/

void fc_max_no_masc(int tipo,struct imagen *bas,struct imagen *ict,double *factor)
{  
   int MM,MMM,V,U,W;
   int i,j,num,NINT=100,P,jmax;
   double nint=100.,max,N,sx,ave;
   double *x,*h,*l;
   double a,b,c,den,new_f,new_h,ix;
   double sx1,sx2,sx3,sx4,sy,sxy,sx2y,n_p;
  /* void determinant();*/
   FILE *file;
tipo=1;
   
V=bas->nfil;
U=bas->ncol;
W=bas->nima;

MM=bas->nfil*bas->ncol;
MMM=bas->nfil*bas->ncol*bas->nima;

P=10*NINT-1;					/*longitud de l'histograma*/
   
   x=dvector(0,MMM-1);
   l=dvector(0,P-1);
   h=dvector(0,P-1);
   
   *factor=0.;
   
   for (i=0;i<P;i++) h[i]=0.;	/*inicialització variables histograma*/
   for (i=0;i<P;i++) l[i]=0.; 
   
   for (i=0;i<MMM;i++) {		/*perfil de quocients d'on es farà l'estima inicial del factor*/
       if (bas->datos[i]>0) {
          /* x[i]=(q2[i]/q1[i]);*/
	  x[i]=	ict->datos[i]/bas->datos[i];
           if (x[i]>0.) {
              num=(int)floor(x[i]*nint+.5);					
              if (num<P) h[num]=h[num]+1.; 
 /*               else  h[P]=h[P]+1.;*/	/*eliminant l'"else" evitem falsos pics q produirien un factor erroni*/
           }
       } 
   } 
   for (i=2;i<P-2;i++) l[i]=(h[i-2]+h[i-1]+h[i]+h[i+1]+h[i+2])/5.; /*histograma suavitzat*/
   
   max=0.;jmax=0;
   
   for (j=0;j<P;j++)  {		/*cerca del màxim*/
       if (max<l[j]) {
           max=l[j];
           jmax=j;
       }
   }

   N=0.;sx=0.;ave=0.;
   
   for (j=jmax-NINT/20;j<=jmax+NINT/20;j++) {		/*promig?*/
       N += h[j];						
       sx += j*h[j]/nint;					
   }
   
   ave=sx/N;
   //*factor=ave;											/*factor d'escala com a promig*/
   printf("\n jmax=%d ave=%f N=%f",jmax,ave,N);

   //for(i=0;i<1000;i++)h[i]=-i*i+1000.*i+10.;
   //jmax=500;
	sx1=0.;sx2=0.;sx3=0.;sx4=0.;sy=0.;sxy=0.;sx2y=0.;		/*ajust  paràbola  y=ax2+bx+c */
	for(i=jmax-NINT/10;i<jmax+NINT/10;i++){
		ix=(float)(i);
		sx1+=ix;
		sx2+=ix*ix;
		sx3+=ix*ix*ix;
		sx4+=ix*ix*ix*ix;
		sy+=h[i];
		sxy+=ix*h[i];
	        sx2y+=ix*ix*h[i];	}
		n_p=nint/5.;
	
	determinant(sx4,sx3,sx2,sx3,sx2,sx1,sx2,sx1,n_p,&den);
	determinant(sx2y,sx3,sx2,sxy,sx2,sx1,sy,sx1,n_p,&a);
	determinant(sx4,sx2y,sx2,sx3,sxy,sx1,sx2,sy,n_p,&b);
	determinant(sx4,sx3,sx2y,sx3,sx2,sxy,sx2,sx1,sy,&c);
	
	a/=den;
	b/=den;
	c/=den;

	new_f=-(b)/(2.*a)/nint;
	new_h=-b*b/(4.*a)+c;

	*factor=new_f;								/*nou factor q prové de l'ajust de la paràbola*/
	//printf("\n a=%f b=%f c=%f new_f=%f new_h=%f",a,b,c,new_f,new_h);
	printf("\n  new_f=%f new_h=%f\n",new_f,new_h);

    free_dvector(x,0,MMM-1);
   
    switch (tipo) 
       {
         
         case 0: break;
         
        
         case 1: file=fopen("cociente.res","w");
                /* if (file==NULL) error_io(106);*/
                 for (j=0;j<P;j++) fprintf(file,"%f %f %f\n",j/100.,h[j],l[j]);
                 fprintf(file,"%f %f",jmax/100.,ave);
                 fclose(file);
               
                 break;
       }   
 free_dvector(l,0,P-1);
 free_dvector(h,0,P-1);  
}
   


/******************** calcula el factor entre dos imagenes aplicando una mascara *****************************/

void fc_max_masc(int tipo,struct imagen *bas,struct imagen *ict,struct imagen *masc, double *factor)

{  
 
   int MM,MMM,V,U,W;
   int i,j,num,NINT=100,P,jmax;
   double nint=100.,max,N,sx,ave;
   double *x,*h,*l;
   double a,b,c,den,new_f,new_h,ix;
   double sx1,sx2,sx3,sx4,sy,sxy,sx2y,n_p;
   /*void determinant();*/
   FILE *file;
   
V=bas->nfil;
U=bas->ncol;
W=bas->nima;

MM=bas->nfil*bas->ncol;
MMM=bas->nfil*bas->ncol*bas->nima;

P=10*NINT-1;					/*longitud de l'histograma*/
   
   x=dvector(0,MMM-1);
   l=dvector(0,P-1);
   h=dvector(0,P-1);
   
   *factor=0.;
   
   for (i=0;i<P;i++) h[i]=0.;	/*inicialització variables histograma*/
   for (i=0;i<P;i++) l[i]=0.; 
   
   for (i=0;i<MMM;i++) {		/*perfil de quocients d'on es farà l'estima inicial del factor*/
       if (masc->datos[i]*bas->datos[i]>0) {
           /*x[i]=(q2[i]/q1[i]);*/
	   
	   x[i]=ict->datos[i]/bas->datos[i];	
           if (x[i]>0.) {
              num=(int)floor(x[i]*nint+.5);					
              if (num<P) h[num]=h[num]+1.; 
 /*               else  h[P]=h[P]+1.;*/	/*eliminant l'"else" evitem falsos pics q produirien un factor erroni*/
           }
       } 
   } 
   for (i=2;i<P-2;i++) l[i]=(h[i-2]+h[i-1]+h[i]+h[i+1]+h[i+2])/5.; /*histograma suavitzat*/
   
   max=0.;jmax=0;
   
   for (j=0;j<P;j++)  {		/*cerca del màxim*/
       if (max<l[j]) {
           max=l[j];
           jmax=j;
       }
   }

   N=0.;sx=0.;ave=0.;
   
   for (j=jmax-NINT/20;j<=jmax+NINT/20;j++) {		/*promig?*/
       N += h[j];						
       sx += j*h[j]/nint;					
   }
   
   ave=sx/N;
   //*factor=ave;											/*factor d'escala com a promig*/
   printf("\n jmax=%d ave=%f N=%f",jmax,ave,N);

   //for(i=0;i<1000;i++)h[i]=-i*i+1000.*i+10.;
   //jmax=500;
	sx1=0.;sx2=0.;sx3=0.;sx4=0.;sy=0.;sxy=0.;sx2y=0.;		/*ajust  paràbola  y=ax2+bx+c */
	for(i=jmax-NINT/10;i<jmax+NINT/10;i++){
		ix=(float)(i);
		sx1+=ix;
		sx2+=ix*ix;
		sx3+=ix*ix*ix;
		sx4+=ix*ix*ix*ix;
		sy+=h[i];
		sxy+=ix*h[i];
	    sx2y+=ix*ix*h[i];	}
		n_p=nint/5.;
	
	determinant(sx4,sx3,sx2,sx3,sx2,sx1,sx2,sx1,n_p,&den);
	determinant(sx2y,sx3,sx2,sxy,sx2,sx1,sy,sx1,n_p,&a);
	determinant(sx4,sx2y,sx2,sx3,sxy,sx1,sx2,sy,n_p,&b);
	determinant(sx4,sx3,sx2y,sx3,sx2,sxy,sx2,sx1,sy,&c);
	
	a/=den;
	b/=den;
	c/=den;

	new_f=-(b)/(2.*a)/nint;
	new_h=-b*b/(4.*a)+c;

	*factor=new_f;								/*nou factor q prové de l'ajust de la paràbola*/
	//printf("\n a=%f b=%f c=%f new_f=%f new_h=%f",a,b,c,new_f,new_h);
	printf("\n  new_f=%f new_h=%f",new_f,new_h);

    free_dvector(x,0,MMM-1);
   
    switch (tipo) 
       {
         
         case 0: break;
         
        
         case 1: file=fopen("cociente.res","w");
                /* if (file==NULL) error_io(106);*/
                 for (j=0;j<P;j++) fprintf(file,"%f %f %f\n",j/100.,h[j],l[j]);
                 fprintf(file,"%f %f",jmax/100.,ave);
                 fclose(file);
               
                 break;
       }   
 free_dvector(l,0,P-1);
 free_dvector(h,0,P-1);  
}
  

/**************************************calcular determinant*****************************************/
void determinant(double a1,double a2,double a3,double a4,double a5,double a6,double a7,double a8,double a9,double *det)

{
	*det=a1*a5*a9+a4*a8*a3+a2*a6*a7-a3*a5*a7-a1*a6*a8-a2*a4*a9;
}
