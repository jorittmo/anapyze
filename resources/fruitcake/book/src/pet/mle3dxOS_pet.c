#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <util/nrutil.h>
#include <util/mp_i4.h>
#include <pet/mp_i4_pet.h>
#include <pet/mle3dxOS_pet.h>
#include <pet/proyectores_pet_hdr.h>
#include <util/convolucion.h>
#include <util/inpout.h>
//#include <util/inpout_double.h>
#define N_INPUTS 15

void error_rec3d_OS(int nerr,char *text);

/*============================= MLExOS =========================================*/

void mlemxOS_pet(float *q,int Nfil,int Ntalls,int Nvox,float *le,float *lc,int Nang,int Nbp,int NbOS,float *n_contes,int OS,char *fmatriu_pes,int NITER,int norma,char *fitxer_ima,char tipout[3])
{
char fitxer_ima_NITER[200]; // nombre del fichero de reconstrucción individualizado por iteracción
int iteracio;		/* index indica la iteracio actual */
int NangOS;		/* nombre d'angles de projeccio a cada subset */  
int Nitems_max;		/* nombre maxims d'elements a les matrius de pesos. Per alloc */                
float *q2;		/* imatge auxiliar */
float *vp;		/* vector pesos (suma per columnes de la matriu mp) */
sparse_i4 mp;		/* variable de la matriu de pesos */
FILE *file;
NangOS=Nang/OS;

// printf("Nfil: %d\nNtalls: %d\n Nvox: %d\n",Nfil,Ntalls,Nvox);
// printf("Nang: %d\n Nbp: %d\n  NboOS: %d\n OS: %d\n NITER: %d\n norma: %d\n",Nang,Nbp,NbOS,OS,NITER,norma);

// recorre las matrices de cada subset y guarda el Nitems_max 
check_nitems_max_mp_pet(fmatriu_pes,OS,&Nitems_max);

// reserva memoria para las tres componentes de la matriz de pesos
if((mp.ar=(float *)calloc(Nitems_max,sizeof(float)))==NULL) error_rec3d_OS(23,"ar");  
if((mp.ja=(int *)calloc(Nitems_max,sizeof(int)))==NULL) error_rec3d_OS(23,"ja");
if((mp.ia=(int *)calloc(NbOS+1,sizeof(int)))==NULL) error_rec3d_OS(23,"ia");    

// reserva de memoria de les variables auxiliares
if((q2=(float *)calloc(Nvox,sizeof(float)))==NULL) error_rec3d_OS(24,"q2");  
if((vp=(float *)calloc(Nvox,sizeof(float)))==NULL) error_rec3d_OS(24,"vp");  

// imagen inicial, un cilindro
// nuevo parametro para hacer cilindro inicial más pequeño
omple_cilindre_mes_petit(q,Nfil,Ntalls,*n_contes/Nvox,3);      

// ejecución de rec_mlexOS_pet para cada iteracción
for(iteracio=0; iteracio<NITER; iteracio ++)
{    
	printf("\n* * * iter: %d * * *\n",iteracio); 
	// reconstrucción   
	rec_mlexOS_pet(q,q2,le,lc,vp,n_contes,OS,Nvox,NbOS,NangOS,norma,fmatriu_pes,&mp);
	
	// Se abre un fichero por cada iteraccion
	sprintf(fitxer_ima_NITER,"%s.iter.%d",fitxer_ima,iteracio);
	if((file=fopen(fitxer_ima_NITER,"w"))==NULL) error_rec3d_OS(31,fitxer_ima);
	output_ima(Nvox,1,q,file,tipout);
	fclose(file);
}

free_mpes(&mp); free(q2); free(vp);
}












/*============================= REC_MLExOS ===================================*/
/* modificada para PET 3D subsets - Pablo/Carles nov2004 */

void rec_mlexOS_pet(float *q,float *q2,float *le,float *lc,float *vp,float *n_contes,int OS,int Nvox,int NbOS,int NangOS,int norma,char *fmatriu_pes,sparse_i4 *mp)
{ 
register int i,k;
int index_le;				//índice sobre los sinogramas experimentales
float n_contes2[1],factor; 		// cuentas calculadas y factor de corrección
char fmatriu_pesxOS[140];		// nombre de la matriz de pesos que se actualiza a cada subset
int *tk,Nvox_kernel; // tamano del kernel PSF para convolucion
float sigma=0.61;  // sigma de la PSF
float *qkernel; // imagen q contiene el kernel
float *cq,*cq2;
float *cvp; // vp convolucionado


if((cq=(float *)calloc(NbOS,sizeof(float)))==NULL) error_rec3d_OS(24,"cq");  
if((cq2=(float *)calloc(NbOS,sizeof(float)))==NULL) error_rec3d_OS(24,"cq2");  

if((cvp=(float *)calloc(Nvox,sizeof(float)))==NULL) error_rec3d_OS(24,"cvp");  
tk=(int *)calloc(1,sizeof(int));

for(k=0;k<OS;k++){

   printf("\nOS::%d\t",k);
   
   // lectura de la matriz de pesos correspondiente al subset -k- 
   sprintf(fmatriu_pesxOS,"%s.%d",fmatriu_pes,k);
   lee_pesos_xOS(mp,fmatriu_pesxOS);
    
   // se obtiene el vector de pesos de la matriz (es la imagen de sensibilidad)
   vector_pesosxOS(vp,mp,Nvox);

   // CONVOLUCION DE VP CON PSF   
   tam_kernel(tk,&sigma);
   Nvox_kernel=(*tk)*(*tk)*1;
   if((qkernel=(float *)calloc(Nvox_kernel,sizeof(float)))==NULL) error_rec3d_OS(24,"vp");  
   gaussian_kernel_q(sigma, qkernel,*tk,*tk,1);
   convolucion_q(vp,cvp,qkernel,64,64,69,*tk,*tk,1);

   // se inicializa una matriz auxiliar 
   for(i=0;i<Nvox;i++) q2[i]=0;  
  
   // se projecta hacia adelante - forward projector
   convolucion_q(q,cq,qkernel,64,64,69,*tk,*tk,1);
   projeccio(lc,cq,-1,mp,NbOS);
   // CONVOLUCION DE Q CON PSF PERO ANTES DEL FORWARD
   

   //sprintf(nombre,"fwd_proj.subset.%d",k);
   //if((file_in=fopen(nombre,"w"))==NULL) error_rec3d_OS(31,"");
   //output_ima(NbOS,1,lc,file_in,"fl");
   //fclose(file_in);
   
   //for(ii=0;ii<NbOS;ii++) printf("%f**",lc[ii]);
   
   // Guarda imagen calculada Xn(i) y su proyección A(i,j)*Xn(i)
   //sprintf(nombre,"fwd_proj.subset.%d",k);
   //if((file_in=fopen(nombre,"w"))==NULL) error_rec3d_OS(31,"");
   //output_ima(NbOS,1,lc,file_in,"fl");
   //fclose(file_in);
   //sprintf(nombre,"object.subset.%d",k);
   //if((file_in=fopen(nombre,"w"))==NULL) error_rec3d_OS(31,"");
   //output_ima(Nvox,1,q,file_in,"fl");
   //fclose(file_in);

   // se calcula el cociente de proyecciones     
   for(i=0,index_le=k*NbOS ; i<NbOS ; i++,index_le++)
   {
	if(lc[i]==0) lc[i]=1.;  //evitando dividir por cero
	lc[i]= le[index_le]/lc[i];    //index_le diferente del index de lc para leer el subset correspondiente
   }

   // se proyecta hacia atrás las proyecciones calculadas (lc)
   per_mt(lc,q2,mp,Nvox);
   convolucion_q(q2,cq2,qkernel,64,64,69,*tk,*tk,1);
   // CONVOLUCION DE Q2 CON PSF DESPUES DEL BACKPROJECTOR

   // q2 es la imagen obtenida hasta el momento, se actualiza dividiendo por 
   // la imagen de sensibilidad.
   for(i=0;i<Nvox;i++) if(cvp[i]!= 0)  q[i] *= cq2[i]/cvp[i];
          
   // normaliza al número de cuentas. 
   if(norma)
   { 
      printf("*N*");
      projeccio(lc,q,-1,mp,NbOS);     
      numeroc(lc,n_contes2,NbOS);
      factor=n_contes[k]/n_contes2[0];
      for(i=0;i<Nvox;i++) q[i]*=factor;
   }
   }
}






/* ============================= NUMERO DE COMPTES PER SUBSET ==============*/
/* modificada para PET 3d subsets Pablo y Carles nov2004 */

void numerocpxOS(float *le,float *n_comptes,int OS,int NbOS)
{
int i,j,k;
float nct=0.;

// inicializa vector de cuentas
for(k=0;k<OS;k++) n_comptes[k]=0;

for(k=0,j=0;k<OS;k++){
    for(i=0;i<NbOS;i++,j++) n_comptes[k]+=le[j];
    nct+= n_comptes[k];
    }
//printf("Reconstruyendo sinogramas con OSEM3D con cuentas totales: %f\n",nct);
}

/* ============================= CHECK NITEMS MAX MPES ==============*/
// Esta funcion lee el num de elementos de cada mp asociada a cada subset y escribe el maximo 
void check_nitems_max_mp_pet(char *fmatriu_pes,int OS,int *Nitems_max)
{
char fmatriu_pesxOS[240];
int p,aux;
FILE *ini;

*Nitems_max=0;

for(p=0;p<OS;p++)
{
	/* ..... generacio del nom de la matriu de pesos .... */
	sprintf(fmatriu_pesxOS,"%s.%d",fmatriu_pes,p);
	
	if((ini = fopen(fmatriu_pesxOS,"r"))==NULL) error_mp(11,fmatriu_pesxOS);
	
	/* se salta dos valores hasta leer el numero de elementos de la matriz que es el tercero 
	Al introducir nuevos campos en sparse ahora hay cambios pues hay que saltar más */
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);	
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	fread (&aux,sizeof(size_t),1,ini);
	//fread (&aux,sizeof(size_t),1,ini);
	fclose(ini);

	if(aux > *Nitems_max) *Nitems_max=aux;
}	
} 
/*=========================================== ERROR =========================*/

void error_rec3d_OS(int nerr,char *text)
{
printf("\n\nERROR EN EL PROGRAMA rec3d\n");
switch(nerr){
    case 23: printf("\n\nProblemas de calloc en mle3dxOS al reservar memoria para el elemento %s de struct mp\n",text); break;
    case 24: printf("\n\nProblemas de calloc en mle3dxOS al reservar memoria para el vector auxiliar %s \n",text); break;
    case 31: printf("\n\nError al abrir el fichero de resultados");
             printf("\nEl fichero  %s  no se puede abrir\n",text); break;
    }
exit(0);
}  
