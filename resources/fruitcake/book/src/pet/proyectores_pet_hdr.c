#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <time.h>
#define EPSILON 0.0001
#include <util/utiles_vectores.h>
#include <util/inpout.h>
#include <util/nrutil.h>
#include <pet/proyectores_pet_hdr.h>
#include <pet/d_punto_lor.h>
#include <pet/STIR.h>
#include <util/mp_i4.h>
#include <pet/mp_i4_pet.h>        
//#define PESMINIM 0.01
#define maximo(a,b) ((a)>=(b) ? (a):(b))
#define minimo(a,b) ((a)<=(b) ? (a):(b))
#define SIGNE(a) (a<-EPSILON?-1:(a>EPSILON?1:0))


// Global variables
extern int Nvox,Nbt,Nbp,Nang;
extern float *q;
int MM,PN,M,N,P;


/* Subrutina que calcula los puntos de corte entre una LOR cualquiera y la malla 3D donde situamos el objeto
la salida siempre son dos puntos de tres componentes paguiar, Enero 2004*/
/*
Los parametros de entrada hacen referencia a la recta respecto la cual voy a calcular
los puntos de corte de entrada y de salida de la LOR en el objeto.
*/

/* Esta es la ecuacion de la recta que describe la LOR(n,teta,phi ó z1-z2):

x=-n*seno + lambda*coseno
y=n*coseno + lambda*seno
z=z_med + lambda*cotan_phi (cuando phi es 90º tengo plano horizontal, con phi' sería ..*tan_phi' y 0º para el plano horizontal)

*/









/********** 



compara 






***********/

void compara(float *dx,float *dy,float *dz,int *cas)
{
if(*dx==*dy)*dx+=EPSILON;
if(*dx==*dz)*dx+=EPSILON;
if(*dy==*dz)*dy+=EPSILON;

if(*dx<*dy){
  if(*dx<*dz) *cas=3; //dx es el menor
  else *cas=1; //dz es el menor
  }
else{
  if (*dy<*dz) *cas=2; //dy es el menor
  else *cas=1;
  }
}









/***************************** 







Calcula los ray_sum 




****************************************************************************/

/* Pequeña modificacion de calcul_ate() en SPECT para que devuelva el valor de la proyeccion en lugar de la atenuacion */
float proy_pet3d(struct coord3d *c1,struct coord3d *c2,struct imagen *ima,float dmax)
{
float ray_sum,adx,vx,vy,vz,dpb,dant,dx,dy,dz,next_x,next_y,next_z;
float A,B,C,AR;  //,a lado del FOV;
float tpix,lvox;
int ip,jp,kp,ind_map;
int incrx,incry,incrz;
int cas,esf;
int i,j,k;
float limite; //valor epsilon para comparar valores y dar signos
float xpix,ypix,zpix,xbin,ybin,zbin;
int Nfil,Ntalls,Nfild2,Npix;

limite=0.0001;

xpix=c1->x;
ypix=c1->y;
zpix=c1->z;
xbin=c2->x;
ybin=c2->y;
zbin=c2->z;
Nfil=ima->nfil;
Ntalls=ima->nima;
tpix=ima->tpix/10.; /* hdr en mm ahora necesito cm */
lvox=ima->tcorte/10.; /* hdr en mm ahora necesito cm */
Nfild2=0.5*ima->nfil;
Npix=Nfil*Nfil;

//calculo de vectores unitarios en cada dirección
//signo es -1 si vx es negativo, 1 si vx es positivo y 0 si vx es cero
vx=xbin-xpix;
incrx=signo(vx,limite);
vy=ybin-ypix;
incry=signo(vy,limite);
vz=zbin-zpix;
incrz=signo(vz,limite);
dpb=sqrt(vx*vx+vy*vy+vz*vz); 
vx/=dpb; //
vy/=dpb;
vz/=dpb;

// Se comprueba si la LOR interseca a la esfera de radio a contenida en el FOV lambda=-B+-sqrt(B2-4AC)/2A
// pretende evitar calculos en las esquinas del FOV
A=(vx*vx)+(vy*vy);
B=2*(vx*xpix+vy*ypix);
C=(xpix*xpix)+(ypix*ypix)-(Nfild2*tpix*Nfild2*tpix);
AR=B*B-(4*A*C);  // argumento de la raiz

if(AR<0)
{	
	esf=0; //no entra en el while y entonces adx=0 y ray_sum=1
}
else
{
	esf=1; //entrará en el while
}



//Indices que corresponden a la coordenada c1 (punto inicial)
ptos_corte_en_ijk(&i,&j,&k,c1,lvox,tpix,Nfild2);

// Índices del voxel siguiente en la dirección de la LOR
if(incrx<=0) ip=i;
else ip=i+1;
if(incry<=0) jp=j;
else jp=j+1;
if(incrz<=0) kp=k;
else kp=k+1;
next_x=(-Nfild2+ip)*tpix;
next_y=(-Nfild2+jp)*tpix;
next_z=(kp)*lvox;   

// Calculo el paso inicial
if(fabs(vx)>EPSILON) dx=(next_x-xpix)/vx;
else dx=dmax;
if(fabs(vy)>EPSILON) dy=(next_y-ypix)/vy;
else dy=dmax;
if(fabs(vz)>EPSILON) dz=(next_z-zpix)/vz;
else dz=dmax;

dant=0;
adx=0;
ind_map=k*Npix+i*Nfil+j; //Npix son los cortes, Nfil son filas y columnas

while(ip>=0 && ip<=Nfil && jp>=0 && jp<=Nfil && kp>=0 && kp<=Ntalls && esf==1)
{
  //  devuelve el valor "cas" segun valores de dx dy dz
  compara(&dx,&dy,&dz,&cas);

  switch(cas){
     case 1: // dz es el menor
		adx+=(dz-dant)*ima->datos[ind_map];
          dant=dz;
          kp+=incrz;
          ind_map+=incrz*Npix;
          next_z+=lvox*incrz;
          dz=(next_z-zpix)/vz;
          break;
     case 2: // dy es el menor
		adx+=(dy-dant)*ima->datos[ind_map];
          dant=dy;
          jp+=incry;
          ind_map+=incry;
          next_y+=tpix*incry;
          dy=(next_y-ypix)/vy;
	  	break;
     case 3: // dx es el menor
          adx+=(dx-dant)*ima->datos[ind_map];
		dant=dx;
          ip+=incrx;
          ind_map+=incrx*Nfil;
          next_x+=tpix*incrx;
          dx=(next_x-xpix)/vx;
          break;
      default:
          printf("\nError al comparar diferenciales en proyectores_pet.h\n");
	  exit(0);
          }
}
ray_sum=exp(adx);
return(ray_sum);
}






/***************************** 




Calcula los ray_sum 





****************************************************************************/

// calcula rayo a va asignando 'pesos' a medida que avanza por el rayo, en este caso tiene en cuenta los detect planos (Siddon)
int ray_tracing_lor(struct coord3d *c1,struct coord3d *c2,sparse_i4 *mp, int Nfil, int Ntallsima, float tpix, float lvox, float dmax, int ind_proj, int Nb_os, int fac_comp, float sum)
{
float vx,vy,vz,dpb,dant,dx,dy,dz,next_x,next_y,next_z;
int ip,jp,kp,ind_map;
int incrx,incry,incrz;
int cas;
int i,j,k,mp_ne_ini;;
float limite; //valor epsilon para comparar valores y dar signos
float xpix,ypix,zpix,xbin,ybin,zbin,ri;
int Nfild2,Npix,out_cyl;

limite=0.0001;

xpix=c1->x;
ypix=c1->y;
zpix=c1->z;
xbin=c2->x;
ybin=c2->y;
zbin=c2->z;
Nfild2=0.5*Nfil;
Npix=Nfil*Nfil;
//printf("\ncoord1:%f,%f,%f coord2:%f,%f,%f",c1->x,c1->y,c1->z,c2->x,c2->y,c2->z);
// control para los puntos sobre las esquinas de la malla
if(xbin==xpix && ybin==xbin && zbin==zpix) 
{
	xpix+=10.*tpix;
	ypix+=10.*tpix;
	zpix+=tpix;
}

//calculo de vectores unitarios en cada dirección
//signo es -1 si vx es negativo, 1 si vx es positivo y 0 si vx es cero
vx=xbin-xpix;
vy=ybin-ypix;
vz=zbin-zpix;
incrx=signo(vx,limite);
incry=signo(vy,limite);
incrz=signo(vz,limite);
dpb=sqrt(vx*vx+vy*vy+vz*vz); 
vx/=dpb; 
vy/=dpb;
vz/=dpb;


//Indices que corresponden a la coordenada c1 (punto inicial) para ray-tracing
ptos_corte_en_ijk_rt(&i,&j,&k,c1,lvox,tpix,Nfild2);
//printf("\ncoord1:%f,%f,%f indexes:%d,%d,%d",c1->x,c1->y,c1->z,i,j,k);
// Índices del voxel siguiente en la dirección de la LOR
if(incrx<=0) ip=i;
else ip=i+1;
if(incry<=0) jp=j;
else jp=j+1;
if(incrz<=0) kp=k;
else kp=k+1;

next_x=((ip-Nfild2)*tpix)+(0.5*tpix);        
next_y=((jp-Nfild2)*tpix)+(0.5*tpix);
next_z=((kp)*lvox)+(0.5*lvox);   
//next_z=((kp)*lvox);   

// Calculo el paso inicial
if(fabs(vx)>EPSILON) dx=(next_x-xpix)/vx;
else dx=dmax;
if(fabs(vy)>EPSILON) dy=(next_y-ypix)/vy;
else dy=dmax;
if(fabs(vz)>EPSILON) dz=(next_z-zpix)/vz;
else dz=dmax;

dant=0;
ind_map=k*Npix+i*Nfil+j; //Npix son los cortes, Nfil son filas y columnas
mp_ne_ini=mp->ne;
if(fac_comp==666) mp->ne+=ind_map;
out_cyl=0;
while(ip>=0 && ip<=Nfil && jp>=0 && jp<=Nfil && kp>=0 && kp<=Ntallsima)
{
  
  //  devuelve el valor "cas" segun valores de "dx-dy-dz"
  compara(&dx,&dy,&dz,&cas);
  
  //comprueba que esté dentro del cilindro-FOV
  ri=sqrt(((next_x-(0.5*tpix))*(next_x-(0.5*tpix)))+((next_y-(0.5*tpix))*(next_y-(0.5*tpix))));
  if(ri>=(Nfild2-2.5)*tpix) out_cyl=1;
  
  //imprime valores
  
  switch(cas)
  {
     case 1: // dz es el menor
	 
		  if(fac_comp==666)
		  {	 	
			mp->ar[mp->ne]=sum*fabs(dz-dant);
			if(out_cyl==1) mp->ar[mp->ne]=0.;
			mp->ja[mp->ne]=ind_map;
			if(incrz>=0) mp->ne+=Npix;
			if(incrz<0) mp->ne-=Npix;
		  }
		  
		  if(fac_comp!=666 && out_cyl==0 && dz!=dant)
		  {
		  
			mp->ar[mp->ne]=sum*fabs(dz-dant);
			mp->ja[mp->ne]=ind_map;
			mp->ne++;
		  }
		            
		  dant=dz;
          kp+=incrz;
		  ind_map+=incrz*Npix;
          next_z+=lvox*incrz;
          dz=(next_z-zpix)/vz;
		  break;
     
	 case 2: // dy es el menor
	 
	 	  if(fac_comp==666)
		  {
			mp->ar[mp->ne]=sum*fabs(dy-dant);
			if(out_cyl==1) mp->ar[mp->ne]=0.;
			mp->ja[mp->ne]=ind_map;
			if(incry>=0) mp->ne++;
			if(incry<0) mp->ne--;	
		  }
		  
		  if(fac_comp!=666 && out_cyl==0 && dy!=dant)
		  {
			mp->ar[mp->ne]=sum*fabs(dy-dant);
			mp->ja[mp->ne]=ind_map;
			mp->ne++;
		  }
		  
		  dant=dy;
          jp+=incry;
          ind_map+=incry;
          next_y+=tpix*incry;
          dy=(next_y-ypix)/vy;
		  break;

     case 3: // dx es el menor
	 	  
		  if(fac_comp==666)
		  {
			mp->ar[mp->ne]=sum*fabs(dx-dant);			
			if(out_cyl==1) mp->ar[mp->ne]=0.;
			mp->ja[mp->ne]=ind_map;	  
			if(incrx>=0) mp->ne+=Nfil;
			if(incrx<0) mp->ne-=Nfil;          
		  }
		  
		  if(fac_comp!=666 && out_cyl==0 && dx!=dant)
		  {
			mp->ar[mp->ne]=sum*fabs(dx-dant);
			mp->ja[mp->ne]=ind_map;	  
//printf("\n%d/mp.ar:%f mp.ja:%d  (%d, %d %d) (%d, %d %d)",mp->ne,mp->ar[mp->ne],mp->ja[mp->ne],i,j,k,ip,jp,kp);		  			
			mp->ne++;
		  }
		  
		  dant=dx;
          ip+=incrx;
          ind_map+=incrx*Nfil;
          next_x+=tpix*incrx;
          dx=(next_x-xpix)/vx;
		  break;
      
	  default:
          printf("\nError al comparar diferenciales en proyectores_pet.h\n");
	  exit(0);
   }
   
   out_cyl=0;

}

if(fac_comp==666) mp->ne+=((Npix*Ntallsima)-(mp->ne-mp_ne_ini));
return(0);
}



int check_if_exist(sparse_i4 *mp, int mp_ne_ini,int indice,int indice_actual,int *s)
{
	int i;
	i=indice_actual;
	while(i>=mp_ne_ini)
	{
		if(mp->ja[i]==indice) 
		{
			*s=i;
		}
		i--;
	}
	
	return 0;
}


void	set_value_mp(sparse_i4 *mp,float sum,float dplor,int indice, float lim_exp)
{
	//printf("dplor:%f pes:%f factor:%f indice:%d\n",dplor,1.-(dplor/lim_exp),sum,indice);
	mp->ar[mp->ne]=sum*(1.-(dplor/lim_exp));
	mp->ja[mp->ne]=indice;
	mp->ne++;
}


void	expansion_X(int npix_exp2d, float next_x, float next_y, float next_z, float tpix, float lvox, struct lor_params *lor, int Nfil, int Ntallsima, int ip, int ind_map, int Npix, sparse_i4 *mp, float sum,float lim_exp,int *mask_ima)
{
	int nex,indice,Nfild2,ii,jj,kk;
	float ri,dplor;
	Nfild2=0.5*Nfil;
	nex=npix_exp2d;
	struct coord3d cp,cc;
	
	while(nex>=0)
	{
		//printf("\texp_X-nex:%d: x:%f y:%f z:%f indice:%d\n",nex,next_x-(nex*tpix),next_y,next_z,ind_map-(nex*Nfil)); 
		ri=sqrt(((next_x-(nex*tpix))*(next_x-(nex*tpix)))+((next_y-(0.*tpix))*(next_y-(0.*tpix))));
  		if(ri<=(Nfild2+0.5)*tpix && ip>=nex && ip<Nfil-nex) 
		{
			indice=ind_map-(nex*Nfil);
			if(mask_ima[indice]==0)
			{
				cp.x=next_x-(nex*tpix);
				cp.y=next_y;
				cp.z=next_z;
				ptos_corte_en_ijk_pWu_fast(&ii,&jj,&kk,&cp,lvox,tpix,Nfild2);
				ptos_centro_voxel_pWu_fast(&cc,ii,jj,kk,lvox,tpix,Nfild2);
				dplor=d_punto_lor(cc.x,cc.y,cc.z,lor->n,lor->seno,lor->coseno,lor->phi,lor->z_med,lor->mod_vlor);
				//dplor=d_punto_lor(next_x-(nex*tpix),next_y,next_z,lor->n,lor->seno,lor->coseno,lor->phi,lor->z_med,lor->mod_vlor);
				/*printf("dplor:%f\n",dplor);*/  
				if(dplor<=lim_exp) 
				{
					set_value_mp(mp,sum,dplor,indice,lim_exp);
					//printf("xpix:%f ypix:%f zpix:%f n:%f phi:%f z_med:%f \n",next_x-(nex*tpix),next_y,next_z,lor->n,lor->phi,lor->z_med);  
				}
				mask_ima[indice]=1;
			}
		}
		//printf("\texp_X-nex:%d: x:%f y:%f z:%f indice:%d\n",nex,next_x+(nex*tpix),next_y,next_z,ind_map+(nex*Nfil)); 				
		ri=sqrt(((next_x+(nex*tpix))*(next_x+(nex*tpix)))+((next_y-(0.*tpix))*(next_y-(0.*tpix))));
  		if(ri<=(Nfild2+0.5)*tpix && ip>=nex && ip<Nfil-nex) 
		{
			indice=ind_map+(nex*Nfil);
			if(mask_ima[indice]==0)
			{
				cp.x=next_x+(nex*tpix);
				cp.y=next_y;
				cp.z=next_z;
				ptos_corte_en_ijk_pWu_fast(&ii,&jj,&kk,&cp,lvox,tpix,Nfild2);
				ptos_centro_voxel_pWu_fast(&cc,ii,jj,kk,lvox,tpix,Nfild2);
				dplor=d_punto_lor(cc.x,cc.y,cc.z,lor->n,lor->seno,lor->coseno,lor->phi,lor->z_med,lor->mod_vlor);
				//dplor=d_punto_lor(next_x+(nex*tpix),next_y,next_z,lor->n,lor->seno,lor->coseno,lor->phi,lor->z_med,lor->mod_vlor);
				if(dplor<=lim_exp) {
					set_value_mp(mp,sum,dplor,indice,lim_exp);
					//printf("xpix:%f ypix:%f zpix:%f n:%f phi:%f z_med:%f \n",next_x+(nex*tpix),next_y,next_z,lor->n,lor->phi,lor->z_med);  
				}
				mask_ima[indice]=1;
			}
		}
		nex--;
	}
}

void	expansion_Y(int npix_exp2d, float next_x, float next_y, float next_z, float tpix, float lvox, struct lor_params *lor, int Nfil, int Ntallsima, int jp, int ind_map, int Npix, sparse_i4 *mp, float sum,float lim_exp,int *mask_ima)
{
	int nex,indice,Nfild2,ii,jj,kk;
	float ri,dplor;
	struct coord3d cp,cc;
	
	Nfild2=0.5*Nfil;
	nex=npix_exp2d;
	
	while(nex>=0)
	{
		//printf("\texp_Y-nex:%d: x:%f y:%f z:%f indice:%d\n",nex,next_x,next_y-(nex*tpix),next_z,ind_map-nex);  
		ri=sqrt((next_x)*(next_x)+((next_y-(nex*tpix))*(next_y-(nex*tpix))));
		//printf("\nINICIO: %f<=%f %d>=%d %d<%d\n",ri,(Nfild2+0.5)*tpix,jp,nex,jp,Nfil-nex);
  		if(ri<=(Nfild2+0.5)*tpix && jp>=nex && jp<Nfil-nex) 
		{
			indice=ind_map-nex;
			//printf("\nindice:%d mask:%d  ",indice,mask_ima[indice]);
			if(mask_ima[indice]==0)
			{
				cp.x=next_x;
				cp.y=next_y-(nex*tpix);
				cp.z=next_z;
				ptos_corte_en_ijk_pWu_fast(&ii,&jj,&kk,&cp,lvox,tpix,Nfild2);
				ptos_centro_voxel_pWu_fast(&cc,ii,jj,kk,lvox,tpix,Nfild2);
				dplor=d_punto_lor(cc.x,cc.y,cc.z,lor->n,lor->seno,lor->coseno,lor->phi,lor->z_med,lor->mod_vlor);
				//dplor=d_punto_lor(next_x,next_y-(nex*tpix),next_z,lor->n,lor->seno,lor->coseno,lor->phi,lor->z_med,lor->mod_vlor);
				//printf("  dplor:%f  ",dplor);
				if(dplor<=lim_exp) {
					set_value_mp(mp,sum,dplor,indice,lim_exp);
					//printf("YES !!! xpix:%f ypix:%f zpix:%f n:%f phi:%f z_med:%f",next_x,next_y+(nex*tpix),next_z,lor->n,lor->phi,lor->z_med);  
				}
				mask_ima[indice]=1;
			}
		}
	     //printf("\texp_Y-nex:%d: x:%f y:%f z:%f indice:%d\n",nex,next_x,next_y+(nex*tpix),next_z,ind_map+nex);  
		ri=sqrt(((next_x)*(next_x))+((next_y+(nex*tpix))*(next_y+(nex*tpix))));
  		if(ri<=(Nfild2+0.5)*tpix && jp>=nex && jp<Nfil-nex) 
		{
			indice=ind_map+nex;
			//printf("\nindice:%d mask:%d  ",indice,mask_ima[indice]);
			if(mask_ima[indice]==0)
			{
				cp.x=next_x;
				cp.y=next_y+(nex*tpix);
				cp.z=next_z;
				ptos_corte_en_ijk_pWu_fast(&ii,&jj,&kk,&cp,lvox,tpix,Nfild2);
				ptos_centro_voxel_pWu_fast(&cc,ii,jj,kk,lvox,tpix,Nfild2);
				dplor=d_punto_lor(cc.x,cc.y,cc.z,lor->n,lor->seno,lor->coseno,lor->phi,lor->z_med,lor->mod_vlor);
				//dplor=d_punto_lor(next_x,next_y+(nex*tpix),next_z,lor->n,lor->seno,lor->coseno,lor->phi,lor->z_med,lor->mod_vlor);
				if(dplor<=lim_exp) {
					set_value_mp(mp,sum,dplor,indice,lim_exp);
					//printf("YES!!! xpix:%f ypix:%f zpix:%f n:%f phi:%f z_med:%f",next_x,next_y+(nex*tpix),next_z,lor->n,lor->phi,lor->z_med); 
				}
				mask_ima[indice]=1;
			}
		}
		nex--;
	}

}


/****************************************************************************/

// calcula rayo a va asignando 'pesos' a medida que avanza por el rayo, según la distancia perpendicular a la LOR.
int ray_tracing_lor_pWu(struct coord3d *c1,struct coord3d *c2,sparse_i4 *mp, int Nfil, int Ntallsima, float tpix, float lvox, float tbin, float dmax, int ind_proj, int Nb_os, int fac_comp, float sum,struct lor_params *lor, int npix_exp2d,struct imagen *ima, int desf)
{
float vx,vy,vz,dpb,dant,dx,dy,dz,next_x,next_y,next_z;
int ind_map;
int incrx,incry,incrz;
int cas,nex;
int i,j,k;
float limite; //valor epsilon para comparar valores y dar signos
float xpix,ypix,zpix,xbin,ybin,zbin;
int Nfild2,Npix,out_cyl;
float lim_exp;
int *mask_ima;


limite=0.0001;
lim_exp=2*tbin;

if(npix_exp2d==0) lim_exp=2.5*tbin;

xpix=c1->x;
ypix=c1->y;
zpix=c1->z;
xbin=c2->x;
ybin=c2->y;
zbin=c2->z;
Nfild2=0.5*Nfil;
Npix=Nfil*Nfil;
//printf("\ncoord1:%f,%f,%f coord2:%f,%f,%f",c1->x,c1->y,c1->z,c2->x,c2->y,c2->z);
// control para los puntos sobre las esquinas de la malla
// if(xbin==xpix && ybin==xbin && zbin==zpix) 
// {
// 	xpix+=10.*tpix;
// 	ypix+=10.*tpix;
// 	zpix+=tpix;
// }

//calculo de vectores unitarios en cada dirección
//signo es -1 si vx es negativo, 1 si vx es positivo y 0 si vx es cero
vx=xbin-xpix;
vy=ybin-ypix;
vz=zbin-zpix;
incrx=signo(vx,limite);
incry=signo(vy,limite); 
incrz=signo(vz,limite);
dpb=sqrt(vx*vx+vy*vy+vz*vz); 
vx/=dpb; 
vy/=dpb;
vz/=dpb;

//Indices que corresponden a la coordenada c1 (punto inicial) para ray-tracing
ptos_corte_en_ijk_pWu_fast(&i,&j,&k,c1,lvox,tpix,Nfild2);

if(i==Nfil) i--;
if(j==Nfil) j--;
if(i==-1) i++;
if(j==-1) j++;

next_x=c1->x;        
next_y=c1->y;
next_z=c1->z; 



// Calculo el paso inicial
if(fabs(vx)>EPSILON) dx=(next_x-xpix)/vx;
else dx=dmax;
if(fabs(vy)>EPSILON) dy=(next_y-ypix)/vy;
else dy=dmax;
if(fabs(vz)>EPSILON) dz=(next_z-zpix)/vz;
else dz=dmax;

dant=-1;
ind_map=k*Npix+i*Nfil+j; //Npix son los cortes, Nfil son filas y columnas
//ima->datos[(i*Nfil+j)+(desf)*Npix]+=1.;
out_cyl=0;
//printf("\nEntra?: %d,%d,%d (%f,%f,%f) LOR n:%f\n****************\n\n",i,j,k,next_x,next_y,next_z,lor->n);
//printf("\n\n\nPunto inicial: %d,%d,%d    Punto siguiente: %d,%d,%d (%d)\n\n\n",i,j,k,ip,jp,kp,ind_map);
if((mask_ima = (int *)calloc(Npix*Ntallsima,sizeof(int)))==NULL) error_fesmp(53,Npix*Ntallsima*sizeof(int)/1048576);


//printf("\n***********************\n");
while(i>=0 && i<Nfil && j>=0 && j<Nfil && k>=0 && k<Ntallsima)
{
  //  devuelve el valor "cas" segun valores de "dx-dy-dz"
 compara(&dx,&dy,&dz,&cas);
 //printf("+++");

  switch(cas)
  {
     case 1: // dz es el menor
	 
		  if(dz!=dant)
		  {
			//printf("\n********************\nZ:next: %f %f %f ijk-(%d,%d,%d) \n",next_x,next_y,next_z,i,j,k);  			
			nex=0;
			while(nex<=npix_exp2d && k<Ntallsima-nex && k>=0)
			{
				//printf("NEX pos:%d\n",nex);
				expansion_X(npix_exp2d, next_x, next_y, next_z+(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, i, ind_map+(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				expansion_Y(npix_exp2d, next_x, next_y, next_z+(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, j, ind_map+(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				nex++;
			}
				
			nex=0;
			while(nex<=npix_exp2d && k>=nex && k<Ntallsima)
			{	 	
				//printf("NEX neg:%d\n",nex);
				expansion_X(npix_exp2d, next_x, next_y, next_z-(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, i, ind_map-(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				expansion_Y(npix_exp2d, next_x, next_y, next_z-(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, j, ind_map-(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				nex++;
			}
			
		  }         
		dant=dz;
		k+=incrz;
		ind_map+=incrz*Npix;
          next_z+=lvox*incrz;
          dz=(next_z-zpix)/vz;
		break;
     
	 case 2: // dy es el menor
	 	  
		  if(dy!=dant)
		  { //printf("\n********************\nY:next: %f %f %f ijk (%d,%d,%d)\n",next_x,next_y,next_z,i,j,k);
			nex=0;
			while(nex<=npix_exp2d && k<Ntallsima-nex && k>=0)
			{
				//printf("NEX pos:%d\n",nex);
				expansion_X(npix_exp2d, next_x, next_y, next_z+(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, i, ind_map+(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				expansion_Y(npix_exp2d, next_x, next_y, next_z+(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, j, ind_map+(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				nex++;
			}
				
			nex=0;
			while(nex<=npix_exp2d && k>=nex && k<Ntallsima)
			{	
				//printf("NEX neg:%d\n",nex);
				expansion_X(npix_exp2d, next_x, next_y, next_z-(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, i, ind_map-(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				expansion_Y(npix_exp2d, next_x, next_y, next_z-(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, j, ind_map-(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				nex++;
			}
			
		  }			
		dant=dy;
          j+=incry;
		ind_map+=incry;
          next_y+=tpix*incry;
          dy=(next_y-ypix)/vy;
		break;

     case 3: // dx es el menor
	 	  
		  if(dx!=dant)
		  { //printf("\n********************\nX:next: %f %f %f ijk (%d,%d,%d)\n",next_x,next_y,next_z,i,j,k);
			nex=0;

			while(nex<=npix_exp2d && k<Ntallsima-nex && k>=0)
			{
				//printf("NEX neg:%d\n",nex);
				expansion_X(npix_exp2d, next_x, next_y, next_z+(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, i, ind_map+(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				expansion_Y(npix_exp2d, next_x, next_y, next_z+(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, j, ind_map+(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				nex++;
			}
				
			nex=0;
			//printf("\n");
			while(nex<=npix_exp2d && k>=nex && k<Ntallsima)
			{				
				//printf("NEX neg:%d\n",nex);
				expansion_X(npix_exp2d, next_x, next_y, next_z-(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, i, ind_map-(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				expansion_Y(npix_exp2d, next_x, next_y, next_z-(nex*lvox), tpix, lvox, lor, Nfil, Ntallsima, j, ind_map-(nex*Npix), Npix, mp, sum, lim_exp,mask_ima);
				nex++;
			}
			//printf("\n\n\n");
				
		  }
		dant=dx;
          i+=incrx;
		ind_map+=incrx*Nfil;
          next_x+=tpix*incrx;
          dx=(next_x-xpix)/vx;;
		break;
      
	  default:
          printf("\nError al comparar diferenciales en proyectores_pet.h\n");
	  exit(0);
   }
   out_cyl=0;
}
free(mask_ima);
return(0);
}





















/***************************** 





Calcula los ray_sum con subsample





****************************************************************************/

// calcula rayo a va asignando 'pesos' a medida que avanza por el rayo, en este caso tiene en cuenta los detect planos
int ray_tracing_lor_sub(struct coord3d *c1,struct coord3d *c2,sparse_i4 *mp, int Nfil, int Ntallsima, float tpix, float lvox, float dmax, int ind_proj, int Nb_os, int fac_comp, float sum)
{
float vx,vy,vz,dpb,dant,dx,dy,dz,next_x,next_y,next_z;
int ip,jp,kp,ind_map;
int incrx,incry,incrz;
int cas,n,s;
int i,j,k;
float limite; //valor epsilon para comparar valores y dar signos
float xpix,ypix,zpix,xbin,ybin,zbin,ri;
int Nfild2,Npix,out_cyl;

limite=0.0001;

xpix=c1->x;
ypix=c1->y;
zpix=c1->z;
xbin=c2->x;
ybin=c2->y;
zbin=c2->z;
Nfild2=0.5*Nfil;
Npix=Nfil*Nfil;
//printf("\ncoord1:%f,%f,%f coord2:%f,%f,%f",c1->x,c1->y,c1->z,c2->x,c2->y,c2->z);
// control para los puntos sobre las esquinas de la malla
if(xbin==xpix && ybin==xbin && zbin==zpix) 
{
	xpix+=10.*tpix;
	ypix+=10.*tpix;
	zpix+=tpix;
}

//calculo de vectores unitarios en cada dirección
//signo es -1 si vx es negativo, 1 si vx es positivo y 0 si vx es cero
vx=xbin-xpix;
vy=ybin-ypix;
vz=zbin-zpix;
incrx=signo(vx,limite);
incry=signo(vy,limite);
incrz=signo(vz,limite);
dpb=sqrt(vx*vx+vy*vy+vz*vz); 
vx/=dpb; 
vy/=dpb;
vz/=dpb;


//Indices que corresponden a la coordenada c1 (punto inicial) para ray-tracing
ptos_corte_en_ijk_rt(&i,&j,&k,c1,lvox,tpix,Nfild2);
//printf("\ncoord1:%f,%f,%f indexes:%d,%d,%d",c1->x,c1->y,c1->z,i,j,k);
// Índices del voxel siguiente en la dirección de la LOR
if(incrx<=0) ip=i;
else ip=i+1;
if(incry<=0) jp=j;
else jp=j+1;
if(incrz<=0) kp=k;
else kp=k+1;

next_x=((ip-Nfild2)*tpix)+(0.5*tpix);        
next_y=((jp-Nfild2)*tpix)+(0.5*tpix);
next_z=((kp)*lvox)+(0.5*lvox);   
//next_z=((kp)*lvox);   

// Calculo el paso inicial
if(fabs(vx)>EPSILON) dx=(next_x-xpix)/vx;
else dx=dmax;
if(fabs(vy)>EPSILON) dy=(next_y-ypix)/vy;
else dy=dmax;
if(fabs(vz)>EPSILON) dz=(next_z-zpix)/vz;
else dz=dmax;

dant=0;
ind_map=k*Npix+i*Nfil+j; //Npix son los cortes, Nfil son filas y columnas
if(fac_comp==666) mp->ne+=ind_map;
out_cyl=0;
//printf("\nxpix:%f ypix:%f xbin:%f ybin:%f",xpix,ypix,xbin,ybin);
  
while(ip>=0 && ip<=Nfil && jp>=0 && jp<=Nfil && kp>=0 && kp<=Ntallsima)
{
  
  //  devuelve el valor "cas" segun valores de "dx-dy-dz"
  compara(&dx,&dy,&dz,&cas);
  
  //comprueba que esté dentro del cilindro-FOV
  ri=sqrt(((next_x-(0.5*tpix))*(next_x-(0.5*tpix)))+((next_y-(0.5*tpix))*(next_y-(0.5*tpix))));
  if(ri>=(Nfild2-2.5)*tpix) out_cyl=1;
  
  //imprime valores
  
  //printf("\ni:%d j:%d k:%d -----------> ip:%d jp:%d kp:%d (ind_map:%d)(lor:%d)",i,j,k,ip,jp,kp,ind_map,mp->ia[ind_proj]);
  switch(cas)
  {
     case 1: // dz es el menor
	 
		  if(fac_comp==666)
		  {	 	
			mp->ar[mp->ne]=sum*fabs(dz-dant);
			if(out_cyl==1) mp->ar[mp->ne]=0.;
			mp->ja[mp->ne]=ind_map;
			if(incrz>=0) mp->ne+=Npix;
			if(incrz<0) mp->ne-=Npix;
		  }
		  
		  if(fac_comp!=666 && out_cyl==0 && dz!=dant)
		  {
		    n=mp->ne-1;
			s=-1;
		  	while(n>=mp->ia[ind_proj])
			{
				if(mp->ja[n]==ind_map) s=n;
				n--;
			}
			
			if(s!=-1)
			{
				mp->ar[s]+=sum*fabs(dz-dant);
			}
			
			if(s==-1)
			{
				mp->ar[mp->ne]=sum*fabs(dz-dant);
				mp->ja[mp->ne]=ind_map;
				mp->ne++;
			}
			
		   }
		            
		  dant=dz;
          kp+=incrz;
		  ind_map+=incrz*Npix;
          next_z+=lvox*incrz;
          dz=(next_z-zpix)/vz;
		  break;
     
	 case 2: // dy es el menor
	 
	 	  if(fac_comp==666)
		  {
			mp->ar[mp->ne]=sum*fabs(dy-dant);
			if(out_cyl==1) mp->ar[mp->ne]=0.;
			mp->ja[mp->ne]=ind_map;
			if(incry>=0) mp->ne++;
			if(incry<0) mp->ne--;	
		  }
		  
		  if(fac_comp!=666 && out_cyl==0 && dy!=dant)
		  {
		    n=mp->ne-1;
			s=-1;
		  	while(n>=mp->ia[ind_proj])
			{
				if(ind_map==mp->ja[n]) s=n;
				n--;
			}
			
			if(s!=-1)
			{
				mp->ar[s]+=sum*fabs(dy-dant);
			}
			
			if(s==-1)
			{
				mp->ar[mp->ne]=sum*fabs(dy-dant);
				mp->ja[mp->ne]=ind_map;
				mp->ne++;
			}
			
		   }
		  
		  dant=dy;
          jp+=incry;
          ind_map+=incry;
          next_y+=tpix*incry;
          dy=(next_y-ypix)/vy;
		  break;

     case 3: // dx es el menor
	 	  
		  if(fac_comp==666)
		  {
			mp->ar[mp->ne]=sum*fabs(dx-dant);
			if(out_cyl==1) mp->ar[mp->ne]=0.;
			mp->ja[mp->ne]=ind_map;	  
			if(incrx>=0) mp->ne+=Nfil;
			if(incrx<0) mp->ne-=Nfil;          
		  }
		  
		  if(fac_comp!=666 && out_cyl==0 && dx!=dant)
		  {
		    n=mp->ne-1;
			s=-1;
		  	while(n>=mp->ia[ind_proj])
			{
				//printf("\nBuscando: n:%d %d vs %d",n,ind_map,mp->ja[n]);
				
				if(ind_map==mp->ja[n]) s=n;
				n--;
			}
			
			if(s!=-1)
			{
				//printf("\tsuma en %d",ind_map);
				mp->ar[s]+=sum*fabs(dx-dant);
			}
			
			if(s==-1)
			{
				//printf("\tnuevo");
				mp->ar[mp->ne]=sum*fabs(dx-dant);
				mp->ja[mp->ne]=ind_map;
				mp->ne++;
			}
			
		   }
		  
		  dant=dx;
          ip+=incrx;
          ind_map+=incrx*Nfil;
          next_x+=tpix*incrx;
          dx=(next_x-xpix)/vx;
		  break;
      
	  default:
          printf("\nError al comparar diferenciales en proyectores_pet.h\n");
	  exit(0);
   }
   
   out_cyl=0;

}

return(0);
}












/***************************** 






Calcula puntos de corte 





*****************************************************************************************************/

void calcul_ptos_corte(float a,float n,float seno,float coseno,float tan_phi,float z_med,struct coord3d *c1,struct coord3d *c2)
{
	float *x,*y,*z,lambda;
	float *modulos;
	int i;
	modulos=vector(0,3);  // funcion de util/nrutil.h de Fruitcake
	x=vector(0,3);
	y=vector(0,3);
	z=vector(0,3);
	/*
	Origen de angulos en eje X, si coseno es 0 entonces angulo 90 direccion del eje y
	*/

	//Calculo de puntos de corte para plano perpendicular al eje x positivo
	if(coseno<=EPSILON && coseno>=-EPSILON){ //direccion del eje y
		x[0]=1000*a;
		y[0]=1000*a;
		z[0]=1000*a;//simula pto corte en infinito(para que no se incluya al comparar)
		}
		else{
		lambda=(a+(n*seno))/coseno;
		x[0]=a;
		y[0]=(n*coseno)+(lambda*seno);
		z[0]=z_med+(lambda*tan_phi);
		}
		modulos[0]=modulo_vector(x[0],y[0],0);
	
	//Calculo de puntos de corte para plano perpendicular al eje y positivo
	if(seno<=EPSILON && seno>=-EPSILON){ // direccion del eje x
		x[1]=1000*a;
		y[1]=1000*a;
		z[1]=1000*a;//simula pto corte en infinito(para que no se incluya al comparar)
		}
		else{
		lambda=(a-(n*coseno))/seno;
		x[1]=(-n*seno)+(lambda*coseno);
		y[1]=a;
		z[1]=z_med+(lambda*tan_phi);
		}
		modulos[1]=modulo_vector(x[1],y[1],0);
	
	//Calculo de puntos de corte para plano perpendicular al eje x negativo
	if(coseno<=EPSILON && coseno>=-EPSILON){ //direccion del eje y
		x[2]=1000*a;
		y[2]=1000*a;
		z[2]=1000*a;//simula pto corte en infinito(para que no se incluya al comparar)
		}
		else{
		lambda=(-a+(n*seno))/coseno;
		x[2]=-a;
		y[2]=(n*coseno)+(lambda*seno);
		z[2]=z_med+(lambda*tan_phi);
		}
		modulos[2]=modulo_vector(x[2],y[2],0);
	
	//Calculo de puntos de corte para plano perpendicular al eje y negativo
	if(seno<=EPSILON && seno>=-EPSILON){ //direccion del eje x
		x[3]=1000*a;
		y[3]=1000*a;
		z[3]=1000*a;//simula pto corte en infinito(para que no se incluya al comparar)
		}
		else{
		lambda=(-a-(n*coseno))/seno;
		x[3]=(-n*seno)+(lambda*coseno);
		y[3]=-a;
		z[3]=z_med+(lambda*tan_phi);
		}
		modulos[3]=modulo_vector(x[3],y[3],0);
	//printf("----------------------------------------------------\n");	
	
	// pone las componentes fuera del cuadrado con valores enormes para que se desprecien en la ordenacion
	for(i=0;i<=3;i++){
	if(x[i]>a || x[i]<-a) x[i]=1000*x[i];
	if(y[i]>a || y[i]<-a) y[i]=1000*y[i];
	}
	// escoge las coordenadas de los dos modulos mas pequeños
	//Se ordenan los vectores x,y,z segun los valores del vector modulos
	ordena_cuatro_vectores(4,modulos,x,y,z);	

	free_vector(modulos,0,3);

	c1->x=x[0];
	c1->y=y[0];
	c1->z=z[0];
	c2->x=x[1];
	c2->y=y[1];
	c2->z=z[1];

	free_vector(x,0,3);
	free_vector(y,0,3);
	free_vector(z,0,3);
}

/************Transforma los ptos de corte calculados con la funcion anterior en los índices del voxel más próximo***************/
void ptos_corte_en_ijk(int *i,int *j,int *k,struct coord3d *c,float lvox,float tpix,int Nfild2)
{
float limite;
limite=0.0001;
//Asignacion al indice k del voxel
*k=(c->z/lvox);

//Asignacion a los índices i y j del voxel
*i=(c->x/tpix)+Nfild2;
*j=(c->y/tpix)+Nfild2;
}


/************Transforma los ptos de corte calculados con la funcion anterior en los índices del voxel más próximo***************/
void ptos_corte_en_ijk_rt(int *i,int *j,int *k,struct coord3d *c,float lvox,float tpix,int Nfild2)
{
float limite;
limite=0.0001;
//Asignacion al indice k del voxel
*k=(c->z-(0.49*lvox))/lvox;
//*k=(c->z)/lvox;
//Asignacion a los índices i y j del voxel
*i=((c->x-(0.5*tpix))/tpix)+Nfild2;
*j=((c->y-(0.5*tpix))/tpix)+Nfild2;
}

/************Transforma los ptos de corte calculados con la funcion anterior en los índices del voxel más próximo***************/
void ptos_corte_en_ijk_pWu_fast(int *i,int *j,int *k,struct coord3d *c,float lvox,float tpix,int Nfild2)
{
float limite;
limite=0.0001;

// *k=(int)(c->z/lvox-0.5);
// *i=(int)(c->x/tpix+0.5)+Nfild2;
// *j=(int)(c->y/tpix+0.5)+Nfild2;

*k=rint((c->z-lvox)/lvox);
*i=rint((c->x/tpix)-0.5+Nfild2);
*j=rint((c->y/tpix)-0.5+Nfild2);
//printf("i %d,j %d, k %d\n",*i,*j,*k);
}

/************Transforma los ptos de corte calculados con la funcion anterior en los índices del voxel más próximo***************/
void ptos_centro_voxel_pWu_fast(struct coord3d *cc,int i,int j,int k,float lvox,float tpix,int Nfild2)
{
// cc->z=(k*lvox)+lvox;
// cc->x=(i*tpix)+0.5-(Nfild2*tpix);
// cc->y=(j*tpix)+0.5-(Nfild2*tpix);
cc->z=(k*lvox)+lvox;
cc->x=(i+0.5-Nfild2)*tpix;
cc->y=(j+0.5-Nfild2)*tpix;

}










// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////





// Función que obtiene la matriz del sistema rPET, Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 
// Preparada para proyectar factores de Atenuación, para ello se elimina el efecto del factor SUM, aplicado para
// tener en cuenta que una LOR puede aparecer en varias posiciones del detector, lo cual no es necesario incluirlo para
// calcular los factores de atenuación 3D.

// Modificado para ATT en Mayo2007





// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int build_srm_rpet_ray_tracing_att(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset)
{
char nom_fitxer[128],nom_fitxer_sb[128],nom_fitxer2[128];
int in,ind_proj,l,k2,k1,Ncol,nsb,Nb_os;
double NITEMS_INI,NITEMS;
float z1,z2,z_med,n,angle;
float n0,angle0,z0,inc_zp,inc_ang,inc_n,d2,sum;
float coseno,seno,phi,factor,mod_vlor,zmin,zmax,tan_phi;
float epsilon=1e-4,dmax;
sparse_i4 mp;
struct coord3d coord1;
struct coord3d coord2;

// params definition
strcpy(nom_fitxer,"sresmx");
strcpy(nom_fitxer2,"sresmx");
sprintf(nom_fitxer2,"%s.pesos",nom_fitxer);
Ncol=Nfil;
mp.n_fil=Nfil; // image pixels in X-Y-directions 
mp.n_col=Ncol; // image pixels in X-Y-directions 
mp.Nbins=(2*ndet)-1; // transaxial detectors 
mp.tbin=tdet/2.;
mp.tpix=mp.Nbins*mp.tbin/Nfil; // image pixel size in X-Y-directions in mm
mp.Ntallsima=2*ndet_z-1; // image pixel in Z-direction
if(axial_diff==0) mp.Ntallsima=ndet_z;  
mp.lvox=tdet_z/2.; // image pixel size in Z-direction in mm 
if(axial_diff==0) mp.lvox=tdet_z;
mp.Nang=Nang;
mp.ndet_z=ndet_z;  
mp.tdet_z=tdet_z;
mp.ddet=ddet;
mp.Nsubsets=Nsubsets;
mp.axial_diff=axial_diff;
mp.inv=0;
Nbt=mp.Nbins*Nang*num_imagenes(ndet_z,axial_diff);	// sinogram elements 
Nb_os=Nbt/Nsubsets; // sinogram elements per subset 
Nvox=Nfil*Ncol*mp.Ntallsima; // image elements 
d2=0.5*Nfil*mp.tpix; //medio FOV
if(fac_comp!=666) NITEMS_INI=(Nb_os/fac_comp)*Nvox; // matrix elements/fac_comp  
if(fac_comp==666) NITEMS_INI=Nb_os*Nvox; // no sparse matrix elements  
NITEMS=NITEMS_INI; // sparse matrix elements estimation  


// no sparse warning
if(fac_comp==666) printf("\n**\nWARNING! No sparse format was selected. You could have memory problems...\n**\n"); 


// Origin coords 
z0=tdet_z/2.; // mm 
angle0=0.; // degrees 
n0=((-(float)mp.Nbins/2.)+0.5)*mp.tbin; // mm



// steps
inc_zp=tdet_z; // axial increment in mm
inc_ang=180./Nang; // rotating in 180 degrees in degrees
inc_n=mp.tbin; // transaxial increment in mm 


// Display for selected params
printf("\n**************************************************************");
printf("\nVolumen de reconstrucción: %dx%dx%d (pixel:%f mm y corte:%f mm)\n",Nfil,Ncol,mp.Ntallsima,mp.tpix,mp.lvox);
printf("\nCaracterísticas del escáner: Detectores planos: %d bins (tamaño:%f mm) x %d ángulos\n",mp.Nbins,mp.tbin,Nang);
printf("Maximum axial difference is %d (Sinogramas totales: %d)\n",axial_diff,num_imagenes(ndet_z,axial_diff));
printf("%d detectores axiales de %f mm. Detectores separados: %f mm\n\n\n",ndet_z,tdet_z,ddet);
printf("\n**************************************************************");


// memory allocation
printf("\n\nInitial memory estimation to allocate per subset %1.0f Mb\n\n",NITEMS*sizeof(float)/1048576.);
if((mp.ar = (float *)calloc(NITEMS,sizeof(float)))==NULL) error_fesmp(53,NITEMS*sizeof(float)/1048576);
if((mp.ja = (int *)calloc(NITEMS,sizeof(int)))==NULL) error_fesmp(53,NITEMS*sizeof(int)/1048576);
mp.ia = (int *)calloc(Nbt+1,sizeof(int));


// SUBSET control  
for(nsb=0;nsb<Nsubsets;nsb++)
{
if(opc_subset==nsb || opc_subset==-1)
	{
		//inicializa mp y genera nombre matriz correspondiente al subset
		mp.ne=0;
		sprintf(nom_fitxer_sb,"%s.%d",nom_fitxer,nsb);
		printf("\n************\n\nSubset: %d\n",nsb);
			
		// AXIAL LOR control
		for(k1=0,z1=z0,ind_proj=0;k1<ndet_z;z1+=inc_zp,k1++)
		{
		for(k2=0,z2=z0;k2<ndet_z;z2+=inc_zp,k2++)
		{
			
			if(abs(k1-k2)>axial_diff)
			{
				continue;
			}
		
			// n-angular independent LOR params definition
			z_med=fabs(z1+z2)/2;  
			zmax=maximo(z1,z2)+epsilon;
			zmin=minimo(z1,z2)-epsilon;
			printf("\nSinogram: z1=%d  z2=%d (memory used: %f Mb)",k1,k2,6.*mp.ne/1048576.);
			phi=atan((z2-z1)/ddet); 
			tan_phi=tan(phi);
			factor=cos(phi);
			mod_vlor=sqrt((fabs(z2-z1)*fabs(z2-z1))+(ddet*ddet));
			
			// ANGULAR LOR control
			for(l=nsb,angle=angle0+(nsb*inc_ang);l<Nang;l+=Nsubsets,angle+=inc_ang*Nsubsets)
			{
				// angular LOR params definition
				seno=sin(angle*M_PI/180.);
				coseno=cos(angle*M_PI/180.);
				
				// TRANSAXIAL LOR control
				for(in=0,n=n0;in<mp.Nbins;in++,ind_proj++,n+=inc_n)
				{
				
					// El factor sum ahora es siempre 1 (Mayo 2007 Pablo Aguiar)
					sum=1.;
					
					//calculate non-null elements only
					mp.ia[ind_proj]=mp.ne; 
					calcul_ptos_corte(d2, n, seno, coseno, tan_phi,z_med, &coord1, &coord2);
					dmax=10000*maximo(fabs(coord1.x-coord2.x),fabs(coord1.y-coord2.y));
					ray_tracing_lor(&coord1,&coord2, &mp, Nfil, mp.Ntallsima, mp.tpix, mp.lvox, dmax, ind_proj, Nb_os, fac_comp, sum);
	
					// redimensiona memoria si es necesario
					if(mp.ne >= NITEMS-Nvox)
					{
						NITEMS+=NITEMS_INI;
						if((mp.ar = (float *)realloc(mp.ar,NITEMS*sizeof(float)))==NULL) 
						{ 
							error_fesmp(53,NITEMS*sizeof(float)/1048576.);
						}
						if((mp.ja = (int *)realloc(mp.ja,NITEMS*sizeof(int)))==NULL)
						{
							error_fesmp(53,NITEMS*sizeof(int)/1048576.);
						}
						printf("\nrealloc for %1.0f Mb\n",NITEMS*sizeof(float)/1048576.);
					}
						
				} // transaxial LOR control
			} // angular LOR control
		} // axial LOR control
		} // axial LOR control
	
	
		// Writting matrix for each subset
		mp.ia[Nb_os]=mp.ne; printf("\n\n\nNo-null elements: %d\n",mp.ne);

		// for sparse format
		if(fac_comp!=666)
		{ 
			guarda_pesos_xOS(&mp,nom_fitxer_sb);
		}
		

	}
}
return 0;
}









// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////





// Función que obtiene la matriz del sistema rPET, Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 





// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int build_srm_rpet_ray_tracing(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset)
{
FILE *mat;
char nom_fitxer[128],nom_fitxer_sb[128],nom_fitxer_sb2[128],nom_fitxer2[128];
int in,ind_proj,l,k2,k1,Ncol,nsb,Nb_os;
double NITEMS_INI,NITEMS;
float z1,z2,z_med,n,angle;
float n0,n_i,angle0,z0,inc_zp,inc_ang,ang_i,inc_n,d2,sum;
float coseno,seno,phi,factor,mod_vlor,zmin,zmax,tan_phi;
float epsilon=1e-4,dmax;
sparse_i4 mp;
struct coord3d coord1;
struct coord3d coord2;

// params definition
strcpy(nom_fitxer,"sresmx");
strcpy(nom_fitxer2,"sresmx");
sprintf(nom_fitxer2,"%s.pesos",nom_fitxer);
Ncol=Nfil;
mp.n_fil=Nfil; // image pixels in X-Y-directions 
mp.n_col=Ncol; // image pixels in X-Y-directions 
mp.Nbins=(2*ndet)-1; // transaxial detectors 
mp.tbin=tdet/2.;
mp.tpix=mp.Nbins*mp.tbin/Nfil; // image pixel size in X-Y-directions in mm
mp.Ntallsima=2*ndet_z-1; // image pixel in Z-direction
if(axial_diff==0) mp.Ntallsima=ndet_z;  
mp.lvox=tdet_z/2.; // image pixel size in Z-direction in mm 
if(axial_diff==0) mp.lvox=tdet_z;
mp.Nang=Nang;
mp.ndet_z=ndet_z;  
mp.tdet_z=tdet_z;
mp.ddet=ddet;
mp.Nsubsets=Nsubsets;
mp.axial_diff=axial_diff;
mp.inv=0;
Nbt=mp.Nbins*Nang*num_imagenes(ndet_z,axial_diff);	// sinogram elements 
Nb_os=Nbt/Nsubsets; // sinogram elements per subset 
Nvox=Nfil*Ncol*mp.Ntallsima; // image elements 
d2=0.5*Nfil*mp.tpix; //medio FOV
if(fac_comp!=666) NITEMS_INI=(Nb_os/fac_comp)*Nvox; // matrix elements/fac_comp  
if(fac_comp==666) NITEMS_INI=(Nb_os/2)*Nvox; // no sparse matrix elements  
NITEMS=NITEMS_INI; // sparse matrix elements estimation  


// no sparse warning
if(fac_comp==666) printf("\n**\nWARNING! No sparse format was selected. You could have memory problems...\n**\n"); 


// Origin coords 
z0=tdet_z/2.; // mm 
angle0=0.; // degrees 
n0=((-(float)mp.Nbins/2.)+0.5)*mp.tbin; // mm



// steps
inc_zp=tdet_z; // axial increment in mm
inc_ang=180./Nang; // rotating in 180 degrees in degrees
inc_n=mp.tbin; // transaxial increment in mm 


// Display for selected params
printf("\n**************************************************************");
printf("\nVolumen de reconstrucción: %dx%dx%d (pixel:%f mm y corte:%f mm)\n",Nfil,Ncol,mp.Ntallsima,mp.tpix,mp.lvox);
printf("\nCaracterísticas del escáner: Detectores planos: %d bins (tamaño:%f mm) x %d ángulos\n",mp.Nbins,mp.tbin,Nang);
printf("Maximum axial difference is %d (Sinogramas totales: %d)\n",axial_diff,num_imagenes(ndet_z,axial_diff));
printf("%d detectores axiales de %f mm. Detectores separados: %f mm\n\n\n",ndet_z,tdet_z,ddet);
printf("\n**************************************************************");


// memory allocation
printf("\n\nInitial memory estimation to allocate per subset %1.0f Mb\n\n",NITEMS*sizeof(float)/1048576.);
if((mp.ar = (float *)calloc(NITEMS,sizeof(float)))==NULL) error_fesmp(53,NITEMS*sizeof(float)/1048576);
if((mp.ja = (int *)calloc(NITEMS,sizeof(int)))==NULL) error_fesmp(53,NITEMS*sizeof(int)/1048576);
mp.ia = (int *)calloc(Nbt+1,sizeof(int));


// SUBSET control  
for(nsb=0;nsb<Nsubsets;nsb++)
{
if(opc_subset==nsb || opc_subset==-1)
	{
		//inicializa mp y genera nombre matriz correspondiente al subset
		mp.ne=0;
		sprintf(nom_fitxer_sb,"%s.%d",nom_fitxer,nsb);
		printf("\n************\n\nSubset: %d\n",nsb);
			
		// AXIAL LOR control
		for(k1=0,z1=z0,ind_proj=0;k1<ndet_z;z1+=inc_zp,k1++)
		{
		for(k2=0,z2=z0;k2<ndet_z;z2+=inc_zp,k2++)
		{
			
			if(abs(k1-k2)>axial_diff)
      //   if(k1>14 || k1<16 || k2>14 || k2<16)
			{
				continue;
			}
		
			// n-angular independent LOR params definition
			z_med=fabs(z1+z2)/2;  
			zmax=maximo(z1,z2)+epsilon;
			zmin=minimo(z1,z2)-epsilon;
			printf("\nSinogram: z1=%d  z2=%d (memory used: %f Mb)",k1,k2,6.*mp.ne/1048576.);
			phi=atan((z2-z1)/ddet); 
			tan_phi=tan(phi);
			factor=cos(phi);
			mod_vlor=sqrt((fabs(z2-z1)*fabs(z2-z1))+(ddet*ddet));
			
			// ANGULAR LOR control
         for(l=nsb,angle=angle0+(nsb*inc_ang);l<Nang;l+=Nsubsets,angle+=inc_ang*Nsubsets)						
         {
				// angular LOR params definition
				seno=sin(angle*M_PI/180.);
				coseno=cos(angle*M_PI/180.);
				
				// TRANSAXIAL LOR control
				for(in=0,n=n0;in<mp.Nbins;in++,ind_proj++,n+=inc_n)
				{
				
					// calcula contribución de LOR en otras posiciones de los detectores
					n_i=n;
					sum=0.;
					while(n_i<=0 && n_i>=n0) 
					{
						ang_i=atan((n_i-n)/0.5*ddet);
						n_i-=inc_n;
						sum+=cos(ang_i*M_PI/180.);
					}
					n_i=n;
					while(n_i>0 && n_i<=-n0) 
					{
						ang_i=atan((n_i-n)/0.5*ddet);
						n_i+=inc_n;
						sum+=cos(ang_i*M_PI/180.);
					}
					sum=1.;
					
					//calculate non-null elements only
					mp.ia[ind_proj]=mp.ne; 
					calcul_ptos_corte(d2, n, seno, coseno, tan_phi,z_med, &coord1, &coord2);
					dmax=10000*maximo(fabs(coord1.x-coord2.x),fabs(coord1.y-coord2.y));
					ray_tracing_lor(&coord1,&coord2, &mp, Nfil, mp.Ntallsima, mp.tpix, mp.lvox, dmax, ind_proj, Nb_os, fac_comp, sum);
	
					// redimensiona memoria si es necesario
					if(mp.ne >= NITEMS-Nvox)
					{
						NITEMS+=NITEMS_INI;
						if((mp.ar = (float *)realloc(mp.ar,NITEMS*sizeof(float)))==NULL) 
						{ 
							error_fesmp(53,NITEMS*sizeof(float)/1048576.);
						}
						if((mp.ja = (int *)realloc(mp.ja,NITEMS*sizeof(int)))==NULL)
						{
							error_fesmp(53,NITEMS*sizeof(int)/1048576.);
						}
						printf("\nrealloc for %1.0f Mb\n",NITEMS*sizeof(float)/1048576.);
					}
						
				} // transaxial LOR control
			} // angular LOR control
		} // axial LOR control
		} // axial LOR control
	
	
		// Writting matrix for each subset
		mp.ia[Nb_os]=mp.ne; printf("\n\n\nNo-null elements: %d\n",mp.ne);

		// for sparse format
		if(fac_comp!=666)
		{ 
			guarda_pesos_xOS(&mp,nom_fitxer_sb);
		}
		
		// for no sparse format
		if(fac_comp==666)
		{
			printf("\n\nOutput: No-Sparse format\n");
			sprintf(nom_fitxer_sb2,"%s.no_sparse",nom_fitxer_sb);
			if(( mat = fopen(nom_fitxer_sb2,"wb"))==NULL) error_fesmp(30,0);
			fwrite (mp.ar,sizeof(float),(size_t)mp.ne,mat);
			fclose (mat);	
		}
	}
}
return 0;
}













// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////




// Función que obtiene la matriz del sistema rPET, Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 
// modificado para sub-sample bin (4points)




// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int build_srm_rpet_ray_tracing_sub(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset)
{
FILE *mat;
char nom_fitxer[128],nom_fitxer_sb[128],nom_fitxer_sb2[128],nom_fitxer2[128];
int in,ind_proj,l,k2,k1,Ncol,nsb,Nb_os;
double NITEMS_INI,NITEMS;
float z1,z2,z_med,n,angle;
float n0,n_i,angle0,z0,inc_zp,inc_ang,ang_i,inc_n,d2,sum;
float coseno,seno,phi,factor,mod_vlor,zmin,zmax,tan_phi;
float epsilon=1e-4,dmax;
float nss,nss0,z1ss,z2ss,z1ss0,z2ss0,shift_ang,angle_ss,ass;
int kn,ka,inv;
sparse_i4 mp;
struct coord3d coord1;
struct coord3d coord2;
int a;


// params definition
strcpy(nom_fitxer,"sresmx");
strcpy(nom_fitxer2,"sresmx");
sprintf(nom_fitxer2,"%s.pesos",nom_fitxer);
Ncol=Nfil;
mp.n_fil=Nfil; // image pixels in X-Y-directions 
mp.n_col=Ncol; // image pixels in X-Y-directions 
mp.Nbins=(2*ndet)-1; // transaxial detectors 
mp.tbin=tdet/2.;
mp.tpix=mp.Nbins*mp.tbin/Nfil; // image pixel size in X-Y-directions in mm
mp.Ntallsima=2*ndet_z-1; // image pixel in Z-direction
if(axial_diff==0) mp.Ntallsima=ndet_z;  
mp.lvox=tdet_z/2.; // image pixel size in Z-direction in mm 
if(axial_diff==0) mp.lvox=tdet_z;
mp.Nang=Nang;
mp.ndet_z=ndet_z;  
mp.tdet_z=tdet_z;
mp.ddet=ddet;
mp.Nsubsets=Nsubsets;
mp.axial_diff=axial_diff;
mp.inv=0;
Nbt=mp.Nbins*Nang*num_imagenes(ndet_z,axial_diff);	// sinogram elements 
Nb_os=Nbt/Nsubsets; // sinogram elements per subset 
Nvox=Nfil*Ncol*mp.Ntallsima; // image elements 
d2=0.5*Nfil*mp.tpix; //medio FOV
if(fac_comp!=666) NITEMS_INI=(Nb_os/fac_comp)*Nvox; // matrix elements/fac_comp  
if(fac_comp==666) NITEMS_INI=Nb_os*Nvox; // no sparse matrix elements  
NITEMS=NITEMS_INI; // sparse matrix elements estimation

// no sparse warning
if(fac_comp==666) printf("\n**\nWARNING! No sparse format was selected. You could have memory problems...\n**\n"); 


// Origin coords 
z0=tdet_z/2.; // mm 
angle0=0.; // degrees 
n0=((-(float)mp.Nbins/2.)+0.5)*mp.tbin; // mm



// steps
inc_zp=tdet_z; // axial increment in mm
inc_ang=180./Nang; // rotating in 180 degrees in degrees
inc_n=mp.tbin; // transaxial increment in mm 
shift_ang=atan(0.5*inc_n/ddet); // shift angle for sub-sample

// Display for selected params
printf("\n**************************************************************");
printf("\nVolumen de reconstrucción: %dx%dx%d (pixel:%f mm y corte:%f mm)\n",Nfil,Ncol,mp.Ntallsima,mp.tpix,mp.lvox);
printf("\nCaracterísticas del escáner: Detectores planos: %d bins (tamaño:%f mm) x %d ángulos\n",mp.Nbins,mp.tbin,Nang);
printf("Maximum axial difference is %d (Sinogramas totales: %d)\n",axial_diff,num_imagenes(ndet_z,axial_diff));
printf("%d detectores axiales de %f mm. Detectores separados: %f mm\n\n\n",ndet_z,tdet_z,ddet);
printf("\n**************************************************************");


// memory allocation
printf("\n\nInitial memory estimation to allocate per subset %1.0f Mb\n\n",NITEMS*sizeof(float)/1048576.);
if((mp.ar = (float *)calloc(NITEMS,sizeof(float)))==NULL) error_fesmp(53,NITEMS*sizeof(float)/1048576);
if((mp.ja = (int *)calloc(NITEMS,sizeof(int)))==NULL) error_fesmp(53,NITEMS*sizeof(int)/1048576);
mp.ia = (int *)calloc(Nbt+1,sizeof(int));


// SUBSET control  
for(nsb=0;nsb<Nsubsets;nsb++)
{
if(opc_subset==nsb || opc_subset==-1)
	{
		//inicializa mp y genera nombre matriz correspondiente al subset
		mp.ne=0;
		sprintf(nom_fitxer_sb,"%s.%d",nom_fitxer,nsb);
		printf("\n************\n\nSubset: %d\n",nsb);
			
		// AXIAL LOR control
		for(k1=0,z1=z0,ind_proj=0;k1<ndet_z;z1+=inc_zp,k1++)
		{
		for(k2=0,z2=z0;k2<ndet_z;z2+=inc_zp,k2++)
		{
			
			if(abs(k1-k2)>axial_diff)
			{
				continue;
			}
		
			// n-angular independent LOR params definition
			printf("\nSinogram: z1=%d  z2=%d (memory used: %f Mb)",k1,k2,6.*mp.ne/1048576.);
			
			// ANGULAR LOR control
			for(l=nsb,angle=angle0+(nsb*inc_ang);l<Nang;l+=Nsubsets,angle+=inc_ang*Nsubsets)
			{
				
				// TRANSAXIAL LOR control
				for(in=0,n=n0;in<mp.Nbins;in++,ind_proj++,n+=inc_n)
				{
					// calcula contribución de LOR en otras posiciones de los detectores
					n_i=n;
					sum=0.;
					while(n_i<=0 && n_i>=n0) 
					{
						ang_i=atan((n_i-n)/0.5*ddet);
						n_i-=inc_n;
						sum+=cos(ang_i*M_PI/180.);
					}
					n_i=n;
					while(n_i>0 && n_i<=-n0) 
					{
						ang_i=atan((n_i-n)/0.5*ddet);
						n_i+=inc_n;
						sum+=cos(ang_i*M_PI/180.);
					}
					
					//bucles para z1,z2 y n que no aumentan el valor de ind_proj
					// es decir formarán parte de la misma LOR
					nss0=n-(inc_n/4.);
					z1ss0=z1-inc_zp/4.;
					z2ss0=z2-inc_zp/4.	;
					//printf("**************\n: z1:%f,z2:%f,n:%f\n**********************\n",z1,z2,n);
					
					z1ss=z1;
					z2ss=z2;
					nss=n;
					a=0;
					mp.ia[ind_proj]=mp.ne;
					//printf("NEW Subsampl\n");					
					//printf("*********\nSubsample z1:%f (entre %f y %f)\n",z1ss,z1ss0,z1+inc_zp/4);
  					for(z1ss=z1ss0;z1ss<=z1+0.001+inc_zp/4.;z1ss+=inc_zp/2.)
  					{
  					for(z2ss=z2ss0;z2ss<=z2+0.001+inc_zp/4.;z2ss+=inc_zp/2.)
  					{
 					for(kn=0,nss=nss0;kn<2;nss+=inc_n/2.,kn++)
					{
 					for(ka=0,ass=0.;ka<2;ass+=shift_ang,ka++)
					{
					if(kn==1) ass=ass*(-1);
					angle_ss=angle+(ass*180./M_PI);
					seno=sin(angle_ss*M_PI/180.);
					coseno=cos(angle_ss*M_PI/180.);
					//printf("\nangle_ss:%f",angle_ss);
						z_med=fabs(z1ss+z2ss)/2;  
						zmax=maximo(z1ss,z2ss)+epsilon;
						zmin=minimo(z1ss,z2ss)-epsilon;
						phi=atan((z2ss-z1ss)/ddet); 
						tan_phi=tan(phi);
						factor=cos(phi);
						mod_vlor=sqrt((fabs(z2ss-z1ss)*fabs(z2ss-z1ss))+(ddet*ddet));
						
						//calculate non-null elements only 
						calcul_ptos_corte(d2, nss, seno, coseno, tan_phi,z_med, &coord1, &coord2);
						dmax=10000*maximo(fabs(coord1.x-coord2.x),fabs(coord1.y-coord2.y));
						ray_tracing_lor_sub(&coord1,&coord2, &mp, Nfil, mp.Ntallsima, mp.tpix, mp.lvox, dmax, ind_proj, Nb_os, fac_comp, sum);
		
						// redimensiona memoria si es necesario
						if(mp.ne >= NITEMS-Nvox)
						{
							NITEMS+=NITEMS_INI;
							if((mp.ar = (float *)realloc(mp.ar,NITEMS*sizeof(float)))==NULL) 
							{ 
								error_fesmp(53,NITEMS*sizeof(float)/1048576.);
							}
							if((mp.ja = (int *)realloc(mp.ja,NITEMS*sizeof(int)))==NULL)
							{
								error_fesmp(53,NITEMS*sizeof(int)/1048576.);
							}
							printf("\nrealloc for %1.0f Mb\n",NITEMS*sizeof(float)/1048576.);
						}
					}
					}
					}
					}
				} // transaxial LOR control
			} // angular LOR control
		} // axial LOR control
		} // axial LOR control
	
	
		// Writting matrix for each subset
		mp.ia[Nb_os]=mp.ne; printf("\n\n\nNo-null elements: %d\n",mp.ne);
		inv=0;
		// for sparse format
		if(fac_comp!=666)
		{
			guarda_pesos_xOS(&mp,nom_fitxer_sb);
		}
		
		// for no sparse format
		if(fac_comp==666)
		{
			printf("\n\nOutput: No-Sparse format\n");
			sprintf(nom_fitxer_sb2,"%s.no_sparse",nom_fitxer_sb);
			if(( mat = fopen(nom_fitxer_sb2,"wb"))==NULL) error_fesmp(30,0);
			fwrite (mp.ar,sizeof(float),(size_t)mp.ne,mat);
			fclose (mat);	
		}
	}
}
return 0;
}














// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////




// Función que obtiene la matriz del sistema rPET, Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 
// modificado para sub-sample bin sin cruzados 



// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int build_srm_rpet_ray_tracing_sub_simple(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset)
{
FILE *mat;
char nom_fitxer[128],nom_fitxer_sb[128],nom_fitxer_sb2[128],nom_fitxer2[128];
int in,ind_proj,l,k2,k1,Ncol,nsb,Nb_os;
double NITEMS_INI,NITEMS;
float z1,z2,z_med,n,angle;
float n0,n_i,angle0,z0,inc_zp,inc_ang,ang_i,inc_n,d2,sum;
float coseno,seno,phi,factor,mod_vlor,zmin,zmax,tan_phi;
float epsilon=1e-4,dmax;
float nss,nss0,z1ss,z2ss,z1ss0,z2ss0,shift_ang,angle_ss;
int kn,inv;
sparse_i4 mp;
struct coord3d coord1;
struct coord3d coord2;
int a;


// params definition
strcpy(nom_fitxer,"sresmx");
strcpy(nom_fitxer2,"sresmx");
sprintf(nom_fitxer2,"%s.pesos",nom_fitxer);
Ncol=Nfil;
mp.n_fil=Nfil; // image pixels in X-Y-directions 
mp.n_col=Ncol; // image pixels in X-Y-directions 
mp.Nbins=(2*ndet)-1; // transaxial detectors 
mp.tbin=tdet/2.;
mp.tpix=mp.Nbins*mp.tbin/Nfil; // image pixel size in X-Y-directions in mm
mp.Ntallsima=2*ndet_z-1; // image pixel in Z-direction
if(axial_diff==0) mp.Ntallsima=ndet_z;  
mp.lvox=tdet_z/2.; // image pixel size in Z-direction in mm 
if(axial_diff==0) mp.lvox=tdet_z;
mp.Nang=Nang;
mp.ndet_z=ndet_z;  
mp.tdet_z=tdet_z;
mp.ddet=ddet;
mp.Nsubsets=Nsubsets;
mp.axial_diff=axial_diff;
mp.inv=0;
Nbt=mp.Nbins*Nang*num_imagenes(ndet_z,axial_diff);	// sinogram elements 
Nb_os=Nbt/Nsubsets; // sinogram elements per subset 
Nvox=Nfil*Ncol*mp.Ntallsima; // image elements 
d2=0.5*Nfil*mp.tpix; //medio FOV
if(fac_comp!=666) NITEMS_INI=(Nb_os/fac_comp)*Nvox; // matrix elements/fac_comp  
if(fac_comp==666) NITEMS_INI=Nb_os*Nvox; // no sparse matrix elements  
NITEMS=NITEMS_INI; // sparse matrix elements estimation

// no sparse warning
if(fac_comp==666) printf("\n**\nWARNING! No sparse format was selected. You could have memory problems...\n**\n"); 


// Origin coords 
z0=tdet_z/2.; // mm 
angle0=0.; // degrees 
n0=((-(float)mp.Nbins/2.)+0.5)*mp.tbin; // mm



// steps
inc_zp=tdet_z; // axial increment in mm
inc_ang=180./Nang; // rotating in 180 degrees in degrees
inc_n=mp.tbin; // transaxial increment in mm 
shift_ang=atan(0.5*inc_n/ddet); // shift angle for sub-sample

// Display for selected params
printf("\n**************************************************************");
printf("\nVolumen de reconstrucción: %dx%dx%d (pixel:%f mm y corte:%f mm)\n",Nfil,Ncol,mp.Ntallsima,mp.tpix,mp.lvox);
printf("\nCaracterísticas del escáner: Detectores planos: %d bins (tamaño:%f mm) x %d ángulos\n",mp.Nbins,mp.tbin,Nang);
printf("Maximum axial difference is %d (Sinogramas totales: %d)\n",axial_diff,num_imagenes(ndet_z,axial_diff));
printf("%d detectores axiales de %f mm. Detectores separados: %f mm\n\n\n",ndet_z,tdet_z,ddet);
printf("\n**************************************************************");


// memory allocation
printf("\n\nInitial memory estimation to allocate per subset %1.0f Mb\n\n",NITEMS*sizeof(float)/1048576.);
if((mp.ar = (float *)calloc(NITEMS,sizeof(float)))==NULL) error_fesmp(53,NITEMS*sizeof(float)/1048576);
if((mp.ja = (int *)calloc(NITEMS,sizeof(int)))==NULL) error_fesmp(53,NITEMS*sizeof(int)/1048576);
mp.ia = (int *)calloc(Nbt+1,sizeof(int));


// SUBSET control  
for(nsb=0;nsb<Nsubsets;nsb++)
{
if(opc_subset==nsb || opc_subset==-1)
	{
		//inicializa mp y genera nombre matriz correspondiente al subset
		mp.ne=0;
		sprintf(nom_fitxer_sb,"%s.%d",nom_fitxer,nsb);
		printf("\n************\n\nSubset: %d\n",nsb);
			
		// AXIAL LOR control
		for(k1=0,z1=z0,ind_proj=0;k1<ndet_z;z1+=inc_zp,k1++)
		{
		for(k2=0,z2=z0;k2<ndet_z;z2+=inc_zp,k2++)
		{
			
			if(abs(k1-k2)>axial_diff)
			{
				continue;
			}
		
			// n-angular independent LOR params definition
			printf("\nSinogram: z1=%d  z2=%d (memory used: %f Mb)",k1,k2,6.*mp.ne/1048576.);
			
			// ANGULAR LOR control
			for(l=nsb,angle=angle0+(nsb*inc_ang);l<Nang;l+=Nsubsets,angle+=inc_ang*Nsubsets)
			{
				
				// TRANSAXIAL LOR control
				for(in=0,n=n0;in<mp.Nbins;in++,ind_proj++,n+=inc_n)
				{
					// calcula contribución de LOR en otras posiciones de los detectores
					n_i=n;
					sum=0.;
					while(n_i<=0 && n_i>=n0) 
					{
						ang_i=atan((n_i-n)/0.5*ddet);
						n_i-=inc_n;
						sum+=cos(ang_i*M_PI/180.);
					}
					n_i=n;
					while(n_i>0 && n_i<=-n0) 
					{
						ang_i=atan((n_i-n)/0.5*ddet);
						n_i+=inc_n;
						sum+=cos(ang_i*M_PI/180.);
					}
				
				
					
					//bucles para z1,z2 y n que no aumentan el valor de ind_proj
					// es decir formarán parte de la misma LOR
					nss0=n-(inc_n/4.);
					z1ss0=z1-inc_zp/4.;
					z2ss0=z2-inc_zp/4.	;
					//printf("**************\n: z1:%f,z2:%f,n:%f\n**********************\n",z1,z2,n);
					
					z1ss=z1;
					z2ss=z2;
					nss=n;
					a=0;
					mp.ia[ind_proj]=mp.ne;
					//printf("NEW Subsampl\n");					
					//printf("*********\nSubsample z1:%f (entre %f y %f)\n",z1ss,z1ss0,z1+inc_zp/4);
  					for(z1ss=z1ss0,z2ss=z2ss0;z1ss<=z1+0.001+inc_zp/4.;z1ss+=inc_zp/2.,z2ss+=inc_zp/2.)
  					{

 					for(kn=0,nss=nss0;kn<2;nss+=inc_n/2.,kn++)
					{
					angle_ss=angle;
					seno=sin(angle_ss*M_PI/180.);
					coseno=cos(angle_ss*M_PI/180.);
					//printf("\nangle_ss:%f",angle_ss);
						z_med=fabs(z1ss+z2ss)/2;  
						zmax=maximo(z1ss,z2ss)+epsilon;
						zmin=minimo(z1ss,z2ss)-epsilon;
						phi=atan((z2ss-z1ss)/ddet); 
						tan_phi=tan(phi);
						factor=cos(phi);
						mod_vlor=sqrt((fabs(z2ss-z1ss)*fabs(z2ss-z1ss))+(ddet*ddet));
						
						//calculate non-null elements only 
						calcul_ptos_corte(d2, nss, seno, coseno, tan_phi,z_med, &coord1, &coord2);
						dmax=10000*maximo(fabs(coord1.x-coord2.x),fabs(coord1.y-coord2.y));
						ray_tracing_lor_sub(&coord1,&coord2, &mp, Nfil, mp.Ntallsima, mp.tpix, mp.lvox, dmax, ind_proj, Nb_os, fac_comp, sum);
		
						// redimensiona memoria si es necesario
						if(mp.ne >= NITEMS-Nvox)
						{
							NITEMS+=NITEMS_INI;
							if((mp.ar = (float *)realloc(mp.ar,NITEMS*sizeof(float)))==NULL) 
							{ 
								error_fesmp(53,NITEMS*sizeof(float)/1048576.);
							}
							if((mp.ja = (int *)realloc(mp.ja,NITEMS*sizeof(int)))==NULL)
							{
								error_fesmp(53,NITEMS*sizeof(int)/1048576.);
							}
							printf("\nrealloc for %1.0f Mb\n",NITEMS*sizeof(float)/1048576.);
						}
					}
					}
				} // transaxial LOR control
			} // angular LOR control
		} // axial LOR control
		} // axial LOR control
	
	
		// Writting matrix for each subset
		mp.ia[Nb_os]=mp.ne; printf("\n\n\nNo-null elements: %d\n",mp.ne);
		inv=0;
		// for sparse format
		if(fac_comp!=666)
		{
			guarda_pesos_xOS(&mp,nom_fitxer_sb);
		}
		
		// for no sparse format
		if(fac_comp==666)
		{
			printf("\n\nOutput: No-Sparse format\n");
			sprintf(nom_fitxer_sb2,"%s.no_sparse",nom_fitxer_sb);
			if(( mat = fopen(nom_fitxer_sb2,"wb"))==NULL) error_fesmp(30,0);
			fwrite (mp.ar,sizeof(float),(size_t)mp.ne,mat);
			fclose (mat);	
		}
	}
}
return 0;
}









// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////




// Función que obtiene la matriz del sistema rPET, Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 
// modificado para sub-sample bin sin cruzados 



// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int build_srm_rpet_ray_tracing_sub_grid(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset)
{
FILE *mat;
char nom_fitxer[128],nom_fitxer_sb[128],nom_fitxer_sb2[128],nom_fitxer2[128];
int in,ind_proj,l,k2,k1,Ncol,nsb,Nb_os;
double NITEMS_INI,NITEMS;
float z1,z2,z_med,n,angle;
float n0,n_i,angle0,z0,inc_zp,inc_ang,ang_i,inc_n,d2,sum;
float coseno,seno,phi,factor,mod_vlor,zmin,zmax,tan_phi;
float epsilon=1e-4,dmax;
float nss,nss0,z1ss,z2ss,z1ss0,z2ss0,shift_ang,angle_ss;
int kn,inv;
sparse_i4 mp;
struct coord3d coord1;
struct coord3d coord2;
int a;


// params definition
strcpy(nom_fitxer,"sresmx");
strcpy(nom_fitxer2,"sresmx");
sprintf(nom_fitxer2,"%s.pesos",nom_fitxer);
Ncol=Nfil;
mp.n_fil=Nfil; // image pixels in X-Y-directions 
mp.n_col=Ncol; // image pixels in X-Y-directions 
mp.Nbins=(2*ndet)-1; // transaxial detectors 
mp.tbin=tdet/2.;
mp.tpix=mp.Nbins*mp.tbin/Nfil; // image pixel size in X-Y-directions in mm
mp.Ntallsima=2*ndet_z-1; // image pixel in Z-direction
if(axial_diff==0) mp.Ntallsima=ndet_z;  
mp.lvox=tdet_z/2.; // image pixel size in Z-direction in mm 
if(axial_diff==0) mp.lvox=tdet_z;
mp.Nang=Nang;
mp.ndet_z=ndet_z;  
mp.tdet_z=tdet_z;
mp.ddet=ddet;
mp.Nsubsets=Nsubsets;
mp.axial_diff=axial_diff;
mp.inv=0;
Nbt=mp.Nbins*Nang*num_imagenes(ndet_z,axial_diff);	// sinogram elements 
Nb_os=Nbt/Nsubsets; // sinogram elements per subset 
Nvox=Nfil*Ncol*mp.Ntallsima; // image elements 
d2=0.5*Nfil*mp.tpix; //medio FOV
if(fac_comp!=666) NITEMS_INI=(Nb_os/fac_comp)*Nvox; // matrix elements/fac_comp  
if(fac_comp==666) NITEMS_INI=Nb_os*Nvox; // no sparse matrix elements  
NITEMS=NITEMS_INI; // sparse matrix elements estimation

// no sparse warning
if(fac_comp==666) printf("\n**\nWARNING! No sparse format was selected. You could have memory problems...\n**\n"); 


// Origin coords 
z0=tdet_z/2.; // mm 
angle0=0.; // degrees 
n0=((-(float)mp.Nbins/2.)+0.5)*mp.tbin; // mm



// steps
inc_zp=tdet_z; // axial increment in mm
inc_ang=180./Nang; // rotating in 180 degrees in degrees
inc_n=mp.tbin; // transaxial increment in mm 
shift_ang=atan(0.5*inc_n/ddet); // shift angle for sub-sample

// Display for selected params
printf("\n**************************************************************");
printf("\nVolumen de reconstrucción: %dx%dx%d (pixel:%f mm y corte:%f mm)\n",Nfil,Ncol,mp.Ntallsima,mp.tpix,mp.lvox);
printf("\nCaracterísticas del escáner: Detectores planos: %d bins (tamaño:%f mm) x %d ángulos\n",mp.Nbins,mp.tbin,Nang);
printf("Maximum axial difference is %d (Sinogramas totales: %d)\n",axial_diff,num_imagenes(ndet_z,axial_diff));
printf("%d detectores axiales de %f mm. Detectores separados: %f mm\n\n\n",ndet_z,tdet_z,ddet);
printf("\n**************************************************************");


// memory allocation
printf("\n\nInitial memory estimation to allocate per subset %1.0f Mb\n\n",NITEMS*sizeof(float)/1048576.);
if((mp.ar = (float *)calloc(NITEMS,sizeof(float)))==NULL) error_fesmp(53,NITEMS*sizeof(float)/1048576);
if((mp.ja = (int *)calloc(NITEMS,sizeof(int)))==NULL) error_fesmp(53,NITEMS*sizeof(int)/1048576);
mp.ia = (int *)calloc(Nbt+1,sizeof(int));


// SUBSET control  
for(nsb=0;nsb<Nsubsets;nsb++)
{
if(opc_subset==nsb || opc_subset==-1)
	{
		//inicializa mp y genera nombre matriz correspondiente al subset
		mp.ne=0;
		sprintf(nom_fitxer_sb,"%s.%d",nom_fitxer,nsb);
		printf("\n************\n\nSubset: %d\n",nsb);
			
		// AXIAL LOR control
		for(k1=0,z1=z0,ind_proj=0;k1<ndet_z;z1+=inc_zp,k1++)
		{
		for(k2=0,z2=z0;k2<ndet_z;z2+=inc_zp,k2++)
		{
			
			if(abs(k1-k2)>axial_diff)
			{
				continue;
			}
		
			// n-angular independent LOR params definition
			printf("\nSinogram: z1=%d  z2=%d (memory used: %f Mb)",k1,k2,6.*mp.ne/1048576.);
			
			// ANGULAR LOR control
			for(l=nsb,angle=angle0+(nsb*inc_ang);l<Nang;l+=Nsubsets,angle+=inc_ang*Nsubsets)
			{
				
				// TRANSAXIAL LOR control
				for(in=0,n=n0;in<mp.Nbins;in++,ind_proj++,n+=inc_n)
				{
					// calcula contribución de LOR en otras posiciones de los detectores
					n_i=n;
					sum=0.;
					while(n_i<=0 && n_i>=n0) 
					{
						ang_i=atan((n_i-n)/0.5*ddet);
						n_i-=inc_n;
						sum+=cos(ang_i*M_PI/180.);
					}
					n_i=n;
					while(n_i>0 && n_i<=-n0) 
					{
						ang_i=atan((n_i-n)/0.5*ddet);
						n_i+=inc_n;
						sum+=cos(ang_i*M_PI/180.);
					}
				
				
					
					//bucles para z1,z2 y n que no aumentan el valor de ind_proj
					// es decir formarán parte de la misma LOR
					nss0=n-(inc_n/4.);
					z1ss0=z1-inc_zp/4.;
					z2ss0=z2-inc_zp/4.	;
					//printf("**************\n: z1:%f,z2:%f,n:%f\n**********************\n",z1,z2,n);
					
					z1ss=z1;
					z2ss=z2;
					nss=n;
					a=0;
					mp.ia[ind_proj]=mp.ne;
					//printf("NEW Subsampl\n");					
					//printf("*********\nSubsample z1:%f (entre %f y %f)\n",z1ss,z1ss0,z1+inc_zp/4);
  					for(z1ss=z1ss0,z2ss=z2ss0;z1ss<=z1+0.001+inc_zp/4.;z1ss+=inc_zp/4.,z2ss+=inc_zp/4.)
  					{

 					for(kn=0,nss=nss0;kn<3;nss+=inc_n/4.,kn++)
					{
					angle_ss=angle;
					seno=sin(angle_ss*M_PI/180.);
					coseno=cos(angle_ss*M_PI/180.);
					//printf("\nangle_ss:%f",angle_ss);
						z_med=fabs(z1ss+z2ss)/2;  
						zmax=maximo(z1ss,z2ss)+epsilon;
						zmin=minimo(z1ss,z2ss)-epsilon;
						phi=atan((z2ss-z1ss)/ddet); 
						tan_phi=tan(phi);
						factor=cos(phi);
						mod_vlor=sqrt((fabs(z2ss-z1ss)*fabs(z2ss-z1ss))+(ddet*ddet));
						
						//calculate non-null elements only 
						calcul_ptos_corte(d2, nss, seno, coseno, tan_phi,z_med, &coord1, &coord2);
						dmax=10000*maximo(fabs(coord1.x-coord2.x),fabs(coord1.y-coord2.y));
						ray_tracing_lor_sub(&coord1,&coord2, &mp, Nfil, mp.Ntallsima, mp.tpix, mp.lvox, dmax, ind_proj, Nb_os, fac_comp, sum);
		
						// redimensiona memoria si es necesario
						if(mp.ne >= NITEMS-Nvox)
						{
							NITEMS+=NITEMS_INI;
							if((mp.ar = (float *)realloc(mp.ar,NITEMS*sizeof(float)))==NULL) 
							{ 
								error_fesmp(53,NITEMS*sizeof(float)/1048576.);
							}
							if((mp.ja = (int *)realloc(mp.ja,NITEMS*sizeof(int)))==NULL)
							{
								error_fesmp(53,NITEMS*sizeof(int)/1048576.);
							}
							printf("\nrealloc for %1.0f Mb\n",NITEMS*sizeof(float)/1048576.);
						}
					}
					}
				} // transaxial LOR control
			} // angular LOR control
		} // axial LOR control
		} // axial LOR control
	
	
		// Writting matrix for each subset
		mp.ia[Nb_os]=mp.ne; printf("\n\n\nNo-null elements: %d\n",mp.ne);
		inv=0;
		// for sparse format
		if(fac_comp!=666)
		{
			guarda_pesos_xOS(&mp,nom_fitxer_sb);
		}
		
		// for no sparse format
		if(fac_comp==666)
		{
			printf("\n\nOutput: No-Sparse format\n");
			sprintf(nom_fitxer_sb2,"%s.no_sparse",nom_fitxer_sb);
			if(( mat = fopen(nom_fitxer_sb2,"wb"))==NULL) error_fesmp(30,0);
			fwrite (mp.ar,sizeof(float),(size_t)mp.ne,mat);
			fclose (mat);	
		}
	}
}
return 0;
}







// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////




// Función que obtiene la matriz del sistema rPET por el Método de WU, Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 




// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int build_srm_rpet_distance(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset,float pesminim)
{
FILE *mat;
char nom_fitxer[128],nom_fitxer_sb[128],nom_fitxer_sb2[128],nom_fitxer2[128];
int i,j,k,in,ind_proj,l,k2,k1,Ncol,Ntallsima,Nbins,ind_ima,nsb,Nb_os;
double NITEMS_INI,NITEMS;
float z1,z2,z_med,n,n_i,angle,ang_i;
float n0,angle0,xpix0,ypix0,zpix0,z0,inc_zp,inc_ang,inc_n,inc_zi,inc_x,inc_y,xpix,ypix,zpix;
float coseno,seno,phi,factor,mod_vlor,pes,tbin,tpix,lvox,dplor,zmin,zmax,sum;
float epsilon=1e-4;
time_t time1,time2;
int ctr,ctr2,npl;
float dist_cen;
sparse_i4 mp;
int aaa;
// params definition
strcpy(nom_fitxer,"sresmx");
strcpy(nom_fitxer2,"sresmx");
sprintf(nom_fitxer2,"%s.pesos",nom_fitxer);
Ncol=Nfil; // image pixels in X-Y-directions 
Nbins=(2*ndet)-1; // transaxial detectors 
tbin=tdet/2.;
tpix=Nbins*tbin/Nfil; // image pixel size in X-Y-directions in mm
Ntallsima=2*ndet_z-1; // image pixel in Z-direction
if(axial_diff==0) Ntallsima=ndet_z;  
lvox=tdet_z/2.; // image pixel size in Z-direction in mm 
if(axial_diff==0) lvox=tdet_z;  
Nbt=Nbins*Nang*num_imagenes(ndet_z,axial_diff);	// sinogram elements 
Nb_os=Nbt/Nsubsets; // sinogram elements per subset 
Nvox=Nfil*Ncol*Ntallsima; // image elements 
if(fac_comp!=666) NITEMS_INI=(Nb_os/fac_comp)*Nvox; // matrix elements/fac_comp  
if(fac_comp==666) NITEMS_INI=Nb_os*Nvox; // no sparse matrix elements  
NITEMS=NITEMS_INI; // sparse matrix elements estimation  


// no sparse warning
if(fac_comp==666) printf("\n**\nWARNING! No sparse format was selected. You could have memory problems...\n**\n"); 


// Origin coords 
z0=tdet_z/2.; // mm 
angle0=0.; // degrees 
n0=((-(float)Nbins/2.)+0.5)*tbin; // mm
zpix0=z0; // mm 
xpix0=((-(float)Nfil/2.)+0.5)*tpix; // mm
ypix0=((-(float)Nfil/2.)+0.5)*tpix; // mm


// steps
inc_zp=tdet_z; // axial increment in mm
inc_ang=180./Nang; // rotating in 180 degrees in degrees
inc_n=tbin; // transaxial increment in mm 
inc_zi=inc_zp/2.; // image axial increment in mm
if(axial_diff==0) inc_zi=inc_zp;
inc_x=inc_y=tpix; // XY image increment in mm 


// Display for selected params
printf("\n**************************************************************");
printf("\nVolumen de reconstrucción: %dx%dx%d (pixel:%f mm y corte:%f mm)\n",Nfil,Ncol,Ntallsima,tpix,lvox);
printf("\nCaracterísticas del escáner: Detectores planos: %d bins (tamaño:%f mm) x %d ángulos\n",Nbins,tbin,Nang);
printf("Maximum axial difference is %d (Sinogramas totales: %d)\n",axial_diff,num_imagenes(ndet_z,axial_diff));
printf("%d detectores axiales de %f mm. Detectores separados: %f mm\n\n\n",ndet_z,tdet_z,ddet);
printf("\n**************************************************************");


// memory allocation
printf("\n\nInitial memory estimation to allocate per subset %1.0f Mb\n\n",NITEMS*sizeof(float)/1048576.);
if((mp.ar = (float *)calloc(NITEMS,sizeof(float)))==NULL) error_fesmp(53,NITEMS*sizeof(float)/1048576);
if((mp.ja = (int *)calloc(NITEMS,sizeof(int)))==NULL) error_fesmp(53,NITEMS*sizeof(int)/1048576);
mp.ia = (int *)calloc(Nbt+1,sizeof(int));


// SUBSET control  
for(nsb=0;nsb<Nsubsets;nsb++)
{
	//inicializa mp y genera nombre matriz correspondiente al subset
	mp.ne=0;
	sprintf(nom_fitxer_sb,"%s.%d",nom_fitxer,nsb);
	printf("\n************\n\nSubset: %d\n",nsb);
	ctr=1;
	ctr2=0;
	
	// AXIAL LOR control
	for(k1=0,z1=z0,ind_proj=0;k1<ndet_z;z1+=inc_zp,k1++)
	{
	for(k2=0,z2=z0;k2<ndet_z;z2+=inc_zp,k2++)
	{
		
		if(abs(k1-k2)>axial_diff)
		{
			continue;
		}
			
		//Control de tiempo para finalizar
		if(ctr==1) 
		{
			time1=time(NULL);
			ctr2++;
			if(ctr2==ndet_z) ctr=2;
		}
		else 
		{
			time2=time(NULL);
			npl=((Nsubsets-nsb+1)*(ndet_z))-k1;
						
			if(difftime(time2,time1)<2.) 
			{
			ctr=0;
			}
			else
			{
			printf("\tTime to finish %1.1f minutes ",npl*difftime(time2,time1)/60.);
			ctr=1;
			ctr2=0;
		}
		}
	
		// n-angular independent LOR params definition
		z_med=fabs(z1+z2)/2;  
		zmax=maximo(z1,z2)+epsilon;
		zmin=minimo(z1,z2)-epsilon;
		printf("\nSinogram: z1=%d  z2=%d (memory used: %f Mb)",k1,k2,6.*mp.ne/1048576.);
		phi=atan((z2-z1)/ddet); 
		factor=cos(phi);
		mod_vlor=sqrt((fabs(z2-z1)*fabs(z2-z1))+(ddet*ddet));
		
		// ANGULAR LOR control
		for(l=nsb,angle=angle0+(nsb*inc_ang);l<Nang;l+=Nsubsets,angle+=inc_ang*Nsubsets)
		{
			// angular LOR params definition
			seno=sin(angle*M_PI/180.);
			coseno=cos(angle*M_PI/180.);
			
			// TRANSAXIAL LOR control
			for(in=0,n=n0;in<Nbins;in++,ind_proj++,n+=inc_n)
			{
				//non-null elements
				mp.ia[ind_proj]=mp.ne; 

				// Z-Y-X IMAGE control
				for(k=0,zpix=zpix0,ind_ima=0;k<Ntallsima;k++,zpix+=inc_zi)
				{
				
					// No se recorren los cortes de la imagen que quedan fuera de los 
					// límites axiales de la L.O.R. Aceleración zfast, paguiar nov2006
					if(fac_comp!=666)
					{
						if(zpix<zmin || zpix>zmax)
						{
							ind_ima+=Nfil*Ncol;
							continue;
						}
						
					}
					
					aaa=0;
// 					if(k1==5 && k2==5) aaa=1;
// 					if(aaa!=1) continue;
// 					
					for(i=0,xpix=xpix0;i<Nfil;i++,xpix+=inc_x)
					{
					for(j=0,ypix=ypix0;j<Ncol;j++,ypix+=inc_y,ind_ima++)
					{
						dist_cen=sqrt((xpix*xpix)+(ypix*ypix));
						// Perpendicular distance from image point to LOR and "pes"
						dplor=d_punto_lor(xpix,ypix,zpix,n,seno,coseno,phi,z_med,mod_vlor);
						if(dplor>tbin || dist_cen>(0.5*Nfil*tpix)) pes=0.;
						else pes=calcul_pes(tbin,dplor,factor);
						
						// PESMINIM discrimitation
						if(pes>pesminim)
						{
							mp.ja[mp.ne]=ind_ima; 
							mp.ar[mp.ne]=0.;
							n_i=n;
							sum=0.;
							//printf("\nxpix:%f (%d) ypix:%f (%d) zpix:%f (%d) n:%f",xpix,i,ypix,j,zpix,k,n);      
							
							// Adding equivalent angular positions
  							while(n_i<=0 && n_i>=n0) 
  							{
  								ang_i=atan((n_i-n)/0.5*ddet);
  								mp.ar[mp.ne]+=pes*cos(ang_i*M_PI/180.);
  								n_i-=inc_n;
  								sum+=cos(ang_i*M_PI/180.);
  							}
  							n_i=n;
  							while(n_i>0 && n_i<=-n0) 
  							{
  								ang_i=atan((n_i-n)/0.5*ddet);
  								mp.ar[mp.ne]+=pes*cos(ang_i*M_PI/180.);
  								n_i+=inc_n;
  								sum+=cos(ang_i*M_PI/180.);
  							}
/*							mp.ar[mp.ne]=pes;*/
							mp.ne++;
							
							// Realloc memory if it is needed 
							if(mp.ne >= NITEMS)
							{
								NITEMS+=NITEMS_INI;
								if((mp.ar = (float *)realloc(mp.ar,NITEMS*sizeof(float)))==NULL) 
								{ 
									error_fesmp(53,NITEMS*sizeof(float)/1048576.);
								}
								if((mp.ja = (int *)realloc(mp.ja,NITEMS*sizeof(int)))==NULL)
								{
									error_fesmp(53,NITEMS*sizeof(int)/1048576.);
								}
								printf("\nrealloc for %1.0f Mb\n",NITEMS*sizeof(float)/1048576.);
							}
						}
						
						// writting null values if no sparse format was selected (666)
						else if(pes<=pesminim && fac_comp==666)
						{
							mp.ja[mp.ne]=ind_ima;
							mp.ar[mp.ne]=0.;
							mp.ne++;
							
							// Realloc memory if it is needed 
							if(mp.ne >= NITEMS)
							{
								NITEMS+=NITEMS_INI;
								if((mp.ar = (float *)realloc(mp.ar,NITEMS*sizeof(float)))==NULL)
								{
									error_fesmp(53,NITEMS*sizeof(float));
								}
								if((mp.ja = (int *)realloc(mp.ja,NITEMS*sizeof(int)))==NULL) 
								{
									error_fesmp(53,NITEMS*sizeof(int));
								}
								printf("\nrealloc for %1.0f Mb\n",NITEMS*sizeof(float)/1048576.);
							} 
						}
					} // x-image control
					} // y-image control
				} // z-image control
			} // transaxial LOR control
		} // angular LOR control
	} // axial LOR control
	} // axial LOR control


	// Writting matrix for each subset
	mp.ia[Nb_os]=mp.ne; printf("\n\n\nNo-null elements: %d\n",mp.ne);
	
	
	//relleno campos de la matriz
	mp.n_fil=Nfil;
	mp.n_col=Ncol;
	mp.Ntallsima=Ntallsima;
	mp.tpix=tpix;
	mp.lvox=lvox;
	mp.ndet_z=ndet_z;
	mp.tdet_z=tdet_z;
	mp.Nbins=Nbins;
	mp.Nang=Nang;
	mp.tbin=tbin;
	mp.ddet=ddet;
	mp.Nsubsets=Nsubsets;
	mp.axial_diff=axial_diff;
	mp.inv=0;
	

	// for sparse format
	if(fac_comp!=666)
	{
		guarda_pesos_xOS(&mp,nom_fitxer_sb);
	}
	
	// for no sparse format
	if(fac_comp==666)
	{
		printf("\n\nOutput: No-Sparse format\n");
		sprintf(nom_fitxer_sb2,"%s.no_sparse",nom_fitxer_sb);
		if(( mat = fopen(nom_fitxer_sb2,"wb"))==NULL) error_fesmp(30,0);
		fwrite (mp.ar,sizeof(float),(size_t)mp.ne,mat);
		fclose (mat);	
	}
}
return 0;
}








// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////





// Función que obtiene la matriz del sistema rPET, Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 





// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int build_srm_rpet_ray_tracing_pWu(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int npix_exp2d,int fac_comp,int opc_subset)
{
FILE *mat;
char nom_fitxer[128],nom_fitxer_sb[128],nom_fitxer_sb2[128],nom_fitxer2[128];
int in,ind_proj,l,k2,k1,Ncol,nsb,Nb_os;
double NITEMS_INI,NITEMS;
float z1,z2,angle;
float n0,n_i,angle0,z0,inc_zp,inc_ang,ang_i,inc_n,d2,sum;
float factor,zmin,zmax,tan_phi;
float epsilon=1e-4,dmax;
float xcon,ycon,zcon,rcon1,rcon2,xini,yini,zini;
sparse_i4 mp;
struct coord3d coord1;
struct coord3d coord2;
struct lor_params lor;
struct imagen ima;
int desf,con,con2;
desf=0;con2=0;
// params definition
strcpy(nom_fitxer,"sresmx");
strcpy(nom_fitxer2,"sresmx");
sprintf(nom_fitxer2,"%s.pesos",nom_fitxer);
Ncol=Nfil;
mp.n_fil=Nfil; // image pixels in X-Y-directions 
mp.n_col=Ncol; // image pixels in X-Y-directions 
mp.Nbins=(2*ndet)-1; // transaxial detectors 
mp.tbin=tdet/2.;
mp.tpix=mp.Nbins*mp.tbin/Nfil; // image pixel size in X-Y-directions in mm
mp.Ntallsima=2*ndet_z-1; // image pixel in Z-direction
if(axial_diff==0) mp.Ntallsima=ndet_z;  
mp.lvox=0.5*tdet_z; // image pixel size in Z-direction in mm 
if(axial_diff==0) mp.lvox=tdet_z;
mp.Nang=Nang;
mp.ndet_z=ndet_z;  
mp.tdet_z=tdet_z;
mp.ddet=ddet;
mp.Nsubsets=Nsubsets;
mp.axial_diff=axial_diff;
mp.inv=0;
mp.sym=0;		
Nbt=mp.Nbins*Nang*num_imagenes(ndet_z,axial_diff);	// sinogram elements 
Nb_os=Nbt/Nsubsets; // sinogram elements per subset 
Nvox=Nfil*Ncol*mp.Ntallsima; // image elements 
d2=0.5*Nfil*mp.tpix; //medio FOV
if(fac_comp!=666) NITEMS_INI=(Nb_os/fac_comp)*Nvox; // matrix elements/fac_comp  
if(fac_comp==666) NITEMS_INI=Nb_os*Nvox; // no sparse matrix elements  
NITEMS=NITEMS_INI; // sparse matrix elements estimation  


// no sparse warning
if(fac_comp==666) printf("\n**\nWARNING! No sparse format was selected. You could have memory problems...\n**\n"); 


// steps
inc_zp=tdet_z; // axial increment in mm
inc_ang=180./Nang; // rotating in 180 degrees in degrees
inc_n=mp.tbin; // transaxial increment in mm 

// Origin coords 
z0=tdet_z/2.; // mm 
angle0=0.; // degrees 
n0=-0.5*((float)mp.Nbins-1.)*mp.tbin; // mm


// Display for selected params
printf("\n**************************************************************");
printf("\nVolumen de reconstrucción: %dx%dx%d (pixel:%f mm y corte:%f mm)\n",Nfil,Ncol,mp.Ntallsima,mp.tpix,mp.lvox);
printf("\nCaracterísticas del escáner: Detectores planos: %d bins (tamaño:%f mm) x %d ángulos\n",mp.Nbins,mp.tbin,Nang);
printf("Maximum axial difference is %d (Sinogramas totales: %d)\n",axial_diff,num_imagenes(ndet_z,axial_diff));
printf("%d detectores axiales de %f mm. Detectores separados: %f mm\n\n\n",ndet_z,tdet_z,ddet);
printf("\n**************************************************************");


// memory allocation
printf("\n\nInitial memory estimation to allocate per subset %1.0f Mb\n\n",NITEMS*sizeof(float)/1048576.);
if((mp.ar = (float *)calloc(NITEMS,sizeof(float)))==NULL) error_fesmp(53,NITEMS*sizeof(float)/1048576);
if((mp.ja = (int *)calloc(NITEMS,sizeof(int)))==NULL) error_fesmp(53,NITEMS*sizeof(int)/1048576);
mp.ia = (int *)calloc(Nbt+1,sizeof(int));

//lee_template_hdr("ptos_corte.hdr",&ima);
// SUBSET control  
for(nsb=0;nsb<Nsubsets;nsb++)
{
if(opc_subset==nsb || opc_subset==-1)
	{
		//inicializa mp y genera nombre matriz correspondiente al subset
		mp.ne=0;
		sprintf(nom_fitxer_sb,"%s.%d",nom_fitxer,nsb);
		printf("\n*************************\n\nSubset: %d *******************************\n",nsb);
			
		// AXIAL LOR control
		for(k1=0,z1=z0,ind_proj=0;k1<ndet_z;z1+=inc_zp,k1++)
		{
			//printf("\nInd_proj:%d\tmp.ia:%d\n",ind_proj,mp.ia[ind_proj-1]);
		for(k2=0,z2=z0;k2<ndet_z;z2+=inc_zp,k2++)
		{
			
			if(abs(k1-k2)>axial_diff)
			{
				continue;
			}
		
 
			// n-angular independent LOR params definition
			lor.z_med=fabs(z1+z2)/2;  
			zmax=maximo(z1,z2)+epsilon;
			zmin=minimo(z1,z2)-epsilon;
			printf("\nSinogram: k1=%d (%f) k2=%d (%f) (memory used: %f Mb) ",k1,z1,k2,z2,6.*mp.ne/1048576.);
			lor.phi=atan((z2-z1)/ddet); 
			tan_phi=tan(lor.phi);
			factor=cos(lor.phi);
			lor.mod_vlor=sqrt((fabs(z2-z1)*fabs(z2-z1))+(ddet*ddet));
			
			// ANGULAR LOR control
			for(l=nsb,angle=angle0+(nsb*inc_ang);l<Nang;l+=Nsubsets,angle+=inc_ang*Nsubsets)
			{
				// angular LOR params definition
				lor.seno=sin(angle*M_PI/180.);
				lor.coseno=cos(angle*M_PI/180.);
				con=0;
				
				// TRANSAXIAL LOR control
				for(in=0,lor.n=n0;in<mp.Nbins;in++,ind_proj++,lor.n+=inc_n)
				{
					
					//printf("mp.ia[%d]:%d\n",ind_proj-1,mp.ia[ind_proj-1]);
					if(ind_proj>=2)
					{
						
						if(mp.ia[ind_proj-1]-mp.ia[ind_proj-2]<0 || mp.ia[ind_proj-1]-mp.ia[ind_proj-2]>mp.n_col*mp.n_fil*mp.Ntallsima) 
						{
							printf("\n\n\n¡Error inesperado (hay una LOR con %d elementos no nulos)!\n\n\n\n\n",mp.ia[ind_proj-1]-mp.ia[ind_proj-2]);
							exit(0);
						}
					}
				
					// calcula contribución de LOR en otras posiciones de los detectores
					n_i=lor.n;
					sum=0.;
 					while(n_i<=0 && n_i>=n0) 
 					{
 						ang_i=atan((n_i-lor.n)/0.5*ddet);
 						n_i-=inc_n;
 						sum+=cos(ang_i*M_PI/180.);
 					}
 					n_i=lor.n;
 					while(n_i>0 && n_i<=-n0) 
 					{
 						ang_i=atan((n_i-lor.n)/0.5*ddet);
 						n_i+=inc_n;
 						sum+=cos(ang_i*M_PI/180.);
 					}
					//sum=1.;
/*					if(((in==29 && k1==0 && k2==34) || (in==29 && k1==34 && k2==0) || (in==29 && k1==17 && k2==17)) || ((in==29 && k1==10 && k2==24) || (in==29 && k1==24 && k2==10)))*/ 
					
					
// 					if((in==29 && k1==0 && k2==34) || (in==29 && k1==34 && k2==0) || (in==29 && k1==10 && k2==24) || (in==29 && k1==24 && k2==10))
// 					{
					//calculate non-null elements only
					mp.ia[ind_proj]=mp.ne; 
					calcul_ptos_corte(d2, lor.n, -lor.seno, lor.coseno, -tan_phi,lor.z_med, &coord1, &coord2);
					//printf("\nP1: %f,%f,%f --- P2: %f,%f,%f  LOR: n:%f z_med:%f",coord1.x,coord1.y,coord1.z,coord2.x,coord2.y,coord2.z,lor.n,lor.z_med);
					dmax=10000*maximo(fabs(coord1.x-coord2.x),fabs(coord1.y-coord2.y));
					 //if(k1==5 && k2==5){
						 
						if(con==0) 
						{ 
							xcon=coord1.x; ycon=coord1.y; zcon=coord1.z;
							if(con2!=0) { xcon=xini; ycon=yini; zcon=zini;}
						}
						
						//printf("Cerca? %f, %f, %f OR %f, %f, %f desde %f %f %f \n",coord1.x,coord1.y,coord1.z,coord2.x,coord2.y,coord2.z,xcon,ycon,zcon);		
						rcon1=((xcon-coord1.x)*(xcon-coord1.x))+((ycon-coord1.y)*(ycon-coord1.y))+((zcon-coord1.z)*(zcon-coord1.z));
						rcon2=((xcon-coord2.x)*(xcon-coord2.x))+((ycon-coord2.y)*(ycon-coord2.y))+((zcon-coord2.z)*(zcon-coord2.z));
						if(rcon1>rcon2)
						{
							xcon=coord2.x; ycon=coord2.y; zcon=coord2.z;
							coord2.x=coord1.x; coord2.y=coord1.y; coord2.z=coord1.z;
							coord1.x=xcon; coord1.y=ycon; coord1.z=zcon;
						}
						else
						{							
							xcon=coord1.x; ycon=coord1.y; zcon=coord1.z;
						}
						//printf("%f, %f, %f\n",coord1.x,coord1.y,coord1.z);						
						if(con==0) 
						{
							xini=coord1.x; yini=coord1.y; zini=coord1.z; //printf("*****\ninicial %f, %f, %f\n\n",xini,yini,zini);
							con=1;
					 	}
						con2=1;
						
    					ray_tracing_lor_pWu(&coord1,&coord2, &mp, Nfil, mp.Ntallsima, mp.tpix, mp.lvox, mp.tbin, dmax, ind_proj, Nb_os, fac_comp, sum,&lor, npix_exp2d,&ima,desf);
					// }
			
					// redimensiona memoria si es necesario
					if(mp.ne >= NITEMS-Nvox)
					{
						NITEMS+=NITEMS_INI;
						if((mp.ar = (float *)realloc(mp.ar,NITEMS*sizeof(float)))==NULL) 
						{ 
							error_fesmp(53,NITEMS*sizeof(float)/1048576.);
						}
						if((mp.ja = (int *)realloc(mp.ja,NITEMS*sizeof(int)))==NULL)
						{
							error_fesmp(53,NITEMS*sizeof(int)/1048576.);
						}
						printf("\nrealloc for %1.0f Mb\n",NITEMS*sizeof(float)/1048576.);
					}
						
				} // transaxial LOR control
			} // angular LOR control

		} // axial LOR control
		} // axial LOR control
	
	
		// Writting matrix for each subset
		mp.ia[Nb_os]=mp.ne; printf("\n\n\nNo-null elements: %d\n",mp.ne);

		// for sparse format
		if(fac_comp!=666)
		{ 
			guarda_pesos_xOS(&mp,nom_fitxer_sb);
		}
		
		// for no sparse format
		if(fac_comp==666)
		{
			printf("\n\nOutput: No-Sparse format\n");
			sprintf(nom_fitxer_sb2,"%s.no_sparse",nom_fitxer_sb);
			if(( mat = fopen(nom_fitxer_sb2,"wb"))==NULL) error_fesmp(30,0);
			fwrite (mp.ar,sizeof(float),(size_t)mp.ne,mat);
			fclose (mat);	
		}
	}
desf++;
}

guarda_imagen_hdr("ptos_corte.hdr",&ima);
return 0;
}



















// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////





// Función que obtiene la matriz del sistema rPET, Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 

// Incluye sym en OZ Junio2007



// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int build_srm_rpet_pWu_symZ(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int npix_exp2d,int fac_comp,int opc_subset)
{
char nom_fitxer[128],nom_fitxer_sb[128],nom_fitxer2[128];
int in,ind_proj,l,k2,k1,Ncol,nsb,Nb_os,ni;
double NITEMS_INI,NITEMS;
float z1,z2,angle;
float n0,n_i,angle0,z0,inc_zp,inc_ang,ang_i,inc_n,d2,sum;
float factor,zmin,zmax,tan_phi;
float epsilon=1e-4,dmax;
float xcon,ycon,zcon,rcon1,rcon2,xini,yini,zini;
sparse_i4 mp;
struct coord3d coord1;
struct coord3d coord2;
struct lor_params lor;
struct imagen ima;
int desf,con,con2;
desf=0;con2=0;
// params definition
strcpy(nom_fitxer,"sresmx");
strcpy(nom_fitxer2,"sresmx");
sprintf(nom_fitxer2,"%s.pesos",nom_fitxer);
Ncol=Nfil;
mp.n_fil=Nfil; // image pixels in X-Y-directions 
mp.n_col=Ncol; // image pixels in X-Y-directions 
mp.Nbins=(2*ndet)-1; // transaxial detectors 
mp.tbin=tdet/2.;
mp.tpix=mp.Nbins*mp.tbin/Nfil; // image pixel size in X-Y-directions in mm
mp.Ntallsima=2*ndet_z-1; // image pixel in Z-direction
if(axial_diff==0) mp.Ntallsima=ndet_z;  
mp.lvox=0.5*tdet_z; // image pixel size in Z-direction in mm 
if(axial_diff==0) mp.lvox=tdet_z;
mp.Nang=Nang;
mp.ndet_z=ndet_z;  
mp.tdet_z=tdet_z;
mp.ddet=ddet;
mp.Nsubsets=Nsubsets;
mp.axial_diff=axial_diff;
mp.inv=0;
mp.sym=111;
/*if(mp.sym==111) ni=((num_imagenes(ndet_z,axial_diff)-mp.ndet_z)/2)+mp.ndet_z;*/
ni=num_imagenes(ndet_z,axial_diff);
Nbt=mp.Nbins*Nang*ni;	// sinogram elements 
Nb_os=Nbt/Nsubsets; // sinogram elements per subset 
Nvox=Nfil*Ncol*mp.Ntallsima; // image elements 
d2=0.5*Nfil*mp.tpix; //medio FOV
if(fac_comp!=666) NITEMS_INI=(Nb_os/fac_comp)*Nvox; // matrix elements/fac_comp  
if(fac_comp==666) NITEMS_INI=Nb_os*Nvox; // no sparse matrix elements  
NITEMS=NITEMS_INI; // sparse matrix elements estimation  


// no sparse warning
if(fac_comp==666) printf("\n**\nWARNING! No sparse format was selected. You could have memory problems...\n**\n"); 


// steps
inc_zp=tdet_z; // axial increment in mm
inc_ang=180./Nang; // rotating in 180 degrees in degrees
inc_n=mp.tbin; // transaxial increment in mm 

// Origin coords 
z0=tdet_z/2.; // mm 
angle0=0.; // degrees 
n0=-0.5*((float)mp.Nbins-1.)*mp.tbin; // mm


// Display for selected params
printf("\n**************************************************************");
printf("\nVolumen de reconstrucción: %dx%dx%d (pixel:%f mm y corte:%f mm)\n",Nfil,Ncol,mp.Ntallsima,mp.tpix,mp.lvox);
printf("\nCaracterísticas del escáner: Detectores planos: %d bins (tamaño:%f mm) x %d ángulos\n",mp.Nbins,mp.tbin,Nang);
printf("Maximum axial difference is %d (Sinogramas totales: %d)\n",axial_diff,num_imagenes(ndet_z,axial_diff));
printf("%d detectores axiales de %f mm. Detectores separados: %f mm\n\n\n",ndet_z,tdet_z,ddet);
printf("\n**************************************************************");


// memory allocation
printf("\n\nInitial memory estimation to allocate per subset %1.0f Mb\n\n",NITEMS*sizeof(float)/1048576.);
if((mp.ar = (float *)calloc(NITEMS,sizeof(float)))==NULL) error_fesmp(53,NITEMS*sizeof(float)/1048576);
if((mp.ja = (int *)calloc(NITEMS,sizeof(int)))==NULL) error_fesmp(53,NITEMS*sizeof(int)/1048576);
mp.ia = (int *)calloc(Nbt+1,sizeof(int));

//lee_template_hdr("ptos_corte.hdr",&ima);
// SUBSET control  
for(nsb=0;nsb<Nsubsets;nsb++)
{
if(opc_subset==nsb || opc_subset==-1)
	{
		//inicializa mp y genera nombre matriz correspondiente al subset
		mp.ne=0;
		sprintf(nom_fitxer_sb,"%s.%d",nom_fitxer,nsb);
		printf("\n*************************\n\nSubset: %d *******************************\n",nsb);
			
		// AXIAL LOR control
		for(k1=0,z1=z0,ind_proj=0;k1<ndet_z;z1+=inc_zp,k1++)
		{
		for(k2=0,z2=z0;k2<ndet_z;z2+=inc_zp,k2++)
		{
			
			if(abs(k1-k2)>axial_diff)
			{
				continue;
			}
			
		
		
 
			// n-angular independent LOR params definition
			lor.z_med=fabs(z1+z2)/2;  
			zmax=maximo(z1,z2)+epsilon;
			zmin=minimo(z1,z2)-epsilon;
			printf("\nSinogram: k1=%d (%f) k2=%d (%f) (memory used: %f Mb) ",k1,z1,k2,z2,6.*mp.ne/1048576.);
			lor.phi=atan((z2-z1)/ddet); 
			tan_phi=tan(lor.phi);
			factor=cos(lor.phi);
			lor.mod_vlor=sqrt((fabs(z2-z1)*fabs(z2-z1))+(ddet*ddet));
			
			// ANGULAR LOR control
			for(l=nsb,angle=angle0+(nsb*inc_ang);l<Nang;l+=Nsubsets,angle+=inc_ang*Nsubsets)
			{
				// angular LOR params definition
				lor.seno=sin(angle*M_PI/180.);
				lor.coseno=cos(angle*M_PI/180.);
				con=0;
				
				// TRANSAXIAL LOR control
				for(in=0,lor.n=n0;in<mp.Nbins;in++,ind_proj++,lor.n+=inc_n)
				{
					
					// Control por si hay errores en el índice mp.ia (tamanho de la LOR)
					//printf("mp.ia[%d]:%d\n",ind_proj,mp.ia[ind_proj]);
					if(ind_proj>=2)
					{
						
						if(mp.ia[ind_proj-1]-mp.ia[ind_proj-2]<0 || mp.ia[ind_proj-1]-mp.ia[ind_proj-2]>mp.n_col*mp.n_fil*mp.Ntallsima) 
						{
							printf("\n\n\n¡Error inesperado (hay una LOR con %d elementos no nulos)!\n\n\n\n\n",mp.ia[ind_proj-1]-mp.ia[ind_proj-2]);
							exit(0);
						}
					}
					
 					if(k1>k2)
 					{
  						mp.ia[ind_proj]=mp.ne;
  						continue;
 					}
						

					// calcula contribución de LOR en otras posiciones de los detectores
					n_i=lor.n;
					sum=0.;
 					while(n_i<=0 && n_i>=n0) 
 					{
 						ang_i=atan((n_i-lor.n)/0.5*ddet);
 						n_i-=inc_n;
 						sum+=cos(ang_i*M_PI/180.);
 					}
 					n_i=lor.n;
 					while(n_i>0 && n_i<=-n0) 
 					{
 						ang_i=atan((n_i-lor.n)/0.5*ddet);
 						n_i+=inc_n;
 						sum+=cos(ang_i*M_PI/180.);
 					}
					//sum=1.;
/*					if(((in==29 && k1==0 && k2==34) || (in==29 && k1==34 && k2==0) || (in==29 && k1==17 && k2==17)) || ((in==29 && k1==10 && k2==24) || (in==29 && k1==24 && k2==10)))*/ 
					
					
// 					if((in==29 && k1==0 && k2==34) || (in==29 && k1==34 && k2==0) || (in==29 && k1==10 && k2==24) || (in==29 && k1==24 && k2==10))
// 					{
					//calculate non-null elements only
					mp.ia[ind_proj]=mp.ne; 
					calcul_ptos_corte(d2, lor.n, -lor.seno, lor.coseno, -tan_phi,lor.z_med, &coord1, &coord2);
					//printf("\nP1: %f,%f,%f --- P2: %f,%f,%f  LOR: n:%f z_med:%f",coord1.x,coord1.y,coord1.z,coord2.x,coord2.y,coord2.z,lor.n,lor.z_med);
					dmax=10000*maximo(fabs(coord1.x-coord2.x),fabs(coord1.y-coord2.y));
					 //if(k1==5 && k2==5){
						 
						if(con==0) 
						{ 
							xcon=coord1.x; ycon=coord1.y; zcon=coord1.z;
							if(con2!=0) { xcon=xini; ycon=yini; zcon=zini;}
						}
						
						//printf("Cerca? %f, %f, %f OR %f, %f, %f desde %f %f %f \n",coord1.x,coord1.y,coord1.z,coord2.x,coord2.y,coord2.z,xcon,ycon,zcon);		
						rcon1=((xcon-coord1.x)*(xcon-coord1.x))+((ycon-coord1.y)*(ycon-coord1.y))+((zcon-coord1.z)*(zcon-coord1.z));
						rcon2=((xcon-coord2.x)*(xcon-coord2.x))+((ycon-coord2.y)*(ycon-coord2.y))+((zcon-coord2.z)*(zcon-coord2.z));
						if(rcon1>rcon2)
						{
							xcon=coord2.x; ycon=coord2.y; zcon=coord2.z;
							coord2.x=coord1.x; coord2.y=coord1.y; coord2.z=coord1.z;
							coord1.x=xcon; coord1.y=ycon; coord1.z=zcon;
						}
						else
						{							
							xcon=coord1.x; ycon=coord1.y; zcon=coord1.z;
						}
						//printf("%f, %f, %f\n",coord1.x,coord1.y,coord1.z);						
						if(con==0) 
						{
							xini=coord1.x; yini=coord1.y; zini=coord1.z; //printf("*****\ninicial %f, %f, %f\n\n",xini,yini,zini);
							con=1;
					 	}
						con2=1;
						
    					ray_tracing_lor_pWu(&coord1,&coord2, &mp, Nfil, mp.Ntallsima, mp.tpix, mp.lvox, mp.tbin, dmax, ind_proj, Nb_os, fac_comp, sum,&lor, npix_exp2d,&ima,desf);
					// }
			
					// redimensiona memoria si es necesario
					if(mp.ne >= NITEMS-Nvox)
					{
						NITEMS+=NITEMS_INI;
						if((mp.ar = (float *)realloc(mp.ar,NITEMS*sizeof(float)))==NULL) 
						{ 
							error_fesmp(53,NITEMS*sizeof(float)/1048576.);
						}
						if((mp.ja = (int *)realloc(mp.ja,NITEMS*sizeof(int)))==NULL)
						{
							error_fesmp(53,NITEMS*sizeof(int)/1048576.);
						}
						printf("\nrealloc for %1.0f Mb\n",NITEMS*sizeof(float)/1048576.);
					}
						
				} // transaxial LOR control
			} // angular LOR control

		} // axial LOR control
		} // axial LOR control
	
	
		// Writting matrix for each subset
		mp.ia[Nb_os]=mp.ne; printf("\n\n\nNo-null elements: %d\n",mp.ne);
/*		for(i=0;i<100;i++) printf("mp.ar: %f\n",mp.ar[i]);*/
		guarda_pesos_xOS(&mp,nom_fitxer_sb);
	}
desf++;
}

guarda_imagen_hdr("ptos_corte.hdr",&ima);
return 0;
}



















// +++++++ 



//Error function 


// +++++++

void error_fesmp(int nerr, int arg)
{
switch(nerr)
{
case 30: printf("\n\nError opening ouput file\n\n"); break;
case 53: printf("\n\nError trying allocate %d Mb\n\n",arg); break;
default: printf("\n\nError in error function\n\n");
}

exit(0);
}










