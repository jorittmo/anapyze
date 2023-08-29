#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <spect/ate3d_spect.h>
#define SIGNE(a) (a<-EPSILON?-1:(a>EPSILON?1:0))
#define EPSILON 1e-10


/**************************************************************************
        Rutines per al calcul de l'atenuacio 3d.
        CONVENI D'INDEXS: x files, y columnes, z talls.
        ORDRE EN QUE ES MOUEN ELS INDEXS: 1er y-j, 2n x-i, 3er z-k
        ORIGEN DE COORDENADES: Al centre del pla del centre (en z) del primer tall

        Modificacions: mapa d'atenuació entrat de baix a dalt. Origen de
        les z al centre del tall d'abaix de tot (Carles 12-3-03)
***************************************************************************/

void compara(float *dx,float *dy,float *dz,int *cas)
{
if(*dx==*dy)*dx+=EPSILON;
if(*dx==*dz)*dx+=EPSILON;
if(*dy==*dz)*dy+=EPSILON;

if(*dx<*dy){
  if(*dx<*dz) *cas=3;
  else *cas=1;
  }
else{
  if (*dy<*dz) *cas=2;
  else *cas=1;
  }

}

/*======================================== ATENUACIO =========================*/

float calcul_ate(float xpix,float ypix,float zpix,float xbin,float ybin,float zbin,float tpix,float lvox,float *mapa,float dmax,int Nfil,int Nfild2,int Npix,int Ntalls,int i,int j,int k)
{
float ate,mudx,vx,vy,vz,dpb,dant,dx,dy,dz,next_x,next_y,next_z,lvoxd2;
register int ip,jp,kp,incrx,incry,incrz,ind_map;
int cas;

/* Calcul del vector de la recta que uneix el voxel i el bin
   Calcul del increments (signes) dels diferents indexs (i,j,k) al anar del voxel al bin */

vx=xbin-xpix;
incrx=SIGNE(vx);
vy=ybin-ypix;
incry=SIGNE(vy);
vz=zbin-zpix;
incrz=SIGNE(vz);

/* Calcul del vector unitari i de la distancia */

dpb=sqrt(vx*vx+vy*vy+vz*vz);
vx/=dpb;
vy/=dpb;
vz/=dpb;
lvoxd2=lvox/2.;

/* calcul del seguent index (i,j,k) de la malla 0->N cadascuna de les N+1 separacions entre voxels
  que atravesara el raig en el seu cami cap el bin */

if(incrx<=0) ip=i;
else ip=i+1;
if(incry<=0) jp=j;
else jp=j+1;
if(incrz<=0) kp=k;
else kp=k+1;

 /* calcul de les coordenades del seguents plans (x,y,z) de separació entre voxels (malla)
  que atravesara el raig en el seu cami cap el bin */

next_x=(-Nfild2+ip)*tpix;
next_y=(-Nfild2+jp)*tpix;
next_z=((float)kp)*lvox-lvoxd2;

/* limitar la busqueda d'nterseccions dins del domini de la imatge. Evitar buscar interseccions de vectors
   paral.lels o quasi paral.lels  */

if(fabs(vx)>EPSILON) dx=(next_x-xpix)/vx;
else dx=dmax;
if(fabs(vy)>EPSILON) dy=(next_y-ypix)/vy;
else dy=dmax;
if(fabs(vz)>EPSILON) dz=(next_z-zpix)/vz;
else dz=dmax;


/* inicialitzacio de les variables. dant: distancia anterior al voxel, mudx: acumulat de la distancia
pel cada coeficient d'atenuacio, ind_map: index incial del mapa d'atenuacio (posició del voxel dins el mapa) */

dant=0;
mudx=0;
ind_map=k*Npix+i*Nfil+j;

/* bucle mentre estem dins del vomun (mapa d'atenuacio)
   compara les distàncies del voxel al proper pla x, y i z de la malla del mapa d'atenuació
   La que es menor marca el canvi d'index dins del mapa d'atenuacio.
   L'espai recorregut dins el voxel corresponent ve determinat per la diferencia entre
   aquesta distancia i la de l'anterior pla atravessat.
   Per evitar els rajos passant per les interseccions (dos indexs o tres variant a l'hora)
   s'hi i afegeix un valor EPSILON a una de les distancies en cas de coincidencia  */

while(ip>=0 && ip<=Nfil && jp>=0 && jp<=Nfil && kp>=0 && kp<=Ntalls){
  compara(&dx,&dy,&dz,&cas);
  switch(cas){
     case 1:
          mudx+=(dz-dant)* *(mapa+ind_map);
          dant=dz;
          kp+=incrz;
          ind_map+=incrz*Npix;
          next_z+=lvox*incrz;
          dz=(next_z-zpix)/vz;
          break;
     case 2:
          mudx+=(dy-dant)* *(mapa+ind_map);
          dant=dy;
          jp+=incry;
          ind_map+=incry;
          next_y+=tpix*incry;
          dy=(next_y-ypix)/vy;
	  break;
     case 3:
          mudx+=(dx-dant)* *(mapa+ind_map);
          dant=dx;
          ip+=incrx;
          ind_map+=incrx*Nfil;
          next_x+=tpix*incrx;
          dx=(next_x-xpix)/vx;
          break;
      default:
          printf("\nError al resultat de comparar d a ate3D.h\n");
	  exit(0);
          }
    }

ate=exp(-mudx);
return(ate);

}

