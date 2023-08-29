/*#define MAX_SIGMES 3*/ /*VALOR INICIAL*/

/*PRUEBA DEL 8 FEBRERO DE 2007*/
/*VAMOS A DEFINIR MAX_SIGMAS COMO 2*/

#define MAX_SIGMES 2
#define EPSILON 1e-10
#define K 0.3989422863
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <spect/pes3d.h>

/*MODIFICACION: Cris mayo 2006*/

/*1.Las constantes (direccion fan beam y paralela) de los colimadores ELSCINT y GSK.La nomencaltura de las constantes es coherente con la que se utiliza en la definicion de los colimadores en simset v2.6)
2.Colimador case1, se incluye el detector
3.Colimador case2, no se incluye el detector, simulacion simset_frey*/

/*paralelo_elscint_L40mm, paralelo_elscint_L50mm; estos nombres quieren decir que los parametros del PSF_num son los del colimador fan beam de la elscint salvo la focal; es para hacer unas pruebas deterministas y ver cuanto se recupera según la PSF que se este considerando*/

enum llista_colimadors {
                 fan_beam_elscint,
                 fan_beam_elscint_simset,
                 paralel3_tesi_carles,
                 paralel4_tesi_carles,
                 paralel5_LEHR_siemens_ecam,
		 fan_beam_gsk_prism3000,
		 fan_beam_elscint_PSFinrinseca,
		 paralelo_elscint_L40mm,
                 paralelo_elscint_L50mm,
		 paralelo_ECAM_L4cm,
		 fan_beam_elscint_L24mm,
		 paralelo_ECAM_sinPSF,
                 no_colimador};

/*aixo es la funcio d'entrada del colimador*/
void parametres_colimador(tipo_colimador *COL)
{
printf("\n parametres col.limador: %d",COL->num);
switch(COL->num){
    case 1:  /*case ELSCINT (fan beam):*/
        COL->F=35.5;
        COL->L=4.;
        COL->A_Y=0.3369; /* constante direccion fan beam*/
        COL->A_Z=0.3369; /*constante direccion paralela*/
        COL->B=0.8;
        COL->w=0.0866;
        COL->sigma_int=0.17;
        COL->do_fanbeam=1;
        break;
   case 2: /*case ELSCINT (fan_beam_simset):*/
        COL->F=35.5;
        COL->L=4.;
        COL->A_Y=0.3369; /*constante direccion fan beam*/
        COL->A_Z=0.3369; /*constante direccion paralela*/
        COL->B=0.;
        COL->w=0.0866;
        COL->sigma_int=0.17;
        COL->do_fanbeam=1;
        break;
   case 3: /*paralel3:*/
        COL->A=0.0275;
        COL->B=0.2;
        COL->do_fanbeam=0;
        break;
   case 4: /*paralel4:*/
        COL->A=0.0172;
        COL->B=0.2;
        COL->do_fanbeam=0;
        break;
   case 5: /*paralel5: (Es la ECAM parametros de la regresion de los datos calculados con PENELOPE (27/07/2006))*/
        COL->A=0.0167;
        COL->B=0.1405;
        COL->do_fanbeam=0;
        break;
   case 6: /*fan_beam_gsk_prism3000: resultats Albert treball 002 GSK_project 21/2/06*/
        COL->F=65.0;
        COL->L=2.7;
        COL->A_Y=0.3575; /*constante direccion fan beam*/
        COL->A_Z=0.3360; /*constante direccion paralela*/
        COL->B=0.0;
        COL->w=0.0866;
        COL->sigma_int=0.17;
        COL->do_fanbeam=1;
        break;
   case 7: /*fan_beam (ELSCINT): constantes muy pequeñas para que PSF sea del orden de la intrinseca*/
        COL->F=35.5;
        COL->L=4.0;
        COL->A_Y=0.00001; /*constante direccion fan beam*/
        COL->A_Z=0.00001; /*constante direccion paralela*/
        COL->B=0.8;
        COL->w=0.0866;
        COL->sigma_int=0.17;
        COL->do_fanbeam=1;
        break;	
   case 8: /*paralelo: (Los parametros son los mismos que el colimador fan beam de la elscinst, F=355000000 mm y la L=40 mm*/ /*noviembre 2006 para hacer unas pruebas deterministas*/
        COL->A=0.0165;
        COL->B=0.093;
        COL->do_fanbeam=0;
        break;
   case 9: /*paralelo: (Los parametros son los mismos que el colimador fan beam de la elscinst, F=355000000 mm y la L=50 mm*/ /*noviembre 2006 para hacer unas pruebas deterministas*/
        COL->A=0.0132;
        COL->B=0.0976;
        COL->do_fanbeam=0;
        break;
   case 10: /*paralelo: (Los parametros son los mismos que el colimador paralelo de la ECAM, pero se pone una L que es 40 mm (que es la L del colimador fan beam de la elscinst*/ /*noviembre 2006 para hacer unas pruebas deterministas*/
        COL->A=0.0101;
        COL->B=0.0998;
        COL->do_fanbeam=0;
        break;
   case 11:  /*case ELSCINT (fan beam) con una L de 2,405 cm y que es la L del colimador paralelo de la ECAM:*/
        COL->F=35.5;
        COL->L=2.405;
        COL->A_Y=0.3369; /* constante direccion fan beam*/
        COL->A_Z=0.3369; /*constante direccion paralela*/
        COL->B=0.8;
        COL->w=0.0866;
        COL->sigma_int=0.17;
        COL->do_fanbeam=1;
	break; 
   case 12: /*paralelo ECAM: PSF es del orden de la intriseca y constante con z*/
        COL->A=0.0;
        COL->B=0.158;
        COL->do_fanbeam=0;
        break;

   default:
        printf("\n\nrutina pes3d.h: Error al nombre de colimador\n\n");
        exit(0);
   }

}

/*========================================= PES ==============================*/
/*calcul pes fan beam eix x o eix fan beam*/
float calcul_pes_fb(float dlat,float dpp,float *gauss,float costheta,float tpix,float tpixd2,float tbind2,tipo_colimador COL)
{
float denom,xc,pes,sigma;

/*formula abast maxim: z i focal referides al front plane*/
denom=sqrt(COL.L*COL.L*(COL.F-dpp)*(COL.F-dpp)-COL.w*COL.w*(COL.L+2.*dpp)*(COL.L+2.*dpp));
xc=COL.A_Y*(dpp+COL.L+COL.B)*COL.w*(2.*COL.F+COL.L)/costheta/denom;
sigma=sqrt(COL.sigma_int*COL.sigma_int+xc*xc);

pes=calcul_pes(dlat/sigma,tbind2/sigma,gauss);

/*eficiencia del col.limador fan-beam no es constant! valors relatius*/
/*pes*=cos(theta)*cos(theta)*(FOCAL+L)/(FOCAL-z);*/
pes*=costheta*costheta*(COL.F-5.)/(COL.F-dpp);	/*canvi 17.02.01, relativa a z=5cm */

return(pes);

}

/*========================================= PES fb z ==============================*/
/*calcul pes fan beam eix z o eix parallel*/
float calcul_pes_fb_z(float dlat,float dpp,float *gauss,float tpix,float tpixd2,float tbind2,tipo_colimador COL)
{
float xc,sigma,pes;

xc=2.*COL.A_Z*COL.w*(dpp+COL.L+COL.B)/COL.L;
sigma=sqrt(COL.sigma_int*COL.sigma_int+xc*xc);

pes=calcul_pes(dlat/sigma,tbind2/sigma,gauss);

/*pes*=costheta*costheta*(COL.F-5.)/(COL.F-dpp);*/

return(pes);

}

/*========================================= PES PAR ==============================*/
/*calcul pes parallel en x i z*/
float calcul_pes_par(float dlat,float dpp,float *gauss,float tpix,float tpixd2,float tbind2,tipo_colimador COL)
{
float sigma,pes;

sigma=COL.A*dpp+COL.B;
/*la distancia del punt al pixel dpp va en cm*/
/*sigma es en cm!!!!*/
/*A no te unitats pero la multipliques per cm i B en cm*/

pes=calcul_pes(dlat/sigma,tbind2/sigma,gauss);

return(pes);

}

/*========================================= PES ==============================*/
float calcul_pes(float dlat_norm,float tbind2_norm,float *gauss)
{
float pes,x,y;
int j1,j2;

if(dlat_norm>=tbind2_norm){
  x=dlat_norm-tbind2_norm;
  if(x>MAX_SIGMES) return 0.;
  y=dlat_norm+tbind2_norm;
  j1=floor(x/DX_GAUSS+EPSILON);
  j2=floor(y/DX_GAUSS+EPSILON);
  pes=fabs(gauss[j2]-gauss[j1]);
  }
else{
  x=tbind2_norm-dlat_norm;
  y=dlat_norm+tbind2_norm;
  j2=floor(y/DX_GAUSS+EPSILON);
  j1=floor(x/DX_GAUSS+EPSILON);
  pes=fabs(gauss[j1]+gauss[j2]-1.);
  }

return(pes);

}

/*==================================== GAUSSIANA =============================*/

void calcul_gaussiana(float *gauss)
{
int i;
float x,g,anterior;

gauss[0]=anterior=0.5;
x=-DX_GAUSS/2.;
for(i=1;i<PUNTS_GAUSS;i++){
  x+=DX_GAUSS;
  g=K*exp(-x*x/2.);
  gauss[i]=anterior+g*DX_GAUSS;
  anterior=gauss[i];
  }
}


// Funciones que incluyen nmax_sg febrero 2007 (Cris, Carles, Judith, Albert)

/*========================================= PES ==============================*/
/*calcul pes fan beam eix x o eix fan beam*/
float calcul_pes_fb_g(float dlat,float dpp,float *gauss,float costheta,float tpix,float tpixd2,float tbind2,tipo_colimador COL,float nmax_sg)
{
float denom,xc,pes,sigma;

/*formula abast maxim: z i focal referides al front plane*/
denom=sqrt(COL.L*COL.L*(COL.F-dpp)*(COL.F-dpp)-COL.w*COL.w*(COL.L+2.*dpp)*(COL.L+2.*dpp));
xc=COL.A_Y*(dpp+COL.L+COL.B)*COL.w*(2.*COL.F+COL.L)/costheta/denom;
sigma=sqrt(COL.sigma_int*COL.sigma_int+xc*xc);

pes=calcul_pes_g(dlat/sigma,tbind2/sigma,gauss,nmax_sg);

/*eficiencia del col.limador fan-beam no es constant! valors relatius*/
/*pes*=cos(theta)*cos(theta)*(FOCAL+L)/(FOCAL-z);*/
pes*=costheta*costheta*(COL.F-5.)/(COL.F-dpp);	/*canvi 17.02.01, relativa a z=5cm */

return(pes);

}

/*========================================= PES fb z ==============================*/
/*calcul pes fan beam eix z o eix parallel*/
float calcul_pes_fb_z_g(float dlat,float dpp,float *gauss,float tpix,float tpixd2,float tbind2,tipo_colimador COL,float nmax_sg)
{
float xc,sigma,pes;

xc=2.*COL.A_Z*COL.w*(dpp+COL.L+COL.B)/COL.L;
sigma=sqrt(COL.sigma_int*COL.sigma_int+xc*xc);

pes=calcul_pes_g(dlat/sigma,tbind2/sigma,gauss,nmax_sg);

/*pes*=costheta*costheta*(COL.F-5.)/(COL.F-dpp);*/

return(pes);

}

/*========================================= PES PAR ==============================*/
/*calcul pes parallel en x i z*/
float calcul_pes_par_g(float dlat,float dpp,float *gauss,float tpix,float tpixd2,float tbind2,tipo_colimador COL,float nmax_sg)
{
float sigma,pes;

sigma=COL.A*dpp+COL.B;
/*la distancia del punt al pixel dpp va en cm*/
/*sigma es en cm!!!!*/
/*A no te unitats pero la multipliques per cm i B en cm*/

pes=calcul_pes_g(dlat/sigma,tbind2/sigma,gauss,nmax_sg);

return(pes);

}

/*========================================= PES ==============================*/

float calcul_pes_g(float dlat_norm,float tbind2_norm,float *gauss,float nmax_sg)
{
float pes,x,y;
int j1,j2;

if(dlat_norm>=tbind2_norm){
  x=dlat_norm-tbind2_norm; /*7 febrero 2007 tiene  que estar dentro de nmax_sg el extremo superior del bin. La primera que vez que cambiamos esto, pusimos que la y es la que tiene que caer dentro de nmax_sig; esto hacia que las matrices calculadas tuvieran un tamaño menor*/
  y=dlat_norm+tbind2_norm;
  
  if(x>nmax_sg) return(0.); /*7 febrero 2007, ponemos el control de si hay que calcular o no en el valor de la distancia x. Asi esta igual que la calcul_pes*/

  j1=floor(x/DX_GAUSS_G+EPSILON);
  j2=floor(y/DX_GAUSS_G+EPSILON);
  pes=fabs(gauss[j2]-gauss[j1]);
  }
else{
  x=tbind2_norm-dlat_norm;
  y=dlat_norm+tbind2_norm;
  j2=floor(y/DX_GAUSS_G+EPSILON);
  j1=floor(x/DX_GAUSS_G+EPSILON);
  pes=fabs(gauss[j1]+gauss[j2]-1.);
  }

return(pes);

}

/*==================================== GAUSSIANA =============================*/
void calcul_gaussiana_g(float *gauss,float nmax_sg,int npunts_g,float *alturagaussiana)
{
int i;
float x,g,anterior;

gauss[0]=anterior=0.5;
x=-DX_GAUSS_G/2.;

for(i=1;i<npunts_g;i++)
{
  x+=DX_GAUSS_G;
  g=K*exp(-x*x/2.);
  gauss[i]=anterior+g*DX_GAUSS_G;
  anterior=gauss[i];
}
 *alturagaussiana=g;
}



