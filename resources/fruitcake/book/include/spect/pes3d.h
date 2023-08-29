#ifndef __PES3D_H
#define __PES3D_H
#define DX_GAUSS 0.05
#define DX_GAUSS_G 0.015

#define PUNTS_GAUSS 200

/*prototipus del pes3d.c llibreria de definicio de colimadors del femsp3d*/

typedef struct
{
  float A,A_Y,A_Z,B,F,w,L,sigma_int;
  int num,do_fanbeam;
}tipo_colimador;

void parametres_colimador(tipo_colimador *COL);

float calcul_pes_fb(float dlat,float dpp,float *gauss,float costheta,float tpix,float tpixd2,float tbind2,tipo_colimador COL);

float calcul_pes_fb_z(float dlat,float dpp,float *gauss,float tpix,float tpixd2,float tbind2,tipo_colimador COL);

float calcul_pes_par(float dlat,float dpp,float *gauss,float tpix,float tpixd2,float tbind2,tipo_colimador COL);

float calcul_pes(float dlat_norm,float tbind2_norm,float *gauss);

void calcul_gaussiana(float *gauss);


// Funciones que incluyen nmax_sg febrero2007 (Cris, Carles, Judith, Albert)

float calcul_pes_fb_g(float dlat,float dpp,float *gauss,float costheta,float tpix,float tpixd2,float tbind2,tipo_colimador COL,float nmax_sg);

float calcul_pes_fb_z_g(float dlat,float dpp,float *gauss,float tpix,float tpixd2,float tbind2,tipo_colimador COL,float nmax_sg);

float calcul_pes_par_g(float dlat,float dpp,float *gauss,float tpix,float tpixd2,float tbind2,tipo_colimador COL,float nmax_sg);

float calcul_pes_g(float dlat_norm,float tbind2_norm,float *gauss,float nmax_sg);

void calcul_gaussiana_g(float *gauss,float nmax_sg,int npunts_g,float *pesminim);


#endif //__PES3D_H
