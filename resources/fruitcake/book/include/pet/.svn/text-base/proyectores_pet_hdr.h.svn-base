#ifndef __PROYECTORES_PET_H
#define __PROYECTORES_PET_H

#include <util/utiles_vectores.h>
#include <pet/d_punto_lor.h>
#include <util/mp_i4.h>
#include <util/inpout.h>

/********** compara ***********/
void compara(float *dx,float *dy,float *dz,int *cas);

// igual que la anterior pero proyecta un objeto sin la EXP(-MU*X)
// se puede utilizar para obtener factores de atenuación (ver corr_ate_pet3d_hdr)
float proy_pet3d(struct coord3d *c1,struct coord3d *c2,struct imagen *ima,float dmax);

// Proyecta un LOR por ray tracing y asigna pesos desde LOI para construir matriz de pesos 
int ray_tracing_lor(struct coord3d *c1,struct coord3d *c2,sparse_i4 *mp, int Nfil, int Ntallsima, float tpix, float lvox, float dmax, int ind_proj, int Nb_os, int fac_comp, float sum);

//funcion empleada en la siguiente
int check_if_exist(sparse_i4 *mp, int mp_ne_ini,int indice,int indice_actual, int *s);

// coloca valor en la matriz de pesos
void	set_value_mp(sparse_i4 *mp,float sum,float dplor, int indice, float lim_exp);

// Expande en X-Y
void	expansion_X(int npix_exp2d, float next_x, float next_y, float next_z, float tpix, float lvox, struct lor_params *lor, int Nfil, int Ntallsima, int ip, int ind_map, int Npix, sparse_i4 *mp, float sum, float tbin,int *mask_ima);
void	expansion_Y(int npix_exp2d, float next_x, float next_y, float next_z, float tpix, float lvox, struct lor_params *lor, int Nfil, int Ntallsima, int jp, int ind_map, int Npix, sparse_i4 *mp, float sum, float tbin,int *mask_ima);


// Proyecta un LOR por ray tracing y asigna pesos desde LOI para construir matriz de pesos 
int ray_tracing_lor_pWu(struct coord3d *c1,struct coord3d *c2,sparse_i4 *mp, int Nfil, int Ntallsima, float tpix, float lvox, float tbin, float dmax, int ind_proj, int Nb_os, int fac_comp, float sum,struct lor_params *lor, int npix_exp2d, struct imagen *ima,int desf);

// específico para sub-sample
int ray_tracing_lor_sub(struct coord3d *c1,struct coord3d *c2,sparse_i4 *mp, int Nfil, int Ntallsima, float tpix, float lvox, float dmax, int ind_proj, int Nb_os, int fac_comp, float sum);

/***************************** Calcula puntos de corte *****************************************************************************************************/
void calcul_ptos_corte(float a,float n,float seno,float coseno,float tan_phi,float z_med,struct coord3d *c1,struct coord3d *c2);

/***********Transforma los ptos de corte calculados con la funcion anterior en los índices del voxel más próximo***************/
void ptos_corte_en_ijk(int *i,int *j,int *k,struct coord3d *c,float lvox,float tpix,int Nfild2);


/***********Transforma los ptos de corte calculados con la funcion anterior en los índices del voxel más próximo***************/
// cambios para ray tracing
void ptos_corte_en_ijk_rt(int *i,int *j,int *k,struct coord3d *c,float lvox,float tpix,int Nfild2);

// cambios para pWu fast
void ptos_corte_en_ijk_pWu_fast(int *i,int *j,int *k,struct coord3d *c,float lvox,float tpix,int Nfild2);

// calcula coordenadas del centro del voxel
void ptos_centro_voxel_pWu_fast(struct coord3d *cc,int i,int j,int k,float lvox,float tpix,int Nfild2);


/* Función que obtiene la matriz del sistema rPET útil para proyectar FACTORES DE ATENUACIÓN, 
Modificado por Pablo Aguiar - Mayo 2007*/
int build_srm_rpet_ray_tracing_att(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset);


/* Función que obtiene la matriz del sistema rPET, 
Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 */
int build_srm_rpet_ray_tracing(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset);

/* Función que obtiene la matriz del sistema rPET, 
Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 */
int build_srm_rpet_ray_tracing_sub(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset);

/* Función que obtiene la matriz del sistema rPET, 
Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 */
int build_srm_rpet_ray_tracing_sub_simple(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset);

/* Función que obtiene la matriz del sistema rPET, 
Pablo Aguiar - Valencia Octubre-Noviembre-Diciembre 2006 */
int build_srm_rpet_ray_tracing_sub_grid(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset);

// Funcion que obtiene matriz por el metodo de la distancia (pseudo-WU)
int build_srm_rpet_distance(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int fac_comp,int opc_subset,float pesminim);

// Implementacion del pseudo-Wu fast
int build_srm_rpet_ray_tracing_pWu(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int npix_exp2d,int fac_comp,int opc_subset);

// Implementacion del pseudo-Wu fast
int build_srm_rpet_pWu_symZ(int Nfil,int ndet_z,float tdet_z,int Nang,int ndet,float tdet,float ddet,int Nsubsets,int axial_diff,int npix_exp2d,int fac_comp,int opc_subset);

/* funcion error de build_srm_rpet_ray_tracing */
void error_fesmp(int nerr, int arg);

#endif /*__PROYECTORES_PET_H*/
