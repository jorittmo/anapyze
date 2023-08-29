
#ifndef __MP_I4__PET_H
#define __MP_I4__PET_H

#include <string.h>
#include <util/mp_i4.h>
#include <util/inpout.h>
#include <pet/mle3dxOS_pet.h> 

void fwd_proj(float *sino,float *ima,sparse_i4 *mp);

void fwd_proj_attenuation(float *sino,float *ima,sparse_i4 *mp);

int fwdproj_att(char *matriz_in,char *objeto_hdr);

void back_proj(float *lc,float *q2,sparse_i4 *mp);

void back_proj_sb(float *lc,float *q2,sparse_i4 *mp, int sb);

void lee_parametros_pesos(char *matriu_pes,sparse_i4 *mp);

void lee_pesos_xOS(sparse_i4 *mp,char *nom_fitxer);

void guarda_pesos_xOS(sparse_i4 *mp,char *nom_fitxer);

void integral_matrix_bins(sparse_i4 *mp, float  *vector);

void ima_sensibilidad(struct imagen *sens,sparse_i4 *mp);

void ima_sensibilidad_mp(struct imagen *sens,char *matriu_pes);

void error_mp_pet(int nerr,char *text);

#endif //__MP_I4_PET_H


