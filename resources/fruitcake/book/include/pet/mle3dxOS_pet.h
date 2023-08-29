#ifndef __MLE3DXOS_PET_H
#define __MLE3DXOS_PET_H

#include <util/mp_i4.h>

void mlemxOS_pet(float *q,int Nfil,int Ntalls,int Nvox,float *le,float *lc,int Nang,int Nbp,int NbOS,float *n_contes,int OS,char *fmatriu_pes,int NITER,int norma,char *fitxer_ima,char tipout[3]);

void rec_mlexOS_pet(float *q,float *q2,float *le,float *lc,float *vp,float *n_contes,int OS,int Nvox,int NbOS,int NangOS,int norma,char *fmatriu_pes,sparse_i4 *mp);

void numerocpxOS(float *le,float *n_comptes,int OS,int NbOS);

void check_nitems_max_mp_pet(char *fmatriu_pes,int OS,int *Nitems_max);


void error_rec(int nerr,char *text);

#endif //__MLE3DXOS_PET_H
