#ifndef __MLE3DXOS_H
#define __MLE3DXOS_H

#include <util/mp_i4.h>
void mlemxOS(float *q,int Nfil,int Ntalls,int Nvox,float *le,float *lc,int Nang,int Nbp,int NbOS,float *n_contes,
	int OS,char *fmatriu_pes,int NITER,int norma,char *fitxer_ima,char *tipout,float *le_scat);
void mlemxOS_mascara(float *q,int Nfil,int Ntalls,int Nvox,float *le,float *lc,int Nang,int Nbp,int NbOS,float *n_contes,
	int OS,char *fmatriu_pes,int NITER,int norma,char *fitxer_ima,char *tipout,float *le_scat,float num_pix_masc);
void rec_mlexOS(float *q,float *q2,float *le,float *lc,float *vp,float *n_contes,int OS,int Nvox,int NbOS,
	int NangOS,int norma,char *fmatriu_pes,sparse_i4 *mp, float *le_scat);
void numerocpxOS(float *le,float *n_comptes,int OS,int NbOS);
void check_nitems_max_mp(char *fmatriu_pes,int OS,int *Nitems_max);
void error_mle3dxOS(int nerr,const char *text);

#endif 
