#ifndef __ATE3D_SPECT_H
#define __ATE3D_SPECT_H

/**************************************************************************
        Rutines per al calcul de l'atenuacio 3d.
        CONVENI D'INDEXS: x files, y columnes, z talls.
        ORDRE EN QUE ES MOUEN ELS INDEXS: 1er y-j, 2n x-i, 3er z-k
        ORIGEN DE COORDENADES: Al centre del pla del centre (en z) del primer tall

        Modificacions: mapa d'atenuació entrat de baix a dalt. Origen de
        les z al centre del tall d'abaix de tot (Carles 12-3-03)
***************************************************************************/

void compara(float *dx,float *dy,float *dz,int *cas);

float calcul_ate(float xpix,float ypix,float zpix,float xbin,float ybin,float zbin,float tpix,float lvox,float *mapa,float dmax,int Nfil,int Nfild2,int Npix,int Ntalls,int i,int j,int k);

#endif //__ATE3D_SPECT_H
