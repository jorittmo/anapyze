#ifndef __IO_IMA2_H
#define __IO_IMA2_H

/*******************************************************************************
*                                                                              * 
*          AQUESTA RUTINA LLEGEIX UNA IMATGE O COLECCIO D'IMATGES              *
*           EN QUALSEVOL FORMAT I HO RETORNA A UN PUNTER DOUBLE.               *
*                                                                              * 
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <util/io_ima2.h>

void error_io2(int nerr);
void input_imatge(int Mfil,int Mcol,int N,double *dq,int tipo,char *nom_fitxer);
void output_imatge(int Mfil,int Mcol,int N,double *dq,char tipo,char *nom_fitxer,double offset,double factor_escala);
void append_imatge(int Mfil,int Mcol,int N,double *dq,char tipo,char *nom_fitxer,double offset,double factor_escala);
void grabar_imatge(int Mfil,int Mcol,int N,double *dq,char tipo,FILE *file,double offset,double factor_escala);

#endif /*__IO_IMA2_H */
