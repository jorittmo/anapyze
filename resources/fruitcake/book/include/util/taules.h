#ifndef TAULES__H
#define TAULES__H

#include <stdio.h>
#include <stdlib.h>

/*=================================== LLEGIR FILA INT =====================*/
void llegir_fila_int(char *fitxer,int *fila,int ncol);
/*=================================== LLEGIR FILA FLOAT =====================*/
void llegir_fila_fl(char *fitxer,float *fila,int ncol);
/*=================================== LLEGIR TAULA FLOAT =====================*/
void llegir_taula_fl(char *fitxer,float **taula,int nfil,int ncol);
/*=================================== SUBTAULA_FILA ==========================*/
void subtaula_fila_fl(float **taula1,float **taula2,int nfil1,int nfil2,int ncol,int primera,int salt);
/*=================================== EXTREU_FILA FLOAT ======================*/
int extreu_fila_fl(float **taula,float *fila,int nfil,int ncol,int n);
/*=================================== SUBTAULA_COL ===========================*/
void subtaula_col_fl(float **taula1,float **taula2,int nfil,int ncol1,int ncol2,int *col);
/*=================================== ESCRIU TAULA FLOAT =====================*/
void escriu_taula_fl(char *fitxer,float **taula,int nfil,int ncol);
/*=================================== ESCRIU TAULA FLOAT TRASPOSTA ===========*/
void escriu_taula_fl_trasp(char *fitxer,float **taula,int nfil,int ncol);
/*=================================== ESCRIU DUES TAULES FLOAT ===============*/
void escriu_2_taules_fl(char *fitxer,float **taula1,int nfil,int ncol1,float **taula2,int ncol2,int pcol_t2);
/*=================================== ESCRIU FILA FLOAT ======================*/
void escriu_fila_fl(char *fitxer,float *fila,int ncol);
/*=================================== ESCRIU FILA FLOAT FIX3 =================*/
void escriu_fila_fl_fix3(char *fitxer,float *fila,int ncol);
/*=================================== ESCRIU FILA INT ======================*/
void escriu_fila_int(char *fitxer,int *fila,int ncol);
/*=========================================== ERROR =========================*/
void error_taules(int nerr,char *text);

#endif /*TAULES__H*/
