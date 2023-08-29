#ifndef CRV__H
#define CRV__H

#include <stdio.h>
#include <math.h>
#include <util/inpout.h>
#define MO 714025
#define IA 1366
#define IC 150889

extern int MM,PN,M,N,P;

void proj_segones(float *le,float *lep,float *les,char CV[80],int inisor,char tipoin[3],float nc[1]);

void cros_val(float *lcp,float *cs,float *lep,float *les,double *acpl,double *adpl,double *acsl,double *adsl,FILE *file2,int iteracio);

void rifa_nou(float *le,float *lp,float *ls,int i0);

void rifa(float *le,float *lp,float *ls,int i0);

void rifa1(float *le,float *lp,int i0);

double logfact(double num);

float ran2(long *idum);

#undef MO
#undef IA
#undef IC

#endif //__CRV_H
