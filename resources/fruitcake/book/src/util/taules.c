#include <stdio.h>
#include <stdlib.h>
#include <util/taules.h>

/*=================================== LLEGIR FILA INT =====================*/

void llegir_fila_int(char *fitxer,int *fila,int ncol)
{
register int j;
FILE *file;

if((file=fopen(fitxer,"r"))==NULL) error_taules(51,fitxer);
for(j=0;j<ncol;j++) fscanf(file,"\n%d",fila+j);
fclose(file);
}

/*=================================== LLEGIR FILA FLOAT =====================*/

void llegir_fila_fl(char *fitxer,float *fila,int ncol)
{
register int j;
FILE *file;

if((file=fopen(fitxer,"r"))==NULL) error_taules(51,fitxer);
for(j=0;j<ncol;j++) fscanf(file,"\n%f",fila+j);
fclose(file);
}

/*=================================== LLEGIR TAULA FLOAT =====================*/

void llegir_taula_fl(char *fitxer,float **taula,int nfil,int ncol)
{
register int i,j;
FILE *file;

if((file=fopen(fitxer,"r"))==NULL) error_taules(61,fitxer);
for(i=0;i<nfil;i++) for(j=0;j<ncol;j++) fscanf(file,"\n%f",&taula[i][j]);
fclose(file);
}

/*=================================== SUBTAULA_FILA ==========================*/

void subtaula_fila_fl(float **taula1,float **taula2,int nfil1,int nfil2,int ncol,int primera,int salt)
{
register int i,j,k;

if((primera+salt*(nfil2-1))>=nfil1) error_taules(42,"");
for(j=0,k=primera;j<nfil2;j++,k+=salt) for(i=0;i<ncol;i++) taula2[j][i]=taula1[k][i];
}

/*=================================== EXTREU_FILA FLOAT ======================*/

int extreu_fila_fl(float **taula,float *fila,int nfil,int ncol,int n)
{
register int i;

if(n>nfil) error_taules(45,"");
for(i=0;i<ncol;i++) fila[i]=taula[n][i];
return 1;
}

/*=================================== SUBTAULA_COL ===========================*/

void subtaula_col_fl(float **taula1,float **taula2,int nfil,int ncol1,int ncol2,int *col)
{
register int i,j;

for(j=0;j<nfil;j++) for(i=0;i<=ncol2;i++) taula2[j][i]=taula1[j][col[i]];
}

/*=================================== ESCRIU TAULA FLOAT =====================*/

void escriu_taula_fl(char *fitxer,float **taula,int nfil,int ncol)
{
register int i,j;
FILE *file;

if((file=fopen(fitxer,"w"))==NULL) error_taules(71,fitxer);
for(i=0;i<nfil;i++){
  for(j=0;j<ncol;j++) fprintf(file,"%f\t",taula[i][j]);
  fprintf(file,"%c\n",' ');
  }
fclose(file);
}

/*=================================== ESCRIU TAULA FLOAT TRASPOSTA ===========*/

void escriu_taula_fl_trasp(char *fitxer,float **taula,int nfil,int ncol)
{
register int i,j;
FILE *file;

if((file=fopen(fitxer,"w"))==NULL) error_taules(71,fitxer);
 for(j=0;j<ncol;j++){
  for(i=0;i<nfil;i++) fprintf(file,"%f\t",taula[i][j]);
  fprintf(file,"%c\n",' ');
  }
fclose(file);
}

/*=================================== ESCRIU DUES TAULES FLOAT ===============*/

void escriu_2_taules_fl(char *fitxer,float **taula1,int nfil,int ncol1,float **taula2,int ncol2,int pcol_t2)
{
register int i,j;
FILE *file;

/* pcol_t2=primera columna taula 2 (pcol_t2=1 permet no tornar a posar les x) */

if((file=fopen(fitxer,"w"))==NULL) error_taules(71,fitxer);
for(i=0;i<nfil;i++){
  for(j=0;j<ncol1;j++) fprintf(file,"%f\t",taula1[i][j]);
  for(j=pcol_t2;j<ncol2;j++) fprintf(file,"%f\t",taula2[i][j]);
  fprintf(file,"%c\n",' ');
  }
fclose(file);
}

/*=================================== ESCRIU FILA FLOAT ======================*/

void escriu_fila_fl(char *fitxer,float *fila,int ncol)
{
register int j;
FILE *file;

if((file=fopen(fitxer,"w"))==NULL) error_taules(81,fitxer);
for(j=0;j<ncol;j++) fprintf(file,"%f\t",fila[j]);
fprintf(file,"%c\n",' ');

fclose(file);
}

/*=================================== ESCRIU FILA FLOAT FIX3 =================*/

void escriu_fila_fl_fix3(char *fitxer,float *fila,int ncol)
{
register int j;
FILE *file;

if((file=fopen(fitxer,"w"))==NULL) error_taules(81,fitxer);
for(j=0;j<ncol;j++) fprintf(file,"%.3f\t",fila[j]);
fprintf(file,"%c\n",' ');

fclose(file);
}

/*=================================== ESCRIU FILA INT ======================*/

void escriu_fila_int(char *fitxer,int *fila,int ncol)
{
register int j;
FILE *file;

if((file=fopen(fitxer,"w"))==NULL) error_taules(81,fitxer);
for(j=0;j<ncol;j++) fprintf(file,"%d\t",fila[j]);
fprintf(file,"%c\n",' ');

fclose(file);
}

/*=========================================== ERROR =========================*/

void error_taules(int nerr,char *text)
{

switch(nerr){
    case 42:  printf("\n\nERROR al seleccionar taula:\n\t(primera+salt*nfil2) s'excedeix del nombre total de files disponibles\n"); break;
    case 45:  printf("\n\nERROR al seleccionar fila al separar:\n\tn s'excedeix del nombre total de files disponibles\n"); break;
    case 51:  printf("\nError al obrir el fitxer de la fila: el fitxer %s no existeix o no es pot obrir\n",text); break;
    case 61:  printf("\nError al obrir el fitxer de la taula: el fitxer %s no existeix o no es pot obrir\n",text); break;
    case 71:  printf("\nError al obrir el fitxer d'escriptura per la taula:\n el fitxer %s no es pot obrir\n",text); break;
    case 81:  printf("\nError al obrir el fitxer d'escriptura per la fila:\n el fitxer %s no es pot obrir\n",text); break;
    default: printf("\n\nError del numero d'error en la funcio error_taules()\n"); 
    }
    
exit(0);
}

