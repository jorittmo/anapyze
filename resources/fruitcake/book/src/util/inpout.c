#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <util/inpout.h>
#include <util/nrutil.h>
#include <util/dbh.h>
#include <util/manejo_strings.h>
#define EOS '\0'
#define maxi(a,b) ((a)>(b) ? (a):(b))
#define mini(a,b) ((a)<(b) ? (a):(b))

char tokensep[]=" \t\n,";
/*===================================== TOKEN DADES FITXER IMATGES ===========*/

void token_dades_fitxer_ima(char *fitxerdades,char *nomfitxer,int *nfil, int *ncol, int *nima,char *tipoin)
{
char buf[240],*token;
FILE *file;

if((file=fopen(fitxerdades,"r"))==NULL) error_io(102,fitxerdades);
token=(char *)calloc(240,sizeof(char));

if(fgets(buf,80,file)==NULL) error_io(33,fitxerdades);
token=strtok(buf,tokensep);
strcpy(nomfitxer,token);
if(fgets(buf,80,file)==NULL) error_io(33,fitxerdades);
token=strtok(buf,tokensep);
*nfil=atoi(token);
if(fgets(buf,80,file)==NULL) error_io(33,fitxerdades);
token=strtok(buf,tokensep);
*ncol=atoi(token);
if(fgets(buf,80,file)==NULL) error_io(33,fitxerdades);
token=strtok(buf,tokensep);
*nima=atoi(token);
if(fgets(buf,80,file)==NULL) error_io(33,fitxerdades);
token=strtok(buf,tokensep);
strncpy(tipoin,token,2);

fclose(file);
}

/*=================================== LLEGIR FITXER ==========================*/

void input_fitxer(int nfil,int ncol,int nima,float *le,char *nom_fitxer,char tipoin[2])
{

int nitems,j;
unsigned short int *sq;
unsigned char *cq;
FILE *file;

nitems=nfil*ncol*nima;

if( (file=fopen(nom_fitxer,"r")) ==NULL) error_io(102,nom_fitxer);

if(strncmp(tipoin,"1b",2)==0 || strncmp(tipoin,"1B",2)==0){
       cq=ucvector(0,nitems-1);
       fread(cq,1,nitems,file);
       for(j=0;j<nitems;j++) *(le+j)=(float) *(cq+j);
       free_ucvector(cq,0,nitems-1);
       }
else if(strncmp(tipoin,"2b",2)==0){  
       sq=usvector(0,nitems-1);
       fread(sq,2,nitems,file);
       for(j=0;j<nitems;j++) *(le+j)=(float) *(sq+j);
       free_usvector(sq,0,nitems-1);
       }
else if(strncmp(tipoin,"fl",2)==0) fread(le,4,nitems,file);
else error_io(103,tipoin);
            
fclose(file);
}            

/*=================================== LLEGIR FITXER con offset paguiar05==========================*/

void input_fitxer_offset(int nfil,int ncol,int nima,float *le,char *nom_fitxer,char tipoin[2],int offset)
{

int nitems,j;
unsigned short int *sq;
short int *sqq;
signed int *iq;
double *dq;
unsigned char *cq;
FILE *file;

nitems=nfil*ncol*nima;

if( (file=fopen(nom_fitxer,"r")) ==NULL) error_io(102,nom_fitxer);

fseek(file,offset,SEEK_SET);

if(strncmp(tipoin,"1b",2)==0 || strncmp(tipoin,"1B",2)==0){
       cq=ucvector(0,nitems-1);
       fread(cq,1,nitems,file);
       for(j=0;j<nitems;j++) *(le+j)=(float) *(cq+j);
       free_ucvector(cq,0,nitems-1);
       }
else if(strncmp(tipoin,"2b",2)==0){  
       sq=usvector(0,nitems-1);
       fread(sq,2,nitems,file);
       for(j=0;j<nitems;j++) *(le+j)=(float) *(sq+j);
       free_usvector(sq,0,nitems-1);
       }
else if(strncmp(tipoin,"i2",2)==0){  
	  sqq=sivector(0,nitems-1);
	  fread(sqq,2,nitems,file);
	  for(j=0;j<nitems;j++) *(le+j)=(float) *(sqq+j);
	  free_sivector(sqq,0,nitems-1);
	  }
else if(strncmp(tipoin,"i4",2)==0){  
	  iq=ivector(0,nitems-1);
	  fread(iq,4,nitems,file);
	  for(j=0;j<nitems;j++) *(le+j)=(float) *(iq+j);
	  free_ivector(iq,0,nitems-1);
	  }
else if(strncmp(tipoin,"fl",2)==0) fread(le,4,nitems,file);

else if(strncmp(tipoin,"db",2)==0) {
	  dq=dvector(0,nitems-1);
	  fread(dq,8,nitems,file);
	  for(j=0;j<nitems;j++) *(le+j)=(float) *(dq+j);
	 }
          

else error_io(103,tipoin);
fclose(file);
}

/*=================================== LLEGIR FITXER SOBRE UNSiGNED CHAR =============*/

void input_fitxer_uchar(int nfil,int ncol,int nima,unsigned char *ue,char *nom_fitxer,char tipoin[2])
{
int nitems;
FILE *file;

nitems=nfil*ncol*nima;

if( (file=fopen(nom_fitxer,"r")) ==NULL) error_io(102,nom_fitxer);

if(strncmp(tipoin,"1b",2)==0)fread(ue,1,nitems,file);
else error_io(108,tipoin);

fclose(file);
}        
/*==================== LLEGIR LLESCA O SUMA DE LLESQUES D'UNA SEQUENCIA ======*/

void input_llesca(int nfil,int ncol,int nima,float *le,int p,int sumar,char *nom_fitxer,char tipoin[2])
{

int i,k,incol,itima,pncol,j,nitems,tima;
unsigned short int *sq;
unsigned char *cq;
float *fq;
FILE *file;

tima=nfil*ncol;
nitems=tima*nima;
pncol=p*ncol;
for(i=0;i<nitems;i++) *(le+i)=0;
 
if( (file=fopen(nom_fitxer,"r")) ==NULL) error_io(105,nom_fitxer);

if(strncmp(tipoin,"1b",2)==0){
       cq=ucvector(0,nitems-1);
       fread(cq,1,nitems,file);
       p--;
       for (k=0;k<sumar;k++){
          pncol=p*ncol;
          for(i=0;i<nima;i++){
            incol=i*ncol;
            itima=i*tima;
            for(j=0;j<ncol;j++) *(le+incol+j)+=(float) *(cq+itima+pncol+j);
            }
          p++;
          }  
       free_ucvector(cq,0,nitems-1);
       }
else if(strncmp(tipoin,"2b",2)==0){  
       sq=usvector(0,nitems-1);
       fread(sq,2,nitems,file);
       p--;
       for(k=0;k<sumar;k++){
          pncol=p*ncol;
          for(i=0;i<nima;i++){
             incol=i*ncol;
             itima=i*tima;
             for(j=0;j<ncol;j++) *(le+incol+j)+=(float) *(sq+itima+pncol+j);
             }
          p++;
          }
       free_usvector(sq,0,nitems-1);
       }
else if(strncmp(tipoin,"fl",2)==0){
       fq=vector(0,nitems-1);
       fread(fq,4,nitems,file);
       p--;
       for(k=0;k<sumar;k++){
          pncol=p*ncol;
          for(i=0;i<nima;i++){
             incol=i*ncol;
             itima=i*tima;
             for(j=0;j<ncol;j++) *(le+incol+j)+= *(fq+itima+pncol+j);
             }
          p++;
          }  
       free_vector(fq,0,nitems-1);
       }
else error_io(104,tipoin);
            
fclose(file);
}                 

/*=============================== INPUT_SELECCIO_IMA_DE_SEQ ==================*/

void input_selecc_ima_de_seq(char *fitxer,char tinp[3],float *seqima,int M,int N,int primera,int nima,int salt)
{
int MM;   
float *seqima2;

MM=M*M;
seqima2=vector(0,MM*N-1);
input_fitxer(M,M,N,seqima2,fitxer,tinp);
if((primera+salt*(nima-1))>=N) error_io(41,"");
subseq(seqima,seqima2,MM,primera,nima,salt);
free(seqima2);
}   

/*========================================= SUBSEQUENCIA =====================*/

void subseq(float *seqima,float *seqima2,int MM,int primera,int nima,int salt)
{
int k,i,j,kMM,jMM,MMs;   

MMs=MM*salt;
for(k=primera,j=0,kMM=MM*primera,jMM=0; k<nima ; k++,kMM+=MMs,j++,jMM+=MM)
                     for(i=0;i<MM;i++) *(seqima+jMM+i)=*(seqima2+kMM+i);
}

/*====================================== E S C R I T U R A ===================*/

void output_ima(int n_fil,int n_col,float *q,FILE *file,char tipout[3])
{

float *qN;
int i,n_items;
unsigned short int *cq;
unsigned char *uq;

n_items=n_fil*n_col;

if(strncmp(tipout,"1b",2)==0){
     uq=(unsigned char *)malloc(n_items);
     qN=(float *)malloc(n_items*4);
     normalitzacio_imatge(q,qN,n_items);
     for(i=0;i<n_items;i++) *(uq+i)=(unsigned char) *(qN+i);
     fwrite(uq,1,n_items,file);
     free(qN);
     free(uq);
     }
else if(strncmp(tipout,"2b",2)==0){
     cq=(unsigned short int *)malloc(n_items*2);
     for(i=0;i<n_items;i++) *(cq+i)=(unsigned short int) *(q+i);
     fwrite(cq,2,n_items,file);
     free(cq);
     }
else if(strncmp(tipout,"fl",2)==0)  fwrite(q,sizeof(float),n_items,file);
else error_io(101,tipout);

}

/*====================================== E S C R I T U R A ===================*/

void output_ima3(int n_fil,int n_col,int SL,float *q,FILE *file,char *tipout)
{

float *qN;
int i,n_items;
unsigned short int *cq;
short int *cqq;
int *iq;
unsigned char *uq;

n_items=n_fil*n_col*SL;

if(strncmp(tipout,"1b",2)==0){ /* Normalizacion a 255 */
     uq=(unsigned char *)malloc(n_items*sizeof(unsigned char));
     qN=(float *)malloc(n_items*sizeof(float));
     normalitzacio_imatge(q,qN,n_items);
     for(i=0;i<n_items;i++) *(uq+i)=(unsigned char) *(qN+i);
     fwrite(uq,sizeof(unsigned char),n_items,file);
     free(qN);
     free(uq);
     }
else if(strncmp(tipout,"1B",2)==0){ /* No se normalizan a 255 */
     uq=(unsigned char *)malloc(n_items*sizeof(unsigned char));
     for(i=0;i<n_items;i++) *(uq+i)=(unsigned char) *(q+i);
     fwrite(uq,sizeof(unsigned char),n_items,file);
     free(uq);
     }
else if(strncmp(tipout,"2b",2)==0){
     cq=(unsigned short int *)malloc(n_items*sizeof(unsigned short int));
     for(i=0;i<n_items;i++) *(cq+i)=(unsigned short int) *(q+i);
     fwrite(cq,sizeof(unsigned short int),n_items,file);
     free(cq);
     }
else if(strncmp(tipout,"i2",2)==0){
	cqq=(short int *)malloc(n_items*sizeof(short int));
	for(i=0;i<n_items;i++) *(cqq+i)=(short int) *(q+i);
	fwrite(cqq,sizeof(short int),n_items,file);
	free(cqq);
	}
else if(strncmp(tipout,"i4",2)==0){
	iq=(int *)malloc(n_items*sizeof(int));
	for(i=0;i<n_items;i++) *(iq+i)=(short int) *(q+i);
	fwrite(iq,sizeof(int),n_items,file);
	free(iq);
	}
else if(strncmp(tipout,"fl",2)==0)  fwrite(q,sizeof(float),n_items,file);
else error_io(101,tipout);

}

/*====================================== E S C R I T U R A - N F =============*/

void output_ima_nf(int n_fil,int n_col,float *q,char *nom_fitxer,char tipout[3])
{

float *qN;
FILE *file;
int i,n_items;
unsigned short int *cq;
unsigned char *uq;

if((file=fopen(nom_fitxer,"w"))==NULL) error_io(106,nom_fitxer);

n_items=n_fil*n_col;

if(strncmp(tipout,"1b",2)==0){
        uq=(unsigned char *)malloc(n_items);
        qN=(float *)malloc(n_items*4);
        normalitzacio_imatge(q,qN,n_items);
        for(i=0;i<n_items;i++) *(uq+i)=(unsigned char) *(qN+i);
        fwrite(uq,1,n_items,file);
        free(qN);
        free(uq);
        }
else if(strncmp(tipout,"2b",2)==0){
        cq=(unsigned short int *)malloc(n_items*2);
        for(i=0;i<n_items;i++) *(cq+i)=(unsigned short int) *(q+i);
        fwrite(cq,2,n_items,file);
        free(cq);
        }
else if(strncmp(tipout,"fl",2)==0) fwrite(q,sizeof(float),n_items,file);
else error_io(101,tipout); 

fclose(file);            
}

/*=================================== E S C R I T U R A - S E Q- N F ==========*/

void output_seq_nf(int n_fil,int n_col,int nima,float *q,char *nom_fitxer,char tipout[3])
{
float *qN;
FILE *file;
register int i,n_items,k,kMM,t_ima;
unsigned short int *cq;
unsigned char *uq;

if((file=fopen(nom_fitxer,"w"))==NULL) error_io(106,nom_fitxer);

n_items=n_fil*n_col*nima;

if(strncmp(tipout,"1b",2)==0){
        t_ima=n_fil*n_col;
        uq=(unsigned char *)malloc(t_ima);
        qN=(float *)malloc(t_ima*4);
        for(k=0,kMM=0;k<nima;k++,kMM+=t_ima){
           normalitzacio_imatge(q+kMM,qN,t_ima);
           for(i=0;i<t_ima;i++) *(uq+i)=(unsigned char) *(qN+i);
           fwrite(uq,1,t_ima,file);
           }
        free(qN);
        free(uq);
        }
else if(strncmp(tipout,"2b",2)==0){
        cq=(unsigned short int *)malloc(n_items*sizeof(short int));
        for(i=0;i<n_items;i++) *(cq+i)=(unsigned short int) *(q+i);
        fwrite(cq,2,n_items,file);
        free(cq);
        }
else if(strncmp(tipout,"fl",2)==0) fwrite(q,sizeof(float),n_items,file);
else error_io(101,tipout); 

fclose(file);            
}

/*========================== NORMALITZACIO DE LA IMATGE A UN VALOR MIG========*/

void normalitzacio_imatge_vm(float *q,float *qout,int nitems,float vm)
{   
register int i;
float factor,nc;

numeroc(q,&nc,nitems);
factor=vm/nc;
for (i=0;i<nitems;i++){
  *(qout+i)= *(q+i)*factor;
  if(*(qout+i)>255.) *(qout+i)=255.;
  }
}

/*=================================== NORMALITZACIO DE LA IMATGE AL MAXIM ====*/
/* modificada para reescalar tambien los negativos julio06 */
void normalitzacio_imatge(float *q,float *qout,int nitems)
{   
register int i;
float mx,mn,factor;

/*reescala por minimo negativo */
mn=10;
for(i=0;i<nitems;i++) mn=mini(mn, *(q+i) );
for (i=0;i<nitems;i++) *(q+i)= *(q+i)-mn;

/* reescala por maximo */
mx=-10;
for(i=0;i<nitems;i++) mx=maxi(mx, *(q+i) );
factor=255./mx;
for (i=0;i<nitems;i++) *(qout+i)= *(q+i)*factor;

}

/*======================================= NUMERO DE CONTES ===================*/

void numeroc(float *l,float *res,int nitems)
{

register int i;

 *res=0;
 for(i=0;i<nitems;i++) *res+=*(l+i);

}

/*======================================= NUMERO DE CONTES PRINT =============*/

void numerocp(float *l,float *res,int nitems)
{

register int i;

 *res=0;
 for(i=0;i<nitems;i++) {
	 *res+=*(l+i);
 }
printf("%f\n",*res);

}

/*======================================= DISPERSIO ===================*/

void dispersio(float *l,int vm,float *res,int nitems)
{

register int i;

 *res=0;
 for(i=0;i<nitems;i++) *res+=(*(l+i)-vm)*(*(l+i)-vm);
 *res=sqrt(*res/(nitems-1));
}

/*============================== OMPLE =======================================*/

void omple(float *ima,int nitems,float valor)
{
register int i;

 for(i=0;i<nitems;i++) *(ima+i)=valor;

}

/*============================== OMPLE CILINDRE ==========================*/

void omple_cilindre(float *ima,int nfil,int ntalls,float valor)
{
register int i,j,k,ind,nfild2;
float x,y,R2;

nfild2=nfil/2;
R2=(nfild2)*(nfild2);
for(k=0,ind=0;k<ntalls;k++){
  for(i=0,x=(float)(-nfild2)+0.5;i<nfil;i++,x+=1.){
    for(j=0,y=(float)(-nfild2)+0.5;j<nfil;j++,ind++,y+=1.){
      if((x*x+y*y)<R2) ima[ind]=valor;
      else ima[ind]=0.;
      }
    }
  }
}

/*============================== OMPLE CILINDRE ==========================*/

void omple_cilindre_mes_petit(float *ima,int nfil,int ntalls,float valor, int rest)
{
register int i,j,k,ind,nfild2;
float x,y,R2;

nfild2=nfil/2;
R2=(nfild2-rest)*(nfild2-rest);
for(k=0,ind=0;k<ntalls;k++){
  for(i=0,x=(float)(-nfild2)+0.5;i<nfil;i++,x+=1.){
    for(j=0,y=(float)(-nfild2)+0.5;j<nfil;j++,ind++,y+=1.){
      if((x*x+y*y)<R2) ima[ind]=valor;
      else ima[ind]=0.;
      }
    }
  }
}
/*============================== OMPLE MASCARA : ccrespo jul 2006 ==========================*/
void omple_mascara_rec3d(float *ima,int nfil,int ntalls,float valor)
{
register int i,nitems;

nitems=nfil*nfil*ntalls;
for(i=0;i<nitems;i++){
if(ima[i] != 0.) {
        ima[i]=valor;
}
}
}
/*==============================================  ITOA  ======================*/

char *itoa(int n,char *s)
{
int i,sign,low,hi;
char c;

if((sign=n)<0) n=-n;
i=0;
do{
  s[i++]=n%10+'0';
  }while((n/=10)>0);
if(sign<0) s[i++]='-';
s[i]=EOS;

for(low=0,hi=i-1;low<hi;low++,hi--){
    c=s[low];
    s[low]=s[hi];
    s[hi]=c;
    }
return(s);
}

/*======================================== Lee imagen en formato Analyze, paguiar ========================== */   
void lee_imagen_hdr(const char *nom_fitxer_hdr, struct imagen *ima)
{
char nom_fitxer_img[128];
short int bitpix,datatype;
int nitems;
/*int i,j1,j2,j3,j4,u1,u2,u3,u4; Variables para Big/Little */
struct dsr hdr; /* Estructura de datos generada por load_hdr(),- dbh.h */

/* Se genera el nombre para *.img */
memset(nom_fitxer_img,0,128);
strncpy(nom_fitxer_img,nom_fitxer_hdr,strlen(nom_fitxer_hdr)-3);
strcat(nom_fitxer_img,"img");

/* Funci�n de libreria dbh.h que lee *.hdr y genera la struct hdr */
//load_hdr(nom_fitxer_hdr,&hdr);  // OJO: Ver que pasa cuando no existe el *.hdr, falta funcion de error!!!??
if(load_hdr(nom_fitxer_hdr,&hdr)!=0) error_io(102,(char *)nom_fitxer_hdr);

/* OJO: dim[2] son las filas en el mricro y dim[1] son las columnas corregido Junio 05*/
ima->nfil=hdr.dime.dim[2];
ima->ncol=hdr.dime.dim[1];
ima->nima=hdr.dime.dim[3];
ima->tpix=hdr.dime.pixdim[1];
ima->tpiy=hdr.dime.pixdim[2];
ima->tcorte=hdr.dime.pixdim[3];
ima->offset=hdr.dime.vox_offset;
bitpix=hdr.dime.bitpix;
datatype=hdr.dime.datatype;
nitems=ima->nfil*ima->ncol*ima->nima;
ima->datos=vector(0,nitems-1);

/* Se convierten los valores de dbh.h al formato de input_fitxer() Falta ver si funcionan signed y unsigned */

	if(bitpix==8 && datatype==2) {
			strcpy(ima->tipo,"1b");/*lectura como char*/
			input_fitxer_offset(ima->nfil,ima->ncol,ima->nima,ima->datos,nom_fitxer_img,ima->tipo,ima->offset);	
			}
	else if(bitpix==16 && datatype==3){ /* Extension local para nosotros de Analyze 7.5 Julio2006 */
			strcpy(ima->tipo,"2b"); /* lectura como unsigned short */
			input_fitxer_offset(ima->nfil,ima->ncol,ima->nima,ima->datos,nom_fitxer_img,ima->tipo,ima->offset);
			}
	else if(bitpix==16 && datatype==4){ 
			strcpy(ima->tipo,"i2"); /* lectura como signed short */
			input_fitxer_offset(ima->nfil,ima->ncol,ima->nima,ima->datos,nom_fitxer_img,ima->tipo,ima->offset);
			}
	else if(bitpix==16 && datatype!=4 && datatype!=3){ 
			strcpy(ima->tipo,"i2"); /* lectura como signed short */
			//printf("Warning: Error in datatype (default is i2)\n");
			input_fitxer_offset(ima->nfil,ima->ncol,ima->nima,ima->datos,nom_fitxer_img,ima->tipo,ima->offset);
			}
	else if(bitpix==32 && datatype==8){ 
			strcpy(ima->tipo,"i4"); /* lectura como signed int */
			input_fitxer_offset(ima->nfil,ima->ncol,ima->nima,ima->datos,nom_fitxer_img,ima->tipo,ima->offset);
			}						
	else if(bitpix==32 && datatype==16) {
			strcpy(ima->tipo,"fl");
			input_fitxer_offset(ima->nfil,ima->ncol,ima->nima,ima->datos,nom_fitxer_img,ima->tipo,ima->offset);
			}
	else if(bitpix==64 && datatype==64) {
			strcpy(ima->tipo,"db");
			input_fitxer_offset(ima->nfil,ima->ncol,ima->nima,ima->datos,nom_fitxer_img,ima->tipo,ima->offset);
			}
}

/*======================================== Lee imagen en formato Analyze, paguiar ========================== */   
void lee_template_hdr(const char *nom_fitxer_hdr, struct imagen *ima)
{
short int bitpix,datatype;
int nitems;
/*int i,j1,j2,j3,j4,u1,u2,u3,u4; Variables para Big/Little */
struct dsr hdr; /* Estructura de datos generada por load_hdr(),- dbh.h */

/* Funci�n de libreria dbh.h que lee *.hdr y genera la struct hdr */
//load_hdr(nom_fitxer_hdr,&hdr);  // OJO: Ver que pasa cuando no existe el *.hdr, falta funcion de error!!!??
if(load_hdr(nom_fitxer_hdr,&hdr)!=0) error_io(102,(char *)nom_fitxer_hdr);

/* OJO: dim[2] son las filas en el mricro y dim[1] son las columnas corregido Junio 05*/
ima->nfil=hdr.dime.dim[2];
ima->ncol=hdr.dime.dim[1];
ima->nima=hdr.dime.dim[3];
ima->tpix=hdr.dime.pixdim[1];
ima->tcorte=hdr.dime.pixdim[3];
ima->offset=hdr.dime.vox_offset;
bitpix=hdr.dime.bitpix;
datatype=hdr.dime.datatype;
nitems=ima->nfil*ima->ncol*ima->nima;
ima->datos=vector(0,nitems-1);

bitpix=hdr.dime.bitpix;

	if(bitpix==8 && datatype==2) {
		strcpy(ima->tipo,"1b");
		}
	else if(bitpix==16 && datatype==3){
		strcpy(ima->tipo,"2b");
		}
	else if(bitpix==16 && datatype==4){
		strcpy(ima->tipo,"i2");
		}
	else if(bitpix==32 && datatype==8){
		strcpy(ima->tipo,"i4");
		}
	else if(bitpix==32 && datatype==16) {
		strcpy(ima->tipo,"fl");
		}
	else if(bitpix==64 && datatype==64) {
		strcpy(ima->tipo,"db");
		}

}



/*======================================== Lee imagen en formato Analyze, paguiar ========================== */   
void lee_imagen_hv(const char *nom_fitxer_hv, struct imagen *ima)
{
	char **palabras;
	char nom_fitxer_v[180];
	int *ntokens,i,nitems;
	
	/* Se genera el nombre para *.v */
	memset(nom_fitxer_v,0,128);
	strncpy(nom_fitxer_v,nom_fitxer_hv,strlen(nom_fitxer_hv)-2);
	strcat(nom_fitxer_v,"v");
	
	ntokens=malloc(sizeof (int));
	cuenta_tokens_ascii(nom_fitxer_hv,ntokens);
	palabras=cmatrix(0,*ntokens,0,100); /* cada token tiene maximo 100 caracteres */
	separa_tokens_ascii(nom_fitxer_hv,palabras);
	
	for(i=0;i<*ntokens;i++)
	{
	
		if(strcmp(palabras[i],"size")==0 && strcmp(palabras[i+1],"[1]")==0 && strcmp(palabras[i+2],":=")==0) 
			ima->nfil=atoi(palabras[i+3]);
		if(strcmp(palabras[i],"size")==0 && strcmp(palabras[i+1],"[2]")==0 && strcmp(palabras[i+2],":=")==0) 
			ima->ncol=atoi(palabras[i+3]);
		if(strcmp(palabras[i],"size")==0 && strcmp(palabras[i+1],"[3]")==0 && strcmp(palabras[i+2],":=")==0) 
			ima->nima=atoi(palabras[i+3]);
		if(strcmp(palabras[i],"(mm")==0 && strcmp(palabras[i+2],"[1]")==0 && strcmp(palabras[i+3],":=")==0) 
			ima->tpix=atof(palabras[i+4]);
		if(strcmp(palabras[i],"(mm")==0 && strcmp(palabras[i+2],"[3]")==0 && strcmp(palabras[i+3],":=")==0) 
			ima->tcorte=atof(palabras[i+4]);
		if(strcmp(palabras[i],"data")==0 && strcmp(palabras[i+1],"offset")==0 && strcmp(palabras[i+2],"in")==0) 
			ima->offset=atof(palabras[i+5]);
		if(strcmp(palabras[i],"!number")==0 && strcmp(palabras[i+1],"format")==0 && strcmp(palabras[i+3],"float")==0) 
			strncpy(ima->tipo,"fl\0",3);
		if(strcmp(palabras[i],"bytes")==0 && strcmp(palabras[i+2],"pixel")==0 && strcmp(palabras[i+4],"4")==0) 
			strncpy(ima->tipo,"fl",2);
	}
	
	free_cmatrix(palabras,0,*ntokens,0,100);
	
	nitems=ima->nfil*ima->ncol*ima->nima;
	ima->datos=vector(0,nitems-1);
	
	if(strncmp(ima->tipo,"fl",2)==0) 
	{
		input_fitxer_offset(ima->nfil,ima->ncol,ima->nima,ima->datos,nom_fitxer_v,ima->tipo,ima->offset);
	}
	else
	{
		printf("Error! Tipo de dato de header interfile no soportado\n\n");
		exit(1);
	}
}




/*======================================== Lee imagen en formato cab del CT del GIR-USC ========================== */   
void lee_imagen_cab(const char *nom_fitxer_cab, struct imagen *ima)
{
double *data;
FILE *file;
int n,i,nitems;
	
if((file = fopen(nom_fitxer_cab,"r"))==NULL) error_io(102,(char*)nom_fitxer_cab);	

data=malloc(sizeof(int));
fread (data,sizeof(double),1,file);
ima->ncol=(int)*data;

fread (data,sizeof(double),1,file);
ima->nfil=(int)*data;

fread (data,sizeof(double),1,file);
ima->nima=(int)*data;

fread (data,sizeof(double),1,file);
ima->tpix=(float)*data;

fread (data,sizeof(double),1,file);
ima->tcorte=(float)*data;
if(ima->tcorte==0) ima->tcorte=1.;

fread (data,sizeof(double),1,file);
ima->Dfo=(float)*data;

fread (data,sizeof(double),1,file);
ima->Dfd=(float)*data;

ima->offset=0.;
strncpy(ima->tipo,"fl\0",3);

// To read data using values (hard code)
ima->ncol=498;
ima->nfil=725;
ima->nima=202;
ima->tpix=0.389;
ima->tcorte=(float)*data;
ima->tcorte=0.5;
ima->Dfo=390;
ima->Dfd=1180;
ima->offset=0.;
strncpy(ima->tipo,"fl\0",3);

// To read data using header values
nitems=ima->nfil*ima->ncol*ima->nima;
ima->datos=vector(0,nitems-1);
data=dvector(0,(ima->nfil*ima->ncol)-1);
for(n=0;n<ima->nima;n++)
{
    fread(data,8,ima->nfil*ima->ncol,file);
    for(i=0;i<ima->nfil*ima->ncol;i++) ima->datos[i]=(float)data[i];
}

}








/*======================================== Guarda imagen en formato Analyze, paguiar ========================== */   
/* Falta escribir defaults */
void guarda_imagen_hdr(const char *nom_fitxer_hdr,struct imagen *ima)
{
/*int i,j1,j2,j3,j4,u1,u2,u3,u4; Variables Big/Endian */
int nitems;
char nom_fitxer_img[128];
struct dsr hdr; /* Estructura de datos generada por load_hdr(),- dbh.h */
FILE *file;

memset(nom_fitxer_img,0,128);
/* Se Genera el nombre para *.img */
strncpy(nom_fitxer_img,nom_fitxer_hdr,strlen(nom_fitxer_hdr)-3);
strcat(nom_fitxer_img,"img");
if((file=fopen(nom_fitxer_img,"wb"))==NULL) error_io(106,nom_fitxer_img);
memset(&hdr,0,sizeof(hdr));

/* Se genera la estructura hdr */
hdr.hk.sizeof_hdr=sizeof(hdr);
hdr.hk.extents=0;
hdr.hk.regular='r';

hdr.dime.dim[0]=4;

/* OJO: dim[2] son las filas en el mricro y dim[1] son las columnas corregido Junio 05*/
hdr.dime.dim[2]=ima->nfil;
hdr.dime.dim[1]=ima->ncol;
hdr.dime.dim[3]=ima->nima;
nitems=ima->nfil*ima->ncol*ima->nima;
hdr.dime.dim[4]=1;
hdr.hist.orient=0;
strcpy(hdr.dime.vox_units,"");
strcpy(hdr.dime.cal_units,"");
hdr.dime.cal_max=0.0;
hdr.dime.cal_min=0.0;

hdr.dime.pixdim[5]=0.0;
hdr.dime.pixdim[6]=0.0;
hdr.dime.pixdim[7]=0.0;
hdr.dime.pixdim[8]=0.0;
hdr.dime.pixdim[0]=1.0;
hdr.dime.pixdim[1]=ima->tpix;
hdr.dime.pixdim[2]=ima->tpiy;
hdr.dime.pixdim[3]=ima->tcorte;
hdr.dime.pixdim[4]=0.0;
//hdr.dime.vox_offset=ima->offset;
hdr.dime.vox_offset=0; /* en ppio es asi pq no se guarda nunca con offset ene-11 */
hdr.dime.funused1=1.0; /* scale factor */
hdr.dime.funused2=0.0;
hdr.dime.funused3=0.0;
hdr.dime.cal_min=0.0;
hdr.dime.cal_max=0.0;


if(strncmp(ima->tipo,"1b",2)==0) { /* 1-Byte normalizando a escala 255 */
		hdr.dime.bitpix=8;
		hdr.dime.datatype=DT_UNSIGNED_CHAR;
		hdr.dime.glmin=0;
		hdr.dime.glmax=255;
		printf("\nOJO: Se normaliza maximo a 255!\n" );
		output_ima3(ima->nfil,ima->ncol,ima->nima,ima->datos,file,"1b");
		}
else if(strncmp(ima->tipo,"1B",2)==0) { /* 1-Byte sin normalizar a escala 255 */
		hdr.dime.bitpix=8;
		hdr.dime.datatype=DT_UNSIGNED_CHAR;
		hdr.dime.glmin=0;
		hdr.dime.glmax=255;
		output_ima3(ima->nfil,ima->ncol,ima->nima,ima->datos,file,"1B");
		}
		else if(strncmp(ima->tipo,"2b",2)==0){ /* tipo definido como extension a Analyze 7.5 en dbh.h Julio 2006 */
		hdr.dime.bitpix=16;
		hdr.dime.datatype=DT_UNSIGNED_SHORT;
		hdr.dime.glmin=0;
		hdr.dime.glmax=65536;
		output_ima3(ima->nfil,ima->ncol,ima->nima,ima->datos,file,"2b");
			}
else if(strncmp(ima->tipo,"i2",2)==0){  
		hdr.dime.bitpix=16;
		hdr.dime.datatype=DT_SIGNED_SHORT;
		hdr.dime.glmin=-32768;
		hdr.dime.glmax=32767;
		output_ima3(ima->nfil,ima->ncol,ima->nima,ima->datos,file,"i2");
			}
else if(strncmp(ima->tipo,"i4",2)==0){  
		hdr.dime.bitpix=32;
		hdr.dime.datatype=DT_SIGNED_INT;
/*		hdr.dime.glmin=-2147483647;
		hdr.dime.glmax=2147483648;*/
		output_ima3(ima->nfil,ima->ncol,ima->nima,ima->datos,file,"i4");
			}			
else if(strncmp(ima->tipo,"fl",2)==0) {
		hdr.dime.bitpix=32;
		hdr.dime.datatype=DT_FLOAT;
		output_ima3(ima->nfil,ima->ncol,ima->nima,ima->datos,file,"fl");
		}

/* Funci�n de libreria dbh.h que lee *.hdr y genera la struct hdr */
save_hdr(nom_fitxer_hdr,&hdr);
fclose(file);
}

void print_hdr(const char *nom_fitxer_hdr)
{
struct dsr hdr; /* Estructura de datos generada por load_hdr(),- dbh.h */
char tipo[2];

/* Funci�n de libreria dbh.h que lee *.hdr y genera la struct hdr */
//load_hdr(nom_fitxer_hdr,&hdr);  // OJO: Ver que pasa cuando no existe el *.hdr, falta funcion de error!!!??
if(load_hdr(nom_fitxer_hdr,&hdr)!=0) error_io(102,(char *)nom_fitxer_hdr);


/* Se convierten los valores de dbh.h al formato de input_fitxer() Falta ver si funcionan signed y unsigned */
	if(hdr.dime.bitpix==8) {
			strcpy(tipo,"1b");	
			}
	else if(hdr.dime.bitpix==16){
			strcpy(tipo,"2b");
			}
	else if(hdr.dime.bitpix==32) {
			strcpy(tipo,"fl");
			}
			
printf("Parametros contenidos en el header: \n");	
printf("Nombre del archivo:%s\n",nom_fitxer_hdr);
printf("\n");
printf("Pixels en direccion Z (nslices): %d\n",hdr.dime.dim[3]);	
printf("Pixels en direccion X (xbins): %d\n",hdr.dime.dim[1]);
printf("Pixels en direccion Y (ybins): %d\n",hdr.dime.dim[2]);
printf("\n");
printf("Tamano pixel en Z (zpix): %f\n",hdr.dime.pixdim[3]);	
printf("Tamano pixel en X (xpix): %f\n",hdr.dime.pixdim[1]);
printf("Tamano pixel en Y (ypix): %f\n",hdr.dime.pixdim[1]);
printf("\n");
printf("Tipo de datos es: %s\n",tipo);
printf("\n");
printf("\n");
}





/* guarda imagen interfile (HS de sinograma-proyeccion) a partir una estructura imagen*/
void guarda_imagen_hs(char *nombre_archivo_hs,struct imagen *ima)
{
char nombre_archivo_s[128]; /* nombre archivo para guardar la imagen */
int i; /* numero bytes por pixel en if */
FILE *file,*file2;

memset(nombre_archivo_s,0,128);
strncpy(nombre_archivo_s,nombre_archivo_hs,strlen(nombre_archivo_hs)-2);
strcat(nombre_archivo_s,"s");
if((file=fopen(nombre_archivo_hs,"wb"))==NULL) error_io(106,nombre_archivo_hs);
if((file2=fopen(nombre_archivo_s,"wb"))==NULL) error_io(106,nombre_archivo_s);



fprintf(file,"!INTERFILE  :=\n");
fprintf(file,"name of data file := %s\n",nombre_archivo_s);
fprintf(file,"!originating system := Userdefined\n");
fprintf(file,"!GENERAL DATA :=\n");
fprintf(file,"!GENERAL IMAGE DATA :=\n");
fprintf(file,"!type of data := PET\n");
fprintf(file,"imagedata byte order := LITTLEENDIAN\n");
fprintf(file,"!PET STUDY (General) :=\n");
fprintf(file,"!PET data type := Emission\n");
fprintf(file,"applied corrections := {arc correction}\n");
if(strncmp(ima->tipo,"1b",2)==0)
{
printf("ERROR: Todavia no implementada funcion para 1b\n");
exit(1);
}
if(strncmp(ima->tipo,"2b",2)==0)
{
printf("ERROR: Todavia no implementada funcion para 2b\n");
exit(1);
}
if(strncmp(ima->tipo,"fl",2)==0)
{
fprintf(file,"!number format := float\n");
fprintf(file,"!number of bytes per pixel := 4\n");
fprintf(file,"\n");output_ima3(ima->nfil,ima->ncol,ima->nima,ima->datos,file2,"fl");
}

fprintf(file,"number of dimensions := 4\n");
fprintf(file,"matrix axis label [4] := segment\n");
fprintf(file,"!matrix size [4] := %d\n",ima->nsegments);
fprintf(file,"matrix axis label [2] := view\n");
fprintf(file,"!matrix size [2] := %d\n",ima->nfil);
fprintf(file,"matrix axis label [3] := axial coordinate\n");
fprintf(file,"!matrix size [3] := {");
i=0;
while(ima->axial_coord[i]!=0)
{
fprintf(file,"%d ",ima->axial_coord[i]);
i++;
}
fprintf(file,"}\n");
fprintf(file,"matrix axis label [1] := tangential coordinate\n");
fprintf(file,"!matrix size [1] := %d\n",ima->ncol);
fprintf(file,"minimum ring difference per segment := { ");
i=0;
while(ima->min_rd[i]!=0)
{
printf("min_rd:%d\n",ima->min_rd[i]);
fprintf(file,"%d ",ima->min_rd[i]);
i++;
}
fprintf(file,"}\n");
fprintf(file,"maximum ring difference per segment := { ");
i=0;
while(ima->max_rd[i]!=0)
{
fprintf(file,"%d ",ima->max_rd[i]);
i++;
}
fprintf(file,"}\n");
fprintf(file,"Scanner parameters:=\n");
fprintf(file,"Scanner type := Userdefined\n");
fprintf(file,"number of rings := %d\n",ima->nrings);
fprintf(file,"number of detectors per ring := 576\n");
fprintf(file,"inner ring diameter (cm) := %f\n",ima->dring-0.5);
fprintf(file,"average depth of interaction (cm) := 0.25\n");
fprintf(file,"distance between rings (cm) := %f\n",ima->tcorte/10.);
fprintf(file,"default bin size (cm) := %f\n",ima->tpix/10.);/*
fprintf(file,"Number of blocks per bucket in transaxial direction         := 3\n");
fprintf(file,"Number of blocks per bucket in axial direction              := 4\n");
fprintf(file,"Number of crystals per block in axial direction             := 6\n");
fprintf(file,"Number of crystals per block in transaxial direction        := 8\n");*/
fprintf(file,"end scanner parameters:=\n");
fprintf(file,"effective central bin size (cm) := %f\n",ima->tpix/10.);
fprintf(file,"data offset in bytes[1] := 0\n");
fprintf(file,"number of time frames := 1\n");
fprintf(file,"!END OF INTERFILE :=\n");
}





/* guarda imagen interfile (HV de una imagen) a partir una estructura imagen*/
void guarda_imagen_hv(char *nombre_archivo_hv,struct imagen *ima)
{
char nombre_archivo_v[128]; /* nombre archivo para guardar la imagen */
FILE *file,*file2;

memset(nombre_archivo_v,0,128);
strncpy(nombre_archivo_v,nombre_archivo_hv,strlen(nombre_archivo_hv)-2);
strcat(nombre_archivo_v,"v");
if((file=fopen(nombre_archivo_hv,"wb"))==NULL) error_io(106,nombre_archivo_hv);
if((file2=fopen(nombre_archivo_v,"wb"))==NULL) error_io(106,nombre_archivo_v);

fprintf(file,"!INTERFILE  :=\n");
fprintf(file,"name of data file := %s\n",nombre_archivo_v);
fprintf(file,"!GENERAL DATA :=\n");
fprintf(file,"!GENERAL IMAGE DATA :=\n");
fprintf(file,"!type of data := PET\n");
fprintf(file,"imagedata byte order := LITTLEENDIAN\n");
fprintf(file,"!PET STUDY (General) :=\n");
fprintf(file,"!PET data type := Image\n");

if(strncmp(ima->tipo,"1b",2)==0)
{
printf("ERROR: Todavia no implementada funcion para 1b\n");
exit(1);
}
if(strncmp(ima->tipo,"2b",2)==0)
{
printf("ERROR: Todavia no implementada funcion para 2b\n");
exit(1);
}
if(strncmp(ima->tipo,"fl",2)==0)
{
fprintf(file,"!number format := float\n");
fprintf(file,"!number of bytes per pixel := 4\n");
output_ima3(ima->nfil,ima->ncol,ima->nima,ima->datos,file2,"fl");
}

/* dimensions */
fprintf(file,"number of dimensions := 3\n");
fprintf(file,"matrix axis label [1] := x\n");
fprintf(file,"!matrix size [1] := %d\n",ima->nfil);
fprintf(file,"scaling factor (mm/pixel) [1] := %f\n",ima->tpix);
fprintf(file,"matrix axis label [2] := y\n");
fprintf(file,"!matrix size [2] := %d\n",ima->ncol);
fprintf(file,"scaling factor (mm/pixel) [2] := %f\n",ima->tpix);
fprintf(file,"matrix axis label [3] := z\n");
fprintf(file,"!matrix size [3] := %d\n",ima->nima);
fprintf(file,"scaling factor (mm/pixel) [3] := %f\n",ima->tcorte);

fprintf(file,"data offset in bytes[1] := 0\n");
fprintf(file,"number of time frames := 1\n");
fprintf(file,"!END OF INTERFILE :=\n");
}


/*=========================================== ERROR =========================*/

void error_io(int nerr,char *text)
{

switch(nerr){
    case 32:  printf("\nFalten dades al fitxer de dades d'un fitxer: %s\n",text);
              printf("\nLes dades que hi ha d'haver son:\n\tNom del fitxer\n\tNumero de files\n\tNumero de columnes\n\ttipus de fitxer (1b-2b-fl)\n");break;
    case 33:  printf("\nFalten dades al fitxer de dades d'un fitxer: %s\n",text);
              printf("\nLes dades que hi ha d'haver son:\n\tNom del fitxer\n\tNumero de files\n\tNumero de columnes\n\ttipus de fitxer (1b-2b-fl)\n\tNumero d'imatges\n");break;
    case 41:  printf("\n\nERROR al seleccionar sequencia:\n\t(primera+salt*nima) s'excedeix del nombre total d'imatges disponibles\n"); break;
    case 101: printf("\n\nERROR en el format d'escritura. Format: %s  erroni\n",text); break;
    case 102: printf("\n\n%s file not found\n",text); break;
    case 103: printf("\n\nError en el format del fitxer de lectura. Format: %s  erroni\n",text); break;
    case 104: printf("\n\nError en el format del fitxer de la sequencia. Format: %s  erroni\n",text);break;
    case 105: printf("\n\nError al obrir el fitxer de la sequencia. El fitxer %s no existeix\n",text); break;
    case 106: printf("\n\nError al abrir el fichero de escritura, el fichero %s no existe\n",text); break;
    default: printf("\n\nError del numero d'error en la funcio error_io()\n"); 
    }
    
exit(0);
} 

