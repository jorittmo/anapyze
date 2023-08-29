#include <stdlib.h>
#include <string.h>
#include <pet/indexOS_pet.h>
#include <util/inpout.h>

#define maxim(a,b) ((a)>=(b)?(a):(b))
#define minim(a,b) ((a)<=(b)?(a):(b))
#define abs(a) ((a)>=0?(a):(-a))

/*================================== INDEXS ==================================*/

void calcul_indexs(int *indexs,int OS,int NangdOS)
{
int i,j,k,im,m,n,*ple,*iOS,*a,*sa,*dif;

iOS=(int *)calloc(OS,sizeof(int));
ple=(int *)calloc(OS,sizeof(int));
a=(int *)calloc(OS,sizeof(int));
sa=(int *)calloc(OS,sizeof(int));
dif=(int *)calloc(OS+1,sizeof(int));

for(i=0;i<OS;i++) iOS[i]=ple[i]=0;
ple[0]=1;

calcul_dif(dif,OS);

im=0;
for(k=1;k<OS;k++){
  for(i=1;i<OS;i++){
    if(!ple[i]){
      j=i-im;
      a[i]=dif[abs(j)];
      }
    }
  for(i=0;i<OS;i++){
      a[i]*=(1-ple[i]);
      sa[i]*=(1-ple[i]);
      sa[i]+=a[i];
      }
  m=n=0;
  for(i=0;i<OS;i++) m=maxim(m,sa[i]);
  for(i=1;i<OS;i++) if(!ple[i]) m=minim(m,sa[i]);
  for(i=1;i<OS;i++) if(sa[i]==m) n=maxim(n,a[i]);
  for(i=OS-1;i>0;i--) if(sa[i]==m) if(a[i]<=n){
       n=a[i];
       im=i;
       }
  iOS[k]=im;
  ple[im]=1;
  }

for(i=0;i<OS;i++) for(j=0;j<NangdOS;j++) indexs[j+i*NangdOS]=iOS[i]+OS*j;

free(a);free(sa);free(dif);free(ple);
}

/*.........................................................................*/

void calcul_dif(int *dif,int OS)
{
int i,OS2;
OS2=OS*OS;

for(i=1;i<=OS;i++) dif[i]+=2*(i*(i-OS))+OS2;
}

/* ============================ ORDENAR PROJECCIONS PMPxOS ======================
    Aquesta rutina agafa unes projeccions en l'ordre dels angles dels subsets i la
    torna ordenada en angles correlatius.................					   */

void ordenar_proj_pmp3dxOS(float *le,float *le2,int Nang,int Nbt,int Nbp,int OS)
{
int i,index_1,index_2,k,*ordre_proj,NangOS;

NangOS=Nang/OS;
ordre_proj=(int*)calloc(Nang,sizeof(int));  

/* .... copia de projeccio a vector auxiliar ... */
for(i=0;i<Nbt;i++) le2[i]=le[i];

/* .... calcul de l'ordre dels angles de projeccio als subsets ... */          
calcul_indexs(ordre_proj,OS,NangOS);  

/* .... per cada angle, ordenar segons angles correlatius ... */
for(k=0,index_1=0 ; k<Nang ;k++)
   for(i=0,index_2=ordre_proj[k]*Nbp; i<Nbp ; i++,index_1++,index_2++) le[index_2]=le2[index_1];
      
free(ordre_proj);
}

/* ============================ REORDENAR PROJECCIONS ============================
    Aquesta rutina agafa unes projeccions ordenada en angles correlatius i la
    torna ordenada en l'ordre d'angles determinats pels subsets.................	   */
    // nov-04 modificado por paguiar para PET
    
    
void reordenar_projxOS_pet(float *le,float *le2,int Nang,int Nbins,int Nsino,int OS)
{
int i,k,kr,kw,l,NangOS,NangOStot,OSxNbins;
NangOS=Nang/OS;  //número de angulos de cada subset
kw=0;
NangOStot=NangOS*Nsino; //número de angulos totales en cada subset
OSxNbins=OS*Nbins; //paso para leer solo angulos de un mismo subset

// Recorre todos los subsets
for(k=0;k<OS;k++)
{
      // kr recorre sinogramas leyendo filas (angulos) de forma discontinua,
      // es decir, saltando tantas filas como subsets se hayan seleccionado.
     for(l=0,kr=k*Nbins;l<NangOStot;l++,kr+=OSxNbins)
     {   
     	for (i=0;i<Nbins;i++,kw++)
		{
     		le2[kw]=le[kr+i];
     	}
     }
}
}

/* ========================= POSAR INDEX AL NOM DE LA MATRIU DE PESOS =================*/

void generar_nom_mpes(char *mpes_nou,char *mpes_original,int k,int OS)
{
char p[3];	/* variable auxiliar per a convertir l'ordre de subset en string*/

strcpy(mpes_nou,mpes_original);
if(OS>1){
     strcat(mpes_nou,".OS");
     strcat(mpes_nou,itoa(k,p));
     } 
}
