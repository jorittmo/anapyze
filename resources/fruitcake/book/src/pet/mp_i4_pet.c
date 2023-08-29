#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pet/mp_i4_pet.h>
#include <pet/STIR.h>


/*============================= FWD_PROJ ====================================*/

void fwd_proj(float *lc,float *q,sparse_i4 *mp)
{
int i,j,NbOS;
int j_min,j_max;

NbOS=mp->Nbins*mp->Nang*num_imagenes(mp->ndet_z,mp->axial_diff)/mp->Nsubsets;
 
for(i=0;i<NbOS;i++)
{
    lc[i]=0.;

    // Explicación: ¿Como se calcula un punto de los sinogramas? 
	// Primero se lee el valor correspondiente de mp.ia[i] y mp.ia[i+1]
	// para saber cuantos valores no nulos hay en la imagen. Se encuentra el intervalo desde j_min hasta j_max 
	// y por lo tanto se obtienen los valores correspondientes en mp.ar, que a su vez se multiplican por el valor de la 
	// imagen en ese mismo punto.
	
    j_min=mp->ia[i];
    j_max=mp->ia[i+1];
	for(j=j_min;j<j_max;j++)  
    {
		lc[i]+= (q[mp->ja[j]])*(mp->ar[j]);
		//printf("mu:%f  X  dx:%f\n",q[mp->ja[j]],mp->ar[j]);
    }
    //printf("valor:%f  ",lc[i]);
    lc[i]=exp(-lc[i]);
    //printf("exp(-valor):%f\n\n\n",lc[i]);	
}

}




/*============================= FWD_PROJ ====================================*/

void fwd_proj_attenuation(float *lc,float *q,sparse_i4 *mp)
{
int i,j,NbOS;
int j_min,j_max;

NbOS=mp->Nbins*mp->Nang*num_imagenes(mp->ndet_z,mp->axial_diff)/mp->Nsubsets;
 
for(i=0;i<NbOS;i++)
{
    lc[i]=0.;
    
    // Explicación: ¿Como se calcula un punto de los sinogramas? 
	// Primero se lee el valor correspondiente de mp.ia[i] y mp.ia[i+1]
	// para saber cuantos valores no nulos hay en la imagen. Se encuentra el intervalo desde j_min hasta j_max 
	// y por lo tanto se obtienen los valores correspondientes en mp.ar, que a su vez se multiplican por el valor de la 
	// imagen en ese mismo punto.
	
    j_min=mp->ia[i];
    j_max=mp->ia[i+1];
	for(j=j_min;j<j_max;j++)  
    {
		lc[i]+= exp(-(q[mp->ja[j]])*(mp->ar[j]));
    }    
}

}

// Utiliza la función anterior para proyectar un objeto y producir un sinograma de factores de atenuación (que se ha de multiplicar X emisión)

int fwdproj_att(char *matriz_in,char *objeto_hdr)
{
	char matriz[180],proy_objeto_hdr[180];
	int i,j,k,i_read,i_write,NbOS,Nvox,sb,Nitems_max;
	sparse_i4 mp;
	struct imagen ima,sino,sino_global;
	
	// relleno estructuras con valores del objeto, después se modifican 
	lee_imagen_hdr(objeto_hdr,&ima);
	lee_imagen_hdr(objeto_hdr,&sino);
	lee_imagen_hdr(objeto_hdr,&sino_global);
	
	// Lectura de la matriz de pesos inicial para obtener parámetros	
	sb=0;
	sprintf(matriz,"%s.%d",matriz_in,sb);
	lee_parametros_pesos(matriz,&mp);
	NbOS=mp.Nbins*mp.Nang*num_imagenes(mp.ndet_z,mp.axial_diff)/mp.Nsubsets;
	Nvox=mp.n_fil*mp.n_col*mp.Ntallsima;
	check_nitems_max_mp_pet(matriz_in,mp.Nsubsets,&Nitems_max);
	mp.ne=Nitems_max;
	
	
	// controlando que coincidan las dimensiones de la matriz y los sinogramas
	if(mp.n_fil!=ima.nfil)
	{
		printf("\n\n\nERROR\nObjeto:%d y Matrices:%d\n",ima.nfil,mp.n_fil);
		return 10;
	}
	if(mp.n_col!=ima.ncol) 
	{
		printf("\n\n\nERROR\nObjeto:%d y Matrices:%d\n",ima.ncol,mp.n_col);
		return 10;
	}
	if(mp.Ntallsima!=ima.nima)
	{
		printf("\n\n\nERROR\nObjeto:%d y Matrices:%d\n",ima.nima,mp.Ntallsima);
		return 10;
	}
	
	// genera salida parcial (para un subset)
	sino.nfil=mp.Nang/mp.Nsubsets;
	sino.ncol=mp.Nbins;
	sino.nima=num_imagenes(mp.ndet_z,mp.axial_diff);
	sino.tpix=mp.tbin;
	sino.tcorte=mp.tdet_z;
	sino.offset=0;
	
	// genera salida global (todos los subsets)
	sino_global.nfil=mp.Nang;
	sino_global.ncol=mp.Nbins;
	sino_global.nima=num_imagenes(mp.ndet_z,mp.axial_diff);
	sino_global.tpix=mp.tbin;
	sino_global.tcorte=mp.tdet_z;
	sino_global.offset=0;
	
	// reserva memoria para las tres componentes de la matriz de pesos
	if((mp.ar=(float *)calloc(mp.ne,sizeof(float)))==NULL) return 5;  
	if((mp.ja=(int *)calloc(mp.ne,sizeof(int)))==NULL) return 5;
	if((mp.ia=(int *)calloc(NbOS+1,sizeof(int)))==NULL) return 5;
	
	// reserva memoria para sinogramas
	if((sino.datos=(float *)calloc(NbOS,sizeof(float)))==NULL) return 5;   
	if((sino_global.datos=(float *)calloc(NbOS*mp.Nsubsets,sizeof(float)))==NULL) return 5;
	
	
	for(sb=0;sb<mp.Nsubsets;sb++)
	{
		
		// Lectura de la matriz de pesos
		sprintf(matriz,"%s.%d",matriz_in,sb);
		lee_pesos_xOS(&mp,matriz);
		
		// mensaje
		printf("Proyectando %s\n",matriz);
		
		// forward projector - se proyecta el objeto
		//fwd_proj_attenuation(sino.datos,ima.datos,&mp);
		fwd_proj(sino.datos,ima.datos,&mp);
		
		for(k=0;k<sino.nima;k++)
		{
			for(j=0;j<sino.nfil;j++)
			{
				//printf("Salta: %d bloques, %d cortes y %d filas\n",sb,k,(mp.Nsubsets*j)+sb);
				for(i=0;i<sino.ncol;i++)
				{
					i_read=(k*sino.nfil*sino.ncol)+(j*sino.ncol)+i;
					i_write=(k*sino_global.nfil*sino.ncol)+(((mp.Nsubsets*j)+sb)*sino.ncol)+i;
					sino_global.datos[i_write]=sino.datos[i_read];
				}
			}
		}
		
	}

	sprintf(proy_objeto_hdr,"proy_global_%s",objeto_hdr);
	guarda_imagen_hdr(proy_objeto_hdr,&sino_global);
	return 0;
}


/*=========================== BACK_PROJ =====================*/

void back_proj(float *lc,float *q2,sparse_i4 *mp)
{
register int i,fila,npixels;

npixels=mp->n_fil*mp->n_col*mp->Ntallsima;

for(i=0;i<npixels;i++) q2[i]=0.;
fila=0;

for(i=0;i<mp->ne;i++)
{
   // se recorren todos los puntos con valores no nulos. Los puntos de la imagen 
   // se calculan a partir del valor del sinograma en la L.O.R. que corresponde y 
   // multiplicado por el peso para el punto imagen que se calcula. Para saber a qué 
   // L.O.R. corresponde cada peso se utiliza el 'if' que aumenta el valor de 'fila' cada 
   // vez que el índice 'i' es igual al mp.ia[fila+1]. 
   
   
   if(i>= mp->ia[fila+1]) fila++;
   q2[mp->ja[i]]+=mp->ar[i]*lc[fila];
}
}

/*=========================== BACK_PROJ =====================*/

void back_proj_sb(float *lc,float *q2,sparse_i4 *mp,int sb)
{    
register int i,fila,npixels;

npixels=mp->n_fil*mp->n_col*mp->Ntallsima;

for(i=0;i<npixels;i++) q2[i]=0.;
fila=0;

for(i=0;i<mp->ne;i++)
{
	// se recorren todos los puntos con valores no nulos. Los puntos de la imagen 
	// se calculan a partir del valor del sinograma en la L.O.R. que corresponde y 
	// multiplicado por el peso para el punto imagen que se calcula. Para saber a qué 
	// L.O.R. corresponde cada peso se utiliza el 'if' que aumenta el valor de 'fila' cada 
	// vez que el índice 'i' es igual al mp.ia[fila+1]. 
	
	
	if(i>= mp->ia[fila+1]) fila++;
	q2[mp->ja[i]]+=mp->ar[i]*lc[(sb*mp->Nbins)+fila];
}
}


/* =================== LEE PARAMS PESOS ======================= */
/* lectura de parámetros de la matriz */
void lee_parametros_pesos(char *matriu_pes_ini, sparse_i4 *mp_ini)
{
FILE *ini;
if((ini = fopen(matriu_pes_ini,"r"))==NULL) error_mp_pet(11,matriu_pes_ini);
fread (&mp_ini->n_fil,sizeof(size_t),1,ini);
fread (&mp_ini->n_col,sizeof(size_t),1,ini);
fread (&mp_ini->tpix,sizeof(size_t),1,ini);
fread (&mp_ini->Ntallsima,sizeof(size_t),1,ini);
fread (&mp_ini->lvox,sizeof(size_t),1,ini);
fread (&mp_ini->ndet_z,sizeof(size_t),1,ini);
fread (&mp_ini->tdet_z,sizeof(size_t),1,ini);
fread (&mp_ini->Nang,sizeof(size_t),1,ini);
fread (&mp_ini->Nbins,sizeof(size_t),1,ini);
fread (&mp_ini->tbin,sizeof(size_t),1,ini);
fread (&mp_ini->ddet,sizeof(size_t),1,ini);
fread (&mp_ini->Nsubsets,sizeof(size_t),1,ini);
fread (&mp_ini->axial_diff,sizeof(size_t),1,ini);
fread (&mp_ini->inv,sizeof(size_t),1,ini);
fread (&mp_ini->sym,sizeof(size_t),1,ini);
fread (&mp_ini->ne,sizeof(size_t),1,ini);
fclose(ini);
/* lo cierro y luego se vuelve a abrir en el calculo propio */

}


/*===========================  PESOSxOS ====================================*/

void lee_pesos_xOS(sparse_i4 *mp,char *nom_fitxer)
{
FILE *mat;
int NbOS;
int ni;
if((mat = fopen(nom_fitxer,"r"))==NULL) error_mp_pet(11,nom_fitxer);

printf("Reading SPARSE  %s ...\n",nom_fitxer);

fread (&mp->n_fil,sizeof(size_t),1,mat);
fread (&mp->n_col,sizeof(size_t),1,mat);
fread (&mp->tpix,sizeof(size_t),1,mat);
fread (&mp->Ntallsima,sizeof(size_t),1,mat);
fread (&mp->lvox,sizeof(size_t),1,mat);
fread (&mp->ndet_z,sizeof(size_t),1,mat);
fread (&mp->tdet_z,sizeof(size_t),1,mat);
fread (&mp->Nang,sizeof(size_t),1,mat);
fread (&mp->Nbins,sizeof(size_t),1,mat);
fread (&mp->tbin,sizeof(size_t),1,mat);
fread (&mp->ddet,sizeof(size_t),1,mat);
fread (&mp->Nsubsets,sizeof(size_t),1,mat);
fread (&mp->axial_diff,sizeof(size_t),1,mat);
fread (&mp->inv,sizeof(size_t),1,mat);
fread (&mp->sym,sizeof(size_t),1,mat);
fread (&mp->ne,sizeof(size_t),1,mat);

printf("Numer of elements  %d ...\n",mp->ne);

ni=num_imagenes(mp->ndet_z,mp->axial_diff);
// if(mp->sym==111) ni=((num_imagenes(mp->ndet_z,mp->axial_diff)-mp->ndet_z)/2)+mp->ndet_z;
printf("\nNumber of images to read:%d\n",ni);
if(mp->inv==0) 
{
	NbOS=(mp->Nbins*mp->Nang*ni)/mp->Nsubsets;
}
else if(mp->inv==1) 
{
	printf("ERROR!, Lectura no definida para matriz invertida");
	exit(0);
}
else 
{
	printf("ERROR!, Campo INV de mp no definido");
	exit(0);
}



fread((int *)mp->ia,sizeof(int),(size_t)(NbOS+1),mat);
fread ((int*)mp->ja,sizeof(int),mp->ne,mat);
fread((float*)mp->ar,sizeof(float),mp->ne,mat);
printf("Elements readed: %d en mp.ia, %d en mp.ja, %d en mp.ar\n",NbOS+1,mp->ne,mp->ne);
printf("\t\tExamples, mp.ia:%d, mp.ja:%d, mp.ar:%f\n",mp->ia[0],mp->ja[0],mp->ar[0]);
printf("\t\tExamples, mp.ia:%d, mp.ja:%d, mp.ar:%f\n",mp->ia[1],mp->ja[1],mp->ar[1]);
printf("\t\tExamples, mp.ia:%d, mp.ja:%d, mp.ar:%f\n",mp->ia[2],mp->ja[2],mp->ar[2]);
printf("\t\tExamples, mp.ia:%d, mp.ja:%d, mp.ar:%f\n",mp->ia[3],mp->ja[3],mp->ar[3]);
printf("   done:\n");
fclose (mat);  
}



/*===========================  PESOSxOS ====================================*/

void guarda_pesos_xOS(sparse_i4 *mp,char *nom_fitxer)
{
FILE *mat;
int Nb_OS,ni;

printf("Saving SPARSE  %s ...",nom_fitxer);

if(!mp->sym) mp->sym=0;

if(( mat = fopen(nom_fitxer,"wb"))==NULL) error_mp_pet(11,nom_fitxer);
fwrite (&mp->n_fil,sizeof(size_t),1,mat); // integer 
fwrite (&mp->n_col,sizeof(size_t),1,mat); // integer 
fwrite (&mp->tpix,sizeof(size_t),1,mat); // float 
fwrite (&mp->Ntallsima,sizeof(size_t),1,mat); // integer 
fwrite (&mp->lvox,sizeof(size_t),1,mat); // float 
fwrite (&mp->ndet_z,sizeof(size_t),1,mat); // integer 
fwrite (&mp->tdet_z,sizeof(size_t),1,mat); // float 
fwrite (&mp->Nang,sizeof(size_t),1,mat); // integer 
fwrite (&mp->Nbins,sizeof(size_t),1,mat); // integer 
fwrite (&mp->tbin,sizeof(size_t),1,mat); // float 
fwrite (&mp->ddet,sizeof(size_t),1,mat); // float 
fwrite (&mp->Nsubsets,sizeof(size_t),1,mat); // integer 
fwrite (&mp->axial_diff,sizeof(size_t),1,mat); // integer
fwrite (&mp->inv,sizeof(size_t),1,mat); // integer
fwrite (&mp->sym,sizeof(size_t),1,mat); // integer
fwrite (&mp->ne,sizeof(size_t),1,mat); // integer	

printf("Numer of elements  %d ...\n",mp->ne);

ni=num_imagenes(mp->ndet_z,mp->axial_diff);
// if(mp->sym==111) ni=((num_imagenes(mp->ndet_z,mp->axial_diff)-mp->ndet_z)/2)+mp->ndet_z;

Nb_OS=(mp->Nbins*mp->Nang*ni)/mp->Nsubsets; //bins number per subset (formato más común)
if (mp->inv==1) Nb_OS=mp->n_fil*mp->n_col*mp->Ntallsima; //voxel number

fwrite ((int*)mp->ia,sizeof(int),(size_t)(Nb_OS+1+1),mat);		
fwrite ((int*)mp->ja,sizeof(int),(size_t)mp->ne,mat);
fwrite ((float*)mp->ar,sizeof(float),(size_t)mp->ne,mat);
fclose (mat);
printf("Elements saved: %d en mp.ia, %d en mp.ja, %d en mp.ar\n",Nb_OS+1,mp->ne,mp->ne);
printf("\t\tExamples, mp.ia:%d, mp.ja:%d, mp.ar:%f\n",mp->ia[0],mp->ja[0],mp->ar[0]);
printf("\t\tExamples, mp.ia:%d, mp.ja:%d, mp.ar:%f\n",mp->ia[1],mp->ja[1],mp->ar[1]);
printf("\t\tExamples, mp.ia:%d, mp.ja:%d, mp.ar:%f\n",mp->ia[2],mp->ja[2],mp->ar[2]);
printf("\t\tExamples, mp.ia:%d, mp.ja:%d, mp.ar:%f\n",mp->ia[3],mp->ja[3],mp->ar[3]);
printf("   done:\n");
}

/*====================================== Obtiene imagen sensibilidad =====================*/

void ima_sensibilidad(struct imagen *sens,sparse_i4 *mp)
{
int i;

// inicializa vector  
for(i=0;i<mp->n_fil*mp->n_col*mp->Ntallsima;i++) 
{
	sens->datos[i]=0.;
}
// se recorre el vector mp.ja que contiene posiciones de la imagen con contribución no nula y se va sumando
// su contribución (contenida en la posición correspondiente de mp.ar)
for(i=0;i<mp->ne;i++) 
{
	sens->datos[mp->ja[i]]+=mp->ar[i];
}
}


/*====================================== Integra SRM para todos los bins =====================*/

void integral_matrix_bins(sparse_i4 *mp, float  *vector)
{
int i;

// inicializa vector  
for(i=0;i<mp->n_fil*mp->n_col*mp->Ntallsima;i++) 
{
	vector[i]=0.;
}
// se recorre el vector mp.ja que contiene posiciones de la imagen con contribución no nula y se va sumando
// su contribución (contenida en la posición correspondiente de mp.ar)
for(i=0;i<mp->ne;i++) 
{
	vector[mp->ja[i]]+=mp->ar[i];
}
}

/*====================================== Obtiene imagen sensibilidad =====================*/

void ima_sensibilidad_mp(struct imagen *sens,char *matriu_pes)
{
sparse_i4 mp;

lee_parametros_pesos(matriu_pes, &mp);
if((mp.ar=(float *)calloc(mp.ne,sizeof(float)))==NULL) error_mp_pet(5, "mp.ar");  
if((mp.ja=(int *)calloc(mp.ne,sizeof(int)))==NULL) error_mp_pet(5, "mp.ja");
if((mp.ia=(int *)calloc((mp.n_fil*mp.n_col*mp.Ntallsima)+1,sizeof(int)))==NULL) error_mp_pet(5, "mp.ia"); 
lee_pesos_xOS(&mp,matriu_pes);

ima_sensibilidad(sens,&mp);
printf("\n********************************************************************************\n");
}



/*=========================================== ERROR =========================*/

void error_mp_pet(int nerr,char *text)
{
switch(nerr){
	case 5: printf("\n\nError reserva memoria: %s",text);break;
    case 11: printf("\n\nError al obrir fitxer de la matriu de pesos");
             printf("\nLa matriu de pesos  %s  no existeix\n",text); break;
    default: printf("\n\nError del numero d'error en la funcio error_mp_pet()"); 
    }
    
exit(0);
}    
