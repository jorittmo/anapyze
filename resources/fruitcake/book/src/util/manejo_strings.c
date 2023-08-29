/*Funciones para manejar strings paguiar marzo 2006*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <util/manejo_strings.h>
#include <util/nrutil.h>
#include <util/inpout.h>

/* en ntokens se guarda el numero de tokens en el archivo ascii que se le pasa  paguiar marzo 2006*/

void cuenta_tokens_ascii(const char *nombre_fichero,int *ntokens)
{
	FILE *file;
	char buf[20000];
	char *token;
	int i;

	i=0;
// se abre el archivo que contiene el Header
	if((file=fopen(nombre_fichero,"r"))==NULL) 
	{
		printf("Error al abrir el archivo ASCII\n");
		exit(1);
	}

	while((fgets(buf,20000,file))!=NULL)
	{
		token=strtok(buf,"\t \n ");
		//printf("token primero:%s\n",token);
		if(token!=NULL) 
		{
			i++;
			//printf("Suma: %d\n",i);
		}
		
		while(token!=NULL)
		{
			token=strtok(NULL,"\t \n ");	
			//printf("token segundo:%s\n",token);

			if(token!=NULL) 
			{
				i++;				
				//printf("Suma: %d\n",i);
			}
		}
	}

*ntokens=i; /*numero de tokens separados*/
}


/*separa un archivo ascii en tokens (de \t y \n) ascii: nombre_fichero char *token es vector donde se almacenaran los tokens,*/
void separa_tokens_ascii(const char *nombre_fichero,char **palabras)
{
FILE *file;
char buf[20000];
char *token;
int i;

i=0;
// se abre el archivo que contiene el Header
if((file=fopen(nombre_fichero,"r"))==NULL) 
{
	printf("Error al abrir el archivo ASCII\n");
	exit(1);
}

while((fgets(buf,20000,file))!=NULL)
{
	token=strtok(buf,"\t \n ");
		
	if(token!=NULL) 
	{
		strcpy(palabras[i],token);
		//printf("%d) token:%s\n",i,token);
		i++;
	}

	while(token!=NULL)
	{
		token=strtok(NULL,"\t \n ");	
		if(token!=NULL) 
		{
			strcpy(palabras[i],token);
			//printf("%d) token:%s\n\n",i,token);
			i++;

		}
	}

}
}

/* en ntokens se guarda el numero de tokens en el string que se le pasa  paguiar mayo 2006*/



void cuenta_tokens(char *buf,int *ntokens,const char *separador)
{
	char *token;
	int i;

	i=0;

	token=strtok(buf,separador);	
	if(token!=NULL) 
	{
		i++;
	}
	while(token!=NULL)
	{
		token=strtok(NULL,separador);	

		if(token!=NULL) 
		{
			i++;
		}
	}
	*ntokens=i;
}



void separa_tokens(char *buf,char **palabras,const char *separador)
{
char *token;
int i;

i=0;
token=strtok(buf,separador);
strcpy(palabras[i],token);	
if(token!=NULL) 
{
	i++;
}
while(token!=NULL)
{
	token=strtok(NULL,separador);	

	if(token!=NULL) 
	{
		strcpy(palabras[i],token);
		i++;
	}
}
}
