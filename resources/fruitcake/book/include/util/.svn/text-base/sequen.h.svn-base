/*******************************************************************************
*                                                                              *
*          AQUESTA LLIBRERIA CONTE UNA SERIE DE RUTINES PER CONVERTIR          *
*          SEQUENCIES D'IMATGES EN UNA SOLA IMATGE, PER SEPARAR UNA IMATGE     *
*          D'UNA SEQUENCIA O UNA LLESCA O TALL D'UNA SEQUENCIA D'IMATGES       *
*          DE FITXERS DE QUALSEVOL FORMAT. FORMATS 1b, 2b, i4, fl, db.         *
*									       *
*******************************************************************************/
#ifndef SEQUEN__H
#define SEQUEN__H

void separa_imatge_1b(unsigned char *seq,unsigned char *q,int NFIL,int NCOL,int NIMA,int n);
void separa_imatge_2b(unsigned short int *seq,unsigned short int *q,int NFIL,int NCOL,int NIMA,int n);
void separa_imatge_i4(int *seq,int *q,int NFIL,int NCOL,int NIMA,int n);
void separa_imatge_fl(float *seq,float *q,int NFIL,int NCOL,int NIMA,int n);
void separa_imatge_db(double *seq,double *q,int NFIL,int NCOL,int NIMA,int n);
void separa_llesca_1b(unsigned char *seq,unsigned char *q,int NFIL,int NCOL,int NIMA,int n);
void separa_llesca_2b(short int *seq,short int *q,int NFIL,int NCOL,int NIMA,int n);
void separa_llesca_i4(int *seq,int *q,int NFIL,int NCOL,int NIMA,int n);
void separa_llesca_fl(float *seq,float *q,int NFIL,int NCOL,int NIMA,int n);
void separa_llesca_db(double *seq,double *q,int NFIL,int NCOL,int NIMA,int n);
void sequencia_a_ima_1b(unsigned char *seq,unsigned char *q,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay); 
void sequencia_a_ima_2b(short int *seq,short int *q,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay);
void sequencia_a_ima_i4(int *seq,int *q,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay); 
void sequencia_a_ima_fl(float *seq,float *q,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay); 
void sequencia_a_ima_db(double *seq,double *q,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay); 
void ima_a_seq_fl(float *q,float *seq,int NFIL,int NCOL,int NIMA,int cada_n,int primera,int nimax,int nimay);
void trasposa_matriu_fl(float *q1,float *q2,int Nfil,int Ncol);
void copia_vectors_fl(float *seq,float *seq2,int Nitems);

#endif //SEQUEN__H
