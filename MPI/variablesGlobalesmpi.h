#ifndef _VARIABLES_GLOBALES_H_
#define _VARIABLES_GLOBALES_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <assert.h>
#include <ctime>
#include <vector>
#include <stddef.h>
#include <sys/time.h>
#include <unistd.h>
#include <vector>
#include <time.h>
#include <string.h>
#include <cstdlib>
#include <mpi.h>
#include <stdarg.h>





#define G3D(V,X,Y,Z)  V[(Z) * ((ii) * (jj)) + (Y) * (ii) + (X)]
#define S3D(V,X,Y,Z,S)  V[(Z) * ((ii) * (jj)) + (Y) * (ii) + (X)]=S
#define CREATEM3D(ii,jj,kk) std::vector<double>((ii)*(jj)*(kk))
#define VECTOR3D std::vector<double>

using namespace std;


double f0 = (181+2)/181;
double f1 = (217+2)/217;

//Cuando hacés crecer por un factor f > 1

int ii =  f0 * 181;
int jj =  f1 * 217;
int kk=   f0* 181 ;

int io = f0  * 32;
int jo = f1  * 98;
int ko = f0  * 85;


//Variables globales porque se me hace mas comodo :B
ofstream info;
ofstream aux;
ofstream datos;
double dt = 0.1;//dias
double h = 1;//mm, malla de 181x217x181 mm (18x22x18 cm)
//int ii = 181+2; //En la imagen final se ven como columnas x) 
//int jj = 217+2; //En la imagen final se ven como filas y)
//int kk = 181+2; // (slices,z)
double C0 = 0; //concentracion inicial de celulas tumorales (por mm3)
double C_max = 1e8; //concentracion maxima de celulas (por mm3)
double C_mig = 1e7; //concentracion de celulas a la que comienza la migracion (por mm3)
double IM = 5; //indice mitotico
int nn = 3000; //corrida de 500 dias
int max_iter = 100;
int dia = 1;
int migracion = 0; // 0 para tumor benigno y 1 para tumor maligno
int diagnostico = 3188; //Para diagnostico a un diametro de 18.26 mm
//int io = 32;
//int jo = 98;
//int ko = 85;
int cantidad1 = 0; //para deteccion de celulas tumorales
int cantidad2 = 0; //para diagnostico
int cantidad3 = 0; //para letalidad



int myints1[] = {508,509,510,514,515,651,652,657,694,695,711,712,713,714,715};
vector<int> callo (myints1, myints1 + sizeof(myints1) / sizeof(int) );


int myints2[] = {273,359,361,624,633,638,639,641};
vector<int> tracto_opt (myints2, myints2 + sizeof(myints2) / sizeof(int) );

int myints3[] = {5,6,71,72,215,216,440,459,502,503,574,576};
vector<int> tallo (myints3, myints3 + sizeof(myints3) / sizeof(int) );



int myints4[] = {1,2,3,4,61,62,63,64,65,66,67,68,80,81,82,83,84,85,86,87,88,89,114,115,116,117,118,119,120,157,158,159,160,161,162,163,164,165,211,212,321,322,564};
vector<int> cerebelo (myints4, myints4 + sizeof(myints4) / sizeof(int) );



int myints5[] = {349,350,351,435,436,464,465,492,493,494,495,496,571,598,601,692,693,379,382,596,597};
vector<int> estriado (myints5, myints5 + sizeof(myints5) / sizeof(int) );


int myints6[] = {452,455,456,457,593,594,595,341,342,453,454,379,382,452,455,456,457,501,640,181,182,183,219,356,358,449,450,513};
vector<int> globo (myints6, myints6+ sizeof(myints6) / sizeof(int) );


int myints7[] = {5,6,71,72,215,216,341,342,343,344,354,355,437,438,440,453,454,498,499,574};
vector<int> medula (myints7, myints7 + sizeof(myints7) / sizeof(int) );



VECTOR3D cerebro = CREATEM3D(ii,jj,kk);
VECTOR3D talairach = CREATEM3D(ii,jj,kk);

VECTOR3D p = CREATEM3D(ii,jj,kk);
VECTOR3D D = CREATEM3D(ii,jj,kk);
VECTOR3D C = CREATEM3D(ii,jj,kk);
VECTOR3D P_optimizado = CREATEM3D(ii,jj,kk);
VECTOR3D M_optimizado = CREATEM3D(ii,jj,kk);
VECTOR3D C_k2 = CREATEM3D(ii,jj,kk);
VECTOR3D C_k1 = CREATEM3D(ii,jj,kk);
VECTOR3D P = CREATEM3D(ii,jj,kk);

VECTOR3D M = CREATEM3D(ii,jj,kk);


//Variables MPI
int world_rank;
int world_size;
int chunkSize = 1;
int scatterSize = chunkSize - 1;
int workingSize = chunkSize + 1;
//FUNCIONES MPI

void screenMessage(const char *fmt, ...) {
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if (world_rank == 0){
        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
    }
}

void doHaloTransfer(VECTOR3D &slice, int rank, int numProcs, int workingSliceSize){
     //printf("[%u] Entering Halo\n", rank);
    
//MPI_Sendrecv(&slice[(ii  * jj) * (workingSliceSize - 2)], ii  * jj, MPI_DOUBLE, rank + 1, 0,
  //        &slice[(ii  * jj) * (workingSliceSize-1)], ii  * jj, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


     if (rank % 2 == 0){//CASO PARES
        //Z0-Z1-Z2-Z3-Z4  Z3-Z4-Z5-Z6-Z7  Z6-Z7-Z8-Z9-ZA  Z9-ZA-ZB-ZC ....
        // SENDRECV Z4(0) <- Z4(1)
        if (rank > 0){
		MPI_Sendrecv(&slice[ii  * jj], ii  * jj, MPI_DOUBLE, rank - 1, 0, &slice[0], ii  * jj, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // SENDRECV Z4(0) <- Z4(1)
        if (rank < numProcs - 1){

		MPI_Sendrecv(&slice[(ii  * jj) * (workingSliceSize - 2)], ii  * jj, MPI_DOUBLE, rank + 1, 0, &slice[(ii  * jj) * (workingSliceSize-1)], ii  * jj, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }else{// CASO IMPARES

        // SENDRECV Z3(1) <- Z3(0)
        if (rank > 0){
	MPI_Sendrecv(&slice[ii  * jj], ii  * jj, MPI_DOUBLE, rank - 1, 0,&slice[0], ii  * jj, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }

        // SENDRECV Z4(0) <- Z4(1)
        if (rank < numProcs - 1){
	MPI_Sendrecv(&slice[(ii  * jj) * (workingSliceSize - 2)], ii  * jj, MPI_DOUBLE, rank + 1, 0,
                &slice[(ii  * jj) * (workingSliceSize-1)], ii  * jj, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }

    }
     //printf("[%u] Exit Halo\n", rank);
    MPI_Barrier(MPI_COMM_WORLD);

}




//FUNCIONES


bool pertenece(vector<int> v, int val){
	
	for(int i=0;i<v.size();i++){
		
			if(v[i]==val) return true;
		
		}
	return false;
}

void copyMatrix(VECTOR3D &mat1, VECTOR3D &mat2){

	for(int k=0;k<chunkSize;k++){
		for(int i=0;i<ii;i++){
			for(int j=0;j<jj;j++){

				S3D(mat1,i,j,k,G3D(mat2,i,j,k));



			}
		}

	}


}

double restaMax(VECTOR3D &mat1, VECTOR3D &mat2){

	double max=0;
	double m1=0;
	double m2=0;
	double abs=0;
	for(int k=0;k<chunkSize;k++){
		for(int i=0;i<ii;i++){
			for(int j=0;j<jj;j++){

				m2=G3D(mat2,i,j,k);
				m1=G3D(mat1,i,j,k);
				if(m1>=m2){
					abs=m1-m2;

				}else{
					abs=m2-m1;

				}
				if(max<abs){
					max=abs;
				}

			}
		}

	}

	return max;


}

void imprimir_matriz(VECTOR3D &mat){

for(int i=0;i<ii;i++){

	for(int j=0;j<jj;j++){

		for(int k=0;k<kk;k++){

			if(!G3D(mat,i,j,k)==0)
			cout<<"valor: "<<G3D(mat,i,j,k)<<endl;

		}

	}

}
	

}
void ReadDifussionData(string dataFile, int tamX, int tamY, int tamZ, int originX, int originY, int originZ, VECTOR3D &difusionMat){

    printf ("Reading difussion data file %s, section (%u,%u,%u),(%u,%u,%u)\n", dataFile.c_str(), originX, originY, originZ, originX + tamX, originY+tamY, originZ+tamZ);

    
    ifstream file ( dataFile.c_str() );
    string value;
    getline ( file, value);
    int x = 0;
    int y = 0;
    int z = 0;

    while (  x < originX + tamX){
                getline ( file, value, ',' );
                x = atoi(value.c_str()) - 1;
                getline ( file, value, ',' );


                y = atoi(value.c_str()) - 1;
                getline ( file, value, ',' );
        
                z = atoi(value.c_str()) - 1 ;
                getline ( file, value );
        
                double s = atof(value.c_str());
		//cout<<s<<endl;
              //printf("%u,%u,%u\n", x,y,z);
               // if (originX <= x && x - originX < tamX  && originY <= y && y - originY < tamY && originZ <= z && z - originZ < tamZ){

                    S3D(difusionMat,x,y,z,s);
            	//}
    }

    printf ("Difussion data read OK\n");
}

void TransformDifusion(){
	
for(int k=0;k<kk;k++){
    for(int j=0;j<jj;j++){
        for(int i=0;i<ii;i++){
			 
				if((G3D(cerebro,i,j,k)>110) && G3D(cerebro,i,j,k)<=225){ //sustancia blanca

					S3D(p,i,j,k,0.107); //proliferacion neta, 1/dia
					
						if(pertenece(callo, G3D(talairach,i,j,k))){ //cuerpo calloso
							S3D(p,i,j,k,0.306); //un 20% mas
						}else if (pertenece(tracto_opt,G3D(talairach,i,j,k))){ //tracto optico
							S3D(D,i,j,k,0.306); //un 20% mas
						}else if(pertenece(tallo,G3D(talairach,i,j,k))){ //tallo cerebral, medula, protuberancia, mesensefalo
							S3D(D,i,j,k,0.204); //un 20% menos
						}else{
							S3D(D,i,j,k,0.255); //igracion neta, mm2/dia
						}
							
				}else{

					if((G3D(cerebro,i,j,k)>=75) && G3D(cerebro,i,j,k)<=110){ //sustancia gris
						S3D(p,i,j,k,0.107); 
							if(pertenece(cerebelo,G3D(talairach,i,j,k))){//cerebelo
								S3D(D,i,j,k,0);
							}else if(pertenece(estriado,G3D(talairach,i,j,k))){//nucleo estriado (caudado + putamen)
								S3D(D,i,j,k,0.0408); //un 20% menos
							}else if(pertenece(globo,G3D(talairach,i,j,k))){//globo palido, sustancia nigra, nucleo subtalamico, nucleo lentiforme, amigdala, claustrum
								S3D(D,i,j,k,0.0408); //un 20% menos
							}else if(pertenece(medula,G3D(talairach,i,j,k))){ //tallo cerebral, medula, protuberancia, mesensefalo
								S3D(D,i,j,k,0.0408); //un 20% menos
							}else{
								S3D(D,i,j,k,0.051);
							}
					            
					}else{
					
					S3D(p,i,j,k, 0);
					S3D(D,i,j,k,0);
					
					}	
				
				}
			}
		}
	}

}

void inicializarCondiciones(){
	
//Origen del tumor en:
int io = 32;
int jo = 98;
int ko = 85;
for(int k=0;k<kk;k++){
   for(int j=0;j<jj;j++){
       for(int i=0;i<ii;i++){
       S3D(C,i,j,k,C0);
       }
	}
}	  
     	S3D(C,io,jo,ko,1);

int dia = 1;
int migracion = 0; // 0 para tumor benigno y 1 para tumor maligno
int diagnostico = 3188; //Para diagnostico a un diametro de 18.26 mm

//Optimizaciones:
for(int k=1;k<kk-1;k++){ 
    for(int j=1;j<jj-1;j++){
        for(int i=1;i<ii-1;i++){
            S3D(P_optimizado,i,j,k,dt * G3D(p,i,j,k));
            S3D(M_optimizado,i,j,k,dt * G3D(D,i,j,k) / (h*h));
        }
    }
}

	
}

/*
void guardar_datos(int n,double error,int cantidad1,int cantidad2,int cantidad3,vector<int>& B,vector<int>& B_R, vector<int>& B_T){
	//para ver el dia actual
	time_t tiempo = time(0);
        struct tm *tlocal = localtime(&tiempo);
        char output[128];
        strftime(output,128,"%d/%m/%y %H:%M:%S",tlocal);
//	if (n%500==0){
//		grabar_matriz(C);
  //  }
    // Grabar info e informar por pantalla:
    if (G3D(C,io,jo,ko) > C_mig && migracion == 0){
            cout<< "comienza migracion"<<endl;
            migracion = 1; 
		}
    if(error > 10)
            cout<< "error = " << error << endl;
  	
    if (n % 10 ==0){

        printf("%s\n",output);
	cout<<"iteracion: "<< dia <<endl;
	info<< output;
	info<<" dia: "<< dia/10<<" ";
	info<< " error: " << error <<endl;
	
	}

        if (migracion == 1){
		info<<"En migración \n"; 
        }
        if (cantidad2 >= diagnostico){
		info<<"diagnosticado"<<endl; 
        }
        if (cantidad3 >= 179594){ //Para area letal (esfera de 70 mm de diametro) ponderada por areas vitales
		info<<"muerte del paciente"<<endl; 
        }

        for (int a =0;a<47;a++){
            if (B[a]>=1&& B_T[a]>0){ //si esa area de Brodmann no es nula
		
                B_R[a]=B[a]*100/B_T[a]; //la relativizo con respecto a la total del Talairach

		info<< "Brodmann: " << a << " = " << B[a]<<endl;
		info<< "Porcentaje: " <<  B_R[a]<<endl;	

            }
        }
        dia++; 
    
}
*/
void iteracion_de_convergencia(int n,int cantidad1,int cantidad2,int cantidad3,VECTOR3D &C,VECTOR3D &C_k1,VECTOR3D &C_k2){

	double max_error=0;
	if (dia <= 1500){
	max_error = 1;
	}else{
	max_error = 10;
	}
	int iter=0;
	double error = 1000;


	copyMatrix(C_k2,C);
	//while((iter < max_iter) && (error > max_error)){
		copyMatrix(C_k1,C_k2);
		for (int k=1;k<chunkSize+1;k++){
			for (int j=1;j<jj-1;j++){
				for (int i=1;i<ii-1;i++){

					//S3D(M,i,j,k,0);
					S3D(M,i,j,k,(G3D(C_k1,i+1,j,k)+G3D(C_k1,i-1,j,k)+G3D(C_k1,i,j+1,k)+G3D(C_k1,i,j-1,k)+G3D(C_k1,i,j,k+1)+G3D(C_k1,i,j,k-1)-6*G3D(C_k1,i,j,k)));

					S3D(C_k2,i,j,k, G3D(M,i,j,k));
					if (G3D(C_k2,i,j,k) > C_max){
						S3D(C_k2,i,j,k, C_max);
					}
					if (G3D(C_k2,i,j,k) < 0.00001){
						S3D(C_k2,i,j,k,0);
					}
					
				}
			}
		}
		error=restaMax(C_k1,C_k2);
		iter++;

	//}
    // Actualizar malla
	copyMatrix(C,C_k2);
	/*if( world_rank==0){
		//error = restaMax(C_k1,C_k2);
		guardar_datos(n,error,cantidad1,cantidad2,cantidad3,B,B_R, B_T);
	
	}*/





}

#endif
