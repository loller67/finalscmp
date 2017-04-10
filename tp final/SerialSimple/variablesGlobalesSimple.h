#ifndef _VARIABLES_GLOBALES_H_
#define _VARIABLES_GLOBALES_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <math.h>
#include <assert.h>
#include <chrono>
#include <ctime>
#include <vector>
#include <stddef.h>
#include <sys/time.h>
#include <unistd.h>
#include <vector>
#include <time.h>
#include <string.h>

#define G3D(V,X,Y,Z)  V[(Z) * ((ii) * (jj)) + (Y) * (ii) + (X)]
#define S3D(V,X,Y,Z,S)  V[(Z) * ((ii) * (jj)) + (Y) * (ii) + (X)]=S
#define CREATEM3D(ii,jj,kk) std::vector<double>((ii)*(jj)*(kk))
#define VECTOR3D std::vector<double>
using namespace std;

//Variables globales
int cantidad1=0;
int cantidad2=0;
int cantidad3=0;
ofstream info;
ofstream datos;
double dt = 0.1;//dias
double h = 1;//mm, malla de 181x217x181 mm (18x22x18 cm)
int ii = 181; //En la imagen final se ven como columnas x) 
int jj = 217; //En la imagen final se ven como filas y)
int kk = 181; // (slices,z)
double C0 = 0; //concentracion inicial de celulas tumorales (por mm3)
double C_max = 1e8; //concentracion maxima de celulas (por mm3)
double C_mig = 1e7; //concentracion de celulas a la que comienza la migracion (por mm3)
double IM = 5; //indice mitotico
int nn = 3000; //corrida de 500 dias
int max_iter = 100;
int dia = 1;
int migracion = 0; // 0 para tumor benigno y 1 para tumor maligno
int diagnostico = 3188; //Para diagnostico a un diametro de 18.26 mm
int io = 32;
int jo = 98;
int ko = 85;
double d=0.051*dt/(h*h);
double p=0.107*dt;

  vector<int> x(ii);
  vector<int> y(jj);
  vector<int> z(kk);
 

vector<int> B(47, 0);// vector que guarda areas de Brodmann absolutas
vector<int> B_R(47, 0); //vector que guarda areas de Brodmann relativizadas

VECTOR3D C = CREATEM3D(ii,jj,kk);
VECTOR3D C_k2 = CREATEM3D(ii,jj,kk);
VECTOR3D C_k1 = CREATEM3D(ii,jj,kk);


//FUNCIONES


//guarda en mat1 el contenido de mat2
void copyMatrix(VECTOR3D &mat1, VECTOR3D &mat2){

	for(int i=0;i<ii;i++){
		for(int j=0;j<jj;j++){
			for(int k=0;k<kk;k++){

				S3D(mat1,i,j,k,G3D(mat2,i,j,k));


			}
		}

	}

}

//calcula el valor de la resta maxima de todas las dimensiones de las dos matrices
double restaMax(VECTOR3D &mat1, VECTOR3D &mat2){

	double max=0;
	double m1=0;
	double m2=0;
	double abs=0;
	for(int i=0;i<ii;i++){
		for(int j=0;j<jj;j++){
			for(int k=0;k<kk;k++){

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

void cargar_vectores(){
	
	x[0]=0;
	y[0]=0;
	z[0]=0;

	for(int i=1;i<ii;i++){
		x[i]=x[i-1]+h;
	}
	for(int j=1;j<jj;j++){
		y[j]=y[j-1]+h;
	}
	for(int k=1;k<kk;k++){
    z[k]=z[k-1]+h;
	}

}
// imprime valores de la matriz, usado para debuggear
void imprimir_matriz(VECTOR3D &mat){
for(int i=0;i<ii;i++){

	for(int j=0;j<jj;j++){

		for(int k=0;k<kk;k++){


			cout<<"valor: "<<G3D(mat,i,j,k)<<endl;

		}

	}

}
	

}
//funcion para obtener la fecha
std::string GetLocalTime() {
    auto now(std::chrono::system_clock::now());
    auto seconds_since_epoch(
                             std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()));
    
    // Construct time_t using 'seconds_since_epoch' rather than 'now' since it is
    // implementation-defined whether the value is rounded or truncated.
    std::time_t now_t(
                      std::chrono::system_clock::to_time_t(
                                                           std::chrono::system_clock::time_point(seconds_since_epoch)));
    
    char temp[10];
    if (!std::strftime(temp, 10, "%H:%M:%S.", std::localtime(&now_t)))
        return "";
    
    return std::string(temp) +
    std::to_string((now.time_since_epoch() - seconds_since_epoch).count());
}


// guarda los valores de la matriz en el archivo .vtk
void dumpMatrixToVtk(VECTOR3D &mat, string fileId){
    cout << "Dumping to VTK..." << endl;
    
    long numDataPoints = ii * jj * kk;
    
    std::filebuf fbmc;
    fbmc.open ("./" + fileId + ".vtk",std::ios::out);
    std::ostream osmc(&fbmc);
    
    osmc << "# vtk DataFile Version 2.0" << endl;
    osmc << "Comment goes here" << endl;
    osmc << "ASCII" << endl;
    
    osmc << "DATASET STRUCTURED_POINTS" << endl;
    osmc << "DIMENSIONS    " <<  ii << " " << jj << " " << kk << endl;
    
    osmc << "ORIGIN    45.000   50.000   45.000" << endl;
    osmc << "SPACING    1.000   1.000   1.000" << endl;
    
    osmc << "POINT_DATA   " << numDataPoints << endl;
    osmc << "SCALARS scalars double" << endl;
    osmc << "LOOKUP_TABLE default" << endl;
    
    for (int i = 0; i < mat.size();i++)
        osmc << mat[i] << endl;
    
    
    fbmc.close();
}
//funcion para levantar las matrices de disco
void ReadDifussionData(string dataFile, int tamX, int tamY, int tamZ, int originX, int originY, int originZ, VECTOR3D &difusionMat){

    printf ("Reading difussion data file %s, section (%u,%u,%u),(%u,%u,%u)\n", dataFile.c_str(), originX, originY, originZ, originX + tamX, originY+tamY, originZ+tamZ);

    
    ifstream file ( dataFile );
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
//imprimir_matriz(difusionMat);
}



//funcion que setea condiciones iniciales
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
	
}


void grabar_matriz(VECTOR3D &mat){
//imprime en datos en algun orden la matriz
	
for(int i=0;i<ii;i++){
	for(int j=0;j<jj;j++){
		for(int k=0;k<kk;k++){
			//if(G3D(mat,i,j,k)!=0)
			datos<<G3D(mat,i,j,k)<<endl;
}
}
}

}

//funcion para guardar datos en disco
void guardar_datos(int n,double error){
	//para ver el dia actual
	time_t tiempo = time(0);
        struct tm *tlocal = localtime(&tiempo);
        char output[128];
        strftime(output,128,"%d/%m/%y %H:%M:%S",tlocal);
	//if (n%500==0){
	//	grabar_matriz(C);
    //}
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
	
	

		if (migracion == 1){
			info<<"En migración \n"; 
		}
		if (cantidad2 >= diagnostico){
			info<<"diagnosticado"<<endl; 
		}
		if (cantidad3 >= 179594){ //Para area letal (esfera de 70 mm de diametro) ponderada por areas vitales
			info<<"muerte del paciente"<<endl; 
		}

	}


        dia++; 
    
}





#endif
