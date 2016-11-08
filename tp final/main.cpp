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
#define G3D(V,X,Y,Z)  V[(Z) * ((ii) * (jj)) + (Y) * (ii) + (X)]
#define S3D(V,X,Y,Z,S)  V[(Z) * ((ii) * (jj)) + (Y) * (ii) + (X)]=S
#define CREATEM3D(ii,jj,kk) std::vector<double>((ii)*(jj)*(kk))
#define VECTOR3D std::vector<double>


using namespace std;


//CONSTANTES

/*
%Paciente 1
%Glioma_1i
%Corrijo ubicacion del tumor y nivel de deteccion del tumor para
%diagnostico y letalidad
*/
double dt = 0.1;//dias
double h = 1;//mm, malla de 181x217x181 mm (18x22x18 cm)
int ii = 181; //En la imagen final se ven como columnas x) 
int jj = 217; //En la imagen final se ven como filas y)
int kk = 181; // (slices,z)
double C0 = 0; //concentracion inicial de celulas tumorales (por mm3)
double C_max = 1e8; //concentracion maxima de celulas (por mm3)
double C_mig = 1e7; //concentracion de celulas a la que comienza la migracion (por mm3)
double IM = 5; //indice mitotico
int nn = 5000; //corrida de 500 dias
int max_iter = 100;


vector<int> callo({508,509,510,514,515,651,652,657,694,695,711,712,713,714,715});

vector<int>  tracto_opt({273,359,361,624,633,638,639,641});//8
vector<int>  tallo({5,6,71,72,215,216,440,459,502,503,574,576});//12
vector<int>  cerebelo({1,2,3,4,61,62,63,64,65,66,67,68,80,81,82,83,84,85,86,87,88,89,114,115,116,117,118,119,120,157,158,159,160,161,162,163,164,165,211,212,321,322,564});//43
vector<int>  estriado({349,350,351,435,436,464,465,492,493,494,495,496,571,598,601,692,693,379,382,596,597}); 
vector<int>  globo({452,455,456,457,593,594,595,341,342,453,454,379,382,452,455,456,457,501,640,181,182,183,219,356,358,449,450,513});
vector<int>  medula({5,6,71,72,215,216,341,342,343,344,354,355,437,438,440,453,454,498,499,574});
vector<int> x(ii);
vector<int> y(jj);
vector<int> z(kk);
/*%Hago mapa de coeficientes
disp(datestr(now,'HH:MM AM'));
cerebro=importdata('Cerebro.mat');
talairach=importdata('Talairach_conv.mat');
*/

/*inicializacion de matrices*/
VECTOR3D cerebro = CREATEM3D(ii,jj,kk);
VECTOR3D talairach = CREATEM3D(ii,jj,kk);

VECTOR3D p = CREATEM3D(ii,jj,kk);
VECTOR3D D = CREATEM3D(ii,jj,kk);
VECTOR3D C = CREATEM3D(ii,jj,kk);
VECTOR3D B = CREATEM3D(ii,jj,kk);
VECTOR3D P_optimizado = CREATEM3D(ii,jj,kk);
VECTOR3D M_optimizado = CREATEM3D(ii,jj,kk);


//FUNCIONES

bool pertenece(vector<int> v, int val){
	
	for(int i=0;i<v.size();i++){
		
			if(v[i]==val) return true;
		
		}
	return false;
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

void ReadDifussionData(string dataFile, int tamX, int tamY, int tamZ, int originX, int originY, int originZ, VECTOR3D &difusionMat){

    printf ("Reading difussion data file %s, section (%u,%u,%u),(%u,%u,%u)\n", dataFile.c_str(), originX, originY, originZ, originX + tamX, originY+tamY, originZ+tamZ);

    
    ifstream file ( dataFile );
    string value;
    getline ( file, value);
    int x = 0;
    int y = 0;
    int z = 0;
    
 
    while ( x < originX + tamX){
                getline ( file, value, ',' );
                x = atoi(value.c_str()) - 1;
                getline ( file, value, ',' );
        
                y = atoi(value.c_str()) - 1;
                getline ( file, value, ',' );
        
                z = atoi(value.c_str()) - 1 ;
                getline ( file, value );
        
                double s = atof(value.c_str());
        
                //printf("%u,%u,%u\n", x,y,z);
                if (originX <= x && x - originX < tamX  && originY <= y && y - originY < tamY && originZ <= z && z - originZ < tamZ)
                    S3D(difusionMat,x - originX,y - originY,z - originZ,s);
            
    }
    printf ("Difussion data read OK\n");

}

void TransformDifusion(){
	
for(int k=0;k<kk;k++){
    for(int j=0;j<jj;j++){
        for(int i=0;i<ii;i++){	
			 
				if((G3D(cerebro,i,j,k)>110) && G3D(cerebro,i,j,k)<=225){ //sustancia blanca
					S3D(p,i,j,k,0.107); //proliferacion neta, 1/dia
					
						if(pertenece(callo, G3D(talairach,i,j,k))) //cuerpo calloso
							S3D(p,i,j,k,0.306); //un 20% mas
				
							
						if (pertenece(tracto_opt,G3D(talairach,i,j,k))) //tracto optico
							S3D(D,i,j,k,0.306); //un 20% mas
							
						if(pertenece(tallo,G3D(talairach,i,j,k))) //tallo cerebral, medula, protuberancia, mesensefalo
							S3D(D,i,j,k,0.204); //un 20% menos

							S3D(D,i,j,k,0.255); //igracion neta, mm2/dia
							
				}else{
					if((G3D(cerebro,i,j,k)>=75) && G3D(cerebro,i,j,k)<=110){ //sustancia gris
						S3D(p,i,j,k,0.107); 
							if(pertenece(cerebelo,G3D(talairach,i,j,k)))//cerebelo
								S3D(D,i,j,k,0);
								
							if(pertenece(estriado,G3D(talairach,i,j,k)))//nucleo estriado (caudado + putamen)
								S3D(D,i,j,k,0.0408); //un 20% menos
								
							if(pertenece(globo,G3D(talairach,i,j,k)))//globo palido, sustancia nigra, nucleo subtalamico, nucleo lentiforme, amigdala, claustrum
								S3D(D,i,j,k,0.0408); //un 20% menos
								
							if(pertenece(medula,G3D(talairach,i,j,k))) //tallo cerebral, medula, protuberancia, mesensefalo
								S3D(D,i,j,k,0.0408); //un 20% menos
								
								
								S3D(D,i,j,k,0.051);
					            
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
     
//C = ones(ii,jj,kk)*C0;
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


int main(){
	
	


    ReadDifussionData("./Cerebro.csv", ii, jj, kk, 0, 0, 0, D);//lee del archivo a matriz
    ReadDifussionData("./Talairach.csv", ii, jj, kk, 0, 0, 0, talairach);//lee del archivo a matriz

    
	
	TransformDifusion();//inicializa valores de la matriz

	inicializarCondiciones();

	
	return 0;
	}
