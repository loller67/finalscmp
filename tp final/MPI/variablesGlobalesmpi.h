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
#include <mpi.h>
#include <stdarg.h>

#define G3D(V,X,Y,Z)  V[(Z) * ((ii) * (jj)) + (Y) * (ii) + (X)]
#define S3D(V,X,Y,Z,S)  V[(Z) * ((ii) * (jj)) + (Y) * (ii) + (X)]=S
#define CREATEM3D(ii,jj,kk) std::vector<double>((ii)*(jj)*(kk))
#define VECTOR3D std::vector<double>

using namespace std;

//Variables globales porque se me hace mas comodo :B

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
int cantidad1 = 0; //para deteccion de celulas tumorales
int cantidad2 = 0; //para diagnostico
int cantidad3 = 0; //para letalidad


  vector<int> callo({508,509,510,514,515,651,652,657,694,695,711,712,713,714,715});
  vector<int>  tracto_opt({273,359,361,624,633,638,639,641});//8
  vector<int>  tallo({5,6,71,72,215,216,440,459,502,503,574,576});//12
  vector<int>  cerebelo({1,2,3,4,61,62,63,64,65,66,67,68,80,81,82,83,84,85,86,87,88,89,114,115,116,117,118,119,120,157,158,159,160,161,162,163,164,165,211,212,321,322,564});//43
  vector<int>  estriado({349,350,351,435,436,464,465,492,493,494,495,496,571,598,601,692,693,379,382,596,597}); 
  vector<int>  globo({452,455,456,457,593,594,595,341,342,453,454,379,382,452,455,456,457,501,640,181,182,183,219,356,358,449,450,513});
  vector<int>  medula({5,6,71,72,215,216,341,342,343,344,354,355,437,438,440,453,454,498,499,574});
vector<int> aux;

  vector<int> x(ii);
  vector<int> y(jj);
  vector<int> z(kk);
  vector<int> v1({893,894});
  vector<int> v2({891,892,978});
  vector<int> v3({895,896,897,898,938,939,940,1019,1088,1091,1092});
  vector<int> v4({804,805,806,808,835,836,1003,1089,1090,1093,1099,1100,1104,1105});
  vector<int> v5({1064,1067,1068,1070,1083,1086,1095,1096,1097,1098});
  vector<int> v6({682,683,809,810,1009,1010,1011,1012,1022,1023,1024,1075,1076,1077,1078,1094});
  vector<int> v7({952,953,968,969,971,972,1039,1040,1043,1044,1065,1066,1080,1081,1087,1101,1102});
  vector<int> v8({1025,1026,1027,1028,1033,1034,1035,1036});
  vector<int> v9({813,814,815,816,841,899,900,901,902,903,943,944,981,982,983,984,1004});
  vector<int> v10({405,406,407,408,410,411,412,413,471,483,484,516,517,754,755});
  vector<int> v11({102,103,109,142,145,152,153,191,192,197,198,235,236,300,301});
//vector<int> v12
  vector<int> v13({288,292,294,295,296,297,365,366,367,368,369,376,400,444,447,589,590,733,734,735});
//vector<int> v14
//vector<int> v15
//vector<int> v16
  vector<int> v17({242,245,303,305,475,477,540});
  vector<int> v18({203,205,247,249,250,254,304,306,307,474,478,487,488,537,631,632,758,760,826,827,857,858});
  vector<int> v19({258,259,263,264,265,266,315,316,414,415,417,418,421,423,425,426,427,428,479,480,566,568,699,700,818,819,820,843,844,848,849,986,987});
  vector<int> v20({8,13,15,16,18,19,23,26,29,46,53,123,124,135,167,171,172,217,218});
  vector<int> v21({77,78,79,97,275,276,286,352});
  vector<int> v22({442,443,445,446,458,460,587,588,634,764,765,766,767});
  vector<int> v23({706,707,762,763,829,830,934,935});
  vector<int> v24({469,470,941,942,1020,1021});
  vector<int> v25({226,227,234,237,281,282,380,381});
//vector<int> v26
  vector<int> v27({432,433});
  vector<int> v28({94,95,169,170,345});
  vector<int> v29({670,671});
  vector<int> v30({429,430,431,481,482,660,661,663,665,669});
  vector<int> v31({756,757,823,824,856,859,863,905,909,915,916,920,923,932,994,995,998,1016,1017,1045,1049,1050,1053});
  vector<int> v32({395,396,610,612,945,946,947,948,1006,1008,1029,1030,1031,1032});
  vector<int> v33({838,839});
  vector<int> v34({177,178,179,180,279,280,283,285,287});
  vector<int> v35({140,141,272});
  vector<int> v36({75,76,129,130,136,166});
  vector<int> v37({208,209,213,214,267,268,269,270,271,323,324,328,329,330,331,332,333,334,420,422,558,560,627,628});
  vector<int> v38({33,34,35,36,38,41,44,45,59,60,101,221,222});
  vector<int> v39({702,703,704,705,708,709,865,906,907,908,910,917,918,958,959,967,970,974,975,989,990,991,1014,1015});
  vector<int> v40({774,778,781,785,869,872,877,878,881,885,927,928,992,993,999,1000,1061,1062,1082});
  vector<int> v41({672,673,719,720,721,722,723,768,770,771});
  vector<int> v42({674,675,727,729,741,742});
  vector<int> v43({738,739,740,792,793,797});
  vector<int> v44({689,690,691,748,749});
  vector<int> v45({602,603,653,655,656});
  vector<int> v46({607,608,658,696,697});
  vector<int> v47({146,147,185,188,298,299,302,392,393,399,401,402,409,466,467,468,511,604});

vector <int> B_T({1200,4371,7001,8085,4237,39416,24752,15300,19996,24206,15580,0,13017,0,0,0,4793,24266,25120,11476,14008,12039,2720,8317,2592,0,309,2125,836,4030,11256,9586,196,1809,1274,3024,8910,10291,9508,22172,3032,2492,1676,3180,4156,7903,11145}); //vector con areas de Brodman totales del Talairach

vector<int> B(47, 0);// vector que guarda areas de Brodmann absolutas
vector<int> B_R(47, 0); //vector que guarda areas de Brodmann relativizadas

VECTOR3D cerebro = CREATEM3D(ii,jj,kk);
VECTOR3D talairach = CREATEM3D(ii,jj,kk);

VECTOR3D p = CREATEM3D(ii,jj,kk);
VECTOR3D D = CREATEM3D(ii,jj,kk);
VECTOR3D C = CREATEM3D(ii,jj,kk);
VECTOR3D P_optimizado = CREATEM3D(ii,jj,kk);
VECTOR3D M_optimizado = CREATEM3D(ii,jj,kk);
VECTOR3D P = CREATEM3D(ii,jj,kk);
VECTOR3D M = CREATEM3D(ii,jj,kk);
VECTOR3D C_k2 = CREATEM3D(ii,jj,kk);
VECTOR3D C_k1 = CREATEM3D(ii,jj,kk);
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

void buscar_areas_Broodman(int i,int j,int k,vector<int>& B){ //busco areas de Brodmann
	
if(pertenece(v1,G3D(talairach,i,j,k)))                               B[1]=B[1]+1;
if(pertenece(v2,G3D(talairach,i,j,k)))                               B[2]=B[2]+1;
if(pertenece(v3,G3D(talairach,i,j,k)))                               B[3]=B[3]+1;
if(pertenece(v4,G3D(talairach,i,j,k)))                               B[4]=B[4]+1;
if(pertenece(v5,G3D(talairach,i,j,k)))                               B[5]=B[5]+1;
if(pertenece(v6,G3D(talairach,i,j,k)))                               B[6]=B[6]+1;
if(pertenece(v7,G3D(talairach,i,j,k)))                               B[7]=B[7]+1;
if(pertenece(v8,G3D(talairach,i,j,k)))                               B[8]=B[8]+1;
if(pertenece(v9,G3D(talairach,i,j,k)))                               B[9]=B[9]+1;
if(pertenece(v10,G3D(talairach,i,j,k)))                              B[10]=B[10]+1;
if(pertenece(v11,G3D(talairach,i,j,k)))                              B[11]=B[11]+1;
if(pertenece(v13,G3D(talairach,i,j,k)))                              B[13]=B[13]+1;
if(pertenece(v17,G3D(talairach,i,j,k)))                              B[17]=B[17]+1;
if(pertenece(v18,G3D(talairach,i,j,k)))                              B[18]=B[18]+1;
if(pertenece(v19,G3D(talairach,i,j,k)))                              B[19]=B[19]+1;
if(pertenece(v20,G3D(talairach,i,j,k)))                              B[20]=B[20]+1;
if(pertenece(v21,G3D(talairach,i,j,k)))                              B[21]=B[21]+1;
if(pertenece(v22,G3D(talairach,i,j,k)))                              B[22]=B[22]+1;
if(pertenece(v23,G3D(talairach,i,j,k)))                              B[23]=B[23]+1;
if(pertenece(v24,G3D(talairach,i,j,k)))                              B[24]=B[24]+1;
if(pertenece(v25,G3D(talairach,i,j,k)))                              B[25]=B[25]+1;
if(pertenece(v27,G3D(talairach,i,j,k)))                              B[27]=B[27]+1;
if(pertenece(v28,G3D(talairach,i,j,k)))                              B[28]=B[28]+1;
if(pertenece(v29,G3D(talairach,i,j,k)))                              B[29]=B[29]+1;
if(pertenece(v30,G3D(talairach,i,j,k)))                              B[30]=B[30]+1;
if(pertenece(v31,G3D(talairach,i,j,k)))                              B[31]=B[31]+1;
if(pertenece(v32,G3D(talairach,i,j,k)))                              B[32]=B[32]+1;
if(pertenece(v33,G3D(talairach,i,j,k)))                              B[33]=B[33]+1;
if(pertenece(v34,G3D(talairach,i,j,k)))                              B[34]=B[34]+1;
if(pertenece(v35,G3D(talairach,i,j,k)))                              B[35]=B[35]+1;
if(pertenece(v36,G3D(talairach,i,j,k)))                              B[36]=B[36]+1;
if(pertenece(v37,G3D(talairach,i,j,k)))                              B[37]=B[37]+1;
if(pertenece(v38,G3D(talairach,i,j,k)))                              B[38]=B[38]+1;
if(pertenece(v39,G3D(talairach,i,j,k)))                              B[39]=B[39]+1;
if(pertenece(v40,G3D(talairach,i,j,k)))                              B[40]=B[40]+1;
if(pertenece(v41,G3D(talairach,i,j,k)))                              B[41]=B[41]+1;
if(pertenece(v42,G3D(talairach,i,j,k)))                              B[42]=B[42]+1;
if(pertenece(v43,G3D(talairach,i,j,k)))                              B[43]=B[43]+1;
if(pertenece(v44,G3D(talairach,i,j,k)))                              B[44]=B[44]+1;
if(pertenece(v45,G3D(talairach,i,j,k)))                              B[45]=B[45]+1;
if(pertenece(v46,G3D(talairach,i,j,k)))                              B[46]=B[46]+1;
if(pertenece(v47,G3D(talairach,i,j,k)))                              B[47]=B[47]+1;
							
	
	
	}

void grabar_matriz(VECTOR3D &mat){
//imprime en datos en algun orden la matriz
	
for(int i=0;i<ii;i++){
	for(int j=0;j<jj;j++){
		for(int k=0;k<kk;k++){
			if(G3D(mat,i,j,k)!=0)
			datos<<G3D(mat,i,j,k)<<endl;
}
}
}

}


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
	cout<<"dia: "<< dia <<endl;
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

void iteracion_de_convergencia(int n,int cantidad1,int cantidad2,int cantidad3,vector<int>& B,vector<int>& B_R, vector<int>& B_T,VECTOR3D &C,VECTOR3D &C_k1,VECTOR3D &C_k2){
	double max_error=0;
	if (dia <= 1500){
	max_error = 1;
	}else{
	max_error = 10;
	}
	int iter=0;
	double error = 1000;


	copyMatrix(C_k2,C);
	while((iter < max_iter) && (error > max_error)){
		copyMatrix(C_k1,C_k2);
		for (int k=1;k<chunkSize;k++){
			for (int j=1;j<jj-1;j++){
				for (int i=1;i<ii-1;i++){

					//S3D(M,i,j,k,0);
					S3D(M,i,j,k, G3D(M_optimizado,i,j,k) * (G3D(C_k1,i+1,j,k)+G3D(C_k1,i-1,j,k)+G3D(C_k1,i,j+1,k)+G3D(C_k1,i,j-1,k)+G3D(C_k1,i,j,k+1)+G3D(C_k1,i,j,k-1)-6*G3D(C_k1,i,j,k)));
						

					S3D(P,i,j,k, G3D(P_optimizado,i,j,k) * G3D(C_k1,i,j,k) * (1 - G3D(C_k1,i,j,k) / C_max));
					S3D(C_k2,i,j,k, G3D(C,i,j,k) + G3D(P,i,j,k) + G3D(M,i,j,k));
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

	}
    // Actualizar malla
	copyMatrix(C,C_k2);
	if( (world_size<=2 && world_rank==0) || ko>=chunkSize*world_rank && ko<=chunkSize*(world_rank+1)-1 ){
		//error = restaMax(C,CK2_slice);
		guardar_datos(n,error,cantidad1,cantidad2,cantidad3,B,B_R, B_T);
	
	}

}

#endif
