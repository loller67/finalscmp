#include "variablesGlobales.h"
//DEJO ITERACION TEMPORAL ACA PARA NO HACER TANTO QUILOMBO EN EL .h
void iteracion_temporal(){
	
	for(int n=0;n<nn;n++){
		
		//Calculo del volumen tumoral y chequeo de areas de Brodmann
		int cantidad1 = 0; //para deteccion de celulas tumorales
		int cantidad2 = 0; //para diagnostico
		int cantidad3 = 0; //para letalidad
		vector <int> B_T({1200,4371,7001,8085,4237,39416,24752,15300,19996,24206,15580,0,13017,0,0,0,4793,24266,25120,11476,14008,12039,2720,8317,2592,0,309,2125,836,4030,11256,9586,196,1809,1274,3024,8910,10291,9508,22172,3032,2492,1676,3180,4156,7903,11145}); //vector con areas de Brodman totales del Talairach
		vector<int> B(47, 0);// vector que guarda areas de Brodmann absolutas
		vector<int> B_R(47, 0); //vector que guarda areas de Brodmann relativizadas
		for(int k=0;k<kk;k++){
			for(int j=0;j<jj;j++){
				for(int i=0;i<ii;i++){
					if (G3D(C,i,j,k) >= 1){
						cantidad1++;
					}
					if(G3D(C,i,j,k) >= 1e7){
						cantidad2=cantidad2+1;
						if ((i>80)&&(i<102)&&(j>90)&&(j<117)&&(k>61)&&(k<71)){ //area del foramen del tentorio
							cantidad3=cantidad3+2;
						}else{
							if(pertenece(cerebelo,G3D(talairach,i,j,k)))
									cantidad3=cantidad3+2;
							if(pertenece(tallo,G3D(talairach,i,j,k)))//tallo cerebral, medula, protuberancia, mesensefalo
									cantidad3=cantidad3+2;
							if(!pertenece(cerebelo,G3D(talairach,i,j,k)&& !pertenece(tallo,G3D(talairach,i,j,k))))
									cantidad3=cantidad3+1;
							}
						}
					
					 if (G3D(D,i,j,k)==0.051){ //si estoy en corteza (de cerebro o cerebelo)
							buscar_areas_Broodman(i,j,k,B);
						}
	   
				} 
			}  
		}   
		
		
		iteracion_de_convergencia(n,cantidad1,cantidad2,cantidad3,B,B_R,B_T);

	}
	

}


int main(){
	
	


    ReadDifussionData("Cerebro.csv", ii, jj, kk, 0, 0, 0, D);//lee del archivo a matriz
    ReadDifussionData("Talairach.csv", ii, jj, kk, 0, 0, 0, talairach);//lee del archivo a matriz

    
	
	TransformDifusion();//inicializa valores de la matriz

	inicializarCondiciones();
	/*
	info = fopen('Info_1i.txt','a+');
	fprintf(info,'Simulacion i del paciente 1 \n');
	status = fclose(info);
	*/

	iteracion_temporal();
	cout<< "PROCESO TERMINADO"<<endl;

	return 0;
	}
