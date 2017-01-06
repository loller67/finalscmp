#include "variablesGlobales.h"
//DEJO ITERACION TEMPORAL ACA PARA NO HACER TANTO QUILOMBO EN EL .h
void iteracion_temporal(){

	for(int n=0;n<nn;n++){
		//cout << "n: "<<n<<endl;
		//Calculo del volumen tumoral y chequeo de areas de Brodmann
/*		if(n%200==0){
			dumpMatrixToVtk(C, "tumor_" + to_string(n));   
		}
*/		for(int k=0;k<kk;k++){
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
	
	

	ReadDifussionData("./Cerebro.csv", 0, 0, 0, ii-1, jj-1, kk-1, cerebro);//lee del archivo a matriz
	ReadDifussionData("./Talaraich.csv", 0, 0, 0, ii-1, jj-1, kk-1, talairach);//lee del archivo a matriz
	printf ("Preprocessing difusion Matrix\n");
	dumpMatrixToVtk(cerebro, "cerebro difusion");
	printf ("Preprocessing talairach Matrix\n");
        dumpMatrixToVtk(talairach, "talairach");




    	printf ("Difusion\n");
	TransformDifusion();//inicializa valores de la matriz
	info.open("info.txt");
	datos.open("datos.txt");
	printf ("Preprocessing initial brain Matrix\n");
	inicializarCondiciones();

	
	
  	info << "Simulacion del paciente \n" ;



	//A PARTIR DE ACA SE CUENTA EL TIEMPO, ESTO SE PUEDE MODIFICAR PARA QUE SE CUENTE EN OTRO LUGAR
	double t1,t2,elapsed;
    	struct timeval tp;
   	int rtn;

    	rtn=gettimeofday(&tp, NULL);
    	t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;  
  
	printf ("Ejecutando iteracion temporal\n");
	iteracion_temporal();

	info.close();
	datos.close();
	//fin toma tiempos
	rtn=gettimeofday(&tp, NULL);
        t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        elapsed=t2-t1;
        printf("Tiempo empleado: %g\n",elapsed);

	printf ("Preprocessing result Matrix\n");
	dumpMatrixToVtk(C, "cerebro_out");   
	cout<< "PROCESO TERMINADO"<<endl;
	return 0;
	}
