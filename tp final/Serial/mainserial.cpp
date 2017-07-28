#include "variablesGlobales.h"
//DEJO ITERACION TEMPORAL ACA PARA NO HACER TANTO QUILOMBO EN EL .h
void iteracion_temporal(){

	for(int n=0;n<nn;n++){

		//descomentar si quieren vtk del tumor temporal
		//if(n%500==0){
		//	dumpMatrixToVtk(C, "tumor_" + to_string(n));   
		//}
		//Calculo del volumen tumoral y chequeo de areas de Brodmann
	if (!(G3D(C,io,jo,ko) < C_mig)){//si no estoy en migracion
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
						
					
					 	if (G3D(D,i,j,k)==0.051){ //si estoy en corteza (de cerebro o cerebelo)
							buscar_areas_Broodman(i,j,k,B);
						}
	   				}
				} 
			}  
		}   
	}
		iteracion_de_convergencia(n,cantidad1,cantidad2,cantidad3,B,B_R,B_T);

	}
	

}



void blur(VECTOR3D& m, int i, int j, int k){
    
    for (int x = 1; x < i; x++){
        for (int y = 1; y < j; y++){
            for (int z = 1; z < k; z++){
                double val = (G3D(m,x,y,z)+G3D(m,x+1,y,z)+G3D(m,x-1,y,z)+G3D(m,x,y+1,z)+G3D(m,x,y-1,z)+G3D(m,x,y,z+1)+G3D(m,x,y,z-1))/8.0;
                S3D( m, x,  y,  z,  val);
            }
        }
    }
    
}

void extender (VECTOR3D& m,VECTOR3D& res, int n){//precondicion, res tiene que estar cargado con todos ceros y ya con la dimension dada y n debe ser el crecimiento

// dentro de cada matriz al hacer el crecimiento siempre se crea n/2 filas y una n/2 columnas de CEROS para luego correr el blur 

	for(int i=0;i<ii;i++){
		for(int j=0; j<jj;j++){
			for(int k=0; k<kk;k++){
			
			//aca falta armar la logica del "crecimiento de la matriz", deberia quedar algo asi
			/*
			1 n vacios/2 2 n vacios/2 3 
			n vacios/2 
			4 n vacios/2 5 n vacios/2 6
			n vacios/2 
			7 n vacios/2 8 n vacios/2 9 

			*/
			S3D(res,n/2*i,n/2*j,n/2*k,G3D(m,i,j,k));


			}


		}

	}

}


int main(){
	
	

	ReadDifussionData("./Cerebro.csv", 0, 0, 0, 181-1, 217-1, 181-1, cerebro);//lee del archivo a matriz
	ReadDifussionData("./Talaraich.csv", 0, 0, 0, 181-1, 217-1, 181-1, talairach);//lee del archivo a matriz
	printf ("Preprocessing difusion Matrix\n");
	//dumpMatrixToVtk(cerebro, "cerebro difusion");
	printf ("Preprocessing talairach Matrix\n");
        //dumpMatrixToVtk(talairach, "talairach");
	extender(talairach_aux,talairach,2);
        extender(cerebro_aux,cerebro,2);
        blur(cerebro,ii-1,jj-1,kk-1);
        blur(talairach, ii-1,jj-1,kk-1);




    	printf ("Difusion\n");
	TransformDifusion();//inicializa valores de la matriz
        //dumpMatrixToVtk(D, "matriz D");
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
	dumpMatrixToVtk(C, "tumor_out");   
	cout<< "PROCESO TERMINADO"<<endl;
	return 0;
	}
