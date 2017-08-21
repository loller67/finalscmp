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
	   				}
				} 
			}  
		}   
	}
		iteracion_de_convergencia(n,cantidad1,cantidad2,cantidad3);

	}
	

}



void blur(VECTOR3D& m, int i, int j, int k){
  

for(int n= 0;n<1;n++ ){
		for (int x = 1; x < i; x++){
		        for (int y = 1; y < j; y++){
		                    for (int z = 1; z < k; z++){
		                                    double val = (G3D(m,x,y,z)+G3D(m,x+1,y,z)+G3D(m,x-1,y,z)+G3D(m,x,y+1,z)+G3D(m,x,y-1,z)+G3D(m,x,y,z+1)+G3D(m,x,y,z-1))/8.0;
		                                                    S3D( m, x,  y,  z,  val);
		                            }
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
	
	
	ReadDifussionData("./Cerebro.csv", 0, 0, 0, 181-1, 217-1, 181-1, cerebro_aux);//lee del archivo a matriz
	ReadDifussionData("./Talaraich.csv", 0, 0, 0, 181-1, 217-1, 181-1, talairach_aux);//lee del archivo a matriz
	//extender(talairach_aux,talairach,2);
	//extender(cerebro_aux,cerebro,2);
	//printf ("Blur cerebro\n");
	//blur(cerebro,ii-1,jj-1,kk-1);
	//printf ("Blur Tala\n");
	//blur(talairach, ii-1,jj-1,kk-1);
	//dumpMatrixToVtk(cerebro, "cerebro difusion");
	//dumpMatrixToVtk(talairach, "el taladro papaaa");
   	printf ("Difusion\n");
	TransformDifusion();//inicializa valores de la matriz
	printf ("Preprocessing initial brain Matrix\n");
	inicializarCondiciones();
	printf ("Ejecutando iteracion temporal\n");
	time_t sec1,sec2;
	time(&sec1); // tiempo en segundos
	// para el momento de aplicar extender y blur, tambien seteamos extendi en true para no hacer lio con el ii jj kk
	iteracion_temporal();
	// proceso cuya duracion quieres medir
	time(&sec2);
	cout << sec2-sec1 << endl;
	cout << "VERSION SERIAL TAMAÑO " << ii <<jj <<kk<<endl; 
	cout<< "PROCESO TERMINADO"<<endl;
	return 0;
	}
