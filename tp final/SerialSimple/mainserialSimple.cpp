#include "variablesGlobalesSimple.h"
//DEJO ITERACION TEMPORAL ACA PARA NO HACER TANTO QUILOMBO EN EL .h
void iteracion_temporal(){

	for(int n=0;n<nn;n++){

	double max_error=0;
	int i,j,k;
	if (dia <= 1500){
        max_error = 1;
	}else{
        max_error = 10;
    }
    int iter=0;
    double error = 1000;
    double p_aux;
    double d_aux;
    copyMatrix(C_k2,C);
               copyMatrix(C_k1,C_k2);
	       for (int k=1;k<kk-1;k++){
		   for (int j=1;j<jj-1;j++){
		       for (int i=1;i<ii-1;i++){
				d_aux=d*(G3D(C_k1,i+1,j,k)+G3D(C_k1,i-1,j,k)+G3D(C_k1,i,j+1,k)+G3D(C_k1,i,j-1,k)+G3D(C_k1,i,j,k+1)+G3D(C_k1,i,j,k-1)-6*G3D(C_k1,i,j,k));
				p_aux=p* G3D(C_k1,i,j,k) * (1 - G3D(C_k1,i,j,k) / C_max);		


				S3D(C_k2,i,j,k, G3D(C,i,j,k) +p_aux+d_aux);
				if (G3D(C_k2,i,j,k) > C_max){
					S3D(C_k2,i,j,k, C_max);
				}
				if (G3D(C_k2,i,j,k) < 0.00001){
					S3D(C_k2,i,j,k,0);
				}
			}
		   }
	       }
  
       //Calcular error y actualizar
       error = restaMax(C_k1,C_k2);
       iter++;
	
   // }
    // Actualizar malla

    copyMatrix(C,C_k2);
    guardar_datos(n,error);

	}
	

}


int main(){
	
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
