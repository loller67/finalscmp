#include "variablesGlobalesopenmp.h"
//DEJO ITERACION TEMPORAL ACA PARA NO HACER TANTO QUILOMBO EN EL .h
void iteracion_temporal(){
	int i,j,k;
	int plus_i,plus_j,plus_k,minus_i,minus_j,minus_k,m,p,m_ant,c,p_nuevo,c_ant;
	for(int n=0;n<nn;n++){
		if(n %100==0){
			cout<<n<<endl;
		}
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
			S3D(C_k1,io,jo,ko,G3D(C_k2,io,jo,ko));
		       //Calcular dominio
			if (G3D(C,io,jo,ko) < C_mig){ //proliferacion
				S3D(M,io,jo,ko,0);
				S3D(P,io,jo,ko, G3D(P_optimizado,io,jo,ko) * G3D(C_k1,io,jo,ko) * (1 - G3D(C_k1,io,jo,ko) / C_max));
				S3D(C_k2,io,jo,ko, G3D(C,io,jo,ko) + G3D(P,io,jo,ko) + G3D(M,io,jo,ko));
				if (G3D(C_k2,io,jo,ko) > C_max){
					S3D(C_k2,io,jo,ko, C_max);
				}
				if (G3D(C_k2,io,jo,ko) < 0.00001){
					S3D(C_k2,io,jo,ko,0);
				}
				if(G3D(C_k1,io,jo,ko)>=G3D(C_k2,io,jo,ko)){// calculo valor absoluto del error
					error=G3D(C_k1,io,jo,ko)-G3D(C_k2,io,jo,ko);

				}else{
					error=G3D(C_k2,io,jo,ko)-G3D(C_k1,io,jo,ko);
				}
		       		iter++;
			}else{//proliferacion y migracion

			       copyMatrix(C_k1,C_k2);
		#pragma omp parallel for collapse(3)  schedule(static) private(plus_i,plus_j,plus_k,minus_i,minus_j,minus_k,m,p,m_ant,c,p_nuevo) num_threads(threads)  
			       for (k=1;k<kk-1;k++){
				   for (j=1;j<jj-1;j++){
				       for (i=1;i<ii-1;i++){

						m_ant=G3D(M_optimizado,i,j,k);
						plus_i= G3D(C_k1,i+1,j,k);
						minus_i=G3D(C_k1,i-1,j,k);
						plus_j=G3D(C_k1,i,j+1,k);
						minus_j=G3D(C_k1,i,j-1,k);
						plus_k=G3D(C_k1,i,j,k+1);
						minus_k=G3D(C_k1,i,j,k-1);
						c=G3D(C_k1,i,j,k);
						p=G3D(P_optimizado,i,j,k);
						S3D(M,i,j,k, m_ant* (plus_i+minus_i+plus_j+minus_j+plus_k+minus_k-6*c));
						m =G3D(M,i,j,k);
						c_ant= G3D(C,i,j,k);
						S3D(P,i,j,k, p * c * (1 - c / C_max));
						p_nuevo= G3D(P,i,j,k);
								
								S3D(C_k2,i,j,k, c_ant + G3D(P,i,j,k) + m);
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
			}
		    }
		    // Actualizar malla
		    
		    copyMatrix(C,C_k2);



			}
	

}


void blurMatrix3d(VECTOR3D& m, int i, int j, int k){
    
    for (int x = 1; x < i; x++){
        for (int y = 1; y < j; y++){
            for (int z = 1; z < k; z++){
                double val = (G3D(m,x,y,z)+G3D(m,x+1,y,z)+G3D(m,x-1,y,z)+G3D(m,x,y+1,z)+G3D(m,x,y-1,z)+G3D(m,x,y,z+1)+G3D(m,x,y,z-1))/8.0;
                S3D( m, x,  y,  z,  val);
            }
        }
    }
    
}


int main(){
	
	ReadDifussionData("./Cerebro.csv", 0, 0, 0, ii-1, jj-1, kk-1, cerebro);//lee del archivo a matriz
	ReadDifussionData("./Talaraich.csv", 0, 0, 0, ii-1, jj-1, kk-1, talairach);//lee del archivo a matriz
	printf ("Difusion\n");
	TransformDifusion();//inicializa valores de la matriz
	printf ("Preprocessing initial brain Matrix\n");
	inicializarCondiciones();
	printf ("Ejecutando iteracion temporal\n");
	time_t sec1,sec2;
	time(&sec1); // tiempo en segundos
	iteracion_temporal();
	// proceso cuya duracion quieres medir
	time(&sec2);
	cout << sec2-sec1 << endl;
	cout<< "PROCESO TERMINADO"<<endl;

	return 0;
	}
