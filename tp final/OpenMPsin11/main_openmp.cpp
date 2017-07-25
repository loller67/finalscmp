#include "variablesGlobalesopenmp.h"
//DEJO ITERACION TEMPORAL ACA PARA NO HACER TANTO QUILOMBO EN EL .h
void iteracion_temporal(){

	if(extendi == true){
		ii = iii;
		jj=jjj;
		kk= kkk;

	}



	for(int n=0;n<nn;n++){
		if(n % 100 ==0){

		cout<<n<<endl;

		}
		double max_error=0;
		int i,j,k;
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
		#pragma omp parallel for collapse(3)  schedule(static) num_threads(threads)  
			       for (int k=1;k<kk-1;k++){
				   for (int j=1;j<jj-1;j++){
				       for (int i=1;i<ii-1;i++){

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
		  
		       //Calcular error y actualizar
		       error = restaMax(C_k1,C_k2);
		       iter++;
			}
		    }
		    // Actualizar malla
		    
		    copyMatrix(C,C_k2);



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
	
	ReadDifussionData("./Cerebro.csv", 0, 0, 0, ii-1, jj-1, kk-1, cerebro);//lee del archivo a matriz
	ReadDifussionData("./Talaraich.csv", 0, 0, 0, ii-1, jj-1, kk-1, talairach);//lee del archivo a matriz
   	printf ("Difusion\n");
	TransformDifusion();//inicializa valores de la matriz
	printf ("Preprocessing initial brain Matrix\n");
	inicializarCondiciones();

	// para el momento de aplicar extender y blur, tambien seteamos extendi en true para no hacer lio con el ii jj kk
	extender(C,C_EXT,2);
	blur(C_EXT,iii-1,jjj-1,kkk-2);

	printf ("Ejecutando iteracion temporal\n");
	iteracion_temporal();
	cout<< "PROCESO TERMINADO"<<endl;
	return 0;
	}
