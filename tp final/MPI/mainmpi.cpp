#include "variablesGlobalesmpi.h"
//DEJO ITERACION TEMPORAL ACA PARA NO HACER TANTO QUILOMBO EN EL .h
void iteracion_temporal(VECTOR3D &C_slice,VECTOR3D &CK1_slice,VECTOR3D &CK2_slice){

	for(int n=0;n<nn;n++){
		//cout << "n: "<<n<<endl;
		//Calculo del volumen tumoral y chequeo de areas de Brodmann

			//if(n%200==0){
			//	dumpMatrixToVtk(C, "tumor_" + to_string(n));   
			//}

			/*if ((G3D(C,io,jo,ko) < C_mig) ){
			for(int k=0;k<kk;k++){
				for(int j=0;j<jj;j++){
					for(int i=0;i<ii;i++){
						if (G3D(C_slice,i,j,k) >= 1){
						}
						if(G3D(C_slice,i,j,k) >= 1e7){
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
			}*/	
			
        //screenMessage("Waiting for processes to border exchange...\n");
        MPI_Barrier(MPI_COMM_WORLD);
        // Intercambio de bordes Cn entre procesos vecinos
        //screenMessage("Starting border exchange...\n");
        	
        doHaloTransfer(C_slice, world_rank, world_size, workingSize);
        
        MPI_Barrier(MPI_COMM_WORLD);
	iteracion_de_convergencia(n,cantidad1,cantidad2,cantidad3,B,B_R,B_T,C_slice,CK1_slice,CK2_slice);

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
 // Initialize the MPI environment
     MPI_Init(NULL, NULL);
      // Find out rank, size
      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 //if (world_rank == 0) {
	
	ReadDifussionData("./Cerebro.csv", 0, 0, 0, ii-1, jj-1, kk-1, cerebro);//lee del archivo a matriz
	ReadDifussionData("./Talaraich.csv", 0, 0, 0, ii-1, jj-1, kk-1, talairach);//lee del archivo a matriz
	printf ("Preprocessing difusion Matrix\n");
	//dumpMatrixToVtk(cerebro, "cerebro difusion");
	printf ("Preprocessing talairach Matrix\n");
      //  dumpMatrixToVtk(talairach, "talairach");


	printf ("Difusion\n");
	TransformDifusion();//inicializa valores de la matriz
        //dumpMatrixToVtk(D, "matriz D");
	printf ("Preprocessing initial brain Matrix\n");
	inicializarCondiciones();
   	info.open("info.txt");
	datos.open("datos.txt");
	aux.open("aux.txt");

//}
    MPI_Barrier(MPI_COMM_WORLD);

    chunkSize = (kk / world_size) ;
    scatterSize = chunkSize - 1;
    workingSize = chunkSize + 1;

    // Capas a procesar (el contenido )
    VECTOR3D CK2_slice        = CREATEM3D(ii, jj, workingSize);
    VECTOR3D CK1_slice        = CREATEM3D(ii, jj, workingSize);
    VECTOR3D C_slice        = CREATEM3D(ii, jj, workingSize);
    screenMessage("Procs: %u, chunkSize: %u, X_SIZE: %u, Y_SIZE: %u, Z_SIZE: %u\n", world_size, chunkSize, ii, jj, kk);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Scatter(&C[0],chunkSize* (ii*jj),MPI_DOUBLE,&C_slice[0],chunkSize * (ii*jj),MPI_DOUBLE,0,MPI_COMM_WORLD);

 // Paso capas internas a los procesos (arranco de Z1)
    
    screenMessage("Scatering CK2\n");
    MPI_Scatter(&C_k1[0],chunkSize* (ii*jj),MPI_DOUBLE,&CK1_slice[0],chunkSize * (ii*jj),MPI_DOUBLE,0,MPI_COMM_WORLD);
    
	screenMessage("Scatering CK1\n");
	MPI_Scatter(&C_k2[0],chunkSize* (ii*jj),MPI_DOUBLE,&CK2_slice[0],chunkSize* (ii*jj) ,MPI_DOUBLE,0,MPI_COMM_WORLD);
 
/*	if(world_size>2){
	S3D(C_slice,io,jo,ko,1);

	S3D(C_slice,io,jo,ko-chunkSize-1,1);
	}
*/
	MPI_Barrier(MPI_COMM_WORLD);
   
   //cout<<G3D(C_slice,io,jo,ko-chunkSize-1)<<"rank " << world_rank <<endl;
    	TransformDifusion();//inicializa valores de la matriz
        //dumpMatrixToVtk(D, "matriz D");
	printf ("Preprocessing initial brain Matrix\n");
	inicializarCondiciones();


//A PARTIR DE ACA SE CUENTA EL TIEMPO, ESTO SE PUEDE MODIFICAR PARA QUE SE CUENTE EN OTRO LUGAR
	double t1,t2,elapsed;
    	struct timeval tp;
   	int rtn;
	if (world_rank == 0){
  	info << "Simulacion del paciente \n" ;

	

 	
        	rtn=gettimeofday(&tp, NULL);
        	t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;    
   	}
  


	printf ("Ejecutando iteracion temporal\n");
	
	iteracion_temporal(C_slice,CK1_slice,CK2_slice);

	info.close();
	datos.close();

	MPI_Barrier(MPI_COMM_WORLD);

	screenMessage("End of loop\n");
    
	MPI_Gather(&C_slice[ii*jj],scatterSize * (ii*jj),MPI_DOUBLE,&C[ii*jj],scatterSize* (ii*jj),MPI_DOUBLE,0,MPI_COMM_WORLD);
    
	screenMessage("End of data gathering\n");

        if (world_rank == 0){
		//fin toma tiempos
		rtn=gettimeofday(&tp, NULL);
		t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
		elapsed=t2-t1;
		
		printf("Tiempo empleado: %g\n",elapsed);

		printf ("Preprocessing result Matrix\n");
		dumpMatrixToVtk(C, "tumor_out");   
		cout<< "PROCESO TERMINADO"<<endl;
        }   

	MPI_Finalize();   
	return 0;
	}
