#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void swap( float ** x, float ** y){

		float *temp = *x;
    *x = *y;
    *y = temp;
}

int main(int argc,char*argv[]){
	int N;
	float *V, *V2;
	int i, iteraciones = 0;
	int converge =0;
	float aux;

 //Controla los argumentos al programa
	if ((argc != 2) || ((N = atoi(argv[1])) <= 0) )
  {
    printf("\nUsar: %s n\n  n: Dimension de la matriz (nxn X nxn)\n", argv[0]);
    exit(1);
  }

 //Aloca memoria para las matrices
  V=(float*)malloc(sizeof(float)*N);
	V2=(float*)malloc(sizeof(float)*N);
	
	//inicializacion del arreglo
	for (i = 0; i < N; i++){
		V[i] = (float)rand()/(float)(RAND_MAX); //funciona en MPI?
		printf (" V[%d] = %f \n",i, V[i]); 
	}
	//loop principal
	while (!converge){
	
		//reduccion
		for (i = 1; i < N - 1; i++){
			V2[i] = (V[i-1] + V[i] + V[i+1])/3;
		}
		
		V2[0]= (V[0] + V[1])/2;
		V2[N-1]= (V[N-1] + V[N-2])/2;
		
		converge = 1;
		i = 1;
		aux = V2[0];
		//chequeo de convergencia
		while ((i < N) && (converge)){ 
			if ((fabs(aux - V2[i])) > 0.01){//si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
					//printf ("%f - " , fabs(V2[0]-V2[i]));
					converge = 0;
			}
			i++;
		}
		
		swap(&V, &V2);// no es necesario, puedo asignar simplemente v2 a v, porque v2 se sobreescribe
		iteraciones++;
		printf ("iteracion = %d\n", iteraciones);
	}
	
	
	for (i = 0; i < N; i++){
		printf (" V2[%d] = %f \n",i, V2[i]); 
	}
	printf ("iteraciones = %d", iteraciones);
}
