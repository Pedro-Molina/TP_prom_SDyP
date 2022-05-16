#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>

int T, N;
float *V, *V2;
int *converge;
int *iteraciones;

void swap( float ** x, float ** y){
	float *temp = *x;
    *x = *y;
    *y = temp;
}

double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

void *function(void *arg) {
	int tid = *(int *) arg;

	int block = N/T;
	int begin = block*tid;
	int end = begin+block;

	float *Vlocal, *V2local;
	Vlocal = V;
	V2local = V2;

	while (!converge[tid]) {

		if (tid == 0) {
			begin++;
			V2[0]= (Vlocal[0] + Vlocal[1])/2;
		} else if (tid == T-1) {
			end--;
			V2[N-1]= (Vlocal[N-1] + Vlocal[N-2])/2;
		}

		//reduccion
		for (int i = begin; i < end; i++){
			V2[i] = (Vlocal[i-1] + Vlocal[i] + Vlocal[i+1])/3;
		}

		//chequeo de convergencia
		int i = begin;
		float aux = V2[0];

		converge[tid] = 1;
		while ((i < end) && (converge)){ 
			if (fabs(aux-V2[i]) > 0.01){//si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
				converge[tid] = 0;
			}
			i++;
		}

		iteraciones[tid]++;
		swap(&V, &V2);

	}

	pthread_exit(NULL);
}

int main(int argc, const char *argv[])
{

	if (argc < 3) {
		printf("Faltan argumentos. programa N T\n");
		return 1;
	}

	N = argv[1];
	T = argv[2];

	//Aloca memoria para las matrices
	V=(float*)malloc(sizeof(float)*N);
	V2=(float*)malloc(sizeof(float)*N);

	converge = (int*) malloc(sizeof(int)*T);
	iteraciones = (int*) malloc(sizeof(int)*T);
	for (int i = 0; i < T; i++) {
		converge[i] = 0;
		iteraciones[i] = 0;
	}

	//inicializacion del arreglo
	for (int i = 0; i < N; i++){
		V[i] = (float)rand()/(float)(RAND_MAX); //funciona en MPI?
		printf (" V[%d] = %f \n",i, V[i]); 
	}

	pthread_t myThreads[T];
	int thread_ids[T];

	double timetick = dwalltime();

	for (int id = 0; id < T; id++) {
		thread_ids[id] = id;
		pthread_create(&myThreads[id], NULL, &function, (void*) &thread_ids[id]);
	}



	for (int id = 0; id < T; id++) {
		pthread_join(&myThreads[id], NULL);
	}

	printf("Tiempo en segundos %f\n", dwalltime() - timetick);

	return 0;
}