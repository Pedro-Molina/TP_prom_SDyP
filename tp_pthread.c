#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>

int T, N;
double *V, *V2;
int *converge;
int *iteraciones;

pthread_barrier_t barrera;

void swap( double ** x, double ** y){
	double *temp = *x;
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

	double *Vlocal, *V2local;
	Vlocal = V;
	V2local = V2;

	printf("Thread %d arranca. begin = %d, end = %d\n", tid, begin, end);

	while (!converge[tid]) {

		/*
		if (tid == 0) {
			begin++;
			V2local[0] = (Vlocal[0] + Vlocal[1])/2;
		} else if (tid == T-1) {
			end--;
			V2local[N-1] = (Vlocal[N-1] + Vlocal[N-2])/2;
		}
		*/

		//reduccion
		for (int i = begin+1; i < end-1; i++){
			V2local[i] = (Vlocal[i-1] + Vlocal[i] + Vlocal[i+1])/3;
		}
		V2local[begin] = (Vlocal[begin] + Vlocal[begin+1]) / 2;
		V2local[end-1] = (Vlocal[end-1] + Vlocal[end-2]) / 2;

		//chequeo de convergencia
		int i = begin;
		double aux = V2[0];

		converge[tid] = 1;
		while ((i < end) && (converge[tid])){ 
			if (fabs(aux - V2local[i]) > 0.05){//si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
				converge[tid] = 0;
				//printf("la resta me dio %f\n", fabs(aux - V2local[i]));
			}
			i++;
		}

		iteraciones[tid]++;

		swap(&Vlocal, &V2local);

		//if (iteraciones[tid] % 5000 == 0) printf("Thread %d termina iteracion %d\n", tid, iteraciones[tid]);
	}
	printf("THREAD %d TERMINA\n", tid);

	pthread_barrier_wait(&barrera);

	pthread_exit(NULL);
}

int main(int argc, const char *argv[])
{

	if (argc < 3) {
		printf("Faltan argumentos. programa N T\n");
		return 1;
	}

	N = atoi(argv[1]);
	T = atoi(argv[2]);

	//Aloca memoria para las matrices
	V=(double*)malloc(sizeof(double)*N);
	V2=(double*)malloc(sizeof(double)*N);

	converge = (int*) malloc(sizeof(int)*T);
	iteraciones = (int*) malloc(sizeof(int)*T);
	for (int i = 0; i < T; i++) {
		converge[i] = 0;
		iteraciones[i] = 0;
	}

	//inicializacion del arreglo
	for (int i = 0; i < N; i++){
		V[i] = (double)rand()/(double)(RAND_MAX); //funciona en MPI?
		printf (" V[%d] = %f \n",i, V[i]); 
	}

	pthread_t myThreads[T];
	int thread_ids[T];


	double timetick = dwalltime();

	pthread_barrier_init(&barrera, NULL, T+1);	// barrera de T+1 threads (se cuenta el main)

	for (int id = 0; id < T; id++) {
		thread_ids[id] = id;
		pthread_create(&myThreads[id], NULL, &function, (void*) &thread_ids[id]);
	}

	pthread_barrier_wait(&barrera);

	// cuando terminan todos los threads, calculo el maximo de iteraciones

	int max_iter = -1;
	for (int i = 0; i < T; i++) {
		if (iteraciones[i] > max_iter) max_iter = iteraciones[i];
	}

	for (int id = 0; id < T; id++) {
		pthread_join(myThreads[id], NULL);
	}

	printf("Tiempo en segundos %f\n", dwalltime() - timetick);

	for (int i = 0; i < N; i++){
		printf (" V2[%d] = %f \n",i, V2[i]); 
	}
	printf ("iteraciones = %d", max_iter);

	return 0;
}