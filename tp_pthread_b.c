#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>

int T, N;
double *V, *V2;
int *converge;
int *iteraciones;

int trabajando;
int esperando;

pthread_barrier_t barrera;
pthread_mutex_t m1, m2, mC;
pthread_cond_t cond1;

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

void swap_casero(double * V, double * V2, int begin, int end){
	double aux;
	for (int i = begin ; i< end ; i++){
		aux = V[i];
		V[i] = V2[i];
		V2[i] = aux;
}

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

		if (tid == 0) {
			begin++;//mall
			V2local[0] = (Vlocal[0] + Vlocal[1])/2;
		} else if (tid == T-1) {
			end--;
			V2local[N-1] = (Vlocal[N-1] + Vlocal[N-2])/2;
		}
		//reduccion
		for (int i = begin; i < end; i++){
			V2local[i] = (Vlocal[i-1] + Vlocal[i] + Vlocal[i+1])/3;
		}
		//2local[begin] = (Vlocal[begin] + Vlocal[begin+1]) / 2;
		//V2local[end-1] = (Vlocal[end-1] + Vlocal[end-2]) / 2;

		//chequeo de convergencia
		int i = begin;
		double aux = V2[0];

		converge[tid] = 1;
		while ((i < end) && (converge[tid])){ 
			if (fabs(aux - V2local[i]) > 0.01){//si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
				converge[tid] = 0;
			}
			i++;
		}

		iteraciones[tid]++;

		//swap_casero(Vlocal, V2local, begin, end);
		swap(&V, &V2);

		// esto es una especie de barrera para esperar a los threads que están trabajando
		pthread_mutex_lock(&m1);
		esperando++;
		printf("Thread %d dice: esperando = %d, trabajando = %d\n", tid, esperando, trabajando);
		if (esperando == trabajando) {//que pasa si se modifica trabajndo despues de que el ultimo hilo haya hecho la comapracion.
			
			esperando = 0;
			pthread_mutex_unlock(&m1);

			pthread_mutex_lock(&mC);
			pthread_cond_broadcast(&cond1);
			pthread_mutex_unlock(&mC);
		} else {
			
			
			pthread_mutex_unlock(&m1);

			pthread_mutex_lock(&mC);
			pthread_cond_wait(&cond1, &mC);//libera el mutx automaticamente
			pthread_mutex_unlock(&mC);
		}

		printf("Thread %d, iteracion %d\n", tid, iteraciones[tid]);
	}
	printf("THREAD %d TERMINA\n", tid);
	pthread_mutex_lock(&m2);
	trabajando--;
	pthread_mutex_unlock(&m2);

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
	trabajando = T; esperando = 0;

	//inicializacion del arreglo
	for (int i = 0; i < N; i++){
		V[i] = (double)rand()/(double)(RAND_MAX); //funciona en MPI?
		printf (" V[%d] = %f \n",i, V[i]); 
	}

	pthread_t myThreads[T];
	int thread_ids[T];


	double timetick = dwalltime();

	pthread_barrier_init(&barrera, NULL, T+1);	// barrera de T+1 threads (se cuenta el main)
	pthread_mutex_init(&m1, NULL);
	pthread_mutex_init(&m2, NULL);
	pthread_mutex_init(&mC, NULL);

	for (int id = 0; id < T; id++) {
		thread_ids[id] = id;
		pthread_create(&myThreads[id], NULL, &function, (void*) &thread_ids[id]);
	}

	pthread_barrier_wait(&barrera);

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

	pthread_barrier_destroy(&barrera);
	pthread_mutex_destroy(&m1);
	pthread_mutex_destroy(&m2);
	pthread_mutex_destroy(&mC);

	return 0;
}
