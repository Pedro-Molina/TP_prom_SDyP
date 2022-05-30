#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>

int T, N;
double *V, *V2;
int *converge;
int iteraciones = 0;

int trabajando;
int esperando;

pthread_barrier_t barrera1, barrera2;
int convergenciaGlobal = 0;

void swap(double **x, double **y)
{
	double *temp = *x;
	*x = *y;
	*y = temp;
}

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

void *function(void *arg)
{
	int tid = *(int *)arg;

	int block = N / T;
	int begin = block * tid;
	int end = begin + block;
	int endConvergencia = end;

	if (tid == 0) begin++;
	if (tid == (T - 1))
	{
		end--;
	} 

	if (tid == 0)
	{
		V2[0] = (V[0] + V[1]) / 2;
	}

	pthread_barrier_wait(&barrera1); // los hilos se deben esperar para que todos tengan el valor de V[0] en la primera iteracion

	while (!convergenciaGlobal) // los hilos deben seguir recalculando los valores si no convergio todo el vector, como en el secuencial, no se puede determinar la convergencia por partes
	{
		
		if (tid == (T - 1))
		{
			V2[N - 1] = (V[N - 1] + V[N - 2]) / 2;
		}
		// reduccion
		for (int i = begin; i < end; i++)
		{
			V2[i] = (V[i - 1] + V[i] + V[i + 1]) / 3;
		}

		// chequeo de convergencia
		int i = begin; // el primero no compara el primer valor (no pasa nada)
		double aux = V2[0];

		converge[tid] = 1;
		while ((i < endConvergencia) && (converge[tid])) 
		{
			if (fabs(aux - V2[i]) > 0.01) // si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
			{
				converge[tid] = 0;
			}
			i++;
		}
	
		pthread_barrier_wait(&barrera1);

		if (tid == 0)
		{
			// chequear convergencia
			convergenciaGlobal = 1;
			int i = 0;
			while ((i < T) && (convergenciaGlobal))
			{
				if (converge[i++] == 0)
				{
					convergenciaGlobal = 0;
				}
			}

			swap(&V, &V2);

			// calculo V[0] para que en la proxima iteracion todos tengan el valor actualizado
			V2[0] = (V[0] + V[1]) / 2;

			iteraciones ++;
		}

		pthread_barrier_wait(&barrera2);
	}

	pthread_exit(NULL);
}

int main(int argc, const char *argv[])
{

	if (argc < 3)
	{
		printf("Faltan argumentos. programa N T\n");
		return 1;
	}

	N = atoi(argv[1]);
	T = atoi(argv[2]);

	// Aloca memoria para los vectores
	V = (double *) malloc(sizeof(double) * N);
	V2 = (double *) malloc(sizeof(double) * N);

	converge = (int *)malloc(sizeof(int) * T);

	for (int i = 0; i < T; i++)
	{
		converge[i] = 0;
	}

	// inicializacion del arreglo
	for (int i = 0; i < N; i++)
	{
		V[i] = (double)rand() / (double)(RAND_MAX); 
	}

	pthread_t myThreads[T];
	int thread_ids[T];

	double timetick = dwalltime();

	pthread_barrier_init(&barrera1, NULL, T); // barrera de T threads
	pthread_barrier_init(&barrera2, NULL, T); // barrera de T threads

	for (int id = 0; id < T; id++)
	{
		thread_ids[id] = id;
		pthread_create(&myThreads[id], NULL, &function, (void *)&thread_ids[id]);
	}
	

	for (int id = 0; id < T; id++)
	{
		pthread_join(myThreads[id], NULL);
	}

	printf("Tiempo en segundos %f\n", dwalltime() - timetick);

	printf("Iteraciones = %d\n", iteraciones);

	pthread_barrier_destroy(&barrera1);
	pthread_barrier_destroy(&barrera2);

	return 0;
}