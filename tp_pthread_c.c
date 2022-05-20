//FUNCIONA PERO MIENTRAS SE INCEMENTAN LOS HILOS INCREMENTA EL TIEMPO DE EJECUCION, CAPAZ MEJORA IMPLEMENTANDO BARRERAS ENTRE HILOS (??)

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

pthread_barrier_t barrera1, barrera2;
int convergencia = 0;

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
	//int endConvergencia = end;
	if (tid == 0) begin++;//se modificaba en cada loop MAL
	if (tid == (T - 1))
	{
		end--;
	} 

	// double *Vlocal, *V2local;
	//  Vlocal = V;//esta asugnacion se hace una sola vez cuando el main swapwea v2 y v1 no se ve reflejado en los valores locales.
	//  V2local = V2;

	//printf("Thread %d arranca. begin = %d, end = %d\n", tid, begin, end);

	while (!convergencia) // los hilos deben seguir recalculando los valores si no convergio todo el vector, como en el secuencial, no se puede determinar la convergencia por partes
	{
		if (tid == 0)
		{
			V2[0] = (V[0] + V[1]) / 2;
		}
		else if (tid == (T - 1))
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
		while ((i < end) && (converge[tid])) // el ultimo no compara el ultimo valor pero funciona igual (??)
		{
			if (fabs(aux - V2[i]) > 0.01) // si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
			{
				converge[tid] = 0;
				// printf("valor abs = %f", fabs(aux - V2[i]));
			}
			i++;
		}

		iteraciones[tid]++;

		// swap(&V, &V2);
		// printf("Thread %d, iteracion %d\n", tid, iteraciones[tid]);
		// if (iteraciones[tid] % 100000 == 0)
		// printf("Thread %d termina iteracion %d\n", tid, iteraciones[tid]);
		
		pthread_barrier_wait(&barrera1);
		// pasada la barrera, terminaron todos los trabajadores y el main chequea la convergencia
		pthread_barrier_wait(&barrera2);
		// pasado el chequeo de convergencia, seguimos trabajando
	}

	//while (!convergencia)
	//{
		//pthread_barrier_wait(&barrera1);
		//pthread_barrier_wait(&barrera2);
	//}

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

	// Aloca memoria para las matrices
	V = (double *)malloc(sizeof(double) * N);
	V2 = (double *)malloc(sizeof(double) * N);

	converge = (int *)malloc(sizeof(int) * T);
	iteraciones = (int *)malloc(sizeof(int) * T);
	for (int i = 0; i < T; i++)
	{
		converge[i] = 0;
		iteraciones[i] = 0;
	}
	// trabajando = T; esperando = 0;

	// inicializacion del arreglo
	for (int i = 0; i < N; i++)
	{
		V[i] = (double)rand() / (double)(RAND_MAX); // funciona en MPI?
		//printf(" V[%d] = %f \n", i, V[i]);
	}

	pthread_t myThreads[T];
	int thread_ids[T];

	double timetick = dwalltime();

	pthread_barrier_init(&barrera1, NULL, T + 1); // barrera de T+1 threads (se cuenta el main)
	pthread_barrier_init(&barrera2, NULL, T + 1); // barrera de T+1 threads (se cuenta el main)

	for (int id = 0; id < T; id++)
	{
		thread_ids[id] = id;
		pthread_create(&myThreads[id], NULL, &function, (void *)&thread_ids[id]);
	}

	convergencia = 0;
	while (!convergencia)
	{
		pthread_barrier_wait(&barrera1);

		// chequear convergencia
		convergencia = 1;
		int i = 0;
		while ((i < T) && (convergencia))
		{
			if (converge[i++] == 0)
			{
				convergencia = 0;
			}
		}

		//printf("Swap \n");
		swap(&V, &V2);

		pthread_barrier_wait(&barrera2);
	}
	int max_iter = -1;
	for (int i = 0; i < T; i++)
	{
		if (iteraciones[i] > max_iter)
			max_iter = iteraciones[i];
	}

	for (int id = 0; id < T; id++)
	{
		pthread_join(myThreads[id], NULL);
	}

	printf("Tiempo en segundos %f\n", dwalltime() - timetick);

	/*for (int i = 0; i < N; i++)
	{
		printf(" V2[%d] = %f \n", i, V2[i]);
	}*/
	printf("iteraciones = %d\n", max_iter);


	/*for (int i = 0; i < T; i++)
	{
		printf(" converge[%d] = %d \n", i, converge[i]);
	}*/

	pthread_barrier_destroy(&barrera1);
	pthread_barrier_destroy(&barrera2);

	return 0;
}