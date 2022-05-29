#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>

int N, T;
double *M, *M2;

int *converge;
int convergenciaGlobal = 0, iteraciones =0;
pthread_barrier_t barrera1, barrera2;

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
    int i, j;
	int block = N/ T;
	int begin = block * tid;
	int end = begin + block;
    double suma;
    int fila;

	int endConvergencia = end;

	if (tid == (T - 1)) end--;
	if (tid == 0)
	{
		M2[0] = (M[0] + M[1] + M[N] + M[N + 1]) / 4.0; 
		begin++;
	}

    while (!convergenciaGlobal) 
    {
        if (tid == 0)
        {
            M2[N - 1] = (M[N - 1] + M[N - 2] + M[2 * N - 1] + M[2 * N - 2]) / 4.0; // superior derecha

            // borde superior
            for (j = 1; j < N-1; j++) {
                suma = 0;
                for (i = 0; i <= 1; i++) {
                    fila = i*N;
                    suma += M[fila+j-1] + M[fila+j] + M[fila+j+1];
                }
                M2[j] = suma / 6.0;
            }
        }
        else if (tid == (T - 1))
        {
            M2[(N-1)*N] = (M[(N-2)*N] + M[(N-2)*N+1] + M[(N-1)*N] + M[(N-1)*N+1]) / 4.0;	// inferior izq
            M2[N*N-1] = (M[N*N-1] + M[N*N-2] + M[(N-1)*N-1] + M[(N-1)*N-2]) / 4.0;			// inferior der

            // borde inferior
            for (j = 1; j < N-1; j++) {
                suma = 0;
                for (i = N-2; i <= N-1; i++) {
                    fila = i*N;
                    suma += M[fila+j-1] + M[fila+j] + M[fila+j+1];
                }
                M2[(N-1)*N+j] = suma / 6.0;
            }
        }

        for (i = begin; i < end; i++){
			for (j = 1; j < N - 1; j++) {
				suma = 0;
				for (int k = i-1; k <= i+1; k++) { 
					fila = k*N;
					suma += M[fila+j-1] + M[fila+j] + M[fila+j+1];
				}
				M2[i*N+j] = suma / 9.0;
			}
		}

        // borde izquierdo
        for (i = begin; i < end; i++) {
			suma = 0;
			for (j = 0; j <= 1; j++) {
				suma += M[(i-1)*N+j] + M[i*N+j] + M[(i+1)*N+j];
			}
			M2[i*N] = suma / 6.0;
		}

		// borde derecho
		for (i = begin; i < end; i++) {
			suma = 0;
			for (j = N-2; j <= N-1; j++) {
				suma += M[(i-1)*N+j] + M[i*N+j] + M[(i+1)*N+j];
			}
			M2[i*N+N-1] = suma / 6.0;
		}

        converge[tid] = 1;
		i = begin; j = 0;
		double aux = M2[0];// esperar a que el primero escriba?
		
		//chequeo de convergencia
		while ((i < endConvergencia) && (converge[tid])) {// correccion de end de convergencia
			while ((j < N) && (converge[tid])) {
				if (fabs(aux - M2[i*N+j]) > 0.01){	//si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
						converge[tid] = 0;
				}
				j++;
			}
			j = 0;
			i++;
		}

		pthread_barrier_wait(&barrera1);
		// pasada la barrera, terminaron todos los trabajadores y el main chequea la convergencia
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

			swap(&M, &M2);

			M2[0] = (M[0] + M[1] + M[N] + M[N + 1]) / 4.0;

			iteraciones ++;
		}
		pthread_barrier_wait(&barrera2);
		// pasado el chequeo de convergencia, seguimos trabajando

    }

    pthread_exit(NULL);
}

int main(int argc, const char *argv[])
{

    int i, j;

    if (argc < 3)
	{
		printf("Faltan argumentos. programa N T\n");
		return 1;
	}

	N = atoi(argv[1]);
	T = atoi(argv[2]);

    // Aloca memoria para las matrices
    M = (double *)malloc(sizeof(double) * N * N);
    M2 = (double *)malloc(sizeof(double) * N * N);

    for (i = 0; i < N; i++){
		for (j = 0; j < N; j++) {
			M[i*N+j] = (double)rand()/(double)(RAND_MAX); 
			//printf ("M[%d] = %f \n", i, M[i*N+j]); 
		}
	}
	
    converge = (int *)malloc(sizeof(int) * T);

	for (int i = 0; i < T; i++)
	{
		converge[i] = 0;
	}

    pthread_t myThreads[T];
	int thread_ids[T];

	double timetick = dwalltime();

    pthread_barrier_init(&barrera1, NULL, T); // barrera de T+1 threads (se cuenta el main)
	pthread_barrier_init(&barrera2, NULL, T); // barrera de T+1 threads (se cuenta el main)


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

    /*for (i = 0; i < N; i++){
		for (j = 0; j < N; j++) {
			//M[i*N+j] = (double)rand()/(double)(RAND_MAX); //funciona en MPI?
			printf ("M[%d] = %f ", i, M[i*N+j]); 
		}
        printf("\n");
	}*/

	pthread_barrier_destroy(&barrera1);
	pthread_barrier_destroy(&barrera2);

	return 0;
}