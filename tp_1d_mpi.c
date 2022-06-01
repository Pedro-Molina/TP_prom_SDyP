#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>

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

void root_process(int size, int N) {
	double *V, *V2, *Vaux;
	int block_size = N / size;
	int converge, iteraciones = 0;
	int convergeGlobal = 0;
	double timetick, aux;

	// Aloca memoria para los vectores
	V = (double *) malloc(sizeof(double) * N);
	V2 = (double *) malloc(sizeof(double) * (block_size +1));
	Vaux = V;	// necesario para el gather

	// Inicializacion del arreglo
	for (int i = 0; i < N; i++)
	{
		V[i] = (double)rand() / (double)(RAND_MAX);
	}

	MPI_Barrier(MPI_COMM_WORLD); //barrera para excluir el tiempo de alocacion de memoria del tiempo del algoritmo
	timetick = dwalltime();

	// Enviar los bloques a cada proceso
	MPI_Scatter(V, block_size, MPI_DOUBLE, V, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	while (!convergeGlobal) 
	{
		MPI_Request request;
		MPI_Status status;

		// envio el valor al vecino derecho
		MPI_Isend(V+block_size-1, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &request);
		// recibo el valor del vecino derecho
		MPI_Irecv(V+block_size, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);
		
		// Reduccion
		for (int i = 1; i < block_size; i++)
		{
			V2[i] = (V[i - 1] + V[i] + V[i + 1]) / 3.0;
		}
		V2[0] = (V[0] + V[1]) / 2.0;
	
		// broadcast de V2[0]
		MPI_Bcast(V2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
		// Chequeo de convergencia
		converge = 1;
		int i = 1;
		aux = V2[0];
		while ((i < block_size) && (converge)) 
		{
			if (fabs(aux - V2[i]) > 0.01) // si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
			{
				converge = 0;
			}
			i++;
		}
	
		iteraciones++;

		swap(&V, &V2);

		// Chequeo de convergencia global
		MPI_Allreduce(&converge, &convergeGlobal, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
	}

	MPI_Gather(V, block_size, MPI_DOUBLE, Vaux, block_size , MPI_DOUBLE, 0, MPI_COMM_WORLD); 

	printf("Tiempo en segundos: %f\n", dwalltime() - timetick);
  	printf("Iteraciones: %d\n", iteraciones);

  	free(V);
	free(V2);
  
}

void worker_process(int rank, int size, int N) {
	double *V, *V2;
	int block_size = N / size + 2;
	int endConvergencia = block_size - 1;
	int converge, iteraciones = 0;
	int convergeGlobal = 0;
	double aux;

	// Aloca memoria para los vectores
	if (rank == size-1) {
		block_size--;
	}
	V = (double *) malloc(sizeof(double) * block_size);
	V2 = (double *) malloc(sizeof(double) * block_size);

	MPI_Barrier(MPI_COMM_WORLD); //barrera para excluir el tiempo de alocacion de memoria del tiempo del algoritmo

	// Recibir el bloque
	MPI_Scatter(V+1, N / size, MPI_DOUBLE, V+1, N / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	while (!convergeGlobal) 
	{
		MPI_Request request;
		MPI_Status status;

		if (rank != size-1) {
			// envio mis valores a los vecinos
			MPI_Isend(V+1, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			MPI_Isend(V+block_size-2, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &request);

			// recibo el valor del vecino izquierdo
			MPI_Irecv(V, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
			// recibo el valor del vecino derecho
			MPI_Irecv(V+block_size-1, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);

		} else {
			// envio el valor al vecino izquierdo
			MPI_Isend(V+1, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			// recibo el valor del vecino izquierdo
			MPI_Irecv(V, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
		}

		// Reduccion
		for (int i = 1; i < block_size-1; i++)
		{
			V2[i] = (V[i - 1] + V[i] + V[i + 1]) / 3.0;
		}
		if (rank == size-1) {
			V2[block_size - 1] = (V[block_size - 1] + V[block_size - 2]) / 2.0;
		}

		// Recibo V2[0] en aux
		MPI_Bcast(&aux, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// Chequeo de convergencia
		converge = 1;
		int i = 1;
		while ((i < endConvergencia) && (converge)) 
		{
			if (fabs(aux - V2[i]) > 0.01) // si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
			{
				converge = 0;
			}
			i++;
		}

		swap(&V, &V2);

		// Chequeo de convergencia global
		MPI_Allreduce(&converge, &convergeGlobal, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

	}
	MPI_Gather(V+1, N / size, MPI_DOUBLE, V, N / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	free(V);
	free(V2);
}


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	int N;
	if ((argc != 2) || ((N = atoi(argv[1])) <= 0))
	{
		printf("\nUsar: %s n\n  n: matriz de NxN\n", argv[0]);
		exit(1);
	}

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		root_process(size, N);
	} else {
		worker_process(rank, size, N);
	}

	MPI_Finalize();

	return 0;
}
