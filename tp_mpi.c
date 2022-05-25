//FUNCIONA PERO MIENTRAS SE INCEMENTAN LOS HILOS INCREMENTA EL TIEMPO DE EJECUCION, CAPAZ MEJORA IMPLEMENTANDO BARRERAS ENTRE HILOS (??)

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>

#define N 512

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

void root_process(int size) {
	double *V, *V2;
	int block_size = N / size;
	int converge, iteraciones = 0;
	int convergeGlobal = 0;
	double timetick;

	printf("Root\n");

	// Aloca memoria para los vectores
	V = (double *) malloc(sizeof(double) * N);
	V2 = (double *) malloc(sizeof(double) * block_size+1);

	// Inicializacion del arreglo
	for (int i = 0; i < N; i++)
	{
		V[i] = (double)rand() / (double)(RAND_MAX); // funciona en MPI?
		printf("V[%d] = %f\n", i, V[i]);
	}

	timetick = dwalltime();

	// Enviar los bloques a cada proceso
  for (int id = 1; id < size-1; id++) {
    MPI_Send(V+id*block_size-1, block_size+2, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
  }
  MPI_Send(V+(size-1)*block_size-1, block_size+1, MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD);

  printf("Enviado todo\n");


  while (!convergeGlobal) 
  {
  	// Reduccion
  	for (int i = 1; i < block_size; i++)
  	{
  		V2[i] = (V[i - 1] + V[i] + V[i + 1]) / 3.0;
  	}
  	V2[0] = (V[0] + V[1]) / 2.0;
  
   	MPI_Bcast(V2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  	// Chequeo de convergencia
  	converge = 1;
  	int i = 1;
  	while ((i < block_size) && (converge)) // el ultimo no compara el ultimo valor pero funciona igual (??)
  	{
  		if (fabs(V2[0] - V2[i]) > 0.01) // si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
  		{
  			converge = 0;
  		}
  		i++;
  	}
  
  	iteraciones++;

  	// Valores compartidos
  	MPI_Request request;
		MPI_Status status;

  	MPI_Isend(V2+block_size-1, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &request);
  	// recibo el valor del vecino derecho
		MPI_Irecv(V2+block_size, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);
  
		swap(&V, &V2);

  	// Chequeo de convergencia global
  	MPI_Allreduce(&converge, &convergeGlobal, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  }

  int max_iteraciones;
  MPI_Reduce(&iteraciones, &max_iteraciones, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);


  printf("Iteraciones: %d\n", max_iteraciones);

	//MPI_Gather(V, block_size, MPI_DOUBLE, V, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	printf("Tiempo en segundos: %f\n", dwalltime() - timetick);
  /*for (int i = 0; i < N; i++) {
  	printf("V[%d] = %f\n", i, V[i]);
	}*/
}

void worker_process(int rank, int size) {
	double *V, *V2;
	int block_size = N / size + 2;
	int converge, iteraciones = 0;
	int convergeGlobal = 0;
	double aux;

	printf("Hola %d\n", rank);

	// Aloca memoria para los vectores
	if (rank == size-1) {
		block_size--;
	}
	V = (double *) malloc(sizeof(double) * block_size);
	V2 = (double *) malloc(sizeof(double) * block_size);

	// Recibir el bloque
  MPI_Recv(V, block_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  printf("%d: Recibi mi vector\n", rank);

  while (!convergeGlobal) 
  {
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
		while ((i < block_size-1) && (converge)) // el ultimo no compara el ultimo valor pero funciona igual (??)
		{
			if (fabs(aux - V2[i]) > 0.01) // si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
			{
				converge = 0;
			}
			i++;
		}

		iteraciones++;

		MPI_Request request;
		MPI_Status status;

		if (rank != size-1) {
			// envio mis valores a los vecinos
			MPI_Isend(V2+1, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			MPI_Isend(V2+block_size-2, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &request);

			// recibo el valor del vecino izquierdo
			MPI_Irecv(V2, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
			// recibo el valor del vecino derecho
			MPI_Irecv(V2+block_size-1, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);

		} else {
			MPI_Isend(V2+1, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			// recibo el valor del vecino izquierdo
			MPI_Irecv(V2, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
		}

		swap(&V, &V2);

		// Chequeo de convergencia global
		MPI_Allreduce(&converge, &convergeGlobal, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

	}

	int max_iteraciones;
  MPI_Reduce(&iteraciones, &max_iteraciones, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  MPI_Gather(V+1, N / size, MPI_DOUBLE, V, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}



int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0) {
    root_process(size);
  } else {
    worker_process(rank, size);
  }
 
  MPI_Finalize();
  
  return 0;
}
