#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>

#define N 1024

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

void root_process(int size)
{
    double *M, *M2, *Maux;
    int block_size = (N * N) / size;
    int converge, iteraciones = 0;
    int convergeGlobal = 0;
    double timetick,suma;
    int fila,i,j,k;
    int begin = 1;
    int end = N/size;

    // Aloca memoria para los vectores
    M = (double *)malloc(sizeof(double) * N * N);
    M2 = (double *)malloc(sizeof(double) * (block_size + N) );
    Maux = M;
    // Inicializacion del arreglo
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            M[i * N + j] = (double)rand() / (double)(RAND_MAX); // funciona en MPI?
        }
    }
    MPI_Barrier(MPI_COMM_WORLD); //barrera para exlcluir el tiempo de alocacion de memoria del tiempo del algoritmo
    timetick = dwalltime();

    // Enviar los bloques a cada proceso

    //MPI_Scatter(message, N/nProcs, MPI_CHAR, part, N/nProcs, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Scatter(M, block_size, MPI_DOUBLE, M2, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    while (!convergeGlobal) 
    {
        MPI_Request request;
		MPI_Status status;

        MPI_Isend(M+block_size-N, N, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &request);
        // recibo el valor del vecino derecho
		MPI_Irecv(M+block_size, N, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);

        // Reduccion
        M2[0] = (M[0] + M[1] + M[N] + M[N + 1]) / 4.0;                         // superior izq
        M2[N - 1] = (M[N - 1] + M[N - 2] + M[2 * N - 1] + M[2 * N - 2]) / 4.0; // superior derecha

        // borde superior
        for (j = 1; j < N - 1; j++)
        {
            suma = 0;
            for (i = 0; i <= 1; i++)
            {
                fila = i * N;
                suma += M[fila + j - 1] + M[fila + j] + M[fila + j + 1];
            }
            M2[j] = suma / 6.0;
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
        for (i = begin; i < end; i++)
        {
            suma = 0;
            for (j = 0; j <= 1; j++)
            {
                suma += M[(i - 1) * N + j] + M[i * N + j] + M[(i + 1) * N + j];
            }
            M2[i * N] = suma / 6.0;
        }

        // borde derecho
        for (i = begin; i < end; i++)
        {
            suma = 0;
            for (j = N - 2; j <= N - 1; j++)
            {
                suma += M[(i - 1) * N + j] + M[i * N + j] + M[(i + 1) * N + j];
            }
            M2[i * N + N - 1] = suma / 6.0;
        }

        MPI_Bcast(M2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        
        i = 0;
        j = 1;
        converge =1;
        double aux = M2[0];
        while ((i < end) && (converge)) { 
			while ((j < N) && (converge)) {
				if (fabs(aux - M2[i*N+j]) > 0.01){	//si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
					converge = 0;
				}
				j++;
			}
			j = 0;
			i++;
		}

        iteraciones++;

        
		swap(&M, &M2);

        // Chequeo de convergencia global
        MPI_Allreduce(&converge, &convergeGlobal, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    }

    //MPI_Gather(part, N/nProcs, MPI_CHAR, message, N/nProcs, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Gather(M, block_size, MPI_DOUBLE, Maux, block_size , MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("Tiempo en segundos: %f\n", dwalltime() - timetick);
    printf("Iteraciones: %d\n", iteraciones);

    free(M);
	free(M2);

}

void worker_process(int rank, int size)
{
    double *M, *M2;
    int block_size = (N * N) / size + 2 * N;
    int converge;
    int convergeGlobal = 0;
    double timetick, aux, suma;
    int begin = 1;
    int end = N/size +1;
    int endConvergencia = end;
    int fila,i,j,k;

    // Aloca memoria para los vectores
    if (rank == size - 1)
    {
        block_size -= N;
        end--;
    }
    M = (double *)malloc(sizeof(double) * block_size);
    M2 = (double *)malloc(sizeof(double) * block_size);

    MPI_Barrier(MPI_COMM_WORLD); //barrera para exlcluir el tiempo de alocacion de memoria del tiempo del algoritmo
    // Recibir el bloque
    //MPI_Scatter(message, N/nProcs, MPI_CHAR, part, N/nProcs, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Scatter(M+N,(N * N) / size, MPI_DOUBLE, M+N, (N * N) / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    while (!convergeGlobal) 
  {
	  // Reduccion
        MPI_Request request;
		MPI_Status status;

        if (rank != size-1) {
			// envio mis valores a los vecinos
			MPI_Isend(M+N, N, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			MPI_Isend(M+block_size-2*N, N, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &request);

			// recibo el valor del vecino izquierdo
			MPI_Irecv(M, N, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
			// recibo el valor del vecino derecho
			MPI_Irecv(M+block_size-N, N, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);

		} else {
			MPI_Isend(M+N, N, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			// recibo el valor del vecino izquierdo
			MPI_Irecv(M, N, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
		}

        if (rank == size-1) {

            M2[(end)*N] = (M[(end-1)*N] + M[(end-1)*N+1] + M[(end)*N] + M[(end)*N+1]) / 4.0;	// inferior izq
            M2[(end)*N + N-1] = (M2[(end)*N + N-1] + M[(end)*N + N-2] + M[(end-1)*N+ N-1] + M[(end-1)*N+ N-2]) / 4.0;				// inferior der

            // borde inferior
            for (j = 1; j < N-1; j++) {
                suma = 0;
                for (i = end-1; i <= end; i++) {
                    fila = i*N;
                    suma += M[fila+j-1] + M[fila+j] + M[fila+j+1];
                }
                M2[(end)*N+j] = suma / 6.0;
            }
        }
	    for (i = begin; i < end ; i++){
			for (j = 1; j < N - 1; j++) {
				suma = 0;
				for (int k = i-1; k <= i+1; k++) { 
					fila = k*N;
					suma += M[fila+j-1] + M[fila+j] + M[fila+j+1];
				}
				M2[i*N+j] = suma / 9.0;
			}
		}

        // borde izquierdo //mejora: recorrer por filas
        for (i = begin; i < end; i++)
        {
            suma = 0;
            for (j = 0; j <= 1; j++)
            {
                suma += M[(i - 1) * N + j] + M[i * N + j] + M[(i + 1) * N + j];
            }
            M2[i * N] = suma / 6.0;
        }

        // borde derecho
        for (i = begin; i < end; i++)
        {
            suma = 0;
            for (j = N - 2; j <= N - 1; j++)
            {
                suma += M[(i - 1) * N + j] + M[i * N + j] + M[(i + 1) * N + j];
            }
            M2[i * N + N - 1] = suma / 6.0;
        }

		// Recibo M2[0] en aux
		MPI_Bcast(&aux, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        i = 1;
        j = 0;
        converge =1;
        
        while ((i < endConvergencia) && (converge)) { 
			while ((j < N) && (converge)) {
				if (fabs(aux - M2[i*N+j]) > 0.01){	//si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
						converge = 0;
				}
				j++;
			}
			j = 0;
			i++;
		}
        
		swap(&M, &M2);

        // Chequeo de convergencia global
		MPI_Allreduce(&converge, &convergeGlobal, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        
	}

	MPI_Gather(M +N, (N * N) / size, MPI_DOUBLE, M+N, (N * N) / size , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    free(M);
	free(M2);

}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
    {
        root_process(size);
    }
    else
    {
        worker_process(rank, size);
    }

    MPI_Finalize();

    return 0;
}
