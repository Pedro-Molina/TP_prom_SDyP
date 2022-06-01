#include <stdio.h>
#include <stdlib.h>
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

int main(int argc, char *argv[])
{
	int N;
	double *M, *M2;
	int i, j, k, iteraciones = 0;
	int fila;
	double suma, aux;
	int converge = 0;

	// Controla los argumentos al programa
	if ((argc != 2) || ((N = atoi(argv[1])) <= 0))
	{
		printf("\nUsar: %s n\n  n: matriz de NxN\n", argv[0]);
		exit(1);
	}

	int begin = 1;
	int end = N-1;

	// Aloca memoria para las matrices
	M = (double *)malloc(sizeof(double) * N * N);
	M2 = (double *)malloc(sizeof(double) * N * N);

	// Inicializacion de la matriz
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			M[i * N + j] = (double)rand() / (double)(RAND_MAX); // funciona en MPI?
		}
	}
	double timetick = dwalltime();

	// loop principal
	while (!converge)
	{
		// esquinas
		M2[0] = (M[0] + M[1] + M[N] + M[N + 1]) / 4.0;														 // superior izq
		M2[N - 1] = (M[N - 1] + M[N - 2] + M[2 * N - 1] + M[2 * N - 2]) / 4.0;								 // superior derecha
		M2[(N - 1) * N] = (M[(N - 2) * N] + M[(N - 2) * N + 1] + M[(N - 1) * N] + M[(N - 1) * N + 1]) / 4.0; // inferior izq
		M2[N * N - 1] = (M[N * N - 1] + M[N * N - 2] + M[(N - 1) * N - 1] + M[(N - 1) * N - 2]) / 4.0;		 // inferior der


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

		// borde inferior
		for (j = 1; j < N - 1; j++)
		{
			suma = 0;
			for (i = N - 2; i <= N - 1; i++)
			{
				fila = i * N;
				suma += M[fila + j - 1] + M[fila + j] + M[fila + j + 1];
			}
			M2[(N - 1) * N + j] = suma / 6.0;
		}
		
		
		// reduccion
		for (i = begin; i < end; i++)
		{
			for (j = 1; j < N - 1; j++)
			{
				suma = 0;
				for (k = i - 1; k <= i + 1; k++)
				{
					fila = k * N;
					suma += M[fila + j - 1] + M[fila + j] + M[fila + j + 1];
				}
				M2[i * N + j] = suma / 9.0;
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

		converge = 1;
		i = 0;
		j = 1;
		aux = M2[0];

		// Chequeo de convergencia
		while ((i < N) && (converge))
		{
			while ((j < N) && (converge))
			{
				if (fabs(aux - M2[i*N+j]) > 0.01)
				{ // si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
					converge = 0;
				}
				j++;
			}
			j = 0;
			i++;
		}

		swap(&M, &M2);
		iteraciones++;
	}

	printf("Tiempo en segundos %f\n", dwalltime() - timetick);
	printf("Iteraciones = %d\n", iteraciones);

	free(M);
	free(M2);
	return 0;
}
