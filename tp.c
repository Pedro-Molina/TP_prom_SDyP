#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

void swap( double ** x, double ** y)
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
	double *V, *V2;
	int i, iteraciones = 0;
	int converge =0;
	double timetick, aux;

 	// Controla los argumentos al programa
	if ((argc != 2) || ((N = atoi(argv[1])) <= 0))
  {
    printf("\nUsar: %s n\n  n: Longitud del vector\n", argv[0]);
    exit(1);
  }

 	// Aloca memoria para los arreglos
  V = (double*) malloc(sizeof(double)*N);
	V2 = (double*) malloc(sizeof(double)*N);
	
	// Inicializacion del arreglo
	for (i = 0; i < N; i++)
	{
		V[i] = (double) rand() / (double)(RAND_MAX); //funciona en MPI?
	}
	
	timetick = dwalltime();
	// Loop principal
	while (!converge)
	{
	
		// Reduccion
		for (i = 1; i < N - 1; i++)
		{
			V2[i] = (V[i-1] + V[i] + V[i+1]) / 3.0;
		}
		
		// Valores de los extremos
		V2[0] = (V[0] + V[1]) / 2.0;
		V2[N-1] = (V[N-1] + V[N-2]) / 2.0;
		
		converge = 1;
		i = 1;
		aux = V2[0];

		// Chequeo de convergencia
		while ((i < N) && (converge))
		{ 
			if ((fabs(aux - V2[i])) > 0.01){ // si la diferencia en mayor a 0.01 el arreglo no llego a la convergencia
					converge = 0;
			}
			i++;
		}
		
		swap(&V, &V2);

		iteraciones++;
	}
	
	printf("Tiempo en segundos: %f\n", dwalltime() - timetick);
	
	printf ("Iteraciones = %d\n", iteraciones);

	/*for (int i = 0; i < N; i++) {
  	printf("V[%d] = %f\n", i, V[i]);
	}*/

}
