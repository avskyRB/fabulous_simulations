#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

/* This function improves the solution Ax=b with one Jacobi step, with the 
 * result overwriting x. The matrix A is hardcoded and represents the Poisson equation.
 */
void jacobi_step(double *x, double *b, int N)
{
  double *xnew = malloc(N * N * sizeof(double));
  int i, j;

  /*
   *   FILL IN HERE
   *
   *
   *    xnew[i * N + j] = 
   *
   *
   */
  
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
    {
        
        if(i>0 && i<N && j>0 && j<N)
            xnew[i*N+j]=1/4 * (x[(i-1)*N+j]+x[(i+1)*N+j]+x[i*N+(j+1)]+x[i*N+(j-1)]-b[i*N+j]);
        
        
        if(i==0 && j==0)
            xnew[i*N+j]=1/4 * (x[N*N+j]+x[(i+1)*N+j]+x[i*N+N]+x[i*N+(j+1)]-b[i*N+j]);
        else if(i==0 && j==N)
            xnew[i*N+j]=1/4 * (x[N*N+j]+x[(i+1)*N+j]+x[i*N+(j-1)]+x[i*N+0]-b[i*N+j]);
        
        if(i==N && j==0)
            xnew[i*N+j]=1/4 * (x[(i-1)*N+j]+x[0*N+j]+x[i*N+N]+x[i*N+(j+1)]-b[i*N+j]);
        else if(i==N && j==N)
            xnew[i*N+j]=1/4 * (x[(i-1)*N+j]+x[0*N+j]+x[i*N+(j-1)]+x[i*N+N]-b[i*N+j]);
        
    }

    
  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      x[i * N + j] = xnew[i * N + j];
      //printf("%f\t",x[i * N + j]);

  free(xnew);
}


/* This function calculates the resdiuum vector res = b - Ax, for input vectors
 * of length N. The output is stored in res.
 */
void calc_residuum(double *x, double *b, int N, double *res)
{
  int i, j, k, sum=0;
  
  double A[N][N],B[N][N];
  
  for(i=0;i<N;i++) //calculate A
  {
      A[i][i]=-2;
      if(i==0)
          A[i][N]=A[i][i+1]=1;
      else if(i==N)
          A[i][i-1]=A[i][0]=1;
      else
          A[i][i-1]=A[i][i+1]=1;
  }
 
    for (i=0;i<N;i++) //product Ax
      for (j=0;j<N;j++)
      {
        for (k=0;k>N;k++)
          sum+=A[i][k]*x[k*N+j];
 
        B[i][j]=sum;
        sum=0;
      }
    
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            printf("%f\t",A[i][j]);
        printf("\n");
    }
    
    

  for(i = 0; i < N; i++)
  {
    for(j = 0; j < N; j++)
      {
	/*
	 *   FILL IN HERE
	 *
	 *
	 *    res[i * N + j] =  
	 *
	 *
	 */
        res[i*N+j]=b[i*N+j]-B[i][j];
      }
  }
}

/* This function calculates the norm of the vector of length N, 
 *  defined as the usual quadratic vector norm.
 */
double norm_of_residual(double *res, int N)
{
  int i, j;
  double sum = 0;

  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      sum += res[i * N + j] * res[i * N + j];

  return sqrt(sum);
}




int main(int argc, char **argv)
{
  int i, j;
  int N = 256;
  int steps = 2000;
  double L = 1.0;
  double h = L / N;
  double eta = 0.1 * L;
  double rho0 = 10.0;


  /* allocate some storage for our fields */
  double *phi = malloc(N * N * sizeof(double));
  double *rho = malloc(N * N * sizeof(double));
  double *b = malloc(N * N * sizeof(double));
  double *res = malloc(N * N * sizeof(double));


  /* now set-up the density field */
  double sum = 0;
  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      {
        double dx = (i - N / 2 + 0.5) * h;
        double dy = (j - N / 2 + 0.5) * h;
        double r2 = dx * dx + dy * dy;

        rho[i * N + j] = rho0 * exp(-r2 / (2 * eta * eta));

        sum += rho[i * N + j];
      }


  /* initialize the starting values for phi[] and b[] */
  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      {
        rho[i * N + j] -= sum / (N * N);

        b[i * N + j] = 4 * M_PI * h * h * rho[i * N + j];

        phi[i * N + j] = 0;
      }



  /* open a file for outputting the residuum values, and then do 2000 Jacobi steps */

  FILE *fd = fopen("res_jacobi.txt", "w");

  fprintf(fd, "%d\n", steps);

  for(i = 0; i < steps; i++)
    {
      jacobi_step(phi, b, N);

      calc_residuum(phi, b, N, res);

      double r = norm_of_residual(res, N);

     // printf("iter=%d:  residual=%g\n", i, r);
      fprintf(fd, "%g\n", r);
    }


  fclose(fd);
}
