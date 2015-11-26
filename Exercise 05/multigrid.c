#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>



/* This function improves the solution Ax=b with one Gauss-Seidel step, with the 
 * result overwriting x. The matrix A is hardcoded and represents the Poisson equation.
 */
void gaussseidel_step(double *x, double *b, int N)
{
  /*
   *   FILL IN HERE
   *
   *
   *    x[i * N + j] = 
   *
   *
   */
  
  int i,j;
    double *xnew = malloc(N * N * sizeof(double));
    
    for(i = 0; i < N; i++) //inizializing xnew
        for(j = 0; j < N; j++)
            xnew[i*N+j]=0;
            
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
            if(i>0 && i<N && j>0 && j<N)
                xnew[i*N + j] = 0.25 * (xnew[(i-1)*N + j] + x[(i+1)*N + j] + x[i*N + (j+1)] + xnew[i*N + (j-1)] - b[i*N + j]);
                
            if(i==0 && j==0)
                xnew[i*N + j] = 0.25 * (xnew[(N-1)*N + j] + x[(i+1)*N + j] + x[i*N + (j+1)] + xnew[i*N + (N-1)] - b[i*N + j]);

            else if(i==0 && j==N-1)
                xnew[i*N + j]= 0.25 * (xnew[(N-1)*N + j] + x[(i+1)*N + j] + x[i*N + (0)] + xnew[i*N + (j-1)] - b[i*N + j]);
        
            if(i==N-1 && j==0)
                xnew[i*N + j]= 0.25 * (xnew[(i-1)*N + j] + x[(0)*N + j] + x[i*N + (j+1)] + xnew[i*N + (N-1)] - b[i*N + j]);

            else if(i==N-1 && j==N-1)
                xnew[i*N + j]= 0.25 * (xnew[(i-1)*N + j] + x[(0)*N + j] + x[i*N + (0)] + xnew[i*N + (j-1)] - b[i*N + j]);  
        }
        
            
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            x[i * N + j] = xnew[i * N + j];

    free(xnew);
  
}


/* This function calculates the resdiuum vector res = b - Ax, for input vectors
 * of length N. The output is stored in res.
 */
void calc_residuum(double *x, double *b, int N, double *res)
{
  /*
   *   FILL IN HERE
   *
   *
   *    res[i * N + j] = 
   *
   *
   */
  
  int i, j;

  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      {
        if(i>0 && i<N && j>0 && j<N)
            res[i*N+j]=b[i*N+j]-(x[(i-1)*N+j]+x[(i+1)*N+j]+x[i*N+(j-1)]+x[i*N+(j+1)]-4*x[i*N+j]);
        
        if(i==0 && j==0)
            res[i*N+j]=b[i*N+j]-(x[(N-1)*N+j]+x[(i+1)*N+j]+x[i*N+(N-1)]+x[i*N+(j+1)]-4*x[i*N+j]);
        else if(i==0 && j==N-1)
            res[i*N+j]=b[i*N+j]-(x[(N-1)*N+j]+x[(i+1)*N+j]+x[i*N+(j-1)]+x[i*N+(0)]-4*x[i*N+j]);
        
        if(i==N-1 && j==0)
            res[i*N+j]=b[i*N+j]-(x[(i-1)*N+j]+x[(0)*N+j]+x[i*N+(N-1)]+x[i*N+(j+1)]-4*x[i*N+j]);
        else if(i==N-1 && j==N-1)
            res[i*N+j]=b[i*N+j]-(x[(i-1)*N+j]+x[(0)*N+j]+x[i*N+(j-1)]+x[i*N+(0)]-4*x[i*N+j]);
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


/* This function restricts the NxN mesh stored in 'fine[ ]' to the NNxNN mesh stored in 'coarse[ ]' 
 */
void do_restrict(int N, double *fine, int NN, double *coarse)
{
    int i, j;
    for(i = 0; i < N; i+=2)
      for(j = 0; j < N; j+=2)
      {
          if(i>0 && i<N-1 && j>0 && j<N-1)
            coarse[i/2 * NN + j/2] =  1/16*fine[(i-1)* N + (j-1)]+ 1/8*fine[i*N + (j-1)] + 1/16*fine[(i+1)*N + (j-1)] + 1/8*fine[(i-1)*N + j] + 1/4*fine[i*N + j] + 1/8*fine[(i+1)*N +j] + 1/16*fine[(i-1)*N + j+1] + 1/8*fine[i*N + j+1] + 1/16*fine[(i+1)*N + j+1];
          if(i==0 && j==0)
            coarse[i/2 * NN + j/2] = 1/16*fine[(i+1)*N + (j+1)] +  1/4*fine[i*N + j] + 1/8*fine[(i+1)*N +j] + 1/8*fine[i*N + j+1];
          else if(i==0 && j==N-1)
            coarse[i/2 * NN + j/2] =  1/8*fine[i*N + (j-1)] + 1/16*fine[(i+1)*N + (j-1)] + 1/4*fine[i*N + j] + 1/8*fine[(i+1)*N +j];
          if(i==N-1 && j==0)
            coarse[i/2 * NN + j/2] =  1/8*fine[(i-1)*N + j] + 1/4*fine[i*N + j]+ 1/16*fine[(i-1)*N + j+1] + 1/8*fine[i*N + j+1];
          else if(i==N-1 && j==N-1)
            coarse[i/2 * NN + j/2] =  1/16*fine[(i-1)* N + (j-1)]+ 1/8*fine[i*N + (j-1)]+ 1/8*fine[(i-1)*N + j] + 1/4*fine[i*N + j];

      }
}

/* This function interpolates the the NNxNN mesh stored in 'coarse[ ]' to NxN mesh stored in 'fine[ ]' 
 */
void do_prolong(int NN, double *coarse, int N, double *fine)
{
  /*
   *   FILL IN HERE
   *
   *
   *   fine[i * N + j] += 
   *
   *
   */
}


/* This function carries out a V-cycle using the multigrid approach, down to N=4.
 * 
 */
void do_v_cycle(double *x, double *b, int N)
{
  gaussseidel_step(x, b, N);

  if(N > 4)
    {
      int i, j, NN = N / 2;

      /* allocate some storage for the residual and error on the current and a coarsened mesh */
      double *res = malloc(N * N * sizeof(double));
      double *err = malloc(N * N * sizeof(double));
      double *res_coarse = malloc(NN * NN * sizeof(double));
      double *err_coarse = malloc(NN * NN * sizeof(double));

      /* calculate the residuum */
      calc_residuum(x, b, N, res);

      /* restrict the residuum */
      do_restrict(N, res, NN, res_coarse);

      /* now multiply the residuum on the coarser mesh (our b)
       * with a factor of 4 to take care of the (h) -> (2h) change, and the fact that we do not change A 
       */
      for(i = 0; i < NN; i++)
        for(j = 0; j < NN; j++)
	  res_coarse[i * NN + j] *= 4;
 

      /* set the coarse error to zero */
      memset(err_coarse, 0, NN * NN * sizeof(double));

      /* call the V-cycle again, recursively */
      do_v_cycle(err_coarse, res_coarse, NN);

      /* now interpolate the error to the finer mesh */
      do_prolong(NN, err_coarse, N, err);

      /* finally add the error to the current solution */
      for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
	  x[i * N + j] += err[i * N + j];
 

      /* free the temporary storage again */
      free(err_coarse);
      free(res_coarse);
      free(err);
      free(res);
    }

  gaussseidel_step(x, b, N);
}



int main(int argc, char **argv)
{
  int i, j;
  int N = 256;
  int NN = N/2;
  int steps = 2000;
  double L = 1.0;
  double h = L / N;
  double eta = 0.1 * L;
  double rho0 = 10.0;
  double *test = malloc(NN * NN * sizeof(double));


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

  
  FILE *fd = fopen("res_multigrid.txt", "w");

  fprintf(fd, "%d\n", steps);

  for(i = 0; i < steps; i++)
    {
      do_v_cycle(phi, b, N);

      calc_residuum(phi, b, N, res);

      double r = norm_of_residual(res, N);

      printf("iter=%d:  residual=%g\n", i, r);

      fprintf(fd, "%g\n", r);
    }

  fclose(fd);

}

