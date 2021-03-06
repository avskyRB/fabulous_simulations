#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#define sigma 3.4*pow(10,-10)
#define kb 1.38*pow(10,-23)
#define eps 1.65*pow(10,-21)



/* this is our data type for storing the information for the particles
 */
typedef struct
{
  double pos[3];
  double vel[3];
  double acc[3];
  double pot;
}
particle;



/* auxiliary function to create a Gaussian random deviate
 */
double gaussian_rnd(void)
{
  double x, y;

  do
    {
      x = drand48();
      y = drand48();
    }
  while(x == 0);

  return sqrt(-2.0 * log(x)) * cos(2 * M_PI * y);
}


/* This function initializes our particle set.
 */
void initialize(particle * p, int nperdim, double boxsize, double temp, double epsilon_in_Kelvin)
{
  int i, j, k, n = 0;

  for(i = 0; i < nperdim; i++)
    for(j = 0; j < nperdim; j++)
      for(k = 0; k < nperdim; k++)
        {
	  /*
	   *   FILL IN HERE
	   *   
	   *    p[n].pos[0] = ........
	   *    p[n].pos[1] = ........
	   *    p[n].pos[2] = ........
	   *
	   *    p[n].vel[0] = .......
	   *    .......
	   */

          p[n].pos[0]=5.0*i+boxsize/sigma;
          p[n].pos[1]=5.0*j+boxsize/sigma;
          p[n].pos[2]=5.0*k+boxsize/sigma;
          
          p[n].vel[0]=gaussian_rnd();
          p[n].vel[1]=gaussian_rnd();
          p[n].vel[2]=gaussian_rnd();
          
          n++;
        }
}


/* This function updates the velocities by applying the accelerations for the given time interval.
 */
void kick(particle *p, int ntot, double dt)
{
  int i, k;

  for(i = 0; i < ntot; i++)
    for(k = 0; k < 3; k++)
      p[i].vel[k] += p[i].acc[k] * dt;
}



/* This function drifts the particles with their velocities for the given time interval.
 * Afterwards, the particles are mapped periodically back to the box if needed.
 */
void drift(particle * p, int ntot, double boxsize, double dt)
{
    int n,i;
    int modi[ntot][3],modf[ntot][3];
  /*
   *   FILL IN HERE
   *   
   *    p[n].pos[0] += ........
   *    ........
   */
  
 /* for(n=0;n<ntot;n++)
      for(i=0;i<3;i++)
        modi[n][i]=p[n].pos[i]/(5*sigma);*/
  
  for(n=0;n<ntot;n++)
      for(i=0;i<3;i++)
      {
        p[n].pos[i]+=p[n].vel[i]*dt;
        //modf[n][i]=p[n].pos[i]/(5*sigma);
        
        /*if(p[n].pos[i]>modi[n][i])
            p[n].pos[i]-=(5*sigma)*modf[n][i];
        if(p[n].pos[i]<modi[n][i])
            p[n].pos[i]+=(5*sigma)*modf[n][i];*/
        
        if(p[n].pos[i]>boxsize/sigma)
            p[n].pos[i]=fabs(p[n].pos[i]-boxsize/sigma);
        if(p[n].pos[i]<boxsize/sigma)
            p[n].pos[i]=boxsize/sigma-p[n].pos[i];
      }
}



/* This function calculates the potentials and forces for all particles. For simplicity,
 * we do this by going through all particle pairs i-j, and then adding the contributions both to i and j.
 */
void calc_forces(particle * p, int ntot, double boxsize, double rcut)
{
  int i, j, k;
  double rcut2 = rcut * rcut;
  double r2, r, r6, r12, dr[3], acc[3], pot;


  /* first, set all the accelerations and potentials to zero */
  for(i = 0; i < ntot; i++)
    {
      p[i].pot = 0;

      for(k = 0; k < 3; k++)
        p[i].acc[k] = 0;
    }


  /* sum over all distinct pairs */
  for(i = 0; i < ntot; i++)
    for(j = i + 1; j < ntot; j++)
      {
        for(k = 0, r2 = 0; k < 3; k++)
          {
            dr[k] = p[i].pos[k] - p[j].pos[k];

            if(dr[k] > 0.5 * boxsize)
              dr[k] -= boxsize;

            if(dr[k] < -0.5 * boxsize)
              dr[k] += boxsize;

            r2 += dr[k] * dr[k];
          }


	/*
	 *   FILL IN HERE
	 *   
	 *    p[i].pot += 
         *    p[j].pot += 
	 *   
         *    p[i].acc[k] += ...
         *    p[j].acc[k] += ...
	 *    ........
	 */
        
            if(sqrt(r2)>rcut)
                p[i].pot=p[j].pot=0;
        
            p[i].pot+=4*(pow(1/sqrt(r2),12.0)-pow(1/sqrt(r2),6.0));
            p[j].pot+=4*(pow(1/sqrt(r2),12.0)-pow(1/sqrt(r2),6.0));
        
            for(k=0;k<3;k++)
            {
                p[i].acc[k]+=48/pow(r2,6.5)-24/pow(r2,3.5);//if you consider two particles one gets an acceleration,
                p[j].acc[k]-=48/pow(r2,6.5)-24/pow(r2,3.5);//the other gets the opposite

                /* APPARENTLY INTEGERS WERE NOT THE THING TO USE HERE*/
            }
      }
}


/* This function calculates the total kinetic and total potential energy, averaged per particle.
 * It also returns the instantanous kinetic temperature.
 */
void calc_energies(particle * p, int ntot, double epsilon_in_Kelvin, double *ekin, double *epot, double *temp)
{
    int n,i;
  /*
   *   FILL IN HERE
   *   
   *     *ekin = ....
   *     *epot = ....
   *     *temp = ....
   */
  for(n=0;n<ntot;n++)
  {
      *ekin+=0.5*(pow(p[n].vel[0],2)+pow(p[n].vel[1],2)+pow(p[n].vel[2],2));
      *epot+=p[n].pot;
      *temp=pow(p[n].vel[0],2)+pow(p[n].vel[1],2)+pow(p[n].vel[2],2);
  }
  *temp *= 2/(3*kb) * 1/ntot * 0.5;
}



/* This function rescales the velocities by a prescribed factor 'fac'
 */
void rescale_velocities(particle * p, int ntot, double fac)
{
    int n,i;
  /*
   *   FILL IN HERE
   *   
   *      p[i].vel[k] *= ...
   *   
   */
  for(n=0;n<ntot;n++)
      for(i=0;i<3;i++)
          p[n].vel[i]*=fac;
}


/*
 * main driver routine
 */
int main(int argc, char **argv)
{
  double epsilon_in_Kelvin = 120.0;     /* energy scale in Lennard Jones potential       */
  int N1d = 8;                          /* particles per dimension                       */
  int N = N1d * N1d * N1d;              /* total particle number                         */
  double mean_spacing = 5.0;            /* mean spacing in dimensionaless internal units */
  double boxsize = N1d * mean_spacing;  /* dimensionless boxsize                         */
  double target_temp_in_Kelvin = 80.0;  /* target temperature                            */
  int nsteps = 60000;                   /* number of steps to take                       */
  double dt = 0.01;                     /* timestep size                                 */

  double ekin, epot, temp;
  int step;

  /* allocate storage for our particles */
  particle *p = malloc(N * sizeof(particle));

  /* let's initialize the particles */
  initialize(p, N1d, boxsize, target_temp_in_Kelvin, epsilon_in_Kelvin);

  /* calculate the forces at t=0 */
  calc_forces(p, N, boxsize, 10.0);

  /* create an output file */
  char fname[100];
  sprintf(fname, "output_%d.txt", (int) target_temp_in_Kelvin);
  FILE *fd = fopen(fname, "w");

  /* measure energies at beginning, and output this to the file and screen */
  calc_energies(p, N, epsilon_in_Kelvin, &ekin, &epot, &temp);
  fprintf(fd, "%6d   %10g   %10g   %10g       %10g\n", 0, ekin, epot, ekin + epot, temp);
  printf("%6d   %10g   %10g   %10g       %10g\n", 0, ekin, epot, ekin + epot, temp);


  /* now we carry out time integration steps using the leapfrog */
  for(step = 0; step < nsteps; step++)
    {
      kick(p, N, 0.5 * dt);

      drift(p, N, boxsize, dt);

      calc_forces(p, N, boxsize, 10.0);

      kick(p, N, 0.5 * dt);

      /* every 100-th step, we call the rescaling routine */
      if((step % 100) == 0)
        {
	  /*
           *   FILL IN HERE
           *   double fac = ....
           */      
          double fac=sqrt(target_temp_in_Kelvin/temp);
     	  rescale_velocities(p, N, fac);
        }

      /* measure energies and output this to the file and screen */
      calc_energies(p, N, epsilon_in_Kelvin, &ekin, &epot, &temp);
      fprintf(fd, "%6d   %10g   %10g   %10g       %10g\n", step + 1, ekin, epot, ekin + epot, temp);
      printf("%6d   %10g   %10g   %10g       %10g\n", step + 1, ekin, epot, ekin + epot, temp);
    }

  fclose(fd);

  return 0;
}
