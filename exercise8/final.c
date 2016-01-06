#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double cs = 2.0;
double u[3];
double u_medium;
double v_medium;
double k[3][3]; 
double alfa[3];
double lambda[3];
int N;
int M;
double tstep;
double xspacing;
double yspacing;
double Lx;
double Ly;
double Tmax;
double n;
void init()
{

	tstep = 0.01;
	N = 60;
	M = 30;
	Lx = 3.0;
	Ly = 1.5;
	xspacing = Lx/N;
	yspacing = Ly/M;
	Tmax = 1.5;
	tstep = 0.01;
}

double * calculate_flux(double uu,double vv,double rrho)
{
	double *flux = malloc(3 * sizeof(double));
	flux[0] = rrho*uu;
	flux[1] = rrho*uu*uu+rrho*cs*cs;
	flux[2] = rrho*uu*vv;
	return flux;
}

double * riemann_roe_isothermal(double *fluxl, double *fluxr, double rhol, double ul, double vl, double rhor, double ur, double vr)
{
	u_medium = (sqrt(rhol)*ul + sqrt(rhor)*ur)/(sqrt(rhol)+sqrt(rhor));
	v_medium = (sqrt(rhol)*vl + sqrt(rhor)*vr)/(sqrt(rhol)+sqrt(rhor));
	k[0][0]=1.0;
	k[1][0] = u_medium - cs;
	k[2][0] = v_medium;

	k[0][1]=1.0;
	k[1][1] = u_medium + cs;
	k[2][1] = v_medium;

	k[0][2]=0.0;
	k[1][2] = 0.0;
	k[2][2] = 1.0;
		
	alfa[0] = (((u_medium + cs)*u[0])-u[1])/(2*cs);
	alfa[1] = (-((u_medium - cs)*u[0])+u[1])/(2*cs);
	alfa[2] = u[2]-(v_medium*u[0]);

	lambda[0] = u_medium - cs;
	lambda[1] = u_medium + cs;
	lambda[2] = u_medium;
	
	double *flux = malloc(3 * sizeof(double));
	for (int i=0; i<3 ; i++)
	{
		for (int j=0; j<3 ; j++)
		{
			flux[i] += (alfa[j]*fabs(lambda[j])*k[i][j]);	
		}
		flux[i] = 0.5*(fluxl[i]+fluxr[i]) - 0.5*flux[i];
	}
    return flux;
}

void sweep(	double state[3][N][M][2], double rho[N][M][2], double u[N][M][2], double v[N][M][2], double tstep, double xspacing, double cs, int N, int M, int offsetleft[2], int offsetright[2])
// ONE TIMESTEP
{
// vector de tres en cada celda por cada timestep,solo 1 timestep
	// initial state??
	
	for(int i=0;i<N;i++)
	{
		for (int j=0;j<M;j++)
		{
	
			rho[i][j][0] = state[0][i][j][0];
			u[i][j][0] = state[1][i][j][0]/rho[i][j][0];
			v[i][j][0] = state[2][i][j][0]/rho[i][j][0];
		}
	}
	for(int i=0;i<N;i++)
	{
		for (int j=0;j<M;j++)
		{
	
		int left[2];
		int right[2];
		int cell[2];
		cell[0] = i;
		cell[1] = j;
		for (int jj = 0; jj<2; jj++) {left[jj] = cell[jj] - offsetleft[jj]; right[jj] = cell[jj] - offsetright[jj];}
		
		if (left[0] == -1) left[0] = N-1;
		if (right[0] == N) right[0] = 0;
		if (left[1] == -1) left[1] = M-1;
		if (right[1] == M) right[1] = 0;
		double uleft = u[left[0]][left[1]][0];
		double uright = u[right[0]][right[1]][0];
		double vleft = v[left[0]][left[1]][0];
		double vright = v[right[0]][right[1]][0];
		double rholeft = rho[left[0]][left[1]][0];
		double rhoright = rho[right[0]][right[1]][0];
		//printf("leftx es: %d \n", left[0]);
		//printf("lefty es: %d \n", left[1]);
		//printf("rholeft es: %f \n", rholeft);
		double *fluxl = calculate_flux(uleft,vleft,rholeft);
		double *fluxr = calculate_flux(u[i][j][0],v[i][j][0],rho[i][j][0]);
		double *fluxR = riemann_roe_isothermal(fluxl,fluxr,rholeft,uleft,vleft,rho[i][j][0],u[i][j][0],v[i][j][0]);
		double *fluxl2 = calculate_flux(u[i][j][0],v[i][j][0],rho[i][j][0]);
		double *fluxr2 = calculate_flux(uright,vright,rhoright);
		double *fluxL = riemann_roe_isothermal(fluxl2,fluxr2,rho[i][j][0],u[i][j][0],v[i][j][0],rhoright,uright,vright);
		//printf("fluxl es: %f \n", fluxl[1]);
		//printf("fluxR es: %f \n", fluxr[1]);

		for (int k=0; k<3; k++)
			{
				state[k][i][j][1] = state[k][i][j][0] + (tstep/xspacing)*(fluxL-fluxR);
			}
		//printf("density next timestep es: %f\n", state[0][i][j][1]);
		rho[i][j][1] = state[0][i][j][1];
		u[i][j][1] = state[1][i][j][1]/rho[i][j][1];
		v[i][j][1] = state[2][i][j][1]/rho[i][j][1];
		}
	}

}


int main()
{
	init();
	double rho[N][M][2];
	double u[N][M][2];
	double v[N][M][2];
	double state[3][N][M][2];
	int offsetrightx[2] = {1,0};
	int offsetleftx[2] = {-1,0};
	int offsetrighty[2] = {0,1};
	int offsetlefty[2] = {0,-1};
	// u and v are fields !!!!????? el maximo de ellos?

/*	for(i=1;i<r;i++)
	{ 
    for(j=1;j<c;j++)
    	{ 
        	if(a[i][j]>t)
            t=a[i][j];
    }
} */
	//if(xspacing>yspacing && u>v) tstep = ccfl * yspacing / (cs + u );
	//elif(xspacing<yspacing && u>v) tstep = ccfl * xspacing / (cs + u );
	//elif(xspacing>yspacing && v>u) tstep = ccfl * yspacing / (cs + v );
	//elif(yspacing>xspacing && v>u) tstep = ccfl * xspacing / (cs + v );

	//INITIAL DENSITY;
	for(int i=0;i<N;i++)
	{
		for (int j=0;j<M;j++)
		{
			u[i][j][0] = 0;
			v[i][j][0] = 0;
			if( (abs((i*(Lx/N))-Lx/2) < Lx/4) && (fabs(j*(Ly/M)-Ly/2)<Ly/4) )
				rho[i][j][0] = 4.0;
			else
				rho[i][j][0]=1.0;
			state[0][i][j][0] = rho[i][j][0];
			state[1][i][j][0] = 0;
			state[2][i][j][0] = 0;
		}
	}


	for(double i=0.0; i<=Tmax ; i=i+tstep)
	{	
		
		sweep(state,rho,u,v,tstep,xspacing,cs,N,M,offsetrightx,offsetleftx);
		for(int k=0;k<3;k++)
		{
		for(int i=0;i<N;i++)
		{
			for (int j=0;j<M;j++)
			{
				state[k][i][j][0] = state[k][i][j][1];

			}
		}
		}
		sweep(state,rho,u,v,tstep,yspacing,cs,N,M,offsetrighty,offsetlefty);
				for(int k=0;k<3;k++)
		{
		for(int i=0;i<N;i++)
		{
			for (int j=0;j<M;j++)
			{
				state[k][i][j][0] = state[k][i][j][1];

			}
		}
		}
	
	}

	for(int i = 0; i < N; i++) 
	{
		for(int j = 0; j < M; j++) 
		{
        printf("%f ", v[i][j][1]);
    	}	
    	printf("\n\r ");
    }





	FILE *f = fopen("result.txt", "w");
	if (f == NULL)
	{
	    printf("Error opening file!\n");
	    exit(1);
	}


	/* print integers and floats */
	for(int j = 0; j < M; j++)
	{
		fprintf(f, "%f  %f\n", j*yspacing, rho[M/2][j][1]);
	}




	fclose(f);


	return 0;
}

