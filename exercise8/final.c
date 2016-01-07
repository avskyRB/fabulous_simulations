#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double cs = 2.0;
double dif[3];
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
double ccfl;

void init()
{

	N = 10;
	M = 5;
	Lx = 3.0/6.0;
	Ly = 1.5/6.0;
	xspacing = Lx/N;
	yspacing = Ly/M;
	Tmax = 1.5;
	ccfl = 0.4;
}



double * calculate_flux(double uu,double vv,double rrho)
{
	double *fluxx = malloc(3 * sizeof(double));
	fluxx[0] = rrho*uu;
	fluxx[1] = rrho*uu*uu+rrho*cs*cs;
	fluxx[2] = rrho*uu*vv;
	return fluxx;
}

double * riemann_roe_isothermal(double *fluxl, double *fluxr, double rhol, double ul, double vl, double rhor, double ur, double vr)
{
	dif[0] = rhor-rhol;
	dif[1] = rhor*ur-rhol*ul;
	dif[2] = rhor*vr-rhol*vl;
	//PROBLEMA : NEGATIVE DENSITY :(
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
		
	alfa[0] = (((u_medium + cs)*dif[0])-dif[1])/(2*cs);
	alfa[1] = (-((u_medium - cs)*dif[0])+dif[1])/(2*cs);
	alfa[2] = dif[2]-(v_medium*dif[0]);

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
	for(int i=0;i<N;i++)
	{
		for (int j=0;j<M;j++)
		{
			for (int k=0; k<3; k++)
				{
					state[k][i][j][1]=0;
				}
			rho[i][j][1] = 0;
			u[i][j][1] = 0;
			v[i][j][1] = 0;
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
		for (int jj = 0; jj<2; jj++) {left[jj] = cell[jj] + offsetleft[jj]; right[jj] = cell[jj] + offsetright[jj];}
		
		// BOUNDARY CONDITIONS
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

		double *fluxl = calculate_flux(uleft,vleft,rholeft);
		double *fluxr = calculate_flux(u[i][j][0],v[i][j][0],rho[i][j][0]);
		double *fluxL = riemann_roe_isothermal(fluxl,fluxr,rholeft,uleft,vleft,rho[i][j][0],u[i][j][0],v[i][j][0]);
		
		double *fluxl2 = calculate_flux(u[i][j][0],v[i][j][0],rho[i][j][0]);
		double *fluxr2 = calculate_flux(uright,vright,rhoright);
		double *fluxR = riemann_roe_isothermal(fluxl2,fluxr2,rho[i][j][0],u[i][j][0],v[i][j][0],rhoright,uright,vright);

		for (int k=0; k<3; k++)
			{
				state[k][i][j][1] = state[k][i][j][0] + (tstep/xspacing)*(fluxL[k]-fluxR[k]);

			}
		rho[i][j][1] = state[0][i][j][1];
		u[i][j][1] = state[1][i][j][1]/rho[i][j][1];
		v[i][j][1] = state[2][i][j][1]/rho[i][j][1];	
		}
	}


	for(int i=0;i<N;i++)
	{
		for(int j=0; j<M;j++)
		{
			printf("%f ", rho[i][j][1]);
		}
		printf("\n");
	}

}

void calculatetstep(double state[3][N][M][2], double maxu, double maxv)
{
			for(int i=0;i<N;i++)
				{ 
    				for(int j=0;j<M;j++)
    				{ 
        				if((state[1][i][j][0]/state[0][i][j][0])>maxu)
           					 maxu=state[1][i][j][0]/state[0][i][j][0];
           				if((state[2][i][j][0]/state[0][i][j][0])>maxv)
           					 maxv=state[2][i][j][0]/state[0][i][j][0];
    				}
				} 
			if(xspacing>yspacing && maxu>maxv); tstep = ccfl * (yspacing / (cs + maxu ));
			if(xspacing<yspacing && maxu>maxv); tstep = ccfl * (xspacing / (cs + maxu ));
			if(xspacing>yspacing && maxv>maxu); tstep = ccfl * (yspacing / (cs + maxv ));
			if(xspacing<yspacing && maxv>maxu); tstep = ccfl * (xspacing / (cs + maxv ));

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
	double t=0;
	double maxu = 0;
	double maxv = 0;

	//INITIAL DENSITY;
	for(int i=0;i<N;i++)
	{
		for (int j=0;j<M;j++)
		{
			u[i][j][0] = 0;
			v[i][j][0] = 0;
			if( (fabs((i*(Lx/N))-Lx/2) < Lx/4) && (fabs(j*(Ly/M)-Ly/2)<Ly/4) )
				rho[i][j][0] = 4.0;
			else
				rho[i][j][0]=1.0;
			state[0][i][j][0] = rho[i][j][0];
			state[1][i][j][0] = u[i][j][0]*rho[i][j][0];
			state[2][i][j][0] = v[i][j][0]*rho[i][j][0];
		}
	}
	

	while(t<0.4)
		{
			calculatetstep(state, maxu, maxv);
			printf("%f\n",tstep);
			printf("******************************************************************sweep en x\n");	
			sweep(state,rho,u,v,tstep,xspacing,cs,N,M,offsetleftx,offsetrightx);

			for(int i=0;i<N;i++)
				{
					for (int j=0;j<M;j++)
						{
							for(int k=0;k<3;k++)
								{
									state[k][i][j][0] = state[k][i][j][1];
								}
								rho[i][j][0] = rho[i][j][1];
								u[i][j][0] = u[i][j][1];
								v[i][j][0] = v[i][j][1];

								
						}
				}
			printf("********************************************************************sweep en y\n");

			/*for(int i=0; i<N;i++)
			{
				for (int j=0; j<M; j++)
				{
					printf("%f ", state[0][i][j][1]);
				}
				printf("\n");
			}*/
			
			sweep(state,rho,v,u,tstep,yspacing,cs,N,M,offsetlefty,offsetrighty);
		
			for(int i=0;i<N;i++)
				{
					for (int j=0;j<M;j++)
						{
							for(int k=0;k<3;k++)
								{
									state[k][i][j][0] = state[k][i][j][1];
								}
								rho[i][j][0] = rho[i][j][1];
								u[i][j][0] = u[i][j][1];
								v[i][j][0] = v[i][j][1];

								
						}
				}

	
			t = t + tstep;
		}



	FILE *f = fopen("result.txt", "w");
	if (f == NULL)
	{
	    printf("Error opening file!\n");
	    exit(1);
	}

	for(int j = 0; j < M; j++)
	{
		fprintf(f, "%f  %f\n", j*yspacing, rho[N/2 - 1][j][1]);
	}

	fclose(f);


	return 0;
}

