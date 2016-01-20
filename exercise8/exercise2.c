#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double rhol = 2.0;
double ul = -1.0;
double vl = -2.0;
double rhor= 1.0;
double ur = -1.0;
double vr = 2.0;
double cs = 2.0;
double u[3];
double u_medium;
double v_medium;
double k[3][3]; 
double alfa[3];
double lambda[3];
double flux[3];

void init()
{
		u[0] = rhor-rhol;
		u[1] = rhor*ur-rhol*ul;
		u[2] = rhor*vr-rhol*vl;
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

}

double * calculate_flux(double uu,double vv,double rrho)
{
	double *flux = malloc(3 * sizeof(double));
	flux[0] = rrho*uu;
	flux[1] = rrho*uu*uu+rrho*cs*cs;
	flux[2] = rrho*uu*vv;
	return flux;
}

void riemann_roe_isothermal(double *fluxl, double *fluxr)
{
	for (int i=0; i<3 ; i++)
	{
		for (int j=0; j<3 ; j++)

		{
			flux[i] += (alfa[j]*fabs(lambda[j])*k[i][j]);	
		}
		flux[i] = 0.5*(fluxl[i]+fluxr[i]) - 0.5*flux[i];
	}
	for(int i = 0; i < 3; i++) 
	{
       	printf("%f ", flux[i]);
    }
}


int main()
{
	init();
	double *fluxl = calculate_flux(ul,vl,rhol);
	double *fluxr = calculate_flux(ur,vr,rhor);
	riemann_roe_isothermal(fluxl, fluxr);
	free(fluxl);
	free (fluxr);
	return 0;
}

