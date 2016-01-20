#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double randd(double min, double max) 
{
    double range = (max - min); 
	double div = (double)RAND_MAX / range;
    return min + ((double)rand() / div);
}

double p(double x)
{
	return exp(-pow(x+2*cos(x)*cos(x),2.0));
}

int main()
{
	FILE *ptr;

	ptr=fopen("metropolis.txt","w");
	
	int i;
	double x0=0,x1,p0,p1,r;

	srand(time(NULL));
	
	for(i=0;i<1000000;i++)
	{
		x0=randd(-4,2);
		x1=x0+randd(-1,1);
	
		if(p(x1)/p(x0)>1)
			fprintf(ptr,"%f\n",p(x1));
		else if(p(x1)/p(x0)<=1)
			fprintf(ptr,"%f\n",p(x0));
	}

	fclose(ptr);
	
	return 0;
}
		
