#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double p(double x)
{
	return 0.078/(pow(x-2,4)+pow(sin(x-3),8));
}

double randd(double min, double max) 
{
    double range = (max - min); 
	double div = (double)RAND_MAX / range;
    return min + ((double)rand() / div);
}

int main()
{
	FILE *ptr;
	ptr=fopen("rej01.txt","w");
	
	srand(time(NULL));
	
	int i=0;
		
	while(i<1000000)
	{
		double x=randd(0,5),y=randd(0,5);
		
		if(y<=p(x))
		{
			fprintf(ptr,"%f\n",x);
			i++;
		}
	}
	
	fclose(ptr);
	
	return 0;
}
