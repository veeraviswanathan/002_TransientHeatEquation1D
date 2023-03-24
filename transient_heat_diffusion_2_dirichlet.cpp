#include<stdio.h>
#include<stdlib.h>
#include<array>
#include<cmath>


double L2norm(double *old, double *newarr);
int main()
{
	//declare domain variables
	double end_time=200 ; // seconds
	double time_now=0.0;
	double dt=1e-5; // time step
	int n_nodes = 20; // no of grid points
	double length = 0.05;// length of domain (m)
	double Tl = 600, Tr = 273; //left and right boundary condition (K)
	double Tinit = 300; // initial temperature for domain
	double k = 2e-3; // Diffusivity (m2/s)
	double* T = nullptr; // Temperature array
	double* Told = nullptr; // Temperature array old
	double* S = nullptr; // Source array
	double residual = 1e20;
	double delta_x = (double)length / n_nodes;
	double rhs_coefficient=0.0;
	double Fr, Fl;
	FILE *fp;

	fp=fopen("user_heat_diffusion_transient_explicit.out","w");
	fprintf(fp,"%s \t ","Time");
	for (int kk=0;kk<n_nodes;kk++)
	{
		fprintf(fp,"%s[%d]\t","T",kk);
	}
	fprintf(fp,"\n");

	fclose(fp);


	//memory allocation
	T 		= (double*)calloc(n_nodes, sizeof(double));
	Told 	= (double*)calloc(n_nodes, sizeof(double));
	S 		= (double*)calloc(n_nodes, sizeof(double));

	// declare index variables
	int i;
	int iter = 0;

	// initialisation block

	rhs_coefficient=k*dt/(delta_x*delta_x);
	for (i = 0;i < n_nodes;i++)
	{
		T[i] = Tinit;
		Told[i] = T[i];
		S[i] = 0 * (dt);
	}
	//S[int(n_nodes/2)]= 5000 * (delta_x * delta_x / k);	
	
	//Gauss Siedel iterative method
	fp=fopen("user_heat_diffusion_transient_explicit.out","a+");
	while(time_now<end_time)
	{
		iter++;
		time_now=time_now+iter*dt;

		// Assign old array
		for (i = 0;i < n_nodes;i++)
		{			
			Told[i] = T[i];
			printf("T: %e ",T[i]);
		}

		T[0] = rhs_coefficient*(Told[1]+2.0*Tl)+(1.0-3.0*rhs_coefficient)*Told[0]+dt*S[0];	
		for (i = 1;i < n_nodes - 1;i++)
		{
			T[i] = rhs_coefficient*(Told[i+1]+T[i-1])+(1.0-2.0*rhs_coefficient)*Told[i]+dt*S[i];
			
		}
		printf("\n");
		T[n_nodes - 1] = rhs_coefficient*(T[n_nodes-2]+2.0*Tr)+(1.0-3.0*rhs_coefficient)*Told[n_nodes - 1]+dt*S[n_nodes - 1];

		//Flux at Right boundary
		Fr = (Tr - T[n_nodes-1]) / delta_x;
		//Flux at left boundary
		Fl = (Tl - T[0]) / delta_x;


		//Write output file
		//if(fmod(time_now,2.0)<5e-3)
		if(iter%50==0)
		{
			fprintf(fp,"%e \t",time_now);
			for (i = 0;i < n_nodes;i++)
				{
					fprintf(fp,"%e \t",T[i]);
				}
			fprintf(fp,"\n");
		}	
		//printf("Iteration: %d Time: %e Temperature midpoint: %e Flux Right: %e Flux left: %e\n", iter,time_now, T[int(n_nodes/2)], Fr,Fl);	
		//printf("rhscoeff: %e Iteration: %d Time: %e Temperature midpoint: %e Flux Right: %e Flux left: %e\n",rhs_coefficient, iter,time_now, T[0], Fr,Fl);	
		
	}
	fclose(fp);
	return(0);
}