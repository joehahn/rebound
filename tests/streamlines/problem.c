#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"

extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;
int N_init = 5;										// number of particles
int Nsteps_per_output = 100;						// steps per output
int Noutput = 250;						            // number of outputs (plus one more)

void problem_init(int argc, char* argv[]){

	// Integrator parameters
	OMEGA = 1.0;
	dt = 1.0e-3;										// timestep
	tmax = dt*Nsteps_per_output*Noutput;				// simulation stop time
//#ifdef OPENGL
//	display_rotate_z = 90;								// display rotation angle (deg)
//	display_rotate_x = 0;	
//#endif

	//Particle parameters
	double x0_max = 10;									// particles' initial -x0_max < x < x0_max
	double y0_max = 50;									// particles' initial y range
    double z0_max = 0.0;								// particles' initial z range
	boxsize = 1;
	root_nx	= 20;
	root_ny = 100;
	root_nz = 10;
	nghostx = 0;
	nghosty = 1;										// ghost boxes along y direction only
	nghostz = 0;
	double radius = 0.0;								// particle radius
	double density = 0.5;								// particle density in cgs when radius is in cm
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = 0.001;
	init_box();

	// Initial conditions
	for (int i=0; i<N_init; i++){
		struct particle p;
		p.x = tools_uniform(-x0_max, x0_max);
		p.y = tools_uniform(-y0_max, y0_max);
		p.z = tools_uniform(-z0_max, z0_max);
		p.vx = 0;
		p.vy = -1.5*p.x*OMEGA;
		p.vz = 0;
		p.ax = 0;
		p.ay = 0;
		p.az = 0;
		p.m = 4*M_PI*density*radius*radius*radius/3;
		p.r = radius;
	    //pt.id = i;
		particles_add(p);
	}
}

void problem_inloop(){

}

void problem_output(){
	// Screen and file output
	if (output_check(Nsteps_per_output*dt)) output_timing();
	if (output_check(Nsteps_per_output*dt)) output_append_ascii("rebound.txt");
}

void problem_finish(){
    // Output final state
	if (output_check(dt/2)) output_timing();
	if (output_check(dt/2)) output_append_ascii("rebound.txt");
	
	/*FILE* of = fopen("error.txt","a+"); 
	double error= N_init-N;
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	double error_limit = 0.1;
	int pass = (error>error_limit?0:1);
	fprintf(of,"%d\t%e\t%e\t%e\t%d\t%e\t",N,dt,error,error_limit,pass,timing_final-timing_initial);
	fprintf(of,"N_init = %d",N_init);
	fprintf(of,"\n");
	fclose(of);*/
}
