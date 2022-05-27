#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "utility.h"
#include "toml.h"

#define POW2(x) ((x)*(x))
#define POW3(x) ((x)*(x)*(x))
#define ABS(x)  (((x) > 0.0) ? (x) : (-(x)))

struct particle{
	long    id;
	double  m;
	double  x;
	double  v;
	double  rho;
	double  e;
    double  h;
	double  dv;
	double  drho;
	double  de;
};
typedef struct particle particle_t;

void   *safe_malloc(size_t, size_t);
double dkernel_near(double, double);
double dkernel_far(double, double);
particle_t create_particle(long, double, double, double, double, double, long);
particle_t *uniform_system(long, long, double, double, double);
void central_blast(particle_t *, long, double, double);
void update_particle(particle_t *, double);
void get_deriv(particle_t *, particle_t *, long, long);
void integrate_euler(particle_t *, long, long, long, double);
double mean(double *, long);
double var(double *, long);

int main(){
    FILE* fp;
    char errbuf[200];
	fp = fopen("input.toml", "r");
	toml_table_t* params = toml_parse_file(fp, errbuf, sizeof(errbuf));
	fclose(fp);

	long NPART    = (long)   toml_int_in(params, "NPART").u.i;
    long NSTEPS   = (long)   toml_int_in(params, "NSTEPS").u.i;
    long NSPH     = (long)   toml_int_in(params, "NSPH").u.i;
    long NRUNS    = (long)   toml_int_in(params, "NRUNS").u.i;
	double TSTART = (double) toml_double_in(params, "TSTART").u.d;
	double TEND   = (double) toml_double_in(params, "TEND").u.d;

	toml_free(params);

    long i;
    double *t_eu;
    struct 	timespec start, end;
    double dt = (TEND - TSTART) / ((double) NSTEPS - 1.0);
    t_eu = (double *) safe_malloc(NRUNS, sizeof (double));

    for (i = 0; i < NRUNS; i++){
        clock_gettime(CLOCK_MONOTONIC, &start);
        particle_t *system;
        system = uniform_system(NPART, NSPH, 1.0, 0.0, 1.0);
        central_blast(system, NPART, 1.0, 1.0e-5);
        integrate_euler(system, NPART, NSPH, NSTEPS, dt);
        clock_gettime(CLOCK_MONOTONIC, &end);
        free(system);
        t_eu[i] = (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec)/1.e9;
    }
    
    printf("%ld\t%.8e\t%.8e\n", NPART, mean(t_eu, NRUNS), sqrt(var(t_eu, NRUNS)/NRUNS));
    free(t_eu);

    return 0;
}

void *safe_malloc(size_t count, size_t size){
	void *ptr;
	ptr = malloc(count * size);

	if (ptr == NULL) {
    	fprintf(stderr,"[ERROR]: Not enough memory (%ld B)\n", (long) (count * size));
    	exit(EXIT_FAILURE);
  	}
  
	return ptr;
}

double dkernel_near(double r, double invh){
	return 0.6666666666666666*POW2(invh)*(-3.0*(r*invh) + 2.25*POW2(r*invh));
}

double dkernel_far(double r, double invh){
	return -0.5*POW2(invh)*POW2(2.0 - (r*invh));
}

particle_t create_particle(long id, double m, double x, double v, double rho, double e, long NSPH){
	particle_t p;
	p.id  = id;
	p.m   = m;
	p.x   = x;
	p.v   = v;
	p.rho = rho;
	p.e   = e;
	p.h   = 0.0;
	return p;
}

particle_t *uniform_system(long NPART, long NSPH, double m, double xmin, double xmax){
	long   		i;
	double 		dx;
	double 		rho;
	particle_t *system;

	dx     = (xmax - xmin)/((double) NPART - 1.0);
	rho    = m/dx;
	system = (particle_t *) safe_malloc(NPART, sizeof (particle_t));

	#pragma omp parallel for default(none) private(i) shared(NPART, m, xmin, dx, rho, system, NSPH) schedule(static)
	for ( i = 0; i < NPART; i++ ){
		system[i] = create_particle(i, m, xmin + i*dx, 0.0, rho, 0.0, NSPH);
	}

	return system;
}

void central_blast(particle_t *system, long NPART, double e_blast, double e_bg){
	long i;

	#pragma omp parallel for default(none) private(i) shared(NPART, system, e_bg) schedule(static)
	for ( i = 0; i < NPART; i++ ){
		system[i].e = e_bg;
	}
	system[(NPART >> 1) - 1].e = e_blast;
	return;
}

void update_particle(particle_t *p, double dt){
	p -> x   += dt * (p -> v);
	p -> v   += dt * (p -> dv);
	p -> rho += dt * (p -> drho);
	p -> e   += dt * (p -> de);
	return;
}

void get_deriv(particle_t *p, particle_t *system, long NPART, long NSPH){
	
	long    i;
	long    imo;
	long   *indexs;
	double *distances;
	double *signs;
	double  td_sph;
	double  invd_sph;
	double  temp0;
	double  temp1;
	double  temp2;
	double  temp3;

	indexs    = (long *)   safe_malloc(NPART, sizeof (long));
	distances = (double *) safe_malloc(NPART, sizeof (double));
	signs     = (double *) safe_malloc(NPART, sizeof (double));

	for ( i = 0; i < NPART; i++ ){
		distances[i]      = (p -> x) - system[i].x;
		signs[i]          = SIGN(distances[i]);
		distances[i]      = ABS(distances[i]);
	}

	indexx(NPART, distances-1, indexs-1);

	p -> h = distances[indexs[NSPH] - 1];

	td_sph    = 2.0*(p -> h);
	invd_sph  = 1.0/(p -> h);
	temp0     = (p -> e)/(p -> rho);

	p -> dv   = 0.0;
	p -> drho = 0.0;
	p -> de   = 0.0;

	for ( i = 1; i <= NSPH; i++ ){
		imo = indexs[i] - 1;

		temp1 = system[imo].e/system[imo].rho + temp0;
		temp2 = system[imo].v - (p -> v);
		temp3 = signs[imo]*dkernel_near(distances[imo], invd_sph);

		p -> dv   -= system[imo].m*temp1*temp3;
		p -> drho -= system[imo].m*temp2*temp3;
		p -> de   -= system[imo].m*temp1*temp2*temp3;
	}
	for ( i = (NSPH+1); i < NPART; i++ ){
		imo = indexs[i] - 1;
		if ( distances[imo] > td_sph ){ break; }

		temp1 = system[imo].e/system[imo].rho + temp0;
		temp2 = system[imo].v - (p -> v);
		temp3 = signs[imo]*dkernel_far(distances[imo], invd_sph);

		p -> dv   -= system[imo].m*temp1*temp3;
		p -> drho -= system[imo].m*temp2*temp3;
		p -> de   -= system[imo].m*temp1*temp2*temp3;
	}

	p -> dv *= 0.4;
	p -> de *= 0.2;

	free(distances);
	free(signs);
	free(indexs);
	
	return;
}

void integrate_euler(particle_t *system, long NPART, long NSPH, long NSTEPS, double dt){
	long ns;
	long i;
	long j;

	for ( ns = 1; ns <= NSTEPS; ns++ ){

		#pragma omp parallel for default(none) private(i) shared(system, NPART, NSPH, dt) schedule(dynamic)
		for ( i = 0; i < NPART; i++ ){
			get_deriv(&system[i], system, NPART, NSPH);
		}
		#pragma omp parallel for default(none) private(i) shared(system, NPART, NSPH, dt) schedule(static)
		for ( i = 0; i < NPART; i++ ){
			update_particle(&system[i], dt);
		}
	}

	return;
}

double mean(double *arr, long n){
	long i;
	double m = 0.0;
	for (i = 0; i < n; i++){
		m += arr[i];
	}
	return m/n;
}

double var(double *arr, long n){
	long i;
	double s;
	double m = mean(arr, n);
	double v = 0.0;
	for (i = 0; i < n; i++){
		s  = arr[i] - m;
		v += s*s;
	}
	return v/(n-1);
}