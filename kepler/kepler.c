#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "toml.h"

#define POW3(x) ((x)*(x)*(x))
#define mGRAV -39.48009132 // AU^3 * yr^{-2} * M_{sun}^{-1}

struct phasespace{
    double x;
    double y;
    double vx;
    double vy;
};
typedef struct phasespace phsp_t;

void   *safe_malloc(size_t, size_t);
phsp_t create_phsp(double, double, double, double);
phsp_t gravity(phsp_t);
phsp_t *integrate(phsp_t (*f)(phsp_t), void (*integrator)(phsp_t (*f)(phsp_t), phsp_t *, long, double), phsp_t, double, long);
void   rk4(phsp_t (*f)(phsp_t), phsp_t *, long, double);
void   leapfrog(phsp_t (*f)(phsp_t), phsp_t *, long, double);
double mean(double *, long);
double var(double *, long);

int main(){
    FILE* fp;
    char errbuf[200];
	fp = fopen("input.toml", "r");
	toml_table_t* params = toml_parse_file(fp, errbuf, sizeof(errbuf));
	fclose(fp);

	long NSTEPS  = (long)   toml_int_in(params, "NSTEPS").u.i;
    long NRUNS   = (long)   toml_int_in(params, "NRUNS").u.i;
	double TIME0 = (double) toml_double_in(params, "TIME0").u.d;
	double TIME1 = (double) toml_double_in(params, "TIME1").u.d;
	double X0    = (double) toml_double_in(params, "X0").u.d;
	double Y0    = (double) toml_double_in(params, "Y0").u.d;
    double VX0   = (double) toml_double_in(params, "X0").u.d;
	double VY0   = (double) toml_double_in(params, "Y0").u.d;

	toml_free(params);

    long i;
    double *t_rk;
	double *t_lp;
    struct timespec start, end;
    phsp_t initial_cond;
    double dt = (TIME1 - TIME0)/((double) NSTEPS - 1.0);
    t_rk = (double *) safe_malloc(NRUNS, sizeof (double));
	t_lp = (double *) safe_malloc(NRUNS, sizeof (double));
    initial_cond= create_phsp(X0, Y0, VX0, VY0);

    for (i = 0; i < NRUNS; i++){
        clock_gettime(CLOCK_MONOTONIC, &start);
        phsp_t *rk_traj = integrate(gravity, rk4, initial_cond, dt, NSTEPS);
        clock_gettime(CLOCK_MONOTONIC, &end);
        t_rk[i] = (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec)/1.e9;
        free(rk_traj);

        clock_gettime(CLOCK_MONOTONIC, &start);
        phsp_t *lp_traj = integrate(gravity, leapfrog, initial_cond, dt, NSTEPS);
        clock_gettime(CLOCK_MONOTONIC, &end);
        t_lp[i] = (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec)/1.e9;
        free(lp_traj);
    }

    printf("%ld\t%.8e\t%.8e\t%.8e\t%.8e\n", NSTEPS, mean(t_rk, NRUNS), mean(t_lp, NRUNS), sqrt(var(t_rk, NRUNS)/NRUNS), sqrt(var(t_lp, NRUNS)/NRUNS));

    free(t_rk);
    free(t_lp);

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

phsp_t create_phsp(double x, double y, double vx, double vy){
    phsp_t phsp;
    phsp.x  = x;
    phsp.y  = y;
    phsp.vx = vx;
    phsp.vy = vy;
    return phsp;
}

phsp_t gravity(phsp_t actual){
    phsp_t derivative;
    double grmt = mGRAV/POW3(sqrt(actual.x*actual.x + actual.y*actual.y));
    derivative = create_phsp(actual.vx, actual.vy, grmt*actual.x, grmt*actual.y);
    return derivative;
}

phsp_t *integrate(phsp_t (*f)(phsp_t), void (*integrator)(phsp_t (*f)(phsp_t), phsp_t *, long, double), phsp_t initial_cond, double dt, long NSTEPS){
    long i;
    phsp_t *solution;

    solution = (phsp_t *) safe_malloc(NSTEPS, sizeof (phsp_t));
    solution[0].x  = initial_cond.x;
    solution[0].y  = initial_cond.y;
    solution[0].vx = initial_cond.vx;
    solution[0].vy = initial_cond.vy;

    for( i = 1; i < NSTEPS; i++ ){
        integrator(f, solution, i, dt);
    }
    return solution;
}

void rk4(phsp_t (*f)(phsp_t), phsp_t *solution, long i, double dt){
    int im1 = i-1;
    double hdt = dt/2.;
    double sdt = dt/6.;
    phsp_t k1;
    phsp_t k2;
    phsp_t k3;
    phsp_t k4;

    k1 = (*f)(solution[im1]);
    solution[i].x  = solution[im1].x  + hdt*k1.x;
    solution[i].y  = solution[im1].y  + hdt*k1.y;
    solution[i].vx = solution[im1].vx + hdt*k1.vx;
    solution[i].vy = solution[im1].vy + hdt*k1.vy;
    k2 = (*f)(solution[i]);
    solution[i].x  = solution[im1].x  + hdt*k2.x;
    solution[i].y  = solution[im1].y  + hdt*k2.y;
    solution[i].vx = solution[im1].vx + hdt*k2.vx;
    solution[i].vy = solution[im1].vy + hdt*k2.vy;
    k3 = (*f)(solution[i]);
    solution[i].x  = solution[im1].x  + dt*k3.x;
    solution[i].y  = solution[im1].y  + dt*k3.y;
    solution[i].vx = solution[im1].vx + dt*k3.vx;
    solution[i].vy = solution[im1].vy + dt*k3.vy;
    k4 = (*f)(solution[i]);
    solution[i].x  = solution[im1].x  + sdt*(k1.x  + 2*k2.x  + 2*k3.x  + k4.x);
    solution[i].y  = solution[im1].y  + sdt*(k1.y  + 2*k2.y  + 2*k3.y  + k4.y);
    solution[i].vx = solution[im1].vx + sdt*(k1.vx + 2*k2.vx + 2*k3.vx + k4.vx);
    solution[i].vy = solution[im1].vy + sdt*(k1.vy + 2*k2.vy + 2*k3.vy + k4.vy);

    return;
}

void leapfrog(phsp_t (*f)(phsp_t), phsp_t *solution, long i, double dt){
    int im1 = i-1;
    double hdt = dt/2.;
    phsp_t half;

    half = (*f)(solution[im1]);
    solution[i].vx = solution[im1].vx + hdt*half.vx;
    solution[i].vy = solution[im1].vy + hdt*half.vy;
    solution[i].x  = solution[im1].x  +  dt*solution[i].vx;
    solution[i].y  = solution[im1].y  +  dt*solution[i].vy;
    half = (*f)(solution[i]);
    solution[i].vx = solution[i].vx + hdt*half.vx;
    solution[i].vy = solution[i].vy + hdt*half.vy;
    
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