#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "toml.h"

void *safe_malloc(size_t, size_t);
void compute(int (*f)(double, double, long), int**, double*, double*, long, long);
int isconvergent_real(double, double, long);
int isconvergent_complex(double, double, long);
double mean(double *, long);
double var(double *, long);

int main(){
	FILE* fp;
    char errbuf[200];
	fp = fopen("input.toml", "r");
	toml_table_t* params = toml_parse_file(fp, errbuf, sizeof(errbuf));
	fclose(fp);

	long NPTS_DIM   = (long)   toml_int_in(params, "NPTS_DIM").u.i;
	long MAXITER    = (long)   toml_int_in(params, "MAXITER").u.i;
	long NRUNS      = (long)   toml_int_in(params, "NRUNS").u.i;
	double REAL_MIN = (double) toml_double_in(params, "REAL_MIN").u.d;
	double IMAG_MIN = (double) toml_double_in(params, "IMAG_MIN").u.d;
	double REAL_MAX = (double) toml_double_in(params, "REAL_MAX").u.d;
	double IMAG_MAX = (double) toml_double_in(params, "IMAG_MAX").u.d;

	toml_free(params);

	double dc_real = (REAL_MAX - REAL_MIN)/(NPTS_DIM - 1.);
	double dc_imag = (IMAG_MAX - IMAG_MIN)/(NPTS_DIM - 1.);
	
	long 	i;
	struct 	timespec start, end;
	double	*t_real;
	double	*t_complex;
	double 	*c_real;
	double 	*c_imag;
	int 	**converged;

	t_real    = (double *) safe_malloc(NRUNS, sizeof (double));
	t_complex = (double *) safe_malloc(NRUNS, sizeof (double));
	c_real    = (double *) safe_malloc(NPTS_DIM, sizeof (double));
	c_imag    = (double *) safe_malloc(NPTS_DIM, sizeof (double));
	converged = (int **)   safe_malloc(NPTS_DIM, sizeof (int**));
	for (i = 0; i < NPTS_DIM; i++){
		c_real[i]    = REAL_MIN + i*dc_real;
		c_imag[i]    = IMAG_MIN + i*dc_imag;
		converged[i] = (int *) safe_malloc(NPTS_DIM, sizeof (int *));
	}

	for (i = 0; i < NRUNS; i++){
		clock_gettime(CLOCK_MONOTONIC, &start);
		compute(isconvergent_real, converged, c_real, c_imag, NPTS_DIM, MAXITER);
		clock_gettime(CLOCK_MONOTONIC, &end);
		t_real[i] = (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec)/1.e9;

		clock_gettime(CLOCK_MONOTONIC, &start);
		compute(isconvergent_complex, converged, c_real, c_imag, NPTS_DIM, MAXITER);
		clock_gettime(CLOCK_MONOTONIC, &end);
		t_complex[i] = (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec)/1.e9;
	}

	printf("%ld\t%.8e\t%.8e\t%.8e\t%.8e\n", NPTS_DIM, mean(t_real, NRUNS), mean(t_complex, NRUNS), sqrt(var(t_real, NRUNS)/NRUNS), sqrt(var(t_complex, NRUNS)/NRUNS));

	free(t_real);
    free(t_complex);
	free(c_real);
	free(c_imag);
	for( i = 0; i < NPTS_DIM; i++ ){ free(converged[i]); }
	free(converged);

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

void compute(int (*f)(double, double, long), int** converged, double* c_real, double* c_imag, long NPTS_DIM, long MAXITER){
	int i, j;

	#pragma omp parallel for default(none) private(i, j) shared(f, converged, c_real, c_imag, NPTS_DIM, MAXITER) schedule(dynamic)
	for( i = 0; i < NPTS_DIM; i++ ){
		for( j = 0; j < NPTS_DIM; j++ ){ converged[i][j] = f(c_real[i], c_imag[j], MAXITER); }
	}
	return;
}

int isconvergent_real(double c_real, double c_imag, long MAXITER){
	int i;
	double z_real = 0.;
	double z_imag = 0.;
	double zr2, zi2;
	
	for( i = 0; i < MAXITER; i++ ){
		z_real = z_real*z_real - z_imag*z_imag + c_real;
		z_imag = 2*z_real*z_imag + c_imag;
		if( (z_real*z_real + z_imag*z_imag) > 4 ){ break; }
	}
	return i;
}

int isconvergent_complex(double c_real, double c_imag, long MAXITER){
	int i;
	double complex c = c_real + c_imag * I;
	double complex z = 0. + 0. * I;
	
	for( i = 0; i < MAXITER; i++ ){
		z = z*z + c;
		if( creal(z)*creal(z) + cimag(z)*cimag(z) > 4 ){ break; }
	}
	return i;
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