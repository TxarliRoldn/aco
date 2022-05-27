#!/bin/bash

# Set environmental variables
export OMP_NUM_THREADS=$(nproc --all)
export JULIA_NUM_THREADS=$(nproc --all)
export NUMBA_NUM_THREADS=$(nproc --all)

# Write headers
printf "NPTS_DIM\tt_real\tt_complex\terr_t_real\terr_t_complex\n" > ccfile
printf "NPTS_DIM\tt_real\tt_complex\terr_t_real\terr_t_complex\n" > cofile
printf "NPTS_DIM\tt_real\tt_complex\terr_t_real\terr_t_complex\n" > jlfile
printf "NPTS_DIM\tt_real\tt_complex\terr_t_real\terr_t_complex\n" > pyfile
printf "NPTS_DIM\tt_real\tt_complex\terr_t_real\terr_t_complex\n" > ppfile

for NPTS_DIM in 10 20 50 100 200 500 1000 5000
do
    echo $NPTS_DIM

#     # Generate input
    sed "s/PATATA/$NPTS_DIM/g" template.toml > input.toml

#     # Launch benchmarks
    ./mandelbrot.out     >> ccfile
    ./mandelbrot_O3.out  >> cofile
    julia mandelbrot.jl  >> jlfile
    python mandelbrot.py >> pyfile
done

for NPTS_DIM in 10 20 50 100 200 500
do
    echo $NPTS_DIM

    # Generate input
    sed "s/PATATA/$NPTS_DIM/g" template.toml > input.toml

    # Launch benchmarks
    python mandelbrot_pure.py >> ppfile
done

rm input.toml