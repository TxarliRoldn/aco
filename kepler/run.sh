#!/bin/bash

# Set environmental variables
export OMP_NUM_THREADS=$(nproc --all)
export JULIA_NUM_THREADS=$(nproc --all)
export NUMBA_NUM_THREADS=$(nproc --all)

# Write headers
printf "NSTEPS\tt_rk4\tt_leapfrog\terr_t_rk4\terr_t_leapfrog\n" > ccfile
printf "NSTEPS\tt_rk4\tt_leapfrog\terr_t_rk4\terr_t_leapfrog\n" > cofile
printf "NSTEPS\tt_rk4\tt_leapfrog\terr_t_rk4\terr_t_leapfrog\n" > jlfile
printf "NSTEPS\tt_rk4\tt_leapfrog\terr_t_rk4\terr_t_leapfrog\n" > pyfile
printf "NSTEPS\tt_rk4\tt_leapfrog\terr_t_rk4\terr_t_leapfrog\n" > ppfile

for NSTEPS in 1000 2000 5000 10000 20000 50000 100000
do
    echo $NSTEPS

    # Generate input
    sed "s/PATATA/$NSTEPS/g" template.toml > input.toml

    # Launch benchmarks
    ./kepler.out     >> ccfile
    ./kepler_O3.out  >> cofile
    julia kepler.jl  >> jlfile
    python kepler.py >> pyfile
done

for NSTEPS in 1000 2000 5000 10000 20000
do
    echo $NSTEPS

    # Generate input
    sed "s/PATATA/$NSTEPS/g" template.toml > input.toml

    # Launch benchmarks
    python kepler_pure.py >> ppfile
done

rm input.toml