#!/bin/bash

# Set environmental variables
export OMP_NUM_THREADS=$(nproc --all)
export JULIA_NUM_THREADS=$(nproc --all)
export NUMBA_NUM_THREADS=$(nproc --all)

# Write headers
printf "NPART\tt_eu\terr_t_eu\n" > ccfile
printf "NPART\tt_eu\terr_t_eu\n" > cofile
printf "NPART\tt_eu\terr_t_eu\n" > jlfile
printf "NPART\tt_eu\terr_t_eu\n" > pyfile
printf "NPART\tt_eu\terr_t_eu\n" > ppfile

for NPART in 50 100 200 500 1000 2000
do
    echo $NPART

    # Generate input
    sed "s/PATATA/$NPART/g" template.toml > input.toml

    # Launch benchmarks
    ./sph.out     >> ccfile
    ./sph_O3.out  >> cofile
    julia sph.jl  >> jlfile
    python sph.py >> pyfile
done

for NPART in 50 100 200 500
do
    echo $NPART

    # Generate input
    sed "s/PATATA/$NPART/g" template.toml > input.toml

    # Launch benchmarks
    python sph_pure.py >> ppfile
done

rm input.toml