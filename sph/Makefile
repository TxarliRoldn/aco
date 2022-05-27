main: toml.o utility.o sph.c
	gcc -o sph.out sph.c utility.o toml.o -fopenmp -lm
	gcc -o sph_O3.out sph.c utility.o toml.o -fopenmp -lm -O3

toml.o: toml.c
	gcc -c toml.c -std=c99 -Wall -Wextra -fpic -O2 -DNDEBUG

utility.o: utility.c
	gcc -c utility.c -O3 -lm

clean:
	rm -f *.o sph.out sph_O3.out