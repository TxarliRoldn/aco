main: toml.o kepler.c
	gcc -o kepler.out kepler.c toml.o -lm
	gcc -o kepler_O3.out kepler.c toml.o -lm -O3

toml.o: toml.c
	gcc -c toml.c -std=c99 -Wall -Wextra -fpic -O2 -DNDEBUG

clean:
	rm -f *.o kepler.out kepler_O3.out