main: toml.o mandelbrot.c
	gcc -o mandelbrot.out mandelbrot.c toml.o -fopenmp -lm
	gcc -o mandelbrot_O3.out mandelbrot.c toml.o -fopenmp -lm -O3

toml.o: toml.c
	gcc -c toml.c -std=c99 -Wall -Wextra -fpic -O2 -DNDEBUG

clean:
	rm -f *.o mandelbrot.out mandelbrot_O3.out