all: main.o lodepng.o
	gcc main.o -o main

lodepng.o: lodepng.c lodepng.h
	gcc -c lodepng.c

clean:
	rm -f *.o
