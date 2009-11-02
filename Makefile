all:
	g++ -o fractals -O2 -lglut fractal.cpp

clean:
	-rm *~ fractals