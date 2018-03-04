CC=g++
CFLAGS=-lxcb -std=c++11 -lpthread

particles: particles.cpp
	$(CC) particles.cpp -o particles $(CFLAGS)

clean:
	rm particles
