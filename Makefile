CC=g++
CFLAGS=-l xcb

particles: particles.cpp
	$(CC) particles.cpp -o particles $(CFLAGS)

clean:
	rm particles
