CC = g++
CFLAGS = -O3 --std=c++11

Optimization.o : Optimization.cpp Optimization.h
	$(CC) -c Optimization.cpp -o Optimization.o $(CFLAGS)
Optimization.exe : Optimization.o ../../Random/random.o
	$(CC) Optimization.o ../../Random/random.o -o Optimization.exe $(CFLAGS)

clean :
	rm *.o Optimization.exe
