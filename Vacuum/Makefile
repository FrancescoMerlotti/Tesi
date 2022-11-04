CC = g++
CFLAGS = -O3 --std=c++11

Vacuum.o : Vacuum.cpp Vacuum.h
	$(CC) -c Vacuum.cpp -o Vacuum.o $(CFLAGS)

Vacuum.exe : Vacuum.o ../Random/random.o
	$(CC) Vacuum.o ../Random/random.o -o Vacuum.exe $(CFLAGS)

clean :
	rm *.o

output_cl :
	rm graph_vacuum/*.dat
