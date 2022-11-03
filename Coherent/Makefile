CC = g++
CFLAGS = -O3 --std=c++11

Coherent.o : Coherent.cpp Coherent.h
	$(CC) -c Coherent.cpp -o Coherent.o $(CFLAGS)

Coherent.exe : Coherent.o ../Random/random.o
	$(CC) Coherent.o ../Random/random.o -o Coherent.exe $(CFLAGS)

clean :
	rm *.o