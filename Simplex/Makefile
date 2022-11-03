CC = g++
CFLAGS = -O3 --std=c++11

Simplex.o : Simplex.cpp Simplex.h
	$(CC) -c Simplex.cpp -o Simplex.o $(CFLAGS)

Simplex.exe : Simplex.o ../Random/random.o
	$(CC) Simplex.o ../Random/random.o -o Simplex.exe $(CFLAGS)

clean :
	rm *.o
