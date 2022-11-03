CC = g++
CFLAGS = -O3 --std=c++11

Thermal.o : Thermal.cpp Thermal.h
	$(CC) -c Thermal.cpp -o Thermal.o $(CFLAGS)

Thermal.exe : Thermal.o ../Random/random.o
	$(CC) Thermal.o ../Random/random.o -o Thermal.exe $(CFLAGS)

clean :
	rm *.o Thermal.exe
