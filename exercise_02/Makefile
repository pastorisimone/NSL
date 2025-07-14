CC = g++               
CFLAGS = -Wall -O3 -std=c++11 -I../../RNG

all: esercizio_02_1.exe esercizio_02_2.exe

esercizio_02_1.exe: esercizio_02_1.o ../ParallelRandomNumberGenerator/random.o
	$(CC) esercizio_02_1.o ../ParallelRandomNumberGenerator/random.o -o esercizio_02_1.exe $(CFLAGS)

esercizio_02_2.exe: esercizio_02_2.o ../ParallelRandomNumberGenerator/random.o
	$(CC) esercizio_02_2.o ../ParallelRandomNumberGenerator/random.o -o esercizio_02_2.exe $(CFLAGS)

esercizio_02_1.o: esercizio_02_1.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c esercizio_02_1.cpp -o esercizio_02_1.o $(CFLAGS)

esercizio_02_2.o: esercizio_02_2.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c esercizio_02_2.cpp -o esercizio_02_2.o $(CFLAGS)

../ParallelRandomNumberGenerator/random.o: ../ParallelRandomNumberGenerator/random.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c ../ParallelRandomNumberGenerator/random.cpp -o ../ParallelRandomNumberGenerator/random.o $(CFLAGS)

clean:
	rm -f *.o esercizio_02_1.exe esercizio_02_2.exe ../ParallelRandomNumberGenerator/*.o 

run_1: esercizio_02_1.exe
	./esercizio_02_1.exe

run_2: esercizio_02_2.exe
	./esercizio_02_2.exe