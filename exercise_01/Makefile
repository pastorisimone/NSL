CC = g++               
CFLAGS = -Wall -O3 -std=c++11 -I../../RNG

all: esercizio_01_1.exe esercizio_01_2.exe esercizio_01_3.exe

esercizio_01_1.exe: esercizio_01_1.o ../ParallelRandomNumberGenerator/random.o
	$(CC) esercizio_01_1.o ../ParallelRandomNumberGenerator/random.o -o esercizio_01_1.exe $(CFLAGS)

esercizio_01_2.exe: esercizio_01_2.o ../ParallelRandomNumberGenerator/random.o
	$(CC) esercizio_01_2.o ../ParallelRandomNumberGenerator/random.o -o esercizio_01_2.exe $(CFLAGS)

esercizio_01_3.exe: esercizio_01_3.o ../ParallelRandomNumberGenerator/random.o
	$(CC) esercizio_01_3.o ../ParallelRandomNumberGenerator/random.o -o esercizio_01_3.exe $(CFLAGS)

esercizio_01_1.o: esercizio_01_1.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c esercizio_01_1.cpp -o esercizio_01_1.o $(CFLAGS)

esercizio_01_2.o: esercizio_01_2.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c esercizio_01_2.cpp -o esercizio_01_2.o $(CFLAGS)

esercizio_01_3.o: esercizio_01_3.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c esercizio_01_3.cpp -o esercizio_01_3.o $(CFLAGS)

../ParallelRandomNumberGenerator/random.o: ../ParallelRandomNumberGenerator/random.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c ../ParallelRandomNumberGenerator/random.cpp -o ../ParallelRandomNumberGenerator/random.o $(CFLAGS)

clean:
	rm -f *.o esercizio_01_1.exe esercizio_01_2.exe esercizio_01_3.exe ../ParallelRandomNumberGenerator/*.o 

run_1: esercizio_01_1.exe
	./esercizio_01_1.exe

run_2: esercizio_01_2.exe
	./esercizio_01_2.exe

run_3: esercizio_01_3.exe
	./esercizio_01_3.exe