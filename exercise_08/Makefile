CC = g++               
CFLAGS = -Wall -O3 -std=c++11 -I../../RNG

esercizio_08.exe: esercizio_08.o ../ParallelRandomNumberGenerator/random.o
	$(CC) esercizio_08.o ../ParallelRandomNumberGenerator/random.o -o esercizio_08.exe $(CFLAGS)

esercizio_08.o: esercizio_08.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c esercizio_08.cpp -o esercizio_08.o $(CFLAGS)

../ParallelRandomNumberGenerator/random.o: ../ParallelRandomNumberGenerator/random.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c ../ParallelRandomNumberGenerator/random.cpp -o ../ParallelRandomNumberGenerator/random.o $(CFLAGS)

clean:
	rm -f *.o esercizio_08.exe ../ParallelRandomNumberGenerator/*.o 

run: esercizio_08.exe
	./esercizio_08.exe