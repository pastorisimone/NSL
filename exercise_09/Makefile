CC = g++               
CFLAGS = -Wall -O3 -std=c++11 -I../../RNG

esercizio_09.exe: esercizio_09.o ../ParallelRandomNumberGenerator/random.o
	$(CC) esercizio_09.o ../ParallelRandomNumberGenerator/random.o -o esercizio_09.exe $(CFLAGS)

esercizio_09.o: esercizio_09.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c esercizio_09.cpp -o esercizio_09.o $(CFLAGS)

../ParallelRandomNumberGenerator/random.o: ../ParallelRandomNumberGenerator/random.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c ../ParallelRandomNumberGenerator/random.cpp -o ../ParallelRandomNumberGenerator/random.o $(CFLAGS)

clean:
	rm -f *.o esercizio_09.exe ../ParallelRandomNumberGenerator/*.o 

run: esercizio_09.exe
	./esercizio_09.exe