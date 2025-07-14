CC = g++               
CFLAGS = -Wall -O3 -std=c++11 -I../../RNG

esercizio_05.exe: esercizio_05.o ../ParallelRandomNumberGenerator/random.o
	$(CC) esercizio_05.o ../ParallelRandomNumberGenerator/random.o -o esercizio_05.exe $(CFLAGS)

esercizio_05.o: esercizio_05.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c esercizio_05.cpp -o esercizio_05.o $(CFLAGS)

../ParallelRandomNumberGenerator/random.o: ../ParallelRandomNumberGenerator/random.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c ../ParallelRandomNumberGenerator/random.cpp -o ../ParallelRandomNumberGenerator/random.o $(CFLAGS)

clean:
	rm -f *.o esercizio_05.exe ../ParallelRandomNumberGenerator/*.o 

run: esercizio_05.exe
	./esercizio_05.exe