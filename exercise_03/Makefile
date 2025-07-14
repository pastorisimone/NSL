CC = g++               
CFLAGS = -Wall -O3 -std=c++11 -I../../RNG

esercizio_03_1.exe: esercizio_03_1.o ../ParallelRandomNumberGenerator/random.o
	$(CC) esercizio_03_1.o ../ParallelRandomNumberGenerator/random.o -o esercizio_03_1.exe $(CFLAGS)

esercizio_03_1.o: esercizio_03_1.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c esercizio_03_1.cpp -o esercizio_03_1.o $(CFLAGS)

../ParallelRandomNumberGenerator/random.o: ../ParallelRandomNumberGenerator/random.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c ../ParallelRandomNumberGenerator/random.cpp -o ../ParallelRandomNumberGenerator/random.o $(CFLAGS)

clean:
	rm -f *.o esercizio_03_1.exe ../ParallelRandomNumberGenerator/*.o 

run: esercizio_03_1.exe
	./esercizio_03_1.exe