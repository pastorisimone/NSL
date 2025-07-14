CC = mpic++
CFLAGS = -Wall -O3 -std=c++11 -I../ParallelRandomNumberGenerator
EXEC = esercizio_10.exe
SRCS = esercizio_10.cpp esercizio_09_lib.cpp ../ParallelRandomNumberGenerator/random.cpp
OBJS = esercizio_10.o esercizio_09_lib.o ../ParallelRandomNumberGenerator/random.o

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) -o $(EXEC) $(CFLAGS)

esercizio_10.o: esercizio_10.cpp esercizio_09_lib.h ../ParallelRandomNumberGenerator/random.h
	$(CC) -c esercizio_10.cpp -o esercizio_10.o $(CFLAGS)

esercizio_09_lib.o: esercizio_09_lib.cpp esercizio_09_lib.h
	$(CC) -c esercizio_09_lib.cpp -o esercizio_09_lib.o $(CFLAGS)

../ParallelRandomNumberGenerator/random.o: ../ParallelRandomNumberGenerator/random.cpp ../ParallelRandomNumberGenerator/random.h
	$(CC) -c ../ParallelRandomNumberGenerator/random.cpp -o ../ParallelRandomNumberGenerator/random.o $(CFLAGS)

clean:
	rm -f *.o $(EXEC) ../ParallelRandomNumberGenerator/*.o

run: $(EXEC)
	mpirun -np 8 ./$(EXEC)
