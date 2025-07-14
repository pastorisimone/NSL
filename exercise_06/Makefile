CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo

SRC_DIR = NSL_SIMULATOR/SOURCE
BIN_DIR = NSL_SIMULATOR
OUTPUT_DIR = OUTPUT
OBJ = $(SRC_DIR)/system.o $(SRC_DIR)/particle.o $(SRC_DIR)/random.o $(SRC_DIR)/NSL_SIMULATOR.o
EXE = $(BIN_DIR)/simulator.exe

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $(EXE) $(AFLAGS)

$(SRC_DIR)/NSL_SIMULATOR.o: $(SRC_DIR)/NSL_SIMULATOR.cpp $(SRC_DIR)/system.h $(SRC_DIR)/particle.h
	$(CC) -c $(SRC_DIR)/NSL_SIMULATOR.cpp -o $(SRC_DIR)/NSL_SIMULATOR.o $(CFLAGS)

$(SRC_DIR)/system.o: $(SRC_DIR)/system.cpp $(SRC_DIR)/system.h $(SRC_DIR)/random.h $(SRC_DIR)/particle.h
	$(CC) -c $(SRC_DIR)/system.cpp -o $(SRC_DIR)/system.o $(CFLAGS)

$(SRC_DIR)/particle.o: $(SRC_DIR)/particle.cpp $(SRC_DIR)/particle.h
	$(CC) -c $(SRC_DIR)/particle.cpp -o $(SRC_DIR)/particle.o $(CFLAGS)

$(SRC_DIR)/random.o: $(SRC_DIR)/random.cpp $(SRC_DIR)/random.h
	$(CC) -c $(SRC_DIR)/random.cpp -o $(SRC_DIR)/random.o $(CFLAGS)

run: all
	cd $(SRC_DIR) && ../../$(EXE)

clean:
	rm -f $(SRC_DIR)/*.o $(EXE)

remove:
	rm -f $(OUTPUT_DIR)/*.*
	rm -f $(OUTPUT_DIR)/CONFIG/*.*