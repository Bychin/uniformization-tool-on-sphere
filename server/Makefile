CC = g++

BUILD_DIR = ./build
BIN_DIR = ./bin
INC_DIR = ./include

CFLAGS = -I$(INC_DIR) -std=c++17

all:
	make \
	configure \
	server

configure:
	mkdir -p $(BIN_DIR)
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BIN_DIR)
	rm -rf $(BUILD_DIR)

server: server.o cfg.o
	$(CC) $(CFLAGS) $(BUILD_DIR)/server.o $(BUILD_DIR)/cfg.o -o $(BIN_DIR)/server

server.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/server.o server.cpp

cfg.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/cfg.o util/cfg.cpp
