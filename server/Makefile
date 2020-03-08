CC = g++

BOOST_DIR = /usr/local/lib/boost_1_72_0
BOOST_INC = -I$(BOOST_DIR)

BUILD_DIR = ./build
BIN_DIR = ./bin
INC_DIR = ./include

CFLAGS = -I$(INC_DIR) $(BOOST_INC) -std=c++17

SERVER_OBJ = server.o cfg.o distribution.o

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

server: $(SERVER_OBJ)
	$(CC) $(CFLAGS) $(addprefix $(BUILD_DIR)/,$(SERVER_OBJ)) -o $(BIN_DIR)/server

server.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/server.o server.cpp

cfg.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/cfg.o util/cfg.cpp

distribution.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/distribution.o distribution.cpp
