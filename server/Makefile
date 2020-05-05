CC = g++

BOOST_DIR = /usr/local/lib/boost_1_72_0
BOOST_INC = -I$(BOOST_DIR)

BUILD_DIR = ./build
BIN_DIR = ./bin
INC_DIR = ./include

CFLAGS = -I$(INC_DIR) $(BOOST_INC) -std=c++17

SERVER_OBJ = server.o cfg.o angular_gauss.o

TEST_SERVER_OBJ = test_angular_gauss_distribution.o
TEST_CFLAGS = -I./

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

test:
	make \
	configure \
	test_angular_gauss_distribution

	$(BIN_DIR)/test_angular_gauss_distribution

# individual executables

test_angular_gauss_distribution: distribution.o test_angular_gauss_distribution.o
	$(CC) $(CFLAGS) $(TEST_CFLAGS) $(BUILD_DIR)/angular_gauss.o $(BUILD_DIR)/test_angular_gauss_distribution.o \
	-o $(BIN_DIR)/test_angular_gauss_distribution

.PHONY: all configure clean server test

# compiled objects

server.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/server.o server.cpp

cfg.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/cfg.o util/cfg.cpp

angular_gauss.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/angular_gauss.o distributions/angular_gauss.cpp

test_angular_gauss_distribution.o:
	$(CC) $(CFLAGS) $(TEST_CFLAGS) -c -o $(BUILD_DIR)/test_angular_gauss_distribution.o test/test_angular_gauss_distribution.cpp
