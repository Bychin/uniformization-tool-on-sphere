CC = g++

BOOST_DIR = /usr/local/lib/boost_1_72_0
BOOST_INC = -I$(BOOST_DIR)

BUILD_DIR = ./build
BIN_DIR = ./bin
INC_DIR = ./include

CFLAGS = -I$(INC_DIR) $(BOOST_INC) -std=c++17

SERVER_OBJ = server.o cfg.o angular_gauss.o classic_grid.o spiral_grid.o \
	isoline_api.o

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

test_angular_gauss_distribution: angular_gauss.o test_angular_gauss_distribution.o
	$(CC) $(CFLAGS) $(TEST_CFLAGS) $(BUILD_DIR)/angular_gauss.o $(BUILD_DIR)/test_angular_gauss_distribution.o \
	-o $(BIN_DIR)/test_angular_gauss_distribution

.PHONY: all configure clean server test

# compiled objects

server.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ server.cpp

cfg.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ util/cfg.cpp

angular_gauss.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ distributions/angular_gauss.cpp

classic_grid.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ grids/classic_grid.cpp

spiral_grid.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ grids/spiral_grid.cpp

isoline_api.o:
	$(CC) $(CFLAGS) -I./ -c -o $(BUILD_DIR)/$@ api/isoline.cpp
	# TODO -I./

test_angular_gauss_distribution.o:
	$(CC) $(CFLAGS) $(TEST_CFLAGS) -c -o $(BUILD_DIR)/$@ test/test_angular_gauss_distribution.cpp
