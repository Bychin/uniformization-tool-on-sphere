CXXC = g++

BOOST_DIR = /usr/local/lib/boost_1_72_0
BOOST_INC = -I$(BOOST_DIR)

BUILD_DIR = ./build
BIN_DIR = ./bin
INC_DIR = ./include

CXXFLAGS = -I$(INC_DIR) -I./ $(BOOST_INC) -std=c++17 -O2

SERVER_OBJ = server.o cfg.o util.o angular_gauss.o classic_grid.o \
	spiral_grid.o isoline_api.o stats_api.o statistics.o


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
	$(CXXC) $(CXXFLAGS) $(addprefix $(BUILD_DIR)/,$(SERVER_OBJ)) -o $(BIN_DIR)/server


.PHONY: all configure clean server test test_angular_gauss_distribution


# compiled objects

server.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ server.cpp

statistics.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ include/psvt/statistics.cpp

cfg.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ util/cfg.cpp

util.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ util/util.cpp

angular_gauss.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ distributions/angular_gauss.cpp

classic_grid.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ grids/classic_grid.cpp

spiral_grid.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ grids/spiral_grid.cpp

isoline_api.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ api/isoline.cpp

stats_api.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ api/stats.cpp
