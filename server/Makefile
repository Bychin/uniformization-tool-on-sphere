CXXC = g++
CC = gcc

FFTW_DIR = /usr/local
FFTW_INC = -I$(FFTW_DIR)/include
FFTW_LIB = -L$(FFTW_DIR)/lib -lfftw3

BOOST_DIR = /usr/local/lib/boost_1_72_0
BOOST_INC = -I$(BOOST_DIR)

S2KIT_DIR = ./include/S2kit/src
S2KIT_INC = -I./include/S2kit/include/

BUILD_DIR = ./build
BIN_DIR = ./bin
INC_DIR = ./include

CFLAGS = $(FFTW_INC) $(S2KIT_INC) -std=c11 -O2 -m64 -fPIC -U__STRICT_ANSI__
LDFLAGS = -lm -m64 -fPIC

CXXFLAGS = -I$(INC_DIR) -I./ $(S2KIT_INC) $(BOOST_INC) -std=c++17 -O2

SERVER_OBJ = server.o cfg.o util.o angular_gauss.o classic_grid.o \
	spiral_grid.o isoline_api.o stats_api.o statistics.o von_mises_fisher.o \
	fisher_kappa.o

# S2Kit objects
DLT_SEMI_SRC = chebyshev_nodes.c pmm.c pml.c cospml.c seminaive.c l2_norms.c \
	vector_funcs.c weights.c
DLT_SEMI_OBJ = $(addsuffix .o,$(basename $(DLT_SEMI_SRC)))
FST_SEMI_OBJ = $(DLT_SEMI_OBJ) naive.o s2kit_util.o
CONV_SEMI_OBJ_MEMO = $(FST_SEMI_OBJ) FST_semi_memo.o


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

server: $(SERVER_OBJ) $(CONV_SEMI_OBJ_MEMO)
	$(CXXC) $(CXXFLAGS) $(addprefix $(BUILD_DIR)/,$(SERVER_OBJ)) \
		$(addprefix $(BUILD_DIR)/,$(CONV_SEMI_OBJ_MEMO)) $(FFTW_LIB) $(LDFLAGS) \
		-o $(BIN_DIR)/server


.PHONY: all configure clean server


# compiled objects

server.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ server.cpp

cfg.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ util/cfg.cpp

util.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ util/util.cpp

angular_gauss.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ distributions/angular_gauss.cpp

von_mises_fisher.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ distributions/von_mises_fisher.cpp

classic_grid.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ grids/classic_grid.cpp

spiral_grid.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ grids/spiral_grid.cpp

isoline_api.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ api/isoline.cpp

stats_api.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ api/stats.cpp


# PSVT compiled objects

fisher_kappa.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ $(INC_DIR)/psvt/fisher_kappa.cpp

statistics.o:
	$(CXXC) $(CXXFLAGS) -c -o $(BUILD_DIR)/$@ $(INC_DIR)/psvt/statistics.cpp


# S2Kit compiled objects

FST_semi_fly.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/FST_semi_fly.c

FST_semi_memo.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/FST_semi_memo.c

cospml.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/legendre_polynomials/cospml.c

pml.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/legendre_polynomials/pml.c

pmm.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/legendre_polynomials/pmm.c

l2_norms.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/legendre_polynomials/util/l2_norms.c

vector_funcs.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/legendre_polynomials/util/vector_funcs.c

naive.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/legendre_transform/naive.c

seminaive.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/legendre_transform/seminaive.c

weights.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/legendre_transform/weights.c

chebyshev_nodes.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/util/chebyshev_nodes.c

s2kit_util.o:
	$(CC) $(CFLAGS) -c -o $(BUILD_DIR)/$@ $(S2KIT_DIR)/util/util.c
