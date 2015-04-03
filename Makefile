CC=gcc
OPTS=-g -Wall -O0
INCLUDE_DIR=/usr/local/include
LIBRARY_DIR=/usr/local/lib
LIBS=-lgsl -lgslcblas -lpthread -lncurses -lm
INCLUDE_DIR_LOCAL=./include
EXE_DIR_LOCAL=./examples
SRC_DIR_LOCAL=./src

EXE_N=main_HR main_HR_LowMemory main_ROSSLER main_ROSSLER_LowMemory main_LEECH_LowMemory main_HH_Channel_Noise_LowMemory main_HH_Channel_Noise_LowMemory_MT main_HH_Channel_Noise_Orio_LowMemory main_Leech_Channel_Noise_Orio_LowMemory main_Leech_Channel_Noise_Orio_LowMemory_MT

EXE=$(patsubst %,$(EXE_DIR_LOCAL)/%,$(EXE_N))

OBJ_N=Solver.o Stack.o HR.o Rossler.o Leech.o HH_Channel_Noise.o HH_Channel_Noise_Orio.o Leech_Channel_Noise_Orio.o
OBJ=$(patsubst %,$(SRC_DIR_LOCAL)/%,$(OBJ_N))

all: $(EXE)

$(EXE_DIR_LOCAL)/%: $(EXE_DIR_LOCAL)/%.c $(OBJ)
	$(CC) $(OPTS) -o $@ $^ -I$(INCLUDE_DIR) -I$(INCLUDE_DIR_LOCAL) -L$(LIBRARY_DIR) $(LIBS)

$(SRC_DIR_LOCAL)/%.o: $(SRC_DIR_LOCAL)/%.c $(INCLUDE_DIR_LOCAL)/%.h
	$(CC) $(OPTS) -o $@ -c $< -I$(INCLUDE_DIR) -I$(INCLUDE_DIR_LOCAL) -L$(LIBRARY_DIR) $(LIBS)
