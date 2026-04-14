CC=gcc
CC_FLAGS=-g -Wall
CC_LIBS=-lm -lfftw3

# Parallel shared-memory build: make OPENMP=1
ifeq ($(OPENMP),1)
  CC_FLAGS += -fopenmp
  CC_LIBS += -fopenmp
endif

SRC_DIR=src
HDR_DIR=include/
OBJ_DIR=obj

# source and object files
SRC_FILES=$(wildcard $(SRC_DIR)/*.c)
OBJ_FILES=$(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILES))

BIN_FILE=cnavier

all: $(OBJ_DIR) $(BIN_FILE)

$(BIN_FILE): $(OBJ_FILES)
	$(CC) $(CC_FLAGS) $^ -I$(HDR_DIR) -o $@ $(CC_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CC_FLAGS) -c $^ -I$(HDR_DIR) -o $@ $(LFLAGS)

$(OBJ_DIR):
	mkdir $@

clean:
	rm -rf $(BIN_FILE) $(OBJ_DIR) $(TBN_DIR) output/*.vtk
