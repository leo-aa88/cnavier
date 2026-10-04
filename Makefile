CC=gcc
CC_FLAGS=-g -Wall
CC_LIBS=-lm -lfftw3

# Parallel shared-memory build: make OPENMP=1
ifeq ($(OPENMP),1)
  CC_FLAGS += -fopenmp
  CC_LIBS += -fopenmp
endif

# GPU build: make CUDA=1 (needs the NVIDIA CUDA toolkit). nvcc is taken from
# PATH, falling back to /usr/local/cuda; set NVCC or CUDA_HOME to override.
ifeq ($(CUDA),1)
  NVCC ?= $(or $(shell command -v nvcc 2>/dev/null),/usr/local/cuda/bin/nvcc)
  CUDA_HOME ?= $(patsubst %/bin/nvcc,%,$(NVCC))
  NVCC_FLAGS ?= -O2
  CC_FLAGS += -DUSE_CUDA
  ifneq ($(wildcard $(CUDA_HOME)/lib64/libcudart*),)
    CC_LIBS += -L$(CUDA_HOME)/lib64 -Wl,-rpath,$(CUDA_HOME)/lib64
  endif
  CC_LIBS += -lcufft -lcudart -lstdc++
endif

SRC_DIR=src
HDR_DIR=include/
OBJ_DIR=obj
TEST_DIR=tests

# source and object files
SRC_FILES=$(wildcard $(SRC_DIR)/*.c)
OBJ_FILES=$(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILES))
ifeq ($(CUDA),1)
  OBJ_FILES += $(patsubst $(SRC_DIR)/%.cu, $(OBJ_DIR)/%.o, $(wildcard $(SRC_DIR)/*.cu))
endif

BIN_FILE=cnavier
TEST_BIN=test_cnavier

# Records the build flags so that switching OPENMP/CUDA rebuilds every object
CONFIG=$(OBJ_DIR)/config

all: $(OBJ_DIR) $(BIN_FILE)

$(BIN_FILE): $(OBJ_FILES)
	$(CC) $(CC_FLAGS) $^ -I$(HDR_DIR) -o $@ $(CC_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(CONFIG)
	$(CC) $(CC_FLAGS) -c $< -I$(HDR_DIR) -o $@ $(LFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu $(CONFIG)
	$(NVCC) $(NVCC_FLAGS) -c $< -I$(HDR_DIR) -o $@

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.c $(CONFIG)
	$(CC) $(CC_FLAGS) -c $< -I$(HDR_DIR) -o $@

$(OBJ_DIR):
	mkdir $@

$(CONFIG): FORCE | $(OBJ_DIR)
	@echo '$(CC) $(CC_FLAGS) $(NVCC) $(NVCC_FLAGS)' | cmp -s - $@ || echo '$(CC) $(CC_FLAGS) $(NVCC) $(NVCC_FLAGS)' > $@

# Tests: make test (add CUDA=1 to also check the GPU backend against the CPU)
$(TEST_BIN): $(OBJ_DIR)/test_solver.o $(filter-out $(OBJ_DIR)/main.o, $(OBJ_FILES))
	$(CC) $(CC_FLAGS) $^ -I$(HDR_DIR) -o $@ $(CC_LIBS)

test: $(OBJ_DIR) $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -rf $(BIN_FILE) $(TEST_BIN) $(OBJ_DIR) $(TBN_DIR) output/*.vtk

FORCE:

.PHONY: all test clean FORCE
