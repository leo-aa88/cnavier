CC=gcc
CC_FLAGS=-g -O2 -Wall
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

# Stricter warnings, as in CI: make WERROR=1
ifeq ($(WERROR),1)
  CC_FLAGS += -Wextra -Werror
  NVCC_WERROR = -Werror all-warnings -Xcompiler -Wall,-Wextra,-Werror
endif

# AddressSanitizer and UndefinedBehaviorSanitizer: make SANITIZE=1 (CPU builds)
ifeq ($(SANITIZE),1)
  CC_FLAGS += -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer
  CC_LIBS += -fsanitize=address,undefined
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
CONV_BIN=convergence_study

# Records the build flags so that switching OPENMP/CUDA rebuilds every object
CONFIG=$(OBJ_DIR)/config

all: $(OBJ_DIR) $(BIN_FILE)

$(BIN_FILE): $(OBJ_FILES)
	$(CC) $(CC_FLAGS) $^ -I$(HDR_DIR) -o $@ $(CC_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(CONFIG)
	$(CC) $(CC_FLAGS) -c $< -I$(HDR_DIR) -o $@ $(LFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu $(CONFIG)
	$(NVCC) $(NVCC_FLAGS) $(NVCC_WERROR) -c $< -I$(HDR_DIR) -o $@

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.c $(CONFIG)
	$(CC) $(CC_FLAGS) -c $< -I$(HDR_DIR) -o $@

$(OBJ_DIR):
	mkdir $@

$(CONFIG): FORCE | $(OBJ_DIR)
	@echo '$(CC) $(CC_FLAGS) $(NVCC) $(NVCC_FLAGS) $(NVCC_WERROR)' | cmp -s - $@ || echo '$(CC) $(CC_FLAGS) $(NVCC) $(NVCC_FLAGS) $(NVCC_WERROR)' > $@

# Tests: make test (add CUDA=1 to also check the GPU backend against the CPU)
$(TEST_BIN): $(OBJ_DIR)/test_solver.o $(OBJ_DIR)/mms.o $(filter-out $(OBJ_DIR)/main.o, $(OBJ_FILES))
	$(CC) $(CC_FLAGS) $^ -I$(HDR_DIR) -o $@ $(CC_LIBS)

test: $(OBJ_DIR) $(TEST_BIN)
	./$(TEST_BIN)

# Spatial convergence study against a manufactured solution (a few minutes;
# add OPENMP=1 to make it faster)
$(CONV_BIN): $(OBJ_DIR)/convergence.o $(OBJ_DIR)/mms.o $(filter-out $(OBJ_DIR)/main.o, $(OBJ_FILES))
	$(CC) $(CC_FLAGS) $^ -I$(HDR_DIR) -o $@ $(CC_LIBS)

convergence: $(OBJ_DIR) $(CONV_BIN)
	./$(CONV_BIN)

# ---------------------------------------------------------------------------
# Checks. CI runs each of these as one step (.github/workflows/ci.yml).
# ---------------------------------------------------------------------------

CLANG_FORMAT ?= clang-format
CLANG_TIDY   ?= clang-tidy
CPPCHECK     ?= cppcheck
VALGRIND     ?= valgrind

FORMAT_FILES = $(wildcard $(SRC_DIR)/*.c $(SRC_DIR)/*.cu $(HDR_DIR)*.h $(TEST_DIR)/*.c)
LINT_FILES   = $(wildcard $(SRC_DIR)/*.c $(TEST_DIR)/*.c)

# A short solver run, in a scratch directory so output/ is left alone
SMOKE_RUN = $(TEST_DIR)/run_in_tmp.sh
SMOKE_ARGS = --tf 0.05 --output-interval 5

VALGRIND_FLAGS = --leak-check=full --show-leak-kinds=all --errors-for-leak-kinds=all \
                 --error-exitcode=1 --quiet

# Rewrite the sources in the style of .clang-format
format:
	$(CLANG_FORMAT) -i $(FORMAT_FILES)

format-check:
	$(CLANG_FORMAT) --dry-run --Werror $(FORMAT_FILES)

cppcheck:
	$(CPPCHECK) --enable=warning,performance,portability --std=c11 --error-exitcode=1 \
	            --inline-suppr --quiet -I$(HDR_DIR) $(LINT_FILES)

# Checks are listed in .clang-tidy. Run twice so that the code behind
# #ifdef USE_CUDA and #ifdef _OPENMP is analysed too. cppcheck explores those
# configurations by itself. Neither tool reads the CUDA source; nvcc's own
# warnings are errors with WERROR=1.
tidy:
	$(CLANG_TIDY) --quiet $(LINT_FILES) -- -I$(HDR_DIR) -std=gnu11
	$(CLANG_TIDY) --quiet $(LINT_FILES) -- -I$(HDR_DIR) -std=gnu11 -DUSE_CUDA -fopenmp

# Test suite and a short solver run under ASan + UBSan
test-asan:
	$(MAKE) SANITIZE=1 $(BIN_FILE) $(TEST_BIN)
	./$(TEST_BIN)
	$(SMOKE_RUN) $(CURDIR)/$(BIN_FILE) $(SMOKE_ARGS)

# Test suite and a short solver run under valgrind; any leak, including
# memory still reachable at exit, is an error. Use a serial CPU build: the
# OpenMP and CUDA runtimes keep memory of their own.
valgrind: $(OBJ_DIR) $(BIN_FILE) $(TEST_BIN)
	$(VALGRIND) $(VALGRIND_FLAGS) ./$(TEST_BIN) > /dev/null
	$(SMOKE_RUN) $(VALGRIND) $(VALGRIND_FLAGS) $(CURDIR)/$(BIN_FILE) $(SMOKE_ARGS) > /dev/null

# Invalid command lines must be refused without writing anything
test-cli: $(OBJ_DIR) $(BIN_FILE)
	$(TEST_DIR)/cli.sh $(CURDIR)/$(BIN_FILE)

# A short run against stored results, and the default case against Ghia et al.
regression: $(OBJ_DIR) $(BIN_FILE)
	$(TEST_DIR)/regression.sh $(CURDIR)/$(BIN_FILE)

clean:
	rm -rf $(BIN_FILE) $(TEST_BIN) $(CONV_BIN) $(OBJ_DIR) $(TBN_DIR) output/*.vtk

FORCE:

.PHONY: all test convergence clean FORCE format format-check cppcheck tidy test-asan valgrind test-cli regression
