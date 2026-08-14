# ----- Compiler -----
CC = clang

# ----- CPU optimization (strong for both Intel and Apple Silicon) -----
CPUFLAGS = -march=native

# ----- Force 20 OpenMP threads -----
export OMP_NUM_THREADS=20
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores

# ----- Compile Flags -----
CFLAGS = -O3 $(CPUFLAGS) \
         -fno-math-errno -fno-trapping-math \
         -DNDEBUG -std=c11 \
         -Wall -Wextra -Wno-unused-parameter \
         -Wshadow -Wpointer-arith -Wstrict-prototypes \
         -funroll-loops \
         -Xpreprocessor -fopenmp \
         -I/opt/homebrew/opt/libomp/include \
         -I../core

LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp -lm

# ----- Directories / Files -----
BIN_DIR = bin
TARGET  = $(BIN_DIR)/all_J_trial_eval

CORE_DIR = ../core
SRCS = main.c all_J_trial_eval.c \
       $(CORE_DIR)/dynamics.c
OBJS = $(SRCS:.c=.o)

# ----- Default target -----
.PHONY: all clean run
all: $(TARGET)

# ----- Create bin directory -----
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

# ----- Linking -----
$(TARGET): $(OBJS) | $(BIN_DIR)
	@echo "🔗 Linking $(TARGET)"
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# ----- Generic compile rule -----
%.o: %.c
	@echo "⚙️  Compiling $<"
	$(CC) $(CFLAGS) -c $< -o $@

# ----- Run -----
run: $(TARGET)
	@echo "🚀 Running with OMP_NUM_THREADS=20"
	OMP_NUM_THREADS=20 $(TARGET) $(ARGS)

# ----- Clean -----
clean:
	@echo "🧹 Cleaning..."
	rm -f $(OBJS) $(TARGET)
