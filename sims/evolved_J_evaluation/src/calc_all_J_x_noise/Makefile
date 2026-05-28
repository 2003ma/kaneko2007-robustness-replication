CC = clang

CPUFLAGS = -march=native

export OMP_NUM_THREADS=20
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores

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

BIN_DIR = bin
TARGET = $(BIN_DIR)/all_J_trial_eval_noise

CORE_DIR = ../core
SRCS = main.c all_J_trial_eval_noise.c \
       $(CORE_DIR)/dynamics.c
OBJS = $(SRCS:.c=.o)

.PHONY: all clean run
all: $(TARGET)

$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

$(TARGET): $(OBJS) | $(BIN_DIR)
	@echo "Linking $(TARGET)"
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.c
	@echo "Compiling $<"
	$(CC) $(CFLAGS) -c $< -o $@

run: $(TARGET)
	OMP_NUM_THREADS=20 $(TARGET) $(ARGS)

clean:
	rm -f $(OBJS) $(TARGET)
