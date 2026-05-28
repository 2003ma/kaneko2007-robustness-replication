# --- Compiler selection ---
CC = clang    
# --- CPU 最適化設定（Apple / Intel 共通で最強） ---
CPUFLAGS = -mcpu=native

# --- OpenMP 強制フル稼働 ---
export OMP_NUM_THREADS=20
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores

# --- Base C flags ---
CFLAGS = -O3 $(CPUFLAGS) \
         -fno-math-errno -fno-trapping-math \
         -DNDEBUG -std=c17 \
         -Wall -Wextra -Wno-unused-parameter \
         -Xpreprocessor -fopenmp \
         -I/opt/homebrew/opt/libomp/include

LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp -lm

# --- Sources / Output ---
SRCDIR = src
MAINDIR = main
SRCS = $(MAINDIR)/main.c $(SRCDIR)/kaneko_robustness.c $(SRCDIR)/save_all_J.c
OBJS = $(SRCS:.c=.o)

BINDIR = bin
BIN = $(BINDIR)/collect_J

# --- Default target ---
.PHONY: all
all: $(BIN)

# --- Link ---
$(BIN): $(OBJS) | $(BINDIR)
	@echo "🔗 Linking $@"
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

# --- Compile rule ---
%.o: %.c
	@echo "⚙️  Compiling $<"
	$(CC) $(CFLAGS) -c $< -o $@

# --- Create bin directory ---
$(BINDIR):
	@mkdir -p $(BINDIR) results

# --- Run target ---
run: $(BIN)
	@echo "🚀 Running with ARGS='$(ARGS)'"
	OMP_NUM_THREADS=20 $(BIN) $(ARGS)

# --- Cleanup ---
.PHONY: clean
clean:
	@echo "🧹 Cleaning..."
	rm -f $(OBJS) $(BIN)
