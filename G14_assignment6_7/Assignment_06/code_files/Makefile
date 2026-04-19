CXX = g++
CXXFLAGS = -O3 -fopenmp -std=c++11 -march=native -Wall
TARGET = parallel_interp
SRCS = main.cpp utils.cpp init.cpp

.PHONY: all clean test benchmark

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

clean:
	rm -f $(TARGET) Mesh.out input.bin

# Quick correctness test
test: $(TARGET)
	@echo "=== Testing Correctness ==="
	OMP_NUM_THREADS=4 ./$(TARGET) Test_input.bin
	@if diff -q Mesh.out Test_Mesh.out > /dev/null; then \
		echo "✓ PASS: Output matches reference"; \
	else \
		echo "✗ FAIL: Output differs from reference"; \
		diff Mesh.out Test_Mesh.out | head -20; \
	fi

# Run full benchmark suite
benchmark: $(TARGET)
	@echo "=== Running Benchmark Suite ==="
	@mkdir -p results
	@for cfg in a b c d e; do \
		if [ -f input_$${cfg}.bin ]; then \
			echo "Config $$cfg:"; \
			for t in 1 2 4 8 16; do \
				OMP_NUM_THREADS=$$t ./$(TARGET) input_$${cfg}.bin 2>&1 | grep "Total interpolation"; \
			done; \
		else \
			echo "Warning: input_$${cfg}.bin not found"; \
		fi \
	done