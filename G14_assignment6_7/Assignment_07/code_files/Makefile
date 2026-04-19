CXX      = g++
CXXFLAGS = -O3 -march=native -funroll-loops
OMPFLAGS = -fopenmp
LIBS     = -lm

# ---------- Serial ----------
SERIAL_SRCS = main.cpp utils.cpp init.cpp
SERIAL_OBJS = $(SERIAL_SRCS:.cpp=.o)
SERIAL_OUT  = main.out

# ---------- Parallel ----------
OMP_SRCS = main_omp.cpp utils_omp.cpp init.cpp
OMP_OBJS = $(OMP_SRCS:.cpp=_omp.o)
OMP_OUT  = main_omp.out

# ---------- Input maker ----------
MAKER_OUT = input_maker.out

.PHONY: all serial parallel maker clean

all: serial parallel maker

serial: $(SERIAL_OUT)

parallel: $(OMP_OUT)

maker: $(MAKER_OUT)

# Serial objects
main.o: main.cpp init.h utils.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

utils.o: utils.cpp utils.h init.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

init.o: init.cpp init.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SERIAL_OUT): main.o utils.o init.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

# Parallel objects (compiled with OpenMP)
main_omp_omp.o: main_omp.cpp init.h utils.h
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -c $< -o $@

utils_omp_omp.o: utils_omp.cpp utils.h init.h
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -c $< -o $@

init_omp.o: init.cpp init.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OMP_OUT): main_omp_omp.o utils_omp_omp.o utils.o init.o
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ $(LIBS) -o $@

$(MAKER_OUT): input_file_maker.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o *.out Mesh.out input.bin
