NUM_NODES=30000
DENSITY=0.05
SEED=0
NUM_THREADS=2
NUM_PROC = 2

CC=g++
CXXFLAGS=-Wall -Werror -g -DNUM_NODES=$(NUM_NODES) -DDENSITY=$(DENSITY) -DSEED=$(SEED) -DNUM_THREADS=$(NUM_THREADS)
LDFLAGS=-lm

TARGETS=serial
SOURCES=$(TARGETS:=.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLES_WIN=$(TARGETS:=.exe) parallel_omp.exe parallel_pthread.exe parallel_mpi.exe parallel_mpi_omp.exe

all: $(TARGETS) parallel_mpi_omp
build: $(TARGETS) parallel_mpi_omp

debug: CXXFLAGS += -DDEBUG
debug: all

$(OBJECTS): $(SOURCES)

$(TARGETS): $(OBJECTS)

parallel_omp: CXXFLAGS += -fopenmp
parallel_omp: LDFLAGS += -fopenmp

parallel_mpi: parallel_mpi.cpp
	mpic++ $(CXXFLAGS) -o $@ parallel_mpi.cpp $(LDFLAGS)

parallel_pthread: LDFLAGS += -lpthread
parallel_pthread:
	$(CC) $(CXXFLAGS) -o $@ parallel_pthread.cpp $(LDFLAGS)

parallel_mpi_omp: parallel_mpi_omp.cpp
	mpic++ $(CXXFLAGS) -fopenmp -o $@ parallel_mpi_omp.cpp $(LDFLAGS)

run-all: run-serial run-omp run-pthread run-mpi run-mpi-omp

run-serial: serial
	./serial

run-omp: parallel_omp
	./parallel_omp

run-pthread: parallel_pthread
	./parallel_pthread

run-mpi: parallel_mpi
	mpirun -np $(NUM_THREADS) --oversubscribe ./parallel_mpi

run-mpi-omp: parallel_mpi_omp
	mpirun -np $(NUM_PROC) --oversubscribe ./parallel_mpi_omp

OUT_FILES=out out_omp out_pthread out_mpi out_mpi_omp

clean:
	rm -f $(TARGETS) parallel_omp parallel_pthread parallel_mpi parallel_mpi_omp *.o $(OUT_FILES)

clean-win:
	del $(EXECUTABLES_WIN) *.o $(OUT_FILES)