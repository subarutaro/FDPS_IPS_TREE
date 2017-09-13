#======================================================================
#   Numerical Libraries and Compilers
#======================================================================

FDPS_LOC = ../../FDPS
FDPS_INC = -I$(FDPS_LOC)/src
CXXFLAGS_COMMON = -std=c++11 -O3 -ffast-math -funroll-loops $(FDPS_INC)
#CXXFLAGS_COMMON = -std=c++11 -O0 -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)

# Debug macros
#IPS_FLAGS += -DIPS_TREE_FORCE_ERROR_CHECK
#IPS_FLAGS += -DIPS_TREE_PSEUDOPARTICLE_QUADRUPOLE_MOMENT_CHECK
#IPS_FLAGS += -DIPS_TREE_PSEUDOPARTICLE_DISTANCE_CHECK

CXXFLAGS_COMMON += $(IPS_FLAGS)

# [1] Serial
CXX = g++
CXXFLAGS = $(CXXFLAGS_COMMON)
CXXFLAGS += -llapack -lblas -lgfortran
#CXXFLAGS+= -DPARTICLE_SIMULATOR_DEBUG_PRINT
# [2] OpenMP
#CXX = g++
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# [3] MPI
#CXX = mpicxx
#MPI_PATH = $(HOME)/software/openmpi/1.8.1
#CXX = $(MPI_PATH)/bin/mpicxx
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -L$(MPI_PATH)/lib
# [4] MPI + OpenMP
#CXX = mpicxx
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp

#----------------------------------------------------------------------
#   Source and Object files
#----------------------------------------------------------------------
SRCS	= main.cpp
HEDS    = pdb_manager.h unit.h user_defined_class.h constraint.h

TARGET = main.out

$(TARGET): Makefile $(SRCS) $(HEDS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

debug:	Makefile $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -DDEBUG_PRINT -o $(TARGET) $(LDFLAGS)


clean:
	rm -f $(TARGET) *.o

