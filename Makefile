# Makefile for locate-cdawg test code
# Note that for this to compile correctly, you must have installed the sdsl library on your system (go to extern/sdsl-lite and ./install.sh)

INCLUDE_FLAGS=-I ~/include -I extern/BWTIL/extern/bitvector/include -I extern/BWTIL/extern -I extern/BWTIL/data_structures -I extern/rlcsa -I include

#LIB_FLAGS: lib directories
LIB_DIR_FLAGS=-L ~/lib

#library flags
LIB_FLAGS=-lsdsl -ldivsufsort -ldivsufsort64

#rlcsa flags:
# Use 64-bit integers in a 64-bit environment.
SIZE_FLAGS = -DMASSIVE_DATA_RLCSA

# Parallelism is supported by either libstdc++ Parallel Mode or MCSTL.
PARALLEL_FLAGS = -DMULTITHREAD_SUPPORT -D_GLIBCXX_PARALLEL -fopenmp
# MCSTL_ROOT = /fs-3/d/jltsiren/suds/mcstl
# PARALLEL_FLAGS = -DMULTITHREAD_SUPPORT -I$(MCSTL_ROOT)/c++ -fopenmp

# Vectors using nibble codes instead of delta codes are faster, but they also
# take up more space.
VECTOR_FLAGS = $(PSI_FLAGS) $(LCP_FLAGS) $(SA_FLAGS)
# PSI_FLAGS = -DUSE_NIBBLE_VECTORS
# LCP_FLAGS = -DSUCCINCT_LCP_VECTOR
# SA_FLAGS = -DSUCCINCT_SA_VECTOR

RLCSA_FLAGS = $(SIZE_FLAGS) $(PARALLEL_FLAGS) $(VECTOR_FLAGS)

CLANG_CXXFLAGS=-Weverything -pedantic -Wno-c++98-compat -Wno-c++98-compat-pedantic $(RLCSA_FLAGS)
GCC_CXXFLAGS=-Wall $(RLCSA_FLAGS)

RLCSA_LIB = extern/rlcsa/librlcsa.a

ifndef STD
	STD=c++11
endif

ifndef OPTFLAGS
	ifneq (yes,$(DEBUG))		
		ARCHFLAGS=-march=native -mtune=generic
		ifneq (yes, $(ASSERTS))
			ASSERTSFLAGS=-DNDEBUG
		else
			ASSERTSFLAGS=
		endif
		OPTFLAGS=-Ofast -fstrict-aliasing -g -ggdb $(ASSERTSFLAGS) $(ARCHFLAGS)
	else
		OPTFLAGS=-O0 -ggdb -g
	endif
endif

CCVER=$(shell $(CXX) --version)

ifneq (,$(findstring g++,$(CXX)))
	MY_CXXFLAGS=$(GCC_CXXFLAGS) $(CXXFLAGS)
else
	MY_CXXFLAGS=$(CLANG_CXXFLAGS) $(CXXFLAGS)
endif


CCVER=$(shell $(CXX) --version)

ifneq (,$(findstring g++,$(CXX)))
	MY_CXXFLAGS=$(GCC_CXXFLAGS) $(CXXFLAGS)
else
	MY_CXXFLAGS=$(CLANG_CXXFLAGS) $(CXXFLAGS)
endif

all: 
	@echo "Compiling with $(CXX)..."
	make clean
	make integer_coding
	make rlcsa-lib
	make -C extern/rlcsa/ build_rlcsa
	make cdawg_utils
	make cdawg_build
	make classic_search_dawg
	make succ_search_dawg

#call this to build all rlcsa dependencies on which lz-rlcsa depends

integer_coding: integer_coding.h
	@echo "Compiling with $(CXX)..."
	$(CXX) integer_coding.c -c 


rlcsa-lib:
	make librlcsa.a -C extern/rlcsa

cdawg_utils: $(RLCSA_LIB) integer_coding cdawg_utils.h dynamic.h
	@echo "Compiling with $(CXX)..."
	$(CXX) 	-std=$(STD) $(MY_CXXFLAGS) $(OPTFLAGS) \
		$(INCLUDE_FLAGS) $(LIB_DIR_FLAGS)  \
		 -c cdawg_utils.cpp  \
		$(LIB_FLAGS) 



classic_search_dawg: $(RLCSA_LIB) integer_coding cdawg_utils dynamic.h
	@echo "Compiling with $(CXX)..."
	$(CXX) 	-std=$(STD) $(MY_CXXFLAGS) $(OPTFLAGS) \
		$(INCLUDE_FLAGS) $(LIB_DIR_FLAGS)  \
		 -o classic_search_dawg.x classic_search_dawg.cpp $(RLCSA_LIB) \
		$(LIB_FLAGS) integer_coding.o cdawg_utils.o

succ_search_dawg: $(RLCSA_LIB) integer_coding cdawg_utils dynamic.h
	@echo "Compiling with $(CXX)..."
	$(CXX) 	-std=$(STD) $(MY_CXXFLAGS) $(OPTFLAGS) \
		$(INCLUDE_FLAGS) $(LIB_DIR_FLAGS)  \
		 -o succ_search_dawg.x succ_search_dawg.cpp $(RLCSA_LIB) \
		$(LIB_FLAGS) integer_coding.o cdawg_utils.o

cdawg_build: $(RLCSA_LIB) integer_coding cdawg_utils dynamic.h
	@echo "Compiling with $(CXX)..."
	$(CXX) 	-std=$(STD) $(MY_CXXFLAGS) $(OPTFLAGS) \
		$(INCLUDE_FLAGS) $(LIB_DIR_FLAGS)  \
		 -o cdawg_build.x cdawg_build.cpp $(RLCSA_LIB) \
		$(LIB_FLAGS) integer_coding.o cdawg_utils.o



clean:
	rm -f cdawg_build.x integer_coding.o cdawg_utils.o succ_search_dawg.x  classic_search_dawg.x cdawg_build.x rlcsa-search.x

clean-rlcsa:
	rm extern/rlcsa/*.o extern/rlcsa/*.a extern/rlcsa/bits/*.o extern/rlcsa/misc/*.o 

clean-all:
	make clean
	make clean-rlcsa
