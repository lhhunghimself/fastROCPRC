CXX=g++
CC=gcc
LD=$(CXX)

MYLIBS=-lm -lc -lgcc -lrt


LOPENMP= -lgomp
OPT_FLAGS= -Ofast -ftree-vectorize -funroll-all-loops -m64 -msse3 -mavx
ifdef DEBUG
 OPT_FLAGS= -ggdb -m64 -msse3 -mavx
endif
FFLAGS=-fimplicit-none

CFLAGS= -I ./ -I ./tclap-1.2.1/include $(MAC_INCLUDES) $(OPT_FLAGS) -std=c++11 -Wno-write-strings -fexpensive-optimizations -fopenmp -ffast-math -finline-functions -frerun-loop-opt -static $(HDRS) $(DEFINES) $(OS)
LIBS= $(LIBDIRS) $(MYLIBS) $(LOPENMP)

# object files needed

OBJ_fastROCPRC = main.o

fastROCPRC :$(OBJ_fastROCPRC)
	$(LD) -o fastROCPRC $(LFLAGS) $(OBJ_fastROCPRC) $(LIBS)
main.o: main.cpp
	$(CXX) $(CFLAGS) -c main.cpp
clean:
	rm *.o
	rm core.* 
# all together now...
default: fastROCPRC
