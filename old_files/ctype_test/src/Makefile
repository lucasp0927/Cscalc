# Directories
BINDIR = ../bin
OBJDIR = ./obj
# Compilers
GCC       := g++
#INCLUDESFLAG  += -I. -I$(CUDA_INSTALL_PATH)/include/
# -I/usr/share/cuda-sdk/C/common/inc
#NVFLAGS =-arch sm_21 -lcublas -L/usr/share/cuda-sdk/C/lib  --compiler-bindir=gcc -m64 --compiler-options '-O3 -fno-strict-aliasing'
#NVFLAGS =-O3 -fno-strict-aliasing
#NVFLAGS_SHARE =-arch sm_21 -shared -Xcompiler -fpic -lcublas -L/usr/share/cuda-sdk/C/lib --compiler-bindir=gcc -m64 --compiler-options '-O3 -fPIC'
FLAGS_SHARE =-shared -O3 -fPIC -lm -fopenmp

all: libmatrixMul.so


libmatrixMul.so: matrixMul.cpp
	$(GCC) $(FLAGS_SHARE) -o $(OBJDIR)/libmatrixMul.so matrixMul.cpp

clean:  tidy
	rm $(BINDIR)/* $(OBJDIR)/*

tidy :
	$(VERBOSE)find . | egrep "#" | xargs rm -f
	$(VERBOSE)find . | egrep "\~" | xargs rm -f
