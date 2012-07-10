# Directories
BINDIR = ../bin
OBJDIR = ./obj
# Compilers
GCC       := g++
ICC       := icc
#INCLUDESFLAG  += -I. -I$(CUDA_INSTALL_PATH)/include/
# -I/usr/share/cuda-sdk/C/common/inc
#NVFLAGS =-arch sm_21 -lcublas -L/usr/share/cuda-sdk/C/lib  --compiler-bindir=gcc -m64 --compiler-options '-O3 -fno-strict-aliasing'
#NVFLAGS =-O3 -fno-strict-aliasing
#NVFLAGS_SHARE =-arch sm_21 -shared -Xcompiler -fpic -lcublas -L/usr/share/cuda-sdk/C/lib --compiler-bindir=gcc -m64 --compiler-options '-O3 -fPIC'
ICCFLAGS_SHARE =-shared -O3 -fPIC -openmp -xHost -fast -unroll-aggressive -wd981 -mkl -Wall #operand order warning
#FLAGS_SHARE =-shared -O3 -fPIC -fopenmp -Wall
FLAGS_SHARE =-shared -O3 -fPIC  -Wall

all: libcumtrapz.so


libcumtrapz.so: cumtrapz.cpp
	$(ICC) $(ICCFLAGS_SHARE) -o $(OBJDIR)/libcumtrapz.so cumtrapz.cpp

clean:  tidy
	rm $(BINDIR)/* $(OBJDIR)/*

tidy :
	$(VERBOSE)find . | egrep "#" | xargs rm -f
	$(VERBOSE)find . | egrep "\~" | xargs rm -f
