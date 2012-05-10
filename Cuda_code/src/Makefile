CUDA_INSTALL_PATH ?= /opt/cuda-toolkit

# Directories
BINDIR = ../bin
OBJDIR = ./obj
# Compilers
NVCC       := $(CUDA_INSTALL_PATH)/bin/nvcc
NVINCLUDES  += -I. -I$(CUDA_INSTALL_PATH)/include/
# -I/usr/share/cuda-sdk/C/common/inc
#NVFLAGS =-arch sm_21 -lcublas -L/usr/share/cuda-sdk/C/lib  --compiler-bindir=gcc -m64 --compiler-options '-O3 -fno-strict-aliasing'
NVFLAGS =-arch sm_21 -lcublas --compiler-bindir=gcc -m64 --compiler-options '-O3 -fno-strict-aliasing'
#NVFLAGS_SHARE =-arch sm_21 -shared -Xcompiler -fpic -lcublas -L/usr/share/cuda-sdk/C/lib --compiler-bindir=gcc -m64 --compiler-options '-O3 -fPIC'
NVFLAGS_SHARE =-arch sm_21 -shared -Xcompiler -fpic -lcublas --compiler-bindir=gcc -m64 --compiler-options '-O3 -fPIC'

all: libmatrixMul.so


libmatrixMul.so: matrixMul.cu
	$(NVCC) $(NVINCLUDES) $(NVFLAGS_SHARE) -o $(OBJDIR)/libmatrixMul.so matrixMul.cu

clean:  tidy
	rm $(BINDIR)/* $(OBJDIR)/*

tidy :
	$(VERBOSE)find . | egrep "#" | xargs rm -f
	$(VERBOSE)find . | egrep "\~" | xargs rm -f
