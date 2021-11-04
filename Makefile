#
# 'make'        build executable file 'main'
# 'make clean'  removes all .o and executable files
#

# if not already defined, set the HDF5 Home
HDF5HOME ?= /home/alex/hdf5

# define the Cpp compiler to use
CXX  ?= g++
NVCC ?= nvcc

# define any compile-time flags
CXXFLAGS	:= -std=c++11 -Wall -Wextra -O3
NVFLAGS     := -O3 --gpu-architecture=sm_35 -Wno-deprecated-gpu-targets

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = -lhdf5

# define output directory
OUTPUT	:= bin

# define build directory
BUILD   := build

# define source directory
SRC		:= src

# define include directory
INCLUDE	:= include $(HDF5HOME)/include

# define lib directory
LIB		:= $(HDF5HOME)/lib

ifeq ($(OS),Windows_NT)
MAIN_CPU	:= main.exe
MAIN_GPU	:= main_gpu.exe
SOURCEDIRS	:= $(SRC)
INCLUDEDIRS	:= $(INCLUDE)
LIBDIRS		:= $(LIB)
FIXPATH = $(subst /,\,$1)
RM			:= del /q /f
MD	:= mkdir
else
MAIN_CPU	:= main
MAIN_GPU	:= main_gpu
SOURCEDIRS	:= $(shell find $(SRC) -type d)
INCLUDEDIRS	:= $(shell find $(INCLUDE) -type d)
LIBDIRS		:= $(shell find $(LIB) -type d)
FIXPATH = $1
RM = rm -f
MD	:= mkdir -p
endif

# define any directories containing header files other than /usr/include
INCLUDES	:= $(patsubst %,-I%, $(INCLUDEDIRS:%/=%))

# define the C libs
LIBS		:= $(patsubst %,-L%, $(LIBDIRS:%/=%))

# define the C source files
SOURCES		:= $(wildcard $(patsubst %,%/*.c*, $(SOURCEDIRS)))

# define the C object files
OBJECTS		:= $(SOURCES:.cpp=.o)
OBJECTS		:= $(OBJECTS:.cu=.o)

OBJECTS_CPU := $(filter-out src/main_gpu.o, $(OBJECTS))
OBJECTS_GPU := $(filter-out src/main.o, $(OBJECTS))

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

OUTPUTMAIN_CPU	:= $(call FIXPATH,$(OUTPUT)/$(MAIN_CPU))
OUTPUTMAIN_GPU	:= $(call FIXPATH,$(OUTPUT)/$(MAIN_GPU))

all: $(OUTPUT) $(MAIN_CPU) $(MAIN_GPU)
	@echo Executing 'all' complete!

cpu: $(OUTPUT) $(MAIN_CPU)
	@echo Executing 'all' complete!

gpu: $(OUTPUT) $(MAIN_GPU)
	@echo Executing 'gpu' complete!

$(OUTPUT):
	$(MD) $(OUTPUT)

$(MAIN_CPU): $(OBJECTS_CPU)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(OUTPUTMAIN_CPU) $(OBJECTS_CPU) $(LFLAGS) $(LIBS)

$(MAIN_GPU): $(OBJECTS_GPU)
	$(NVCC) $(NVFLAGS) $(INCLUDES) -o $(OUTPUTMAIN_GPU) $(OBJECTS_GPU) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

src/main_gpu.o:
	$(NVCC) $(NVFLAGS) $(INCLUDES) -c src/main_gpu.cu -o $@

.PHONY: clean
clean:
	$(RM) $(call FIXPATH,$(OBJECTS))
	@echo Cleanup complete!