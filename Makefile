
EXEC   = GalIC
CONFIG   = Config.sh
BUILD_DIR = build
SRC_DIR = src

PARAMFILE = Model_H1.param
N := 16

###################
#determine SYSTYPE#
###################
ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif

MAKEFILES = Makefile config-makefile
ifeq ($(wildcard Makefile.systype), Makefile.systype)
MAKEFILES += Makefile.systype
endif



PERL	 = /usr/bin/perl
RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) BUILD_DIR=$(BUILD_DIR) make -f config-makefile)
CONFIGVARS := $(shell cat $(BUILD_DIR)/galicconfig.h)



#MPICHLIB = -lmpich
GMPLIB   = -lgmp
GSLLIB   = -lgsl -lgslcblas
MATHLIB  = -lm



# modules for Magny
# module add mvapich2/gcc/64/1.6-qlc

ifeq ($(SYSTYPE),"APHI") 
CC       =   mpicc
CXX      =   mpicxx
#OPTIMIZE =   -g -w -m64 -O3 -msse3
OPTIMIZE =   -g -w -m64 -O3 -march=native
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
else
OPTIMIZE +=  -Wno-unknown-pragmas
endif
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  
FFTW_LIBS=  
GMP_INCL =  
GMP_LIBS =  
MPICHLIB =
HDF5INCL = 
HDF5LIB  = 
#OPT      +=  -DNOCALLSOFSYSTEM
#OPT      +=  -DIMPOSE_PINNING
#OPT      +=  -DUSE_SSE
endif


ifndef LINKER
LINKER = $(CC)
endif



##########################################
#determine the needed object/header files#
##########################################

SUBDIRS = . 

OBJS =   main.o allocate.o  allvars.o  disk.o   grid.o  bulge.o  set_particles.o parallel_sort.o \
	     halo.o  init.o  io.o  mymalloc.o  orbit_response.o  parameters.o  structure.o  system.o  disp_fields.o \
	     forcetree/gravtree.o forcetree/forcetree.o forcetree/forcetree_walk.o domain/peano.o domain/pqueue.o \
	     domain/domain.o domain/domain_balance.o domain/domain_counttogo.o  domain/domain_exchange.o \
	     domain/domain_rearrange.o domain/domain_sort_kernels.o domain/domain_toplevel.o domain/domain_vars.o domain/domain_box.o


INCL += allvars.h proto.h

SUBDIRS += forcetree domain

################################
#determine the needed libraries#
################################


ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
HDF5LIB  = 
endif

ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
THREAD_LIB = 
endif


##########################
#combine compiler options#
##########################

CFLAGS = $(OPTIMIZE) $(OPT) $(HDF5INCL) $(GSL_INCL) $(FFTW_INCL) $(ODE_INCL) $(GMP_INCL) $(MKL_INCL) $(CUDA_INCL) -I$(BUILD_DIR)

LIBS = $(MATHLIB) $(HDF5LIB) $(MPICHLIB) $(GSL_LIBS) $(GSLLIB) $(FFTW_LIB) $(GMP_LIBS) $(GMPLIB) $(ODE_LIB) $(MKL_LIBS) $(THREAD_LIB) $(CUDA_LIBS)


SUBDIRS := $(addprefix $(BUILD_DIR)/,$(SUBDIRS))
OBJS := $(addprefix $(BUILD_DIR)/,$(OBJS)) $(BUILD_DIR)/compile_time_info.o
INCL := $(addprefix $(SRC_DIR)/,$(INCL)) $(BUILD_DIR)/galicconfig.h


################
#create subdirs#
################
RESULT := $(shell mkdir -p $(SUBDIRS)  )



#############
#build rules#
#############

all: $(EXEC)

$(EXEC): $(OBJS)
	$(LINKER) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC)
	mpirun -n $(N) -f hostfile ./$(EXEC) $(PARAMFILE)

bg: $(OBJS)
	$(LINKER) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC)
	mpirun -n $(N) -f hostfile ./$(EXEC) $(PARAMFILE)  1> log.out.txt 2> log.err.txt &

clean:
	rm -f $(OBJS) $(EXEC) lib$(LIBRARY).a
	rm -f $(BUILD_DIR)/compile_time_info.c $(BUILD_DIR)/galicconfig.h

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(INCL) $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/compile_time_info.o: $(BUILD_DIR)/compile_time_info.c $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@
 
