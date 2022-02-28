#
# You can set parameters on the command line:
#
# make CC=mpicc
#
# Or run the configure script to set various parameters. Usually
# defaults are all you need, specially if the dependencies are in
# $HOME/gkylsoft and you are using standard compilers (not building on
# GPUs)

ARCH_FLAGS = -march=native
# Warning flags: -Wall -Wno-unused-variable -Wno-unused-function -Wno-missing-braces
CFLAGS = -O3 -g -ffast-math -fPIC -MMD -MP
LDFLAGS =

# CUDA flags
USING_NVCC =
NVCC_FLAGS = 
CUDA_LIBS =
ifeq ($(CC), nvcc)
       USING_NVCC = yes
       CFLAGS = -O3 -g --forward-unknown-to-host-compiler -MMD -MP -fPIC
       NVCC_FLAGS = -x cu -dc -arch=sm_70 --compiler-options="-fPIC" 
       LDFLAGS += -arch=sm_70
       CUDA_LIBS = -lcublas
endif

# Install prefix
PREFIX = ${HOME}/gkylsoft
# Build directory
ifdef USING_NVCC
	BUILD_DIR ?= cuda-build
else
	BUILD_DIR ?= build
endif
# Directories containing source code
SRC_DIRS := minus zero apps kernels

# Default lapack include and libraries: we prefer linking to static library
LAPACK_INC = ${HOME}/gkylsoft/OpenBLAS/include
LAPACK_LIB_DIR = ${HOME}/gkylsoft/OpenBLAS/lib
LAPACK_LIB = ${HOME}/gkylsoft/OpenBLAS/lib/libopenblas.a

# SuperLU includes and librararies
SUPERLU_INC = ${HOME}/gkylsoft/superlu/include
SUPERLU_LIB_DIR = ${HOME}/gkylsoft/superlu/lib
SUPERLU_LIB = ${HOME}/gkylsoft/superlu/lib/libsuperlu.a

# list of includes from kernels
KERN_INC_DIRS = $(shell find $(SRC_DIRS) -type d)
KERN_INCLUDES = $(addprefix -I,$(KERN_INC_DIRS))

# Include config.mak file (if it exists) to overide defaults above
-include config.mak

# determine OS we are running on
UNAME = $(shell uname)

# On OSX we should use Accelerate framework
ifeq ($(UNAME), Darwin)
	LAPACK_INC = zero # dummy
	LAPACK_LIB = -framework Accelerate
	CFLAGS += -DGKYL_USING_FRAMEWORK_ACCELERATE
	SHFLAGS += -dynamiclib -install_name ${PREFIX}/gkylzero/lib/libgkylzero.so
else
	SHFLAGS += -shared
endif

# all includes
INCLUDES = -Iminus -Iminus/STC/include -Izero -Iapps -Iregression ${KERN_INCLUDES} -I${LAPACK_INC} -I${SUPERLU_INC}

# Build commands for CUDA source
$(BUILD_DIR)/%.cu.o: %.cu
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

# Build commands for C source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(ARCH_FLAGS) $(INCLUDES) -c $< -o $@

# Specialized build commands for kernels when using nvcc
ifdef USING_NVCC

# Unfortunately, due to the limitations of the NVCC compiler to treat
# device code in C files, we need to force compile the kernel code
# using the -x cu flag

$(BUILD_DIR)/kernels/bin_op/%.c.o : kernels/bin_op/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/gyrokinetic/%.c.o : kernels/gyrokinetic/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/lbo/%.c.o : kernels/lbo/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/maxwell/%.c.o : kernels/maxwell/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/vlasov/%.c.o : kernels/vlasov/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/basis/%.c.o : kernels/basis/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

endif

# List of source files
SRCS := $(shell find $(SRC_DIRS) -name *.c)
ifdef USING_NVCC
	SRCS += $(shell find zero -name *.cu)
endif

# List of object files that will be built
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

# List of regression and unit tests
REGS := $(patsubst %.c,${BUILD_DIR}/%,$(wildcard regression/rt_*.c))
UNITS := $(patsubst %.c,${BUILD_DIR}/%,$(wildcard unit/ctest_*.c))

# We need to build CUDA unit-test objects
UNIT_CU_SRCS =
UNIT_CU_OBJS =
ifdef USING_NVCC
	UNIT_CU_SRCS = $(shell find unit -name *.cu)
	UNIT_CU_OBJS = $(UNIT_CU_SRCS:%=$(BUILD_DIR)/%.o)
endif

# Header files
HEADERS := $(wildcard minus/*.h) $(wildcard zero/*.h) $(wildcard apps/*.h) $(wildcard kernels/*/*.h)

# Library name
ifeq ($(CC), nvcc)
	G0LIB = gkylzerocuda
else
	G0LIB = gkylzero
endif

# Object files to compile in library
LIBOBJS = ${OBJS}

# names of libraries to build
G0STLIB = lib${G0LIB}.a
G0SHLIB = lib${G0LIB}.so

# Make targets: libraries, regression tests and unit tests
all: ${BUILD_DIR}/${G0STLIB} ${BUILD_DIR}/${G0SHLIB} ${REGS} ${UNITS}

# Library archive
${BUILD_DIR}/${G0STLIB}: ${LIBOBJS}
	ar -crs $@ ${LIBOBJS}

# Dynamic (shared) library
${BUILD_DIR}/${G0SHLIB}: ${LIBOBJS}
	${CC} ${SHFLAGS} ${LDFLAGS} ${LIBOBJS} ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread -o $@

# Regression tests
${BUILD_DIR}/regression/%: regression/%.c ${BUILD_DIR}/${G0STLIB} regression/rt_arg_parse.h
	$(MKDIR_P) ${BUILD_DIR}/regression
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) ${BUILD_DIR}/${G0STLIB} ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread

# Unit tests
${BUILD_DIR}/unit/%: unit/%.c ${BUILD_DIR}/${G0STLIB} ${UNIT_CU_OBJS}
	$(MKDIR_P) ${BUILD_DIR}/unit
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) ${UNIT_CU_OBJS} ${BUILD_DIR}/${G0STLIB} ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread

.PHONY: check check1 clean partclean install

# Run all unit tests
check: ${UNITS}
	$(foreach unit,${UNITS},$(unit);)

install: all
	$(MKDIR_P) ${PREFIX}/gkylzero/include
	${MKDIR_P} ${PREFIX}/gkylzero/lib
	${MKDIR_P} ${PREFIX}/gkylzero/bin
	${MKDIR_P} ${PREFIX}/gkylzero/share
	cp ${HEADERS} ${PREFIX}/gkylzero/include
	cp -f ${BUILD_DIR}/${G0STLIB} ${PREFIX}/gkylzero/lib
	cp -f ${BUILD_DIR}/${G0SHLIB} ${PREFIX}/gkylzero/lib
	cp -f Makefile.sample ${PREFIX}/gkylzero/share/Makefile
	cp -f regression/rt_arg_parse.h ${PREFIX}/gkylzero/share/rt_arg_parse.h
	cp -f regression/rt_twostream.c ${PREFIX}/gkylzero/share/rt_twostream.c
	cp -f ${BUILD_DIR}/regression/rt_vlasov_kerntm ${PREFIX}/gkylzero/bin/

clean:
	rm -rf ${BUILD_DIR}

# partclean does not delete the kernel object files as they do not
# always need to be rebuilt and are very expensive to build when using
# nvcc
partclean:
	rm -rf ${BUILD_DIR}/${G0STLIB} ${BUILD_DIR}/${G0SHLIB} ${BUILD_DIR}/zero ${BUILD_DIR}/app ${BUILD_DIR}/regression ${BUILD_DIR}/unit

# command to make dir
MKDIR_P ?= mkdir -p	
