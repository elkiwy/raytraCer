SRC := src
OBJ := obj
BUILD := build
CONTAINER := raytracer-memory-test
CC := gcc-10
#CC := gcc

SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))
OBJECTS_WITHOUT_MAIN := $(filter-out $(OBJ)/main.o, $(OBJECTS))
KERNEL_SRC = program.cl


LINUX_EXTRA_LIBS= -L/usr/lib/x86_64-linux-gnu -lm
OPENMP = -fopenmp
OPENCL = -framework OpenCL -arch x86_64 -DUNIX -DDEBUG -DMAC
WARNS = -Wextra -Werror=implicit-int -Werror=incompatible-pointer-types -Werror=int-conversion
LIBS = -L/usr/local/lib $(OPENCL)
FLAGS = -g -Og $(WARNS) $(LIBS) #-fsanitize=address -fsanitize=undefined




####
#
# Main compilation
#
####

#Compile main and link all the objects
build/raytraCer: clean $(SRC)/kernelSource.h $(OBJECTS)
	$(CC) $(OBJECTS) $(FLAGS) $(OPENMP) -o $@

#Create the header file from the OpenCL source
$(SRC)/kernelSource.h:
	echo '#ifndef _KERNELSRC_H' >> $(SRC)/kernelSource.h
	echo '#define _KERNELSRC_H' >> $(SRC)/kernelSource.h
	xxd -i $(SRC)/$(KERNEL_SRC) >> $(SRC)/kernelSource.h; \
	echo '#endif /* _KERNELSRC_H */' >> $(SRC)/kernelSource.h

#Compile all the objects
$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(FLAGS) $(OPENMP) -c $< -o $@




####
#
# Lib compilation
#
####

#Create libtracing.a archive with all the compiled objects
lib: $(OBJECTS_WITHOUT_MAIN)
	rm -f libtracing.a
	ar -cvq libtracing.a $(OBJECTS_WITHOUT_MAIN)

#SymLink all source files to /usr/local/include/raytraCer and copy the lib file
install: lib
	ln -sf "$(shell pwd)/src" /usr/local/include/raytraCer
	cp libtracing.a /usr/local/lib/libtracing.a



####
#
# Utils
#
####

#Run the program
run: build/raytraCer
	./build/raytraCer

#Clean all the objects file and builds
clean:
	rm -f $(SRC)/kernelSource.h
	rm -f $(BUILD)/raytraCer && rm -f $(OBJ)/*.o
