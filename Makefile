SRC := src
OBJ := obj
BUILD := build
CONTAINER := raytracer-memory-test
CC := gcc-10

SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))
OBJECTS_WITHOUT_MAIN := $(filter-out $(OBJ)/main.o, $(OBJECTS))



LINUX_EXTRA_LIBS= -L/usr/lib/x86_64-linux-gnu -lm
OPENMP = #-fopenmp

OPENCL = -framework OpenCL -arch x86_64 -DUNIX -DDEBUG -DMAC
WARNS = -Wextra -Werror=implicit-int -Werror=incompatible-pointer-types -Werror=int-conversion
LIBS = -L/usr/local/lib $(OPENCL)
FLAGS = -g -Og $(WARNS) $(LIBS)  #-fsanitize=address -fsanitize=undefined


build/raytraCer: clean $(OBJECTS)
	$(CC) $(OBJECTS) $(FLAGS) $(OPENMP) -o $@

build/raytraCer_linux: clean $(OBJECTS)
	$(CC) $(OBJECTS) $(FLAGS) $(LINUX_EXTRA_LIBS) $(OPENMP) -o $@



$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(FLAGS) $(OPENMP) -c $< -o $@

clean:
	rm -f $(BUILD)/raytraCer && rm -f $(OBJ)/*.o


run: build/raytraCer
	./build/raytraCer

test_highres: build/raytraCer
	./build/raytraCer -w 1024 -h 1024 -s 512 -o output_highres.png

test_midres: build/raytraCer
	./build/raytraCer -w 512 -h 512 -s 256 -o output_midres.png

test_lowres: build/raytraCer
	./build/raytraCer -w 256 -h 256 -s 128 -o output_lowres.png












lib: $(OBJECTS_WITHOUT_MAIN)
	rm -f libtracing.a
	ar -cvq libtracing.a $(OBJECTS_WITHOUT_MAIN)
	cp libtracing.a /usr/local/lib/libtracing.a


install: lib
	ln -sf "$(shell pwd)/src" /usr/local/include/raytraCer


docker-build:
	cd memory-test && docker build -t $(CONTAINER) .

docker-run:
	docker run -ti -v ~/Documents/raytraCer:/raytracer $(CONTAINER)
docker-run-memory-test:
	docker run -ti -v ~/Documents/raytraCer:/raytracer $(CONTAINER) bash -c "cd raytracer; make valgrind;"

valgrind: build/raytraCer_linux
	valgrind --leak-check=yes ./build/raytraCer_linux -o output.ppm
