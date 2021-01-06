SRC := src
OBJ := obj
BUILD := build
CONTAINER := raytracer-memory-test
CC := gcc-10

SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))


FLAGS = -g -Og -Wextra -Werror=implicit-int -Werror=incompatible-pointer-types -Werror=int-conversion -L/usr/local/lib #-fsanitize=address -fsanitize=undefined

LINUX_EXTRA_LIBS= -L/usr/lib/x86_64-linux-gnu -lm
OPENMP=-fopenmp


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


test: build/raytraCer
	./build/raytraCer output.png





docker-build:
	cd memory-test && docker build -t $(CONTAINER) .

docker-run:
	docker run -ti -v ~/Documents/raytraCer:/raytracer $(CONTAINER)
docker-run-memory-test:
	docker run -ti -v ~/Documents/raytraCer:/raytracer $(CONTAINER) bash -c "cd raytracer; make valgrind;"

valgrind: build/raytraCer_linux
	valgrind --leak-check=yes ./build/raytraCer_linux output.ppm
