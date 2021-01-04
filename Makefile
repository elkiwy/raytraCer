SRC := src
OBJ := obj
BUILD := build
CONTAINER := raytracer-memory-test

SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))


FLAGS = -g -Og -Wshadow -Wextra -Werror=implicit-int -Werror=incompatible-pointer-types -Werror=int-conversion
FLAGS_SANITIZE = -fsanitize=address -fsanitize=undefined
#LIBS = -L /usr/local/lib -l SDL2-2.0.0 -l SDL2_ttf -l SDL2_image -l SDL2_gfx

LINUX_EXTRA_LIBS=-L/usr/lib/x86_64-linux-gnu -lm


build/raytraCer: clean $(OBJECTS)
	gcc $(OBJECTS) $(FLAGS) $(FLAGS_SANITIZE) $(LIBS) $(LINUX_EXTRA_LIBS) -o $@


$(OBJ)/%.o: $(SRC)/%.c
	gcc $(FLAGS) -c $< -o $@

clean:
	rm -f $(BUILD)/raytraCer && rm -f $(OBJ)/*.o

run: build/raytraCer
	./build/raytraCer


test: build/raytraCer
	./build/raytraCer output.ppm



docker-build:
	cd memory-test && docker build -t $(CONTAINER) .

docker-run:
	docker run -ti -v ~/Documents/raytraCer:/raytracer $(CONTAINER)
docker-run-memory-test:
	docker run -ti -v ~/Documents/raytraCer:/raytracer $(CONTAINER) bash -c "cd raytracer; make valgrind;"

valgrind: build/raytraCer
	valgrind --leak-check=yes ./build/raytraCer output.ppm
