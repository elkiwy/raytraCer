SRC := src
OBJ := obj
BUILD := build

SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))

FLAGS = -g -Og -Wshadow -Wextra -Werror=implicit-int -Werror=incompatible-pointer-types -Werror=int-conversion -fsanitize=address -fsanitize=undefined
#LIBS = -L /usr/local/lib -l SDL2-2.0.0 -l SDL2_ttf -l SDL2_image -l SDL2_gfx


build/raytraCer: clean $(OBJECTS)
	gcc $(OBJECTS) $(FLAGS) $(LIBS) -o $@


$(OBJ)/%.o: $(SRC)/%.c
	gcc $(FLAGS) -c $< -o $@

clean:
	rm -f $(BUILD)/raytraCer && rm -f $(OBJ)/*.o

run: build/raytraCer
	./build/raytraCer


test: build/raytraCer
	./build/raytraCer output.ppm
