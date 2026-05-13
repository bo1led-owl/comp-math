SRC = $(wildcard ./*.c)
CC_FLAGS = -std=c23 -lm -Wall -Wextra -Wpedantic -O3 -flto=full -march=native
# CC_FLAGS = -std=c23 -lm -Wall -Wextra -Wpedantic -fsanitize=address,undefined -fsanitize=leak -g

run: build/main
	./build/main

clean:
	rm -r build

build/main: build $(SRC)
	$(CC) $(CC_FLAGS) $(SRC) -o build/main

build:
	mkdir build/

.PHONY: run clean
