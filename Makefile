task1: build build/task1
	./build/task1

build/task%: %.cpp
	clang++ -std=c++23 -O2 $< -o $@

build:
	mkdir build
