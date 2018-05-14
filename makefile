CC=g++
CFLAGS=-std=c++11 -pthread -O3

myprogram: myprogram.cpp mylib.cpp
	$(CC) -o myprogram myprogram.cpp mylib.cpp $(CFLAGS)
