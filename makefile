CC=g++

all: xia2root
#all: xia2root  xia2ev2

xia2root: xia2root.cpp
	$(CC) xia2root.cpp -o xia2root `root-config --cflags --glibs`

xia2ev2: xia2ev2_nopart.cpp
	$(CC) xia2ev2_nopart.cpp -o xia2ev2_nopart
