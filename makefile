CC=g++

all: xia2root  xia2ev2_nopart pixie2root scan pxi-time-order

xia2root: xia2root.cpp
	$(CC) xia2root.cpp -o xia2root `root-config --cflags --glibs`

xia2ev2_nopart: xia2ev2_nopart.cpp
	$(CC) xia2ev2_nopart.cpp -o xia2ev2_nopart

pixie2root: pixie2root.cpp
	$(CC) pixie2root.cpp -o pixie2root `root-config --cflags --glibs`

scan: scan.c
	$(CC) scan.c -o scan

pxi-time-order: pxi-time-order.c
	$(CC) pxi-time-order.c -o pxi-time-order 
