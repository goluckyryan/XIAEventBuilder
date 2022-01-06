CC=g++

#all: xia2root  xia2ev2_nopart pixie2root scan pxi-time-order
#all: xia2root  xia2ev2_nopart pixie2root scan evt2root evt2hist
all: xia2root  to2root evt2root evt2hist pxi-time-order ev22txt EventBuilder

#this is FSU evt to root
xia2root: armory/xia2root.cpp
	$(CC) armory/xia2root.cpp -o xia2root `root-config --cflags --glibs`

#xia2ev2_nopart: armory/xia2ev2_nopart.cpp
#	$(CC) armory/xia2ev2_nopart.cpp -o xia2ev2_nopart

#this is for eventbuild 
to2root: armory/to2root.cpp
	$(CC) armory/to2root.cpp -o to2root `root-config --cflags --glibs`

#this is for online root
evt2root: armory/evt2root.cpp
	$(CC) armory/evt2root.cpp -o evt2root `root-config --cflags --glibs`

#this is for online spectrums
evt2hist: armory/evt2hist.cpp
	$(CC) armory/evt2hist.cpp -o evt2hist `root-config --cflags --glibs`

pxi-time-order: armory/pxi-time-order.c
	$(CC) armory/pxi-time-order.c -o pxi-time-order 
	
ev22txt:  armory/ev22txt.cpp
	$(CC) armory/ev22txt.cpp -o ev22txt

EventBuilder: armory/EventBuilder.cpp
	$(CC) armory/EventBuilder.cpp -o EventBuilder `root-config --cflags --glibs`
