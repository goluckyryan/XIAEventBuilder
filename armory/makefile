CC=g++

#all: xia2root  xia2ev2_nopart pixie2root scan pxi-time-order
#all: xia2root  xia2ev2_nopart pixie2root scan evt2root evt2hist
#all: xia2root  to2root MergeEVT evt2hist pxi-time-order ev22txt EventBuilder
#all: xia2root  to2root MergeEVT pxi-time-order ev22txt EventBuilder
all: to2root MergeEVT ev22txt EventBuilder pxi-time-order EventBuilder2

#this is FSU evt to root
xia2root: ../armory/xia2root.cpp
	$(CC) ../armory/xia2root.cpp -o xia2root `root-config --cflags --glibs`

#xia2ev2_nopart: armory/xia2ev2_nopart.cpp
#	$(CC) armory/xia2ev2_nopart.cpp -o xia2ev2_nopart

#this is for eventbuild 
to2root: ../armory/to2root.cpp
	$(CC) ../armory/to2root.cpp -o to2root `root-config --cflags --glibs`

#this is for online root
MergeEVT: ../armory/MergeEVT.cpp
	$(CC) ../armory/MergeEVT.cpp -o MergeEVT `root-config --cflags --glibs`

#this is for online spectrums
evt2hist: ../armory/evt2hist.cpp
	$(CC) ../armory/evt2hist.cpp -o evt2hist `root-config --cflags --glibs`

pxi-time-order: ../armory/pxi-time-order.c
	$(CC) ../armory/pxi-time-order.c -o pxi-time-order 
	
ev22txt:  ../armory/ev22txt.cpp
	$(CC) ../armory/ev22txt.cpp -o ev22txt

EventBuilder: ../armory/EventBuilder.cpp
	$(CC) ../armory/EventBuilder.cpp -o EventBuilder `root-config --cflags --glibs`

clean:
	-rm xia2root to2root MergeEVT evt2hist pxi-time-order ev22txt EventBuilder
