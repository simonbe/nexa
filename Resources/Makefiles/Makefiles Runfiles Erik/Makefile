CC = mpiCC

Objects = main.o Shapes.o
OutputName = shapes

CFLAGS = -O3
CXX_INCLUDES = -I/usr/include/openmpi -I/usr/local/include -I../network
CXX_LIBS = ../network/libnetwork.a

.SUFFIXES: .cxx .cpp

.cxx.o:
	$(CC) $(CFLAGS) $(CXX_INCLUDES) -c $<

.cpp.o:
	$(CC) $(CFLAGS) $(CXX_INCLUDES) -c $<

#------------ targets --------------------------------------------
# describe how to create the targets - often there will be only one target

.PHONY: all clean

all: network project

project: $(Objects)
	$(CC) $(CFLAGS) $(CXX_INCLUDES) $(Objects) $(CXX_LIBS) -o $(OutputName) 

network:
	cd ../network && $(MAKE)

clean:
	cd ../network && $(MAKE) clean
	rm -f $(Objects) $(OutputName)

