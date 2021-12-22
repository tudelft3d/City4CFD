# makefile for LGPL-licensed LAStools
#
#COPTS    = -g -Wall -Wno-deprecated -DDEBUG 
COPTS     = -O3 -Wall -Wno-deprecated -Wno-strict-aliasing -Wno-write-strings -Wno-unused-result -DNDEBUG 
#COMPILER  = CC
COMPILER  = g++
LINKER  = g++
#BITS     = -64

#LIBS     = -L/usr/lib64
#LIBS     = -L/usr/lib32
#INCLUDE  = -I/usr/include

LASLIBS     = -L../LASlib/lib
LASINCLUDE  = -I../LASzip/src -I../LASlib/inc 

all: laszip las2las lasmerge

laszip: laszip.o geoprojectionconverter.o 
	${LINKER} ${BITS} ${COPTS} laszip.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	cp $@ ../bin

las2las: las2las.o geoprojectionconverter.o
	${LINKER} ${BITS} ${COPTS} las2las.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	cp $@ ../bin

lasmerge: lasmerge.o geoprojectionconverter.o 
	${LINKER} ${BITS} ${COPTS} lasmerge.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	cp $@ ../bin

.cpp.o: 
	${COMPILER} ${BITS} -c ${COPTS} ${INCLUDE} $(LASINCLUDE) $< -o $@

.c.o: 
	${COMPILER} ${BITS} -c ${COPTS} ${INCLUDE} $(LASINCLUDE) $< -o $@

clean:
	rm -rf *.o
	rm -rf laszip las2las lasmerge 

clobber:
	rm -rf *.o
	rm -rf laszip las2las lasmerge
	rm -rf *~
