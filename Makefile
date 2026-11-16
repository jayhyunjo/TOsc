CC=g++
CC+=-DDEBUG -g  
CFLAGS=-c -Wall -m64
LDFLAGS=-fPIC
DIR_SRC = ./src
SOURCES=read_oscillation_v01.cxx $(wildcard $(DIR_SRC)/*.cxx)
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=read_oscillation_v01

#ROOTSYS=/home/xji/data0/software/root_build

CFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS += $(shell $(ROOTSYS)/bin/root-config --libs) 

#CFLAGS += -std=c++11 -I./inc/ -I$(ROOTSYS)/include/
CFLAGS += -I./inc/ -I$(ROOTSYS)/include/

SCAN=scan_appearance_cls
SCANTOYS=scan_appearance_cls_toys
CHECKMAP=check_dchi2_map

all: $(SOURCES) $(EXECUTABLE) $(SCAN) $(SCANTOYS) $(CHECKMAP)

$(EXECUTABLE):$(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) -lMinuit2

$(SCAN): $(SCAN).o $(DIR_SRC)/TOsc.o
	$(CC) -o $@ $(SCAN).o $(DIR_SRC)/TOsc.o $(LDFLAGS) -lMinuit2

$(SCANTOYS): $(SCANTOYS).o $(DIR_SRC)/TOsc.o
	$(CC) -o $@ $(SCANTOYS).o $(DIR_SRC)/TOsc.o $(LDFLAGS) -lMinuit2

#.C.o:
%.o:%.cxx
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(DIR_SRC)/*.o; rm $(EXECUTABLE); rm -f *.pcm *.d *.so
canv:
	rm -f canv*
rroot:
	rm -f *.root

$(CHECKMAP): $(CHECKMAP).o $(DIR_SRC)/TOsc.o
	$(CC) -o $@ $(CHECKMAP).o $(DIR_SRC)/TOsc.o $(LDFLAGS) -lMinuit2
