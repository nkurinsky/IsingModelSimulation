CC=g++ -std=c++0x -O1 -g
CFLAGS=-c -Wall -Wextra
SOURCES=main.cpp lattice.cpp modelstats.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=PottsMC_lm

prefix=~/local

all: $(SOURCES) $(EXECUTABLE) install clean

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm *.o $(EXECUTABLE)

install: ${EXECUTABLE}
	install -m 0755 ${EXECUTABLE} $(prefix)/bin

.PHONY: install
