CC=g++ -std=c++0x
CFLAGS=-c -Wall
SOURCES=main.cpp lattice.cpp modelstats.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=PottsMC

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm *.o $(EXECUTABLE)
