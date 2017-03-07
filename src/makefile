CC=g++
CFLAGS= -c -g -Wall
LDFLAGS=
INCLUDES =
LIBS = 

SOURCES=main.cpp TEDriver.cpp TEntropy.cpp TimingToSequence.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=DTE

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -std=c++11  $< -o $@

clean:
	rm -rf *o $(EXECUTABLE)
