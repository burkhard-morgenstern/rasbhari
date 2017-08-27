CC=g++
CFLAGS=-Wall -O3 -std=c++11 # -g -fno-omit-frame-pointer -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -Wl,--no-as-needed -lprofiler -ltcmalloc -Wl,--as-needed 
SOURCES=src/main.cpp src/rasbimp.cpp src/rasbopt.cpp src/rasbhari.cpp src/rasbcomp.cpp src/speedsens.cpp src/sensmem.cpp src/patternset.cpp src/pattern.cpp
HEADER=src/rasbimp.hpp src/rasbopt.hpp src/rasbhari.hpp src/rasbcomp.hpp src/speedsens.hpp src/sensmem.hpp src/patternset.hpp src/pattern.hpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=rasbhari

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: $(HEADER)
		$(CC) -c $(CFLAGS) $< -o $@
clean:
	find ./src/ -name "*.o" -delete
	find ./ -name $(EXECUTABLE) -delete
