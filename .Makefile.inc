# Makefile for DynamicRPlugin

# Specify the path to fastjet-config as an input parameter
fastjetprefix = /usr/local
fastjetconfig = ${fastjetprefix}/bin/fastjet-config

# Compiler and flags
#CXX = $(shell $(fastjetconfig) --cxx)
CXX = g++ -std=c++11
CXXFLAGS = $(shell $(fastjetconfig) --cxxflags) -fPIC
LDFLAGS = $(shell $(fastjetconfig) --libs) 
CXXFLAGS += -I.

# Source files
SRC_FILES = DynamicRJetPlugin.cc

# Object files
OBJ_FILES = $(SRC_FILES:.cc=.o)

# Library file
LIBRARY = libDynamicRPlugin.so

# Target for building the shared library
all: $(LIBRARY)

# Rule for building object files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule for building the shared library
$(LIBRARY): $(OBJ_FILES)
	$(CXX) -shared -o $@ $(OBJ_FILES) $(LDFLAGS)

# Target for installing the library
install: all
	# Create directories if not exist
	mkdir -p $(fastjetprefix)/lib
	mkdir -p $(fastjetprefix)/include/fastjet

	# Copy library to lib directory
	cp $(LIBRARY) $(fastjetprefix)/lib

	# Copy header files to include directory
	cp fastjet/DynamicRJetPlugin.hh $(fastjetprefix)/include/fastjet

# Clean up
clean:
	rm -f $(OBJ_FILES) $(LIBRARY)

# Phony target to force make to always rebuild
.PHONY: all install clean

