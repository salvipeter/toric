all: toric-test

LIBGEOM=../libgeom

INCLUDES=-I$(LIBGEOM)
LIBS=-L$(LIBGEOM)/release -lgeom

CXXFLAGS=-Wall -pedantic -std=c++20 -O3 -g $(INCLUDES)

toric-test: toric-test.o toric.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)
