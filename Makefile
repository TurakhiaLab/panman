CXX := clang++
CXXFLAGS := -std=c++17
LDLIBS := -lprotobuf -ltbb

SRCS := main.cpp PangenomeMAT.cpp mutation_annotation.pb.cc
HDRS := PangenomeMAT.hpp mutation_annotation.pb.h

all: main

main: $(OBJS) $(HDRS)
	$(CXX) $(CXXFLAGS) -o main $(SRCS) $(LDLIBS)

clean:
	rm -rf main