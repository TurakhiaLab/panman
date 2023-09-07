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

libminimap2.a: 
	cd src/minimap2_src \
	make libminimap2.a \

integration_test: libminimap2.a
	gcc -g -O2 src/main_test.c src/minimap2_src/libminimap2.a src/minimap2_src/integration_test.c -o integration_test -lz \
