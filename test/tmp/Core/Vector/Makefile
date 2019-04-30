CXXFLAGS ?= -std=c++11 -Wall -Wextra -Wpedantic
CXXFLAGS += -I../../../include/Cubism

TestVector: TestVector.cpp
	$(CXX) $(CXXFLAGS) -O0 -o $@ $^

%.s: %.cpp
	$(CXX) $(CXXFLAGS) -c -O3 -DNDEBUG -S -o $@ $<

clean:
	rm -f TestVector *~ *.o *.swp *.s
