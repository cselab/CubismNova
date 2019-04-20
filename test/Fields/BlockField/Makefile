CXXFLAGS ?= -std=c++11 -Wall -Wextra -Wpedantic
CXXFLAGS += -I../../../include/Cubism

TestBlockField: FieldIndexing.o TestBlockField.cpp
	$(CXX) $(CXXFLAGS) -O0 -o $@ $^

FieldIndexing.o: FieldIndexing.cpp
	$(CXX) $(CXXFLAGS) -c -O3 -DNDEBUG -o $@ $<

%.s: %.cpp
	$(CXX) $(CXXFLAGS) -c -O3 -DNDEBUG -S -o $@ $<

clean:
	rm -f TestBlockField *~ *.o *.swp *.s
