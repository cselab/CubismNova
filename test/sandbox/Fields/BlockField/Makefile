CXXFLAGS ?= -std=c++11 -Wall -Wextra -Wpedantic
CXXFLAGS += -I../../../../include/Cubism

TestBlockField: FieldIndexing.o TestBlockField.cpp
	$(CXX) $(CXXFLAGS) -O0 -o $@ $^

TestBlockFieldAVX256: FieldIndexingAVX256.o TestBlockField.cpp
	$(CXX) $(CXXFLAGS) -O0 -o $@ $^

FieldIndexing.o: FieldIndexing.cpp
	$(CXX) $(CXXFLAGS) -c -O3 -DNDEBUG -o $@ $<

FieldIndexingAVX256.o: FieldIndexing.cpp
	$(CXX) $(CXXFLAGS) -c -O3 -DNDEBUG -DCUBISMTEST_INTRIN -m64 -march=native -mtune=native -o $@ $<

%.s: %.cpp
	$(CXX) $(CXXFLAGS) -c -O3 -DNDEBUG -S -o $@ $<

clean:
	rm -f TestBlockField TestBlockFieldAVX256 *~ *.o *.swp *.s
