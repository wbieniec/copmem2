CPP = g++
FLAGS  = -march=native -O3 -funroll-loops
FLAGS1  = -march=native -O3 -funroll-loops -pthread


MemoryFill.o: MemoryFill.src/MemoryFill.cpp
	$(CPP) -c $(FLAGS) MemoryFill.src/MemoryFill.cpp -o MemoryFill.src/MemoryFill.o


MemoryFill: MemoryFill.o
	$(CPP) $(FLAGS) MemoryFill.src/MemoryFill.o -o memoryfill

format.o: fmt/format.cc fmt/format.h fmt/format.h fmt/format-inl.h fmt/core.h
	$(CPP) $(FLAGS) -c fmt/format.cc -o fmt/format.o

common-xxhash.o: common/xxhash.c common/xxhash.h
	$(CPP) $(FLAGS) -c common/xxhash.c -o common/xxhash.o

common-city.o: common/city.cpp
	$(CPP) $(FLAGS) -c common/city.cpp -o common/city.o

common-metrohash64.o: common/metrohash64.cpp
	$(CPP) $(FLAGS) -c common/metrohash64.cpp -o common/metrohash64.o

common-StopWatch.o: common/StopWatch.cpp common/StopWatch.h
	$(CPP) $(FLAGS) -c common/StopWatch.cpp -o common/StopWatch.o

common-StringUtils.o: common/StringUtils.cpp common/StringUtils.h
	$(CPP) $(FLAGS) -c common/StringUtils.cpp -o common/StringUtils.o

copmem2: format.o common-xxhash.o common-xxhash.o common-city.o common-metrohash64.o common-StopWatch.o common-StringUtils.o CopMEM2.src/CopMEM2.cpp
	$(CPP) $(FLAGS1) CopMEM2.src/CopMEM2.cpp fmt/format.o common/StopWatch.o common/metrohash64.o common/city.o common/StringUtils.o -o copmem2


clean:
	rm -R -f *.o copmem memoryfill
	rm -R -f *.o common/*.o fmt/*.o MemoryFill.src/*.o CopMEM.src/*.o CopMEM2.src/*.o copmem copmem2 memoryfill

all: copmem2 MemoryFill
