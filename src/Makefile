CC	= g++
ARGS	= -O3 -pthread -mavx2 -msse4

all : TEST

TEST:  main.o seq.o spec_avx.o spec_nonSIMD.o spec_avxunroll.o sse.o OffLine.o
	$(CC) $(ARGS) $^ -o $@
	rm *.o

main.o: main.cpp OffLine.hpp seq.hpp spec_nonSIMD.hpp spec_avx.hpp spec_avxunroll.hpp sse.hpp
	$(CC) $(ARGS) -c $<

seq.o: seq.cpp seq.hpp OffLine.hpp
	$(CC) $(ARGS) -c $<

spec_nonSIMD.o: spec_nonSIMD.cpp spec_nonSIMD.hpp OffLine.hpp
	$(CC) $(ARGS) -c $<

spec_avx.o: spec_avx.cpp spec_avx.hpp OffLine.hpp
	$(CC) $(ARGS) -c $< 

spec_avxunroll.o: spec_avxunroll.cpp spec_avxunroll.hpp OffLine.hpp
	$(CC) $(ARGS) -c $< 

sse.o: sse.cpp sse.hpp OffLine.hpp
	$(CC) $(ARGS) -c $< 

OffLine.o: OffLine.cpp OffLine.hpp
	$(CC) $(ARGS) -c $< 

clean :
	rm -f *.o 

