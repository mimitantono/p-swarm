COMMON=-g # -pg
COMPILEOPT=-Wall -msse -msse2 -mssse3 -msse4.2 -mtune=core2 -Icityhash 
CXX = g++ 
FLAGS = -c $(INCLUDES)
CXXFLAGS=$(COMPILEOPT) $(COMMON) -O3
OBJS=main.o matrix.o cpu_info.o qgram.o cityhash/city.o \
	util.o db.o scan.o search8.o search16.o property.o \
	clusterresult.o cluster.o seqinfo.o clusterdata.o \
#	shd/vector_filter.o shd/mask.o shd/bit_convert.o shd/popcount.o 
OBJS_TEST=main.s matrix.s cpu_info.s qgram.s cityhash/city.s \
	util.s db.s scan.s search8.s search16.s property.s \
	clusterresult.s cluster.s seqinfo.s clusterdata.s \
#	shd/vector_filter.s shd/mask.s shd/bit_convert.s shd/popcount.s 
LINKFLAGS=-g 
LIBS=-lpthread

DEPS=Makefile main.h cityhash/config.h cityhash/city.h \
#	shd/vector_filter.h shd/mask.h shd/bit_convert.h shd/popcount.h

%.o : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.s : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -D DEBUG -c -S -o $@ $<

all : main test

main :  $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)

test :  $(OBJS_TEST)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS_TEST) $(LIBS)

clean :
	rm -f *.o *.s *~ $(PROG) gmon.out 
	rm -f main test
