COMMON=-g # -pg
COMPILEOPT=-Wall -msse -msse2 -mssse3 -mtune=core2 -Icityhash 
CXX = g++ 
FLAGS = -c $(INCLUDES)
CXXFLAGS=$(COMPILEOPT) $(COMMON) -O3
OBJS=main.o matrix.o cpu_info.o qgram.o cityhash/city.o \
     util.o db.o scan.o search8.o search16.o property.o \
     clusterresult.o cluster.o 
OBJS_TEST=main.s matrix.s cpu_info.s qgram.s cityhash/city.s \
     util.s db.s scan.s search8.s search16.s property.s \
     clusterresult.s cluster.s 
LINKFLAGS=-g 
LIBS=-lpthread

DEPS=Makefile main.h cluster.h cityhash/config.h cityhash/city.h 

.SUFFIXES:.o .cc

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
