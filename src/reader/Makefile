CXX=g++
CXXFLAGS=-g -O3 -Wall -fopenmp  -std=c++17 -D_GLIBCXX_PARALLEL -Wextra -Wshadow -Wfloat-equal -Wundef -Wstrict-overflow=5 -Wwrite-strings -Wcast-qual -Wswitch-default -march=native -Wswitch-enum -Wconversion -Wunreachable-code -Wformat=2
RM=rm -f
SRCS=test.cpp reader.cpp	event.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: reader

reader: $(OBJS)
	$(CXX)	$(LDFLAGS)	-o	read	$(OBJS)	$(LDLIBS)

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~ .depend

include .depend
