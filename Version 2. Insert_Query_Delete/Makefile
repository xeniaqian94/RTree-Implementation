CXX:=g++
CXXFLAGS:=-c
INCLUDES:=
LIBS:=
EXE:=a1

OBJS:=main.o rtree.o rtnode.o boundingbox.o

all: ${EXE}

${EXE}: ${OBJS}
	$(CXX) -o $@ $^ ${LIBS}

%.o: %.cpp
	$(CXX) ${CXXFLAGS} ${INCLUUDES} -o $@ $<

.PHONY: all clean

clean:
	rm -f ${OBJS} ${EXE}
