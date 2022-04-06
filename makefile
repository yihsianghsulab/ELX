EXECUTABLES = eLX
TARGETS = $(EXECUTABLES)

NM := $(patsubst %.cpp,%.o,$(wildcard /newmat11/*.cpp))

AR = ar
CXX = g++
CXXFLAGS = -O3 -Wno-deprecated -I./ 
CPPSRC	= eLCX_main.cpp eLC_model.cpp cdflib.cpp 
OBJS	= eLCX_main.o eLC_model.o cdflib.o
LIBPATH = -L./
LDFLAGS = $(LIBPATH) -lm -lnewmat

.PHONY: test

all: $(EXECUTABLES)

eLX: eLX.o
	$(CXX) -o $@ $(OBJS) libnewmat.a

eLX.o:
	$(CXX) -c $(CPPSRC) $(CXXFLAGS) 

libnewmat.a : $(NM)
	$(AR) rc $@ $^

clean : 
	rm -f *.o ./newmat11/*.o *.a
	rm -f $(EXECUTABLES)
