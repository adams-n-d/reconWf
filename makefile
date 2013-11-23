include /home/nda/code/scicomp/tril11/Makefile.export.Trilinos

CXX=$(Trilinos_CXX_COMPILER)
CFLAGS=-c -g -O2 -Wall
SOURCES=MyLSQR.cpp MatVec.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TPLIBS = -lmatio
EXE=MyLSQR

all: $(SOURCES) $(EXE)
		
$(EXE): $(OBJECTS) 
	$(CXX) $(OBJECTS) $(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) $(LDFLAGS) $(TPLIBS) -o $@

.cpp.o:
	$(CXX) $(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) $(Trilinos_CXX_COMPILER_FLAGS) $(CFLAGS) $(TPLIBS)  $(Trilinos_LIBRARY_DIRS) $< -o $@


clean:
	rm $(EXE) $(OBJECTS)

