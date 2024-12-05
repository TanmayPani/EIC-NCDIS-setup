
ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = c++
endif

#TODO: Set the path to the root-config, pythia8-config and fastjet-config executables
#Not needed if the executables are in the PATH

ROOT_CONFIG = $(HOME)/heplibs/root/bin/root-config
PYTHIA_CONFIG = $(HOME)/heplibs/pythia8310/bin/pythia8-config
FASTJET_CONFIG = $(HOME)/heplibs/fastjet-install/bin/fastjet-config

CXXFLAGS     =
CXXFLAGS     += $(shell $(ROOT_CONFIG) --cflags)
CXXFLAGS     += $(shell $(PYTHIA_CONFIG) --cxxflags)
CXXFLAGS     += $(shell $(FASTJET_CONFIG) --cxxflags)

LDFLAGS      = 
LDFLAGS     += $(shell $(ROOT_CONFIG) --libs)
LDFLAGS     += $(shell $(FASTJET_CONFIG) --libs --plugins)
LDFLAGS     += $(shell $(PYTHIA_CONFIG) --ldflags)

all: eicJetMaker.exe

eicJetMaker.exe: eicJetMaker.o eicBeamShape.o
	$(CXX) eicJetMaker.o eicBeamShape.o -o eicJetMaker.exe $(LDFLAGS) 

eicJetMaker.o : eicJetMaker.cpp eicBeamShape.hh
	$(CXX) $(CXXFLAGS) -c eicJetMaker.cpp -o eicJetMaker.o -I.

eicBeamShape.o: eicBeamShape.cpp
	$(CXX) $(CXXFLAGS) -c eicBeamShape.cpp -o eicBeamShape.o

clean :
	rm -vf *.o *.exe *~
