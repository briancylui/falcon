# Build libFalcon
# Created 17 Jun 2015 HBP - Les Houches
# ----------------------------------------------------------------------------
ifndef ROOTSYS
$(error * Define ROOTSYS)
endif

#ifndef FASTJET
#$(error * define FASTJET. Edit setup.(c)sh then source setup.(c)sh)
#endif
# ----------------------------------------------------------------------------
NAME	:= Falcon
incdir	:= include
srcdir	:= src
libdir	:= lib
bindir	:= bin
$(shell mkdir -p lib; mkdir -p tmp)

# get lists of sources
CINTSRCS:= $(wildcard $(srcdir)/*_dict.cc)

SRCS	:= $(srcdir)/FalconTester.cc
CXXSRCS	:= $(srcdir)/TKDTreeBinning.cxx $(srcdir)/TKDTree.cxx 
#OTHERSRCS:= $(filter-out $(CINTSRCS) $(SRCS),$(wildcard $(srcdir)/*.cc))

# list of dictionaries to be created
DICTIONARIES	:= $(SRCS:.cc=_dict.cc) #$(CXXSRCS:.cxx=_dict.cc)

# get list of objects
OBJECTS		:= $(OTHERSRCS:.cc=.o) $(SRCS:.cc=.o) $(DICTIONARIES:.cc=.o)
CXXOBJECTS	:= $(CXXSRCS:.cxx=.o)

#say := $(shell echo "DICTIONARIES:     $(DICTIONARIES)" >& 2)
#say := $(shell echo "" >& 2)
#say := $(shell echo "SRCS: $(SRCS)" >& 2)
#say := $(shell echo "OBJECTS: $(OBJECTS)" >& 2)
#$(error bye)
# ----------------------------------------------------------------------------
ROOTCINT	:= rootcint

# check for clang++, otherwise use g++
COMPILER	:= $(shell which clang++ >& $(PWD)/.cxx; tail $(PWD)/.cxx)
COMPILER	:= $(shell basename "$(COMPILER)")
ifeq ($(COMPILER),clang++)
CXX		:= clang++
LD		:= clang++
else
CXX		:= g++
LD		:= g++
endif
CPPFLAGS	:= -I. -I$(DELPHES) -I$(DELPHES)/external -I$(incdir)
CXXFLAGS	:= -O -Wall -fPIC -g -ansi -Wshadow -Wextra \
$(shell root-config --cflags) # -fopenmp 
LDFLAGS		:= -g 
# ----------------------------------------------------------------------------
# which operating system?
OS := $(shell uname -s)
ifeq ($(OS),Darwin)
	LDFLAGS += -dynamiclib
	LDEXT	:= .dylib
else
	LDFLAGS	+= -shared
	LDEXT	:= .so
endif	
LDFLAGS += $(shell root-config --ldflags) -L$(DELPHES) # -fopenmp 
LIBS 	:= -lDelphes -lPyROOT -lTMVA -lTMVAGui $(shell root-config --libs --nonew)
LIBRARY	:= $(libdir)/lib$(NAME)$(LDEXT)
# ----------------------------------------------------------------------------
all: $(LIBRARY)

$(LIBRARY)	: $(OBJECTS)
	@echo ""
	@echo "=> Linking shared library $@"
	$(LD) $(LDFLAGS) $^ $(LIBS)  -o $@

$(OBJECTS)	: %.o	: 	%.cc
	@echo ""
	@echo "=> Compiling $<"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(CXXOBJECTS)	: %.o	: 	%.cxx
	@echo ""
	@echo "=> Compiling $<"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(DICTIONARIES)	: $(srcdir)/%_dict.cc	: $(incdir)/%.h
	@echo ""
	@echo "=> Building dictionary $@"
	$(ROOTCINT) -f $@ -c $(CPPFLAGS) $^
	find $(srcdir) -name "*.pcm" -exec mv {} $(libdir) \;

tidy:
	rm -rf $(srcdir)/*_dict*.* $(srcdir)/*.o 

clean:
	rm -rf $(libdir)/* $(srcdir)/*_dict*.* $(srcdir)/*.o
