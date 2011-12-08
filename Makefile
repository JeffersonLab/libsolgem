#------------------------------------------------------------------------------
include Makefile.arch

#------------------------------------------------------------------------------

NAME    := libsolgem

SOLINCLUDE =

#------------------------------------------------------------------------------
# Hall A Analyzer

PLATFORM = $(shell uname -s)-$(shell uname -i)
HOSTNAME = $(shell hostname)

#ifeq ($(HOSTNAME),fvm13)
ANALYZER ?= $(HOME)/ANALYZER
#else
ANALYZER ?= /dev/null
#endif

SOLINCLUDE += -I$(ANALYZER)/src -I$(ANALYZER)/hana_decode

#------------------------------------------------------------------------------
# EVIO

#ifeq ($(HOSTNAME),fvm13)
EVIO ?= ../libevio
SOLINCLUDE += -I$(EVIO)/include
LDFLAGS  += -L$(EVIO)/$(PLATFORM)/lib -levioxx -levio -lz -lexpat
#else
EVIO ?= /dev/null

SOLINCLUDE += -I$(EVIO)/src/libsrc -I$(EVIO)/src/libsrc++
LDFLAGS  += -L$(EVIO)/lib -levioxx -levio -lz -lexpat
#endif

#------------------------------------------------------------------------------

CXXFLAGS += $(SOLINCLUDE)

DICT	= $(NAME)_dict
SRC   = src/TSolAnalyzer.cxx \
        src/TSolEVIOFile.cxx \
        src/TSolEvData.cxx \
        src/TSolGEMChamber.cxx \
        src/TSolGEMCluster.cxx \
        src/TSolGEMData.cxx \
        src/TSolGEMPlane.cxx \
        src/TSolGEMVStrip.cxx \
        src/TSolSimAux.cxx \
        src/TSolSimGEMDigitization.cxx \
        src/TSolWedge.cxx \
        src/TSolSpec.cxx 
        OBJS	= $(SRC:.$(SrcSuf)=.$(ObjSuf)) $(DICT).o
HDR	= $(SRC:.$(SrcSuf)=.h) src/Linkdef.h
ROHDR	= $(SRC:.$(SrcSuf)=.h) src/Linkdef.h

LIBSOLGEM	= libsolgem.so

PROGRAMS	= $(LIBSOLGEM)

all:	$(PROGRAMS)

$(LIBSOLGEM):	$(OBJS)
	$(LD) $(SOFLAGS) $^ -o $@ $(LDFLAGS) 

clean:
	@rm -f $(OBJS) $(PROGRAMS) *dict.*

$(DICT).cxx: $(ROHDR) 
	$(ROOTCINT) -f $@ -c $(SOLINCLUDE) $^ 

%.$(ObjSuf): %.$(SrcSuf)
	$(CXX)  $(CXXFLAGS) -c -o $@ $<
