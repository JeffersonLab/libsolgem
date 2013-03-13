#------------------------------------------------------------------------------
include Makefile.arch

#------------------------------------------------------------------------------

NAME    := libsolgem

SOLINCLUDE := -I$(shell pwd)/src

#------------------------------------------------------------------------------
# Hall A Analyzer

# Analyzer default location,
ANALYZER ?= $(HOME)/ANALYZER
# Possible Analyzer header locations, will be used in the order found
ANAINCDIRS  := $(wildcard $(addprefix $(ANALYZER)/, include src hana_decode hana_scaler))
ifeq ($(strip $(ANAINCDIRS)),)
  $(error No Analyzer header files found. Check $$ANALYZER)
endif


#------------------------------------------------------------------------------
# EVIO

PLATFORM = $(shell uname -s)-$(shell uname -i)
# EVIO default location, as a last resort if $EVIO isn't set in build env
EVIO ?= ../libevio
# Possible EVIO header directories, will be used in the order found
EVIOINC := $(wildcard $(addprefix $(EVIO)/, include src/libsrc src/libsrc++))
# Possible EVIO library locations, the first one found will be used
EVIOLIB := $(firstword $(wildcard $(addprefix $(EVIO)/, $(PLATFORM)/lib lib)))
ifeq ($(strip $(EVIOINC)),)
  $(error No EVIO header files found. Check $$EVIO)
endif
ifeq ($(strip $(EVIOLIB)),)
  $(error No EVIO library directory found. Check $$EVIO)
endif
ifeq (debug,$(findstring debug,$(ROOTBUILD)))
  DBGSFX = _d
endif
SOLINCLUDE += $(addprefix -I, $(EVIOINC) )
LDFLAGS  += -L$(EVIOLIB) -levioxx$(DBGSFX) -levio$(DBGSFX) -lz -lexpat

# Some of the analyzer include dirs conflict with headers in
# EVIO
SOLINCLUDE += $(addprefix -I, $(ANAINCDIRS) )

#------------------------------------------------------------------------------

CXXFLAGS += $(SOLINCLUDE)

DICT	= $(NAME)_dict
SRC   = src/TSolAnalyzer.cxx \
	src/TSolClusters.cxx \
        src/TSolEVIOFile.cxx \
        src/TSolGEMChamber.cxx \
        src/TSolGEMCluster.cxx \
        src/TSolGEMData.cxx \
        src/TSolGEMPlane.cxx \
        src/TSolGEMVStrip.cxx \
        src/TSolSimAux.cxx \
        src/TSolSimGEMDigitization.cxx \
        src/TSolWedge.cxx \
        src/TSolSpec.cxx \
        src/TSolSimEvent.cxx \
        src/TSolSimFile.cxx \
        src/TSolSimDecoder.cxx

OBJS	= $(SRC:.cxx=.$(ObjSuf)) $(DICT).o
HDR	= $(SRC:.cxx=.h) src/Linkdef.h
ROHDR	= $(SRC:.cxx=.h) src/Linkdef.h

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
