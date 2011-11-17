#------------------------------------------------------------------------------
include Makefile.arch

#------------------------------------------------------------------------------

NAME    := libsolgem

SOLINCLUDE =

#------------------------------------------------------------------------------
# Hall A Analyzer

ANALYZER ?= /dev/null

SOLINCLUDE += -I$(ANALYZER)/src -I$(ANALYZER)/hana_decode

#------------------------------------------------------------------------------
# EVIO

EVIO ?= /dev/null

SOLINCLUDE += -I$(EVIO)/src/libsrc -I$(EVIO)/src/libsrc++
LDFLAGS  += -L$(EVIO)/lib -levioxx -levio -lz -lexpat

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
