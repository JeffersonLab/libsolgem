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
LDFLAGS  += -L$(EVIO)/lib -levio -levioxx -lz -lexpat

#------------------------------------------------------------------------------

CXXFLAGS += $(SOLINCLUDE)

DICT	= $(NAME)_dict
SRC	= src/TSolGEMCluster.cxx src/TSolGEMPlane.cxx src/TSolSpec.cxx \
	  src/TSolAnalyzer.cxx
OBJS	= $(SRC:.$(SrcSuf)=.$(ObjSuf)) $(DICT).o
HDR	= $(SRC:.$(SrcSuf)=.h) src/Linkdef.h

LIBSOLGEM	= libsolgem.so

PROGRAMS	= $(LIBSOLGEM)

all:	$(PROGRAMS)

$(LIBSOLGEM):	$(OBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@ 

clean:
	@rm -f $(OBJS) $(PROGRAMS) *dict.*

$(DICT).cxx: $(HDR) 
	$(ROOTCINT) -f $@ -c $(SOLINCLUDE) $^ 

%.$(ObjSuf): %.$(SrcSuf)
	$(CXX)  $(CXXFLAGS) -c -o $@ $<
