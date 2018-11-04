all: DigiPass ReplayMCDig

ifdef DEBUG
  CXXFLAGS    = -g -O0 -Wall -Wextra
else
  CXXFLAGS    = -g -O2 -Wall -Wextra
endif

ifndef ANALYZER
  $(error $$ANALYZER environment variable not defined)
endif

INCDIRS  = $(wildcard $(addprefix $(ANALYZER)/, include src hana_decode))
#INCDIRS     = ../src ../hana_decode
INCDIRS  +=${LIBSBSGEM}/src
INCDIRS  +=${TREESEARCH}
INCLUDES    = $(addprefix -I, $(INCDIRS))

CXX         = $(shell root-config --cxx)
ROOTCFLAGS  = $(shell root-config --cflags)
LD          = $(shell root-config --ld)
LDFLAGS     = $(shell root-config --ldflags)

CXXFLAGS   += $(ROOTCFLAGS) $(INCLUDES)
#ALLINCLUDES  = -I$(shell root-config --incdir) $(INCLUDES)
LIBS        = $(shell root-config --libs) -L${ANALYZER} -lHallA -ldc

# Add EVIO lib, needed by libdc
ifndef EVIO_LIBDIR
  EVIO_LIBDIR = ${ANALYZER}
endif
LIBS += -L$(EVIO_LIBDIR) -levio

# Add EVIO lib, needed by libdc
ifndef DIG_LIBDIR
  DIG_LIBDIR = ${LIBSBSGEM}
endif
LIBS += -L$(DIG_LIBDIR) -lEG -lsolgem

ifndef TS_LIBDIR
  TS_LIBDIR = ${TREESEARCH}
endif
LIBS += -L$(TS_LIBDIR) -lMinuit -lTreeSearch -lTreeSearch -lTreeSearch-SBS


#dbconvert:	dbconvert.o
#		$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

DigiPass:	DigiPass.o
		$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

ReplayMCDig:	ReplayMCDig.o
		$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

clean:
		rm -f DigiPass DigiPass.o
		rm -f ReplayMCDig ReplayMCDig.o

%.o:		%.cxx Makefile
		$(CXX) $(CXXFLAGS) -o $@ -c $<

.PHONY: all clean
