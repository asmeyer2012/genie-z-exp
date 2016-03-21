#
# Makefile for the GENIE event reweighting applications 
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gRwght1Scan gRwghtZExpAxFF gRwghtZExpDirect gRwghtCovariance

all: $(TGT)

gRwght1Scan: FORCE
	$(CXX) $(CXXFLAGS) -c gRwght1Scan.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRwght1Scan.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/grwght1scan

gRwghtZExpAxFF: FORCE
	$(CXX) $(CXXFLAGS) -c gRwghtZExpAxFF.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRwghtZExpAxFF.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/grwghtzexpaxff

gRwghtZExpDirect: FORCE
	$(CXX) $(CXXFLAGS) -c gRwghtZExpDirect.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRwghtZExpDirect.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/grwghtzexpdirect

gRwghtCovariance: FORCE
	$(CXX) $(CXXFLAGS) -c gRwghtCovariance.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRwghtCovariance.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/grwghtcov

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(GENIE_BIN_PATH)/grwght1scan
	$(GENIE_BIN_PATH)/grwghtzexpaxff
	$(GENIE_BIN_PATH)/grwghtzexpdirect
	$(GENIE_BIN_PATH)/grwghtcov

distclean: FORCE
	$(RM)
	$(GENIE_BIN_INSTALLATION_PATH)/grwght1scan
	$(GENIE_BIN_INSTALLATION_PATH)/grwghtzexpaxff
	$(GENIE_BIN_INSTALLATION_PATH)/grwghtzexpdirect
	$(GENIE_BIN_INSTALLATION_PATH)/grwghtcov

FORCE:

# DO NOT DELETE
