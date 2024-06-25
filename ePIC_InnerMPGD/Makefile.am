AUTOMAKE_OPTIONS = foreign

# list of shared libraries to produce
lib_LTLIBRARIES = \
	libePIC_InnerMPGD_Detector.la

AM_CPPFLAGS = \
	-I$(includedir) \
	-I$(OFFLINE_MAIN)/include \
	-I$(ROOTSYS)/include \
	-I$(G4_MAIN)/include \
	-I$(OPT_SPHENIX)/include 

AM_LDFLAGS = \
	-L$(libdir) \
	-L$(OFFLINE_MAIN)/lib 

libePIC_InnerMPGD_Detector_la_LIBADD = \
	-lSubsysReco \
	-lg4detectors \
	-lg4testbench

pkginclude_HEADERS = \
	ePIC_InnerMPGD_Subsystem.h

libePIC_InnerMPGD_Detector_la_SOURCES = \
	ePIC_InnerMPGD_Detector.cc \
	ePIC_InnerMPGD_Subsystem.cc  \
	ePIC_InnerMPGD_SteppingAction.cc


###################################################
# linking test

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
	testexternals_ePIC_InnerMPGD

testexternals_ePIC_InnerMPGD_SOURCES = testexternals.cc
testexternals_ePIC_InnerMPGD_LDADD = libePIC_InnerMPGD_Detector.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >>$@
	echo " return 0;" >> $@
	echo "}" >> $@


clean-local:
	rm -f $(BUILT_SOURCES)



