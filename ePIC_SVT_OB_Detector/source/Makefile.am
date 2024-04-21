AUTOMAKE_OPTIONS = foreign

# list of shared libraries to produce
lib_LTLIBRARIES = \
	libePIC_SVT_OB_Detector.la

AM_CPPFLAGS = \
	-I$(includedir) \
	-I$(OFFLINE_MAIN)/include \
	-I$(ROOTSYS)/include \
	-I$(G4_MAIN)/include \
	-I$(OPT_SPHENIX)/include 

AM_LDFLAGS = \
	-L$(libdir) \
	-L$(OFFLINE_MAIN)/lib 

libePIC_SVT_OB_Detector_la_LIBADD = \
	-lSubsysReco \
	-lg4detectors \
	-lg4testbench

pkginclude_HEADERS = \
	ePIC_SVT_OB_Subsystem.h

libePIC_SVT_OB_Detector_la_SOURCES = \
	ePIC_SVT_OB_Detector.cc \
	ePIC_SVT_OB_Subsystem.cc  \
	ePIC_SVT_OB_SteppingAction.cc


###################################################
# linking test

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
	testexternals_ePIC_SVT_OB

testexternals_ePIC_SVT_OB_SOURCES = testexternals.cc
testexternals_ePIC_SVT_OB_LDADD = libePIC_SVT_OB_Detector.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >>$@
	echo " return 0;" >> $@
	echo "}" >> $@


clean-local:
	rm -f $(BUILT_SOURCES)



