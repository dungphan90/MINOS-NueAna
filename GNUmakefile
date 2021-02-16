# Makefile for NueAna

SUBDIRS = Module Reweight NueAnaTools Display Extrapolation Sensitivity ParticlePID MultiBinAna

 
# Run rootcint

#ROOTCINT:=NO
ROOTCINT:=YES

# library

LIB:=lib$(PACKAGE)

# library contents

LIBCXXFILES:=$(wildcard *.cxx)


############

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
#include SoftRelTools/arch_spec_sigc++.mk
