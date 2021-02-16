#
# GNUmakefile for simple package
#

SUBDIRS = Plotter Tools

ROOTCINT := YES

PACKAGE := NueAnaParticleAna

LIB := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
override LIBLIBS = -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lNueAna
########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
#


