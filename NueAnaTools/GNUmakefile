# Makefile for NueAna Display

# Run rootcint

#ROOTCINTFLAGS := -p
ROOTCINT:=YES
                                                                                
vpath NueAna/NueAnaTools/DCGraph.h $(SRT_PRIVATE_CONTEXT)/include/ \
                                        $(SRT_PUBLIC_CONTEXT)/include/
CINTLIST += NueAna/NueAnaTools/DCGraph.h
                                                                                
vpath NueAna/NueAnaTools/DCEdge.h $(SRT_PRIVATE_CONTEXT)/include/ \
                                        $(SRT_PUBLIC_CONTEXT)/include/
CINTLIST += NueAna/NueAnaTools/DCEdge.h
                                                                                
vpath NueAna/NueAnaTools/DCVertex.h $(SRT_PRIVATE_CONTEXT)/include/ \
                                        $(SRT_PUBLIC_CONTEXT)/include/
CINTLIST += NueAna/NueAnaTools/DCVertex.h
                                                                                
vpath NueAna/NueAnaTools/DCHit.h $(SRT_PRIVATE_CONTEXT)/include/ \
                                        $(SRT_PUBLIC_CONTEXT)/include/
CINTLIST += NueAna/NueAnaTools/DCHit.h
                                                                                


#LIB:=lib$(PACKAGE)Display
PACKAGE:=NueAnaTools
LIB:=lib$(PACKAGE)

# library contents

LIBCXXFILES:=$(wildcard *.cxx)

############

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_sigc++.mk



