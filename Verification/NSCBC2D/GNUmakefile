# GNUmake build for the 2-D NSCBC analysis tool, for those not using CMake.
# Modelled on amrex/Tools/Plotfile/GNUmakefile: fielddump.cpp carries its own
# main(), so the multiple_executables pattern builds it directly.
#
#   make -j              -> fielddump.<machineSuffix>.ex
#   make COMP=llvm -j    -> match whatever you built PeleC with
#
# DIM must match the plotfiles being read.

PELEC_HOME  ?= ../..
AMREX_HOME  ?= $(PELEC_HOME)/Submodules/PelePhysics/Submodules/amrex

DEBUG     = FALSE
DIM       = 2
COMP      = gnu
PRECISION = DOUBLE
USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = FALSE

BL_NO_FORT = TRUE

programs += fielddump

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

multiple_executables = $(addsuffix .$(machineSuffix).ex, $(programs))
default: $(multiple_executables)

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) $(multiple_executables)
