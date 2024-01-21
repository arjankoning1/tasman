# This file is part of TASMAN
# 
#    Author: Arjan Koning Copyright: (C) 2023
#
# TASMAN is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.
#
# Basic installation (verified with the sample cases)
#
# FC="gfortran " FFLAGS=" "
#
# Distribution FC & FFLAGS (options provided by J-C Sublet)
#
# FC="gfortran " FFLAGS=" -Ofast "
# FC="ifort    " FFLAGS=" -Ofast "
# FC="nagfor   " FFLAGS=" -w     "
#
# Development FC & FFLAGS  (options provided by J-C Sublet)
#
# FC="gfortran " FFLAGS=" -Wall -fcheck=all -Og -g -fbacktrace   "
# FC="ifort    " FFLAGS=" -O0 -g -traceback -check all -debug all"
# FC="nagfor   " FFLAGS=" -C=all -O0 -g -gline                   "
#
CWD = $(shell pwd)
PREFIX ?= $(CWD)/../

FC = gfortran
FFLAGS = -w

fsub = $(shell echo [A-z]*.f90)
fsub1 = $(shell echo [A-z]*.f)
obj = $(fsub:.f=.o) $(fsub1:.f=.o)
src = $(fsub) $(fsub1)

all: change tasman install clean

change : 
	../path_change

tasman : $(obj)
	$(FC) $(FFLAGS) -o  tasman  $(obj)

install:
	mkdir -p $(PREFIX)bin
	mv tasman $(PREFIX)bin  

clean:
	rm -f *.o *.mod

