#======================================================================
#
#   file:               Makefile
#
#   purpose:            Makefile for compiling, running and testing
#                       the BEPS model (as provided by Mousong Wu),
#                       it's derivative codes and optimisation experiments.
#
#   created:            2020/02
#   last:               2021/04
#
#   author:             Michael Vossbeck (The Inversion Lab)
#
#======================================================================

#
RM := rm
# Set path to Mkdepends script; assumes that any Makefile including
# this file is in a sibling of the src directory, in which Mkdepends
# resides
Mkdepends := ./Mkdepends

# Set up special characters
null  :=

# Newer makes set the CURDIR variable.
CURDIR := $(shell pwd)

# Determine platform
UNAMES := $(shell uname -s)

HOST := $(shell hostname)

#-- derivative code(s)
FWEXT    := _fw
BWEXT    := _bw
FWVEXT   := _fwv
DVFWSRC  := bepsderiv$(FWEXT).f90
DVBWSRC  := bepsderiv$(BWEXT).f90
DVFWVSRC := bepsderiv$(FWVEXT).f90
DVFWVMOD := diffsizes.f90

# platform: invoke with PLATFORM=iLab to run package on iLab computer
PLATFORM =
ifeq ($(HOST),ilabdesk)
  PLATFORM = iLab
else ifeq ($(HOST),ilabpower)
  PLATFORM = iLab
else ifeq ($(HOST),mvobook)
  PLATFORM = iLab
else ifeq ($(HOST),txkbook2)
  PLATFORM = iLab
endif

# Check for the netcdf library and include directories 
#  NETCDFDIR  := /software/netcdf-3.6.3-intel
#  NETCDFDIR  := /data/home/mousong/BEPS_folder/software
FC         := gfortran
CC         := gcc
ifeq ($(PLATFORM),iLab)
  ifeq ($(FC),gfortran)
    LIB_NETCDF := /public1/home/sc30439/soft/new-netcdf/install/
    INC_NETCDF := /public1/home/sc30439/soft/new-netcdf/install//include
    MOD_NETCDF := /public1/home/sc30439/soft/new-netcdf/install/include
    LDFLAGS     = $(shell $(LIB_NETCDF)/bin/nf-config --flibs)
    LDFLAGS     = $(LIB_NETCDF)/lib64/libnetcdff.a -lnetcdf -ldl -lm
  else ifeq ($(findstring ifort,$(FC)),ifort)
    LIB_NETCDF := /public1/home/sc30439/soft/new-netcdf/install/
    INC_NETCDF := /public1/home/sc30439/soft/new-netcdf/install//include
    MOD_NETCDF := /public1/home/sc30439/soft/new-netcdf/install//include
    LDFLAGS     = $(shell $(LIB_NETCDF)/bin/nf-config --flibs)
    LDFLAGS     = $(LIB_NETCDF)/lib/libnetcdff.a -lnetcdf -ldl -lm
  endif
  MPIDIR     := /public1/soft/intel/2017/compilers_and_libraries_2017.7.259/linux/mpi
  INC_MPI    := $(MPIDIR)/intel64/include/
  LIB_MPI    := $(MPIDIR)/intel64/lib
else
  LIB_NETCDF := /home/hl/netcdflib/lib
  INC_NETCDF := /home/hl/netcdflib/include
  MOD_NETCDF := /home/hl/netcdflib/include
  MPIDIR     := /home/hl/netcdflib
  INC_MPI    := $(MPIDIR)/include
  LIB_MPI    := $(MPIDIR)/lib
  LDFLAGS     = $(shell $(LIB_NETCDF)/../bin/nf-config --flibs)

##  FC         := mpiifort
##  CC         := mpicc
##  NETCDFDIR  := /public1/soft/netcdf/4.4.1-parallel-icc18
##  LIB_NETCDF := $(NETCDFDIR)/lib
##  INC_NETCDF := $(NETCDFDIR)/include
##  MOD_NETCDF := $(NETCDFDOR)/include
##  MPIDIR     := /public1/soft/intel/2018/compilers_and_libr
##  INC_MPI    := $(MPIDIR)/include
##  LIB_MPI    := $(MPIDIR)/lib
  # NETCDFDIR  := /public1/home/sc30439/zhaoyq-netcdf/install 
  # LIB_NETCDF := /public1/home/sc30439/zhaoyq-netcdf/install/lib
  # INC_NETCDF := /public1/home/sc30439/zhaoyq-netcdf/install/include
  # MOD_NETCDF := /public1/home/sc30439/zhaoyq-netcdf/install/include
  # MPIDIR     :=/public1/soft/mpich2/1.0.8
  # INC_MPI    :=$(MPIDIR)/include
  # LIB_MPI    :=$(MPIDIR)/lib
#   NETCDFDIR  := /public1/home/sc30439/zhaoyq-netcdf/install
#   LIB_NETCDF := $(NETCDFDIR)/lib
#   INC_NETCDF := $(NETCDFDIR)/include
#   MOD_NETCDF := $(LIB_NETCDF)

# # Setting MPI
# #  MPIDIR     :=/share/apps/intel/impi/2017.1.132
# #  MPIDIR     :=/software/mpich2-1.5-intel
#   MPIDIR     :=/public1/soft/mpich2/1.0.8
#   INC_MPI    :=$(MPIDIR)/include
#   LIB_MPI    :=$(MPIDIR)/lib
endif

# setting esmf
# ESMFDIR    :=./esmf_wrf_timemgr/unittests
ESMFDIR    :=


# Set if Shared memory multi-processing will be used
ifeq ($(SMP),$(null))
  SMP := FALSE
endif

CPPDEF := cpp

# Set optimization on by default
ifeq ($(OPT),$(null))
  OPT := TRUE
endif

ifeq ($(OPT),TRUE)
  CPPDEF := -DOPT
endif

# Load dependency search path.
dirs := . $(shell cat Filepath)

# Set cpp search path, include netcdf
#-- iLab::want to have unique entries
# cpp_dirs := $(dirs) $(INC_NETCDF) $(MOD_NETCDF) $(INC_MPI)  # $(LIB_MPI)
cpp_dirs := $(shell echo $(dirs) $(INC_NETCDF) $(MOD_NETCDF) $(INC_MPI) | awk -v RS='[[:space:]]+' '!a[$$0]++{printf "%s%s", $$0, RT}')
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Expand any tildes in directory names. Change spaces to colons.
# (the vpath itself is set elsewhere, based on this variable)
vpath_dirs    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir))) 
vpath_dirs    := $(subst $(space),:,$(vpath_dirs)) 

#-- get list of BEPS core source files (already preprocessed, i.e. .f90 files)
BEPS_A_F90_SOURCES := \
aerodynamic_conductance.f90 \
AnGsMod.f90 \
beps_con.f90 \
cos_plant.f90 \
cos_grnd.f90 \
beps_helper.f90 \
beps_par.f90 \
beps_phenology.f90 \
beps_cropMod.f90 \
beps_soilMod.f90 \
beps_time_manager.f90 \
bepstype.f90 \
controlInput_mod.f90 \
bepstypeInit.f90 \
ecoRespMod.f90 \
endrun.f90 \
ESMF_AlarmClockMod.f90 \
ESMF_AlarmMod.f90 \
ESMF_BaseMod.f90 \
ESMF_BaseTimeMod.f90 \
ESMF_CalendarMod.f90 \
ESMF_ClockMod.f90 \
ESMF.f90 \
ESMF_FractionMod.f90 \
ESMF_ShrTimeMod.f90 \
ESMF_Stubs.f90 \
ESMF_TimeIntervalMod.f90 \
ESMF_TimeMod.f90 \
get_CO2.f90 \
get_COS.f90 \
inter_prg.f90 \
MeatMod.f90 \
meteoMod.f90 \
mid_results.f90 \
netRadiationMod.f90 \
outputMod.f90 \
PF_run.f90 \
MC_run.f90 \
rainsnowMod.f90 \
readcoef.f90 \
readparam.f90 \
restart.f90 \
retrive_soilp.f90 \
s_coszs.f90 \
SensibleHeat.f90 \
shr_kind_mod.f90 \
wrf_error_fatal.f90 \
wrf_message.f90


#-- add differentiation/assimilation framework related sources
BEPSFUNC_SRCS := bepsfunc_setup.f90 bepsfunc.f90 mo_bepsfunc_ctl.f90 mo_prior.f90 misfit.f90 cost.f90 obs.f90 devprior.f90 prior.f90
#-- active sources
ASRCS := $(BEPS_A_F90_SOURCES) $(BEPSFUNC_SRCS)
#-- passive sources
PSRCS = 
#-- driver sources
DRV_SOURCES +=  driver.f90 runsimobs.f90 rundfprv.f90 runderivfwv.f90

#-- list of (Fortran) files that require preprocessing
PP_SOURCES := $(shell ls *.F90)



# Architecture-specific flags and rules
#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------
DBG = 1

ifeq ($(UNAMES),Linux)
  CPPDEF  += -DLINUX -DFORTRANUNDERSCORE
  CFLAGS  := $(CPPDEF)
  #-- define basic Fortran compiler options
  ifeq ($(FC),gfortran)
    FFLAGS  := -Wall -Wno-unused -Wno-unused-parameter -Wno-unused-dummy-argument -fdefault-real-8
    FFOPTI  = -O2
  else ifeq ($(findstring ifort,$(FC)),ifort)
    FFLAGS = -warn all -warn noexternals -real-size 64
    FFOPTI = -O2
  else ifeq ($(FC),mpiifort)
    FFLAGS  := -warn all,nounused -real-size 64
    FFOPTI = -O2
  endif

  #-- put here optimisation specific options
  ifeq ($(DBG),1)
    ifeq ($(FC),gfortran)
      FFLAGS += -g -O0 -fcheck=all -fbacktrace -finit-integer=-99999 -finit-real=snan  -ffpe-trap=invalid,zero,overflow,underflow
    else ifeq ($(findstring ifort,$(FC)),ifort)
      FFLAGS += -O0 -g -check all -check noarg_temp_created -traceback -init=arrays,snan
    else ifeq ($(FC),mpiifort)
      FFLAGS += -g -O0
    endif
  else
 #  FFLAGS += -O2 -c -r8 -i4 -132 -convert big_endian -assume byterecl
    FFLAGS += $(FFOPTI)
  endif

  CFLAGS      += -m64 -g
  LDFLAGS     += -m64

  # FFLAGS += -I$(INC_NETCDF) -I$(INC_MPI) -I$(ESMFDIR) -L$(LIB_NETCDF) -L$(LIB_MPI) -L$(ESMFDIR) \
  #            $(CPPDEF) $(cpp_path) -lnetcdff -lnetcdf -lesmf_time
  ifneq ($(ESMFDIR),$(null))
    FFLAGS += -I$(ESMFDIR)
    LDFAGS += -L$(ESMFDIR) -lesmf_time
  endif
  FFLAGS  += $(cpp_path)
  CFLAGS  += $(cpp_path)
  LDFLAGS += -L$(LIB_NETCDF) -L$(LIB_MPI) \
              $(CPPDEF) $(cpp_path) -lnetcdff -lnetcdf
endif


#------------------------------------------------------------------------
# finalise path/compile settings
#------------------------------------------------------------------------
# Set the vpath for all file types EXCEPT .o
# We do this for individual file types rather than generally using
# VPATH, because for .o files, we don't want to use files from a
# different build (e.g., in building the unit tester, we don't want to
# use .o files from the main build)
vpath %.F90 $(vpath_dirs)
vpath %.f90 $(vpath_dirs)
vpath %.c   $(vpath_dirs)
vpath %.h   $(vpath_dirs)

# Append user defined compiler and load flags to Makefile defaults
CFLAGS   += $(USER_CFLAGS)
FFLAGS   += $(USER_FFLAGS)
LDFLAGS  += $(USER_LDFLAGS)

# Set user specified linker
ifneq ($(USER_LINKER),$(null))
  LINKER := $(USER_LINKER)
else
  LINKER := $(FC)
endif




#------------------------------------------------------------------------
#                   n u m e r i c a l   r e c i p e s   l i b r a r y
#------------------------------------------------------------------------
NUMRECDIR := nr
NUMRECLIB := libnr-$(FC).a
$(NUMRECLIB): 
	$(MAKE) -C $(NUMRECDIR) ../$@ FC=$(FC) FFLAGS="$(FFLAGS) $(FFXTRA)"
nrlib: $(NUMRECLIB)
.PHONY:nrlib


#------------------------------------------------------------------------
#                   A D   s u p p o r t   l i b r a r y
#------------------------------------------------------------------------
AD_SUPPORTDIR  = adsupportlib
ifndef AD_SUPPORTLIB
  AD_SUPPORTLIB  = libadsupport.a
endif
AD_SUPPORTOBJS = adBuffer.o adStack.o
adBuffer.o: $(AD_SUPPORTDIR)/adBuffer.f
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<
adStack.o: $(AD_SUPPORTDIR)/adStack.c
	$(CC) -fPIC -I$(AD_SUPPORTDIR) $(CFLAGS) $(CFXTRA) -c -o $@ $<

$(AD_SUPPORTLIB): $(AD_SUPPORTOBJS)
	ar rcs $(AD_SUPPORTLIB) $^
	ranlib $(AD_SUPPORTLIB)
adsuplib: $(AD_SUPPORTLIB)
.PHONY: adsuplib
SCRATCH_LST += $(AD_SUPPORTLIB)


#------------------------------------------------------------------------
#                   b u i l d i n g   d r i v e r s
#------------------------------------------------------------------------
runbepsorg.x: driver.o $(addsuffix .o, $(basename $(BEPS_A_F90_SOURCES)))
	$(LINKER) -o $@ $^ $(LDFLAGS)
SCRATCH_LST += runbepsorg.x

runsimobs.x: runsimobs.o $(addsuffix .o, $(basename $(ASRCS) $(PSRCS)))
	$(LINKER) -o $@ $^ $(LDFLAGS)
SCRATCH_LST += runsimobs.x

rundfprv.x: $(addsuffix .o, $(basename $(ASRCS) $(PSRCS2))) $(DVBWSRC:.f90=.o) rundfprv.o $(AD_SUPPORTLIB) $(NUMRECLIB)
	$(LINKER) -o $@ $^ $(LDFLAGS)
SCRATCH_LST += rundfprv.x

runderivfwv.x: runderivfwv.o $(addsuffix .o, $(basename $(ASRCS) $(PSRCS))) $(DVFWVMOD:.f90=.o) $(DVFWVSRC:.f90=.o)
	$(LINKER) -o $@ $^ $(LDFLAGS)
SCRATCH_LST += runderivfwv.x


#------------------------------------------------------------------------
#                   s p e c i a l   p r e p r o c e s s i n g
#------------------------------------------------------------------------
runderivfwv.f90: runderiv.F90
	$(FC) -cpp -E -P -DDERIV_FWV $(cpp_path) $(CPPDEF) $(CPPXTRA) $< > $@



#------------------------------------------------------------------------
#                   d e p e n d e n c y   b u i l d i n g
#------------------------------------------------------------------------
$(CURDIR)/Srcfiles: runderivfwv.f90 $(PP_SOURCES:.F90=.f90)
	@$(shell ls -1 *.f90 > $@)
DEPFILE := $(CURDIR)/Depends
# dep: $(CURDIR)/Filepath $(ALL_SOURCES)
#-- iLab::remove circular dependency for obs.o generated by Mkdepends utility
dep: $(CURDIR)/Filepath $(CURDIR)/Srcfiles
	$(Mkdepends) $(CURDIR)/Filepath $(CURDIR)/Srcfiles > $(DEPFILE)
	sed -i 's/^obs.o.*/obs.o : obs.f90/g' $(DEPFILE)
	@$(RM) -f $(CURDIR)/Srcfiles
.PHONY: dep


#------------------------------------------------------------------------
#                   c l e a n i n g
#------------------------------------------------------------------------
clean:
	$(RM) -f *~ *.o *.mod

scratch: clean
	$(RM) -f $(DEPFILE) $(BEPS_A_F90_SOURCES) $(SCRATCH_LST)


#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------
.SUFFIXES:
.SUFFIXES: .F90 .f90 .c .o


#-- ILAB-NOTE::preprocessing syntax fixed to gfortran (for now)
.F90.f90:
	gfortran -cpp -E -P $(cpp_path) $(CPPDEF) $(CPPXTRA) $< > $@
	sed -i 's/public operator/public :: operator/g' $@

.f90.o:
	$(FC) -c $(FFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<

-include $(CURDIR)/Depends
