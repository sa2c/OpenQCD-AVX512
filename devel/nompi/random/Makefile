################################################################################
#
# Makefile to compile and link C programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
#
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq
.PHONY: all


# main programs and required modules

MAIN = check1 check2 check3 check4 check5 check6 check7\
       time1 time2

RANDOM = ranlxs ranlxd gauss random_su3

UTILS = utils

SU3FCTS = su3prod su3ren

EXTRAS = chebyshev ks_test pchi_square stat

MODULES = $(RANDOM) $(UTILS) $(SU3FCTS) $(EXTRAS)


# search path for modules

MDIR = ../../../modules

VPATH = $(MDIR)/nompi/extras:$(MDIR)/nompi/utils:\
        $(MDIR)/random:$(MDIR)/su3fcts


# additional include directories

INCPATH = ../../../include/nompi ../../../include


# additional libraries to be included

LIBS = m

LIBPATH =


# scheduling and optimization options

CFLAGS = -std=c89 -pedantic -fstrict-aliasing -fno-inline \
         -Wall -Wno-long-long -Wstrict-prototypes -Werror \
         -O -mno-avx -Dx64 -DPM

LFLAGS =


############################## do not change ###################################

SHELL=/bin/bash
CC=$(GCC)

PGMS= $(MAIN) $(MODULES)

INCDIRS = $(addprefix -I,$(INCPATH))

OBJECTS = $(addsuffix .o,$(MODULES))

LDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS))

-include $(addsuffix .d,$(PGMS))


# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.c Makefile
	@ $(GCC) -MM -ansi $(INCDIRS) $< -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.c Makefile
	$(CC) $< -c $(CFLAGS) $(INCDIRS) -o $@


# rule to link object files

$(MAIN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(LFLAGS) $(LDFLAGS) -o $@


# produce executables

mkxeq: $(MAIN)


# remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# clean directory

clean:
	@ -rm -rf *.d *.o .tmp $(MAIN)
.PHONY: clean

################################################################################
