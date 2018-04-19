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

MAIN = check1 check2 time1 time2

RANDOM = ranlxs ranlxd gauss

UTILS = utils

LINALG = cmatrix cmatrix_dble

MODULES = $(RANDOM) $(UTILS) $(LINALG)


# search path for modules

MDIR = ../../../modules

VPATH = $(MDIR)/nompi/utils:$(MDIR)/random:$(MDIR)/linalg


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
	@ -rm -rf *.d *.o $(MAIN)
.PHONY: clean

################################################################################
