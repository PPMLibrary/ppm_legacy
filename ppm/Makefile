#-------------------------------------------------------------------------
#  File         :  Makefile
#-------------------------------------------------------------------------
#
#  Purpose      :  Compilation
#
#  Remarks      :
#
#  References   :
#
#  Revisions    :
#-------------------------------------------------------------------------
#  $Log: Makefile,v $
#  Revision 1.1.1.1  2007/07/13 10:16:28  ivos
#  CBL version of the PPM library
#
#  Revision 1.7  2006/08/14 14:05:12  ymishin
#  Removed fftw3.f detection from the automatic dependency files
#
#  Revision 1.6  2006/04/25 15:37:03  pchatela
#  Added directory creation
#
#  Revision 1.5  2006/04/07 07:12:35  pchatela
#  Removed the dreaded dependency entry for mpif.h
#  Cleaner removal of comment lines in sed stuff
#
#  Revision 1.4  2006/02/20 17:34:55  pchatela
#  Added flag to keep pre-processed sources
#  for debugging purposes
#
#  Revision 1.3  2005/04/29 11:17:28  pchatela
#  Additional fixes to get compilation going on G5's
#  - Deleting left over spaces in dependecy lines
#  - Moved the -P option to the beginning of the CPP options
#  - Allow for two pre-processors: one for generating the dependencies
#    (this will be typically gnu cpp) and one for the actual code-compilation
#    as some compilers, e.g. xlf, are quite picky about the code they
#    are passed
#
#  
#
#-------------------------------------------------------------------------
#  Parallel Particle Mesh Library (PPM)
#  Institute of Computational Science
#  ETH Zentrum, Hirschengraben 84
#  CH-8092 Zurich, Switzerland
#-------------------------------------------------------------------------

SHELL   := /bin/sh

include Make.gfortran

# Default value for the dependency pre-processor
# = same as code-generating pre-processor
DEPCPP ?= $(CPP)

# Name of the library
LIB_PPM := libppm.a
# Sources dir
SRC_DIR := $(BUILDDIR)/src
# Where to put .o and .d files
OBJ_DIR := $(BUILDDIR)/objects
# Target
TARGET  := $(OBJ_DIR)/$(LIB_PPM)
# Installation directories
INSTALL_DIR := $(BUILDDIR)/lib
MODULES_DIR := $(BUILDDIR)/include

# List the source files
SOURCES :=  $(notdir $(wildcard $(SRC_DIR)/ppm_module_*.f))

# Now list all the paths to check
VPATH := $(SRC_DIR):
VPATH += $(MPIINCLDIR):
VPATH += $(FFTWINCLDIR):
VPATH += $(MODULES_DIR):

# Useful derived macros
INCLS := $(patsubst %,-I%,$(subst :, ,$(VPATH)))
OBJECTS := $(SOURCES:%.f=$(OBJ_DIR)/%.o)
MODULES := $(SOURCES:%.f=$(MODULES_DIR)/%.mod)
MODSRCS := $(SOURCES:%.f=$(MODULES_DIR)/__%.f)
DEPENDENCIES := $(SOURCES:%.f=$(OBJ_DIR)/%.d)

# If you want to keep the pre-compiled sources!
ifdef KEEPCPP
RMCPP := 'touch'
else
RMCPP := 'rm'
endif

$(warning Checking for directories...)
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(MODULES_DIR))
$(shell mkdir -p $(INSTALL_DIR))


.DEFAULT: ;

# Default target is everything
all: $(TARGET)

# final linking
$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)

# How to make dependency files
# Some explanations: 
# 1) we use the preprocessor to find the includes
# 2) we add the dependency file as a target
# 3) find all the INCLUDE statements that are not inside a comment
# (actually we just check for lines starting with !)
# We exclude mpif.h from this search. 
# If make does not find mpif.h (user did not set the MPIINCLDIR correctly), 
# this could lead to an infinite loop for the dependcy files.
# 4) find all the USE statements that are not inside a comment
# 5) append a bogus line
$(OBJ_DIR)/%.d: $(SRC_DIR)/%.f
	@echo '[ quietly making' $@ ']'
	@$(DEPCPP) $(INCLS) -M $< | \
	sed -e 's#$*.o#$(OBJ_DIR)/$*.o $(OBJ_DIR)/$*.d#' \
	    -e 's#$$#\\#' -e's#\\\\#\\#' > $@
	@$(DEPCPP) -P $(DEFINE) $(INCLS) $< > $(OBJ_DIR)/__$*.f
	@grep "INCLUDE " $(OBJ_DIR)/__$*.f | \
        sed -e 's#^[ \t]*##;s#[ \t]*##' \
	    -e  '/^!/d;/mpif.h/d' \
	    -e  '/^!/d;/fftw3.f/d' \
            -e 's#INCLUDE ##;s#$$#\\#' \
            -e 's#"##g' -e "s#'##g" -e 's# ##g' >> $@
	@echo '# end of source dependencies for .o and .d files' >> $@
	@echo ''$(OBJ_DIR)/$*.o ':\' >> $@ 
	@grep "USE " $(OBJ_DIR)/__$*.f | \
	sed -e 's#^[ \t]*##;s#[ \t]*##' \
	    -e '/^!/d' \
	    -e 's#,.*##' | \
	sed -e 's#USE #$(OBJ_DIR)/#' \
	    -e 's# ##g;s#$$#.o\\#' >> $@
	@echo '# end of module dependencies for .o file' >> $@
	@rm $(OBJ_DIR)/__$*.f

# How to make .o and .mod files
# Note how we move the .mod file: we use * and not $* because
# compilers lower-case the file name, thus creating a mismatch
# in certain cases, e.g. ppm_module_neighlist_MkNeighIdx.f
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(CPP) -P $(DEFINE) $(INCLS) $< > $(OBJ_DIR)/__$*.f
	$(FC) $(INCLS) $(FFLAGS) -c -o $@ $(OBJ_DIR)/__$*.f
	@mv $(BUILDDIR)/*.mod $(MODULES_DIR)
	@$(RMCPP) $(OBJ_DIR)/__$*.f
#	@rm $(OBJ_DIR)/__$*.f

# Remove the objects so that everything recompiles
clean:
	@echo '[ Cleaning up ]'
	rm -f $(OBJECTS)
	rm -f $(MODULES)
	rm -f $(TARGET)

# Get rid of the .d's too. Notice this may produce weird behavior during
# the next invocation, since the makefile itself needs the .d's.
reallyclean:
	rm -f $(OBJECTS)
	rm -f $(MODULES)
	rm -f $(TARGET)
	rm -f $(DEPENDENCIES)

# Make new
new: clean all install

# Make install
install: all
	@echo '[ deploying to '$(INSTALL_DIR)']'
	@cp $(TARGET) $(INSTALL_DIR)
	@ranlib $(INSTALL_DIR)/$(LIB_PPM)

# Do not run this makefile until all dependency files are up-to-date
-include $(DEPENDENCIES)

