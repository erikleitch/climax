#-----------------------------------------------------------------------
# Include Climax-specific variables
#-----------------------------------------------------------------------

SHELL = /bin/sh

include Makefile.defs

#-----------------------------------------------------------------------
# Targets
#-----------------------------------------------------------------------

all: make_directives make_dirs libs bins msg

MAKE_OBJS = util fftutil pgutil program cfitsio slalib models datasets

#-----------------------------------------------------------------------                                 
# Add optional targets, depending on compile directives
#-----------------------------------------------------------------------

ifneq ($(MATLAB_PATH),)
  MAKE_OBJS += matlab
endif

ifneq ($(PYTHON_INC_PATH),)
  MAKE_OBJS += python
endif

# Rule for making dox documentation

dox:
	@for dir in $(MAKE_OBJS) ; do \
	echo ''; echo 'Making dox in: ' $$dir; echo '';\
	  if [ -d $$dir ] ; then (cd $$dir; $(MAKE) dox); fi ; \
	done

# Rule for making libraries

libs:
	@for dir in $(MAKE_OBJS) ; do \
	  echo '';\
	  echo '=======================================================================';\
	  echo '                        Making libs in: ' $$dir; \
	  echo '=======================================================================';\
	  echo '';\
	  if [ -d $$dir ] ; then (cd $$dir; $(MAKE) libs); fi ; \
	done

# Rule for making binaries

bins:
	@for dir in $(MAKE_OBJS) ; do \
	  echo '';\
	  echo '=======================================================================';\
	  echo '                        Making bins in: ' $$dir; \
	  echo '=======================================================================';\
	  echo '';\
	  if [ -d $$dir ] ; then (cd $$dir; $(MAKE) bins); fi ; \
	done

# Test programs

test:
	@for dir in $(MAKE_OBJS) ; do \
	  echo 'Making test programs in: ' $$dir ; \
	  if [ -d $$dir/Test ] ; then (cd $$dir/Test; $(MAKE) test); fi ; \
	done

# Clean directives

clean_obj:
	@for dir in $(MAKE_OBJS) ; do \
	  echo ''; echo 'Cleaning: ' $$dir; echo ''; \
	  if [ -d $$dir ] ; then (cd $$dir; $(MAKE) clean); fi ; \
	  if [ -d $$dir/Test ] ; then (cd $$dir/Test; $(MAKE) clean_test); fi ; \
	done

clean_depend:
	@for dir in $(MAKE_OBJS) ; do \
	  echo ''; echo 'Cleaning dependency files in: ' $$dir ; echo ''; \
	  if [ -d $$dir ] ; then (cd $$dir; $(MAKE) clean_depend); fi ; \
	  if [ -d $$dir/Test ] ; then (cd $$dir/Test; $(MAKE) clean_depend); fi ; \
	done

clean_dirs:
	@if [ -d $(LIBDIR) ] ; then \rm -fr $(LIBDIR) ; fi ;
	@if [ -d $(BINDIR) ] ; then \rm -fr $(BINDIR) ; fi ;
	@if [ -d gcp ] ; then  \rm gcp; fi ;

clean: clean_obj clean_depend clean_dirs

# If the directives don't exist, create the directives file from the template

make_make_directives:
	@if [ ! -f Makefile.directives ] ; then cp Makefile.directives.template Makefile.directives ; fi ;

# If any directives have changed since the last compile, touch relevant files

make_directives:
	$(MAKE) make_make_directives
ifeq ($(DIR_HAVE_CHANGED), y)
	@echo ""
	@echo "Directives have changed since the last compile: touching Directives.h"
	@echo ""
	@touch util/Directives.h
	@echo "OLD_COMPILE_WITH_DEBUG = " $(COMPILE_WITH_DEBUG)  >> Makefile.directives.last
	$(MAKE) clean_depend
endif

# Make directories needed for compilation

make_dirs:
	@if [ ! -d $(LIBDIR) ]  ; then  mkdir $(LIBDIR) ; fi ;
	@if [ ! -d $(BINDIR) ]  ; then  mkdir -p $(BINDIR);  fi; 
	@if [ ! -d $(HELPDIR) ] ; then  mkdir -p $(HELPDIR);  fi; 
	@if [ ! -d gcp ] ; then  ln -s . gcp; fi;


ifeq ($(MAC_OSX), 1)
msg:
	@echo ""
	@echo "======================================================================="
	@echo "#                                                                      "
	@echo "# Finished compiling                                                   "
	@echo "#                                                                      "
	@echo "# Before you attempt to run any executable in this code tree,          "
	@echo "# make sure you add                                                    "
	@echo "#                                                                      "
	@echo "#            $(TOP)/lib                                                "
	@echo "#                                                                      "
	@echo "# and                                                                  "
	@echo "#                                                                      "
	@echo "#            $(TOOLSDIR)/lib                                           "
	@echo "#                                                                      "
	@echo "# to your LD_LIBRARY_PATH and DYLD_LIBRARY_PATH environment variables  "
	@echo "#                                                                      "
	@echo "======================================================================="
	@echo ""
else
msg:
	@echo ""
	@echo "======================================================================="
	@echo "#                                                                      "
	@echo "# Finished compiling                                                   "
	@echo "#                                                                      "
	@echo "# Before you attempt to run any executable in this code tree,          "
	@echo "# make sure you add                                                    "
	@echo "#                                                                      "
	@echo "#            $(TOP)/lib                                                "
	@echo "#                                                                      "
	@echo "# and                                                                  "
	@echo "#                                                                      "
	@echo "#            $(TOOLSDIR)/lib                                           "
	@echo "#                                                                      "
	@echo "# to your LD_LIBRARY_PATH environment variable                         "
	@echo "#                                                                      "
	@echo "======================================================================="
	@echo ""
endif
