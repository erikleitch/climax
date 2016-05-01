PGUTIL=tPgUtil
$(PGUTIL): tPgUtil.o 
	$(FC) $(FC_FLAGS) -o $@ tPgUtil.o \
	$(LDPATH) $(LIBS) $(LIBDIR)/lib$(PROJ_PREFIX)PgUtil$(LIBSO_SUFFIX) \
	$(FC_PG_RPATH) $(PGLD) $(X11LD) -lXmu -lpng

BINEXC = $(PGUTIL)

test: generic_test $(BINEXC)

include $(TOP)/Makefile_generic
