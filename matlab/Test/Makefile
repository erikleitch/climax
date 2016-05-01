ENG=tWriteMiriadData
$(ENG): tWriteMiriadData.o 
	$(CC) $(CCINC) -o $@ $(CCFLAGS) $(LDPATH) $(RPATH) $(LIBS) \
	-Xlinker -R -Xlinker /opt/Matlab2009a/R2009a_linux32/unix/bin/glnx86 \
	-L/opt/Matlab2009a/R2009a_linux32/unix/bin/glnx86 -lmat -leng \
	/opt/Matlab2009a/R2009a_linux32/unix/sys/os/glnx86/libstdc++.so.6 \
	$<

#	/opt/Matlab2009a/R2009a_linux32/unix/bin/glnx86/libustdio.so.38 \

RUN=tRunMatlab
$(RUN): tRunMatlab.o 
	$(CC) $(CCINC) -o $@ $(CCFLAGS) $(LDPATH) $(RPATH) $(LIBS) \
	-Xlinker -R -Xlinker /opt/Matlab2009a/R2009a_linux32/unix/bin/glnx86 \
	-L/opt/Matlab2009a/R2009a_linux32/unix/bin/glnx86 -lmat -leng \
	/opt/Matlab2009a/R2009a_linux32/unix/sys/os/glnx86/libstdc++.so.6 \
	$<

#	/opt/Matlab2009a/R2009a_linux32/unix/bin/glnx86/libustdio.so.38 \

BINEXC = $(RUN) $(ENG)

test: generic_test

include $(TOP)/Makefile_generic

