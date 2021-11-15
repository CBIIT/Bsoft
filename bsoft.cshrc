# Environmental variables for Bsoft on Darwin (x86_64)
 
setenv BSOFT .
setenv BPARAM $BSOFT/parameters/
 
if ( $?PATH ) then
	setenv PATH $BSOFT/bin:$PATH
else
	setenv PATH $BSOFT/bin
endif
 
if ( $?LD_LIBRARY_PATH ) then
	setenv LD_LIBRARY_PATH $BSOFT/lib:$LD_LIBRARY_PATH
else
	setenv LD_LIBRARY_PATH $BSOFT/lib
endif
 
if ( $?DYLD_LIBRARY_PATH ) then
	setenv DYLD_LIBRARY_PATH $BSOFT/lib:$DYLD_LIBRARY_PATH
else
	setenv DYLD_LIBRARY_PATH $BSOFT/lib
endif
 
