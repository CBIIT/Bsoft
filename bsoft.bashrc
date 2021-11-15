# Environmental variables for Bsoft on Darwin (x86_64)
 
BSOFT=.
export BSOFT
BPARAM=$BSOFT/parameters/
export BPARAM
 
if [ "$PATH" ]; then
	PATH=$BSOFT/bin:$PATH
else
	PATH=$BSOFT/bin
fi
export PATH
 
if [ "$LD_LIBRARY_PATH" ]; then
	LD_LIBRARY_PATH=$BSOFT/lib:$LD_LIBRARY_PATH
else
	LD_LIBRARY_PATH=$BSOFT/lib
fi
export LD_LIBRARY_PATH
 
if [ ${DYLD_LIBRARY_PATH} ]; then
	DYLD_LIBRARY_PATH=$BSOFT/lib:$DYLD_LIBRARY_PATH
else
	DYLD_LIBRARY_PATH=$BSOFT/lib
fi
export DYLD_LIBRARY_PATH
 
