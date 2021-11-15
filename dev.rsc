
# Functions for Bsoft configuration
# Bernard Heymann
# 20030818 - 20210624

developer_resources()
{
	# Determine the system and developer resources
	SYS=`uname -s | cut -f1 -d"-"`
	OSVER=`uname -v`
	SDK=""
	if [[ $SYS =~ Darwin ]]; then
		SDK='/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk'
#		CSDK=`ls $SDK | head -1`
#		SDK=$SDK/$CSDK
	fi
	
#	echo "SDK = $SDK"

	# Determine the available compiler
	gv=`gcc 2>&1 | cut -f1 -d":"`
	if [[ $gv == "clang" ]]; then
		CC='clang'
		CXX='clang++'
	elif [[ $gv =~ llvm-gcc ]]; then
		CC='llvm-gcc'
		CXX='llvm-g++'
	elif [[ $gv =~ gcc ]]; then
		CC='gcc'
		CXX='g++'
	else
		echo Compilers not found!
		exit 1
	fi
}

check_dependencies()
{
	if [[ $SYS =~ Darwin ]]; then
		DEPTH='-depth 1 -type d'
	else
		DEPTH='-maxdepth 1 -type d'
	fi
	
	for DIR in "$BDIR"; do
		if [ ! -d "$LIBFFTW" ]; then
			LIBFFTW=`find $DIR $DEPTH -name 'fftw*'`
		fi

		if [ ! -d "$LIBTIFF" ]; then
			LIBTIFF=`find $DIR $DEPTH -name 'tiff*'`
		fi

		if [ ! -d "$LIBJPEG" ]; then
			LIBJPEG=`find $DIR $DEPTH -name 'jpeg*'`
		fi

		if [ ! -d "$LIBPNG" ]; then
			LIBPNG=`find $DIR $DEPTH -name 'libpng*'`
		fi

	done

	# Required
	HAVE_FFTW=0
	HAVE_XML=0

	if [ -e $LIBFFTW/api/fftw3.h ]; then
		echo "# FFTW3       found"
		HAVE_FFTW=1
	else
		echo "# No FFTW3 libraries found!"
		echo "# Looking for: $LIBFFTW/api/fftw3.h"
	fi

	if [ -e $SDK/usr/include/libxml2 ]; then
		echo "# libxml2     found"
		HAVE_XML=1
		x=`grep -A 3 xmlStrPrintf $SDK/usr/include/libxml2/libxml/xmlstring.h | grep msg | grep xmlChar`
		if [[ $x =~ xmlChar ]]; then
			HAVE_XML=2
		fi
	else
		echo "# No XML library found!"
		echo "# Looking for: $SDK/usr/include/libxml2"
	fi

	# Optional
	HAVE_TK=0
	HAVE_TIFF=0
	HAVE_PNG=0
	HAVE_JPEG=0
	
	# Tcl/Tk
	if [[ $SYS =~ Darwin ]]; then
#		if [ -e /Library/Frameworks/Tk.framework ]; then
#			TCL=/Library/Frameworks/Tcl.framework
#			TK=/Library/Frameworks/Tk.framework
#		else
			TCL=$SDK/System/Library/Frameworks/Tcl.framework
			TK=$SDK/System/Library/Frameworks/Tk.framework
#		fi
        TCLINC=$TCL/Headers
		TKINC=$TK/Headers
	elif [[ $SYS =~ Linux ]]; then
		TCL=/usr
		TK=/usr
		TCLINC=$TCL/include
		TKINC=$TK/include
		if [ -e /usr/include/tk ]; then
			TCLINC=/usr/include/tcl
			TKINC=/usr/include/tk
		fi
	fi

	if [ -e $TKINC/tk.h ]; then
		echo "# Tcl/Tk      found"
		HAVE_TK=1
	else
		echo "# No Tcl/Tk library found! Bshow will not be compiled!"
		echo "# Looking for: $TKINC/tk.h"
		TKO=""
	fi

	if [ -e $LIBTIFF/libtiff/tiff.h ]; then
		echo "# libtiff     found"
		HAVE_TIFF=1
	else
		echo "# No TIFF library found! TIFF files will not be supported!"
		echo "# Looking for: $LIBTIFF/libtiff/tiff.h"
	fi

	if [ -e $LIBPNG/png.h ]; then
		echo "# libpng      found"
		HAVE_PNG=1
	else
		echo "# No PNG library found! PNG files will not be supported!"
		echo "# Looking for: $LIBPNG/png.h"
	fi

	if [ -e $LIBJPEG/jpeglib.h ]; then
		echo "# libjpeg     found"
		HAVE_JPEG=1
	else
		echo "# No JPEG library found! JPEG files will not be supported!"
		echo "# Looking for: $LIBJPEG/jpeglib.h"
	fi


}

show_configuration()
{
	echo "--------------------------------------------------------------------------"
	today=`date +"%Y%m%d"`
	echo "Date:	$today"
	echo

	echo "System:		$OSVER"
	echo "C compiler:	$CC"
	echo "C++ compiler:	$CXX"
	echo "Developer resources: $SDK"
	echo

	echo "Bsoft source directory: $BSOFT"
	echo "Bsoft installation directory: $BINSTALL"
	echo "Dependencies:"
	echo "	FFTW3: $LIBFFTW"
	echo "	TIFF:  $LIBTIFF"
	echo "	JPEG:  $LIBJPEG"
	echo "	PNG:   $LIBPNG"
    echo "  Tcl:   $TCL"
    echo "  Tk:    $TK"
	echo "--------------------------------------------------------------------------"
	echo
}

detect_parallelization()
{
	echo "--------------------------------------------------------------------------"
	# Parallelization: GCD or OpenMP
	PAR=""
	if [[ $CC =~ clang ]]; then
		if [[ -e "$SDK/usr/include/dispatch" ]]; then
			echo "Parallelization: Grand Central Dispatch detected"
			PAR=GCD
		else
			echo "No parallelization using Grand Central Dispatch available"
			PAR=NONE
		fi
	else
		omp=`echo | cpp -fopenmp -dM | grep -i openmp | cut -f2 -d" "`
		if [ $omp = _OPENMP ]; then
			echo "Parallelization: OpenMP detected"
			PAR=OMP
		else
			echo "No parallelization using OpenMP available"
			PAR=NONE
		fi
	fi
	echo "--------------------------------------------------------------------------"
	echo
}
