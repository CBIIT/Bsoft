# Makefile for Bsoft
# Fri Jul 16 11:19:10 EDT 2021

# Bsoft directory: /Users/bernard/b20
# System: Darwin

SHELL = /bin/sh -x
MAKE = make

# Dependencies
LIBFFTW=/Users/bernard/b20/fftw-3.3.8
LIBTIFF=/Users/bernard/b20/tiff-4.3.0
LIBJPEG=/Users/bernard/b20/jpeg-9d
LIBPNG=/Users/bernard/b20/libpng-1.6.37

# Bsoft
BSOFT=/Users/bernard/b20/bsoft

TARFILE=bsoft2_1_2.tar
TGZFILE=bsoft2_1_2.tgz

all: fftw  tiff png jpeg bsoft

fftw:
	cd $(LIBFFTW) && $(MAKE) && $(MAKE) install

tiff:
	cd $(LIBTIFF) && $(MAKE) && $(MAKE) install

png:
	cd $(LIBPNG) && $(MAKE) && $(MAKE) install

jpeg:
	cd $(LIBJPEG) && $(MAKE) && $(MAKE) install

bsoft:
	cd $(BSOFT) && $(MAKE)

.PHONY : fftw  tiff png jpeg bsoft configure install clean tar

install:
	cd $(BSOFT) && $(MAKE) install

clean:
	cd $(LIBFFTW) && $(MAKE) clean
	cd $(LIBTIFF) && $(MAKE) clean
	cd $(LIBPNG) && $(MAKE) clean
	cd $(LIBJPEG) && $(MAKE) clean
	cd $(BSOFT) && $(MAKE) clean

cleanbsoft:
	cd $(BSOFT) && $(MAKE) clean

tar:
	tar cfv $(TARFILE) bconf bdist.plist bpkg breq.plist btar dev.rsc
	tar rfv $(TARFILE) bsoft/BVERSION bsoft/bsoft_conf bsoft/benv bsoft/bsetup
	tar rfv $(TARFILE) bsoft/include/*.h bsoft/src/*.cpp bsoft/src/*/*.cpp
	tar rfv $(TARFILE) bsoft/tcltk/*.h bsoft/tcltk/*.cpp
	tar rfv $(TARFILE) bsoft/Brun.app bsoft/Bshow.app bsoft/Scripts
	tar rfv $(TARFILE) bsoft/parameters bsoft/doc
	gzip -c $(TARFILE) > $(TGZFILE)

