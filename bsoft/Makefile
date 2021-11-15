# Makefile for Bsoft
# Fri Jul 16 11:19:10 EDT 2021

# Bsoft source directory:
BDIR=/Users/bernard/b20/bsoft

# Bsoft install directory:
BINSTALL=../../b20install

# System: Darwin

SHELL = /bin/sh -x
MAKE = make

# Compiler settings
CC=clang
CXX=clang++

CFLAGS=-O3 -fPIC -std=c++14 -Wall -Wno-sign-compare -arch x86_64 -DHAVE_XML -DHAVE_GCD -DHAVE_TIFF -DHAVE_PNG -DHAVE_JPEG 
INCLUDES=-I./ -I/Users/bernard/b20/bsoft/include -I/Users/bernard/b20/bsoft/radon -I/Users/bernard/b20/bsoft/eer -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/libxml2 -I/Users/bernard/b20/fftw-3.3.8/include -I/Users/bernard/b20/tiff-4.3.0/libtiff -I/Users/bernard/b20/libpng-1.6.37 -I/Users/bernard/b20/jpeg-9d -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Tcl.framework/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Tk.framework/Headers

# Library loader settings
LDFLAGS=-all_load -dynamiclib -arch x86_64 -headerpad_max_install_names -install_name @loader_path/../lib/libbsoft.dylib
TKLDFLAGS=-all_load -dynamiclib -arch x86_64 -framework Tcl -framework Tk -lbsoft -headerpad_max_install_names -install_name @loader_path/../lib/libbsoft.dylib/libbshow.dylib
LINKOPTS=-bind_at_load -fPIC -arch x86_64 -L/Users/bernard/b20/bsoft/lib -lbsoft
LINKLIBS= -lxml2 -llzma  -lpthread -lz
LIBLIST= /Users/bernard/b20/fftw-3.3.8/lib/libfftw3f.a /Users/bernard/b20/fftw-3.3.8/lib/libfftw3f_threads.a /Users/bernard/b20/tiff-4.3.0/lib/libtiff.a /Users/bernard/b20/libpng-1.6.37/lib/libpng.a /Users/bernard/b20/jpeg-9d/lib/libjpeg.a
LIBBSOFT=lib/libbsoft.dylib
LIBTCLTK=/Users/bernard/b20/bsoft/lib/libbshow.dylib
LIBDIR=
TKO=/Users/bernard/b20/bsoft/tcltk/bsoft_tcl.o

# Sources, objects and executables
EXE_SOURCES=$(wildcard src/*.cpp)
LIB_SOURCES=$(wildcard src/*/*.cpp) $(wildcard radon/*.cpp) $(wildcard eer/*.cpp)
TK_SOURCES=$(wildcard tcltk/*.cpp)
EXE_OBJECTS=$(patsubst %.cpp,%.o,$(EXE_SOURCES))
LIB_OBJECTS=$(patsubst %.cpp,%.o,$(LIB_SOURCES))
TK_OBJECTS=$(patsubst %.cpp,%.o,$(TK_SOURCES))
EXECUTABLES=$(patsubst src/%.cpp,bin/%,$(EXE_SOURCES))

TARFILE=bsoft.tar
TGZFILE=bsoft.tgz

all: $(LIB_OBJECTS) $(TK_OBJECTS) $(RDN_OBJECTS) $(LIBBSOFT) $(LIBTCLTK) $(EXE_OBJECTS) $(EXECUTABLES)

.PHONY : clean tar

%.o: %.cpp
	$(CXX) $(INCLUDES) $(CFLAGS) $< -c -o $@

$(LIBBSOFT):
	$(CXX) -o $(LIBBSOFT) $(LDFLAGS) $(TKO) $(LIB_OBJECTS) $(LIBDIR) $(LINKLIBS) $(LIBLIST)

$(EXECUTABLES):
	$(CXX) -o $@ $(patsubst bin/%,src/%.o,$@) $(LINKOPTS) $(LIBDIR) $(LINKLIBS)

$(LIBTCLTK):
	$(CXX) $(LINKOPTS) $(TKLDFLAGS) -o $@ $(TK_OBJECTS) $(LIBDIR) $(LINKLIBS)

install:
	mkdir -p $(BINSTALL)
	cp BVERSION $(BINSTALL)
	cp bsetup $(BINSTALL)
	cp benv $(BINSTALL)
	cp -r $(BDIR)/lib $(BINSTALL)
	cp -r $(BDIR)/bin $(BINSTALL)
	cp -r $(BDIR)/tcltk $(BINSTALL)
	cp -r $(BDIR)/parameters $(BINSTALL)
	cp -r $(BDIR)/Scripts $(BINSTALL)
	cd $(BINSTALL) && ./benv $(BINSTALL)
	cd $(BINSTALL)/bin && ln -vfs ../tcltk/bshow bshow && ln -vfs ../tcltk/brun brun
	cp macinstall $(BINSTALL)
	cd $(BINSTALL)/bin && ln -vfs ../tcltk/bshowX bshowX && ln -vfs ../tcltk/brunX brunX
	cp -r $(BDIR)/Bshow.app $(BINSTALL)
	cp -r $(BDIR)/Brun.app $(BINSTALL)
	cp -r $(BDIR)/doc $(BINSTALL)

clean:
	rm */*.o src/*/*.o bin/* lib/*
	rm -f */blevel* */bmaskmod* src/rwparam/rwmodel_param.cpp

tar:
	tar cfv $(TARFILE) include/*.h tcltk/*.h $(LIB_SOURCES) $(TK_SOURCES) $(EXE_SOURCES)
	gzip -c $(TARFILE) > $(TGZFILE)

