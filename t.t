clang++ -I./ -I/Users/bernard/b20/bsoft/include -I/Users/bernard/b20/bsoft/radon -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/libxml2 -I/Users/bernard/b20/fftw-3.3.8/include -I/Users/bernard/b20/tiff-4.1.0/libtiff -I/Users/bernard/b20/libpng-1.6.37 -I/Users/bernard/b20/jpeg-9d -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Tcl.framework/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Tk.framework/Headers -O3 -fPIC -std=c++11 -Wall -Wno-sign-compare -arch x86_64 -DHAVE_XML -DHAVE_GCD -DHAVE_TIFF -DHAVE_PNG -DHAVE_JPEG  src/model/model_symmetry.cpp -c -o src/model/model_symmetry.o
+ clang++ -I./ -I/Users/bernard/b20/bsoft/include -I/Users/bernard/b20/bsoft/radon -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/libxml2 -I/Users/bernard/b20/fftw-3.3.8/include -I/Users/bernard/b20/tiff-4.1.0/libtiff -I/Users/bernard/b20/libpng-1.6.37 -I/Users/bernard/b20/jpeg-9d -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Tcl.framework/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Tk.framework/Headers -O3 -fPIC -std=c++11 -Wall -Wno-sign-compare -arch x86_64 -DHAVE_XML -DHAVE_GCD -DHAVE_TIFF -DHAVE_PNG -DHAVE_JPEG src/model/model_symmetry.cpp -c -o src/model/model_symmetry.o
clang++ -I./ -I/Users/bernard/b20/bsoft/include -I/Users/bernard/b20/bsoft/radon -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/libxml2 -I/Users/bernard/b20/fftw-3.3.8/include -I/Users/bernard/b20/tiff-4.1.0/libtiff -I/Users/bernard/b20/libpng-1.6.37 -I/Users/bernard/b20/jpeg-9d -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Tcl.framework/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Tk.framework/Headers -O3 -fPIC -std=c++11 -Wall -Wno-sign-compare -arch x86_64 -DHAVE_XML -DHAVE_GCD -DHAVE_TIFF -DHAVE_PNG -DHAVE_JPEG  src/model/model_transform.cpp -c -o src/model/model_transform.o
+ clang++ -I./ -I/Users/bernard/b20/bsoft/include -I/Users/bernard/b20/bsoft/radon -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/libxml2 -I/Users/bernard/b20/fftw-3.3.8/include -I/Users/bernard/b20/tiff-4.1.0/libtiff -I/Users/bernard/b20/libpng-1.6.37 -I/Users/bernard/b20/jpeg-9d -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Tcl.framework/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Tk.framework/Headers -O3 -fPIC -std=c++11 -Wall -Wno-sign-compare -arch x86_64 -DHAVE_XML -DHAVE_GCD -DHAVE_TIFF -DHAVE_PNG -DHAVE_JPEG src/model/model_transform.cpp -c -o src/model/model_transform.o
src/model/model_transform.cpp:102:63: error: expected ')'
                        comp->location((comp->location() - origin) * scale + origin;
                                                                                   ^
src/model/model_transform.cpp:102:18: note: to match this '('
                        comp->location((comp->location() - origin) * scale + origin;
                                      ^
src/model/model_transform.cpp:206:44: error: expected ')'
                comp2->location(comp->location() - origin;;
                                                         ^
src/model/model_transform.cpp:206:18: note: to match this '('
                comp2->location(comp->location() - origin;;
                               ^
src/model/model_transform.cpp:428:9: error: 'loc' is a private member of 'Bcomponent'
                comp->loc *= binf;
                      ^
/Users/bernard/b20/bsoft/include/rwmodel.h:157:18: note: declared private here
        Vector3<float>          loc;                    // Location coordinates (angstroms)
                                ^
3 errors generated.
make: *** [src/model/model_transform.o] Error 1
