/**
@file	rwBCR.cpp
@brief	Functions for reading and writing AFM BCR-STM files
@author Bernard Heymann
@date	Created: 20170214
@date	Modified: 20170214
**/

#include "rwBCR.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a BCR-STM image format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int			error code (<0 means failure).

A 2D image format intended for atomic force microscopy.
The format starts with a 2048 or 4096 byte text header,
composed of tag = value pairs:
	Tag				Value(s)		Remark
	fileformat		bcrstm			16 bit unsigned short
					bcrf			32 bit float
					bcrstm_unicode	16 bit unsigned short
					bcrf_unicode	32 bit float
	headersize		2048, 4096 (unicode)
	xpixels			x-size
	ypixels			y-size
	xlength			x-scanrange
	ylength			y-scanrange
	xunit, yunit, zunit		units for the three axes. If not defined nm will be the default unit.
	xlabel, ylabel, zlabel	labels for the three axes.
	current			tunneling current in nA (optional)
	bias			bias voltage in V (optional). 
	starttime		starting time of the scanning (MM DD YY hh:mm:ss:hh) (optional).
	scanspeed		measured in nm/sec (optional).
	intelmode		1: Little Endian; 0: Big Endian
	forcecurve		1: indicates that the data contain force curves with the approach curve followed by the retraction curve (optional)
	bit2nm			scale factor for scaling the integer height data to nm.
	xoffset, yoffset	physical offset in nm (optional). 
	voidpixels		number of void pixels, if the field is not present the number is set to zero. For the 16 bit integer bcrstm format void pixels should be set equal to 32767. For the 32 bit floating point bcrf format void pixels are set to 3.402823466E+38.

**/
int 	readBCR(Bimage* p, int readdata)
{
    ifstream		fimg;
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;

	p->images(1);
	p->channels(1);
	p->sizeZ(1);

	long			count(0), hsize(BCRSIZE);
    char			aline[MAXLINELEN];
	Bstring			s(" "), tag, val;

	while ( !fimg.eof() && count < hsize && s.length() ) {
		fimg.getline(aline, MAXLINELEN);
		s = aline;
		count += s.length();
//		cout << count << tab << s << endl;
		if ( s.contains("=") ) {
			tag = s.pre('=').no_space();
			val = s.post('=').no_space();
			if ( tag == "fileformat" ) {
				if ( val.contains("bcrstm") ) 
					p->data_type(UShort);
				else
					p->data_type(Float);
			} else if ( tag == "headersize" ) {
				hsize = val.integer();
			} else if ( tag == "xpixels" ) {
				p->sizeX(val.integer());
			} else if ( tag == "ypixels" ) {
				p->sizeY(val.integer());
			}
		}
	}
	s = 0;
	fimg.clear();
	
//	cout << "Count: " << count << endl;
	
	p->data_offset(hsize);

	if ( readdata) {
		fimg.seekg(p->data_offset(), ios_base::beg);
//		cout << "offset: " << p->data_offset() << endl;
		p->data_alloc();
//		cout << "alloc_size: " << p->alloc_size() << endl;
		fimg.read((char *)p->data_pointer(), p->alloc_size());
    	if ( fimg.fail() ) {
			cerr << "Error: No data read!" << endl;
			cerr << fimg.rdstate() << endl;
			if (fimg.bad()) {
        		cerr << "I/O error while reading\n";
    		} else if (fimg.eof()) {
        		cerr << "End of file reached successfully\n";
    		} else if (fimg.fail()) {
        		cerr << "Non-integer data encountered\n";
			}
			return -1;
		}
	}
	
	fimg.close();
	
	return 0;
}

/**
@brief	Writing a BCR-STM image format.
@param	*p			the image structure.
@return	int			error code (<0 means failure).
**/
int 	writeBCR(Bimage* p)
{
    ofstream        fimg;
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;
	
	Bstring			header("fileformat = ");
	
	if ( p->data_type() == Float ) header += "bcrf\n";
	else header += "bcrstm\n";
	header += "headersize = " + Bstring(BCRSIZE, "%d\n");
	header += "xpixels = " + Bstring(p->sizeX(), "%d\n");
	header += "ypixels = " + Bstring(p->sizeY(), "%d\n");
	
//	cout << "header length: " << header.length() << endl;
//	cout << header << "---" << endl;
//	cout << BCRSIZE - header.length() << endl;
	
	header += Bstring(' ', BCRSIZE - header.length());
//	cout << "header length: " << header.length() << endl;

	fimg.write((char *)header.c_str(), BCRSIZE);
	
	header = 0;

	fimg.write((char *)p->data_pointer(), p->alloc_size());
		
	fimg.close();
	
	return 0;
}

