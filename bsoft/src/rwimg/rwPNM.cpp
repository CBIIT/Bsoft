/**
@file	rwPNM.cpp
@brief	Functions for reading and writing PBM, PGM and PPM files
@author Bernard Heymann
@date	Created: 20110317
@date 	Modified: 20120211
**/

#include "rwPNM.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

char*		pnm_next_item(char* aptr)
{
	while ( !isspace(*aptr) ) aptr++;	// Skip any initial characters
	
	while ( isspace(*aptr) ) aptr++;	// Skip all white space
	
	while ( *aptr == '#' ) {			// Skip any comments
		while ( *aptr != '\n' ) aptr++;
		aptr++;
		while ( isspace(*aptr) ) aptr++;
	}
	
	return aptr;
}

unsigned char*	pnm_read_data(ifstream* fimg, Bimage* p)
{
	char			aline[MAXLINELEN], *aptr;
	int				n;
	long			i(0), len, v;
	long			imgsize = p->sizeX()*p->sizeY()*p->sizeZ()*p->channels();
	long			x, y, new_x = (p->sizeX()-1)/8 + 1;

	unsigned char*	data = p->data_alloc();
	unsigned short*	usdata = (unsigned short *) data;
	
	fimg->seekg(p->data_offset(), ios_base::cur);

    while ( fimg->getline(aline, MAXLINELEN) && i < imgsize ) {
		len = strlen(aline);
		for ( aptr = aline; *aptr != '\n' && aptr < aline + len; ) {
			if ( sscanf(aptr, "%ld%n", &v, &n) ) {
				if ( p->data_type() == Bit ) {
					if ( v ) {
						x = i%p->sizeX();
						y = i/p->sizeX();
						data[y*new_x+x/8] |= 0x80 >> x%8;
					}
				} else if ( p->data_type() == UCharacter ) data[i] = v;
				else usdata[i] = v;
				i++;
			}
			aptr += n;
		}
	}

	return data;
}

/**
@brief	Reading a PNM type file.
A simple image format.
	PBM			bit image format.
	PGM			grayscale image format.
	PPM			RGB color image format.
	File format extensions: .pbm, .pgm, .ppm
	The header contains ascii text indicating the image and element sizes.
	Data types:	byte, unsigned short
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@param 	img_select		image selection in multi-image file (-1 = all images).
@return	int 				0, <0 on error.
**/
int			readPNM(Bimage* p, int readdata, int img_select)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPNM: Filename = " << p->file_name() << endl;

	p->data_type(UCharacter);

	p->sizeZ(1);
	p->channels(1);
	
	string			ext(extension(p->file_name()));
	if ( ext.find("pbm") != string::npos ) p->data_type(Bit);
	else if ( ext.find("pgm") != string::npos ) p->compound_type(TSimple);
	else if ( ext.find("ppm") != string::npos ) {
		p->compound_type(TRGB);
		p->channels(3);
	} else {
		cerr << "Error: File with extension " << ext << " not recognized!" << endl;
		return -1;
	}
		
    ifstream*		fimg = new ifstream;
    fimg->open(p->file_name());
    if ( fimg->fail() ) return -1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPNM: File opened" << endl;

    char			aline[MAXLINELEN], *aptr;
	int				magic(0);
	for ( int i=0; i<MAXLINELEN; ++i ) aline[i] = 0;
	
	fimg->seekg(0, ios_base::end);
	int				len = fimg->tellg();
	fimg->seekg(0);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPNM: File length = " << len << endl;
	if ( len > MAXLINELEN ) len = MAXLINELEN;
	
	if ( !fimg->read(aline, len) ) return -2;
	sscanf(aline, "P%d", &magic);
	switch ( magic ) {
		case 1:
		case 4: p->data_type(Bit); p->channels(1); break;
		case 2:
		case 5: p->compound_type(TSimple); p->channels(1); break;
		case 3:
		case 6: p->compound_type(TRGB); p->channels(3); break;
		default: break;
	}
	
	aptr = pnm_next_item(aline);
	p->sizeX(atoi(aptr));
	aptr = pnm_next_item(aptr);
	p->sizeY(atoi(aptr));

	p->page_size(p->size());
	
	if ( p->data_type() > Bit ) {
		aptr = pnm_next_item(aptr);
		if ( atoi(aptr) > UCHAR_MAX ) p->data_type(UShort);
	}
	
	while ( !isspace(*aptr) ) aptr++;	// Skip any initial characters
	p->data_offset(aptr + 1 - aline);
		
//	if ( img_select > -1 ) {
//		if ( img_select >= p->images() ) img_select = p->images() - 1;
//		p->images(1);
//	} else img_select = -1;
	p->images(1);
	
	
	if ( !readdata ) {
		fimg->close();
		return 0;
	}
	
	if ( magic < 4 ) {
		pnm_read_data(fimg, p);
	} else {
		fimg->seekg(p->data_offset(), ios_base::beg);
		long	readsize = p->alloc_size();
		unsigned char*	data = p->data_alloc();
		fimg->read((char *)data, readsize);
	}
	
	fimg->close();
	
	delete fimg;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPNM: File read: type = " << magic << endl;

	return 0;
}

/**
@brief	Writing a PNM type file.
A simple image format.
@param	*p			the image structure.
@return	int 				0, <0 on error.
**/
int 	writePNM(Bimage* p)
{
	string			ext(extension(p->file_name()));
	
	if ( ext.find("pbm") != string::npos ) p->change_type(Bit);
	else if ( ext.find("pgm") != string::npos ) p->color_to_simple();
	
	if ( p->data_type() > UShort ) p->change_type(UShort);
	else if ( p->data_type() > UCharacter ) p->change_type(UCharacter);

    ofstream        fimg(p->file_name());
    if ( fimg.fail() ) return -1;
	
	long			n;
	long 			imgsize = p->sizeX()*p->sizeY()*p->channels()*p->data_type_size();
	if ( p->data_type() == Bit ) imgsize = (p->page_size()[0]/8)*p->sizeY();
	int				max = (int) p->maximum();
	if ( p->data_type() == UShort && max <= UCHAR_MAX ) max = UCHAR_MAX + 1;
	unsigned char*	data = p->data_pointer();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePNM: datatype=" << p->data_type() << endl;
	
	for ( n=0; n<p->images(); n++, data += imgsize ) {
		if ( p->data_type() == Bit ) fimg << "P4" << endl;
		else if ( p->compound_type() == TSimple ) fimg << "P5" << endl;
		else if ( p->compound_type() == TRGB ) fimg << "P6" << endl;
		else cerr << "Error: Compound type " << p->compound_type() << " not supported!" << endl;
		fimg << p->sizeX() << " " << p->sizeY() << endl;
		if ( p->data_type() > Bit ) {
			if ( p->data_type() == UCharacter || p->data_type() == UShort ) fimg << max << endl;
			else cerr << "Error: Data type " << p->data_type() << " not supported!" << endl;
		}
		fimg.write((char *)data, imgsize);
	}
	
	fimg.close();
	
	return 0;
}

