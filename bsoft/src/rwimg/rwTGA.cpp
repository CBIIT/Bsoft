/**
@file	rwTGA.cpp
@brief	Functions for reading and writing Truevision TGA files
@author Bernard Heymann
@date	Created: 20150811
@date 	Modified: 20150815
**/

#include "rwTGA.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Given a pointer to packed TGA color data, and the number of bits used
// per color return an unpacked 32-bit color.
RGB<unsigned char> 	TGA_to_RGB(const unsigned char * pColorData, unsigned int bitsPerColor )
{
    RGB<unsigned char> 	color;
    switch( bitsPerColor )
    {
    case 15:
        // MSB           LSB
        // xRRRRRGG GGGBBBBB 
        color[0] = (pColorData[1] >> 2) & 0x1F;
        color[1] = ((pColorData[1] << 3) & 0x1C) | ((pColorData[0] >> 5) & 0x07);
        color[2] = (pColorData[0] & 0x1F);
  
        // Convert 5-bit channels to 8-bit channels by left shifting by three
        // and repeating the first three bits to cover the range [0,255] with
        // even spacing. For example:
        //   channel input bits:        4 3 2 1 0
        //   channel output bits: 4 3 2 1 0 4 3 2
        color[0] = (color[0] << 3) | (color[0] >> 2);
        color[1] = (color[1] << 3) | (color[1] >> 2);
        color[2] = (color[2] << 3) | (color[2] >> 2);
        break;
  
    case 24:
        // MSB                    LSB
        // RRRRRRRR GGGGGGGG BBBBBBBB 
        color[0] = pColorData[2];
        color[1] = pColorData[1];
        color[2] = pColorData[0];
        break;
	}
	
    return color;
}

// Given a pointer to packed TGA color data, and the number of bits used
// per color return an unpacked 32-bit color.
RGBA<unsigned char> TGA_to_RGBA(const unsigned char * pColorData, unsigned int bitsPerColor )
{
    RGBA<unsigned char> 	color;
    switch( bitsPerColor )
    {
    case 15:
        // MSB           LSB
        // xRRRRRGG GGGBBBBB 
        color[3] = 255;
        color[0] = (pColorData[1] >> 2) & 0x1F;
        color[1] = ((pColorData[1] << 3) & 0x1C) | ((pColorData[0] >> 5) & 0x07);
        color[2] = (pColorData[0] & 0x1F);
  
        // Convert 5-bit channels to 8-bit channels by left shifting by three
        // and repeating the first three bits to cover the range [0,255] with
        // even spacing. For example:
        //   channel input bits:        4 3 2 1 0
        //   channel output bits: 4 3 2 1 0 4 3 2
        color[0] = (color[0] << 3) | (color[0] >> 2);
        color[1] = (color[1] << 3) | (color[1] >> 2);
        color[2] = (color[2] << 3) | (color[2] >> 2);
        break;
  
    case 16:
        // MSB           LSB
        // ARRRRRGG GGGBBBBB
        color[3] = 255 * ((pColorData[1] & 0x80) >> 7);
        color[0] = (pColorData[1] >> 2) & 0x1F;
        color[1] = ((pColorData[1] << 3) & 0x1C) | ((pColorData[0] >> 5) & 0x07);
        color[2] = (pColorData[0] & 0x1F);
  
        // Convert 5-bit channels to 8-bit channels by left shifting by three
        // and repeating the first three bits to cover the range [0,255] with
        // even spacing. For example:
        //   channel input bits:        4 3 2 1 0
        //   channel output bits: 4 3 2 1 0 4 3 2
        color[0] = (color[0] << 3) | (color[0] >> 2);
        color[1] = (color[1] << 3) | (color[1] >> 2);
        color[2] = (color[2] << 3) | (color[2] >> 2);
        break;
  
    case 24:
        // MSB                    LSB
        // RRRRRRRR GGGGGGGG BBBBBBBB 
        color[3] = 255;
        color[0] = pColorData[2];
        color[1] = pColorData[1];
        color[2] = pColorData[0];
        break;
  
    case 32:
        // MSB                             LSB
        // AAAAAAAA RRRRRRRR GGGGGGGG BBBBBBBB 
        color[3] = pColorData[3];
        color[0] = pColorData[2];
        color[1] = pColorData[1];
        color[2] = pColorData[0];
        break;
    }
	
	color[3] = 255 - color[3];
	
//	cout << color << endl;
  
    return color;
}

/**
@brief	Reading a Truevision TGA image file format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int			error code (<0 means failure).
A color image format.
	Header size:				18 bytes (fixed).
	File format extensions:  	.tga
	The identifier is ???.
	Byte order:					Little endian.
	Data types: 				1 = byte, 2 = float, 3 = complex float
								4 = 3-value vector, 5 = view
**/
int 	readTGA(Bimage *p, int readdata)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readTGA: " << p->file_name() << endl;
	
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

	TGAhead*	header = new TGAhead;
	TGAfoot*	footer = new TGAfoot;
	
//	cout << "header size: " << sizeof(TGAhead) << endl;
//	cout << "footer size: " << sizeof(TGAfoot) << endl;
	
	fimg->seekg(-25, fimg->end);
	fimg->read((char *)footer, 25);
	if ( fimg->fail() ) return -2;
	
	if ( strncmp(footer->sig, "TRUEVISION-XFILE", 16) == 0 ) {
		cout << "Footer found" << endl;
		cout << "Extension offset: " << footer->extoff << endl;
		cout << "Developer offset: " << footer->devoff << endl;
	}

	fimg->seekg(0, fimg->beg);
	fimg->read((char *)header, TGASIZE);
	if ( fimg->fail() ) return -2;

	char		id[256];
	fimg->read(id, header->idlen);
	if ( fimg->fail() ) return -3;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readTGA: " << id << endl;
		cout << id << endl;
		cout << "Image type: " << (int)header->imgtype << endl;
		cout << "Width: " << header->width << " height: " << header->height << endl;
		cout << "Pixel depth: " << (int)header->pixeldepth << endl;
		cout << "Image descriptor: " << (int)header->imgdesc << endl;
		cout << "Color map: " << header->colmapindex << tab << header->colmaplength
			<< tab << (int)header->colmapbits << endl;
	}
	
	p->label(id);
	
    // Determine byte order and swap bytes if from little-endian machine
    unsigned char*	b = (unsigned char *) header;
    long			i, sb(0);
    if ( header->imgtype > 11 ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Swapping header byte order for 2-byte types" << endl;
    	sb = 1;
    	for ( i=3; i<6; i+=2 ) swapbytes(b+i, 2);
    	for ( i=8; i<17; i+=2 ) swapbytes(b+i, 2);
		swapbytes(b+3, 2);
    }
	
	
	int				entrylen((header->colmapbits-1)/8+1);
	if ( header->colmaplength ) {
		char			colmap[entrylen*header->colmaplength];
		fimg->read(colmap, entrylen*header->colmaplength);
		if ( fimg->fail() ) return -4;
	}

	// Map the parameters
	p->size(header->width, header->height, 1);
	p->images(1);
	p->channels(1);
	if ( header->imgtype == 3 ) {
		switch ( header->pixeldepth ) {
			case 8: p->data_type(UCharacter); break;
			case 16: p->data_type(UShort); break;
			case 32: p->data_type(UInteger);
			default: p->data_type(UCharacter); break;
		}
	} else {
		p->data_type(UCharacter);
		if ( header->pixeldepth == 15 || header->pixeldepth == 24 ) {
			p->channels(3);
			p->compound_type(TRGB);
		} else {
			p->channels(4);
			p->compound_type(TRGBA);
		}
	}
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readTGA: mode=" << header->imgtype << " datatype=" << p->data_type() << endl;
	p->data_offset(TGASIZE + header->idlen + entrylen*header->colmaplength);
	
	// Allocating the single sub-image and setting its parameters

	tm*			t = p->get_localtime();
	p->set_time(t);
	
	unsigned char*	nudata;
	long			j, bytes((header->pixeldepth-1)/8+1);
	long			ds(header->width*header->height*bytes);
		
	if ( readdata ) {
		nudata = new unsigned char[ds];
		if ( fread_large(nudata, ds, p->data_offset(), fimg) < 0 ) return -1;
		if ( sb && bytes > 1 ) swapbytes(ds, nudata, bytes);
		if ( p->channels() < 3 ) {
			p->data_assign((unsigned char *) nudata);
		} else {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG readTGA: converting color data" << endl;
			p->data_alloc();
			for ( i=j=0; i<header->width*header->height; i++, j+=bytes )
				if ( p->channels() == 3 ) p->set(i, TGA_to_RGB(nudata+j, header->pixeldepth));
				else p->set(i, TGA_to_RGBA(nudata+j, header->pixeldepth));
			delete[] nudata;
		}
		if ( header->imgdesc & 16 ) p->reslice("-xyz");
		if ( header->imgdesc & 32 ) p->reslice("x-yz");
	}
		
	delete header;
	delete footer;
	
	fimg->close();
	delete fimg;
		
	return 0;
}

/**
@brief	Writing a Truevision TGA map image file format.
@param	*p			the image structure.
@return	int			error code (<0 means failure).
**/
int 	writeTGA(Bimage *p)
{
	
    switch ( p->data_type() ) {
		case Bit: case UCharacter: case SCharacter:
			p->change_type(UCharacter); break;
		case UShort: case Short:
			p->change_type(UShort); break;
		case UInteger: case Integer:
			p->change_type(UInteger); break;
		case Float: case Double:
			p->change_type(UInteger); break;
    	default:
			p->change_type(UCharacter); break;
    }
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeTGA: header size=" << sizeof(TGAhead) << endl;
	
	TGAhead*	header = new TGAhead;
	memset(header, 0, sizeof(TGAhead));
	
	// Map the parameters
	if ( p->label().length() > 255 ) header->idlen = (unsigned char)255;
	else header->idlen = p->label().length();
	
	if ( p->compound_type() == TSimple ) {
		header->imgtype = 3;
		header->pixeldepth = p->data_type_size()*8;
	} else {
		header->imgtype = 2;
		header->pixeldepth = p->channels()*8;
	}
	
	header->width = p->size()[0];
	header->height = p->size()[1];
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeTGA: pixeldepth=" << header->pixeldepth << endl;
	
	long			i;
	
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write((char *)header, TGASIZE);
	fimg.write((char *)p->label().c_str(), header->idlen);
	if ( p->data_pointer() ) {
		char*			data = new char[p->alloc_size()];
		unsigned char*	ptr;
		if ( p->compound_type() < TRGB ) {
			for ( i=0, ptr=p->data_pointer(); i<p->alloc_size(); i++, ptr++ )
				data[i] = *ptr;
		} else {
			for ( i=0, ptr=p->data_pointer(); i<p->alloc_size(); ptr+=p->channels() ) {
				data[i++] = ptr[2];
				data[i++] = ptr[1];
				data[i++] = ptr[0];
				if ( p->channels() > 3 ) data[i++] = 255 - ptr[3];
			}
		}
		fimg.write(data, p->alloc_size());
		delete[] data;
	}
	
	fimg.close();
	
	delete header;
		
	return 0;
}

