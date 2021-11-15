/**
@file	rwSPE.cpp
@brief	Functions for reading and writing SPE CCD files
@author Bernard Heymann
@date	Created: 20081105
@date 	Modified: 20120409
**/

#include "rwSPE.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

extern char* month[];

/**
@brief	Reading a SPE image file format.
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@param 	img_select		image selection in multi-image file (-1 = all images).
@return	int					error code (<0 means failure).
A Princeton Instruments CCD image file format.
	Header size:				4100 bytes (fixed).
	File format extension:  	.SPE.
	Byte order determination:	data type field.
	Data types: 				0 = float, 1 = int,
								2 = short, 3 = unsigned short,
								4 = string/char, 5 = double,
								6 = char, 7 = unsigned char/byte.
**/
int			readSPE(Bimage* p, int readdata, int img_select)
{
 	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSPE: header size = " << sizeof(SPEhead) << endl;
	
	SPEhead*	header = new SPEhead;
	
	fimg->read((char *)header, SPESIZE);
	if ( fimg->fail() ) return -2;
	
    // Determine byte order and swap bytes if from little-endian machine
    int     	sb = 0;
    if ( header->datatype > 127 ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Swapping header byte order" << endl;
    	sb = 1;
		swapbytes((unsigned char *)&header->datatype, 2);
		swapbytes((unsigned char *)&header->xdim, 2);
		swapbytes((unsigned char *)&header->ydim, 2);
		swapbytes((unsigned char *)&header->NumFrames, 2);
		swapbytes((unsigned char *)&header->MaxIntensity, 4);
		swapbytes((unsigned char *)&header->MinIntensity, 4);
    }

	p->label(header->Comments);

	tm*			t = p->get_localtime();
	p->set_time(t);
 
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSPE: xdim=" << header->xdim << " ydim=" << header->ydim << " NumFrames=" << header->NumFrames << endl;
   
	// Map the parameters
	long			nimg(header->NumFrames);
	p->size(header->xdim, header->ydim, 1);
	p->channels(1);
	
	if ( img_select > -1 ) {
		if ( img_select >= nimg ) img_select = nimg - 1;
		nimg = 1;
	}
	
	p->images(nimg);

	switch ( header->datatype ) {
		case 0: p->data_type(Float); break;
		case 1: p->data_type(Integer); break;
		case 2: p->data_type(Short); break;
		case 3: p->data_type(UShort); break;
		case 4: p->data_type(SCharacter); break;
		case 5: p->data_type(Double); break;
		case 6: p->data_type(SCharacter); break;
		case 7: p->data_type(UCharacter); break;
		default: p->data_type(UShort); break;
	}

	p->maximum(header->MaxIntensity); 
	p->minimum(header->MinIntensity);	
	
	p->data_offset(SPESIZE);
	
	delete header;

	long			readsize = p->alloc_size();

	if ( readdata ) {
		p->data_alloc();
		fread_large(p->data_pointer(), readsize, p->data_offset(), fimg);
		if ( sb ) swapbytes(readsize, p->data_pointer(), p->data_type_size());
	}
	
	fimg->close();
	delete fimg;
	
	return 0;
}

/**
@brief	Writing a SPE image file format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
A Princeton Instruments CCD image file format.
**/
int			writeSPE(Bimage* p)
{
	p->color_to_simple();
	if ( p->data_type() != UCharacter && p->data_type() != UShort && 
		p->data_type() != Short && p->data_type() != Float ) p->change_type(Float);
	
	SPEhead*	header = new SPEhead;
	memset(header, 0, sizeof(SPEhead));

	strncpy(header->Comments, p->label().c_str(), 400);

	// Map the parameters
    tm* 		t = p->get_localtime();
	snprintf(header->date, 10, "%02d%3s%04d", t->tm_mday, month[t->tm_mon], t->tm_year);
	header->xdim = p->sizeX();
	header->ydim = p->sizeY();
	header->NumFrames = p->images();
	switch ( p->data_type() ) {
		case UCharacter: header->datatype = 7; break;
		case SCharacter: header->datatype = 6; break;
		case UShort: header->datatype = 3; break;
		case Short: header->datatype = 2; break;
		case Integer: header->datatype = 1; break;
		case Float: header->datatype = 0; break;
		case Double: header->datatype = 5; break;
		default: header->datatype = 0; break;
	}
	p->data_offset(SPESIZE);

	header->MaxIntensity = p->maximum(); 
	header->MinIntensity = p->minimum(); 
	
	long 			datasize = p->alloc_size();
	
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write((char *)header, p->data_offset());
	
	fimg.write((char *)p->data_pointer(), datasize);
	
	fimg.close();
	
	delete header;
		
	return 0;
}

