/**
@file	rwBrookhavenSTEM.cpp
@brief	Functions for reading and writing Brookhaven STEM files
@author Bernard Heymann
@date	Created: 20050729
@date	Modified: 20130307
**/

#include "rwBrookhavenSTEM.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

#define	LUTSIZE	256

int			img_stem_lut_apply(Bimage* p, int neg_stain_flag)
{
	if ( !p->data_pointer() || p->images() < 2 ) return -1;
	
	long			i;
	long			imgsize = p->sizeX()*p->sizeY();
	unsigned char*	data = p->data_pointer();
	int				j;
	int				prefac(146);
	int				la_gain(10);
	int				sa_gain(10);
	int				bf_gain(0);
	
	if ( neg_stain_flag ) {
		la_gain = 2;
		sa_gain = 0;
		bf_gain = 1;
	}
	
	double			ila(0), isa(0), ibf(0), it;
	float			la_lut[LUTSIZE], sa_lut[LUTSIZE], bf_lut[LUTSIZE];
	
	// Data originating from the microscope
	for ( i=0; i<LUTSIZE; i++ )
		la_lut[i] = sa_lut[i] = bf_lut[i] = i;
			
	for ( j=0; j<7000 && ( ila < LUTSIZE-1 || isa < LUTSIZE-1 ); j++ ) {
		ila = (1 - exp(-j/2500.0)) * exp(-j/10000.0);
		isa = (1 - exp(-j/1964.0)) * exp(-j/2000.0);
		ibf = exp(-j/1100.0);
		it = ila + isa + ibf;
//		cout << ila << " " << isa << " " << ibf << " " << it << endl;
		ila = prefac*la_gain*ila/it;
		isa = prefac*sa_gain*isa/it;
		ibf = prefac*bf_gain*ibf/it;
//		cout << j << " " << ila << " " << isa << " " << ibf << " " << it << endl;
		if ( ila < LUTSIZE ) la_lut[(int)ila] = prefac*la_gain*j/2500.0;
		if ( isa < LUTSIZE ) sa_lut[(int)isa] = prefac*sa_gain*j/1964.0;
		if ( ibf < LUTSIZE ) bf_lut[(int)ibf] = prefac*bf_gain*(1 - j/1100.0);
	}

	for ( i=0; i<LUTSIZE; i++ ) {
//		cout << i << " " << la_lut[i] - i << " " << sa_lut[i] - i << " " << bf_lut[i] - i << endl;
		if ( la_lut[i] >= LUTSIZE ) la_lut[i] = LUTSIZE - 1;
		if ( sa_lut[i] >= LUTSIZE ) sa_lut[i] = LUTSIZE - 1;
		if ( bf_lut[i] >= LUTSIZE ) bf_lut[i] = LUTSIZE - 1;
		if ( bf_lut[i] < 0 ) bf_lut[i] = 0;
	}

	for ( i=0; i<imgsize; i++ ) data[i] = (int)la_lut[data[i]];
		
	if ( neg_stain_flag ) {
		for ( i=0, j=imgsize; i<imgsize; i++, j++ ) {
			data[j] = (int)bf_lut[data[j]];
			it = data[j] - data[i]/la_gain;
			if ( it < 0 ) it = 0;
			data[i] = (int) it;
		}
	} else
		for ( i=0, j=imgsize; i<imgsize; i++, j++ ) data[j] = (int)sa_lut[data[j]];
	
	return 0;
}

/**
@brief	Reading a BrookhavenSTEM image file format.
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@return	int					error code (<0 means failure).
A 2D image format used for the Brookhaven STEM.
	It has a text header of 4096 bytes (fixed size).
	The default size is 512x512.
	File format extension:  	.dat.
	Data type: 					unsigned char.
	The header specifies the following keywords:
		SCAN:		pixel size in internal units
	The data is packed into two interleaved channels which are unpacked
	into two separate images in this function. The first image contains
	the signal from the LA (large area) detector and the second image
	contains the signal from either the SA (small area) detector, or
	the BF (bright field) detector. The data is corrected for signal
	non-linearity on the fly when read in.
**/
int 	readBrookhavenSTEM(Bimage* p, int readdata)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

	long	i, j;
	char*			header = new char[BrookhavenSTEMSIZE];
	for ( i=0; i<BrookhavenSTEMSIZE; i++ ) header[i] = 0;
	
	fimg->read(header, BrookhavenSTEMSIZE-1);
	if ( fimg->fail() ) return -2;
	
	// Extract the image attributes from the file header
//	p->label(header);

	double			pixel_size(1);
	char*			aptr = strstr(header, "SCAN");
	if ( aptr ) aptr += 5;
	else aptr = header + 72;	// Old format
	sscanf(aptr, "%lf", &pixel_size);
	pixel_size /= 0.0512;	// Convert from internal units to angstrom
	if ( pixel_size < 1e-3 ) pixel_size = 1;
	
	int				neg_stain_flag(0);
	char			subhead[128];
	strncpy(subhead, header+16, 40);
	for ( i=0; i<40; i++ ) subhead[i] = tolower(subhead[i]);
	if ( strstr(subhead, "van") ) neg_stain_flag = 1;
	if ( strstr(subhead, "ur") ) neg_stain_flag = 1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readBrookhavenSTEM: image header parsed (" << neg_stain_flag << ")" << endl;

	p->data_offset(BrookhavenSTEMSIZE);

	fimg->seekg(0, ios::end);
	long  			datasize = (long)fimg->tellg() - p->data_offset();
	
    p->size(512, 512, 1);
    p->channels(1);
	
	p->data_type(UCharacter);
	if ( datasize > 300000 ) p->images(2);
	else p->images(1);

	p->sampling(pixel_size, pixel_size, 1);

	delete[] header;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readBrookhavenSTEM: file_size=" << datasize << " alloc_size=" << p->alloc_size() << endl;

	if ( readdata ) {
		p->data_alloc();
		fread_large(p->data_pointer(), p->alloc_size(), p->data_offset(), fimg);
	}
	
	fimg->close();
	delete fimg;

	// Unpack the interleaved images
	if ( p->images() > 1 && p->data_pointer() ) {
		long	imgsize = p->sizeX()*p->sizeY();
		unsigned char*	data = p->data_pointer();
		unsigned char*	temp = new unsigned char[imgsize];
		for ( i=j=0; i<imgsize; i++, j+=2 ) {
			temp[i] = data[j+1];
			data[i] = data[j];
		}
		for ( i=0, j=imgsize; i<imgsize; i++, j++ )
			data[j] = temp[i];
		delete[] temp;
		if ( readdata > 1 ) img_stem_lut_apply(p, neg_stain_flag);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readBrookhavenSTEM: image read" << endl << endl;
	
	return 0;
}

/**
@brief	Writing a BrookhavenSTEM image file format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
A 2D image format used for the Brookhaven STEM.
**/
int 	writeBrookhavenSTEM(Bimage* p)
{
	p->change_type(UCharacter);
	
	char*	header = new char[BrookhavenSTEMSIZE];
	memset(header, ' ', BrookhavenSTEMSIZE);
	
	strncpy(header, p->label().c_str(), BrookhavenSTEMSIZE);
	
	snprintf(header + 70, BrookhavenSTEMSIZE-70, "SCAN:%6.3f", p->sampling(0)[0]*0.0512);
	
	p->data_offset(BrookhavenSTEMSIZE);
	
	long			i, j, k;
	long   			imgsize = p->sizeX()*p->sizeY()*p->sizeZ();
	long   			datasize = p->data_size();
	unsigned char*	data = p->data_pointer();
	unsigned char*	temp = data;

	// Interleaved data
	if ( p->images() > 1 ) {
		temp = new unsigned char[datasize];
		for ( i=0, j=imgsize, k=0; i<imgsize; i++, j++, k+=2 ) {
			temp[k] = data[i];
			temp[k+1] = data[j];
		}
	}
	
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write((char *)header, p->data_offset());
	if ( p->data_pointer() ) fimg.write((char *)temp, datasize);
	
	fimg.close();

	delete[] header;
	if ( temp != data ) delete[] temp;
		
	return 0;
}

