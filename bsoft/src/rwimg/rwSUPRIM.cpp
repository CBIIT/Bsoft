/**
@file	rwSUPRIM.cpp
@brief	Functions for reading and writing SUPRIM files
@author Bernard Heymann
@date	Created: 19990930
@date 	Modified: 20120409
**/

#include "rwSUPRIM.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

View		view_from_suprim_euler(SUPRIMhead* header)
{
	View		view;
	
	view[0] = cos(header->reg[PHI_BP].f)*sin(header->reg[THE_BP].f);
	view[1] = sin(header->reg[PHI_BP].f)*sin(header->reg[THE_BP].f);
	view[2] = cos(header->reg[THE_BP].f);
	view[3] = angle_set_negPI_to_PI(header->reg[PSI_BP].f + header->reg[PHI_BP].f);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG view_from_imagic_euler: " << view << endl; 
 
	return view;
}

int			suprim_euler_from_view(SUPRIMhead* header, View view)
{
	if ( acos(view[2]) >  1e-14 )
		header->reg[PHI_BP].f = atan2(view[1], view[0]);
	header->reg[THE_BP].f = acos(view[2]);
	header->reg[PSI_BP].f = view[3] - header->reg[PHI_BP].f;
	header->reg[PSI_BP].f = angle_set_negPI_to_PI(header->reg[PSI_BP].f);
	
	return 0;
}


/**
@brief	Reading a SUPRIM map image file format.
A 3D image format used in electron microscopy.
	Header size:				36 bytes basic fields
								512 bytes additional fields called "registers"
								up to 1024 bytes of "trace" records
	File format extensions:  	.spm, .sup
	Byte order determination:	Version and third dimension values
								must be less than 256*256.
	Data types: 				1 = byte, 2 = float, 3 = complex float,
								4 = run-length-encoded float,
								5 = short, 6 = RGB, 7 = int, 8 = double.
	Transform type: 			Standard
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int					error code (<0 means failure).
**/
int 	readSUPRIM(Bimage *p, int readdata)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

	char		version;
	fimg->read(&version, 1);
	if ( fimg->fail() ) {
		error_show("Suprim version number incorrect", __FILE__, __LINE__);
		return -2;
	}
	
	SUPRIMhead*	header = new SUPRIMhead;
	memset(header, 0, sizeof(SUPRIMhead));
	fimg->read((char *)header, SUPRIMSIZE);
	if ( fimg->fail() ) {
		error_show("Suprim header not read", __FILE__, __LINE__);
		return -2;
	}
	
    // Determine byte order and swap bytes if from little-endian machine
    unsigned char*	b = (unsigned char *) header;
    int				i, sb = 0;
    if ( version < 0x40 || abs(header->reg[NSLICES].l) > SWAPTRIG ) {
		if ( verbose & VERB_PROCESS )
    		cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    	sb = 1;
		int 	extent = SUPRIMSIZE;	// Swap whole header
    	for ( i=0; i<extent; i+=4 ) swapbytes(b+i, 4);
    }
    
    // Convert VAX floating point types if necessary
/*    if ( header->min > header->max ) {
		if ( verbose & VERB_PROCESS )
    		cerr << "Warning: Converting VAX floating point" << endl;
    	vax = 1;
		for ( i=20; i<36; i+=4 )		// All floating point parameters
    	    vax2ieee(b+i, sb);
    }
*/	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readSUPRIM: version byte = " << version << endl;
		cout << "DEBUG readSUPRIM: data size = " << header->ncol << "," << header->nrow << "," << header->reg[NSLICES].l << endl;
	}
	
	// Trace records
	unsigned char*	ptr = (unsigned char *) header->trace;
	int*			iptr = (int *) ptr;
//	fread( ptr, 4, 1, fimg );
	fimg->read((char *)ptr, 4);
	if ( sb ) swapbytes(ptr, 4);
	int 	it = sizeof(int);
	int		narg = *iptr + 2;
	if ( verbose & VERB_DEBUG ) cout << "DEBUG readSUPRIM: number of trace arguments = " << narg << endl;
	if ( narg > 256 ) {
		cerr << "Error: The number of trace arguments are too large! (" << narg << ")" << endl;
		cerr << "	First trace bytes: " << ptr[0] << " " << ptr[1] << " " << ptr[2] << " " << ptr[3] << endl;
		return -3;
	}
	for ( i=0; i<narg; i++) {
//		fread( &ptr[it], 4, 1, fimg );
		fimg->read((char *)&ptr[it], 4);
		if ( sb ) swapbytes(&ptr[it], 4);
		iptr = (int *) &ptr[it];
		it += sizeof(int);
//		fread( &ptr[it], *iptr, 1, fimg );
		fimg->read((char *)&ptr[it], *iptr);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readSUPRIM: trace " << i+1 << " (" << *iptr << "): " << &ptr[it] << endl;
		if ( i == 2 ) p->label((char *) &ptr[it]);
		if ( i==2 && (verbose & VERB_DEBUG) )
			cout << "DEBUG readSUPRIM: label: " << p->label() << endl;
		it += *iptr;
	}
	it += sizeof(int);
	p->data_offset(SUPRIMSIZE + it + 1);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSUPRIM: data offset = " << p->data_offset() << endl;
    
	// Map the parameters
	p->size(header->ncol, header->nrow, header->reg[NSLICES].l);
	if ( header->type == 2 ) p->fourier_type(Standard);
	p->compound_type(TSimple);
	p->images(1);
	p->channels(1);
	switch ( header->intern ) {
		case 1: p->data_type(UCharacter); break; 	// BINARY
		case 2: p->data_type(Float); break;
		case 3: p->data_type(Float); p->channels(2); p->compound_type(TComplex); break;
		case 4: p->data_type(Float); break; 	// RLC_DATA
		case 5: p->data_type(Short); break;
		case 6: p->data_type(UCharacter); p->channels(3); p->compound_type(TRGB); break;
		case 7: p->data_type(Integer); break;
		case 8: p->data_type(Double); break;		// DOUBLE
		default: p->data_type(UCharacter); break;
	}
	p->sampling(header->reg[SAMP_DIST].f, header->reg[SAMP_DIST].f, header->reg[SAMP_DIST].f);
	p->minimum(header->min);
	p->maximum(header->max);
	p->average(header->av);
	p->standard_deviation(header->sd);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSUPRIM: data type = " << p->data_type() << endl;
	
//	p->image->origin(header->reg[COL_ORG].f, header->reg[ROW_ORG].f, header->reg[SLICE_ORG].f);
	p->origin(header->reg[COL_ORG].f, header->reg[ROW_ORG].f, header->reg[SLICE_ORG].f);
	p->background(0, header->reg[BACKGROUND].f);
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSUPRIM: background = " << p->background(long(0)) << endl;
	
	p->image->view(view_from_suprim_euler(header));

	delete header;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSUPRIM: Header setup done" << endl;

	long	readsize = p->alloc_size();
		
	if ( readdata ) {
		p->data_alloc();
		fread_large(p->data_pointer(), readsize, p->data_offset(), fimg);
		if ( sb ) swapbytes(readsize, p->data_pointer(), p->data_type_size());
	}
	
	fimg->close();
	delete fimg;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSUPRIM: Data read" << endl;

	return 0;
}

/**
@brief	Writing a SUPRIM map image file format.
A 3D image format used in electron microscopy.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeSUPRIM(Bimage *p)
{
	if ( p->compound_type() == TComplex ) p->change_type(Float);
	
    switch ( p->data_type() ) {
		case Bit: case UCharacter: case SCharacter:
			p->change_type(UCharacter); break;
    	case UShort: case Short:
			p->change_type(Short); break;
    	case UInteger: case Integer:
			p->change_type(Integer); break;
    	default: break;
    }
	
	SUPRIMhead*	header = new SUPRIMhead;
	memset(header, 0, sizeof(SUPRIMhead));
	
	// Map the parameters
	char	version = 0x43;
	header->ncol = p->sizeX();
	header->nrow = p->sizeY();
	header->reg[NSLICES].l = p->sizeZ();
	header->reg[SAMP_DIST].f = p->sampling(0)[0];
	header->type = 1;						// Image
	if ( p->fourier_type() ) header->type = 2;	// Fourier transform
	switch ( p->data_type() ) {
		case UCharacter: header->intern = 1; 		// BINARY
			if ( p->compound_type() == TRGB && p->channels() == 3 )
				header->intern = 6;
			break; 	// RGB
		case Float: header->intern = ( p->compound_type() == TComplex )? 3: 2; break;
//		case Float: header->intern = 4; break; 	// RLC_DATA
		case Short: header->intern = 5; break;
		case Integer: header->intern = 7; break;
		case Double: header->intern = 8; break;	// DOUBLE
		default: header->intern = 1; break;
	}
	header->min = p->minimum();
	header->max = p->maximum();
	header->av = p->average();
	header->sd = p->standard_deviation();
	
	header->reg[COL_ORG].f = p->image->origin()[0];
	header->reg[ROW_ORG].f = p->image->origin()[1];
	header->reg[SLICE_ORG].f = p->image->origin()[2];
	header->reg[BACKGROUND].f = p->background(long(0));
		
	suprim_euler_from_view(header, p->image->view());
			
	// Trace records
	char*			ptr = (char *) header->trace;
	int*			iptr = (int *) ptr;
	iptr[0] = 1;
	iptr[1] = 12;
	memcpy(ptr + 8, "Bsoft2.0", 8);
	iptr = (int *) (ptr + 20);
	*iptr = 24;
	memcpy(ptr + 24, asctime(p->get_localtime())+4, 20);
	iptr = (int *) (ptr + 48);
//	*iptr = 8;
//	memcpy(ptr + 52, "unknown", 7);
//	p->data_offset(SUPRIMSIZE + 65);
	long			llen = p->label().length();
	if ( llen > 800 ) llen = 800;
	*iptr = llen;
	memcpy(ptr + 52, p->label().c_str(), llen);
	p->data_offset(SUPRIMSIZE + llen + 57);
	
	long 			datatypesize = p->data_type_size();
	header->format = datatypesize;
	long 			datasize = p->alloc_size();
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG writeSUPRIM: version byte=" << version << endl;
		cout << "DEBUG writeSUPRIM: size:" << header->ncol << "," << header->nrow << "," << header->reg[NSLICES].l << endl;
		cout << "DEBUG writeSUPRIM: offset=" << p->data_offset() << " data size=" << datasize << endl;
		cout << "DEBUG writeSUPRIM: number of trace arguments = " << iptr[0] << endl;
		cout << "DEBUG writeSUPRIM: trace = ";
		for ( int i=0; i<20; i++ ) cout << i << " " << header->trace[i] << " " << header->trace[i] << endl;
		cout << endl;
	}
    
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write(&version, 1);
	fimg.write((char *)header, p->data_offset() - 1);
	
	if ( p->data_pointer() ) fimg.write((char *)p->data_pointer(), datasize);
	
	fimg.close();
	
	delete header;
	
	return 0;
}

