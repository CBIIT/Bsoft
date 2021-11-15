/**
@file	rwSER.cpp
@brief	Reading FEI SER files
@author Bernard Heymann
@date	Created: 20150130
@date	Modified: 20150630
**/

#include <time.h>

#include "rwSER.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			read_field(ifstream* fimg, char* p, size_t n)
{
	size_t		where = fimg->tellg();
	
	fimg->read(p, n);
	
	if ( fimg->fail() ) {
		cerr << "Reading failed at " << where << endl;
		return -2;
	}
	
	return n;
}

int 		readSER(Bimage* p, int readdata)
{
	ifstream*			fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSER: File opened" << endl;

	int					offset, iv, nx, ny;
	short				datatype;
	double				ox, oy, ux, uy;
	
	fimg->seekg(22);
	read_field(fimg, (char *) &offset, sizeof(int));
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSER: offset: " << offset << endl;

	fimg->seekg(offset);
	read_field(fimg, (char *) &offset, sizeof(int));
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSER: offset: " << offset << endl;

	fimg->seekg(offset);

	read_field(fimg, (char *) &ox, sizeof(double));
	read_field(fimg, (char *) &ux, sizeof(double));
	read_field(fimg, (char *) &iv, sizeof(int));

	read_field(fimg, (char *) &oy, sizeof(double));
	read_field(fimg, (char *) &uy, sizeof(double));
	read_field(fimg, (char *) &iv, sizeof(int));
	
	read_field(fimg, (char *) &datatype, sizeof(short));

	read_field(fimg, (char *) &nx, sizeof(int));
	read_field(fimg, (char *) &ny, sizeof(int));

	p->size(nx, ny, 1);
	p->channels(1);
	p->images(1);

	p->data_offset(offset+50);	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSER: data offset: " << p->data_offset() << endl;
	
	p->sampling(ux*1e10, uy*1e10, 1);
//	p->image->origin(-ox/ux, -oy/uy, 0);
	p->origin(-ox/ux, -oy/uy, 0.0);
	
	switch ( datatype ) {
		case 1: p->data_type(UCharacter); break;
		case 2: p->data_type(UShort); break;
		case 3: p->data_type(UInteger); break;
		case 4: p->data_type(SCharacter); break;
		case 5: p->data_type(Short); break;
		case 6: p->data_type(Integer); break;
		case 7: p->data_type(Float); break;
		case 8: p->data_type(Double); break;
		case 9: p->data_type(Float); p->compound_type(TComplex); break;
		case 10: p->data_type(Double); p->compound_type(TComplex); break;
	}
	
	if ( readdata ) {
		p->data_alloc();
		fread_large(p->data_pointer(), p->alloc_size(), p->data_offset(), fimg);
	}
	
	fimg->close();
	delete fimg;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwSER: Data read" << endl;
	
	return 0;
}

int 		writeSER(Bimage* p)
{
//	cerr << "Writing is not supported for this format!" << endl;

	if ( p->data_type() == ULong ) p->change_type(UInteger);
	if ( p->data_type() == Long ) p->change_type(Integer);
	
	p->data_offset(126);
	
	short			sv;
	int				iv;
//	float			fv;
	double			dv;
	char			str[8] = "Number";
	
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	sv = 0x4949;	// ByteOrder
	fimg.write((char *) &sv, sizeof(short));
	
	sv = 0;			// SeriesID
	fimg.write((char *) &sv, sizeof(short));
	
	sv = 0x0210;	// SeriesVersion
	fimg.write((char *) &sv, sizeof(short));
	
	iv = 0x4122;	// DataTypeID
	fimg.write((char *) &iv, sizeof(int));
	
	iv = 0x4152;	// TagTypeID
	fimg.write((char *) &iv, sizeof(int));
	
	iv = 1;			// TotalNumberElements
	fimg.write((char *) &iv, sizeof(int));
	
	iv = 1;			// ValidNumberElements
	fimg.write((char *) &iv, sizeof(int));
	
	iv = 68;		// OffsetArrayOffset
	fimg.write((char *) &iv, sizeof(int));
	
	iv = 1;			// NumberDimensions
	fimg.write((char *) &iv, sizeof(int));

	iv = 1;			// DimensionSize
	fimg.write((char *) &iv, sizeof(int));

	dv = 0;			// CalibrationOffset
	fimg.write((char *) &dv, sizeof(double));

	dv = 0;			// CalibrationDelta
	fimg.write((char *) &dv, sizeof(double));

	iv = 0;			// CalibrationElement
	fimg.write((char *) &iv, sizeof(int));

	iv = 6;			// DescriptionLength
	fimg.write((char *) &iv, sizeof(int));
	if ( iv ) fimg.write(str, iv);

	iv = 0;			// UnitsLength: 0
	fimg.write((char *) &iv, sizeof(int));
	if ( iv ) fimg.write(str, iv);

	iv = 76;		// DataOffset[0]: 76
	fimg.write((char *) &iv, sizeof(int));

	iv = p->data_offset() + p->alloc_size();	// TagOffset[0]: 67108990
	fimg.write((char *) &iv, sizeof(int));

	dv = -p->image->origin()[0]*p->sampling(0)[0]*1e-10;	// CalibrationOffsetX
	fimg.write((char *) &dv, sizeof(double));

	dv = p->sampling(0)[0]*1e-10;	// CalibrationDeltaX
	fimg.write((char *) &dv, sizeof(double));

	iv = 0;			// CalibrationElementX
	fimg.write((char *) &iv, sizeof(int));

	dv = -p->image->origin()[1]*p->sampling(0)[1]*1e-10;	// CalibrationOffsetY
	fimg.write((char *) &dv, sizeof(double));

	dv = p->sampling(0)[1]*1e-10;	// CalibrationDeltaY
	fimg.write((char *) &dv, sizeof(double));

	iv = 0;			// CalibrationElementY
	fimg.write((char *) &iv, sizeof(int));

	// DataType
	switch ( p->data_type() ) {
		case Bit: sv = 1; break;
		case UCharacter: sv = 1; break;
		case SCharacter: sv = 4; break;
		case UShort: sv = 2; break;
		case Short: sv = 5; break;
		case UInteger: sv = 3; break;
		case Integer: sv = 6; break;
//		case ULong: sv = 3; break;
//		case Long: sv = 4; break;
		case Float: sv = 7; break;
		case Double: sv = 8; break;
		default: break;
	}
	if ( p->compound_type() == TComplex ) sv += 2;
	fimg.write((char *) &sv, sizeof(short));
	
	iv = p->sizeX();	// ArraySizeX
	fimg.write((char *) &iv, sizeof(int));

	iv = p->sizeY();	// ArraySizeY
	fimg.write((char *) &iv, sizeof(int));

	// Data
	fimg.write((char *) p->data_pointer(), p->alloc_size());
	
	sv = 0x4152;		// TagTypeID
	fimg.write((char *) &sv, sizeof(short));
	
	time_t	t = time(NULL);		// Time
	fimg.write((char *) &t, sizeof(time_t));

	fimg.close();
	
	return 0;
}

