/**
@file	rwIMAGIC.cpp
@brief	Functions for reading and writing Image Science's Imagic files
@author Bernard Heymann
@date	Created: 19990424
@date 	Modified: 20180329
**/

#include "rwIMAGIC.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

View		view_from_imagic_euler(IMAGIChead* header)
{
	View		view;
	
	view[0] = cos(header->gamma)*sin(header->beta);
	view[1] = sin(header->gamma)*sin(header->beta);
	view[2] = cos(header->beta);
	view[3] = angle_set_negPI_to_PI(header->alpha + header->gamma);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG view_from_imagic_euler: " << view << endl; 
	
	return view;
}

int			imagic_euler_from_view(IMAGIChead* header, View view)
{
	if ( acos(view[2]) >  1e-14 )
		header->gamma = atan2(view[1], view[0]);
	header->beta  = acos(view[2]);
	header->alpha = view[3] - header->gamma;
	header->alpha = angle_set_negPI_to_PI(header->alpha);
	
	return 0;
}

/**
@brief	Reading an IMAGIC image format.
A 2D file format for the IMAGIC package.
	The header is stored in a separate file with extension ".hed" and
		a fixed size of 1024 bytes per image.
	The image data is stored in a single block in a file with the
		extension ".img".
	Machine stamp (4 bytes or integer):
		VAX/VMS			1 0 0 0 (16777216)
		Little endian	2 2 2 2 (33686018)
		Big endian		4 4 4 4 (67372036)
	Byte order determination:	Year and hour values
								must be less than 256*256.
	Data types: 				PACK = byte, INTG = short, REAL = float,
								RECO,COMP = complex float.
	Transform type: 			Centered (COMP data type)
								RECO is not a transform
	Note that the x and y dimensions are interchanged (actually a display issue).
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@param 	img_select		image selection in multi-image file (-1 = all images).
@return	int					error code (<0 means failure).
**/
int 	readIMAGIC(Bimage* p, int readdata, int img_select)
{
    // get the filename without extension to find the header file
	string			filename(p->file_name());
	string			headername = replace_extension(filename, "hed");
	filename = replace_extension(filename, "img");
	p->file_name(filename);

	if ( verbose & VERB_LABEL )
	    cout << "Reading header file:            " << headername << endl;
	
    // open and read the header file
	ifstream*		fhed = new ifstream(headername);
	if ( fhed->fail() ) return -1;	
	
	IMAGIChead*	header = new IMAGIChead;
	
	fhed->read((char *)header, IMAGICSIZE);
	if ( fhed->fail() ) return -2;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readIMAGIC: sizeof(IMAGIChead)=" << sizeof(IMAGIChead) << " IMAGICSIZE=" << IMAGICSIZE << endl;

    // Determine byte order and swap bytes if from little-endian machine
    unsigned char*	b = (unsigned char *) header;
    int				sb = 0;
    long		i, extent = IMAGICSIZE - 916;  // exclude char bytes from swapping
    if ( ( abs(header->nyear) > SWAPTRIG ) || ( header->ixlp > SWAPTRIG ) ) {
		if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    	sb = 1;
    	for ( i=0; i<extent; i+=4 ) 
    	    if ( i != 56 )	    	    // exclude type string
				swapbytes(b+i, 4);
    }
    
	// Map the parameters
	long 			nimg(0);
	if ( header->i4lp < 1 ) {	// This deals with the old header
		p->size(header->iylp, header->ixlp, 1);
		nimg = header->ifn + 1;
	} else {
		p->size(header->iylp, header->ixlp, header->izlp);
		nimg = header->i4lp;
	}
    
	p->channels(1);

	long			j(0), n(0), imgstart(0), imgend(nimg-1);
	double			min(1e37), max(-1e37), sum(0), ssum(0);
	unsigned char*	hend;
	if ( img_select > -1 ) {
		if ( img_select >= nimg ) img_select = nimg - 1;
		p->images(1);
		imgstart = imgend = img_select;
	} else {
		img_select = -1;
		p->images(nimg);
	}

	if ( strstr(header->type,"PACK") ) p->data_type(UCharacter);
    if ( strstr(header->type,"INTG") ) p->data_type(Short);
    if ( strstr(header->type,"REAL") ) p->data_type(Float);
    if ( strstr(header->type,"RECO") ) {
    	p->data_type(Float);	// Complex data
		p->compound_type(TComplex); p->channels(2);
		p->fourier_type(NoTransform);
	}
    if ( strstr(header->type,"COMP") ) {
    	p->data_type(Float);	// Complex transform data
		p->compound_type(TComplex); p->channels(2);
		p->fourier_type(Centered);
	}
	
	p->sampling(header->resolx, header->resoly, header->resolz);
    
    // Set min-max values and other statistical values
    if ( header->sigma == 0 && header->varian != 0 )
    	header->sigma = sqrt(header->varian);
    if ( header->densmax == 0 && header->densmin == 0 && header->sigma != 0 ) {
    	header->densmin = header->avdens - header->sigma;
    	header->densmax = header->avdens + header->sigma;
    }
    p->minimum(header->densmin);
    p->maximum(header->densmax);
	p->average(header->avdens);
	p->standard_deviation(header->sigma);
	
	p->label(header->history);
	
	tm*			t = p->get_localtime();
	t->tm_mday = header->ndate;
    t->tm_mon = header->nmonth - 1;
    t->tm_year = header->nyear - 1900;
    t->tm_hour = header->nhour;
    t->tm_min = header->nminut;
    t->tm_sec = header->nsec;
	p->set_time(t);
	
	
	// Get the sub-image information
	
	fhed->seekg(0, ios::beg);
	for ( i=n=0; i<nimg; i++ ) {
		fhed->read((char *)header, IMAGICSIZE);
		if ( fhed->fail() ) return -3;
		if ( img_select < 0 || img_select == i ) {
			hend = (unsigned char *) header + extent;
			if ( sb ) for ( b = (unsigned char *) header; b<hend; b+=4 ) swapbytes(b, 4);
			if ( header->sigma == 0 && header->varian != 0 )
				header->sigma = sqrt(header->varian);
			if ( header->densmax == 0 && header->densmin == 0 && header->sigma != 0 ) {
				header->densmin = header->avdens - header->sigma;
				header->densmax = header->avdens + header->sigma;
			}
			p->image[j].minimum(header->densmin);
			p->image[j].maximum(header->densmax);
			p->image[j].average(header->avdens);
			p->image[j].standard_deviation(header->sigma);
//			p->image[j].ox = -header->iyold;
//			p->image[j].oy = -header->ixold;
			p->image[j].origin(-header->nxstart, -header->nystart, -header->nzstart);	// New header
			p->image[j].view(view_from_imagic_euler(header));
			j++;
			if ( min > header->densmin ) min = header->densmin;
			if ( max < header->densmax ) max = header->densmax;
			sum += header->avdens;
			ssum += header->sigma*header->sigma;
			n++;
		}
	}
	
	p->minimum(min);
	p->maximum(max);
	p->average(sum/n);
	p->standard_deviation(sqrt(ssum/n));
	
	fhed->close();
	delete fhed;
	    
	delete header;

	if ( !readdata ) return 0;
	
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	
	
	long	pagesize = p->channels()*p->sizeX()*p->sizeY()*p->sizeZ()*p->data_type_size();
	long	offset(p->data_offset() + imgstart*pagesize);
	unsigned char*	data = NULL;

	if ( readdata ) {
		p->data_alloc();	
		if ( p->compound_type() == TSimple ) {
			for ( i=imgstart, data=p->data_pointer(); i<=imgend; i++, data+=pagesize, offset+=pagesize ) {
				fread_large(data, pagesize, offset, fimg);
			}
			if ( sb ) swapbytes(p->alloc_size(), p->data_pointer(), p->data_type_size());
		} else {
			data = new unsigned char[pagesize];
			for ( i=imgstart; i<=imgend; i++, offset+=pagesize ) {
				fread_large(data, pagesize, offset, fimg);
				if ( sb ) swapbytes(pagesize, data, p->data_type_size());
				p->unpack_transform(i-imgstart, data, Centered);
			}
			delete[] data;
		}
	}
	
	fimg->close();
	delete fimg;
	    
	return 0;
}

/**
@brief	Writing an IMAGIC image format.
A file format for the IMAGIC package.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeIMAGIC(Bimage* p)
{
	p->color_to_simple();
	
	if ( p->compound_type() == TComplex ) p->change_type(Float);
	
    switch ( p->data_type() ) {
		case Bit: case UCharacter: case SCharacter:
			p->change_type(UCharacter); break;
    	case UShort: case Short:
			p->change_type(Short); break;
    	default:
			p->change_type(Float); break;
    }
	
	IMAGIChead*	header = new IMAGIChead;
	memset(header, 0, sizeof(IMAGIChead));
	
	// Machine stamp
	long			i, n;		
	char			machine_stamp[4] = {0,0,0,0};
	int*			rt = (int *) machine_stamp;
	switch ( systype(0) ) {
		case BigIEEE:
			for ( i=0; i<4; i++ ) machine_stamp[i] = 4;
			break;
		case LittleIEEE:
			for ( i=0; i<4; i++ ) machine_stamp[i] = 2;
			break;
		case LittleVAX:
			machine_stamp[0] = 1;
			break;
		default:
			break;
	}
//	header->realtype = *((int *) machine_stamp);
	header->realtype = *rt;
	
    // Fill in the file header
	header->nhfr = 1;
	header->npix2 = p->sizeX()*p->sizeY();
	header->npixel = p->sizeX()*p->sizeY();
    header->iylp = p->sizeX();
    header->ixlp = p->sizeY();
    header->izlp = p->sizeZ();
    header->i4lp = p->images();
	
	header->resolx = p->sampling(0)[0];
	header->resoly = p->sampling(0)[1];
	header->resolz = p->sampling(0)[2];
	
	tm*			t = p->get_localtime();
    header->ndate = t->tm_mday;
    header->nmonth = t->tm_mon + 1;
    header->nyear = t->tm_year + 1900;
    header->nhour = t->tm_hour;
    header->nminut = t->tm_min;
    header->nsec = t->tm_sec;

    switch( p->data_type() ){
        case UCharacter : memcpy(header->type,"PACK", 4); break;
        case Short : memcpy(header->type,"INTG", 4); break;
        case Float : memcpy(header->type,"REAL", 4);
			if ( p->compound_type() == TComplex ) memcpy(header->type,"RECO", 4);
			if ( p->fourier_type() >= Standard ) memcpy(header->type,"COMP", 4);
			break;
        default : memcpy(header->type,"PACK", 4);
    }
    
	header->densmin = p->minimum();
	header->densmax = p->maximum();
	header->avdens = header->oldavd = p->average();
	header->sigma = p->standard_deviation();
	header->varian = p->standard_deviation()*p->standard_deviation();

	strncpy(header->name, p->file_name().c_str(), 80);
	strncpy(header->history, p->label().c_str(), 228);
	
    // get the filename without extension to find the header file
	string			filename(p->file_name());
	string			headername = replace_extension(filename, "hed");
	filename = replace_extension(filename, "img");
	p->file_name(filename);

	if ( verbose & VERB_LABEL )
	    cout << "Writing header file:            " << headername << endl;

	ofstream		fhed(headername);
	if ( fhed.fail() )  return -1;

	for ( n=0; n<p->images(); n++ ) {
		header->imn =	n + 1;			// Image number
		header->ifn =	p->images() - n - 1;	// Number of images following this one
		header->iyold = (int) -(p->image[n].origin()[0] + 0.5);
		header->ixold = (int) -(p->image[n].origin()[1] + 0.5);
		header->nxstart = (int) -(p->image[n].origin()[0] + 0.5);
		header->nystart = (int) -(p->image[n].origin()[1] + 0.5);
		header->nzstart = (int) -(p->image[n].origin()[2] + 0.5);
		header->densmin = p->image[n].minimum();
		header->densmax = p->image[n].maximum();
		header->avdens = header->oldavd = p->image[n].average();
		header->sigma = p->image[n].standard_deviation();
		header->varian = p->image[n].standard_deviation()*p->image[n].standard_deviation();
		imagic_euler_from_view(header, p->image[n].view());
		fhed.write((char *)header, IMAGICSIZE);
	}
	
	fhed.close();
	
	delete header;

	if ( !p->data_pointer() ) return -2;
		
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -3;

	long 			datasize = p->alloc_size();
	unsigned char*	data = NULL;

	if ( p->compound_type() < TComplex ) {
		fimg.write((char *)p->data_pointer(), datasize);
	} else {
		data = new unsigned char[datasize];
		p->pack_transform(data, Centered);
		fimg.write((char *)data, datasize);
		delete[] data;
	}
	
	fimg.close();
	
	return 0;
}

