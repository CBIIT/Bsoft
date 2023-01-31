/**
@file	rwSPIDER.cpp
@brief	Functions for reading and writing SPIDER files
@author Bernard Heymann
@date	Created: 19990410
@date 	Modified: 20120321
**/

#include "utilities.h"
#include "file_util.h"
#include "rwSPIDER.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

extern char* month[];

View		view_from_spider_euler(SPIDERhead* header)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG view_from_spider_euler: " 
			   << header->gamma << "," << header->theta << "," << header->phi << endl;
	
	View	view;
	
	view[0] = cos(header->phi)*sin(header->theta);
	view[1] = sin(header->phi)*sin(header->theta);
	view[2] = cos(header->theta);
	view[3] = angle_set_negPI_to_PI(header->gamma + header->phi);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG view_from_spider_euler: " << view << endl; 
	
	return view;
}

int			spider_euler_from_view(SPIDERhead* header, View view)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG spider_euler_from_view: " << view << endl; 
	
	if ( acos(view[2]) >  1e-14 )
		header->phi = atan2(view[1], view[0]);
	header->theta = acos(view[2]);
	header->gamma = view[3] - header->phi;
	header->gamma = angle_set_negPI_to_PI(header->gamma);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG spider_euler_from_view: " 
			   << header->gamma << "," << header->theta << "," << header->phi << endl;
 	
	return 0;
}

/**
@brief	Reading a SPIDER image file format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@param 	img_select	image selection in multi-image file (-1 = all images).
@return	int			error code (<0 means failure).
A 3D multi-image format used in electron microscopy.
	Header size:				1024 bytes (not same as data offset!).
	Data offset:				sizeof(float)*x_size*ceil(1024/x_size)
	File format extensions:  	.spi
	The identifier is a 4-byte machine stamp:
				1	Big-endian IEEE 	17 17 00 00
                2	VAX 				34 65 00 00
				3	Cray				-
                4	Little-endian IEEE	68 65 00 00
                5	Convex				85 17 00 00
				6	Fijitsu VP			-
				(Note: not always implemented - so unreliable)
	Byte order determination:	File type and third dimension values
								must be less than 256*256.
	Data type: 					only float.
	Transform type: 			Hermitian
								The x-dimension contains the x-size
								of the full transform
	A multi-image file has a global header followed by a header and data
	for each sub-image.
**/
int 	readSPIDER(Bimage* p, int readdata, int img_select)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	
	
	SPIDERhead*	header = new SPIDERhead;
	
	fimg->read((char *)header, SPIDERSIZE);
	if ( fimg->fail() ) return -2;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readSPIDER: Spider header size = " << sizeof(SPIDERhead) << endl;
		cout << "DEBUG readSPIDER: iform=" << header->iform << endl;
	}
	
    // Determine byte order and swap bytes if from different-endian machine
    unsigned char*	b = (unsigned char *) header;
    int				i, j, sb = 0;
    int				extent = SPIDERSIZE - 180;  // exclude char bytes from swapping
    if ( ( fabs(header->nslice) > SWAPTRIG ) || ( fabs(header->iform) > SWAPTRIG ) ||
    	( fabs(header->nslice) < 1 ) ) {
    	if ( verbose & VERB_PROCESS )
			cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    	sb = 1;
    	for ( i=0; i<extent; i+=4 ) swapbytes(b+i, 4);
    }
        
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSPIDER: nsam=" << header->nsam << " nrow=" << header->nrow 
			<< " nslice=" << header->nslice << " istack=" << header->istack << endl;
	
	// Map the parameters
	p->size((long) header->nsam, (long) header->nrow, (long) header->nslice);
	p->channels(1);
	p->data_type(Float);

	long		xstore = p->sizeX();
	if ( header->iform < 0 ) {
		p->compound_type(TComplex);
		p->channels(2);
		xstore = (int) header->nsam/2;
		if ( header->iform > -20 ) {
			p->sizeX(p->sizeX() + (int) (10 + header->iform));	// 2D Fourier transform
		} else {
			p->sizeX(p->sizeX() + (int) (20 + header->iform));	// 3D Fourier transform
		}
	}
	
	long				nimg(1);
	if ( header->istack > 0 ) nimg = header->maxim;
	
	long 				imgstart(0);
	long 				imgend(nimg - 1);
	unsigned char*		hend;

	if ( img_select > -1 ) {
		if ( img_select >= nimg ) img_select = nimg - 1;
		imgstart = img_select;
		imgend = img_select;
		nimg = 1;
	}
	
	p->images(nimg);

	p->page_size(p->size());
	p->data_offset((long) header->labbyt);
	p->minimum(header->fmin);
	p->maximum(header->fmax);
	p->average(header->av);
	p->standard_deviation(header->sig);
	p->sampling(header->scale, header->scale, header->scale);
	
	tm*			t = p->get_localtime();
	sscanf(header->ctim, "%d:%d:%d", &t->tm_hour, &t->tm_min, &t->tm_sec);
    sscanf(header->cdat, "%d", &t->tm_mday);
	for ( i=0; i<12 && strncmp(&(header->cdat[3]),month[i],3); i++ ) ;
	if ( i < 12 ) t->tm_mon = i;
	sscanf(&(header->cdat[7]), "%d", &t->tm_year);
	if ( t->tm_year < 70 ) t->tm_year += 100;
	p->set_time(t);

	p->label(header->ctit);
	
	long 		header_size = p->data_offset();
	long 		image_size = header_size + p->channels()*xstore*p->sizeY()*p->sizeZ()*p->data_type_size();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSPIDER: header_size=" << header_size << " image_size=" << image_size << endl;
	
	
	p->origin(header->xoff, header->yoff, header->zoff);
	p->image->minimum(header->fmin);
	p->image->maximum(header->fmax);
	p->image->average(header->av);
	p->image->standard_deviation(header->sig);
	p->image->view(view_from_spider_euler(header));
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSPIDER: imgstart=" << imgstart << " imgend=" << imgend << endl;

	if ( header->istack > 0 ) {
		p->data_offset(p->data_offset() + header_size);
		for ( i=imgstart, j=0; i<=imgend; i++, j++ ) {
			fimg->seekg(header_size + i*image_size, ios::beg);
			fimg->read((char *)header, SPIDERSIZE);
			if ( fimg->fail() ) {
				cerr << "Error: Reading image " << i+1 << " header failed!" << endl;
				fimg->close();
				delete fimg;
				return -3;
			}
			hend = (unsigned char *) header + extent;
			if ( sb ) for ( b = (unsigned char *) header; b<hend; b+=4 ) swapbytes(b, 4);
			p->image[j].minimum(header->fmin);
			p->image[j].maximum(header->fmax);
			p->image[j].average(header->av);
			p->image[j].standard_deviation(header->sig);
			p->image[j].origin(header->xoff, header->yoff, header->zoff);
			p->image[j].view(view_from_spider_euler(header));
		}
	}

	delete header;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSPIDER: offset=" << p->data_offset() << endl;
		
	long			pagesize = p->channels()*xstore*p->sizeY()*p->sizeZ()*p->data_type_size();
	long			offset(p->data_offset() + imgstart*image_size);
	unsigned char*	data = NULL;

	if ( readdata ) {
		p->data_alloc();	
		if ( p->compound_type() == TSimple ) {
			for ( i=imgstart, data=p->data_pointer(); i<=imgend; i++, data+=pagesize, offset+=image_size ) {
				fread_large(data, pagesize, offset, fimg);
			}
			if ( sb ) swapbytes(p->alloc_size(), p->data_pointer(), p->data_type_size());
		} else {
			data = new unsigned char[pagesize];
			for ( i=imgstart; i<=imgend; i++, offset+=image_size ) {
				fread_large(data, pagesize, offset, fimg);
				if ( sb ) swapbytes(pagesize, data, p->data_type_size());
				p->unpack_transform(i-imgstart, data, Hermitian);
			}
			delete[] data;
		}
	}
	
	fimg->close();
	delete fimg;
	
	return 0;
}

/**
@brief	Writing a SPIDER image file format.
A 3D image format used in electron microscopy.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeSPIDER(Bimage* p)
{
	p->color_to_simple();
	p->change_type(Float);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeSPIDER: Writing Spider file" << endl;

    float		lenbyt = sizeof(float)*p->sizeX();				// Record length (in bytes)
    float		labrec = floor(SPIDERSIZE/lenbyt);			// # header records
    if ( fmod(SPIDERSIZE,lenbyt) != 0 ) labrec++;
    float		labbyt = labrec*lenbyt; 					// Size of header in bytes
	p->data_offset((long) labbyt);

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG writeSPIDER: Header size=" << sizeof(SPIDERhead) << " labbyt=" << labbyt << endl;
		cout << "DEBUG writeSPIDER: labrec=" << labrec << " lenbyt=" << lenbyt << endl;
	}
	
	char*		buf = new char[p->data_offset()];
	SPIDERhead*	header = (SPIDERhead *) buf;
	memset(buf, 0, p->data_offset());
	
	// Map the parameters
    header->lenbyt = lenbyt;									// Record length (in bytes)
    header->labrec = labrec;									// # header records
    header->labbyt = labbyt; 									// Size of header in bytes
    header->irec = labrec + floor((p->sizeX()*p->sizeY()*p->sizeZ()*sizeof(float))/lenbyt + 0.999999);	// Total # records
	header->nsam = p->sizeX();
	header->nrow = p->sizeY();
	header->nslice = p->sizeZ();
	header->scale = p->sampling(0)[0];	// Note this ignores any anisotropy

	// If a transform, then the physical storage in x is only half+1
	long		xstore = p->sizeX();
	if ( p->compound_type() == TComplex ) {
		xstore = p->sizeX()/2 + 1;
		header->nsam = 2*xstore;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeSPIDER: Size: " << header->nsam << "," << header->nrow << "," << header->nslice << endl;

	if ( p->sizeZ() < 2 ) {
		if ( p->compound_type() == TSimple )
    		header->iform = 1;				// 2D image
		else
			header->iform = -12 + (int)p->sizeX()%2;   // 2D Fourier transform
	} else {
		if ( p->compound_type() == TSimple )
    		header->iform = 3;				// 3D volume
		else
			header->iform = -22 + (int)p->sizeX()%2;   // 3D Fourier transform
	}
	header->istack = 0;
	header->imami = 1;
	header->fmin = p->minimum();
	header->fmax = p->maximum();
	header->av = p->average();
	header->sig = p->standard_deviation();
	
	// For multi-image files
	if ( p->images() > 1 ) {
		header->istack = 2;
		header->inuse = -1;
		header->maxim = p->images();
	}
	
	header->xoff = p->image->origin()[0];
	header->yoff = p->image->origin()[1];
	header->zoff = p->image->origin()[2];
	
	spider_euler_from_view(header, p->image->view());
			
    tm* 		t = p->get_localtime();
 	strftime (header->cdat, 12, "%d-%b-%y", t);
 	strftime (header->ctim, 8, "%H:%M:%S", t);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeSPIDER: date=" << header->cdat << " time=" << header->ctim << endl;
	strncpy(header->ctit, p->label().c_str(), 160);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG writeSPIDER: Date and time: " << header->cdat << " " << header->ctim << endl;
		cout << "DEBUG writeSPIDER: Text label: " << header->ctit << endl;
	}

	long			n;
	long			datatypesize = p->channels()*p->data_type_size();
	long 			datasize = xstore*p->sizeY()*p->sizeZ()*datatypesize;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG writeSPIDER: type size=" << datatypesize << " data size=" << datasize << " offset=" << p->data_offset() << endl;
		cout << "DEBUG writeSPIDER: file name=" << p->file_name() << endl;
	}
	
	unsigned char*	data = NULL;
 	unsigned char*	aptr = p->data_pointer();
	if ( p->compound_type() == TComplex ) data = new unsigned char[datasize];

	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeSPIDER: file " << p->file_name() << " opened" << endl;
	
	fimg.write((char *)header, p->data_offset());
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeSPIDER: File header written" << endl;
	
	if ( p->images() == 1 ) {
		if ( p->compound_type() < TComplex ) {
			fimg.write((char *)p->data_pointer(), datasize);
		} else {
			p->pack_transform(data, Hermitian);
			fimg.write((char *)data, datasize);
		}
	} else {
		header->istack = 0;
		header->inuse = 0;
		header->maxim = 0;
		for ( n=0; n<p->images(); n++ ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG writeSPIDER: Writing image number " << n << endl;
			header->imgnum = n + 1;
			header->fmin = p->image[n].minimum();
			header->fmax = p->image[n].maximum();
			header->av = p->image[n].average();
			header->sig = p->image[n].standard_deviation();
			header->xoff = p->image[n].origin()[0];
			header->yoff = p->image[n].origin()[1];
			header->zoff = p->image[n].origin()[2];
			spider_euler_from_view(header, p->image[n].view());
			fimg.write((char *)header, p->data_offset());
			if ( p->compound_type() < TComplex ) {
				fimg.write((char *)aptr, datasize);
			} else {
				p->pack_transform(n, data, Hermitian);
				fimg.write((char *)data, datasize);
			}
			aptr += datasize;
		}
	}
	
	fimg.close();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeSPIDER: File written" << endl;
	
	if ( data ) delete[] data;
	
	delete[] buf;
	
	return 0;
}

