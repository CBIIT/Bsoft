/**
@file	rwPIF.cpp
@brief	Functions for reading and writing PIF files
@author Bernard Heymann
@date	Created: 19991112
@date 	Modified: 20120706
**/

#include "rwPIF.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

extern char* month[];

// The PIF format's internal scale factor
#define		PIF_SCALE	1e-5
double		PIFscale = PIF_SCALE;

// Internal function prototypes  
char*	read_sf(ifstream* fimg, Bimage* p, int mode, int sb);

/**
@brief	Reading a PIF image and structure factor file format.
A 2D and 3D image format used in electron microscopy.
	There are two types of files:
	1.	A typical binary image format with a file header and a header 
		for each sub-image.
	2.	A binary structure factor format with the same file header as 
		for the image format
	File header size:			512 bytes followed by optional 
								768 byte colour table
	Image header size:			512 bytes for each image 
	File format extensions:  	.pif
	Identifier: 				supposed to be the first 8 bytes
	Byte order determination:	Data type and third dimension values
								must be less than 256*256.
	Data types: 				0 = byte, 1 = short, 2 = float,
								3 = complex short, 4 = complex float.
	Transform type: 			List of structure factors.
	Note:						To avoid endianness problems associated 
								with different floating point formats
								(particularly those under VMS), this
								format was intended to store only integers.
								A 16 byte string following the 8th byte
								encodes a floating point conversion scale
								to convert an integer map and some of the
								fields in the header back to floating point.
								However, this format was designed in the
								CCP4/MRC style format, and the floating
								point data type was retained, as well as
								complex types used for transforms.
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@param 	img_select		image selection in multi-image file (-1 = all images).
@return	int					error code (<0 means failure).
**/
int 	readPIF(Bimage* p, int readdata, int img_select)
{
//	double		PIFscale = PIF_SCALE;
	
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;
	
	long			header_size = sizeof(PIFhead);
	long			image_header_size = sizeof(PIFimagehead);
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readPIF: Size of file header: " <<  header_size << endl;
		cout << "DEBUG readPIF: Size of image header: " << image_header_size << endl;
	}	

	PIFhead*	header = new PIFhead;
	
	fimg->read((char *)header, header_size);
	if ( fimg->fail() ) return -2;
	
    // Determine byte order and swap bytes if from little-endian machine
    unsigned char*		b = (unsigned char *) header;
    long				i, j, sb = 0;
    if ( ( abs( header->mode ) > SWAPTRIG ) || ( abs(header->nz) > SWAPTRIG ) ) {
    	if ( verbose & VERB_PROCESS )
			cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    	sb = 1;
		long 			extent = 84; // exclude labels from swapping
    	for ( i=0; i<extent; i+=4 ) {
			if ( i > 20 && i < 32 ) swapbytes(b+i, 4);
			if ( i > 60 ) swapbytes(b+i, 4);
		}
    }
    
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIF: Mode = " << header->mode << endl;
	
	// Map the parameters of the file header
	// All numerical parameters are integer and some must be converted
	// to floating point using a scaling factor
//	long	pad = 0;
	char			realscalefactor[20];
	strncpy(realscalefactor, header->realscalefactor, 16);
	realscalefactor[16] = 0;
	sscanf(realscalefactor, "%lf", &PIFscale);
//	p->scale = PIFscale;
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIF: PIF scale factor: " << realscalefactor << " " << PIFscale << endl;
	p->size(header->nx, header->ny, header->nz);
	p->page_size(p->size());
	p->channels(1);
	if ( header->numimages > 1 ) {
    	if ( verbose & VERB_DEBUG )
			printf( "DEBUG readPIF: Multi-image file with %d images!\n", 
					header->numimages);
//		pad = 512;
	}
	if ( header->numimages > 1 && header->htype < 1 ) {
		cerr << "Error: Images with different sizes!" << endl;
		return -1;
	}

	switch ( header->mode ) {		// Note: There are many more types!
		case 0: 
		case 6: p->data_type(UCharacter); break;
		case 1:
		case 7:
		case 20:
		case 88: p->data_type(Short); break;
		case 2:
		case 21:
		case 22:
		case 46: p->data_type(Integer); break;
		case 9: p->data_type(Float); break;
		case 40: p->data_type(Double); break;
		case 3:
		case 8: p->data_type(Short); p->compound_type(TComplex); p->channels(2); break;
		case 4:
		case 31:
		case 32: p->data_type(Integer); p->compound_type(TComplex); p->channels(2); break;
		case 10: p->data_type(Float); p->compound_type(TComplex); p->channels(2); break;
		case 41: p->data_type(Double); p->compound_type(TComplex); p->channels(2); break;
		default: p->data_type(UCharacter); break;
	}
	
	if ( header->even < 2) p->label(header->label);
	else p->label((char *) &header->even);
	
	p->data_offset(PIFSIZE);
	if ( header->mode == 97 ) {		// Depth-cued ===> colourmap
		PIFcolor*	colormap = new PIFcolor;
		fimg->read((char *)colormap, sizeof(PIFcolor));
		p->data_offset(p->data_offset() + sizeof(PIFcolor));
		delete colormap;
	}
	
	// Map the parameters in the image headers
	long   xstore = p->sizeX();
	if ( p->compound_type() == TComplex ) xstore = p->sizeX()/2 + 1;
	long   image_size = xstore*p->sizeY()*p->sizeZ()*p->channels()*p->data_type_size()
						+ image_header_size;
	PIFimagehead*	imageheader = new PIFimagehead;
	PIFimagehead_new* imageheader_new = (PIFimagehead_new *) imageheader;
	
	long			n(0), imgstart(0);
	long			imgend(header->numimages - 1);
	double			min(1e37), max(-1e37), sum(0), ssum(0);
	Vector3<double>	ori;
	
	if ( img_select > -1 ) {
		if ( img_select >= header->numimages ) img_select = header->numimages - 1;
		imgstart = img_select;
		imgend = img_select;
		p->images(1);
	} else {
		p->images(header->numimages);
	}
	
	tm*			t = p->get_localtime();
	char		mon[4];
	
	fimg->seekg(header_size + imgstart*image_size, ios_base::beg);

	for ( i = imgstart; i<=imgend; i++ ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readPIF: Reading image " << i+1 << " at position " 
					<< header_size + i*image_size << endl;
//		fimg->seekg(header_size + i*image_size, ios_base::beg);
		fimg->read((char *)imageheader, image_header_size);
		if ( fimg->fail() ) {
			cerr << "Error: Reading image " << i+1 << " header failed!" << endl;
			fimg->close();
			delete fimg;
			return -3;
		}
		
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readPIF: time string: " << imageheader->timestamp << endl;
		if ( strlen(imageheader->timestamp) ) {
			j = sscanf(imageheader->timestamp, "%3s %3s %2d %2d:%2d:%2d %4d", 
				mon, mon, &t->tm_mday, &t->tm_hour, &t->tm_min, &t->tm_sec, &t->tm_year);
			if ( j > 6 ) {
				for ( j=0; strncmp(mon, month[j], 3) != 0; j++ ) ;
				if ( j < 12 ) t->tm_mon = j;
				t->tm_year -= 1900;
				p->set_time(t);
			}
		}
		
    	if ( sb ) {
    		b = (unsigned char *) imageheader;
    		for ( j=0; j<368; j+=4 ) {
				if ( j < 116 || j > 248 ) swapbytes(b+j, 4);
			}
		}
		j = ( p->images() > 1 )? i: 0;
		p->background(j, imageheader->bkgnd);
		p->image[j].minimum(PIFscale*imageheader->min);
		p->image[j].maximum(PIFscale*imageheader->max);
		p->image[j].average(PIFscale*imageheader->mean);
		p->image[j].standard_deviation(PIFscale*imageheader->stddev);
		if ( abs(imageheader->nxstart) < 1e6 )
			ori[0] = -imageheader->nxstart;
		if ( abs(imageheader->nystart) < 1e6 )
			ori[1] = -imageheader->nystart;
		if ( abs(imageheader->nzstart) < 1e6 )
			ori[2] = -imageheader->nzstart;
		// Trigger for new image header
		if ( imageheader_new->ctfMode < 10 ) {
			if ( fabs(PIFscale*imageheader_new->origX) < 1e6 )
				ori[0] = PIFscale*imageheader_new->origX;
			if ( fabs(PIFscale*imageheader_new->origY) < 1e6 )
				ori[1] = PIFscale*imageheader_new->origY;
			if ( fabs(PIFscale*imageheader_new->origZ) < 1e6 )
				ori[2] = PIFscale*imageheader_new->origZ;
			p->image[j].view(PIFscale*imageheader_new->view_x, PIFscale*imageheader_new->view_y, 
					PIFscale*imageheader_new->view_z, PIFscale*imageheader_new->view_angle);
		} else {
			if ( fabs(PIFscale*imageheader->xorigin) < 1e6 )
				ori[0] = PIFscale*imageheader->xorigin;
			if ( fabs(PIFscale*imageheader->yorigin) < 1e6 )
				ori[1] = PIFscale*imageheader->yorigin;
			if ( fabs(PIFscale*imageheader->zorigin) < 1e6 )
				ori[2] = PIFscale*imageheader->zorigin;
			p->image[j].view(PIFscale*imageheader->view_x, PIFscale*imageheader->view_y, 
					PIFscale*imageheader->view_z, PIFscale*imageheader->view_angle);
		}
		p->image[j].origin(ori);
		if ( min > p->image[j].minimum() ) min = p->image[j].minimum();
		if ( max < p->image[j].maximum() ) max = p->image[j].maximum();
		sum += p->image[j].average();
		ssum += p->image[j].standard_deviation()*
			p->image[j].standard_deviation();
		n++;
		fimg->seekg(image_size-image_header_size, ios_base::cur);
    }
    
	if ( ( verbose & VERB_DEBUG ) && sb )
		cout << "DEBUG readPIF: Header swapping done" << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIF: min=" << imageheader->min << " max=" << 
			imageheader->max << " mean=" << imageheader->mean << " std=" << imageheader->stddev << endl;
	p->minimum(min);
	p->maximum(max);
	p->average(sum/n);
	p->standard_deviation(sqrt(ssum/n));

	UnitCell	uc(PIFscale*imageheader->xlength, PIFscale*imageheader->ylength,
					PIFscale*imageheader->zlength, PIFscale*imageheader->alpha,
					PIFscale*imageheader->beta, PIFscale*imageheader->gamma);
	p->unit_cell(uc);
	
	p->space_group(imageheader->ispg);
	
	if ( p->space_group() < 2 ) p->space_group(1);
	
	// Pixel size
	Vector3<double>		sam(1,1,1);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIF: stored pixelsize = " << imageheader->pixelsize_x 
			<< " " << imageheader->pixelsize_y 
			<< " " << imageheader->pixelsize_z << endl;
	if ( imageheader->mx > 0 && imageheader->my > 0 && imageheader->mz > 0)
		sam = {uc.a()/imageheader->mx, uc.b()/imageheader->my, uc.c()/imageheader->mz};
	
	// Trigger for new image header
	if ( imageheader_new->ctfMode < 10 ) {
		if ( imageheader_new->pixelSize > 0 && imageheader_new->pixelSize < 1e9 )
			sam[0] = sam[1] = sam[2] = PIFscale*imageheader_new->pixelSize;
	} else {
		if ( imageheader->pixelsize_x > 0 && imageheader->pixelsize_x < 1e9 )
			sam[0] = PIFscale*imageheader->pixelsize_x;
		if ( imageheader->pixelsize_y > 0 && imageheader->pixelsize_y < 1e9 )
			sam[1] = PIFscale*imageheader->pixelsize_y;
		if ( imageheader->pixelsize_z > 0 && imageheader->pixelsize_z < 1e9 )
			sam[2] = PIFscale*imageheader->pixelsize_z;
	}

	if ( imageheader_new->pixelsize_x > 0 && imageheader_new->pixelsize_x < 1e9 )
			sam[0] = PIFscale*imageheader_new->pixelsize_x;
	if ( imageheader_new->pixelsize_y > 0 && imageheader_new->pixelsize_y < 1e9 )
			sam[1] = PIFscale*imageheader_new->pixelsize_y;
	if ( imageheader_new->pixelsize_z > 0 && imageheader_new->pixelsize_z < 1e9 )
			sam[2] = PIFscale*imageheader_new->pixelsize_z;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIF: assigned pixelsize = " << sam << endl; 
	
	p->sampling(sam);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIF: converted pixelsize = " << sam << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIF: Mode = " << header->mode << endl;

	if ( header->mode == 31 || header->mode == 32 ) {
		read_sf( fimg, p, header->mode, sb );
		cerr << "rwPIF: Structure factor I/O not implemented!" << endl << endl;
		return -4;
	}
	
	long	pagesize = p->channels()*xstore*p->sizeY()*p->sizeZ()*p->data_type_size();
	long	offset(p->data_offset() + imgstart*image_size);
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
				p->unpack_transform(i-imgstart, data, CentHerm);
			}
			delete[] data;
		}
	}
	
	fimg->close();
	delete fimg;
	    
	delete header;
	delete imageheader;
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIF: Done" << endl;
	
	return 0;
}

char*	read_sf(ifstream* fimg, Bimage* p, int mode, int sb)
{
	error_show("Note: Structure factor format within PIF is not supported!", __FILE__, __LINE__);
	
	return NULL;
}

/**
@brief	Writing a PIF image file format.
A 2D and 3D image format used in electron microscopy.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writePIF(Bimage* p)
{
//	double		PIFscale = PIFSCALE;
	
	p->color_to_simple();
	
    switch ( p->data_type() ) {
		case Bit:
    	case SCharacter: p->change_type(UCharacter); break;
    	case UShort: p->change_type(Short); break;
    	default: break;
    }
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePIF: Writing PIF file" << endl;

	// Map the parameters to the file header
	long			header_size = sizeof(PIFhead);
	long			image_header_size = sizeof(PIFimagehead);
	PIFhead*		header = new PIFhead;
	memset(header, 0, sizeof(PIFhead));
	
	// Scale factor to convert integers to floating point
	if ( PIFscale <= 0 || PIFscale > PIF_SCALE ) {
//		if ( p->scale > 0 && p->scale < PIF_SCALE ) PIFscale = p->scale;
//		else PIFscale = PIF_SCALE;
		PIFscale = PIF_SCALE;
		if ( PIFscale > (p->maximum() - p->minimum())*PIF_SCALE ) PIFscale = (p->maximum() - p->minimum())*PIF_SCALE;
	}
	if ( fabs(p->minimum()/PIFscale) > INT_MAX ) PIFscale = fabs(p->minimum()/INT_MAX);
	if ( fabs(p->maximum()/PIFscale) > INT_MAX ) PIFscale = fabs(p->maximum()/INT_MAX);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePIF: PIF scale factor: " << PIFscale << endl;
	double		invPIFscale = 1.0L/PIFscale;
	char		astring[32];

	header->file_id[0] = 8;		// Need to change this!
	header->file_id[4] = 8;		// Need to change this!
	snprintf(astring, 16, "   %13e", PIFscale);
	strncpy(header->realscalefactor, astring, 16);
    snprintf(astring, 16, "%16.9E", invPIFscale);
	strncpy(header->floatScaleStr, astring, 16);
	header->numimages = p->images();
	header->endianness = 1; 	// Byte order
	if ( systype(0) >= LittleIEEE ) header->endianness = 0;
	header->htype = 1;			// Flag to indicate all images of same size
	header->nx = p->sizeX();
	header->ny = p->sizeY();
	header->nz = p->sizeZ();
	switch ( p->data_type() ) {	// Many more types!
		case UCharacter:
			header->mode = 0;
			if ( p->images() > 1 ) header->mode = 6;
			break;
		case Short:
			header->mode = 1;
			if ( p->sizeZ() > 1 ) header->mode = 20;
			if ( p->compound_type() == TComplex ) header->mode = 3;
			break;
		case Integer:
			header->mode = 2;
			if ( p->sizeZ() > 1 ) header->mode = 21;
			if ( p->images() > 1 ) header->mode = 46;
			if ( p->compound_type() == TComplex ) header->mode = 4;
			break;
		case Float:
			header->mode = 9;
			if ( p->compound_type() == TComplex ) header->mode = 10;
			break;
		case Double:
			header->mode = 40;
			if ( p->compound_type() == TComplex ) header->mode = 41;
			break;
		default: header->mode = 0; break;
	}
	snprintf(astring, 32, "Bsoft %s", BVERSION);
	memset(header->genprogram, ' ', 32);
	memcpy(header->genprogram, astring, strlen(astring));
	
	// Map the parameters to the image header
	UnitCell	uc = p->unit_cell();

	PIFimagehead*	imageheader = new PIFimagehead;
	PIFimagehead_new* imageheader_new = (PIFimagehead_new *) imageheader;
	memset(imageheader, 0, sizeof(PIFimagehead));
	
	strncpy(imageheader->timestamp, asctime(p->get_localtime()), 30);
	imageheader->mode = header->mode;
	imageheader->nx = p->sizeX();
	imageheader->ny = p->sizeY();
	imageheader->nz = p->sizeZ();
	imageheader->mx = (int) (uc.a()/p->sampling(0)[0] + 0.5);
	imageheader->my = (int) (uc.b()/p->sampling(0)[1] + 0.5);
	imageheader->mz = (int) (uc.c()/p->sampling(0)[2] + 0.5);
	imageheader_new->pixelsize_x = (int) (p->sampling(0)[0]*invPIFscale + 0.5);
	imageheader_new->pixelsize_y = (int) (p->sampling(0)[1]*invPIFscale + 0.5);
	imageheader_new->pixelsize_z = (int) (p->sampling(0)[2]*invPIFscale + 0.5);
	imageheader_new->pixelSize = imageheader_new->pixelsize_x;
	imageheader_new->pixSizUnits = 1;
	imageheader->mapc = 1;
	imageheader->mapr = 2;
	imageheader->maps = 3;
	imageheader->min = (int) (p->minimum()*invPIFscale + 0.5);
	imageheader->max = (int) (p->maximum()*invPIFscale + 0.5);
	imageheader->mean = (int) (p->average()*invPIFscale + 0.5);
	imageheader->stddev = (int) (p->standard_deviation()*invPIFscale + 0.5);
	imageheader->xlength = (int) (uc.a()*invPIFscale + 0.5);
	imageheader->ylength = (int) (uc.b()*invPIFscale + 0.5);
	imageheader->zlength = (int) (uc.c()*invPIFscale + 0.5);
	imageheader->alpha = (int) (uc.alpha()*180.0/(M_PI*PIFscale) + 0.5);
	imageheader->beta  = (int) (uc.beta()*180.0/(M_PI*PIFscale) + 0.5);
	imageheader->gamma = (int) (uc.gamma()*180.0/(M_PI*PIFscale) + 0.5);
	imageheader->ispg = p->space_group();
	imageheader->aoverb = (int) (1.0*invPIFscale + 0.5);
	imageheader->map_abang = imageheader->gamma;
	imageheader->dela = (int) (1.0*invPIFscale + 0.5);
	imageheader->delb = (int) (1.0*invPIFscale + 0.5);
	imageheader->delc = (int) (1.0*invPIFscale + 0.5);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG writePIF: pixelsize=" << PIFscale*imageheader->pixelsize_x << "," << 
			PIFscale*imageheader->pixelsize_y << "," << PIFscale*imageheader->pixelsize_z << endl;
		cout << "DEBUG writePIF: min=" << imageheader->min << " max=" << imageheader->max << 
			" mean=" << imageheader->mean << " std=" << imageheader->stddev << endl;
		cout << "DEBUG writePIF: Writing labels" << endl;
	}
	
	long 	label_size = p->label().length();
	if ( label_size > 255 ) label_size = 255;
	memcpy(header->label, p->label().c_str(), label_size);
	header->label[label_size] = 0;
	if ( label_size > 79 ) label_size = 79;
	memcpy(imageheader->title, p->label().c_str(), label_size);
	imageheader->title[label_size] = 0;

//	imageheader->nsymbt = p->nsym*80;
	
	p->data_offset(PIFSIZE);
	
	long	xstore = p->sizeX();
	if ( p->compound_type() == TComplex ) xstore = p->sizeX()/2 + 1;
	long	datatypesize = p->channels()*p->data_type_size();
	long	datasize = xstore*p->sizeY()*p->sizeZ()*datatypesize;

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG writePIF: header size=" << header_size << " image header size =" << 
			image_header_size << " new image header size = " << sizeof(PIFimagehead_new) << endl;
		cout << "DEBUG writePIF: offset=" << p->data_offset() << " typesize=" << 
			datatypesize << " datasize =" << datasize << endl;
	}	
	
	long   			n;
	unsigned char*	data = NULL;
 	unsigned char*	aptr = p->data_pointer();

	if ( p->compound_type() == TComplex ) data = new unsigned char[datasize];

	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write((char *)header, header_size);
	
	if ( header->mode == 31 || header->mode == 32 ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG writePIF: Writing structure factors not yet implemented" << endl;

	} else {
		for ( n=0; n<p->images(); n++ ) {
			imageheader->bkgnd = (int) (p->background(n) + 0.5);
			imageheader->min = (int) (p->image[n].minimum()*invPIFscale + 0.5);
			imageheader->max = (int) (p->image[n].maximum()*invPIFscale + 0.5);
			imageheader->mean = (int) (p->image[n].average()*invPIFscale + 0.5);
			imageheader->stddev = (int) (p->image[n].standard_deviation()*invPIFscale + 0.5);
			if ( p->image[n].origin()[0] >= 0 )
				imageheader->nxstart = (int) -(p->image[n].origin()[0] + 0.5);
			else
				imageheader->nxstart = (int) -(p->image[n].origin()[0] - 0.5);
			if ( p->image[n].origin()[1] >= 0 )
				imageheader->nystart = (int) -(p->image[n].origin()[1] + 0.5);
			else
				imageheader->nystart = (int) -(p->image[n].origin()[1] - 0.5);
			if ( p->image[n].origin()[2] >= 0 )
				imageheader->nzstart = (int) -(p->image[n].origin()[2] + 0.5);
			else
				imageheader->nzstart = (int) -(p->image[n].origin()[2] - 0.5);
			imageheader_new->origX = (int) (p->image[n].origin()[0]*invPIFscale + 0.5);
			imageheader_new->origY = (int) (p->image[n].origin()[1]*invPIFscale + 0.5);
			imageheader_new->origZ = (int) (p->image[n].origin()[2]*invPIFscale + 0.5);
			imageheader_new->view_x = (int) (p->image[n].view()[0]*invPIFscale + 0.5);
			imageheader_new->view_y = (int) (p->image[n].view()[1]*invPIFscale + 0.5);
			imageheader_new->view_z = (int) (p->image[n].view()[2]*invPIFscale + 0.5);
			imageheader_new->view_angle = (int) (p->image[n].view().angle()*invPIFscale + 0.5);
			fimg.write((char *)imageheader, image_header_size);
			if ( p->compound_type() < TComplex ) {
				fimg.write((char *)aptr, datasize);
			} else {
				p->pack_transform(n, data, CentHerm);
				fimg.write((char *)data, datasize);
			}
			aptr += datasize;
		}
	}
	
	fimg.close();
	
	if ( data ) delete[] data;
	
	delete header;
	delete imageheader;
		
	return 0;
}

