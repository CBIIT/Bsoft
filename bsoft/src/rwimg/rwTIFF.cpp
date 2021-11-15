/**
@file	rwTIFF.cpp
@brief	Functions for reading and writing TIFF files
@author Bernard Heymann
@date	Created: 19990509
@date 	Modified: 20210615
**/

#include "rwTIFF.h"
#include "utilities.h"

#ifdef HAVE_TIFF

#include "tiff.h"
#include "tiffio.h"
#include "tiffiop.h"

#ifdef HAVE_XML
#include "rwxml.h"
#endif

// Tags for EER compression
#define TIFF_COMPRESSION_EER_V0 			65000
#define TIFF_COMPRESSION_EER_V1 			65001
#define TIFFTAG_EER_ACQUISITION_METADATA	65001

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

extern char* month[];

// Initialize EER decodec
int		TIFFInitEER(TIFF* tif, int scheme);

// Custom tags
static const TIFFFieldInfo xtiffFieldInfo[] = {
	{ TIFF_COMPRESSION_EER_V0, TIFF_VARIABLE, TIFF_VARIABLE, TIFF_UNDEFINED, FIELD_CUSTOM, true, true, const_cast<char*>("EER compression v0") },
	{ TIFF_COMPRESSION_EER_V1, TIFF_VARIABLE, TIFF_VARIABLE, TIFF_UNDEFINED, FIELD_CUSTOM, true, true, const_cast<char*>("EER compression v1") },
	{ TIFFTAG_EER_ACQUISITION_METADATA, TIFF_VARIABLE, TIFF_VARIABLE, TIFF_UNDEFINED, FIELD_CUSTOM, true, true, const_cast<char*>("Acquisition metadata") }
};

static void registerTIFFTags(TIFF *tif)
{
//    int error = TIFFMergeFieldInfo(tif, xtiffFieldInfo, 3);
    TIFFMergeFieldInfo(tif, xtiffFieldInfo, 3);
}

/**
@brief	Reading a TIFF image file format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@param 	img_select	image selection in multi-image file (-1 = all images).
@return	int			error code (<0 means failure).
A 2D and 3D image format commonly used for graphics.
	TIFF 6.0 library of Sam Lefler.
	Data types: 		bit, byte, signed char, unsigned short,
							signed short, int, float
	Color models:		gray scale, RGB, RGBA, CMYK
	Indexed color model
**/
int 	readTIFF(Bimage* p, int readdata, int img_select)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readTIFF: Reading image: " << p->file_name() << endl;

	TIFFSetTagExtender(registerTIFFTags);
	TIFFRegisterCODEC(TIFF_COMPRESSION_EER_V1, "EER compression", TIFFInitEER);

    TIFF*       	fimg;
    if ( ( fimg = TIFFOpen(p->file_name().c_str(), "r") ) == NULL ) return -1;

	int 			flags(0);
	if ( verbose & VERB_DEBUG ) {
//		cout << "Magic number:" << fimg->tif_header.tiff_magic << endl;
		flags = TIFFPRINT_STRIPS | TIFFPRINT_CURVES | TIFFPRINT_COLORMAP;
		TIFFPrintDirectory(fimg, stdout, flags);
	}
	
	// Set defaults
	int				notdone(1), nimg(0);
	long   			i, n, iz, iy, j, jz, jy, k;
	unsigned int	x(1), y(1), z(1), px, py, pz;
	unsigned short	compression(COMPRESSION_NONE);
	unsigned short	samplesperpixel(1), bitspersample(0), lutsize(0);
	unsigned short	sampleformat(SAMPLEFORMAT_UINT), photometric(0);
	unsigned short	planarconfig(PLANARCONFIG_CONTIG);
	unsigned int 	orientation(ORIENTATION_BOTLEFT);
	float			xresolution(1), yresolution(1);
	unsigned short* red = NULL;
	unsigned short*	green = NULL;
	unsigned short*	blue = NULL;

	TIFFGetField(fimg, TIFFTAG_COMPRESSION, &compression);
	(*p)["compression"] = compression;

	if ( compression == TIFF_COMPRESSION_EER_V1 ) {
		TIFFInitEER(fimg, compression);
		if ( verbose & VERB_PROCESS ) {
			char*		eer_metadata = NULL;
			int			eer_md_count = 0;
			TIFFGetField(fimg, TIFFTAG_EER_ACQUISITION_METADATA, &eer_md_count, &eer_metadata);
			JSvalue		jmd = json_from_xml(eer_metadata);
			cout << "EER metadata:" << endl << jmd << endl << endl;
		}
	}
	
	// Get the number of images in the file
	if ( img_select > -1 ) {
		while ( TIFFReadDirectory(fimg) ) nimg++;	
		if ( img_select > nimg ) img_select = nimg;
		p->images(1);
		TIFFSetDirectory(fimg, img_select);
		TIFFGetField(fimg, TIFFTAG_IMAGEWIDTH,  &x);
		TIFFGetField(fimg, TIFFTAG_IMAGELENGTH, &y);
		TIFFGetField(fimg, TIFFTAG_IMAGEDEPTH,  &z);
		p->size(x, y, z);
	} else {
		do {
			TIFFGetField(fimg, TIFFTAG_IMAGEWIDTH,  &x);	// Test if image is the same size as the first image
			TIFFGetField(fimg, TIFFTAG_IMAGELENGTH, &y);
			TIFFGetField(fimg, TIFFTAG_IMAGEDEPTH,  &z);
			if ( nimg == 0 ) {
				p->size(x, y, z);
				nimg++;
			} else {
				if ( x != p->sizeX() || y != p->sizeY() || z != p->sizeZ() ) {
					cerr << "Warning: TIFF image " << nimg << " is of a different size: " << x << "," << y << "," << z << endl;
					cerr << " ---> Only " << nimg << " images read!" << endl;
					notdone = 0;
				} else nimg++;
			}
		} while ( notdone && TIFFReadDirectory(fimg) );
		p->images(nimg);
		TIFFSetDirectory(fimg, 0);
	}
	
	tm*			t = p->get_localtime();
	char*		timestring = NULL, mon[4];
	TIFFGetField(fimg, TIFFTAG_DATETIME, &timestring);
	if ( timestring ) {
		if ( verbose & VERB_DEBUG )
			cout << "TIFF timestring: " << timestring << endl;
		if ( isdigit(timestring[0]) ) {
			sscanf(timestring, "%4d:%2d:%2d %2d:%2d:%2d", 
				&t->tm_year, &t->tm_mon, &t->tm_mday, 
				&t->tm_hour, &t->tm_min, &t->tm_sec);
			t->tm_mon -= 1;
		} else {
			sscanf(timestring, "%3s %3s %2d %2d:%2d:%2d %4d", 
				mon, mon, &t->tm_mday, &t->tm_hour, &t->tm_min, &t->tm_sec, &t->tm_year);
			for ( i=0; strncmp(mon, month[i], 3) != 0; i++ ) ;
			if ( i < 12 ) t->tm_mon = i;
		}
		t->tm_year -= 1900;
		p->set_time(t);
	}
		
	TIFFGetField(fimg, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
	
	// Allocate the sub-image array structure
//	cout << "***** TIFF images read = " << p->n << endl;
	
	p->channels(samplesperpixel);
	if ( p->sizeX() < 1 ) p->sizeX(1);
	if ( p->sizeY() < 1 ) p->sizeY(1);
	if ( p->sizeZ() < 1 ) p->sizeZ(1);
	if ( p->channels() < 1 ) p->channels(1);
	p->page_size(p->size());
	if ( TIFFIsTiled(fimg) ) {
		TIFFGetField(fimg, TIFFTAG_TILEWIDTH,  &x);
		TIFFGetField(fimg, TIFFTAG_TILELENGTH, &y);
		TIFFGetField(fimg, TIFFTAG_TILEDEPTH,  &z);
		p->page_size(x, y, z);
	}

	TIFFGetField(fimg, TIFFTAG_XRESOLUTION, &xresolution);
	TIFFGetField(fimg, TIFFTAG_YRESOLUTION, &yresolution);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readTIFF: xresolution=" << xresolution << " yresolution=" << yresolution << endl;
	p->sampling(xresolution, yresolution, 1);
	if ( p->sizeZ() > 1 && fabs(p->sampling(0)[0] - p->sampling(0)[1]) < 0.001 )
		p->sampling(xresolution, yresolution, xresolution);
	
	TIFFGetField(fimg, TIFFTAG_BITSPERSAMPLE, &bitspersample);
	TIFFGetField(fimg, TIFFTAG_PHOTOMETRIC, &photometric);
	TIFFGetField(fimg, TIFFTAG_SAMPLEFORMAT, &sampleformat);
	TIFFGetField(fimg, TIFFTAG_ORIENTATION, &orientation);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readTIFF: bitspersample=" << bitspersample << " tiff datatype=" << sampleformat << endl;
	
	switch ( photometric ) {
		case PHOTOMETRIC_MINISWHITE:
		case PHOTOMETRIC_MINISBLACK:
			p->compound_type(TSimple); break;
		case PHOTOMETRIC_RGB:
		case PHOTOMETRIC_YCBCR:
			p->compound_type(TRGB);
			if ( p->channels() == 4 ) p->compound_type(TRGBA);
			break;
		case PHOTOMETRIC_SEPARATED:
			p->compound_type(TCMYK);
			break;
		case PHOTOMETRIC_PALETTE:
			p->compound_type(TRGB);
			lutsize = 1;
			for ( i=0; i<bitspersample; i++ ) lutsize *= 2;
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG readTIFF: lutsize=" << lutsize << endl;
			break;
		default:
			p->compound_type(TSimple);
	}
	
	switch ( bitspersample ) {
		case 1:
			p->data_type(Bit);
			break;
		case 8:
			p->data_type(UCharacter);
			if ( sampleformat == SAMPLEFORMAT_INT )
				p->data_type(SCharacter);
			break;
		case 16:
			p->data_type(UShort);
			if ( sampleformat == SAMPLEFORMAT_INT )
				p->data_type(Short);
			break;
		case 32:
			p->data_type(UInteger);
			if ( sampleformat == SAMPLEFORMAT_INT )
				p->data_type(Integer);
			else if ( sampleformat == SAMPLEFORMAT_IEEEFP )
				p->data_type(Float);
			break;
		case 64:
			p->data_type(ULong);
			if ( sampleformat == SAMPLEFORMAT_INT )
				p->data_type(Long);
			else if ( sampleformat == SAMPLEFORMAT_IEEEFP )
				p->data_type(Double);
			break;
		default:
			p->data_type(UCharacter);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readTIFF: Datatype=" << p->data_type() << " Compoundtype=" << p->compound_type() << endl;
	
	unsigned short			min(0), max(0);
	double			dmin(0), dmax(0);
	if ( p->data_type() <= Short ) {
		TIFFGetField(fimg, TIFFTAG_MINSAMPLEVALUE, &min);
		TIFFGetField(fimg, TIFFTAG_MAXSAMPLEVALUE, &max);
		p->minimum(min);
		p->maximum(max);
	} else {
		TIFFGetField(fimg, TIFFTAG_SMINSAMPLEVALUE, &dmin);
		TIFFGetField(fimg, TIFFTAG_SMAXSAMPLEVALUE, &dmax);
		p->minimum(dmin);
		p->maximum(dmax);
	}

	char*			label = NULL;
	TIFFGetField(fimg, TIFFTAG_IMAGEDESCRIPTION, &label);
	if ( label ) p->label(label);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readTIFF: Minimum and maximum: " << p->minimum() << " " << p->maximum() << endl;
	
	int 			tilesize;
	long   			datatypesize = p->data_type_size();
	long			elementsize = p->channels()*datatypesize;
	long 			linesize = p->sizeX()*elementsize;
	unsigned char* 	data = NULL;
	char 			*page;
	
	if ( p->data_type() == Bit ) linesize = p->page_size()[0]/8;

	TIFFGetField(fimg, TIFFTAG_PLANARCONFIG, &planarconfig);
	if ( samplesperpixel > 1 && planarconfig == PLANARCONFIG_SEPARATE ) {
		cerr << "Error: Separate planes not supported!" << endl;
//		bexit(-1);
	}
/*
	TIFFDirectory *td = &fimg->tif_dir;
	cout << "fimg pointer: " << fimg << endl;
	cout << "raw data pointer: " << fimg->tif_rawdata << endl;
	cout << "Strip offset entry: " << td->td_stripoffset_p << endl;
	p->data_offset((unsigned char *)td->td_stripoffset_p-(unsigned char *)fimg);
*/
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readTIFF: readdata=" << readdata << " datasize=" << 
			p->data_size() << " linesize=" << linesize << endl;
		
	if ( !readdata ) {
		TIFFClose(fimg);
		if ( lutsize ) {
			p->data_type(UShort);
			p->channels(3);
		}
		return 0;
	}

	data = p->data_alloc();

	unsigned short*	rgbdata = NULL;
	unsigned short	idx;
	long   			datasize = p->images()*p->sizeX()*p->sizeY()*p->sizeZ();
	if ( lutsize ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readTIFF: lutsize=" << lutsize << endl;
		rgbdata = new unsigned short[3*datasize];
	}
	
	n = nimg = 0;
	do {
		if ( lutsize ) {
			TIFFGetField(fimg, TIFFTAG_COLORMAP, &red, &green, &blue);
			if ( verbose & VERB_DEBUG ) {
				cout << "DEBUG readTIFF: lut:" << endl;
				for ( i=0; i<lutsize; i++ )
					cout << i << tab << red[i] << tab << green[i] << tab << blue[i] << endl;
			}
		}
			TIFFGetField(fimg, TIFFTAG_PLANARCONFIG, &planarconfig);
			if ( samplesperpixel > 1 && planarconfig == PLANARCONFIG_SEPARATE ) {
				cerr << "Error: Separate planes not supported!" << endl;
				bexit(-1);
			}
			if ( TIFFIsTiled(fimg) ) {
				TIFFGetField(fimg, TIFFTAG_TILEWIDTH,  &x);
				TIFFGetField(fimg, TIFFTAG_TILELENGTH, &y);
				TIFFGetField(fimg, TIFFTAG_TILEDEPTH,  &z);
				p->page_size(x, y, z);
				tilesize = TIFFTileSize(fimg);
			
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG readTIFF: tile size: " << x << " x " << y 
						<< " x " << z << " = " << tilesize << endl;

				page = new char[tilesize];
				for ( pz=0; pz<p->sizeZ(); pz+=p->page_size()[2] ) {
					for ( py=0; py<p->sizeY(); py+=p->page_size()[1] ) {
						for ( px=0; px<p->sizeX(); px+=p->page_size()[0] ) {
							TIFFReadTile(fimg, (char *) page, px, py, pz, 0);
							for ( z=0; z<p->page_size()[2] && z<p->sizeZ()-pz; z++ ) {
								iz = (nimg*p->sizeZ() + z + pz)*p->sizeY();
								jz = z*p->page_size()[1];
								for ( y=0; y<p->page_size()[1] && y<p->sizeY()-py; y++ ) {
									iy = (iz + y + py)*p->sizeX();
									jy = (jz + y)*p->page_size()[0];
									if ( p->data_type() == Bit ) {
										iy /= 8;
										jy /= 8;
										memcpy(data+iy, page+jy, linesize);
									} else {
										for ( x=0; x<p->page_size()[0] && x<p->sizeX()-px; x++ ) {
											i = (iy + x + px)*elementsize;
											j = (jy + x)*elementsize;
											memcpy(data+i, page+j, elementsize);
										}
									}
								}
							}
						}
					}
				}
				delete[] page;
			
			} else if ( compression == TIFF_COMPRESSION_EER_V1 ) {
				i = nimg*p->image_size();
				fimg->tif_scanlinesize = p->sizeX()*p->sizeY();
				TIFFReadScanline(fimg, data+i, 0, 0);
				if ( verbose )
					cout << "Electrons counted: " << fimg->tif_rawcc << endl;
			} else {

				linesize = TIFFScanlineSize(fimg);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG readTIFF: scanline size=" << linesize << endl;
				
/*
				auto nrOfStrips = TIFFNumberOfStrips(fimg);
				for (uint32_t i = 0; i < nrOfStrips; i++) {
					auto stripSize = TIFFRawStripSize(fimg, i);
					cout << nimg << tab << nrOfStrips << tab << i << tab << stripSize << tab << fimg->tif_scanlinesize << endl;
				}
*/
				i = nimg*p->sizeZ()*p->sizeY()*linesize;
				for ( z=0; z<p->sizeZ(); z++ ) {
					for ( y=0; y<p->sizeY(); y++, i+=linesize ) {
						// Note: Scanline rows are ambiguous for 3D images
						TIFFReadScanline(fimg, data+i, y, 0);
						if ( lutsize ) {
							if ( bitspersample == 8 ) {
								for ( x=0, j=i, k=3*j; x<p->sizeX(); x++, j++ ) {
									idx = data[j];
									rgbdata[k++] = red[idx];
									rgbdata[k++] = green[idx];
									rgbdata[k++] = blue[idx];
								}
							} else {
								unsigned short*	usdata = (unsigned short *) data;
								for ( x=0, j=i/2, k=3*j; x<p->sizeX(); x++, j++ ) {
									idx = usdata[j];
									rgbdata[k++] = red[idx];
									rgbdata[k++] = green[idx];
									rgbdata[k++] = blue[idx];
								}
							}
						}
					}
				}
			}
			p->image[nimg].view(0,0,1,0);
			nimg++;
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG readTIFF: image " << n << " read" << endl;
//		}
		n++;
	} while ( n < p->images() && TIFFReadDirectory(fimg) );
		
	TIFFClose(fimg);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readTIFF: images read: " << p->images() << endl << endl;

	if ( lutsize ) {
		p->data_type(UShort);
		p->channels(3);
		p->data_assign((unsigned char *) rgbdata);
	}

	// Adherance to conventions:
	//		Density is positive-valued displayed as white on black background
	//		Display orientation is with the origin at the bottom left
	Bstring		order("xyz");
	if ( p->data_pointer() ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readTIFF: Tiff photometric order: " << photometric << endl;
		if ( photometric == PHOTOMETRIC_MINISWHITE ) p->invert();
		switch ( orientation ) {
			case ORIENTATION_TOPLEFT:	order = "x-yz"; break;
			case ORIENTATION_TOPRIGHT:	order = "-x-yz"; break;
			case ORIENTATION_BOTRIGHT:	order = "-xyz"; break;
			case ORIENTATION_LEFTTOP:	order = "-yxz"; break;
			case ORIENTATION_RIGHTTOP:	order = "-y-xz"; break;
			case ORIENTATION_RIGHTBOT:	order = "y-xz"; break;
			case ORIENTATION_LEFTBOT:	order = "yxz"; break;
			default:	order = "xyz";
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readTIFF: Tiff image order: " << order << endl;
		if ( orientation != ORIENTATION_BOTLEFT ) p->reslice(order);
	}	

	return 0;
}

/**
@brief	Writing a TIFF image file format.
@param	*p			the image structure.
@param	flags		flags to set output properties.
@return	int			error code (<0 means failure).
A 2D and 3D image format commonly used for graphics.
	TIFF 6.0 library of Sam Lefler.
	flags:
	1		write tiled TIFF images (default scanline)
	2+		write compressed images.
**/
int 	writeTIFF(Bimage* p, int flags)
{
    TIFF			*fimg;
    if ( ( fimg = TIFFOpen(p->file_name().c_str(), "w") ) == NULL ) return -1;
	
	// Map the parameters
	int 			istiled = flags & 1;
//	unsigned short	compression = flags & 62;
//	unsigned short	compression(COMPRESSION_LZW);
	unsigned short	compression(flags & COMPRESSION_LZW);
	long   			i, j, y, z, n;
	unsigned short	min = (unsigned short) p->minimum();
	unsigned short	max = (unsigned short) p->maximum();
	unsigned short	resolution_unit(1);
	float			xresolution = p->sampling(0)[0];
	float			yresolution = p->sampling(0)[1];
	unsigned int  	rowsperstrip = (unsigned int) -1;
	
	long 			datatypesize = p->data_type_size();
	long 			datasize = p->sizeX()*p->sizeY()*p->sizeZ()*p->channels()*datatypesize;
	long 			linesize = p->sizeX()*p->channels()*datatypesize;
	long 			unitsize = 64;
	long 			pagesize, tilerowsize(0);
	char* 			page;
//	unsigned char*	scanline;
	
	if ( p->data_type() == Bit ) {
		linesize = p->page_size()[0]/8;
		datasize = linesize*p->sizeY()*p->sizeZ();
	}
	
	p->page_size((int) (unitsize*floor(p->sizeX()*1.0/unitsize + 0.99999)),
				(int) (unitsize*floor(p->sizeY()*1.0/unitsize + 0.99999)), 1);
	
	long tilex = (int) (unitsize*floor(p->sizeX()*1.0/unitsize + 0.99999));
	long tiley = (int) (unitsize*floor(p->sizeY()*1.0/unitsize + 0.99999));
	long tilez(1);
	
	if ( compression ) (*p)["compression"] = compression;
	
	tm*		t = p->get_localtime();
	char	timestring[20];
	snprintf(timestring, 20, "%4d:%02d:%02d %02d:%02d:%02d",
			t->tm_year+1900, t->tm_mon+1, t->tm_mday, 
			t->tm_hour, t->tm_min, t->tm_sec);
	
	for ( n=0; n<p->images(); n++ ) {
		TIFFSetDirectory(fimg, n);
		TIFFSetField(fimg, TIFFTAG_IMAGEWIDTH, p->sizeX());
		TIFFSetField(fimg, TIFFTAG_IMAGELENGTH, p->sizeY());
		TIFFSetField(fimg, TIFFTAG_IMAGEDEPTH, p->sizeZ());
		TIFFSetField(fimg, TIFFTAG_SAMPLESPERPIXEL, p->channels());
		if ( istiled ) {
			TIFFSetField(fimg, TIFFTAG_TILEWIDTH, tilex);
			TIFFSetField(fimg, TIFFTAG_TILELENGTH, tiley);
			TIFFSetField(fimg, TIFFTAG_TILEDEPTH, tilez);
		} else
			TIFFSetField(fimg, TIFFTAG_ROWSPERSTRIP,
	    		TIFFDefaultStripSize(fimg, rowsperstrip));
		if ( compression )
			TIFFSetField(fimg, TIFFTAG_COMPRESSION, compression);
	
		if ( verbose & VERB_DEBUG ) {
			cout << "DEBUG writeTIFF: Image " << n << ": size " << p->sizeX() << " " 
				<< p->sizeY() << " " << p->sizeZ() << " " << p->channels() << endl;
			cout << "DEBUG writeTIFF: page dimensions: " <<
					p->page_size()[0] << " " << p->page_size()[1] << " " << p->page_size()[2] << endl;
		}
	
		TIFFSetField(fimg, TIFFTAG_RESOLUTIONUNIT, resolution_unit);
		TIFFSetField(fimg, TIFFTAG_XRESOLUTION, xresolution);
		TIFFSetField(fimg, TIFFTAG_YRESOLUTION, yresolution);
	
		TIFFSetField(fimg, TIFFTAG_IMAGEDESCRIPTION, p->label().c_str());
		TIFFSetField(fimg, TIFFTAG_DATETIME, timestring);
	
		if ( p->data_type() <= Short ) {
			TIFFSetField(fimg, TIFFTAG_MINSAMPLEVALUE, min);
			TIFFSetField(fimg, TIFFTAG_MAXSAMPLEVALUE, max);
		} else {
			TIFFSetField(fimg, TIFFTAG_SMINSAMPLEVALUE, p->minimum());
			TIFFSetField(fimg, TIFFTAG_SMAXSAMPLEVALUE, p->maximum());
		}
		
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG writeTIFF: min & max: " << min << " " << max << endl;
	
		switch ( p->data_type() ) {
			case Bit:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 1);
				break;
			case UCharacter:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 8);
				break;
			case SCharacter:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 8);
				break;
			case UShort:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 16);
				break;
			case Short:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 16);
				break;
			case UInteger:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 32);
				break;
			case Integer:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 32);
				break;
			case ULong:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 64);
				break;
			case Long:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 64);
				break;
			case Float:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 32);
				break;
			case Double:
				TIFFSetField(fimg, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 64);
				break;
			default:
				TIFFSetField(fimg, TIFFTAG_BITSPERSAMPLE, 8);
		}
		
		int			extraTypes[1] = {EXTRASAMPLE_ASSOCALPHA};
		switch ( p->compound_type() ) {
			case TSimple:
				TIFFSetField(fimg, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
				break;
			case TVector3:
			case TRGB:
				TIFFSetField(fimg, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
				break;
			case TRGBA:
				TIFFSetField(fimg, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
				TIFFSetField(fimg, TIFFTAG_EXTRASAMPLES, 1, extraTypes);
				break;
			case TCMYK:
				TIFFSetField(fimg, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_SEPARATED);
				TIFFSetField(fimg, TIFFTAG_INKSET, INKSET_CMYK);
				break;
			default:
				TIFFSetField(fimg, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		}

		TIFFSetField(fimg, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(fimg, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
	
		if ( verbose & VERB_DEBUG ) {
			cout << "DEBUG writeTIFF: configuration set to planar" << endl;
			cout << "DEBUG writeTIFF: Typesize=" << datatypesize << " Datasize=" << datasize << endl;
			if ( istiled )
				cout << "DEBUG writeTIFF: Tilesize=" << TIFFTileSize(fimg) << endl;
			else
				cout << "DEBUG writeTIFF: Scanlinesize=" << TIFFScanlineSize(fimg) << endl;
			cout << "DEBUG writeTIFF: Pages: ";
		}
	
		if ( istiled ) {
			pagesize = TIFFTileSize(fimg);
			page = new char[pagesize];
			tilerowsize = TIFFTileRowSize(fimg);
			for ( z=0; z<p->sizeZ(); z++ ) {
				memset(page, 0, pagesize);
				for ( y=0; y<p->sizeY(); y++ ) {
					i = y*tilerowsize;
					j = ((n*p->sizeZ() + z)*p->sizeY() + y)*linesize;
					memcpy(page+i, p->data_pointer()+j, linesize);
				}
				if ( verbose & VERB_DEBUG ) cout << " " << z;
				TIFFWriteTile(fimg, page, 0, 0, z, 0);
			}
			delete[] page;
		} else {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG writeTIFF: linesize=" << linesize << endl;
			for ( z=0; z<p->sizeZ(); z++ ) {
				for ( y=0; y<p->sizeY(); y++ ) {
					j = ((n*p->sizeZ() + z)*p->sizeY() + y)*linesize;
					TIFFWriteScanline(fimg, p->data_pointer()+j, y, 0);
					// Note: Scanline rows are ambiguous for 3D images
				}
				if ( verbose & VERB_DEBUG ) cout << " " << z;
			}
		}
		if ( verbose & VERB_DEBUG ) cout << endl;
		
		TIFFWriteDirectory(fimg);
	} // End of multi-image loop
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeTIFF: data written" << endl;
	
	TIFFClose(fimg);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeTIFF: Memory freed" << endl;

	return 0;
}


#include "tiffiop.h"

static int
EERFixupTags(TIFF* tif)
{
	(void) tif;
	return (1);
}

/*
 * Encode a hunk of pixels.
 */
static int
EEREncode(TIFF* tif, uint8_t* pp, tmsize_t cc, uint16_t s)
{
	(void) pp; (void) cc; (void) s;
	return (1);
}

/* Fixes bitshift of all-1 code */
int			shift_right(unsigned long& bits, long n)
{
	if ( bits < ULLONG_MAX ) {
		bits >>= n;
		return 0;
	}
	
	bits /= 2;
	bits >>= (n-1);

	return 0;
}

/*
 * Decode a hunk of pixels.
	cc is the intended row length after decompression
 */
static int
EERDecode(TIFF* tif, uint8_t* buf, tmsize_t cc, uint16_t s)
{
	static const char module[] = "EERDecode";
	(void) s;
//	cout << cc << tab << tif->tif_rawcc << tab << (8*tif->tif_rawcc)/11 << tab << tif->tif_rawcc%8 << endl;
//	cout << endl << endl << endl;
/*
	ofstream		fraw("test2.raw");
	if ( fraw.fail() )  return -1;
	
	fraw.write((char *)tif->tif_rawcp, tif->tif_rawcc);
	
	fraw.close();
	
	bexit(0);
*/
//	uint8_t* b = buf;

	long			pxbits(7);					// Value size for a pixel
	long			subpx(2);					// Pixel division depth
	long			codebits(pxbits+2*subpx);	// Full size of a code element
	long			maxcode((1<<codebits)-1);	// Maximum code element value covering all bits
	long			maxval((1<<pxbits)-1);		// Maximum number of pixels in a code element
	long			n(0), xbits, rbits(64), npx(0), ne(0);
	unsigned long	code(0);
	long			rawcc(tif->tif_rawcc);
	unsigned long*	rawcp = (unsigned long *)tif->tif_rawcp;
	unsigned long	bits = *rawcp;			// Read the first word
	rawcc -= sizeof(long);
	
	while ( n<cc && rawcc > 0 ) {
		code = bits & maxcode;
//		bits >>= codebits;
		shift_right(bits, codebits);
		rbits -= codebits;
		if ( rbits <= 0 ) {
			bits = *(++rawcp);	// Read the next word
			rawcc -= sizeof(long);
			xbits = -rbits;
			if ( xbits > 0 && xbits <= codebits ) {	// Add the missing part
				code |= ((bits & ((1<<xbits)-1)) << (codebits-xbits));
//				bits >>= xbits;
				shift_right(bits, xbits);
			}
			rbits += 64;
		}
		npx = code & maxval;
		n += npx;
		buf += npx;
		if ( n >= cc ) break;
		if ( npx < maxval ) {		// Count an electron
			*buf = 1;
			buf++;
			n++;
			ne++;
		} else {
//			cout << code << tab << rbits << endl;
//			bitset<64> bs(bits);
//			cout << bits << tab << bs << endl;
			rbits += 2*subpx;
			if ( rbits > 64 ) {
				rbits -= 64;
				rawcp--;
				rawcc += sizeof(long);
			}
			bits = *rawcp;
//			bits >>= (codebits - rbits);
			shift_right(bits, 64 - rbits);
//			bitset<64> bs(bits);
//			bs = bits;
//			cout << bits << tab << bs << endl;
//			cout << (bits & maxval) << endl;
		}
//		cout << npx << tab << n << tab << rbits << tab << ne << endl;
//		cout << ne << tab << npx << endl;
	}
//	cout << n << tab << rawcc << tab << ne << endl;
/*
	for ( rawcc=0; rawcc<4096; ++rawcc )
		cout << " " << int(*(b++));
	cout << endl;
	bexit(0);
*/
	tif->tif_rawcc = ne;
	
	if ( n > cc ) {
		TIFFErrorExt(tif->tif_clientdata, module,
"Not the correct length for scanline %d, expected %ld bytes, decompressed to %ld bytes",
		             tif->tif_row,
		             cc,
		             n);
//		bexit(-1);
		return (0);
	}

	return (1);
}

/*
 * Seek forwards nrows in the current strip.
 */
static int	EERSeek(TIFF* tif, uint32_t nrows)
{
//	tif->tif_rawcp += nrows * tif->tif_scanlinesize;
//	tif->tif_rawcc -= nrows * tif->tif_scanlinesize;
	return (1);
}

/*
 * Initialize EER compression.
 */
int		TIFFInitEER(TIFF* tif, int scheme)
{
	(void) scheme;
//	assert(scheme == TIFF_COMPRESSION_EER_V1);
//	cout << "Initializing EER decodec" << endl;
	tif->tif_fixuptags = EERFixupTags;
	tif->tif_decoderow = EERDecode;
	tif->tif_decodestrip = EERDecode;
	tif->tif_decodetile = EERDecode;
	tif->tif_encoderow = EEREncode;
	tif->tif_encodestrip = EEREncode;
	tif->tif_encodetile = EEREncode;
	tif->tif_seek = EERSeek;
	return (1);
}

#else

int 	readTIFF(Bimage* p, int readdata, int img_select)
{
	cerr << "Error: TIFF files are not supported!" << endl;
	
	return -1;
}

int 	writeTIFF(Bimage* p, int istiled)
{
	cerr << "Error: TIFF files are not supported!" << endl;
	
	return -1;
}


#endif
