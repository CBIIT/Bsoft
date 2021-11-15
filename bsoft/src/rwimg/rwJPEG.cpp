/**
@file	rwJPEG.cpp
@brief	Functions for reading and writing JPEG files
@author Bernard Heymann
@date	Created: 20030510
@date 	Modified: 20210402
**/

#include "rwJPEG.h"
#include "utilities.h"

#ifdef HAVE_JPEG

#include "jpeglib.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a JPEG image file format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int			error code (<0 means failure).
A 2D image format commonly used for graphics.
	JPEG version 6b library.
	The data is read to invert the y-axis to adhere to the top-left
	convention in the JPEG specification.
	Data types: 				byte, RGB
**/
int			readJPEG(Bimage* p, int readdata)
{
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;

	FILE*			infile;			/* source file */
	JSAMPARRAY		buffer;			/* Output row buffer */
	int				row_stride;		/* physical row width in output buffer */

	if ((infile = fopen(p->file_name().c_str(), "rb")) == NULL) return -1;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);

	jpeg_stdio_src(&cinfo, infile);

	jpeg_read_header(&cinfo, TRUE);
  
	p->images(1);
	p->size(cinfo.image_width, cinfo.image_height, 1);
	p->channels(cinfo.num_components);
	if ( cinfo.out_color_space == JCS_RGB ) {
		p->compound_type(TRGB);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readJPEG: Color model: RGB" << endl;
	} else
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readJPEG: Color model: Gray" << endl;
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readJPEG: Image size: " << p->size() << endl;
	
	
	if ( readdata ) {
		unsigned char*	data = p->data_alloc();

		jpeg_start_decompress(&cinfo);

		row_stride = cinfo.image_width * cinfo.num_components;
  
		buffer = (*cinfo.mem->alloc_sarray)
			((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

		while (cinfo.output_scanline < cinfo.image_height) {
			(void) jpeg_read_scanlines(&cinfo, buffer, 1);
			memcpy(&data[(p->sizeY()-cinfo.output_scanline) * row_stride], buffer[0], row_stride);
		}

		jpeg_finish_decompress(&cinfo);
	}
	
	jpeg_destroy_decompress(&cinfo);

	fclose(infile);
	
	return 0;
}

/**
@brief	Writing a JPEG image file format.
@param	*p			the image structure.
@param	quality		level of compression (25-100).
@return	int			error code (<0 means failure).
A 2D image format commonly used for graphics.
	JPEG version 6b library.
	The data is written to invert the y-axis to adhere to the top-left
	convention in the JPEG specification.
**/
int		writeJPEG(Bimage* p, int quality)
{
	if ( quality < 25 ) quality = 100;
	
	p->change_type(UCharacter);
	
	if ( p->sizeZ() > 1 )
		cerr << "Warning: Only the first slice of the image will be written!" << endl;
	
	if ( p->images() > 1 )
		cerr << "Warning: Only the first image will be written!" << endl;
	
	int			err(0);
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;

	FILE*		outfile;				/* target file */
	JSAMPROW	row_pointer[1];			/* pointer to JSAMPLE row[s] */
	int			row_stride;				/* physical row width in image buffer */

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	if ((outfile = fopen(p->file_name().c_str(), "wb")) == NULL) return -1;

	jpeg_stdio_dest(&cinfo, outfile);

	cinfo.image_width = p->sizeX(); 			/* image width and height, in pixels */
	cinfo.image_height = p->sizeY();
	cinfo.input_components = p->channels();		/* # of color components per pixel */
	cinfo.in_color_space = JCS_GRAYSCALE;
	if ( p->compound_type() == TRGB || p->compound_type() == TVector3 )
		cinfo.in_color_space = JCS_RGB;
		
	unsigned char*	data = p->data_pointer();

	jpeg_set_defaults(&cinfo);

	jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

	jpeg_start_compress(&cinfo, TRUE);

	row_stride = p->sizeX()*p->channels();	/* JSAMPLEs per row in image_buffer */

	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer[0] = (JSAMPLE *) &data[(p->sizeY()-cinfo.next_scanline-1) * row_stride];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	jpeg_finish_compress(&cinfo);

	fclose(outfile);

	jpeg_destroy_compress(&cinfo);
	
	return err;
}

#else

int 	readJPEG(Bimage* p, int readdata)
{
	cerr << "Error: JPEG files are not supported!" << endl;
	
	return -1;
}

int 	writeJPEG(Bimage* p)
{
	cerr << "Error: JPEG files are not supported!" << endl;
	
	return -1;
}

#endif
