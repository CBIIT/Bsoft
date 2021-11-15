/**
@file	rwPNG.cpp
@brief	Functions for reading and writing PNG files
@author Bernard Heymann
@date	Created: 20041223
@date 	Modified: 20170317
**/

#include "rwPNG.h"
#include "utilities.h"

#ifdef HAVE_PNG

#include "png.h"
#include "zlib.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading an PNG or text image format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int			error code (<0 means failure).
A 2D image format commonly used on the Web.
	libpng version 1.2.8 - December 3, 2004
	The PNG image format specifies 5 bit depths, with only three spported here:
		1		bitmap
		2		(not supported)
		4		(not supported)
		8		unsigned char
		16		unsigned short
	Color models:
		Gray scale, RGB and RGBA
		(Gray scale with alpha not supported)
**/
int 	readPNG(Bimage* p, int readdata)
{
	FILE*				infile;			/* source file */
    unsigned char		sig[8];

	if ((infile = fopen(p->file_name().c_str(), "rb")) == NULL) return -1;

	// Signature check
    size_t				nr = fread(sig, 1, 8, infile);
    if ( !nr ) return -1;
    if ( png_sig_cmp(sig, 0, 8) ) return -1;   /* bad signature */

    /* could pass pointers to user-defined error handlers instead of NULLs: */
    png_structp			png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if ( !png_ptr ) return -4;   /* out of memory */

    png_infop			info_ptr = png_create_info_struct(png_ptr);
    if ( !info_ptr ) {
        png_destroy_read_struct(&png_ptr, NULL, NULL);
        return -4;   /* out of memory */
    }

    /* setjmp() must be called in every function that calls a PNG-reading
     * libpng function */
    if ( setjmp(png_jmpbuf(png_ptr)) ) {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        return -2;
    }

    png_init_io(png_ptr, infile);
    png_set_sig_bytes(png_ptr, 8);  /* we already read the 8 signature bytes */

    png_read_info(png_ptr, info_ptr);  /* read all PNG info up to image data */

    /* expand palette images to RGB, low-bit-depth grayscale images to 8 bits,
     * transparency chunks to full alpha channel; strip 16-bit-per-sample
     * images to 8 bits per sample; and convert grayscale to RGB[A] */

	int				bit_depth = png_get_bit_depth(png_ptr, info_ptr);
	int				color_type = png_get_color_type(png_ptr, info_ptr);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPNG: color_type=" << color_type << " bit_depth=" << bit_depth << endl;
	
    if ( color_type == PNG_COLOR_TYPE_PALETTE )
        png_set_palette_to_rgb(png_ptr);
    if ( color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8 )
        png_set_expand_gray_1_2_4_to_8(png_ptr);
    if ( png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS) )
        png_set_tRNS_to_alpha(png_ptr);
//    if ( bit_depth == 16 )
//      png_set_strip_16(png_ptr);
//    if ( color_type == PNG_COLOR_TYPE_GRAY ||
//			color_type == PNG_COLOR_TYPE_GRAY_ALPHA )
//      png_set_gray_to_rgb(png_ptr);

	png_set_interlace_handling(png_ptr);
	
    png_read_update_info(png_ptr, info_ptr);

	color_type = png_get_color_type(png_ptr, info_ptr);

    int				rowbytes = png_get_rowbytes(png_ptr, info_ptr);
	
	// Transferring the information
	p->size(png_get_image_width(png_ptr, info_ptr), png_get_image_height(png_ptr, info_ptr), 1);
	p->channels(png_get_channels(png_ptr, info_ptr));
	p->images(1);

	switch ( bit_depth) {
		case 1: p->data_type(Bit); break;
		case 8: p->data_type(UCharacter); break;
		case 16: p->data_type(UShort); break;
		default: p->data_type(UCharacter); break;
	}
	
	switch ( color_type ) {
		case PNG_COLOR_TYPE_GRAY: p->compound_type(TSimple); break;
		case PNG_COLOR_TYPE_RGB: p->compound_type(TRGB); break;
		case PNG_COLOR_TYPE_RGB_ALPHA: p->compound_type(TRGBA); break;
		default: break;
	}
	
	p->sampling(png_get_x_pixels_per_meter(png_ptr, info_ptr),
		png_get_y_pixels_per_meter(png_ptr, info_ptr), 1);
//	p->image->origin(png_get_x_offset_pixels(png_ptr, info_ptr), 
//		png_get_y_offset_pixels(png_ptr, info_ptr), 0);
	p->origin(png_get_x_offset_pixels(png_ptr, info_ptr), 
		png_get_y_offset_pixels(png_ptr, info_ptr), 0.0);

	tm*				t = p->get_localtime();
	png_time*		modtime;
	png_get_tIME(png_ptr, info_ptr, &modtime);
//	t->tm_year = modtime->year - 1900;
//	t->tm_mon = modtime->month - 1;
//	t->tm_mday = modtime->day;
//	t->tm_hour = modtime->hour;
//	t->tm_min = modtime->minute;
//	t->tm_sec = modtime->second;puts("***");
	p->set_time(t);
	
	if ( readdata ) {
		unsigned char*	data = p->data_alloc();

		png_bytep*		row_pointers = new png_bytep[p->sizeY()];
	
	// Row pointers are inverted with respect to y
		int				i, j;
		for ( i=0, j=p->sizeY()-1; i<p->sizeY(); i++, j--, data += rowbytes )
			row_pointers[j] = (png_bytep) data;

/*		if ( verbose & VERB_DEBUG )
			for ( i=0; i<p->sizeY(); i++ )
				cout << "DEBUG readPNG: row_pointers[" << i << "]=" << static_cast<void*>(row_pointers[i]) << endl;
*/
		png_read_image(png_ptr, row_pointers);

//		png_read_end(png_ptr, NULL);

		delete[] row_pointers;
	}
	
	png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
	
	fclose(infile);
	
	return 0;
}

/**
@brief	Writing an PNG or text image format.
@param	*p			the image structure.
@return	int			error code (<0 means failure).
A 2D image format commonly used on the Web.
	libpng version 1.2.8 - December 3, 2004
**/
int 	writePNG(Bimage* p)
{
	if ( p->data_type() == SCharacter ) p->change_type(UCharacter);
	if ( p->data_type() > UShort ) p->change_type(UShort);

    /* could also replace libpng warning-handler (final NULL), but no need: */

    png_structp		png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) return -4;   /* out of memory */

    png_infop		info_ptr = png_create_info_struct(png_ptr);
    if ( !info_ptr ) {
        png_destroy_write_struct(&png_ptr, NULL);
        return -4;   /* out of memory */
    }

    if ( setjmp(png_jmpbuf(png_ptr)) ) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        return -2;
    }

	FILE*			outfile;
	if ((outfile = fopen(p->file_name().c_str(), "wb")) == NULL) return -1;

    /* make sure outfile is (re)opened in BINARY mode */
	png_init_io(png_ptr, outfile);

    /* set the compression levels--in general, always want to leave filtering
     * turned on (except for palette images) and allow all of the filters,
     * which is the default; want 32K zlib window, unless entire image buffer
     * is 16K or smaller (unknown here)--also the default; usually want max
     * compression (NOT the default); and remaining compression flags should
     * be left alone */

	png_set_compression_level(png_ptr, Z_BEST_COMPRESSION);
/*
    >> this is default for no filtering; Z_FILTERED is default otherwise:
    png_set_compression_strategy(png_ptr, Z_DEFAULT_STRATEGY);
    >> these are all defaults:
    png_set_compression_mem_level(png_ptr, 8);
    png_set_compression_window_bits(png_ptr, 15);
    png_set_compression_method(png_ptr, 8);
 */


    /* set the image parameters appropriately */
	int				color_type = PNG_COLOR_TYPE_GRAY;
	switch ( p->compound_type() ) {
		case TSimple: color_type = PNG_COLOR_TYPE_GRAY; break;
		case TRGB: color_type = PNG_COLOR_TYPE_RGB; break;
		case TRGBA: color_type = PNG_COLOR_TYPE_RGB_ALPHA; break;
		default: color_type = PNG_COLOR_TYPE_GRAY; break;
	}
	
	int				bit_depth = 8*p->data_type_size();

//    interlace_type = mainprog_ptr->interlaced? PNG_INTERLACE_ADAM7 :
//                                             PNG_INTERLACE_NONE;
    int				interlace_type = PNG_INTERLACE_NONE;

    png_set_IHDR(png_ptr, info_ptr, p->sizeX(), p->sizeY(),
		bit_depth, color_type, interlace_type,
		PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	
	tm*				t = p->get_localtime();
	png_time		modtime;
	png_convert_from_struct_tm(&modtime, t);
//	png_convert_from_time_t(&modtime, p->time);	// Converts to GM time
	png_set_tIME(png_ptr, info_ptr, &modtime);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePNG: time=" << mktime(t) << endl;
	
/*
//    if (mainprog_ptr->have_time) {
//      png_time  modtime;
//
//        png_convert_from_time_t(&modtime, mainprog_ptr->modtime);
//        png_set_tIME(png_ptr, info_ptr, &modtime);
//    }

    if (mainprog_ptr->have_text) {
        png_text  text[6];
        int  num_text = 0;

        if (mainprog_ptr->have_text & TEXT_TITLE) {
            text[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
            text[num_text].key = "Title";
            text[num_text].text = mainprog_ptr->title;
            ++num_text;
        }
        if (mainprog_ptr->have_text & TEXT_AUTHOR) {
            text[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
            text[num_text].key = "Author";
            text[num_text].text = mainprog_ptr->author;
            ++num_text;
        }
        if (mainprog_ptr->have_text & TEXT_DESC) {
            text[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
            text[num_text].key = "Description";
            text[num_text].text = mainprog_ptr->desc;
            ++num_text;
        }
        if (mainprog_ptr->have_text & TEXT_COPY) {
            text[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
            text[num_text].key = "Copyright";
            text[num_text].text = mainprog_ptr->copyright;
            ++num_text;
        }
        if (mainprog_ptr->have_text & TEXT_EMAIL) {
            text[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
            text[num_text].key = "E-mail";
            text[num_text].text = mainprog_ptr->email;
            ++num_text;
        }
        if (mainprog_ptr->have_text & TEXT_URL) {
            text[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
            text[num_text].key = "URL";
            text[num_text].text = mainprog_ptr->url;
            ++num_text;
        }
        png_set_text(png_ptr, info_ptr, text, num_text);
    }
*/

    /* if we wanted to write any more text info *after* the image data, we
     * would set up text struct(s) here and call png_set_text() again, with
     * just the new data; png_set_tIME() could also go here, but it would
     * have no effect since we already called it above (only one tIME chunk
     * allowed) */


	// Optional significant bit (sBIT) chunk
/*	png_color_8 	sig_bit;
	
	// If we are dealing with a grayscale image then
	sig_bit.gray = bit_depth;
	
	png_set_sBIT(png_ptr, info_ptr, &sig_bit);

	png_set_shift(png_ptr, &sig_bit);  // to scale low-bit-depth values
*/

	/* set up the transformations:  for now, just pack low-bit-depth pixels
	 * into bytes (one, two or four pixels per byte) */
	png_set_packing(png_ptr);

	/* Optionally write comments into the image */
/*	{
//		png_infop 	write_info_ptr;
		png_text 	text_ptr[1];
		
		char key0[]="Title";
		char text0[]=p->label();
		text_ptr[0].key = key0;
		text_ptr[0].text = text0;
		text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
		text_ptr[0].itxt_length = 0;
		text_ptr[0].lang = NULL;
		text_ptr[0].lang_key = NULL;
		
		png_set_text(write_ptr, write_info_ptr, text_ptr, 1);
	}*/

	/* write all chunks up to (but not including) first IDAT */
	png_write_info(png_ptr, info_ptr);

	// Row pointers are inverted with respect to y
	png_bytep*		row_pointers = new png_bytep[p->sizeY()];
	
	int				i, j;
	int				rowbytes = p->sizeX()*p->channels()*p->data_type_size();
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePNG: rowbytes=" << rowbytes << endl;

	unsigned char*	data = p->data_pointer();
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePNG: data=" << static_cast<void*>(data) << endl;
	
	for ( i=0, j=p->sizeY()-1; i<p->sizeY(); i++, j--, data += rowbytes )
		row_pointers[j] = (png_bytep) data;

/*	if ( verbose & VERB_DEBUG )
		for ( i=0; i<p->sizeY(); i++ )
			cout << "DEBUG writePNG: row_pointers[" << i << "]=" << static_cast<void*>(row_pointers[i]) << endl;
*/
    png_write_image(png_ptr, row_pointers);

    png_write_end(png_ptr, NULL);

	png_destroy_write_struct(&png_ptr, &info_ptr);

	fclose(outfile);
	
	delete[] row_pointers;
	
	return 0;
}

#else

int 	readPNG(Bimage* p, int readdata)
{
	cerr << "Error: PNG files are not supported!" << endl;
	
	return -1;
}

int 	writePNG(Bimage* p)
{
	cerr << "Error: PNG files are not supported!" << endl;
	
	return -1;
}

#endif
