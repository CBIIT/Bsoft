/**
@file	Bimage_color.cpp
@brief	Library routines for mangaing color images
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20200511
**/

#include "Bimage.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief A simple image is converted to a color image.

	All three color values in each voxel are the same (i.e., gray).

**/
void		Bimage::simple_to_rgb()
{
	if ( compoundtype == TRGB ) return;
	
	if ( compoundtype != TSimple ) {
		cerr << "Error in simple_to_rgb: Conversion must be from simple to RGB!" << endl;
		return;
	}

	double				v;
	double				scale = 255.0/(max - min);
	double				shift = -min;
	
	long		ds = x*y*z*n;
	RGB<unsigned char>*	rgb_data = new RGB<unsigned char>[ds];
	
	for ( long j=0; j<ds; j++ ) {
		v = ((*this)[j] + shift) * scale;
		rgb_data[j] = RGB<unsigned char>(v, v, v);
	}
	
	data_type(UCharacter);
	compoundtype = TRGB;
	c = 3;
	
	data_assign((unsigned char *) rgb_data);

	statistics();
}

/**
@brief A simple image is converted to a color image.

	All three color values in each voxel are the same (i.e., gray).
	The alpha value is set to 255.

**/
void		Bimage::simple_to_rgba()
{
	if ( compoundtype == TRGBA ) return;
	
	if ( compoundtype != TSimple ) {
		cerr << "Error in simple_to_rgba: Conversion must be from simple to RGBA!" << endl;
		return;
	}
	
	double					v;
	double					scale = 255.0/(max - min);
	double					shift = -min;
	
	long			ds = x*y*z*n;
	RGBA<unsigned char>*	rgba_data = new RGBA<unsigned char>[ds];
	
	for ( long j=0; j<ds; j++ ) {
		v = ((*this)[j] + shift) * scale;
		rgba_data[j] = RGBA<unsigned char>(v, v, v, 255);
	}
	
	data_type(UCharacter);
	compoundtype = TRGBA;
	c = 4;
	
	data_assign((unsigned char *) rgba_data);
	
	statistics();
}

/**
@brief A color image is converted to a simple image.

	The new value for each voxel is the average of the RGB color values.

**/
void		Bimage::color_to_simple()
{
	if ( compoundtype == TSimple ) return;
	if ( compoundtype == TComplex ) return;
	if ( compoundtype == TCMYK ) cmyk_to_rgb();
	if ( compoundtype == TRGBA ) rgba_to_rgb();
	
	long			ds = x*y*z*n;
	float*			nudata = new float[ds];

	for ( long j=0; j<ds; j++ )
		nudata[j] = (rgb(j)).average();
	
	data_type(Float);
	compoundtype = TSimple;
	c = 1;
	
	data_assign((unsigned char *) nudata);
	
	statistics();
}

/**
@brief An alpha channel is added to a RGB color image.

	The alpha channel value is set to 255.

**/
void		Bimage::rgb_to_rgba()
{
	if ( compoundtype == TRGBA ) return;
	
	if ( compoundtype != TRGB ) {
		cerr << "Error in Bimage::rgb_to_rgba: Conversion must be from RGB to RGBA!" << endl;
		return;
	}
	
	long			ds = x*y*z*n;
	RGBA<unsigned char>*	rgba_data = new RGBA<unsigned char>[ds];
	
	for ( long j=0; j<ds; j++ )
		rgba_data[j] = RGBA<unsigned char>((rgb(j)).r(), (rgb(j)).g(), (rgb(j)).b(), 255);
	
	data_type(UCharacter);
	compoundtype = TRGBA;
	c = 4;
	
	data_assign((unsigned char *) rgba_data);
	
	statistics();
}

/**
@brief The alpha channel is delete from an RGBA color image.
**/
void		Bimage::rgba_to_rgb()
{
	if ( compoundtype == TRGB ) return;
	
	if ( compoundtype != TRGBA ) {
		cerr << "Error in Bimage::rgba_to_rgb: Conversion must be from RGBA to RGB!" << endl;
		return;
	}
		
	long			ds = x*y*z*n;
	RGB<unsigned char>*		rgb_data = new RGB<unsigned char>[ds];
		
	for ( long j=0; j<ds; j++ )
		rgb_data[j] = RGB<unsigned char>((rgba(j)).r(), (rgba(j)).g(), (rgba(j)).b());
		
	data_type(UCharacter);
	compoundtype = TRGB;
	c = 3;
		
	data_assign((unsigned char *) rgb_data);
		
	statistics();
}

/**
@brief Converts an RGB color image to CMYK.
**/
void		Bimage::rgb_to_cmyk()
{
	if ( compoundtype == TCMYK ) return;
	
	if ( compoundtype != TRGB ) {
		cerr << "Error in Bimage::rgb_to_cmyk: Conversion must be from RGB to CMYK!" << endl;
		return;
	}
	
	long			ds = x*y*z*n;
	RGB<unsigned char>		rgb_data;
	CMYK<unsigned char>*	cmyk_data = new CMYK<unsigned char>[ds];
	
	for ( long j=0; j<ds; j++ ) {
		rgb_data = RGB<unsigned char>((rgb(j)).r(), (rgb(j)).g(), (rgb(j)).b());
		cmyk_data[j] = CMYK<unsigned char>(rgb_data);
	}
	
	data_type(UCharacter);
	compoundtype = TCMYK;
	c = 4;
	
	data_assign((unsigned char *) cmyk_data);
	
	statistics();
}

/**
@brief Converts a CMYK color image to RGB.
**/
void		Bimage::cmyk_to_rgb()
{
	if ( compoundtype == TRGB ) return;
	
	if ( compoundtype != TCMYK ) {
		cerr << "Error in Bimage::rgb_to_cmyk: Conversion must be from CMYK to RGB!" << endl;
		return;
	}
	
	long			ds = x*y*z*n;
	CMYK<unsigned char>		cmyk_data;
	RGB<unsigned char>*		rgb_data = new RGB<unsigned char>[ds];
	
	for ( long j=0; j<ds; j++ ) {
		cmyk_data = CMYK<unsigned char>((cmyk(j)).c(), (cmyk(j)).m(), (cmyk(j))[1], (cmyk(j)).k());
		rgb_data[j] = RGB<unsigned char>(cmyk_data);
	}
	
	data_type(UCharacter);
	compoundtype = TRGB;
	c = 3;
	
	data_assign((unsigned char *) rgb_data);
	
	statistics();
}

/**
@brief 	Converts a gray-scale image to a single color.
@param 	color		color selection (0=red, 1=green, 2=blue).
@param 	cmin		lower grayscale boundary.
@param 	cmax		upper grayscale boundary.
@param 	flag		sets the type of conversion.
@return int			0.

	A grayscale image is converted to a selected color between the given minimum and maximum.
	The flag parameter bit:
		bit 1		calculate a subtractive color range.
		bit 2		keep the ranges beyond the minimum and maximum as gray

**/
int 		Bimage::one_color(int color, double cmin, double cmax, int flag)
{
	if ( cmin == cmax ) {
		cmin = min;
		cmax = max;
	}
	
	if ( verbose & VERB_PROCESS )
	    cout << "Coloring the image: " << color << endl << endl;
	
    long		   	i, cc;
    double			gscale = 255.0/(max - min);
    double			scale = 255.0/(cmax - cmin);
	double			v;

	RGB<unsigned char>*	nudata = new RGB<unsigned char>[datasize];
	
	for ( i=0; i<datasize; i++ ) {
		v = (*this)[i];
	    if ( ( v >= cmin ) && ( v <= cmax ) ) {
      		v = scale*(v - cmin);
      		if ( v < 0 ) v = 0;
      		if ( v > 255 ) v = 255;
			if ( flag & 1 ) {
				for ( cc = 0; cc < 3; ++cc ) {
					if ( cc == color )
						nudata[i][color] = 255;
					else
						nudata[i][cc] = (unsigned char) (255 - v);
				}
			} else {
				nudata[i][color] = (unsigned char) v;
			}
		} else if ( flag & 2 ) {
			v = gscale*(v - min);
			nudata[i] = RGB<unsigned char>(v,v,v);
		} else if ( v > cmax ) {
			nudata[i][color] = 255;
		}
	}

    data_type(UCharacter);
    compound_type(TRGB);
    channels(3);

    data_assign((unsigned char *) nudata);
    
	return 0;
}

/**
@brief 	Converts a gray-scale image to red.
@param 	cmin		lower grayscale boundary.
@param 	cmax		upper grayscale boundary.
@param 	flag 		sets the type of conversion
@return int			0.

	A grayscale image is converted to red between the given minimum and maximum.

**/
int 		Bimage::color_red(double cmin, double cmax, int flag)
{
	return one_color(0, cmin, cmax, flag);
}

/**
@brief 	Converts a gray-scale image to green.
@param 	cmin		lower grayscale boundary.
@param 	cmax		upper grayscale boundary.
@param 	flag		sets the type of conversion
@return int			0.

	A grayscale image is converted to green between the given minimum and maximum.

**/
int 		Bimage::color_green(double cmin, double cmax, int flag)
{
	return one_color(1, cmin, cmax, flag);
}

/**
@brief 	Converts a gray-scale image to blue.
@param 	cmin		lower grayscale boundary.
@param 	cmax		upper grayscale boundary.
@param 	flag  		 sets the type of conversion
@return int			0.

	A grayscale image is converted to blue between the given minimum and maximum.

**/
int 		Bimage::color_blue(double cmin, double cmax, int flag)
{
	return one_color(2, cmin, cmax, flag);
}

/**
@brief 	Generates a pure color image without intensity.
@return int				0.

	Pure color is defined as:
		        col
		col = --------
		      sum(col)

**/
int 		Bimage::pure_color()
{
	if ( compoundtype != TRGB && compoundtype != TRGBA ) return 0;
	
	if ( verbose & VERB_LABEL )
	    cout << "Generating a pure color image" << endl << endl;
	
    long   		i, j, cc;
	double			sum, fac;
    
	for ( i=0; i<datasize; i+=c ) {
		for ( sum=cc=0, j=i; cc<3; cc++, j++ ) sum += (*this)[j];
		if ( sum ) {
			fac = dtmax/sum;
			for ( cc=0, j=i; cc<3; cc++, j++ )
				set(j, (*this)[j] * fac);
		}
	}
    
	return 0;
}

/**
@brief 	Combines two colored images.
@param 	*p				second image.
@return int				0.

	Each pair of voxels is summed and the result truncated to 255.

**/
int 		Bimage::color_combine(Bimage* p)
{
	if ( compoundtype != TRGB && compoundtype != TRGBA ) return 0;
	if ( p->compound_type() != TRGB && p->compound_type() != TRGBA ) return 0;
	
	if ( verbose & VERB_LABEL )
	    cout << "Combining two color images" << endl << endl;
	
    long		   i, j, cc;
	double			v;
    
	for ( i=0; i<datasize; i+=c ) {
		for ( cc=0, j=i; cc<c; cc++, j++ ) {
			v = (*this)[j] + (*p)[j];
			if ( v > 255 ) v = 255;
			set(j, v);
		}
	}
    
	return 0;
}


/**
@brief 	Colorizes an image to a spectrum.
@param 	cmin			lower grayscale boundary.
@param 	cmax			upper grayscale boundary.
@return Bimage*			color scale image.

	A grayscale image is converted to RGB and colored from blue (low) to 
	red (high). Only the values between the given minimum and maximum is used.

**/
Bimage*		Bimage::color_spectrum(double cmin, double cmax)
{
	change_type(Float);
	
	if ( cmin == cmax ) {
		cmin = min;
		cmax = max;
	}
	
	if ( verbose & VERB_PROCESS )
	    cout << "Coloring the image" << endl << endl;
	
    long				i;
	double				value;
    double				gscale(255.0/(max - min));
	
    RGB<unsigned char>* nudata = new RGB<unsigned char>[datasize];
	
	for ( i=0; i<datasize; i++ ) {
		value = (*this)[i*c];
	    if ( ( value >= cmin ) && ( value <= cmax ) ) {
			nudata[i].spectrum(value, cmin, cmax);
		} else {
			value = gscale*(value - min);
			nudata[i] = RGB<unsigned char>(value, value, value);
		}
	}

    data_type(UCharacter);
    compound_type(TRGB);
    channels(3);

    data_assign((unsigned char *) nudata);

    // Generate a color scale image
	if ( verbose & VERB_PROCESS )
	    cout << "Generating a color scale image" << endl << endl;

    long				j, k;
	RGB<unsigned char>	rgb;
	
	Bimage*				p = new Bimage(UCharacter, TRGB, 256, 20, 1, 1);

    for ( i=0; i<p->sizeX(); i++ ) {
		rgb.spectrum((double)i, 0.0, 255.0);
		for ( j=0, k=i; j<p->sizeY(); j++, k+=p->sizeX() ) p->set(k, rgb);
	}
	
	return p;
}

/**
@brief 	Colorizes an image with blue positive and red negative.
@param 	red_min			beginning of red gradient (most negative).
@param 	white_min 		end of red gradient (fade into white).
@param 	white_max 		beginning of blue gradient (start with white).
@param 	blue_max		end of blue gradient (most positive)
@return Bimage*			color scale image.

	A grayscale image is converted to RGB and colored with blue positive, 
	red negative, and white in between.

**/
Bimage*		Bimage::red_white_blue(double red_min, double white_min,
				double white_max, double blue_max)
{
	if ( verbose & VERB_PROCESS ) {
	    cout << "Generating a red-white-blue image:" << endl;
	    cout << "Red gradient start and end:     " << red_min << " " << white_min << endl;
	    cout << "Blue gradient start and end:    " << white_max << " " << blue_max << endl;
	} else if ( verbose & VERB_LABEL )
	    cout << "Generating a red-white-blue image" << endl << endl;
	
    long				i, j, k;
    int					red, blu, grn;
	double				v;
    double				blue_scale(255.0/(blue_max - white_max));
    double				red_scale(255.0/(white_min - red_min));
	RGB<unsigned char>	rgb;

    RGB<unsigned char>* nudata = new RGB<unsigned char>[datasize];
	
	for ( i=0; i<datasize; i++ ) {
		v = (*this)[i];
		red = grn = blu = 255;						// White part
		if ( v > blue_max ) {						// Blue part
			red = grn = 0;
		} else if ( v < red_min ) {					// Red part
			blu = grn = 0;
		} else if ( v > white_max ) {				// Blue gradient
      		red = (int) (blue_scale*(blue_max - v));
      		grn = (int) (blue_scale*(blue_max - v));
		} else if ( v < white_min ) {				// Red gradient
      		blu = (int) (red_scale*(v - red_min));
      		grn = (int) (red_scale*(v - red_min));
		}
      	if ( red < 0 ) red = 0;
      	if ( red > 255 ) red = 255;
      	if ( grn < 0 ) grn = 0;
      	if ( grn > 255 ) grn = 255;
      	if ( blu < 0 ) blu = 0;
      	if ( blu > 255 ) blu = 255;
		nudata[i] = RGB<unsigned char>(red, grn, blu);
	}

    data_type(UCharacter);
    compound_type(TRGB);
    channels(3);
	
    data_assign((unsigned char *) nudata);

    // Generate a color scale image
	if ( verbose & VERB_PROCESS )
	    cout << "Generating a color scale image: colorscale.tif" << endl << endl;
	
	Bimage* 	p = new Bimage(UCharacter, TRGB, 256, 20, 1, 1);
	double		scale((max - min)/255);
	
//	cout << "min=" << min << tab << "max=" << max << tab << "scale=" << scale << endl;
	
    for ( i=0; i<p->sizeX(); i++ ) {
		red = grn = blu = 255;					// White part
		v = scale*i + min;
		if ( v > blue_max ) {					// Blue part
			red = grn = 0;
		} else if ( v < red_min ) { 			// Red part
			blu = grn = 0;
		} else if ( v > white_max ) {			// Blue gradient
			red = (int) (blue_scale*(blue_max - v));
			grn = (int) (blue_scale*(blue_max - v));
		} else if ( v < white_min ) {			// Red gradient
			blu = (int) (red_scale*(v - red_min));
			grn = (int) (red_scale*(v - red_min));
		}
		if ( red < 0 ) red = 0;
		if ( red > 255 ) red = 255;
		if ( grn < 0 ) grn = 0;
		if ( grn > 255 ) grn = 255;
		if ( blu < 0 ) blu = 0;
		if ( blu > 255 ) blu = 255;
		rgb = RGB<unsigned char>(red, grn, blu);
		for ( j=0, k=i; j<p->sizeY(); j++, k+=p->sizeX() ) p->set(k, rgb);
	}
	
	return p;
}

/**
@brief 	Generates a phase color wheel.
@return int				0.

	The image is filled with a color wheel where the location (x,y) 
	is converted to polar form:
		a = arctan(y/x)
		r = sqrt(x^2+y^2)
	The color is a function of the angle a and saturation is a function
	of the distance from the origin.

**/
int			Bimage::phase_colour_wheel()
{
	long				i, xx, yy;
	double				rx, ry, a, r, rmax(256);

	origin(size()/2);
	RGB<unsigned char>	rgb;
	
	for ( i=yy=0; yy<y; yy++ ) {
		ry = yy - image->origin()[1];
		for ( xx=0; xx<x; xx++, i++ ) {
			rx = xx - image->origin()[0];
			a = atan2(ry, rx);
			r = sqrt(rx*rx + ry*ry);
			if ( r < rmax ) {
				rgb.phase(a, r);
				set(i, rgb);
			} else {
				set(i, RGB<unsigned char>(0,0,0));
			}
		}
	}
	
	return 0;
}


