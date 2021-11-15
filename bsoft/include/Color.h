/**
@file	Color.h
@brief	Classes for color objects
@author Bernard Heymann
@date	Created: 20050112
@date	Modified: 20170420
**/

#include "Bstring.h"

template <typename Type> class RGBA;
template <typename Type> class CMYK;

#ifndef _COLOR_
#define _COLOR_
/************************************************************************
@Object: class RGB
@Description:
	Class for an RGB color object.
@Features:
	The internal variables are an array of 3 numbers.
*************************************************************************/
template <typename Type>
class RGB {
private:
	Type	data[3];
public:
	RGB()	{ for ( int i=0; i<3; i++ ) data[i] = 0; }
	RGB(const RGB& c) { for ( int i=0; i<3; i++ ) data[i] = c.data[i]; }
	RGB(const double d[3]) { for ( int i=0; i<3; i++ ) data[i] = d[i]; }
	RGB(const double rr, const double gg, const double bb) {
		data[0] = (Type) rr; data[1] = (Type) gg; data[2] = (Type) bb; }
	RGB	operator=(const RGB& c) {
		for ( int i=0; i<3; i++ ) data[i] = c.data[i];
		return *this;
	}
	RGB	operator-(const RGB& c) {
		for ( int i=0; i<3; i++ ) data[i] -= c.data[i];
		return *this;
	}
	Type&	operator[](const unsigned int i) {
		if ( i < 3 ) return data[i];
		else return data[0];
	}
	Type	r() { return data[0]; }
	Type	g() { return data[1]; }
	Type	b() { return data[2]; }
	void	red() { data[0] = 1; data[1] = data[2] = 0; }
	template <typename T2> RGB	operator=(RGBA<T2>& c) {
		for ( int i=0; i<3; i++ ) data[i] = c[i];
		return *this;
	}
	template <typename T2> operator RGB<T2>() const {
		return RGB<T2>(data[0], data[1], data[2]);
	}
	RGB(RGBA<Type>& c) { for ( int i=0; i<3; i++ ) data[i] = c[i]; }
	RGB(CMYK<Type>& cmyk) {
		double			r, g, b, k, k1;
		double			max = 255;
		if ( sizeof(Type) > 1 ) max = 1;		
		k = cmyk.k();
		k1 = 1 - k/max;
		r = cmyk.c()*k1 + k;
		g = cmyk.m()*k1 + k;
		b = cmyk[1]*k1 + k;
		if ( r > max ) r = max;
		if ( g > max ) g = max;
		if ( b > max ) b = max;
		data[0] = (Type) (max - r);
		data[1] = (Type) (max - g);
		data[2] = (Type) (max - b);
	}
	RGB(Bstring& hex) {
		int			i;
		for ( i=0; i<3; i++ ) data[i] = 255;		
		if ( hex != "white" )
			sscanf(hex.c_str(), "#%2x%2x%2x", &data[0], &data[1], &data[2]);		
		if ( sizeof(Type) > 1 ) for ( i=0; i<3; i++ ) data[i] /= 255.0;
	}
	Bstring	hex() {
		int			i, s=1;
		Bstring		hex("#");
		if ( sizeof(Type) > 1 ) s = 255;
		for ( i=0; i<3; i++ ) hex += Bstring(s*data[i], "%02x");
		return hex;
	}
	double	average() { return ((double)data[0] + (double)data[1] + (double)data[2])/3; }
	double	rms() { return sqrt(((double)data[0]*data[0] + (double)data[1]*data[1] + (double)data[2]*data[2])/3); }
	void	spectrum(double value, double cmin, double cmax ) {
		if ( fabs(cmin - cmax) < 1e-6 ) {
			cmin -= 1e-6;
			cmax += 1e-6;
		}
		
		double			type_max = ( sizeof(Type) > 1 )? 1: 255;
		double 			red, blu, grn;
		double			mid = 0.5*(cmax + cmin);
		double			scale = 2*type_max/(cmax - cmin);
		
		red = scale*(value - mid);
		if ( red < 0 ) red = 0;
		if ( red > type_max ) red = type_max;
		blu = scale*(mid - value);
		if ( blu < 0 ) blu = 0;
		if ( blu > type_max ) blu = type_max;
		grn = type_max - red - blu;
		if ( grn < 0 ) grn = 0;
		if ( grn > type_max ) grn = type_max;
		
		data[0] = (Type) red;
		data[1] = (Type) grn;
		data[2] = (Type) blu;
	}
	void	phase(double phi, double amp_ratio) {
		if ( amp_ratio > 255 ) amp_ratio = 255;
		if ( sizeof(Type) > 1 ) if ( amp_ratio > 1 ) amp_ratio = 1;		
		double		red, grn, blu;		
		
		phi /= M_PI;
		if ( phi > 0 ) {		// Positive angles
			red = 2 - 3*phi;
			grn = 3*phi;
			blu = 3*phi - 2;
		} else {				// Negative angles
			red = 2 + 3*phi;
			grn = -3*phi - 2;
			blu = -3*phi;
		}
		if ( red < 0 ) red = 0; // Clip to range 
		if ( grn < 0 ) grn = 0;
		if ( blu < 0 ) blu = 0;
		if ( red > 1 ) red = 1;
		if ( grn > 1 ) grn = 1;
		if ( blu > 1 ) blu = 1;
		
		data[0] = (Type) (amp_ratio*red);
		data[1] = (Type) (amp_ratio*grn);
		data[2] = (Type) (amp_ratio*blu);
	}
	void	random_color() {
		data[0] = (Type) (random()*255/RAND_MAX);
		data[1] = (Type) (random()*255/RAND_MAX);
		data[2] = (Type) (random()*255/RAND_MAX);
	}
} ;

template <typename Type>
ostream& operator<<(ostream& output, RGB<Type> rgb) {
	output.setf(ios::fixed, ios::floatfield);
	output << rgb.r() << tab << rgb.g() << tab << rgb.b();
	return output;
}

/************************************************************************
 @Object: class RGBA
 @Description:
 Class for an RGBA color object.
 @Features:
 The internal variables are an array of 3 numbers.
 *************************************************************************/
template <typename Type>
class RGBA {
private:
	Type	data[4];
public:
	RGBA()	{ for ( int i=0; i<4; i++ ) data[i] = 0; }
	RGBA(const RGBA& c) { for ( int i=0; i<4; i++ ) data[i] = c.data[i]; }
	RGBA(const double d[4]) { for ( int i=0; i<4; i++ ) data[i] = d[i]; }
	RGBA(const double rr, const double gg, const double bb) {
		data[0] = (Type) rr; data[1] = (Type) gg; data[2] = (Type) bb; data[3] = 0; }
	RGBA(const double rr, const double gg, const double bb, const double aa) {
		data[0] = (Type) rr; data[1] = (Type) gg; data[2] = (Type) bb; data[3] = (Type) aa; }
	RGBA	operator=(const RGBA& c) {
		for ( int i=0; i<4; i++ ) data[i] = c.data[i];
		return *this;
	}
	template <typename T2> RGBA	operator=(RGB<T2>& c) {
		for ( int i=0; i<3; i++ ) data[i] = c[i];
		return *this;
	}
	template <typename T2> operator RGBA<T2>() const {
		return RGBA<T2>(data[0], data[1], data[2], data[3]);
	}
	RGBA	operator+(const RGBA& c) {
		return RGBA(data[0] + c.data[0], data[1] + c.data[1], data[2] + c.data[2], data[3] + c.data[3]);
	} 
	RGBA	operator/(const double d) {
		return RGBA(data[0]/d, data[1] /d, data[2]/d, data[3]/d);
	} 
	RGBA(RGB<Type>& c) { for ( int i=0; i<3; i++ ) data[i] = c[i]; data[3] = 0; }
	RGBA(CMYK<Type>& cmyk) {
		RGB<Type>	rgb(cmyk);
		for ( int i=0; i<3; i++ ) data[i] = rgb[i];
		data[3] = 0;
	}
	Type&	operator[](const unsigned int i) {
		if ( i < 4 ) return data[i];
		else return data[0];
	}
	Type	r() { return data[0]; }
	Type	g() { return data[1]; }
	Type	b() { return data[2]; }
	Type	a() { return data[3]; }
	void	red() { data[0] = data[3] = 1; data[1] = data[2] = 0; }
	RGBA(Bstring& hex) {
		int				i;
		unsigned int	urgb[3];
		for ( i=0; i<4; i++ ) data[i] = 255;
		if ( hex != "white" )
			sscanf(hex.c_str(), "#%2x%2x%2x", &urgb[0], &urgb[1], &urgb[2]);		
		for ( i=0; i<3; i++ ) data[i] = urgb[i];
		if ( sizeof(Type) > 1 ) for ( i=0; i<4; i++ ) data[i] /= 255.0;
	}
	Bstring	hex() {
		int			i, s=1;
		Bstring		hex("#");
		if ( sizeof(Type) > 1 ) s = 255;
		for ( i=0; i<3; i++ ) hex += Bstring((int)(s*data[i]), "%02x");
		return hex;
	}
	void	spectrum(double value, double cmin, double cmax ) {
		if ( fabs(cmin - cmax) < 1e-6 ) {
			cmin -= 1e-6;
			cmax += 1e-6;
		}
		
		double			type_max = ( sizeof(Type) > 1 )? 1: 255;
		double 			red, blu, grn;
		double			mid = 0.5*(cmax + cmin);
		double			scale = 2*type_max/(cmax - cmin);
		
		red = scale*(value - mid);
		if ( red < 0 ) red = 0;
		if ( red > type_max ) red = type_max;
		blu = scale*(mid - value);
		if ( blu < 0 ) blu = 0;
		if ( blu > type_max ) blu = type_max;
		grn = type_max - red - blu;
		if ( grn < 0 ) grn = 0;
		if ( grn > type_max ) grn = type_max;
		
		data[0] = (Type) red;
		data[1] = (Type) grn;
		data[2] = (Type) blu;
		data[3] = (Type) type_max;
	}
} ;

template <typename Type>
ostream& operator<<(ostream& output, RGBA<Type> rgba) {
	output.setf(ios::fixed, ios::floatfield);
	output << rgba.r() << tab << rgba.g() << tab << rgba.b() << tab << rgba.a();
	return output;
}

/************************************************************************
 @Object: class CMYK
 @Description:
 Class for a CMYK color object.
 @Features:
 The internal variables are an array of 3 numbers.
 *************************************************************************/
template <typename Type>
class CMYK {
private:
	Type	data[4];
public:
	CMYK()	{ for ( int i=0; i<4; i++ ) data[i] = 0; }
	CMYK(const CMYK& col) { for ( int i=0; i<4; i++ ) data[i] = col.data[i]; }
	CMYK(const double d[4]) { for ( int i=0; i<4; i++ ) data[i] = d[i]; }
	CMYK(const vector<double>& d) { for ( int i=0; i<4 && i<d.size(); i++ ) data[i] = d[i]; }
	CMYK(const double cc, const double mm, const double yy, const double kk) {
		data[0] = (Type) cc; data[1] = (Type) mm; data[2] = (Type) yy; data[3] = (Type) kk; }
	CMYK	operator=(const CMYK& col) {
		for ( int i=0; i<4; i++ ) data[i] = col.data[i];
		return *this;
	}
	CMYK(RGB<Type>& rgb) {
		double			r1, g1, b1, k, k1, c, m, y;
		double			max = 255;
		if ( sizeof(Type) > 1 ) max = 1;				
		r1 = max - rgb.r();
		g1 = max - rgb.g();
		b1 = max - rgb.b();
		k = r1;
		if ( k > g1 ) k = g1;
		if ( k > b1 ) k = b1;
		k1 = max/(max - k);
		c = (r1 - k)*k1;
		m = (g1 - k)*k1;
		y = (b1 - k)*k1;
		if ( c > max ) c = max;
		if ( m > max ) m = max;
		if ( y > max ) y = max;
		data[0] = (Type) c; data[1] = (Type) m; data[2] = (Type) y; data[3] = (Type) k;
	}
	Type&	operator[](const unsigned int i) {
		if ( i < 4 ) return data[i];
		return data[0];
	}
	Type	c() { return data[0]; }
	Type	m() { return data[1]; }
	Type	y() { return data[2]; }
	Type	k() { return data[3]; }
} ;

template <typename Type>
ostream& operator<<(ostream& output, CMYK<Type> cmyk) {
	output.setf(ios::fixed, ios::floatfield);
	output << cmyk.c() << tab << cmyk.m() << tab << cmyk[1] << tab << cmyk.k();
	return output;
}

#endif

