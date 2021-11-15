/**
@file	ps_views.cpp
@brief	Postscript output for views
@author Bernard Heymann
@date	Created: 20011127
@date	Modified: 20201125
**/

#include "ps_views.h"
#include "ps_plot.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes

/**
@brief 	Generates postscript plot of views on the unit sphere.
@param 	&filename			output postscript file name.
@param 	&symmetry_string 	symmetry string to print at the top of the page.
@param 	*view				linked list of views.
@param 	flags				flags.
@return int					0.

	The plotting options are determined by the flags argument:
		0 = plot views
		1 = plot numbered views
		2 = plot views with shading according to occurrence.
			The gray level indicates an estimated increase in signal-to-noise ratio.

**/
int 		ps_views(Bstring& filename, Bstring& symmetry_string, View* view, int flags)
{
	if ( verbose & VERB_LABEL )
		cout << "Writing postscript file: " << filename << endl << endl;
	
	ofstream*	fps = ps_open_and_init(filename, symmetry_string, 1, 600, 800);
	
	*fps << "/Helvetica-Bold findfont 20 scalefont setfont" << endl;
	*fps << "40 740 moveto (" << filename << ") show" << endl;
	
	ps_views(fps, symmetry_string, view, flags);

	ps_close(fps);

	return 0;
}

int 		ps_views(ofstream* fps, Bstring& symmetry_string, View* view, int flags)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_views: symmetry=" << symmetry_string << " flags=" << flags << endl;
		
	int 		i, m, left(50), bottom(150), width(500), height(500), order(1);
	View*		v;

	int			h(100);
	int			ix, iy, nx(2*h+1), ny(2*h+1), sx(h), sy(h);
	double		isx(1.0L/sx), isy(1.0L/sy);
	float		snr(0.1);
	float*		count = NULL;

	long		n = count_list((char *) view);
	double		textscale(0.06), textoffset(0.02);
	if ( n > 50 ) textscale = 3.0/n;
	if ( textscale < 5 ) textscale = 5;
	textoffset = textscale/3;
	
	if ( flags == 2 ) {
		count = new float[nx*ny];
		for ( i=0; i<nx*ny; i++ ) count[i] = 0;
		for ( v = view, n=m=0; v; v = v->next ) {
			ix = (int) (v->x()*sx + sx);
			iy = (int) (v->y()*sy + sy);
			i = iy*nx+ix;
			count[i]++;
			if ( m < count[i] ) m = (int) count[i];
			n++;
		}
		snr = 1/sqrt((double)m);
		if ( snr < 0.1 ) snr = 0.1;
	}
	
	*fps << "%%%%Page: Views 1" << endl;
	*fps << "/Helvetica-Bold findfont 14 scalefont setfont" << endl;
	*fps << "40 720 moveto (Symmetry:        " << symmetry_string << ") show" << endl;
	*fps << "40 700 moveto (Number of views: " << n << ") show" << endl;
	*fps << "/Left " << left << " def\n/Bottom " << bottom << " def\n/Height " << height << " def\n/Width " << width << " def" << endl;
	*fps << "/Scale 200 def" << endl;
	
	*fps << "/View [" << endl;
	if ( flags == 2 ) {
		*fps << "%x y count weight" << endl;
		for ( i=iy=0; iy<ny; iy++ )
			for ( ix=0; ix<nx; ix++, i++ )
				if ( count[i] ) *fps << ix*isx - 1 << " " << iy*isy - 1 << " " << count[i] << " " << 1/(1+1/(snr*sqrt(count[i]))) << endl;
		delete[] count;
	} else {
		*fps << "%x y z a" << endl;
		for ( v = view; v; v = v->next )
			*fps << v->x() << " " << v->y() << " " << v->z() << " " << v->angle() << endl;
	}
	*fps << "] def" << endl;
		
	*fps << "/Geom { newpath " << endl;
	*fps << "	1 0 moveto 0 0 1 0 360 arc stroke" << endl;
	if ( symmetry_string[0] == 'T' ) {
		*fps << "	0 0 moveto 1 0 lineto" << endl;
		*fps << "	-0.577 0.577 moveto 0.577 -0.577 lineto stroke" << endl;
		*fps << "	-0.577 -0.577 moveto 0.577 0.577 lineto stroke" << endl;
		*fps << "	0 1 moveto 0 0.39 0.61 90 162 arc stroke" << endl;
		*fps << "	-0.577 0.577 moveto -0.39 0 0.61 108 180 arc stroke" << endl;
		*fps << "	1 0 moveto 0.39 0 0.61 0 72 arc stroke" << endl;
		*fps << "	0.577 0.577 moveto 0 0.39 0.61 18 90 arc stroke" << endl;
		*fps << "	0 -1 moveto 0 -0.39 0.61 -90 -18 arc stroke" << endl;
		*fps << "	0.577 -0.577 moveto 0.39 0 0.61 -72 0 arc stroke" << endl;
		*fps << "	-1 0 moveto -0.39 0 0.61 180 252 arc stroke" << endl;
		*fps << "	-0.577 -0.577 moveto 0 -0.39 0.61 198 270 arc stroke" << endl;
	} else if ( symmetry_string[0] == 'O' ) {
		*fps << "	-1 0 moveto 1 0 lineto" << endl;
		*fps << "	0 -1 moveto 0 1 lineto" << endl;
		*fps << "	0 0 moveto 0.707 0.707 lineto stroke" << endl;
		*fps << "	0 0 moveto 0.707 -0.707 lineto stroke" << endl;
		*fps << "	-0.577 0.577 moveto -0.707 0.707 lineto stroke" << endl;
		*fps << "	-0.577 -0.577 moveto -0.707 -0.707 lineto stroke" << endl;
		*fps << "	0.577 -0.577 moveto -0.637 0 1.344 -25.4 25.4 arc stroke" << endl;
		*fps << "	0.577 0.577 moveto 0 -0.637 1.344 64.6 115.4 arc stroke" << endl;
		*fps << "	-0.577 0.577 moveto 0.637 0 1.344 154.6 205.4 arc stroke" << endl;
		*fps << "	-0.577 -0.577 moveto 0 0.637 1.344 -115.4 -64.6 arc stroke" << endl;
	} else if ( symmetry_string[0] == 'I' ) {
		if ( symmetry_string.contains("I90") ) {
			*fps << "-0.526 0 moveto 0.526 0 lineto" << endl;
			*fps << "0 -0.357 lineto -0.526 0 lineto 0 0.357 lineto 0.526 0 lineto" << endl;
			*fps << "0 0.357 moveto 0 -0.357 lineto" << endl;
			*fps << "-0.851 -0.526 moveto" << endl;
			*fps << "-0.526 0 lineto" << endl;
			*fps << "0 -0.851 lineto" << endl;
			*fps << "-0.851 -0.526 lineto" << endl;
			*fps << "-0.851 0.526 lineto" << endl;
			*fps << "0 0.851 lineto" << endl;
			*fps << "-0.526 0 lineto" << endl;
			*fps << "-0.851 0.526 lineto" << endl;
			*fps << "0.851 -0.526 moveto" << endl;
			*fps << "0.526 0 lineto" << endl;
			*fps << "0 -0.851 lineto" << endl;
			*fps << "0.851 -0.526 lineto" << endl;
			*fps << "0.851 0.526 lineto" << endl;
			*fps << "0 0.851 lineto" << endl;
			*fps << "0.526 0 lineto" << endl;
			*fps << "0.851 0.526 lineto" << endl;
		} else {
			*fps << "0 0 moveto 0.357 0 lineto" << endl;
			*fps << "0 0.526 lineto 0 -0.526 lineto 0.357 0 lineto" << endl;
			*fps << "-0.526 -0.851 moveto" << endl;
			*fps << "0 -0.526 lineto" << endl;
			*fps << "-0.851 0 lineto" << endl;
			*fps << "-0.526 -0.851 lineto" << endl;
			*fps << "0.526 -0.851 lineto" << endl;
			*fps << "0.851 0 lineto" << endl;
			*fps << "0 -0.526 lineto" << endl;
			*fps << "0.526 -0.851 lineto" << endl;
			*fps << "-0.526 0.851 moveto" << endl;
			*fps << "0 0.526 lineto" << endl;
			*fps << "-0.851 0 lineto" << endl;
			*fps << "-0.526 0.851 lineto" << endl;
			*fps << "0.526 0.851 lineto" << endl;
			*fps << "0.851 0 lineto" << endl;
			*fps << "0 0.526 lineto" << endl;
			*fps << "0.526 0.851 lineto" << endl;
		}
	} else {
		sscanf(&symmetry_string[1], "%d", &order);
		*fps << "	0 0 moveto 1 0 lineto" << endl;
		for ( i=0; i<order; i++ )
			*fps << "	0 0 moveto " << cos((i+0.5)*M_PI*2.0/order) << " " 
				<< sin((i+0.5)*M_PI*2.0/order) << " lineto stroke" << endl;
	}
	*fps << "	closepath } def" << endl;

	*fps << "gsave" << endl;
	*fps << "	Width 2 div Left add Height 2 div Bottom add translate" << endl;
	*fps << "	Scale Scale scale" << endl;
	*fps << "	0.005 setlinewidth Geom stroke" << endl;
	*fps << "	/Helvetica findfont " << textscale << " scalefont setfont" << endl;
	if ( flags == 2 ) {
		*fps << "	0 0.1 1 {" << endl;
		*fps << "		/f exch def" << endl;
		*fps << "		0 4 View length 4 sub {" << endl;
		*fps << "			/Index exch def" << endl;
		*fps << "			/x View Index get def" << endl;
		*fps << "			/y View Index 1 add get def" << endl;
		*fps << "			/v View Index 3 add get def" << endl;
		*fps << "			v f ge {" << endl;
		*fps << "				1 v sub setgray" << endl;
		*fps << "				x y 0.01 0 360 arc fill" << endl;
		*fps << "			} if" << endl;
		*fps << "		} for stroke" << endl;
		*fps << "	} for" << endl;
	} else {
		*fps << "	/n 0 def" << endl;
		*fps << "	0 4 View length 4 sub {" << endl;
		*fps << "		/Index exch def" << endl;
		*fps << "		/x View Index get def" << endl;
		*fps << "		/y View Index 1 add get def" << endl;
		*fps << "		/z View Index 2 add get def" << endl;
		*fps << "		z 0 lt { x y 0.02 0 360 arc stroke }" << endl;
		*fps << "		{ x y 0.01 0 360 arc fill } ifelse" << endl;
		if ( flags == 1 ) {
			*fps << "		/n n 1 add def" << endl;
			*fps << "		x y moveto " << textoffset << " " << textoffset <<
				" rmoveto n cvi (xxxx) cvs show stroke" << endl;
		}
		*fps << "	} for stroke" << endl;
	}
	*fps << "grestore" << endl;
	*fps << "showpage" << endl;
	
	return 0;
}

int 		ps_views2(ofstream* fps, string& symmetry_string, list<View2<float>>& view, int flags)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_views: symmetry=" << symmetry_string << " flags=" << flags << endl;
		
	int 		i, m, left(50), bottom(150), width(500), height(500), order(1);

	int			h(100);
	int			ix, iy, nx(2*h+1), ny(2*h+1), sx(h), sy(h);
	double		isx(1.0L/sx), isy(1.0L/sy);
	float		snr(0.1);
	vector<float>	count(nx*ny,0);

	long		n = view.size();
	double		textscale(0.06), textoffset(0.02);
	if ( n > 50 ) textscale = 3.0/n;
	if ( textscale < 5 ) textscale = 5;
	textoffset = textscale/3;
	
	if ( flags == 2 ) {
		n = m = 0;
		for ( auto v = view.begin(); v != view.end(); ++v ) {
			ix = (int) (v->x()*sx + sx);
			iy = (int) (v->y()*sy + sy);
			i = iy*nx+ix;
			count[i]++;
			if ( m < count[i] ) m = (int) count[i];
			n++;
		}
		snr = 1/sqrt((double)m);
		if ( snr < 0.1 ) snr = 0.1;
	}
	
	*fps << "%%%%Page: Views 1" << endl;
	*fps << "/Helvetica-Bold findfont 14 scalefont setfont" << endl;
	*fps << "40 720 moveto (Symmetry:        " << symmetry_string << ") show" << endl;
	*fps << "40 700 moveto (Number of views: " << n << ") show" << endl;
	*fps << "/Left " << left << " def\n/Bottom " << bottom << " def\n/Height " << height << " def\n/Width " << width << " def" << endl;
	*fps << "/Scale 200 def" << endl;
	
	*fps << "/View [" << endl;
	if ( flags == 2 ) {
		*fps << "%x y count weight" << endl;
		for ( i=iy=0; iy<ny; iy++ )
			for ( ix=0; ix<nx; ix++, i++ )
				if ( count[i] ) *fps << ix*isx - 1 << " " << iy*isy - 1 << " " << count[i] << " " << 1/(1+1/(snr*sqrt(count[i]))) << endl;
	} else {
		*fps << "%x y z a" << endl;
		for ( auto v = view.begin(); v != view.end(); ++v )
			*fps << v->x() << " " << v->y() << " " << v->z() << " " << v->angle() << endl;
	}
	*fps << "] def" << endl;
		
	*fps << "/Geom { newpath " << endl;
	*fps << "	1 0 moveto 0 0 1 0 360 arc stroke" << endl;
	if ( symmetry_string[0] == 'T' ) {
		*fps << "	0 0 moveto 1 0 lineto" << endl;
		*fps << "	-0.577 0.577 moveto 0.577 -0.577 lineto stroke" << endl;
		*fps << "	-0.577 -0.577 moveto 0.577 0.577 lineto stroke" << endl;
		*fps << "	0 1 moveto 0 0.39 0.61 90 162 arc stroke" << endl;
		*fps << "	-0.577 0.577 moveto -0.39 0 0.61 108 180 arc stroke" << endl;
		*fps << "	1 0 moveto 0.39 0 0.61 0 72 arc stroke" << endl;
		*fps << "	0.577 0.577 moveto 0 0.39 0.61 18 90 arc stroke" << endl;
		*fps << "	0 -1 moveto 0 -0.39 0.61 -90 -18 arc stroke" << endl;
		*fps << "	0.577 -0.577 moveto 0.39 0 0.61 -72 0 arc stroke" << endl;
		*fps << "	-1 0 moveto -0.39 0 0.61 180 252 arc stroke" << endl;
		*fps << "	-0.577 -0.577 moveto 0 -0.39 0.61 198 270 arc stroke" << endl;
	} else if ( symmetry_string[0] == 'O' ) {
		*fps << "	-1 0 moveto 1 0 lineto" << endl;
		*fps << "	0 -1 moveto 0 1 lineto" << endl;
		*fps << "	0 0 moveto 0.707 0.707 lineto stroke" << endl;
		*fps << "	0 0 moveto 0.707 -0.707 lineto stroke" << endl;
		*fps << "	-0.577 0.577 moveto -0.707 0.707 lineto stroke" << endl;
		*fps << "	-0.577 -0.577 moveto -0.707 -0.707 lineto stroke" << endl;
		*fps << "	0.577 -0.577 moveto -0.637 0 1.344 -25.4 25.4 arc stroke" << endl;
		*fps << "	0.577 0.577 moveto 0 -0.637 1.344 64.6 115.4 arc stroke" << endl;
		*fps << "	-0.577 0.577 moveto 0.637 0 1.344 154.6 205.4 arc stroke" << endl;
		*fps << "	-0.577 -0.577 moveto 0 0.637 1.344 -115.4 -64.6 arc stroke" << endl;
	} else if ( symmetry_string[0] == 'I' ) {
		if ( symmetry_string[1] == '9' ) {
			*fps << "-0.526 0 moveto 0.526 0 lineto" << endl;
			*fps << "0 -0.357 lineto -0.526 0 lineto 0 0.357 lineto 0.526 0 lineto" << endl;
			*fps << "0 0.357 moveto 0 -0.357 lineto" << endl;
			*fps << "-0.851 -0.526 moveto" << endl;
			*fps << "-0.526 0 lineto" << endl;
			*fps << "0 -0.851 lineto" << endl;
			*fps << "-0.851 -0.526 lineto" << endl;
			*fps << "-0.851 0.526 lineto" << endl;
			*fps << "0 0.851 lineto" << endl;
			*fps << "-0.526 0 lineto" << endl;
			*fps << "-0.851 0.526 lineto" << endl;
			*fps << "0.851 -0.526 moveto" << endl;
			*fps << "0.526 0 lineto" << endl;
			*fps << "0 -0.851 lineto" << endl;
			*fps << "0.851 -0.526 lineto" << endl;
			*fps << "0.851 0.526 lineto" << endl;
			*fps << "0 0.851 lineto" << endl;
			*fps << "0.526 0 lineto" << endl;
			*fps << "0.851 0.526 lineto" << endl;
		} else {
			*fps << "0 0 moveto 0.357 0 lineto" << endl;
			*fps << "0 0.526 lineto 0 -0.526 lineto 0.357 0 lineto" << endl;
			*fps << "-0.526 -0.851 moveto" << endl;
			*fps << "0 -0.526 lineto" << endl;
			*fps << "-0.851 0 lineto" << endl;
			*fps << "-0.526 -0.851 lineto" << endl;
			*fps << "0.526 -0.851 lineto" << endl;
			*fps << "0.851 0 lineto" << endl;
			*fps << "0 -0.526 lineto" << endl;
			*fps << "0.526 -0.851 lineto" << endl;
			*fps << "-0.526 0.851 moveto" << endl;
			*fps << "0 0.526 lineto" << endl;
			*fps << "-0.851 0 lineto" << endl;
			*fps << "-0.526 0.851 lineto" << endl;
			*fps << "0.526 0.851 lineto" << endl;
			*fps << "0.851 0 lineto" << endl;
			*fps << "0 0.526 lineto" << endl;
			*fps << "0.526 0.851 lineto" << endl;
		}
	} else {
		sscanf(&symmetry_string[1], "%d", &order);
		*fps << "	0 0 moveto 1 0 lineto" << endl;
		for ( i=0; i<order; i++ )
			*fps << "	0 0 moveto " << cos((i+0.5)*M_PI*2.0/order) << " "
				<< sin((i+0.5)*M_PI*2.0/order) << " lineto stroke" << endl;
	}
	*fps << "	closepath } def" << endl;

	*fps << "gsave" << endl;
	*fps << "	Width 2 div Left add Height 2 div Bottom add translate" << endl;
	*fps << "	Scale Scale scale" << endl;
	*fps << "	0.005 setlinewidth Geom stroke" << endl;
	*fps << "	/Helvetica findfont " << textscale << " scalefont setfont" << endl;
	if ( flags == 2 ) {
		*fps << "	0 0.1 1 {" << endl;
		*fps << "		/f exch def" << endl;
		*fps << "		0 4 View length 4 sub {" << endl;
		*fps << "			/Index exch def" << endl;
		*fps << "			/x View Index get def" << endl;
		*fps << "			/y View Index 1 add get def" << endl;
		*fps << "			/v View Index 3 add get def" << endl;
		*fps << "			v f ge {" << endl;
		*fps << "				1 v sub setgray" << endl;
		*fps << "				x y 0.01 0 360 arc fill" << endl;
		*fps << "			} if" << endl;
		*fps << "		} for stroke" << endl;
		*fps << "	} for" << endl;
	} else {
		*fps << "	/n 0 def" << endl;
		*fps << "	0 4 View length 4 sub {" << endl;
		*fps << "		/Index exch def" << endl;
		*fps << "		/x View Index get def" << endl;
		*fps << "		/y View Index 1 add get def" << endl;
		*fps << "		/z View Index 2 add get def" << endl;
		*fps << "		z 0 lt { x y 0.02 0 360 arc stroke }" << endl;
		*fps << "		{ x y 0.01 0 360 arc fill } ifelse" << endl;
		if ( flags == 1 ) {
			*fps << "		/n n 1 add def" << endl;
			*fps << "		x y moveto " << textoffset << " " << textoffset <<
				" rmoveto n cvi (xxxx) cvs show stroke" << endl;
		}
		*fps << "	} for stroke" << endl;
	}
	*fps << "grestore" << endl;
	*fps << "showpage" << endl;
	
	return 0;
}

/**
@brief 	Generates postscript plot of views projected as phi and theta.
@param 	&filename		output postscript file name.
@param 	*view			linked list of views.
@return int				0.

**/
int 		ps_views(Bstring& filename, View* view)
{
	if ( !view ) return -1;
	
	Bstring			title("Views");
	ofstream*		fps = ps_open_and_init(filename, title, 1, 600, 800);

	ps_views(fps, view);
	
	ps_close(fps);
	
	return 0;
}

/**
@brief 	Generates postscript plot of views projected as phi and theta.
@param 	*fps			output postscript file stream.
@param 	*view			linked list of views.
@return int				0.

**/
int 		ps_views(ofstream* fps, View* view)
{
	int 			i(1), left(50), bottom(200), width(500), height(250);
	float			x, y;
	Euler			euler;
	
	View*			v;
	if ( !view ) return -1;
	
	*fps << "%%Page: " << i << " " << i << endl;
	*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
//	*fps << "50 755 moveto (Views: " << filename << ") show" << endl;
	*fps << "/Data [" << endl << "%x y fom" << endl;
	for ( v = view; v; v = v->next ) {
		euler = Euler(*v);
		x = euler.phi()*fabs(sin(euler.theta())) + M_PI;
		y = euler.theta();
		*fps << x << " " << y << " 0.5" << endl;
	}
	*fps << "] def" << endl;

	ps_phi_theta_plot(fps, left, bottom, width, height, 2);
	
	*fps << "showpage" << endl;
	
	return 0;
}

/**
@brief 	Postscript plot of the distribution of sets of views.
@param 	&filename		output postscript file name.
@param 	&title 			title.
@param 	nv				number of views.
@param 	*views			list of views.
@param 	ns				number of sets.
@param 	*fom			FOM values (nv*ns values).
@return int 				0, error if <0.
**/
int			ps_sets_of_views(Bstring& filename, Bstring& title, int nv, View* views, int ns, double* fom)
{
	int 			i, j, k, left(50), bottom(200), width(500), height(250);
	float			x, y, vx, vy;
	double			fmax;
	View*			v;
	Euler			euler;
	
	ofstream*		fps = ps_open_and_init(filename, title, ns, 600, 800);
	
	for ( i=1, k=0; i<=ns; i++ ) {
		*fps << "%%Page: " << i << " " << i << endl;
		*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
		*fps << "50 755 moveto (" << title << ": " << i << ") show" << endl;
		*fps << "/Data [" << endl << "%x y fom" << endl;
		for ( j=0, v=views, fmax=0; j<nv && v; j++, k++, v = v->next ) {
			euler = Euler(*v);
			x = euler.phi()*fabs(sin(euler.theta())) + M_PI;
			y = euler.theta();
			*fps << x << " " << y << " " << fom[k] << endl;
//			if ( j == i-1 ) {
			if ( fmax < fom[k] ) {
				fmax = fom[k];
				vx = x;
				vy = y;
			}
		}
		*fps << "] def" << endl;
		*fps << "/Best [" << vx << " " << vy << "] def" << endl;
		ps_phi_theta_plot(fps, left, bottom, width, height, 3);
		*fps << "showpage" << endl;
	}

	ps_close(fps);
	
	return 0;
}

/**
@brief 	Generates postscript plot of projected phi and theta.
@param 	*fps				postscript file stream.
@param 	left				left edge of plot.
@param 	bottom				bottom edge of plot.
@param 	width				width of plot.
@param 	height				height of plot.
@param 	ncol				number of columns in the Data table.
@return int					0.

	The plotting options are determined by the flags argument:
		0 = plot views
		1 = plot numbered views
		2 = plot views with shading according to occurrence.
			The gray level indicates an estimated increase in signal-to-noise ratio.

**/
int			ps_phi_theta_plot(ofstream* fps, int left, int bottom, int width, int height, int ncol)
{
	double			x, y, scale(width*1.0/TWOPI);
	double			phi, theta;
	
	*fps << "/Ncol " << ncol << " def" << endl;
	*fps << "/Scale " << scale << " def" << endl;
	*fps << "/Height " << height << " def" << endl << "/Width " << width << " def" << endl;
	*fps << "/Frame { newpath 0 0 moveto Width 0 lineto Width Height lineto 0 Height lineto closepath } def" << endl;

	*fps << left+width/2-10 << " " << bottom-30 << " moveto (Phi) show" << endl;
	*fps << left-30 << " " << bottom+height/2-20 << " moveto gsave 90 rotate (Theta) show grestore" << endl;
	
	*fps << "gsave" << endl;
	*fps << "	" << left << " " << bottom << " translate" << endl;
	*fps << "	Frame stroke" << endl;
//	*fps << "	Frame clip" << endl;
	*fps << "	Scale Scale scale" << endl;
	*fps << "	" << 1/scale << " setlinewidth" << endl;
	*fps << "	0.5 setgray" << endl;
	*fps << "	0 " << M_PI_2 << " moveto " << TWOPI << " " << M_PI_2 << " lineto" << endl;
	for ( phi = -M_PI; phi < M_PI+0.01; phi += M_PI/3.0 ) {
		*fps << "	" << M_PI << " 0 moveto" << endl;
		for ( theta = 0; theta < M_PI+0.01; theta += M_PI/90.0 ) {
			x = phi * fabs(sin(theta)) + M_PI;
			y = theta;
			*fps << "	" << x << " " << y << " lineto" << endl;
		}
	}
	*fps << "	stroke" << endl;
	*fps << "	/Helvetica findfont 64 scalefont setfont" << endl;
	*fps << "	0 setgray" << endl;
	*fps << "	Ncol Ncol Data length Ncol sub {" << endl;
	*fps << "		/Index exch def" << endl;
	*fps << "		/x Data Index get def" << endl;
	*fps << "		/y Data Index 1 add get def" << endl;
	*fps << "		/f Data Index 2 add get def" << endl;
	*fps << "		1 f sub setgray" << endl;
	*fps << "		x y 0.03 0 360 arc fill" << endl;
	*fps << "	} for" << endl;
	*fps << "	0 setgray" << endl;
	*fps << "	//userdict /Best known {" << endl;
	*fps << "		/x Best 0 get def" << endl;
	*fps << "		/y Best 1 get def" << endl;
	*fps << "		x y 0.05 0 360 arc stroke" << endl;
	*fps << "	} if" << endl;
	*fps << "grestore" << endl;
		
	ps_scale(fps, left, bottom, left+width, bottom, -180, 180, 60, 10, M_PI_2, 0, 12, 0);
	
	ps_scale(fps, left, bottom, left, bottom+height, 0, 180, 30, 10, 0, 0, 12, 0);
	
	return 0;
}

