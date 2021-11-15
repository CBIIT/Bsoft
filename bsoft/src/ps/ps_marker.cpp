/**
@file	ps_marker.cpp
@brief	Methods for postscript tools dealing with markers
@author Samuel Payne and Bernard Heymann
@date	Created: 20010725
@date	Modified: 20200312 (BH)
**/

#include "ps_marker.h"
#include "ps_plot.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Postscript plot of the errors associated with fiducial markers.
@param 	*fps		output stream.
@param 	&title 		title.
@param 	*project	project parameter structure.
@return int 				0, error if <0.
**/
int			ps_marker_plots(ofstream* fps, Bstring& title, Bproject* project)
{
	long				i(0), width(0), height(0), left(50), bottom(50), nmg(0), nmark(0);
	double				min_shift(1e10);
	
	Bfield*				field = project->field;
	Bmicrograph*		mg_ref = project->field->mg;
	Bmicrograph*		mg;
	Breconstruction*	rec = project->rec;
	Bmarker*			mark;
	Bimage*				p;
	Vector3<double>		shift;
	
	if ( verbose & VERB_PROCESS )
		cout << "Generating file with marker errors" << endl << endl;
	
	// Find the dimensions of the micrographs from the images
	for ( field=project->field; field; field=field->next ) {
		for ( mg=field->mg; mg; mg=mg->next ) {
//			if ( fabs(mg_ref->tilt_angle) > fabs(mg->tilt_angle) ) mg_ref = mg;
			shift = mg->origin - rec->origin;
			shift[2] = 0;
			if ( shift.length() < min_shift ) {
				min_shift = shift.length();
				mg_ref = mg;
				p = read_img(mg->fmg, 0, 0);
				if ( !p ) {
					error_show("Error in ps_marker_errors", __FILE__, __LINE__);
					cerr << "Micrograph " << mg->fmg << " not found" << endl;
					return -1;
				}
				if ( width  < p->sizeX() ) width  = p->sizeX();
				if ( height < p->sizeY() ) height = p->sizeY();
				delete p;
			}
			nmg++;
		}
	}
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_marker_errors: width=" << width << " height=" << height << endl;

	for ( nmark=0, mark=mg_ref->mark; mark; mark=mark->next, nmark++ ) ;
	
	double		scale(500.0/width);
	if ( scale > 600.0/height ) scale = 600.0/height;
	
	long		plot_width(scale*width);
	long		plot_height(scale*height);
	
	double		mark_rad(2.5/scale);
	double		mark_diam(2*mark_rad);
	double		err_scale(5/scale);

	*fps << setprecision(4) << "/Scale " << scale << " def" << endl;
	*fps << "/Height " << plot_height << " def\n/Width " << plot_width << " def" << endl;
	*fps << "/Frame { newpath 0 0 moveto Width 0 lineto Width Height lineto 0 Height lineto closepath } def" << endl;
	
	ps_define_arrowline(fps);
	
	*fps << "/Axis { newpath " << mark_rad/2 << " setlinewidth [10 10] 0 setdash "
		<< mg_ref->origin[0]*(1 - cos(mg_ref->tilt_axis)) << " " 
		<< mg_ref->origin[1]*(1 - sin(mg_ref->tilt_axis)) << " moveto "
		<< mg_ref->origin[0]*(1 + cos(mg_ref->tilt_axis)) << " "
		<< mg_ref->origin[1]*(1 + sin(mg_ref->tilt_axis)) << " lineto closepath} def" << endl;

	*fps << "/Ori { newpath " << mg_ref->origin[0] << " " << mg_ref->origin[1] << " 10 0 360 arc fill closepath} def" << endl;
	
	for ( field=project->field; field; field=field->next, i++ ) {
		*fps << "%%Page: " << i+1 << " " << i+1 << endl;
		*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
		*fps << "50 755 moveto (" << title << ": " << field->id << ") show" << endl;
		*fps << "50 735 moveto (Marker traces) show" << endl;

		*fps << "/Pos [" << endl << "%mg";
		for ( mark=mg_ref->mark; mark; mark=mark->next )
			*fps << " " << mark->id << "x " << mark->id << "y";
		*fps << endl;
		for ( mg=field->mg; mg; mg=mg->next ) {
			*fps << mg->img_num;
			shift = mg_ref->origin - mg->origin;
			for ( mark=mg->mark; mark; mark=mark->next )
				*fps << " " << mark->loc[0] + shift[0] << " " << mark->loc[1] + shift[1];
			*fps << endl;
		}
		*fps << "] def" << endl;

		*fps << "gsave" << endl;
		*fps << "	" << left << " " << bottom << " translate" << endl;
		*fps << "	Frame stroke" << endl;
		*fps << "	Scale Scale scale" << endl;
		*fps << "	/Helvetica findfont " << width/50 << " scalefont setfont" << endl;
		*fps << "	Ori stroke" << endl;
		*fps << "	Axis stroke" << endl;
		*fps << "	1 setlinewidth [] 0 setdash" << endl;
		*fps << "	/nm 1 def" << endl;
		*fps << "	1 1 " << nmark << " { " << endl;
		*fps << "		/Index exch def" << endl;
		*fps << "		/id Index def" << endl;
		*fps << "		/x Pos nm get def" << endl;
		*fps << "		/y Pos nm 1 add get def" << endl;
		*fps << "		x y moveto " << mark_diam << " " << mark_diam << " rmoveto id cvi (xxxx) cvs show stroke" << endl;
		*fps << "		x y " << mark_rad << " 0 360 arc closepath fill" << endl;
		*fps << "		x y moveto" << endl;
		*fps << "		" << 2*nmark+1 << " " << 2*nmark+1 << " Pos length " << 2*nmark+1 << " sub { " << endl;
		*fps << "			/Index2 exch def" << endl;
		*fps << "			/x Pos Index2 nm add get def" << endl;
		*fps << "			/y Pos Index2 nm add 1 add get def" << endl;
		*fps << "			x y lineto" << endl;
		*fps << "		} for stroke" << endl;
		*fps << "		x y " << mark_rad << " 0 360 arc closepath stroke" << endl;
		*fps << "		/nm nm 2 add def" << endl;
		*fps << "	} for" << endl;
		*fps << "grestore" << endl;
	
		*fps << "showpage" << endl;
	}
	
	for ( field=project->field; field; field=field->next, i++ ) {
		*fps << "%%Page: " << i+1 << " " << i+1 << endl;
		*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
		*fps << "50 755 moveto (" << title << ": " << field->id << ") show" << endl;
		*fps << "50 735 moveto (Marker errors (x" << err_scale << ")) show" << endl;

		*fps << "/Loc [\n%id x y" << endl;
		for ( nmark=0, mark=mg_ref->mark; mark; mark=mark->next, nmark++ )
			*fps << mark->id << " " << mark->loc[0] << " " << mark->loc[1] << endl;
		*fps << "] def" << endl;
		
		*fps << "/Err [" << endl << "%mg";
		for ( mark=mg_ref->mark; mark; mark=mark->next )
			*fps << " " << mark->id << "x " << mark->id << "y";
		*fps  << endl;
		for ( mg=field->mg; mg; mg=mg->next ) {
			*fps << mg->img_num; 
			for ( mark=mg->mark; mark; mark=mark->next )
				*fps << " " << mark->err[0] << " " << mark->err[1];
			*fps << endl;
		}
		*fps << "] def" << endl;

		*fps << "gsave" << endl;
		*fps << "	" << left << " " << bottom << " translate" << endl;
		*fps << "	Frame stroke" << endl;
		*fps << "	Scale Scale scale" << endl;
		*fps << "	/Helvetica findfont " << width/50 << " scalefont setfont" << endl;
		*fps << "	Ori stroke" << endl;
		*fps << "	Axis stroke" << endl;
		*fps << "	1 setlinewidth [] 0 setdash" << endl;
		*fps << "	/nm -1 def" << endl;
		*fps << "	0 3 Loc length 3 sub { " << endl;
		*fps << "		/Index exch def" << endl;
		*fps << "		/id Loc Index get def" << endl;
		*fps << "		/x Loc Index 1 add get def" << endl;
		*fps << "		/y Loc Index 2 add get def" << endl;
		*fps << "		/nm nm 2 add def" << endl;
		*fps << "		x y moveto " << mark_diam << " " << mark_diam << " rmoveto id cvi (xxxx) cvs show stroke" << endl;
		*fps << "		x y " << mark_rad << " 0 360 arc closepath fill" << endl;
		*fps << "		x y moveto" << endl;
		*fps << "		0 " << 2*nmark+1 << " Err length " << 2*nmark+1 << " sub { " << endl;
		*fps << "			/Index2 exch def" << endl;
		*fps << "			/ex Err Index2 nm add get " << err_scale << " mul x add def" << endl;
		*fps << "			/ey Err Index2 nm add 1 add get " << err_scale << " mul y add def" << endl;
		*fps << "			ex ey lineto" << endl;
		*fps << "		} for stroke" << endl;
		*fps << "	} for" << endl;
		*fps << "grestore" << endl;
	
		*fps << "showpage" << endl;
	}
	
	if ( project->rec ) {
		rec = project->rec;
		*fps << "%%Page: " << i+1 << " " << i+1 << endl;
		*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
		*fps << "50 755 moveto (" << title << ": " << rec->id << ") show" << endl;
		*fps << "50 735 moveto (Marker error sum vectors) show" << endl;
	
		*fps << "/Vec [" << endl << "%id x y ex ey" << endl;
		for ( mark=rec->mark; mark; mark=mark->next )
			*fps << mark->id << " " << mark->loc[0] << " " << 
				mark->loc[1] << " " << mark->err[0] << " " << mark->err[1] << endl;
		*fps << "] def" << endl;
		
		*fps << "gsave" << endl;
		*fps << "	" << left << " " << bottom << " translate" << endl;
		*fps << "	Frame stroke" << endl;
		*fps << "	Scale Scale scale" << endl;
		*fps << "	/Helvetica findfont " << width/50 << " scalefont setfont" << endl;
		*fps << "	Ori stroke" << endl;
		*fps << "	Axis stroke" << endl;
		*fps << "	1 setlinewidth [] 0 setdash" << endl;
		*fps << "	0 5 Vec length 5 sub { " << endl;
		*fps << "		/Index exch def" << endl;
		*fps << "		/id Vec Index get def" << endl;
		*fps << "		/x Vec Index 1 add get def" << endl;
		*fps << "		/y Vec Index 2 add get def" << endl;
		*fps << "		/ex Vec Index 3 add get x add def" << endl;
		*fps << "		/ey Vec Index 4 add get y add def" << endl;
		*fps << "		/d ex x sub ey y sub add def" << endl;
		*fps << "		x y moveto " << mark_diam << " " << mark_diam << " rmoveto id cvi (xxxx) cvs show stroke" << endl;
		*fps << "		x y " << mark_rad << " 0 360 arc closepath fill" << endl;
		*fps << "		d 0 ne {" << endl;
		*fps << "			x y ex ey arrowline" << endl;
//		*fps << "			ex ey x y arrowline" << endl;
		*fps << "		} if" << endl;
		*fps << "	} for" << endl;
		*fps << "grestore" << endl;
	
		*fps << "showpage" << endl;
	}
	
	return 0;
}


/**
@author  Samuel Payne & Bernard Heymann
@brief 	Output a postscript file that graphs the deviations of points in a 
	Bmarker from their predicted location (using the transform
	parameters).
@param	filename
@param 	*set1		the first set of point(refset).
@param 	*set2		the second set of point(applyset).
@param 	t			the transform parameters.
@param 	size		size of image or frame around markers.
@param 	err_scale 	scaling factor for the error.
@return int 0
**/
int			ps_marker_errors(Bstring& filename, Bmarker* set1, Bmarker* set2, 
				Transform t, Vector3<long> size, double err_scale)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_marker_errors: Creating " << filename << endl;
	
	Bstring		title("Matched deviations");
	
	ofstream*	fps = ps_open_and_init(filename, title, 1, 600, 800);

	*fps << "/Helvetica-Bold findfont 18 scalefont setfont" << endl;
	*fps << "40 770 moveto (" << filename << ": " << title << ") show" << endl;
	
	ps_marker_errors(fps, set1, set2, t, size, err_scale);
	
	ps_close(fps);
	
	return 0;
}

int			ps_marker_errors(ofstream* fps, Bmarker* set1, Bmarker* set2, 
				Transform t, Vector3<long> size, double err_scale)
{
	if ( !set1 || !set2 ) {
		cerr << "Error: The two sets must have the same number of points!" << endl;
		return -1;
	}
	
	long			num(0);
	double			xdev(0), ydev(0);
	Bmarker			*m1, *m2;
	Vector3<double>	min = set1->loc;
	Vector3<double>	max = set1->loc;
	
	//determine the deviation of every point in set2 from its predicted 
	//coordinates using the tranfsform parameters
	for ( m1=set1, m2=set2; m1; m1=m1->next, m2=m2->next ) {
		min = min.min(m1->loc);
		max = max.max(m1->loc);
		m2->err[0] = t.scale[0]*(m1->loc[0] *cos(t.angle) + 
				m1->loc[1]*sin(t.angle)) + t.trans[0] - m2->loc[0];
		m2->err[1] = t.scale[1]*(-1* m1->loc[0] *sin(t.angle) + 
				m1->loc[1]*cos(t.angle)) + t.trans[1] - m2->loc[1];
		xdev += m2->err[0]*m2->err[0];
		ydev += m2->err[1]*m2->err[1];
		num++;
	}
	xdev = sqrt(xdev/num);
	ydev = sqrt(ydev/num);
	
	long 		w = size[0];
	long 		h = size[1];
	double 		scale = 650.0/h;		
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_marker_errors: xdev " << xdev << " ydev " << ydev << endl;
	
	*fps << "/Helvetica-Bold findfont 10 scalefont setfont" << endl;
	*fps << "40 60 moveto (Shift: " << t.trans[0] << " " << t.trans[1] << ") show" << endl;
	*fps << "40 45 moveto (Scale: " << t.scale[0] << " " << t.scale[1] << ") show" << endl;
	*fps << "40 30 moveto (Angle: " << t.angle*180.0/M_PI << ") show" << endl;
	*fps << "40 15 moveto (Deviation: " << xdev << " " << ydev << " (scaled by " << err_scale << ")) show" << endl;
	*fps << "/x_Origin 100 def \n/y_Origin 100 def \n/Height " << h << " def \n/Width " << w << " def" << endl;
	*fps << "/Scale " << scale << " def" << endl;
	*fps << "/err_scale " << err_scale << " def" << endl;
	*fps << "/Frame { newpath 0 0 moveto 0 Height lineto Width Height lineto Width 0 lineto closepath } def" << endl;
		
	*fps << "/Data [" << endl;
	for ( m2=set2; m2; m2=m2->next )
		*fps << m2->loc[0] << " " << m2->loc[1] << " " << m2->err[0] << " " << m2->err[1] << endl;
	*fps << "] def\n" << endl;
	
	*fps << "300 80 moveto (x axis) show" << endl;
	*fps << "80 300 moveto gsave 90 rotate (y axis) show grestore" << endl;
	*fps << "100 80 moveto 0 (xxxxx) cvs show\n 80 100 moveto 0 (xxxxx) cvs show" << endl;
	*fps << w*scale+80 << " 80 moveto Width (xxxxx) cvs show" << endl;
	*fps << "60 " << h*scale+100 << " moveto Height (xxxxx) cvs show" << endl;
	*fps << "x_Origin y_Origin translate\nScale Scale scale" << endl;
	*fps << "gsave \n\t5 setlinewidth \n\tFrame stroke \n\t1 setlinewidth" << endl;
	*fps << "\t0 4 Data length 1 sub {  \n\t\t/Index exch def" << endl;
	*fps << "\t\tData Index get Data Index 1 add get err_scale 2 div 0 360 arc fill" << endl;
	*fps << "\t\tData Index get Data Index 1 add get moveto" << endl;
	*fps << "\t\tData Index 2 add get err_scale mul Data Index 3 add get err_scale mul rlineto stroke" << endl;
	*fps << "\t} for \ngrestore" << endl;
	
	*fps << "showpage" << endl;	
	
	return 0;
}

/**
@author  Samuel Payne
@brief 	Output a postscript file that draws a line between matching features
	in a pair of coordinate sets.
@param 	*set1		the first set of coordinates
@param 	*set2		second set of coordinates
@param 	&filename			name of file to be written
@return int 0
**/
int 		ps_marker_match(Bmarker* set1, Bmarker* set2, Bstring& filename) 
{ 
	if ( !set1 || !set2 ) {
		cerr << "Error: The two sets must have the same number of points!" << endl;
		return -1;
	}
	
	long				num = 0;
	Bmarker			*m1, *m2;
	Vector3<double>	min = set1->loc;
	Vector3<double>	max = set1->loc;
	
	//determine the deviation of every point in set2 from its predicted 
	//coordinates using the tranfsform parameters
	for ( m1=set1, m2=set2; m1; m1=m1->next, m2=m2->next ) {
		min = min.min(m1->loc);
		max = max.max(m1->loc);
		num++;
	}
	
    long 		w = (long) (max[0] - min[0]);
    long 		h = (long) (max[1] - min[1]);
	
	double 		scale = 540.0/w;		
	
	Bstring		title("Matched Features");
	ofstream*	fps = ps_open_and_init(filename, title, 1, 600, 800);

	*fps << "/Helvetica-Bold findfont 18 scalefont setfont" << endl;
	*fps << "40 770 moveto \n(Matched Features) show" << endl;
	*fps << "/x_Origin 50 def \n/y_Origin 50 def \n/Height " << h << " def \n/Width " << w << " def" << endl;
	*fps << "/Frame { newpath 0 0 moveto 0 Height lineto Width Height lineto Width 0 lineto closepath } def" << endl;
	*fps << "/Dot { newpath 0 0 moveto 0 0 5 0 360 arc fill closepath} def" << endl;
	*fps << "/Ring { newpath 0 0 moveto 0 0 5 0 360 arc stroke closepath} def" << endl;
	*fps << "/Scale " << scale << " def" << endl;
		
	*fps << "/Data [" << endl;
	for ( m1=set1, m2=set2; m1; m1=m1->next, m2=m2->next )
		*fps << m1->loc[0] << " " << m1->loc[1] << " " << m2->loc[0] << " " << m2->loc[1] << endl;
	*fps << "] def\n" << endl;
	
	*fps << "300 20 moveto (x axis) show" << endl;
	*fps << "30 350 moveto gsave 90 rotate (y axis) show grestore" << endl;
	*fps << "50 20 moveto 0 (xxxx) cvs show\n 30 50 moveto 0 (xxxx) cvs show" << endl;
	*fps << w*scale << " 20 moveto " << w << " (xxxx) cvs show" << endl;
	*fps << "10 " << h*scale << " moveto " << h << " (xxxx) cvs show" << endl;
	*fps << "x_Origin y_Origin translate\nScale Scale scale" << endl;
	*fps << "gsave \n\t5 setlinewidth \n\tFrame stroke \n\t1 setlinewidth" << endl;
	*fps << "\t0 4 Data length 1 sub {  \n\t\t/Index exch def" << endl;
	*fps << "\t\tData Index get Data Index 1 add get moveto" << endl;
	*fps << "\t\tData Index 2 add get Data Index 3 add get lineto stroke" << endl;
	*fps << "\t\tgsave \n\t\t\tData Index get Data Index 1 add get translate \n\t\t\tDot stroke" << endl;
	*fps << "\t\tgrestore" << endl;
	*fps << "\t\tgsave \n\t\t\tData Index 2 add get Data Index 3 add get translate \n\t\t\tRing stroke" << endl;
	*fps << "\t\tgrestore \n\t} for \ngrestore" << endl;
	
	*fps << "showpage" << endl;
	
	ps_close(fps);
	
	return 0;
}


/**
@author Jessica Mavadia, Bernard Heymann
@date	20080610 - 20100211
@brief 	Generates a PostScript file that displays X by X, Y by Y, 
	and Z by Z plots for diagnostic purposes.
@param 	&filename	PostScript file name.
@param 	*project	project with two tilt series reconstructions.
@return int					0, <0 on error
**/
int		ps_dual_zcompare(Bstring& filename, Bproject* project)
{
	if ( !project->rec ) return -1;
	if ( !project->rec->next ) return -1;

	Breconstruction*	rec1 = project->rec;
	Breconstruction*	rec2 = project->rec->next;
	Bmarker*			mark1 = NULL;
	Bmarker*			mark2 = NULL;
	Bstring				x_label("X coordinate");
	Bstring				y_label("Y coordinate");
	Bstring				z_label("Z coordinate");
	long				i, npt;
	
	Vector3<double>	cmin, cmax;
	
	for ( npt=0, mark1 = rec1->mark; mark1; mark1 = mark1->next ) if ( mark1->sel ) {
		for ( mark2 = rec2->mark; mark2 && mark2->id != mark1->id; mark2 = mark2->next ) ;
		if ( mark2 && mark2->sel ) npt++;
	}
	
	if ( verbose )
		cout << "Generating plot " << filename << endl << endl;

	Bstring			title("Comparing dual tilt marker positions");
	Bplot*			plot = new Bplot(3, npt, 6);
	plot->title(title);
	
	for ( i=0, mark1 = rec1->mark; mark1; mark1 = mark1->next ) if ( mark1->sel ) {
		for ( mark2 = rec2->mark; mark2 && mark2->id != mark1->id; mark2 = mark2->next ) ;
		if ( mark2 && mark2->sel ) {
			cmin = cmin.min(mark1->loc);
			cmin = cmin.min(mark2->loc);
			cmax = cmax.max(mark1->loc);
			cmax = cmax.max(mark2->loc);
			(*plot)[i]			= mark1->loc[0];
			(*plot)[i+npt]		= mark1->loc[1];
			(*plot)[i+2*npt]	= mark1->loc[2];
			(*plot)[i+3*npt]	= mark2->loc[0];
			(*plot)[i+4*npt]	= mark2->loc[1];
			(*plot)[i+5*npt]	= mark2->loc[2];
			i++;
		}
	}
//	cout << "Pages: " << plot->pages() << endl;
	for ( i=0; i<3; i++ ) { // cout << "Setting up page " << plot->page(i).number() << endl;
		plot->page(i).title(title);
		plot->page(i).columns(2);
		plot->page(i).column(0).axis(1);
		plot->page(i).column(1).axis(3);
		plot->page(i).column(1).type(0);
		plot->page(i).column(1).element_size(3);
		plot->page(i).axis(1).min(cmin[i]);
		plot->page(i).axis(1).max(cmax[i]);
		plot->page(i).axis(3).min(cmin[i]);
		plot->page(i).axis(3).max(cmax[i]);
	}
	plot->page(0).column(0).label(x_label);
	plot->page(0).column(1).label(x_label);
	plot->page(1).column(0).label(y_label);
	plot->page(1).column(1).label(y_label);
	plot->page(2).column(0).label(z_label);
	plot->page(2).column(1).label(z_label);
	plot->page(0).column(0).number(0);
	plot->page(0).column(1).number(3);
	plot->page(1).column(0).number(1);
	plot->page(1).column(1).number(4);
	plot->page(2).column(0).number(2);
	plot->page(2).column(1).number(5);
	
//	cout << "Zcompare setup done " << endl;
	
	ps_plot(filename, plot);
	
	delete plot;
	
	return 0;
}
