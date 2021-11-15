/**
@file	ps_micrograph.cpp
@brief	Postscript output for micrographs
@author Bernard Heymann
@date	Created: 20011127
@date	Modified: 20210706
**/

#include "ps_micrograph.h"
#include "ps_plot.h"
#include "ps_views.h"
#include "rwimg.h"
#include "symmetry.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
int			ps_origins(ofstream* fps, Bstring& title, Bproject* project, int flags);

/**
@brief 	Postscript plot of the micrograph origins.
@param 	&filename	output postscript file name.
@param 	&title 		title.
@param 	*project	project parameter structure.
@return int 			0, error if <0.
**/
int			ps_mg_origins(Bstring& filename, Bstring& title, Bproject* project)
{
	int 			i, width(0), height(0), left=(50), bottom=(50), nmg(0);
	
	Bfield*			field = project->field;
	if ( !field ) return -1;
	
	Bmicrograph*	mg = field->mg;
	if ( !mg ) return -1;
	
	Bimage*			p;
	
	// Find the dimensions of the micrographs from the coordinate sets
	for ( field=project->field; field; field=field->next ) {
		for ( mg=field->mg; mg; mg=mg->next ) {
			p = read_img(mg->fmg, 0, 0);
			if ( !p ) {
				error_show("Error in ps_mg_origins", __FILE__, __LINE__);
				cerr << "Micrograph " << mg->fmg << " not found" << endl;
				return -1;
			}
			if ( width  < p->sizeX() ) width  = p->sizeX();
			if ( height < p->sizeY() ) height = p->sizeY();
			delete p;
			nmg++;
		}
	}
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_mg_origins: width=" << width << " height=" << height << endl;
	
	float		scale(500.0/width);
	if ( scale > 700.0/height ) scale = 700.0/height;
	width *= scale;
	height *= scale;

	
	ofstream*	fps = ps_open_and_init(filename, title, 1, 600, 800);

	*fps << "%%%%Page: " << 1 << " " << 1 << endl;
	*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
	*fps << "50 755 moveto (" << title << ") show" << endl;
	*fps << "/Scale " << scale << " def" << endl;
	*fps << "/Height " << height << " def\n/Width " << width << " def" << endl;
	*fps << "/Frame { newpath 0 0 moveto Width 0 lineto Width Height lineto 0 Height lineto closepath } def" << endl;
	*fps << "/Data [\n%# ox oy" << endl;
	
	for ( field=project->field; field; field=field->next )
		for ( i=1, mg=field->mg; mg; mg=mg->next, ++i )
			*fps << i << " " << mg->origin[0] << " " << mg->origin[1] << endl;
	
	*fps << "] def" << endl;

	*fps << "gsave" << endl;
	*fps << "	" << left << " " << bottom << " translate" << endl;
	*fps << "	Frame stroke" << endl;
	*fps << "	Scale Scale scale" << endl;
	*fps << "	/Helvetica findfont 64 scalefont setfont" << endl;
	*fps << "	/x Data 1 get def" << endl;
	*fps << "	/y Data 2 get def" << endl;
	*fps << "	x y 5 0 360 arc closepath fill" << endl;
	*fps << "	x y moveto" << endl;
	*fps << "	3 3 Data length 3 sub { " << endl;
	*fps << "		/Index exch def" << endl;
	*fps << "		/i Data Index get def" << endl;
	*fps << "		/x Data Index 1 add get def" << endl;
	*fps << "		/y Data Index 2 add get def" << endl;
	*fps << "		x y lineto" << endl;
/*	*fps << "	/x Data 0 get x_Min sub x_Scale mul def" << endl;
	*fps << "	/y Data 1 get y_Min sub y_Scale mul def" << endl;
	*fps << "	x y moveto" << endl;
	*fps << "	3 3 Data length 3 sub {" << endl;
	*fps << "		/x_Index exch def" << endl;
	*fps << "		/y_Index x_Index 1 add def" << endl;
	*fps << "		/x Data x_Index get x_Min sub x_Scale mul def" << endl;
	*fps << "		/y Data y_Index get y_Min sub y_Scale mul def" << endl;
	*fps << "		x y lineto" << endl;*/
	*fps << "	} for stroke" << endl;
	*fps << "grestore" << endl;
	
	*fps << "showpage" << endl;
	
	ps_close(fps);
	
	return 0;
}

/**
@brief 	Postscript plot of the positions of particles on micrographs.
@param 	&filename	output postscript file name.
@param 	&title 		title.
@param 	*project	project parameter structure.
@return int 			0, error if <0.
**/
int			ps_mg_particle_positions(Bstring& filename, Bstring& title, Bproject* project)
{
	int 			i, width(0), height(0), left=(50), bottom=(50), nmg(0);
	
	Bfield*			field = project->field;
	if ( !field ) return -1;
	
	Bmicrograph*	mg = field->mg;
	if ( !mg ) return -1;
	
	Bparticle*		part = mg->part;
	if ( !part ) return -1;
   
	Bimage*			p;
	
	// Find the dimensions of the micrographs from the coordinate sets
	for ( field=project->field; field; field=field->next ) {
		for ( mg=field->mg; mg; mg=mg->next ) {
			p = read_img(mg->fmg, 0, 0);
			if ( !p ) {
				error_show("Error in ps_mg_particle_positions", __FILE__, __LINE__);
				cerr << "Micrograph " << mg->fmg << " not found" << endl;
				return -1;
			}
			if ( width  < p->sizeX() ) width  = p->sizeX();
			if ( height < p->sizeY() ) height = p->sizeY();
			delete p;
			nmg++;
		}
	}
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_mg_particle_positions: width=" << width << " height=" << height << endl;
	
	float		scale(500.0/width);
	if ( scale > 700.0/height ) scale = 700.0/height;
	
	ofstream*	fps = ps_open_and_init(filename, title, nmg, 600, 800);
	
	for ( field=project->field; field; field=field->next ) {
		for ( i=0, mg=field->mg; mg; mg=mg->next, i++ ) {
			*fps << "%%%%Page: " << i+1 << " " << i+1 << endl;
			*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
			*fps << "50 755 moveto (" << title << ": " << mg->fpart << ") show" << endl;
			*fps << "/Scale " << scale << " def" << endl;
			*fps << "/Height 700 def\n/Width 500 def" << endl;
			*fps << "/Frame { newpath 0 0 moveto Width 0 lineto Width Height lineto 0 Height lineto closepath } def" << endl;
	
			*fps << "/Data [\n%x y mag sel" << endl;
			for ( part=mg->part; part; part=part->next ) if ( part->sel )
				*fps << part->loc[0] << " " << part->loc[1] << " " << 
						part->mag << " " << part->sel << endl;
			*fps << "] def" << endl;

			*fps << "gsave" << endl;
			*fps << "	" << left << " " << bottom << " translate" << endl;
			*fps << "	Frame stroke" << endl;
//			*fps << "	Frame clip" << endl;
			*fps << "	Scale Scale scale" << endl;
			*fps << "	/Helvetica findfont 64 scalefont setfont" << endl;
			*fps << "	0 4 Data length 4 sub { " << endl;
			*fps << "		/Index exch def" << endl;
			*fps << "		/x Data Index get def" << endl;
			*fps << "		/y Data Index 1 add get def" << endl;
			*fps << "		/mag Data Index 2 add get def" << endl;
			*fps << "		/sel Data Index 3 add get def" << endl;
			*fps << "		x y moveto 0 30 rmoveto sel cvi (xxxx) cvs show stroke" << endl;
			*fps << "		x y moveto 0 -30 rmoveto mag cvr (x.xxx) cvs show stroke" << endl;
			*fps << "	} for" << endl;
			*fps << "grestore" << endl;
	
			*fps << "showpage" << endl;
		}
	}
	
	ps_close(fps);
	
	return 0;
}

/**
@brief 	Generates asymmetric unit particle view and origin postscript plots.
@param 	&filename			output postscript file name.
@param 	&title				title.
@param 	&symmetry_string 	symmetry string to print at the top of the page.
@param 	*project			parameter structure.
@param 	selection			selection number (-1 selects positives, 0 selects all).
@return int					number of view panels.
**/
int 		ps_particle_views_origins(Bstring& filename, Bstring& title, 
					Bstring& symmetry_string, Bproject* project, int selection)
{
	Bsymmetry	sym(symmetry_string);
	
	View*		view = views_from_project(project, selection);
	
	change_views_to_asymmetric_unit(sym, view);
	
	if ( verbose & VERB_LABEL )
		cout << "Writing postscript file: " << filename << endl << endl;
	
	ofstream*	fps = ps_open_and_init(filename, symmetry_string, 1, 600, 800);
	
	*fps << "/Helvetica-Bold findfont 20 scalefont setfont" << endl;
	*fps << "40 740 moveto (" << filename << ") show" << endl;
	
	ps_views(fps, symmetry_string, view, 2);
	
	kill_list((char *) view, sizeof(View));	
	
	*fps << "/Helvetica-Bold findfont 20 scalefont setfont" << endl;
	*fps << "40 740 moveto (" << filename << ") show" << endl;
	
	ps_origins(fps, title, project, 2);

	ps_close(fps);

	return 0;
}

/**
@brief 	Generates phi-theta views postscript plots.
@param 	&filename			output postscript file name.
@param 	&title				title.
@param 	*project			parameter structure.
@param 	selection			selection number (-1 selects positives, 0 selects all).
@return int					number of view panels.
**/
int 		ps_particle_phi_theta(Bstring& filename, Bstring& title, 
					Bproject* project, int selection)
{
	View*		views = views_from_project(project, selection);
	
	if ( verbose & VERB_LABEL )
		cout << "Writing postscript file: " << filename << endl << endl;
	
	ofstream*	fps = ps_open_and_init(filename, title, 1, 600, 800);
	
	*fps << "/Helvetica-Bold findfont 20 scalefont setfont" << endl;
	*fps << "40 740 moveto (" << filename << ") show" << endl;
	
	ps_views(fps, views);
	
	kill_list((char *) views, sizeof(View));	
	
	*fps << "/Helvetica-Bold findfont 20 scalefont setfont" << endl;
	*fps << "40 740 moveto (" << filename << ") show" << endl;
	
	ps_close(fps);

	return 0;
}


/**
@author  Eduardo Sanz-Garcia and Bernard Heymann
@brief 	Generates postscript plot of origins on a square with the size of the box.
@param 	&filename		output postscript file name.
@param 	&title			title.
@param 	*project		micrograph particle parameter structure.
@param 	flags			flags.
@return int				0, error if <0.

	The plotting options are determined by the flags argument:
		0 = plot origins
		1 = plot numbered origins
		2 = plot origins with shading according to occurrence.
			The gray level indicates the ratio to the maximum occurrence.

**/
int 		ps_origins(Bstring& filename, Bstring& title, Bproject* project, int flags)
{
	if ( verbose & VERB_LABEL )
		cout << "Writing postscript file: " << filename << endl << endl;
	
	ofstream*	fps = ps_open_and_init(filename, title, 1, 600, 800);
	
	*fps << "/Helvetica-Bold findfont 20 scalefont setfont" << endl;
	*fps << "40 740 moveto (" << filename << ") show" << endl;
	
	ps_origins(fps, title, project, flags);

	ps_close(fps);

	return 0;
}

int			ps_origins(ofstream* fps, Bstring& title, Bproject* project, int flags)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_origins: flags=" << flags << endl;
		
	int					left(50), bottom(150), width(500), height(500);
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = project->rec;
	Bparticle*			part = NULL;
   
	//    Discover the dimension of the first particle image
	//    I supposed that the rest have the same dimentions (not always true).
	Bimage*				p = NULL;
	long		nx = 0, ny = 0, nz = 0;
	
	if ( project->select < 1 ) {
		if ( field && field->mg ) {
			nx = (long) field->mg->box_size[0];
			ny = (long) field->mg->box_size[1];
			nz = 1;
		}
		for ( field = project->field; field && !p; field = field->next )
			for ( mg = field->mg; mg && !p; mg = mg->next )
				if ( ! access(mg->fpart.c_str(), R_OK) )
					p = read_img(mg->fpart, 0, 0);   // Read one image header
	} else {
		if ( rec ) {
			nx = (long) rec->box_size[0];
			ny = (long) rec->box_size[1];
			nz = (long) rec->box_size[2];
		}
		for ( rec = project->rec; rec && !p; rec = rec->next ) {
			if ( ! access(rec->fpart.c_str(), R_OK) )
				p = read_img(rec->fpart, 0, 0);   // Read one image header
			for ( part = rec->part; part && !p; part = part->next )
				if ( ! access(part->fpart.c_str(), R_OK) )
					p = read_img(part->fpart, 0, 0);   // Read one image header
		}
	}
	
	if ( p ) {
		nx = p->sizeX();
		ny = p->sizeY();
		nz = p->sizeZ();
		delete p;
	}
	
	if ( nx*ny < 1 ) {
		error_show("Error in ps_origins: zero image size!", __FILE__, __LINE__);
		return -1;
	}
	
	long	n = 0;                   // number of particles
	float			cx = nx/2;
	float			cy = ny/2;
	float			scale = width*1.0/nx;
	if ( verbose & VERB_FULL ) {
		cout << "Dimensions of the image particles: " << nx << " * " << ny << endl;
		cout << "Center of the particle images:     " << cx << " * " << cy << endl;
	}

	long			i, m, ix, iy;
	float*			count = NULL;

	if ( flags == 2 ) {
		count = new float[nx*ny];
		for ( i=0; i<nx*ny; i++ ) count[i] = 0;
		if ( project->select < 1 ) {
			for ( n=m=0, field = project->field; field; field = field->next )
				for ( mg = field->mg; mg; mg = mg->next )
					for ( part = mg->part; part; part = part->next )
						if ( part->sel > 0 ) {
							ix = (long) part->ori[0];
							iy = (long) part->ori[1];
							if ( ix < 0 ) ix += nx;
							if ( iy < 0 ) iy += ny;
							if ( ix >= 0 && iy >= 0 && ix < nx && iy < ny ) {
								i = iy*nx+ix;
								count[i]++;
								if ( m < count[i] ) m = (int) count[i];
								n++;
							}
						}
		} else {
			for ( n=m=0, rec = project->rec; rec; rec = rec->next )
				for ( part = rec->part; part; part = part->next )
					if ( part->sel > 0 ) {
						ix = (long) part->ori[0];
						iy = (long) part->ori[1];
						if ( ix < 0 ) ix += nx;
						if ( iy < 0 ) iy += ny;
						if ( ix >= 0 && iy >= 0 && ix < nx && iy < ny ) {
							i = iy*nx+ix;
							count[i]++;
							if ( m < count[i] ) m = (int) count[i];
							n++;
						}
					}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_origins: number of particles selected = " << n << endl;
		
	*fps << "%%%%Page: Origins 2" << endl;
	*fps << "/Helvetica-Bold findfont 14 scalefont setfont" << endl;
	*fps << "/Left " << left << " def\n/Bottom " << bottom << " def\n/Height " << height << " def\n/Width " << width << " def" << endl;

	*fps << "/Origin [" << endl;
	if ( flags == 2 ) {
		*fps << "%x y count weight" << endl;
		for ( i=iy=0; iy<ny; iy++ )
			for ( ix=0; ix<nx; ix++, i++ )
				if ( count[i] ) *fps << ix << " " << iy << " " << count[i] << " " << count[i]/m << endl;
		delete[] count;
	} else {
		*fps << "%x y count weight" << endl;
		if ( project->select < 1 ) {
			for ( n=0, field = project->field; field; field = field->next )
				for ( mg = field->mg; mg; mg = mg->next )
					for ( part = mg->part; part; part = part->next )
						if ( part->sel > 0 ) {
							*fps << part->ori[0] << " " << part->ori[1] << " 0.0 0.0" << endl;
							n++;
						}
		} else {
			for ( rec = project->rec; rec; rec = rec->next )
				for ( part = rec->part; part; part = part->next )
					if ( part->sel > 0 ) {
						*fps << part->ori[0] << " " << part->ori[1] << " 0.0 0.0" << endl;
						n++;
					}
		}
	}

	*fps << "] def" << endl;

	*fps << "/Dim " << nx << " def" << endl;
	*fps << "/Rad Dim 2 div def" << endl;
	*fps << "/Scale Width Dim div def" << endl;
	*fps << "/Box { newpath " << endl;
	*fps << "	Rad Rad Rad 0 360 arc stroke" << endl;					// Circle
	*fps << "	0 0 moveto 0 Dim rlineto Dim 0 rlineto 0 Dim neg rlineto Dim neg 0 rlineto" << endl;	// Square
	*fps << "	0 Rad moveto Dim 0 rlineto" << endl;					// Cross
	*fps << "	Rad 0 moveto 0 Dim rlineto" << endl;
	*fps << "	closepath } def" << endl;

	*fps << "40 720 moveto (Number of origins:   " << n << ") show" << endl;
	if ( nz < 2 ) *fps << "40 700 moveto (Particle dimensions: " << nx << " x " << ny << ") show" << endl;
	else *fps << "40 700 moveto (Particle dimensions: " << nx << " x " << ny << " x " << nz << ") show" << endl;
	*fps << "40 680 moveto (Maximum count:       " << m << ") show" << endl;
	*fps << "/Helvetica findfont 16 scalefont setfont" << endl;
	*fps << "gsave" << endl;
	*fps << "	Left Bottom translate Scale Scale scale" << endl;
	*fps << "	" << 1/scale << " setlinewidth Box stroke" << endl;
	if ( flags == 2 ) {
		*fps << "	0 0.1 1 {" << endl;
		*fps << "		/f exch def" << endl;
		*fps << "		0 4 Origin length 4 sub {" << endl;
		*fps << "			/Index exch def" << endl;
		*fps << "			/x Origin Index get def" << endl;
		*fps << "			/y Origin Index 1 add get def" << endl;
		*fps << "			/v Origin Index 3 add get def" << endl;
		*fps << "			v f ge {" << endl;
		*fps << "				1 v sub setgray" << endl;
		*fps << "				x y 0.5 0 360 arc fill" << endl;
		*fps << "			} if" << endl;
		*fps << "		} for stroke" << endl;
		*fps << "	} for" << endl;
	} else {
		*fps << "	/n 0 def" << endl;
		*fps << "	0 4 Origin length 4 sub {" << endl;
		*fps << "		/Index exch def" << endl;
		*fps << "		/x Origin Index get def" << endl;
		*fps << "		/y Origin Index 1 add get def" << endl;
		*fps << "		x y 0.5 0 360 arc fill" << endl;
		if ( flags == 1 ) {
			*fps << "		/n n 1 add def" << endl;
			*fps << "		x y moveto 0.02 0.02 rmoveto n cvi (xxxx) cvs show stroke" << endl;
		}
		*fps << "	} for stroke" << endl;
	}
	*fps << "grestore" << endl;
	*fps << "showpage" << endl;
	
	return 0;
}

/**
@brief 	Postscript plot of the particle FOM histogram.
@param 	&filename	output postscript file name.
@param 	*project	project parameter structure.
@return int 				0, error if <0.
**/
int			ps_part_fom_histogram(Bstring& filename, Bproject* project)
{
	if ( !project ) return 0;
	
	int					f, nfom;
	long				i, j, nh(200);
	Bstring				tag;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;

	for ( f=nfom=0; f<NFOM; f++ ) if ( project->fom_tag[f] ) nfom = f + 1;

	Bstring			title("FOM histogram");
	Bplot*			plot = new Bplot(1, nh+1, nfom+1);
	plot->title(title);
	
	plot->page(0).title(title);
	plot->page(0).columns(nfom+1);
	for ( i=0; i<=nfom; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("FOM");
	plot->page(0).column(0).axis(1);
	for ( i=1; i<=nfom; i++ ) {
		tag = get_fom_tag(project->fom_tag[i-1]);
		plot->page(0).column(i).type(1);
//		plot->page(0).column(i).label("Count");
		plot->page(0).column(i).label(tag);
		plot->page(0).column(i).axis(3);
		plot->page(0).column(i).element_size(0.5);
	}
//	plot->page(0).column(nfom).label("Count");
	plot->page(0).column(1).color(0,0,0);
	if ( nfom > 1 ) plot->page(0).column(2).color(1,0,0);
	if ( nfom > 2 ) plot->page(0).column(3).color(0,1,0);
	if ( nfom > 3 ) plot->page(0).column(4).color(0,0,1);
	if ( nfom > 4 ) plot->page(0).column(5).color(1,1,0);
	plot->page(0).axis(1).min(-1);
	plot->page(0).axis(1).max(1);
	plot->page(0).axis(1).inc(0.2);
	plot->page(0).axis(1).label("FOM");
	plot->page(0).axis(3).label("Count");

/*
	for ( f=i=0; f<NFOM; f++ ) if ( project->fom_tag[f] ) {
		tag = get_fom_tag(project->fom_tag[f]);
		switch ( i ) {
			case 0: *fps << "	0 0 0 setcolor" << endl; break;
			case 1: *fps << "	1 0 0 setcolor" << endl; break;
			case 2: *fps << "	0 1 0 setcolor" << endl; break;
			case 3: *fps << "	0 0 1 setcolor" << endl; break;
			case 4: *fps << "	1 1 0 setcolor" << endl; break;
		}
		*fps << left << " " << bottom-100-20*i << " moveto (" << tag << ") show" << endl;
		i++;
	}
*/
	for ( i=0; i<=nh; i++ ) (*plot)[i] = i*2.0/nh - 1;
	
	if ( verbose & VERB_PROCESS ) {
		if ( nfom < 2 )
			cout << "Generating file " << filename << " with a FOM histogram" << endl << endl;
		else
			cout << "Generating file " << filename << " with " << nfom << " FOM histograms" << endl << endl;
	}
	
	for ( f=0, j=plot->rows(); f<nfom; f++, j+=plot->rows() ) if ( project->fom_tag[f] ) {
		if ( project->select < 1 ) {
			for ( field = project->field; field; field = field->next ) {
				for ( mg = field->mg; mg; mg = mg->next ) {
					for ( part = mg->part; part; part = part->next ) if ( part->sel ) {
						i = (long) ((part->fom[f]+1)*nh/2 + 0.5);
						if ( i < 0 ) i = 0;
						if ( i > nh ) i = nh;
						(*plot)[i+j] += 1;
					}
				}
			}
		} else {
			for ( rec = project->rec; rec; rec = rec->next ) {
				for ( part = rec->part; part; part = part->next ) if ( part->sel ) {
					i = (long) ((part->fom[f]+1)*nh/2 + 0.5);
					if ( i < 0 ) i = 0;
					if ( i > nh ) i = nh;
					(*plot)[i+j] += 1;
				}
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_part_fom_histogram: Histogram data transferred" << endl;
	
	ofstream*		fps = ps_open_and_init(filename, plot);
	
	ps_graph(fps, plot, 1);
	
	delete plot;


	long		npart(0);
	if ( project->select < 1 ) npart = project_count_mg_part_selected(project);
	else npart = project_count_rec_part_selected(project);

	if ( nfom > 1 && npart > 0) {
		title = "FOM plot";
		plot = new Bplot(1, npart, nfom);
		plot->title(title);
//		cout << "DEBUG ps_part_fom_histogram: page=" << plot->page(0).number() << endl;
	
		plot->page(0).title(title);
		plot->page(0).columns(nfom);
		for ( i=0; i<nfom; i++ ) plot->page(0).column(i).number(i);
		tag = get_fom_tag(project->fom_tag[0]);
		plot->page(0).column(0).label(tag);
		plot->page(0).column(0).axis(1);
		for ( i=1; i<nfom; i++ ) {
			tag = get_fom_tag(project->fom_tag[i]);
			plot->page(0).column(i).type(0);
			plot->page(0).column(i).label(tag);
			plot->page(0).column(i).axis(3);
			plot->page(0).column(i).element_size(0.9);
		}
		
		plot->page(0).column(1).color(1,0,0);
		if ( nfom > 2 ) plot->page(0).column(2).color(0,1,0);
		if ( nfom > 3 ) plot->page(0).column(3).color(0,0,1);
		if ( nfom > 4 ) plot->page(0).column(4).color(1,1,0);
		
		for ( i=1; i<4; i+=2 ) {
			plot->page(0).axis(i).min(-1);
			plot->page(0).axis(i).max(1);
			plot->page(0).axis(i).inc(0.2);
		}
	
		for ( f=0, j=0; f<nfom; f++, j+=npart ) if ( project->fom_tag[f] ) {
			if ( project->select < 1 ) {
				for ( i=0, field = project->field; field; field = field->next ) {
					for ( mg = field->mg; mg; mg = mg->next ) {
						for ( part = mg->part; part; part = part->next ) if ( part->sel ) {
							(*plot)[i+j] = part->fom[f];
							i++;
						}
					}
				}
			} else {
				for ( i=0, rec = project->rec; rec; rec = rec->next ) {
					for ( part = rec->part; part; part = part->next ) if ( part->sel ) {
						(*plot)[i+j] = part->fom[f];
						i++;
					}
				}
			}
		}

		if ( verbose & VERB_DEBUG )
			cout << "DEBUG ps_part_fom_histogram: FOM data: npart=" << npart << " transferred=" << i << endl;

		ps_graph(fps, plot, 1);
		
		delete plot;
	}

	ps_class_average_fom(fps, project);
	
	ps_close(fps);	
	
	return 0;
}

/**
@brief 	Postscript histogram plot of micrograph defocus parameters.
@param 	&filename	output postscript file name.
@param 	*project	project parameter structure.
@return int 				0, error if <0.
**/
int			ps_defocus_histogram(Bstring& filename, Bproject* project)
{
	if ( !project ) return 0;
	
	int					n(0);
	long				i, nc(2), nh(100);
	double				damax(10), da, dd, a, c, s, ca(0), sa(0), r2(0);
	double				daa(0), das(0), dda(0), dds(0);
	Bfield*				field;
	Bmicrograph*		mg;
//	Breconstruction*	rec;

	Bstring			title("Defocus histogram");
	Bplot*			plot = new Bplot(1, nh+1, nc);
	plot->title(title);

	plot->page(0).title(title);
	plot->page(0).columns(nc);
	for ( i=0; i<nc; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Defocus");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(1);
	plot->page(0).column(1).label("Count");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).element_size(0.5);
	plot->page(0).column(1).color(0,0,0);
	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(damax);
	plot->page(0).axis(1).inc(2);
	plot->page(0).axis(1).label("Defocus");
	plot->page(0).axis(3).label("Count");

	for ( i=0; i<=nh; i++ ) (*plot)[i] = i*damax/nh;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->ctf ) {
			da = mg->ctf->defocus_average() * 1e-4;
			dd = mg->ctf->defocus_deviation() * 1e-4;
			a = mg->ctf->astigmatism_angle();
			c = cos(a);
			s = sin(a);
			ca += c;
			sa += s;
//			r2 += c*c + s*s;
			daa += da;
			dda += dd;
			das += da*da;
			dds += dd*dd;
			n++;
			i = (long) (nh*da/damax + 0.5);
			if ( i < 0 ) i = 0;
			if ( i > nh ) i = nh;
			(*plot)[i+nh] += 1;
		}
	}
	
	daa /= n;
	dda /= n;
	das = das/n - daa*daa;
	dds = dds/n - dda*dda;
	if ( das > 0 ) das = sqrt(das);
	else das = 0;
	if ( dds > 0 ) dds = sqrt(dds);
	else dds = 0;
	
	ca /= n;
	sa /= n;
	a = atan2(sa, ca);
	r2 = ca*ca + sa*sa;
	
	if ( verbose ) {
		cout << "Defocus average:                " << daa << " (" << das << ")" << endl;
		cout << "Defocus deviation:              " << dda << " (" << dds << ")" << endl;
		cout << "Astigmatism angle:              " << a*180/M_PI << " (" << sqrt(-log(r2)) << ")" << endl;
	}

	if ( verbose & VERB_PROCESS )
		cout << "Generating file " << filename << " with defocus histograms" << endl << endl;

	ofstream*		fps = ps_open_and_init(filename, plot);
	
	ps_graph(fps, plot, 1);
	
	delete plot;


	ps_close(fps);	
	
	return 0;
}

/**
@brief 	Postscript plot of micrograph astigmatism parameters.
@param 	&filename	output postscript file name.
@param 	*project	project parameter structure.
@return int 		0, error if <0.
**/
int			ps_astigmatism_plot(Bstring& filename, Bproject* project)
{
	if ( !project ) return 0;
	
	long				i, j, nc(2);
	Bfield*				field;
	Bmicrograph*		mg;

	long				nmg = project_count_micrographs(project);
	
	Bstring				title("Astigmatism plot");
	Bplot*				plot = new Bplot(1, nmg, nc);
	plot->title(title);

	plot->page(0).title(title);
	plot->page(0).columns(nc);
	for ( i=0; i<nc; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Defocus Deviation");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(0);
	plot->page(0).column(1).label("Astigmatism Angle");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).element_size(2);
	plot->page(0).column(1).color(0,0,0);
//	plot->page(0).axis(1).min(0);
//	plot->page(0).axis(1).max(1e4);
//	plot->page(0).axis(1).inc(1e3);
	plot->page(0).axis(1).label("Defocus Deviation");
//	plot->page(0).axis(3).min(-90);
//	plot->page(0).axis(3).max(90);
//	plot->page(0).axis(3).inc(15);
	plot->page(0).axis(3).label("Astigmatism Angle");

	for ( i=0, j=nmg, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->ctf ) {
			(*plot)[i++] = mg->ctf->defocus_deviation();
			(*plot)[j++] = mg->ctf->astigmatism_angle()*180.0/M_PI;
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Generating file " << filename << " with an astigmatism plot" << endl << endl;

	ofstream*		fps = ps_open_and_init(filename, plot);
	
	ps_graph(fps, plot, 1);
	
	delete plot;

	ps_close(fps);	
	
	return 0;
}


/**
@brief 	Postscript histogram plot of class average FOM's.
@param 	*fps		output stream.
@param 	*project	project parameter structure.
@return int 		0, error if <0.
**/
int			ps_class_average_fom(ofstream* fps, Bproject* project)
{
	if ( !project ) return 0;
	if ( !project->class_avg ) return 0;
	
	long			i, nc(0), nmax(0), nt(0);
	double			fommax(0);
	Bparticle*		ca = project->class_avg;
	
	for ( ca = project->class_avg; ca; ca = ca->next ) if ( ca->sel ) {
		if ( nmax < ca->sel ) nmax = ca->sel;
		if ( fommax < ca->fom[0] ) fommax = ca->fom[0];
		nt += ca->sel;
		nc++;
	}
	
	nmax = log(nmax) + 1;

	if ( verbose & VERB_PROCESS )
		cout << "Generating a class average plot" << endl << endl;
	
	Bstring			title("Class average FOMs");
	Bplot*			plot = new Bplot(1, nc, 2);
	plot->title(title);

	plot->page(0).title(title);
	plot->page(0).columns(2);
	for ( i=0; i<2; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("log(number)");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(0);
	plot->page(0).column(1).label("FOM");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).element_size(5);
	plot->page(0).column(1).color(0,0,0);
	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(nmax);
	plot->page(0).axis(1).flags(2);
//	plot->page(0).axis(1).inc(2);
	plot->page(0).axis(1).label("log(number)");
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(fommax);
	plot->page(0).axis(3).label("FOM");

	for ( i=0, ca = project->class_avg; ca; ca = ca->next ) if ( ca->sel ) {
		(*plot)[i] = log(ca->sel);
		(*plot)[i+nc] = ca->fom[0];
		i++;
	}

	Bstring			txt;
	txt = Bstring(nc, "Number of class averages: %ld");
	plot->page(0).add_text(txt);
	txt = Bstring(nt, "Number of particles: %ld");
	plot->page(0).add_text(txt);
			
	ps_graph(fps, plot, 1);
	
	delete plot;

	return 0;
}
