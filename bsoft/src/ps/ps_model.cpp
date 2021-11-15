/**
@file	ps_model.cpp
@brief	Postscript output for models
@author Bernard Heymann
@date	Created: 20090203
@date	Modified: 20190201
**/

#include "ps_model.h"
#include "ps_plot.h"
#include "ps_views.h"
#include "model_views.h"
#include "model_compare.h"
#include "model_util.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Postscript plot of the distribution of model views.
@param 	&filename		output postscript file name.
@param 	*model			model parameter structure.
@param 	combined		flag to show all models combined.
@return int 			0, error if <0.
**/
int			ps_model_views(Bstring& filename, Bmodel* model, int combined)
{
	int 			i(1), left(50), bottom(200), width(500), height(250), nm(1);
	double			x, y;
	Euler			euler;
	
	Bmodel*			mp;
	Bcomponent*		comp = model->comp;
	if ( !comp ) return -1;

	if ( !combined )
		for ( nm=0, mp = model; mp; mp = mp->next ) nm++;

	Bstring			title("Model views");
	
	ofstream*		fps = ps_open_and_init(filename, title, nm, 600, 800);
	
	if ( combined ) {
		*fps << "%%Page: " << i << " " << i << endl;
		*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
		*fps << "50 755 moveto (" << title << ": " << model->identifier() << ") show" << endl;
		*fps << "/Data [" << endl << "%x y fom sel" << endl;
		for ( mp = model; mp; mp = mp->next ) {
			for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
				View	v(comp->view()[0],comp->view()[1],comp->view()[2],comp->view()[3]);
				euler = Euler(v);
				x = euler.phi()*fabs(sin(euler.theta())) + M_PI;
				y = euler.theta();
				*fps << setprecision(4) << x << " " << y << " " << comp->FOM() << " " << comp->select() << endl;
			}
		}
		*fps << "] def" << endl;
		ps_phi_theta_plot(fps, left, bottom, width, height, 4);
		*fps << "showpage" << endl;
	} else {
		for ( i=1, mp = model; mp; mp = mp->next, i++ ) {
			*fps << "%%Page: " << i << " " << i << endl;
			*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
			*fps << "50 755 moveto (" << title << ": " << mp->identifier() << ") show" << endl;
			*fps << "/Data [" << endl << "%x y fom sel" << endl;
			for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
				View	v(comp->view()[0],comp->view()[1],comp->view()[2],comp->view()[3]);
				euler = Euler(v);
				x = euler.phi()*fabs(sin(euler.theta())) + M_PI;
				y = euler.theta();
				*fps << setprecision(4) << x << " " << y << " " << comp->FOM() << " " << comp->select() << endl;
			}
			*fps << "] def" << endl;
			ps_phi_theta_plot(fps, left, bottom, width, height, 4);
			*fps << "showpage" << endl;
		}
	}

	ps_close(fps);
	
	return 0;
}

/**
@brief 	Postscript plot of the distribution of model views.
@param 	&filename			output postscript file name.
@param 	*model				model parameter structure.
@param 	&symmetry_string 	symmetry.
@param 	combined			flag to show all models combined.
@return int 				0, error if <0.
**/
int			ps_model_symmetry_views(Bstring& filename, Bmodel* model, string& symmetry_string, int combined)
{
	int 			i(1), nm(1);
	Bmodel*			mp;
	
	if ( !model->comp ) return -1;

	if ( !combined )
		for ( nm=0, mp = model; mp; mp = mp->next ) nm++;

	Bstring			title("Model views");
	
	ofstream*		fps = ps_open_and_init(filename, title, nm, 600, 800);
	
	list<View2<float>>	view;
	
	if ( combined ) {
		*fps << "%%Page: " << i << " " << i << endl;
		*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
		*fps << "50 755 moveto (" << title << ": " << model->identifier() << ") show" << endl;
		*fps << "/Data [" << endl << "%x y fom sel" << endl;
		view = views_from_models(model);
		ps_views2(fps, symmetry_string, view, 2);
//		kill_list((char *) view, sizeof(View));
		*fps << "showpage" << endl;
	} else {
		for ( i=1, mp = model; mp; mp = mp->next, i++ ) {
			*fps << "%%Page: " << i << " " << i << endl;
			*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
			*fps << "50 755 moveto (" << title << ": " << mp->identifier() << ") show" << endl;
			*fps << "/Data [" << endl << "%x y fom sel" << endl;
			view = views_from_model(model);
			ps_views2(fps, symmetry_string, view, 2);
//			kill_list((char *) view, sizeof(View));
			*fps << "showpage" << endl;
		}
	}

	ps_close(fps);
	
	return 0;
}

/**
@brief 	Postscript plot of the model component FOM histogram.
@param 	&filename	output postscript file name.
@param 	*model		model parameter structure.
@return int 				0, error if <0.
**/
int			ps_model_fom_histogram(Bstring& filename, Bmodel* model)
{
	if ( !model ) return 0;
	
	int				nmod, i, j, nh(100), max(0);
	Bmodel* 		mp;
	Bcomponent*		comp;
	
	for ( nmod=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) nmod++;

	if ( verbose & VERB_PROCESS ) {
		if ( nmod < 2 )
			cout << "Generating file " << filename << " with a FOM histogram" << endl << endl;
		else
			cout << "Generating file " << filename << " with " << nmod << " FOM histograms" << endl << endl;
	}

	Bstring			title("FOM histogram");
	Bplot*			plot = new Bplot(nmod, nh+1, nmod+1);
	plot->title(title);

	for ( i=0; i<=nh; i++ ) (*plot)[i] = i*1.0/nh;

	for ( j=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		j += plot->rows();
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			i = (long) (comp->FOM()*nh + 0.5);
			if ( i < 0 ) i = 0;
			if ( i > nh ) i = nh;
			(*plot)[j+i] += 1;
		}
	}

	for ( i=plot->rows(); i<plot->rows()*plot->columns(); i++ )
		if ( max < (*plot)[i] ) max = (int) (*plot)[i];
	max = 10*(max/10 + 1);
	
	for ( i=0; i<nmod; i++ ) {
		plot->page(i).title(title);
		plot->page(i).columns(2);
		plot->page(i).column(0).number(0);
		plot->page(i).column(1).number(i+1);
		plot->page(i).column(0).label("FOM");
		plot->page(i).column(1).label("Count");
		plot->page(i).column(0).axis(1);
		plot->page(i).column(1).axis(3);
		plot->page(i).column(1).type(1);
		plot->page(i).column(1).element_size(0.5);
		plot->page(i).axis(1).max(1);
		plot->page(i).axis(3).max(max);
		plot->page(i).axis(1).label("FOM");
		plot->page(i).axis(3).label("Count");
	}
	
	ps_plot(filename, plot);
	
	delete plot;
	
	return 0;
}

/**
@brief 	Generates a postscript plot for a model occupancy analysis.
@author	Daniel Nemecek
@author	Bernard Heymann
@param 	*model			model.
@param 	cutoff			coverage cutoff (if 0 don't show).
@param 	bins			number of bins for histograms.
@param 	nfit			number of binomial distributions to fit.
@param 	*distrib		distribution array ((3+nfit)*(ncomp+1)).
@param 	*prob			weight and probability array (2*nfit).
@param 	R				binomial fit residual.
@param 	&filename		postscript file name.
@return int				0, <0 on error.

	The distribution and probability arrays must already be allocated
	and the content calculated with model_occupancy_distribution.

**/
int			ps_model_occupancy(Bmodel* model, double cutoff, int bins, int nfit, 
				vector<double>& distrib, vector<double>& prob, double R, Bstring& filename)
{
	if ( filename.length() < 1 ) return -1;
	
	Bmodel*			mp;
	Bcomponent*		comp;

	long			ncomp = model_maxnum_components(model);
	
	long			i, j;	
	long			pages(3);
	long			width(600);
	long			height(800);
	long			npoint(ncomp+1); 
	long			nfunc(3+nfit);
//	double			x_min(0), x_max(ncomp);
//	double			y_min(0), y_max(0);
	double			y_max(0);
	double			hmin(1e6), hmax(-1e6);
	double			scale, shift;
	Bstring			title("Occupancy Analysis");
	Bstring			label, txt;
	
	// Page 1
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_model_occupancy: first page " << endl;

	Bplot*			plot = new Bplot(1, npoint, nfunc+1);
	plot->title(title);

	title = filename + ": Occupancy";
	
	plot->page(0).title(title);
	plot->page(0).columns(nfunc+1);
	for ( i=0; i<=nfunc; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Occupancy");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).label("Count");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).type(1);
	plot->page(0).column(2).label("StDev");
	plot->page(0).column(2).axis(3);
	plot->page(0).column(2).type(3);
	for ( i=3; i<=nfunc; i++ ) {
		label = Bstring(i-3, "Fit%d");
		plot->page(0).column(i).type(2);
		plot->page(0).column(i).label(label);
		plot->page(0).column(i).axis(3);
		plot->page(0).column(i).element_size(0.5);
	}
	plot->page(0).column(3).label("FitSum");
	plot->page(0).column(4).color(1,0,0);
	if ( nfit > 1 ) plot->page(0).column(5).color(0,1,0);
	if ( nfit > 2 ) plot->page(0).column(6).color(0,0,1);
	if ( nfit > 3 ) plot->page(0).column(7).color(1,1,0);
	
	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(ncomp);
	plot->page(0).axis(1).label("Occupancy");
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(1);
	plot->page(0).axis(3).label("Models (fraction)");
	
	if ( cutoff > 0 ) {
		txt = Bstring(cutoff, "Coverage cutoff: %g");
		plot->page(0).add_text(txt);
	}
	
	if ( nfit > 1 ) {
		txt = "Weights:     ";
		for ( i=0; i<nfit; i++ ) txt += Bstring(prob[i], " %f");
		plot->page(0).add_text(txt);
		txt = "Probability: ";
		for ( i=0; i<nfit; i++ ) txt += Bstring(prob[i+nfit], " %f");
		plot->page(0).add_text(txt);
	}

	txt = Bstring(R, "R = %g");
	plot->page(0).add_text(txt);

	for ( i=0; i<npoint; i++) (*plot)[i] = i;
	for ( i=0, j=npoint; i<npoint*nfunc; i++, j++ ) (*plot)[j] = distrib[i];
	
	ofstream*		fps = ps_open_and_init(filename, title, pages, width, height);

	ps_graph(fps, plot, 1);
	
	delete plot;

	// Page 2
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_model_occupancy: second page " << endl;

	for ( hmin=1e6, hmax=-1e6, mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) {
			if ( comp->density() > hmax ) hmax = comp->density();
			if ( comp->density() < hmin ) hmin = comp->density();
		}
	}
	
	if ( hmin == hmax ) {
		cerr << "Warning: The component densities vary only from " << hmin << " to " << hmax << endl;
		hmax = hmin + 1;
	}
	
	scale = (bins - 1)/(hmax - hmin);
	
	plot = new Bplot(1, bins, 2);

	title = filename + ": Average density";
	
	plot->page(0).title(title);
	plot->page(0).columns(2);
	for ( i=0; i<2; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Density");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).label("Components");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).type(1);
	plot->page(0).column(1).element_size(0.5);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_model_occupancy: bins=" << bins << " scale=" << scale << endl;

	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) {
			i = (long) (scale*(comp->density() - hmin) + 0.5);
			if ( i < bins ) (*plot)[i+bins] += 1;
		}
	}
	
	for ( i=0, j=bins, y_max = 0; i<bins; i++, j++ ) {
		(*plot)[i] = i/scale + hmin;
		if ( y_max < (*plot)[j] ) y_max = (*plot)[j];
	}
	y_max *= 1.05;
	
	plot->page(0).axis(1).min(hmin);
	plot->page(0).axis(1).max(hmax);
	plot->page(0).axis(1).label("Density");
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(y_max);
	plot->page(0).axis(3).label("Components");
	
	ps_graph(fps, plot, 1);
	
	delete plot;

	// Page 3
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_model_occupancy: third page " << endl;

	for ( hmin=1e6, hmax=-1e6, mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next, i++ ) {
			shift = comp->location().distance(comp->force());
			if ( shift > hmax ) hmax = shift;
			if ( shift < hmin ) hmin = shift;
		}
	}
	
	scale = (bins - 1)/(hmax - hmin);
	
	plot = new Bplot(1, bins, 2);

	title = filename + ": Shifts";
	
	plot->page(0).title(title);
	plot->page(0).columns(2);
	for ( i=0; i<2; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Density");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).label("Components");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).type(1);
	plot->page(0).column(1).element_size(0.5);
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) {
			if ( comp->select() ) { // count shifts of occupied locations only
				shift = comp->location().distance(comp->force());
				i = (long) (scale*(shift - hmin) + 0.5);
				if ( i < bins ) (*plot)[i+bins] += 1;
			}
		}
	}
	
	for ( i=0, j=bins, y_max = 0; i<bins; i++, j++ ) {
		(*plot)[i] = i/scale + hmin;
		if ( y_max < (*plot)[j] ) y_max = (*plot)[j];
	}
	y_max *= 1.05;
	
	plot->page(0).axis(1).min(hmin);
	plot->page(0).axis(1).max(hmax);
	plot->page(0).axis(1).label("Shift (A)");
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(y_max);
	plot->page(0).axis(3).label("Components");
	
	ps_graph(fps, plot, 1);
	
	delete plot;

	ps_close(fps);
	
	return 0;
}

