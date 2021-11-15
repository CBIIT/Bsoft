/**
@file	bmatch.cpp
@brief	Program to search an image with multiple projections and orient them.
@author Bernard Heymann
@date	Created: 20170720
@date	Modified: 20190821

clang++ -o bin/bmatch src/bmatch.cpp -I. -I/usr/local/include -I/Users/bernard/b20/bsoft/include -L/Users/bernard/b20/bsoft/lib -lbsoft -DHAVE_FFTW3 -std=c++11 -I/Users/bernard/b20/fftw-3.3.6-pl2/include
**/

#include "rwimg.h"
#include "rwmodel.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bmodel*		img_match(Bimage* p, Bimage* pref, double res_hi, double res_lo, 
				double angle_inc, double threshold, long exclusion, double mask_rad);
Bmodel*		img_match_one(Bimage* p, long img_num, Bimage* pref, double res_hi, double res_lo, 
				double angle_inc, double threshold, long exclusion, double mask_rad, fft_plan planf, fft_plan planb);

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmatch [options] img.mrc cc.map",
"--------------------------------------",
"Searching an image with a set of 2D templates/references",
"Orientations are taken from the input projections.",
" ",
//"Actions:",
//"-variance                Calculate shell variances of input projections.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 2,3.5,1        Sampling before rescaling (angstrom/voxel, a single value sets all three).",
"-resolution 25,200       Resolution limits (angstrom).",
"-angle 2.5               Angular increment (default 3 degrees).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-threshold 0.24          Threshold to find correlation peaks.",
"-exclusion 150           Minimum distance between correlation peaks.",
"-radius 50               Mask radius (slightly bigger than particle).",
" ",
"Input:",
"-reference file.pif      Input reference projections file.",
" ",
"Output:",
"-output newmod.star      New model file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
    
	Vector3<double> origin;						// Origin
	int				set_origin(0);				// Flag for setting the origin
	Vector3<double>	sampling;					// Sampling before rescaling
	double			res_hi(0);					// High resolution limit
	double			res_lo(0); 					// Low resolution limit
	double 			angle_inc(M_PI/60);			// Angular increment
	double			threshold(0.1);				// Threshold for peak detection
	long			exclusion(10);				// Peak detection exclusion size
	double			mask_rad(70);				// Particle masking radius
	Bstring			ref_filename;				// Input projections file name
	Bstring			outfile;					// Output parameter file name
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(res_hi, res_lo) < 1 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
 		if ( curropt->tag == "angle" )
			if ( ( angle_inc = curropt->value.real() ) < 0.1 )
				cerr << "-angle: An angle must be specified!" << endl;
		if ( curropt->tag == "sampling" )
			sampling = curropt->scale();
 		if ( curropt->tag == "threshold" )
			if ( ( threshold = curropt->value.real() ) < 0.001 )
				cerr << "-threshold: A resolution limit must be specified!" << endl;
 		if ( curropt->tag == "exclusion" )
			if ( ( exclusion = curropt->value.integer() ) < 2 )
				cerr << "-exclusion: A minimum distance between peaks must be specified!" << endl;
 		if ( curropt->tag == "radius" )
			if ( ( mask_rad = curropt->value.real() ) < 0.001 )
				cerr << "-radius: A mask radius in pixels must be specified!" << endl;
		if ( curropt->tag == "reference" )
			ref_filename = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);

	double			ti = timer_start();
	
    Bimage*	 		p = NULL;
    Bimage*	 		pref = NULL;
    
	if ( ( p = read_img(argv[optind++], 1, -1) ) == NULL ) {
		cerr << "Error: File not read!" << endl;
		bexit(-1);
	}

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}
	
	if ( ref_filename.length() ) {
		if ( ( pref = read_img(ref_filename, 1, -1) ) == NULL ) {
			cerr << "Error: File not read!" << endl;
			bexit(-1);
		}
	}

	if ( set_origin ) {
		if ( set_origin == 2 ) pref->origin(p->size()/2);
		else pref->origin(origin);
	}
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();
	
	if ( sampling.volume() > 0 ) {
		p->sampling(sampling);
		if ( pref ) pref->sampling(sampling);
	}

	Bmodel*		model = NULL;
	if ( pref )
		model = img_match(p, pref, res_hi, res_lo, angle_inc, threshold, exclusion, mask_rad);

	// Write an output parameter format file if a name is given
    if ( outfile.length() && model ) {
		write_model(outfile, model);
	}

	model_kill(model);
	
	delete p;
	delete pref;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	return 0;
}

Bmodel*		img_match(Bimage* p, Bimage* pref, double res_hi, double res_lo,
				double angle_inc, double threshold, long exclusion, double mask_rad)
{
	Bmodel*			modlist = NULL;
	Bmodel*			model;
	
	fft_plan		planf = fft_setup_plan(p->size(), FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(p->size(), FFTW_BACKWARD, 1);
	
	for ( long n = 0; n < p->images(); ++n ) {
		model = img_match_one(p, n, pref, res_hi, res_lo, 
			angle_inc, threshold, exclusion, mask_rad, planf, planb);
//		model_add(&modlist, model);
		if ( modlist ) modlist->add(model);
		else modlist = model;
	}
	
    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	return modlist;
}

Bimage*		img_reference_prepare(Bimage* pref, long nn, Vector3<long> size)
{
	Bimage*		pref1 = pref->extract(nn);
	
	if ( pref1->size() != size ) {
		pref1->calculate_background();
		pref1->pad(size, FILL_BACKGROUND);
	}
	
//	write_img("t.pif", pref1, 0);
	
	return pref1;
}


/*
*/
Bmodel*		img_match_one(Bimage* p, long img_num, Bimage* pref, double res_hi, double res_lo, 
				double angle_inc, double threshold, long exclusion, double mask_rad, fft_plan planf, fft_plan planb)
{
	if ( angle_inc < M_PI/180 ) angle_inc = M_PI/180;
	
	long			i, j, k, nn;
	double			angle, amin(0), amax(TWOPI), ccmax(0), ccmax1;
	Vector3<long>	coor, rect(2*mask_rad, 2*mask_rad, 1);
	Vector3<double>	start, nucoor;
	View			view;

	Bstring			id(img_num, "%d");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;
	Bcomptype*		ct = NULL;
	
	model->mapfile(p->file_name());
	model->image_number(img_num);

	Bimage*			p1 = p->extract(img_num);

	p1->calculate_background();

	Bimage*			pcc = new Bimage(Float, TSimple, p1->size(), 1);
	Bimage*			pang = pcc->copy();
	Bimage*			pnn = new Bimage(Integer, TSimple, p1->size(), 1);
	pnn->fill(-1);

	if ( verbose ) {
		cout << "Global search for image " << img_num << ":" << endl;
		cout << "Proj\tvx\tvy\tvz\tva\tccmax" << setprecision(4) << endl;
	}
	for ( nn=0; nn<pref->images(); ++nn ) {
		id = Bstring(nn+1, "%d");
//		model_add_type_by_id(model, id);
		model->add_type(id);
		Bimage*		pref1 = img_reference_prepare(pref, nn, p1->size());
		ccmax1 = 0;
		for ( angle = amin; angle <= amax; angle += angle_inc ) {
			Bimage*		prot = pref1->rotate(pref1->size(), angle);
			Bimage*		pcc1 = p1->cross_correlate(prot, res_hi, res_lo, planf, planb);
			if ( !pcc1 ) bexit(-1);
			for ( j=0; j<pcc1->image_size(); ++j ) if ( (*pcc)[j] < (*pcc1)[j] ) {
				pcc->set(j, (*pcc1)[j]);
				pang->set(j, angle);
				pnn->set(j, nn);
			}
			if ( ccmax1 < pcc1->maximum() ) ccmax1 = pcc1->maximum();
			delete pcc1;
			delete prot;
		}
		if ( verbose )
			cout << nn << tab << pref->image[nn].view() << tab << ccmax1 << endl;
		if ( ccmax < ccmax1 ) ccmax = ccmax1;
		delete pref1;
	}
	if ( verbose )
		cout << endl;

	if ( threshold >= ccmax ) {
		threshold = 0.95*ccmax;
		if ( verbose )
			cout << "Adjusting threshold to " << threshold << endl;
	}

	Bimage*			ppeak = pcc->region_peaks(exclusion, threshold, 0, 1);
	write_img("cc.pif", pcc, 0);
	write_img("nn.pif", pnn, 0);
	
	pcc->clear();
	
	for ( i=j=0; i<ppeak->image_size(); ++i ) if ( (*ppeak)[i] ) j++;
	
	if ( verbose )
		cout << "Number of peaks:          " << j << endl << endl;

	threshold *= p1->sizeX()*1.0/pref->sizeX();
	
	if ( verbose ) {
		cout << "Adjusting threshold to " << threshold << endl;
		cout << "Local search for image " << img_num << ":" << endl;
		cout << "x\ty\tz\tvx\tvy\tvz\tva\tccmax" << setprecision(4) << endl;
	}
	rect = Vector3<long>(2*mask_rad, 2*mask_rad, 1);
	angle_inc /= 10;
	for ( i=j=0; i<ppeak->image_size(); ++i ) if ( (*ppeak)[i] ) {
		nn = (*pnn)[i];
		for ( k=0, ct = model->type; k<nn; ++k ) ct = ct->next;
		Bimage*		pref1 = pref->extract(nn);
		pref1->calculate_background();
		coor = vector3_set_PBC(pref1->image->origin() + ppeak->coordinates(i), p->size());
		start[0] = coor[0] - mask_rad;
		start[1] = coor[1] - mask_rad;
		Bimage*		pe = p1->copy();
		pe->calculate_background();
		pe->origin(pref1->image->origin());
		pe->edge(0, 1, rect, start, 5, FILL_BACKGROUND);
//			Bstring		fn(j+1, "t%d.pif");
//			write_img(fn, pe, 0);
		ccmax1 = 0;
		amin = (*pang)[i] - 10*angle_inc;
		amax = (*pang)[i] + 10*angle_inc;
		coor = vector3_set_PBC(ppeak->coordinates(i), p->size());
		nucoor = coor;
		for ( angle = amin; angle <= amax; angle += angle_inc ) {
			Bimage*		prot = pref1->rotate(pref1->size(), angle);
			Bimage*		pcc1 = pe->cross_correlate(prot, res_hi, res_lo);
			if ( !pcc1 ) bexit(-1);
			pcc1->origin(coor);
			pcc1->find_peak(5);
			if ( (*pcc)[i] < pcc1->image->FOM() ) {
				pcc->set(i, pcc1->image->FOM());
				pang->set(i, angle);
				nucoor = pcc1->image->origin();
			}
			if ( ccmax1 < pcc1->image->FOM() ) ccmax1 = pcc1->image->FOM();
			delete pcc1;
			delete prot;
		}
		delete pe;
		delete pref1;
		view = pref->image[nn].view();
		view.angle(view.angle() - (*pang)[i]);
		if ( p->sizeX() > 1 && nucoor[0] >= p->sizeX()/2 ) nucoor[0] -= p->sizeX();
		if ( p->sizeY() > 1 && nucoor[1] >= p->sizeY()/2 ) nucoor[1] -= p->sizeY();
		if ( p->sizeZ() > 1 && nucoor[2] >= p->sizeZ()/2 ) nucoor[2] -= p->sizeZ();
		if ( verbose )
			cout << nucoor << tab << view << tab << ccmax1 << endl;
		if ( ccmax1 >= threshold ) {
//			comp = component_add(&model->comp, ++j);
			if ( comp ) comp = comp->add(++j);
			else model->comp = comp = new Bcomponent(++j);
			comp->type(ct);
			comp->location(nucoor*p->sampling(0));
			view = view.backward();
			comp->view(View2<float>(view[0],view[1],view[2],view[3]));
			comp->FOM(ccmax1);
			comp->radius(10*p->sampling(0)[0]);
		}
	}
	if ( verbose )
		cout << "Number of peaks:          " << j << endl << endl;

	delete ppeak;
	delete pcc;
	delete pang;
	delete pnn;
	delete p1;

	return model;
}


