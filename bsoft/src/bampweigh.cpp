/**
@file	bampweigh.cpp
@brief	Program to filter images.
@author Bernard Heymann
@date	Created: 20040714
@date	Modified: 20191125
**/

#include "rwimg.h"
#include "rwFSC_XML.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bampweigh [options] input.img output.img",
"-----------------------------------------------",
"Weighs the amplitudes of a map using a reference map's amplitudes.",
" ",
"Actions:",
"-invert                  Invert density in the image.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation after filtering.",
"-dose 2.5                Weigh the amplitudes by accumulated dose using the given dose/frame (e/Å2).",
"-normalize 2             Normalize amplitudes (1:amplitude; 2:power).",
"-root                    Weigh by the square root of the amplitudes.",
"-square                  Weigh by the squares of the amplitudes.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; default from input file; a single value can be given).",
"-resolution 15           Resolution limit (angstrom)",
" ",
"Input:",
"-RPS file.txt            Radial power spectrum.",
"-FSC file.xml            FSC curve.",
"-Reference file.map      Reference map.",
"-Mask mask.tif           Mask file to exclude regions in reciprocal space.",
"-mask mask.mrc           Real space mask to apply before weighing.",
" ",
NULL
};

Bplot*		read_rps(Bstring& filename)
{
	if ( verbose & VERB_PROCESS )
		cout << "# Reading radial power spectrum:  " << filename << endl;

    ifstream			f;
    f.open(filename.c_str());
    if ( f.fail() ) return NULL;
	
	int					i, j, ncol(2), nrow(0);
	char				aline[MAXLINELEN];
    double				s, p;
    vector<double>		vs, vp;

	f.getline(aline, MAXLINELEN);

 	while ( !f.eof() ) {
 		f >> s >> p;
 		if ( p ) {
 			vs.push_back(s);
 			vp.push_back(p);
			nrow++;
		}
	}
	
	f.close();

	Bstring				title = "RPS Plot";
	Bplot*				plot = new Bplot(1, nrow, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Spatial Frequency (A)");
	plot->page(0).column(1).label("Power");
//	plot->page(0).axis(1).min(0);
//	plot->page(0).axis(1).max(1/hi_res);
//	plot->page(0).axis(1).inc(0.1/hi_res);
	plot->page(0).axis(1).flags(1);
	plot->page(0).axis(1).label("Spatial Frequency (1/Å)");
//	plot->page(0).axis(3).min(0);
//	plot->page(0).axis(3).max(1);
//	plot->page(0).axis(3).inc(0.1);
	plot->page(0).column(1).type(2);
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).axis(3);

 	for ( i=0, j=nrow; i<nrow; ++i, ++j ) {
		(*plot)[i] = vs[i];
		(*plot)[j] = vp[i];
		if ( verbose & VERB_DEBUG )
			cout << vs[i] << tab << vp[i] << endl;
	}
	
	return plot;
}

int 	main(int argc, char **argv)
{
	// Initialize variables
	int 			setinvert(0);				// Flag to invert density
	double			weigh_dose(0);				// Dose per frame to weigh by dose
	int				normalize(0);				// Flag to normalize amplitudes
	int				root(0);					// Flag to weigh by sqrt of amplitudes
	int				square(0);					// Flag to weigh by squares of amplitudes
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	double			resolution(0);				// Resolution limit
	double			nuavg(0), nustd(0); 		// Values for rescaling
	Bstring			rpsfile;					// Radial power spectrum
	Bstring			fscfile;					// FSC curve
	Bstring			reffile;					// Reference map
    Bstring			maskfile;					// Mask to exclude regions in reciprocal space
    Bstring			realmask;					// Real space mask
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "invert" ) setinvert = 1;
		if ( curropt->tag == "dose" )
			if ( ( weigh_dose = curropt->value.real() ) < 0.001 )
				cerr << "-dose: The dose per frame must be specified!" << endl;
		if ( curropt->tag == "normalize" )
			if ( ( normalize = curropt->value.integer() ) < 1 )
				cerr << "-normalize: A normalization option must be specified!" << endl;
		if ( curropt->tag == "root" ) root = 1;
		if ( curropt->tag == "square" ) square = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.1 )
				cerr << "-resolution: A resolution must be specified!" << endl;
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			else if ( nustd <= 0 )
				cerr << "-rescale: A positive standard deviation must be specified!" << endl;
		}
		if ( curropt->tag == "RPS" )
			rpsfile = curropt->filename();
		if ( curropt->tag == "FSC" )
			fscfile = curropt->filename();
		if ( curropt->tag == "Reference" )
			reffile = curropt->filename();
 		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
 		if ( curropt->tag == "mask" )
			realmask = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
#ifdef HAVE_GCD
	fftwf_init_threads();
	fftwf_plan_with_nthreads(system_processors());
#endif
	if ( verbose )
		cout << "Number of threads:              " << system_processors() << endl;
	
	// Read image file
	Bimage*		p = read_img(argv[optind++], 1, -1);
	if ( p == NULL ) bexit(-1);
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();			// Preserve old type
	
	if ( sam.volume() > 0 ) p->sampling(sam);

	if ( setinvert ) p->invert();
	
	Bplot*		plot = NULL;
	if ( fscfile.length() ) {
		plot = xml_read_fsc(fscfile);
		if ( !plot ) {
			cerr << "Error: The FSC file " << fscfile << " was not read!" << endl;
			bexit(-1);
		}
	} else if ( rpsfile.length() ) {
		plot = read_rps(rpsfile);
		if ( !plot ) {
			cerr << "Error: The RPS file " << rpsfile << " was not read!" << endl;
			bexit(-1);
		}
	}
	
	Bimage*		pref = NULL;
    Bimage*		pmask = NULL;
    Bimage*		prealmask = NULL;
	if ( reffile.length() && ( pref = read_img(reffile, 1, -1) ) ) {
		if ( maskfile.length() )
			pmask = read_img(maskfile, 1, -1);
	}

	if ( realmask.length() ) {
		prealmask = read_img(realmask, 1, -1);
		if ( prealmask ) {
			p->multiply(prealmask);
			delete prealmask;
		}
	}
	
	if ( pref ) {
		p->fspace_weigh(pref, pmask, resolution);
		delete pref;
	} else if ( plot ) {
		if ( fscfile.length() )
			p->fspace_weigh_FSC_curve(plot, resolution);
		else if ( rpsfile.length() )
			p->fspace_weigh_RPS_curve(plot, resolution);
	} else if ( weigh_dose ) {
		p->fspace_weigh_dose(weigh_dose);
	} else if ( normalize > 0 ) {
		p->fspace_normalize_radial(pmask, resolution, (normalize==1));
	} else if ( root ) {
		p->fspace_sqrt_amp();
	} else if ( square ) {
		p->fspace_square_amp();
	} else {
		p->fspace_weigh_C_curve(resolution);
	}

	if ( plot ) delete plot;
	if ( pmask ) delete pmask;

	if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
	
	if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}
	
	delete p;
	
#ifdef HAVE_GCD
	fftwf_cleanup_threads();
#endif

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
