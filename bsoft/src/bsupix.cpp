/**
@file	bsupix.cpp
@brief	Segment images
@author Bernard Heymann
@date	Created: 20160320
@date	Modified: 20210303
**/

#include "rwimg.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

vector<Bsuperpixel>	superpixel_from_json(JSvalue& root);
JSvalue		json_from_superpixel(vector<Bsuperpixel>& seg);

/* Usage assistence */
const char* use[] = {
" ",
"Usage: bsupix [options] input.img output.img",
"--------------------------------------------",
"Segments an image into superpixels.",
" ",
"Actions:",
"-average 7,5,3           Averaging/smoothing filter: kernel size.",
"-variance 11             Calculate a local variance image using the given size kernel.",
"-gaussian 11,2.6         Gaussian smoothing filter: kernel size and sigma.",
"-superpixels 123         Segment into superpixels with the given step size.",
"-impose 3                Assign segment feature:",
"                         1: segment average",
"                         2: lowest neighboring segment average",
"                         3: difference from lowest neighboring segment average",
"                         4: highest neighboring segment average",
"                         5: difference from highest neighboring segment average",
"-colorize                Generate random segment colors.",
"-pattern neigbor         Calculate local or neighbor binary pattern.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-datatype u              Force writing of a new data type.",
" ",
"Parameters for superpixel segmentation:",
"-iterations 50           Number of segmentation iterations (default 10).",
"-bin 3                   Number of binning levels (default 1).",
"-weight 0.5              Color weight (default 0.2).",
"-stop 1.5                Stopping condition (default 1%).",
" ",
//"Parameters for segmentation analysis:",
//" ",
"Input:",
"-maskin mask.mrc         Multi-level input mask.",
"-input segfile.json      Segmentation input file.",
" ",
"Output:",
"-maskout mask.mrc        Multi-level output mask.",
"-output segfile.json     Segmentation output file.",
" ",
NULL
};

int 	main(int argc, char* argv[])
{
    // Initialize variables
	Vector3<long>	average_kernel;			// Average filter kernel size
	long 			var_kernel(0);			// Variance filter kernel size
	long			gauss_kernel(0), sigma(0);	// Gaussian kernel size and sigma
	long			step(0);				// Superpixel step size
	double			colorweight(0.2);		// Color weighting
	long			iterations(10);			// Segmentation iterations
	long			bin(1);					// Number of binning levels
	double			stop(1);				// Stopping condition
	int				impose(0);				// Flag to impose segment colors
	int				colorize(0);			// Generate color mask
	DataType		nudatatype(Unknown_Type);	// Conversion to new type
	Bstring			maskin;					// Multi-level input mask file name
	Bstring			maskout;				// Multi-level output mask file name
	Bstring			segin;					// Segmentation input file name
	Bstring			segout;					// Segmentation output file name

	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "average" ) {
			if ( ( i = curropt->values(average_kernel[0],
					average_kernel[1], average_kernel[2]) ) < 1 )
				cerr << "-average: The kernel edge size must be specified!" << endl;
			else
				if ( i < 2 ) average_kernel[2] = average_kernel[1] = average_kernel[0];
		}
		if ( curropt->tag == "variance" )
			if ( ( var_kernel = curropt->integer() ) < 2 )
				cerr << "-variance: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "gaussian" )
			if ( curropt->values(gauss_kernel, sigma) < 2 )
				cerr << "-gaussian: The kernel edge size and sigma must be specified!" << endl;
		if ( curropt->tag == "superpixels" )
			if ( ( step = curropt->integer() ) < 2 )
				cerr << "-superpixels: A step size greater than 1 must be specified!" << endl;
		if ( curropt->tag == "impose" )
			if ( ( impose = curropt->integer() ) < 1 )
				cerr << "-impose: One option must be specified!" << endl;
		if ( curropt->tag == "colorize" )
			colorize = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "iterations" )
			if ( ( iterations = curropt->integer() ) < 1 )
				cerr << "-iterations: At least one iteration must be specified!" << endl;
		if ( curropt->tag == "bin" )
			if ( ( bin = curropt->integer() ) < 1 )
				cerr << "-bin: At number of binning levels must be specified!" << endl;
		if ( curropt->tag == "weight" )
			if ( ( colorweight = curropt->real() ) < 0.001 )
				cerr << "-weight: A color weight must be specified!" << endl;
		if ( curropt->tag == "stop" )
			if ( ( stop = curropt->real() ) < 0.001 )
				cerr << "-stop: A stopping percentage must be specified!" << endl;
 		if ( curropt->tag == "maskin" )
			maskin = curropt->filename();
 		if ( curropt->tag == "maskout" )
			maskout = curropt->filename();
 		if ( curropt->tag == "input" )
			segin = curropt->filename();
 		if ( curropt->tag == "output" )
			segout = curropt->filename();
   }
	option_kill(option);
        
	double		ti = timer_start();
	
    // Read the input files
	int 		dataflag = 0;
	if ( optind < argc - 1 || maskin.length() ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	
	if ( !p )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	if ( !dataflag ) bexit(0);
	
	if ( p->standard_deviation() <= 0 ) p->statistics();
	
	if ( average_kernel[0] > 0 ) p->filter_average(average_kernel);

	if ( var_kernel > 1 ) p->variance(var_kernel);

	if ( sigma > 0 ) p->filter_gaussian(gauss_kernel, sigma);

	Bimage*				pmask = NULL;
	
	if ( maskin.length() ) {
		pmask = read_img(maskin, 1, -1);
		if ( !p->check_if_same_size(pmask) ) {
			cerr << "Error: The image and mask must be the same size!" << endl;
			bexit(-1);
		}
	}

	vector<Bsuperpixel>	seg;

	if ( segin.length() ) {
		if ( verbose )
			cout << "Reading " << segin << endl;
		JSvalue		sp = JSparser().parse(segin.c_str());
		seg = superpixel_from_json(sp);
/*		if ( pmask )
			vector<long>		vstep = {step,step,step};
			vstep *= 2;
			if ( vstep[2] > z ) vstep[2] = z;
			p->superpixels_update(pmask, vstep, colorweight, seg);
		}*/
	} else if ( step ) {
		if ( bin == 1 )
			seg = p->superpixels(step, colorweight, iterations, stop);
		else
			seg = p->superpixels(step, colorweight, iterations, bin, stop);
		pmask = p->next;
		p->next = NULL;
	}

	if ( segout.length() ) {
		if ( verbose )
			cout << "Writing " << segout << endl << endl;
		JSvalue		sp = json_from_superpixel(seg);
		string		fn(segout.c_str());
		sp.write(fn);
	}

	if ( impose ) {
		if ( seg.size() )
			p->impose_superpixels(pmask, seg, impose);
		else if ( impose == 1 )
			p->level_masked_stats(pmask);
	}
	
	if ( colorize ) pmask->levelmask_colorize();
	
	if ( maskout.length() )
		write_img(maskout, pmask, 0);
	
    // Write an output multilevel mask file if a file name is given
    Bstring		filename(argv[optind]);
    if ( optind < argc ) {
		p->change_type(nudatatype);
    	write_img(filename, p, 0);
    	if ( impose == 1 && p->next ) {
			filename = filename.pre_rev('.') + "_var." + filename.post_rev('.');
			p->next->change_type(nudatatype);
			write_img(filename, p->next, 0);
    	}
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	return 0;
}

vector<Bsuperpixel>	superpixel_from_json(JSvalue& root)
{
	vector<Bsuperpixel>	seg;

	vector<JSvalue>		splist = root.array();
	
	for ( auto it = splist.begin(); it != splist.end(); ++it ) {
		JSvalue&		jv = *it;
		Bsuperpixel		sp;
		vector<double>	coor = jv["coordinates"].array_real();
		vector<double>	avg = jv["channel_averages"].array_real();
		vector<double>	var = jv["channel_variances"].array_real();
		vector<long>	nb = jv["neighbors"].array_integer();
		sp.channels(avg.size());
		sp.index(jv["id"].integer());
		sp.count(jv["count"].integer());
		sp.coordinates(coor);
		sp.channels(avg);
		sp.variances(var);
		sp.weight(jv["weight"].real());
		sp.neighbors(nb);
		seg.push_back(sp);
	}
	
	if ( verbose )
		cout << "Superpixels:                   " << seg.size() << endl << endl;
		
	return seg;
}

JSvalue		jsvalue_from_neighbors(Bsuperpixel& sp)
{
	vector<JSvalue>		vec;

	for ( long j=0; j<NNEIGHBOR && sp.neighbor(j) >= 0; ++j )
		vec.push_back(sp.neighbor(j));
	
	return JSvalue(vec);
}


JSvalue		json_from_superpixel(vector<Bsuperpixel>& seg)
{
	vector<JSvalue>			splist;
	
	for ( auto it = seg.begin(); it != seg.end(); ++it ) {
		map<string, JSvalue>	sp;
		sp["id"] = it->index();
		sp["count"] = it->count();
		sp["coordinates"] = JSvalue(it->coordinates());
		sp["channel_averages"] = JSvalue(it->channels());
		sp["channel_variances"] = JSvalue(it->variances());
		sp["weight"] = it->weight();
		sp["neighbors"] = jsvalue_from_neighbors(*it);
		splist.push_back(sp);
	}
	
	return splist;
}

