/**
@file	balign.cpp
@brief	Aligns two images.
@author Bernard Heymann
@date	Created: 20021111
@date	Modified: 20160603
**/

#include "rwimg.h"
#include "ps_marker.h"
#include "linked_list.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Transform	img_align(Bimage* p1, Bimage* p2, Vector3<long> tile_size, 
				double res_lo, double res_hi, double max_shift, int filter_flag, int refine_flag);

// Usage assistance
const char* use[] = {
" ",
"Usage: balign [options] input1.img input2.img output.img",
"--------------------------------------------------------",
"The second input image is aligned to the first, transformed and output.",
"(Note: This program is intended to align very large micrograph images).",
" ",
"Actions:",
"-add                     Add input images after transforming the second.",
"-filterextremes          Filter micrograph extremes before aligning.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-tile 1024,1024,5        Size of tiles for correlation (default 512,512,1).",
"-resolution 30,500       High and low resolution limits for cross-correlation (default 0.1,1000).",
"-maxshift 55             Maximum allowed shift (default 1/4 of tile).",
" ",
NULL
};

int			main(int argc, char* argv[])
{
	// Initialize all settings
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;				// Units for the three axes (A/pixel)
	int 			add(0);
	Vector3<long> 	tile_size(512,512,1);
	double			res_hi(0.1);
	double			res_lo(1000);
	double			max_shift(0);
	int				filter_flag(0);
	int 			refine_flag(1);
	    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "add" )
			add = 1;
		if ( curropt->tag == "filterextremes" )
			filter_flag = 1;
		if ( curropt->tag == "tile" ) tile_size = curropt->size();
		if ( curropt->tag == "resolution" ) {
    	    if ( curropt->values(res_hi, res_lo) < 1 )
				cerr << "-resolution: A high resolution limit must be specified." << endl;
			else if ( res_hi > res_lo )
				swap(res_hi, res_lo);
        }
		if ( curropt->tag == "maxshift" )
    	    if ( ( max_shift = curropt->value.real() ) < 1 )
				cerr << "-maxshift: A maximum shift in pixels must be specified." << endl;
    }
	option_kill(option);

	double		ti = timer_start();
	
    // Read the input files
	int 		dataflag(0);
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p1 = read_img(argv[optind++], dataflag, -1);
	if ( p1 == NULL ) {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
    Bimage* 	p2 = read_img(argv[optind++], dataflag, -1);
	if ( p2 == NULL ) {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	if ( sam.volume() ) {
		p1->sampling(sam);
		p2->sampling(sam);
	}
	
	Transform		t = img_align(p1, p2, tile_size, res_lo, res_hi, max_shift, filter_flag, refine_flag);

	p2->transform(t.scale, t.origin, t.trans, Matrix3(t.axis, t.angle), FILL_USER, p2->average());

	Bimage*			p = NULL;
	if ( add ) p = *p1 + *p2;
	if ( !p ) p = p2;

    // Write an output file if a file name is given
    if ( argc > optind && strspn(argv[optind],"-") != 1 ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}
	
	if ( p != p2 ) delete p;
	delete p1;
	delete p2;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/**
@brief 	Aligns the second image to the first.
@param 	*p1				first image.
@param 	*p2				second image (transformed)
@param 	tile_size		3-valued vector for the size of sub-images.
@param 	res_lo			low resolution limit for cross-correlation.
@param 	res_hi			high resolution limit for cross-correlation.
@param 	max_shift		maximum shift allowed (default 1/4 of tile).
@param 	filter_flag		flag to filter micrograph extremes.
@param 	refine_flag 	flag to turn on refinement of shift.
@return Transform		structure with shift, scale, rotation angle, and R factor.

	The second image is transformed to fit on the first image.

**/
Transform	img_align(Bimage* p1, Bimage* p2, Vector3<long> tile_size, 
				double res_lo, double res_hi, double max_shift, int filter_flag, int refine_flag)
{
	long				i, iter, maxiter(5);
	Vector3<double>		shift;
	Transform			t, best_t;
	best_t.fom = 1e37;
	
	Vector3<long> 		start1, common_size, step_size;
	
	// Find the largest area covering both images
	common_size = p1->size().min(p2->size());
	common_size = common_size.max(1);
	
	// Ensure the tile fits into the images
	tile_size = tile_size.max(1);
	tile_size = tile_size.min(common_size);
	if ( tile_size[0] >= common_size[0] && common_size[0] > 1 ) {
		tile_size[0] = (int) (0.4*common_size[0]);
//		start1[0] = (common_size[0] - tile_size[0])/2;
	}
	if ( tile_size[1] >= common_size[1] && common_size[1] > 1 ) {
		tile_size[1] = (int) (0.4*common_size[1]);
//		start1[1] = (common_size[1] - tile_size[1])/2;
	}
	if ( tile_size[2] >= common_size[2] && common_size[2] > 1 ) {
		tile_size[2] = (int) (0.4*common_size[2]);
//		start1[2] = (common_size[2] - tile_size[2])/2;
	}

	if ( max_shift < 1 ) max_shift = tile_size[0]/4;
		
	Vector3<double>		half_tile;
	half_tile[0] = (tile_size[0] - 1)/2;
	half_tile[1] = (tile_size[1] - 1)/2;
	half_tile[2] = (tile_size[2] - 1)/2;
	
	Bimage*				pex1 = NULL;
	Bimage*				pex2 = NULL;
	Bmarker*			set1 = NULL;
	Bmarker*			set2 = NULL;
	Bmarker*			m1 = NULL;
	Bmarker*			m2 = NULL;
	Bstring				fn1(p1->file_name()), fn2(p2->file_name());
	Bstring				filename = fn1.pre_rev('.') + "_" + fn2.pre_rev('.') + ".ps";
	
	if ( verbose & VERB_LABEL ) {
		cout << "Aligning images ";
		if ( refine_flag ) cout << "with ";
		else cout << "without ";
		cout << "shift refinement." << endl << endl;
	}
	
	if ( common_size[0] > 2*tile_size[0] )
		start1[0] = (common_size[0]
				- (int) (tile_size[0]*(common_size[0]/tile_size[0])))/2;
	if ( common_size[1] > 2*tile_size[1] )
		start1[1] = (common_size[1]
				- (int) (tile_size[1]*(common_size[1]/tile_size[1])))/2;
	if ( common_size[2] > 2*tile_size[2] )
		start1[2] = (common_size[2]
				- (int) (tile_size[2]*(common_size[2]/tile_size[2])))/2;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Common size:                    " << common_size << endl;
		cout << "Tile size:                      " << tile_size << endl;
		cout << "Maximum shift:                  " << max_shift << endl;
		cout << "Reference start:                " << start1 << endl << endl;
	}
	
	if ( res_hi < 0.1 || res_hi > tile_size[0]*p1->sampling(0)[0] ) res_hi = 4.0*p1->sampling(0)[0];
	
	if ( filter_flag ) {
		p1->filter_extremes();
		p2->filter_extremes();
	}
	
	pex1 = p1->extract_tiles(0, start1, common_size, tile_size, step_size, 0);
	
	vector<Vector3<long>>	tile_start(pex1->images());
	
	for ( i=0; i<pex1->images(); i++ ) {
		m1 = (Bmarker *) add_item((char **) &m1, sizeof(Bmarker));
		m2 = (Bmarker *) add_item((char **) &m2, sizeof(Bmarker));
		if ( !set1 ) set1 = m1;
		if ( !set2 ) set2 = m2;
		m1->sel = m2->sel = 1;
		m1->fom = m2->fom = 1;
		m1->loc = pex1->image[i].origin() + half_tile;
		tile_start[i] = pex1->image[i].origin();
	}
	
	for ( iter=0; iter<maxiter && best_t.fom>0.1; iter++ ) {
		pex2 = p2->extract_tiles(0, tile_start, tile_size);
		// Determines the shift of the second image onto the first
		pex1->find_shift(pex2, NULL, res_hi, res_lo, max_shift, 0, refine_flag);
		for ( i=0, m2=set2; i<pex1->images(); i++, m2=m2->next ) {
			shift = pex1->image[i].origin() - pex2->image[i].origin();
			if ( shift[0] > half_tile[0] ) shift[0] -= tile_size[0];
			if ( shift[1] > half_tile[1] ) shift[1] -= tile_size[1];
			if ( shift[2] > half_tile[2] ) shift[2] -= tile_size[2];
			// To get the location in the second image, the shift must be subtracted
			m2->loc = tile_start[i] + half_tile - shift;
		}
		t = markers_find_rottrans(set1, set2, 0.01);
//		if ( verbose & VERB_PROCESS )
			cout << "Fit " << iter+1 << ":\t" << p1->file_name() << tab << p2->file_name() << tab 
				<< t.trans << tab
				<< t.scale << tab 
				<< t.angle*180.0/M_PI << tab << t.fom << endl;
		if ( best_t.fom > t.fom ) best_t = t;
		else iter = maxiter;
		for ( i=0, m1=set1, m2=set2; i<pex1->images(); i++, m1=m1->next, m2=m2->next ) {
			m2->loc[0] = best_t.scale[0]*(m1->loc[0]*cos(best_t.angle) -
					m1->loc[1]*sin(best_t.angle)) + t.trans[0] + 0.5;
			m2->loc[1] = best_t.scale[1]*(m1->loc[0]*sin(best_t.angle) +
					m1->loc[1]*cos(best_t.angle)) + t.trans[1] + 0.5;
			m2->loc[2] = best_t.scale[2]*m1->loc[2] + t.trans[2];
			tile_start[i] = m2->loc - half_tile;
		}
		delete pex2;
	}
	
	if ( verbose & VERB_RESULT ) 
		cout << "Best fit:\t" <<
				p1->file_name() << tab << p2->file_name() << tab <<
				best_t.trans << tab <<
				best_t.scale << tab <<
				best_t.angle*180.0/M_PI << tab << best_t.fom << endl;
	
	delete pex1;

	ps_marker_errors(filename, set1, set2, t, p1->size(), 20);

	kill_list((char *) set1, sizeof(Bmarker));
	kill_list((char *) set2, sizeof(Bmarker));
	
	return best_t;
}

