/**
@file	beditimg.cpp
@brief	Segment images into density regions
@author Bernard Heymann
@date	Created: 20000901
@date	Modified: 20190506
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: beditimg [options] [input.img] output.img",
"------------------------------------------------",
"Creates a new image or modifies an input image by inserting or deleting regions.",
" ",
"Actions:",
"-create 20,60,30         Create a new image of this size (pixels/voxels).",
"                         Note: No input image is read and the first file name is taken as the output image.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
"-line 14,3,44,25,86,2    Start and end of line (voxels).",
"-box 0,3,44,50,30,44     Start (lower-left-back voxel) and size of box (voxels).",
"-sphere 22,5,39,120      Center and radius from center of sphere (voxels).",
"-shell 32,45,23,20,120   Center and inner and outer radii of shell (voxels).",
"-cylinder 41,76,9,12,23  Center, radius and height of cylinder (voxels).",
"-gaussian 45,55,12,14.5  Center and sigma of gaussian sphere (voxels).",
"-quadratic 1,2,1,2,1,2,1 Create a quadratic image (7 parameters).",
"-chirp 0.01,45           Create a chirp image: frequency scale and shift (angstrom & degrees).",
"-gradient quad           Correct for the gradient in an image: linear or quadratic.",
"-tetrahedron 9,-3,4,...  Extract a tetrahedron defined by four sets of 3-valued vectors.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-origin 60.4,-32.2,44    Origin of image (voxels).",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-images 12               Number of images (default: 1).",
"-fill 0.02               Fill value: average (default), background, or value.",
"-background -1.2         Fill a newly created image with a background value.",
"-edgewidth 10            Gaussian width of edge of region to delete (voxels, <0 inverts selection).",
"-wrap                    Turn wrapping on (default off).",
" ",
"Parameters for placing an image (use with -place):",
"-location 33,12,56       Location to place an input image.",
"-rotate 0,0.4,0.7,52     Rotation applied before placing an input image.",
" ",
"Input:",
"-place img.pif           An image to place within a larger image.",
" ",
NULL
};

int 	main(int argc, char* argv[])
{
    // Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	int				set_origin(0);
	Vector3<double>	origin;						// Origin in file
	Vector3<long>	nusize;						// Size of new file
	long 			nimg(1);					// Number of images to create
	double			nuavg(0), nustd(0); 		// For rescaling
	int 			fill_type(FILL_AVERAGE);	// Use average
	double			fill(0);					// Default fill is average
	double			background(0);				// Default background is zero
	Vector3<double>	linestart, lineend(-1,-1,-1);	// Start and end of line
	Vector3<double>	boxstart;					// Start of box
	Vector3<long>	boxsize;					// Size of box
	Vector3<double>	spherecenter;				// Sphere center
	double			sphereradius(0); 			// Sphere radius
	Vector3<double>	shellcenter;				// Shell center
	double			shellminrad(0); 			// Shell minimum radius
	double			shellmaxrad(0); 			// Shell maximum radius
	Vector3<double>	cylcenter;					// Cylinder center
	double			cylradius(0);				// Cylinder radius
	double			cylheight(0);				// Cylinder height
	double			sigma(0);					// Gaussian sphere width
	double			quad[7] = {0,0,0,0,0,0,0};	// Quadratic parameters
	double			width(0); 					// Gaussian width minimum = sharp edge
	double			chirp_scale(0);				// Chirp frequency scale
	double			chirp_shift(0);				// Chirp frequency shift (degrees)
	int 			correct_gradient(0);		// Linear gradient correction flag
	int 			tetrahedron(0);				// Flag to extract a tetrahedron
	Vector3<double> tet[4];						// Tetrahedron
	int 			setwrap(0);					// Wrapping flag
	Vector3<double>	loc;						// Location for image to place
	double			angle(0);					// Rotation angle
	Vector3<double>	axis(0.0,0.0,1.0);			// Rotation axis
	Bstring			imgfile;					// Image to place

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		long 			i, j;
		vector<double>	d = curropt->value.split_into_doubles(",");
		if ( curropt->tag == "create" )
			nusize = curropt->size();
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "images" )
			if ( ( nimg = curropt->value.integer() ) < 1 )
				cerr << "-images: The number of images must be specified!" << endl;
		if ( curropt->tag == "rescale" )
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "background" )
			if ( ( background = curropt->value.real() ) == 0 )
				cerr << "-background: A background value must be specified!" << endl;
		if ( curropt->tag == "line" )
			curropt->line(linestart, lineend);
		if ( curropt->tag == "box" )
			if ( curropt->box(boxstart, boxsize) < 6 )
				cerr << "-box: A position and size must be specified!" << endl;
		if ( curropt->tag == "edgewidth" ) {
			width = curropt->value.real();
	        if ( fabs(width) < 0.001 )
				cerr << "-edgewidth: An edge width must be specified!" << endl;
		}
		if ( curropt->tag == "sphere" )
			if ( curropt->values(spherecenter[0], 
					spherecenter[1], spherecenter[2], sphereradius) < 4 )
				cerr << "-sphere: A position and radius must be specified!" << endl;
		if ( curropt->tag == "shell" ) {
			if ( d.size() < 5 )
				cerr << "-shell: A position and two radii must be specified!" << endl;
			else {
				for ( i=0; i<3; i++ ) shellcenter[i] = d[i];
				shellminrad = d[3];
				shellmaxrad = d[4];
			}
		}
		if ( curropt->tag == "cylinder" ) {
			if ( d.size() < 5 )
				cerr << "-cylinder: A position, radius and height must be specified!" << endl;
			else {
				for ( i=0; i<3; i++ ) cylcenter[i] = d[i];
				cylradius = d[3];
				cylheight = d[4];
			}
		}
		if ( curropt->tag == "gaussian" ) {
			if ( d.size() < 4 )
				cerr << "-gaussian: 4 parameters must be specified!" << endl;
			else {
				for ( i=0; i<3; i++ ) spherecenter[i] = d[i];
				sigma = d[3];
			}
		}
		if ( curropt->tag == "quadratic" ) {
			if ( d.size() < 7 )
				cerr << "-quadratic: 7 parameters must be specified!" << endl;
			else
				for ( size_t i=0; i<d.size(); i++ ) quad[i] = d[i];
		}
		if ( curropt->tag == "chirp" ) {
			if ( curropt->values(chirp_scale, chirp_shift) < 1 )
				cerr << "-chirp: A frequency scale must be specified!" << endl;
			else
				chirp_shift *= M_PI/180.0;
		}
		if ( curropt->tag == "gradient" ) {
			if ( curropt->value[0] == 'q' ) correct_gradient = 2;
        	else correct_gradient = 1;
		}
		if ( curropt->tag == "tetrahedron" ) {
			if ( d.size() < 12 )
				cerr << "-tetrahedron: Twelve values must be specified!" << endl;
			else {
				for ( i=0; i<4; i++ )
					for ( j=0; j<3; j++ )
						tet[i][j] = d[i];
				tetrahedron = 1;
			}
		}
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "wrap" )
			setwrap = 1;
 		if ( curropt->tag == "location" )
			if ( curropt->values(loc[0], loc[1], loc[2]) < 3 )
				cerr << "-location: All three location coordinates must be specified!" << endl;
		if ( curropt->tag == "rotate" ) {
			if ( curropt->values(axis[0], axis[1], axis[2], angle) < 4 )
				cerr << "-rotate: All three vector elements and an angle must be specified!" << endl;
			else
				angle *= M_PI/180.0;
		}
		if ( curropt->tag == "place" )
			imgfile = curropt->filename();
	}
	option_kill(option);
        
	double		ti = timer_start();
	
    Bimage* 	p = NULL;
	if ( nusize[0] > 0 && nusize[1] > 0 && nusize[2] > 0 ) {
		if ( nudatatype == Unknown_Type ) nudatatype = Float;
		p = new Bimage(nudatatype, TSimple, nusize, 1);
		if ( background ) p->fill(background);
	} else {
    	p = read_img(argv[optind++], 1, -1);
		if ( !p ) {
			cerr << "Error: No input file read or output size specified!" << endl;
			bexit(-1);
		}
		if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();	// Preserve the old type
	}
	
	if ( fill_type == FILL_AVERAGE ) fill = p->average();
	if ( fill_type == FILL_BACKGROUND ) fill = p->background(long(0));
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}
		
	if ( sam.volume() ) p->sampling(sam);
	
	if ( lineend[2] >= 0 )
		p->line(linestart, lineend, width, fill_type, fill);

	if ( boxsize.volume() > 0 )
		p->shape(0, boxsize, boxstart, width, fill_type, fill);

	if ( sphereradius )
		p->sphere(spherecenter, sphereradius, width, fill_type, fill);

	if ( sigma > 0 )
		p->gaussian_sphere(0, spherecenter, sigma, fill);
	
	if ( shellmaxrad > 0 ) {
		if ( setwrap )
			p->shell_wrap(0, shellcenter, shellminrad, shellmaxrad, width, fill_type, fill);
		else
			p->shell(0, shellcenter, shellminrad, shellmaxrad, width, fill_type, fill);
	}
	
	if ( cylradius && cylheight )
		p->cylinder(cylcenter, cylradius, cylheight, width, fill_type, fill);

	if ( quad[6] ) 
    	p->quadric(quad);

	if ( chirp_scale )
    	p->chirp(chirp_scale, chirp_shift);

	if ( imgfile.length() ) {
		Bimage*		pp = read_img(imgfile, 1, 0);
		if ( angle ) pp->rotate(axis, angle);
		p->place(0, pp, loc, 1e30, 1, 0, 0);
		delete pp;
	}

	if ( tetrahedron ) {
		Bimage*		pnu = p->extract_tetrahedron(tet, fill_type, fill);
		delete p;
		p = pnu;
	}
	
	if ( correct_gradient == 1 ) p->gradient_correction();
	else if ( correct_gradient == 2 ) p->quadric_correct(p->quadric_fit());

	if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
	
	if ( nimg > 1 ) {
		Bimage*		pm = p->copy(nimg);
		for ( long n=1; n<nimg; ++n ) pm->replace(n, p);
		delete p;
		p = pm;
	}

    // Write an output file if a file name is given
    if ( optind < argc ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}
	
	delete p;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

