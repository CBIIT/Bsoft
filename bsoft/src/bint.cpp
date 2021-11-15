/**
@file	bint.cpp
@brief	Interpolation of 2D and 3D images.
@author Bernard Heymann
@date	Created: 19990904
@date	Modified: 20190424
**/

#include "rwimg.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bint [options] input.img output.img",
"------------------------------------------",
"Applies geometric transformations requiring interpolation to 2D and 3D images.",
"The transformation is set up so that only a single interpolation is done",
"to preserve information as far as possible.",
" ",
"Actions:",
"-size 10,50,8            New size (pixels).",
"-translate 0,-10,30      Translation after rotating (default 0,0,0).",
"-rotate -2.5,5.5,9,45    Rotate around a vector by an angle (default 0,0,0,0).",
"-toview 0.3,0.5,0.8,33   Rotate to view: vector {xyz} and angle (default 0,0,1,0).",
"-fromview 0.3,0.5,0.8,33 Rotate from view: vector {xyz} and angle (default 0,0,1,0).",
"-toeuler 16,120,85       Euler angle rotation: Phi, theta, psi (around axes z, y, z).",
"-fromeuler 16,120,85     Euler angle rotation: Psi, theta, phi (around axes z, y, z).",
"-scale 1.5,1.5,1.5       Relative scale (changes the units).",
"-enlarge 2,3,1           Enlarge the image by the given integer scale.",
"-integer 2,32            Special integer interpolation by given kernel size (changes units).",
"                         If old origin falls on a voxel, new origin also falls on a voxel.",
"                         Second value forces odd/even dimensions (default all odd).",
"-bin 3,2,1               Binning by given kernel size - one value = isotropic binning.",
"                         (i.e., 2 = binning with a 2x2 kernel, changes the units).",
"-median 3                Median binning (changes the units).",
"-invert                  Invert density in the image.",
"-slices                  Interpret multiple 2D images as z-slices of a single 3D image",
"-images                  Interpret slices of a single 3D image as 2D images",
"-rescale -0.1,5.2        Rescale output data to average and standard deviation.",
"-project                 Project a 3D image along the z-axis after transformation.",
"-Skew +                  + skews and - removes skewing (default off).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.2,0.7,2.0    Sampling before transformations (angstrom/pixel).",
"-newsampling 1.2,0.7,2.0 Sampling set after transformations (angstrom/pixel).",
"-origin 10,-10,20        Origin for rotation (default 0,0,0).",
"-neworigin 10,-10,20     Origin set after transformations.",
"-fill 127                Fill value for resizing (default background).",
"-unitcell 10,23,77,90,90,90 Unit cell parameters.",
" ",
"Output:",
"-compression 2           Compression type: 5=LZW (TIFF only).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);		// Conversion to new type
	Vector3<double> sampling;						// Units for the three axes (A/pixel)
	Vector3<double> newsampling;					// Units for the three axes (A/pixel)
	int 			setinvert(0);					// Flag to invert density
	int				znswitch(0);					// 0=not, 1=n2z, 2=z2n
	int 			setresize(0);
	int 			settransform(0);
	int 			setproject(0);
	int 			setrescale(0);
	double			nuavg(0), nustd(0);				// Rescaling to average and stdev
	int 			rotate(0); 						// Rotation flag
	int 			backward(0); 					// Backward rotation flag
	int 			set_skew(0);					// Skewing flag
	Vector3<long>	bin = {1,1,1};					// No binning
	int 			median_binning(0);				// No median binning
	int 			integer_factor(0);				// No special integer interpolation
	int				intbin_flag(21);				// Odd dimesnions for integer interpolation
	Vector3<long>	newsize;
	Vector3<double>	origin;							// Origin for rotation
	int				set_origin(0);					// Flag to set origin
	int				set_nuorigin(0);				// Flag to set new origin
	Vector3<double>	nuorigin;						// Origin set after transformations
	int 			fill_type(FILL_BACKGROUND);		// Type of fill value
	double			fill(0);						// Fill value for new areas
	Vector3<long>	enlarge;						// Scale for integer enlargement
	Vector3<double> scale(1.0,1.0,1.0);
	Vector3<long>	shift;
	Vector3<double>	translate;
	double			angle(0);						// Rotation angle
	Vector3<double>	axis(0.0,0.0,1.0);				// Rotation axis
	View			view;							// View to rotate to
	Euler			euler;							// Euler angles
	UnitCell		uc(0,0,0,M_PI_2,M_PI_2,M_PI_2);
	int				compression(0);				// Output compression type
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "slices" )
			znswitch = 1;
		if ( curropt->tag == "images" )
			znswitch = 2;
		if ( curropt->tag == "rescale" ) {
        	if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			else
				setrescale = 1;
			if ( nustd <= 0 ) setrescale = 0;
		}
		if ( curropt->tag == "sampling" )
			sampling = curropt->scale();
		if ( curropt->tag == "newsampling" )
			newsampling = curropt->scale();
		if ( curropt->tag == "size" ) {
			newsize = curropt->size();
			settransform = 1;
		}
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "neworigin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_nuorigin = 2;
			} else {
				nuorigin = curropt->origin();
				set_nuorigin = 1;
			}
		}
		if ( curropt->tag == "translate" ) {
			translate = curropt->vector3();
			if ( translate.length() )
				settransform = 1;
			else
				cerr << "-translate: At least one translation value must be specified!" << endl;
		}
		if ( curropt->tag == "rotate" ) {
			if ( curropt->values(axis[0], axis[1], axis[2], angle) < 4 )
				cerr << "-rotate: All three vector elements and an angle must be specified!" << endl;
			else {
				settransform = 1;
				rotate = 1;
				angle *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "toview" ) {
			view = curropt->view();
			settransform = 1;
			rotate = 1;
			backward = 1;
		}
		if ( curropt->tag == "fromview" ) {
			view = curropt->view();
			settransform = 1;
			rotate = 1;
		}
		if ( curropt->tag == "toeuler" ) {
			euler = curropt->euler();
			view = euler.view();
			settransform = 1;
			rotate = 1;
			backward = 1;
		}
		if ( curropt->tag == "fromeuler" ) {
			euler = curropt->euler();
			view = euler.view();
			settransform = 1;
			rotate = 1;
		}
		if ( curropt->tag == "scale" ) {
			scale = curropt->scale();
			if ( scale.volume() )
				settransform = 1;
			else
				cerr << "-scale: At least one value must be specified!" << endl;
		}
		if ( curropt->tag == "enlarge" ) {
			enlarge = curropt->scale();
			if ( enlarge.volume() < 1 )
				cerr << "-enlarge: All three values must be specified!" << endl;
		}
		if ( curropt->tag == "integer" )
			if ( curropt->values(integer_factor, intbin_flag) < 1 )
				cerr << "-integer: An integer must be specified!" << endl;
		if ( curropt->tag == "bin" ) {
			if ( ( i = curropt->values(bin[0], bin[1], bin[2]) ) < 1 )
				cerr << "-bin: At least one bin size must be specified!" << endl;
			else
				if ( i < 2 ) bin[2] = bin[1] = bin[0];
		}
		if ( curropt->tag == "median" )
			if ( ( median_binning = curropt->value.integer() ) < 1 )
				cerr << "-median: A binning size must be specified!" << endl;
		if ( curropt->tag == "project" )
			setproject = 1;
		if ( curropt->tag == "Skew" ) {
			settransform = 1;
			if ( curropt->value.contains("+" ) )
				set_skew = 1;
			else if ( curropt->value.contains("-" ) )
				set_skew = -1;
			else {
				cerr << "-Skew: A skewing direction (+/-) must be specified!" << endl;
				settransform = 0;
			}
		}
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "unitcell" ) {
        	uc = curropt->unit_cell();
			settransform = 1;
        }
		if ( curropt->tag == "compression" )
			if ( ( compression = curropt->integer() ) < 1 )
				cerr << "-compression: A number must be specified!" << endl;
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	// Read image file
	int 		dataflag(0);
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	if ( p == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	if ( znswitch == 1 )
		if ( p->images_to_slices() < 0 ) bexit(-1);
	
	if ( znswitch == 2 ) 
		if ( p->slices_to_images() < 0 ) bexit(-1);

	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();		// Preserve the old type
	else if ( nudatatype > p->data_type() )
		p->change_type(nudatatype);
	
	if ( sampling.volume() > 0 ) p->sampling(sampling);
	
	// Scaling is always positive
	if ( fabs(scale.volume()) < 0.001 ) {
		cerr << "Error: Scales too small!" << endl;
		bexit(-1);
	}
	
	if ( setinvert ) p->invert();

	if ( p->background(long(0)) < 1e-37 ) p->calculate_background();
	if ( fill_type == FILL_AVERAGE ) fill = p->average();
	if ( fill_type == FILL_BACKGROUND ) fill = p->background(long(0));
	
	if ( uc.check() ) p->unit_cell(uc);
	else uc = p->unit_cell();
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}

	if ( bin[0] )
		p->bin(bin);
	else if ( median_binning )
		p->median_bin(median_binning);
	else if ( integer_factor )
		p->integer_interpolation(integer_factor, intbin_flag);
	
	if ( enlarge.volume() )
		p->enlarge(enlarge);
	
	Matrix3		mat(1);
	Bimage*		pnu = NULL;
	
	if ( settransform ) {
	
		if ( set_skew ) {
			mat = uc.skew_rotation(set_skew==-1);
		} else if ( rotate ) {
			if ( view[2] < 1 || view.angle() != 0 ) {
				if ( p->sizeZ() < 2 ) view = View(0, 0, 1, view.angle());
				mat = view.matrix();
				if ( backward )
					mat = mat.transpose();
			} else if ( angle != 0 ) {
				if ( p->sizeZ() < 2 ) {
					axis[0] = axis[1] = 0;
					axis[2] = 1;
				}
				mat = Matrix3(axis, angle);
//				cout << mat << endl;
			}
			view = View(mat.transpose());
			if ( verbose & VERB_FULL )
				cout << "View: " << view << endl;
		}

		if ( fabs(scale[0] - 1) < 1e-10 ) setresize += 1;
		if ( fabs(scale[1] - 1) < 1e-10 ) setresize += 1;
		if ( fabs(scale[2] - 1) < 1e-10 ) setresize += 1;
		shift[0] = (long)floor(translate[0]+0.5);
		shift[1] = (long)floor(translate[1]+0.5);
		shift[2] = (long)floor(translate[2]+0.5);
		if ( fabs(translate[0] - shift[0]) < 1e-10 ) setresize += 1;
		if ( fabs(translate[1] - shift[1]) < 1e-10 ) setresize += 1;
		if ( fabs(translate[2] - shift[2]) < 1e-10 ) setresize += 1;
		
		if ( set_skew || rotate || setresize < 6 ) {
			if ( newsize.volume() < 1 ) newsize = p->size();
			pnu = p->transform(newsize, scale, p->image->origin(), translate, mat, fill_type, fill);
			delete p;
			p = pnu;
		} else {		 //cerr << "minimal!" << endl;
			p->resize(newsize, shift, fill_type, fill);
		}
		
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG bint: Transformations done" << endl;
	}
	
	if ( setrescale ) p->rescale_to_avg_std(nuavg, nustd);

	if ( setproject ) p->project('z', 1);

	if ( newsampling.volume() > 0 ) p->sampling(newsampling);

	if ( set_nuorigin ) {
		if ( set_nuorigin == 2 ) p->origin(p->size()/2);
		else p->origin(nuorigin);
	}
			
	if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, compression);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

