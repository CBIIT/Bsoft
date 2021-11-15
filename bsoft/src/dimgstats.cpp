/**
@file	dimgstats.cpp
@brief	Calculates statistical measures within a set of images or between sets of images
@author David Belnap and Bernard Heymann
@date	Created: 20051213
@date	Modified: 20190208 (BH)
**/

#include "Bstring.h"
#include "rwimg.h"
#include "img_combine.h"
#include "math_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"
#include "Vector3.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
int			img_fom_sqrt(Bimage* p);
Bimage* 	img_students_t_test(int n, Bimage* p, vector<double>& weight);
Bimage* 	img_students_t_test_equal(int set1, int set2, Bimage* p1, Bimage* p2, vector<double>& weight);
Bimage* 	img_students_t_test_unequal(int set1, int set2, Bimage* p1, Bimage* p2, vector<double>& weight);
Bimage* 	img_f_test(int set1, int set2, Bimage* p1, Bimage* p2, vector<double>& weight);

// Usage assistence
const char* use[] = {
" ",
"Usage: dimgstats [options] input.images",
"---------------------------------------",
"Calculates statistical measures within a set of images or between sets of images.",
"Corresponding pixels are compared--i.e. pixel m in image 1 is compared to pixel m",
"in image 2.",
" ",
"Definitions: ",
"   Variance              var = Sum[ (rho - avg_rho)(rho - avg_rho) ] / (n - 1)",
"     numerator portion    Vn = Sum[ (rho - avg_rho)(rho - avg_rho) ]",
"   Standard deviation    std = sqrt( var )",
"   Student's t-test (t)  SD = sqrt[(Vn1 + Vn2)(1/m + 1/n) / (m + n - 2) ]",
"     approx. equal var     t = ( avg1 - avg2 ) / SD",
"     unequal var           t = ( avg1 - avg2 ) / sqrt[ var1/m + var2/n ]",
"   F-test                  F = var1 / var2",
"avg = average (mean), m = number of data points in set 1, n = number in set 2,",
"rho = value of a given pixel at a given (x,y,z) coordinate",
"SD = standard error of the difference of the means, sqrt = square root",
"Sum = summation of.",
" ",
"See Milligan and Flicker (1987) J. Cell Biol. 105:29-39 for a reference",
"to the use of the Student's t-test in 3D electron microscopy.",
" ",
"Any number of input images may be given and may include the wild card \"*\".",
"(However, VMS does not support Unix style usage of wild cards)",
"All images are converted to floating point.",
" ",
"Actions:",
"-Ttest equal             Do a Student's t-test to compare two sets (equal/unequal).",
"-Ftest                   Do a F-test to compare two sets.",
"-average                 Compute an average image for each set.",
"-stdev                   Compute a standard deviation image for each set.",
"-variance                Compute a variance image for each set.",
"-significance            Compute significance-level (probability) image.",
"-rescale -0.1,5.2        Rescale input images to average, standard deviation.",
" ",
"Parameters:",
"-datatype u              Force writing of a new data type",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-sets 10,12              Number of input files in each set (default = single set).",
"-verbose 7               Verbosity of output",
"-weights 0.1,0.3,0.9     Weights for input files (default 1 for every file)",
" ",
"Output:",
"-output name.img         base name for output image(s), suffix added based on action",
" ",
NULL
};



int 		main(int argc, char **argv)
{
	// Initialize variables
	Bstring   		filebase;					// Output file base name
	Bstring*  		file_list  = NULL;			// List of file names for first and second sets
	Bstring*  		file_list2 = NULL;			// List of file names for second set
	Bstring   		outfile;					// Name of output file
	Bstring   		weight_string;				// input weights

	int       		avgflag(0);					// Flag to output average image
	int       		avgvar(0);					// Flag indicating that average and variance need to be calculated
	int       		fflag(0);					// Flag for F-test
	int       		nfiles;						// Number of files entered
	int       		set1(0), set2(0);			// Number of images in sets 1 and 2, for Student's t-test image
	int       		setweight(0);				// Weight flag, default weights are 1
	int       		signif(0);					// Flag to calculate significance-level image
	int       		stdev(0);					// Flag to calculate and output standard deviation
	int       		sttflag[2] = {0,0};			// Flags for Student's t-test equal or unequal variance
	int       		variance(0);				// Flag to output variance

	double			nuavg(0), nustd(0);			// Average and standard deviation for rescaling
	vector<double>	weight, weight2;			// Weights for each image file (default = 1)

	DataType		nudatatype(Unknown_Type);	// Conversion to new data type
	Bimage			*p1 = NULL, *p2 = NULL;		// Average images for sets 1 and 2
	long			i;							// Index
	Vector3<double>	sam;						// Units for the three axes (A/pixel)


	// Command-line options
	int			optind;
	Boption*	option = get_option_list(use, argc, argv, optind);
	Boption*	curropt;
	for ( curropt = option; curropt; curropt = curropt->next )  {
		if ( curropt->tag == "Ttest" )  {
			avgvar  = 1;
			if ( curropt->value[0] == 'e' ) sttflag[0] = 1;   
			if ( curropt->value[0] == 'u' ) sttflag[1] = 1;   
		}
		if ( curropt->tag == "Ftest" )  {
			fflag = 1;   avgvar  = 1;
		}
		if ( curropt->tag == "average" )  {
        	avgflag = 1;   avgvar  = 1;
		}
		if ( curropt->tag == "stdev" )  {
			stdev = 1;   avgvar  = 1;
		}
		if ( curropt->tag == "variance" )  {
			variance = 1;   avgvar  = 1;
		}
		if ( curropt->tag == "significance" )
			signif = 1;
		if ( curropt->tag == "rescale" )
        	if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified" << endl;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "sets" ) {
        	if ( curropt->values(set1, set2) < 2 )
				cerr << "-sets: Two integers must be specified" << endl;
			else if ( set1<=0 || set2<=0 )
				cerr << "-sets: Each set must have at least one member." << endl;
		}
		if ( curropt->tag == "weights" ) {
			weight_string = curropt->value;
			if ( weight_string.length() < 1 )
				cerr << "-weights: Weights for all files must be specified" << endl;
			else
				setweight = 1;
		}
		if ( curropt->tag == "output" )
        	filebase = curropt->filename();
	}
	option_kill(option);

	double		ti = timer_start();

	if ( optind >= argc )  {
		cerr << "No input files!" << endl;
		bexit(-1);
	}
	if ( (sttflag[0] || sttflag[1] || fflag) && (set1 == 0)  )  {
		cerr << "You must designate two image sets (option \"-sets\") for the t-test or F-test!" << endl;
		bexit(-1);
	}


	// Set up image file names
	nfiles = 0;
	while ( optind < argc ) {
		if ( !set1 || ( set1 && nfiles < set1 ) )
			string_add(&file_list, argv[optind++]);
		else
			string_add(&file_list2, argv[optind++]);
		nfiles++;
	}

	if ( set1 < 1 && set2 < 1 )  // If set option unused, assume only one set
		 set1 = nfiles;
	else                         // If two sets, check that number of files and numbers for sets are in agreement
		if ( (set1+set2) != nfiles )  {
			cerr << "The number of files (" << nfiles << ") does not agree with the numbers designated for sets 1 ("
				<< set1 << ") and 2 (" << set2 << ")." << endl;
			bexit(-1);
		}

	// Set weights for each input file, default is 1.0
	i = 0;
	if ( setweight ) {
		weight = weight_string.split_into_doubles(",");
		if ( weight.size() < nfiles ) {
			cerr << "Error: The number of weights must be equal to the number of input files" << endl;
			bexit(-1);
		}
	} else {
		for ( i=0; i<nfiles; i++ ) weight.push_back(1.0);
	}

	if ( set2 ) 
		for ( i=set1; i<nfiles; i++ ) weight2.push_back(weight[i]);

	// Compute average and variance images
	if ( avgvar )  {
		p1 = img_add_weighed(file_list, weight, nuavg, nustd, 3);
		if ( set2 )
			p2 = img_add_weighed(file_list2, weight2, nuavg, nustd, 3);

		if ( filebase.length() )  {   // write image(s) to file
			if ( sam.volume() > 0 )  {
				p1->sampling(sam);
				p2->sampling(sam);
			}
			p1->change_type(nudatatype);
			if ( p2 ) p2->change_type(nudatatype);
			if ( avgflag )  {
				if ( set2 )
					outfile = filebase.pre_rev('.') + "_avg1." + filebase.post_rev('.');
				else
					outfile = filebase.pre_rev('.') + "_avg." + filebase.post_rev('.');
				write_img(outfile, p1, 0);
				if ( p2 )  {
					outfile = filebase.pre_rev('.') + "_avg2." + filebase.post_rev('.');
					write_img(outfile, p2, 0);
				}
			}
			if ( variance )  {
				if ( set2 )
					outfile = filebase.pre_rev('.') + "_var1." + filebase.post_rev('.');
				else
					outfile = filebase.pre_rev('.') + "_var." + filebase.post_rev('.');
				write_img(outfile, p1->next, 0);
				if ( p2 )  {
					outfile = filebase.pre_rev('.') + "_var2." + filebase.post_rev('.');
					write_img(outfile, p2->next, 0);
				}
			}
		}
	}

	string_kill(file_list);

	// F-test image
	if ( fflag )  {

		Bimage*   pf = img_f_test(set1, set2, p1, p2, weight);

		if ( filebase.length() )  {
			outfile = filebase.pre_rev('.') + "_ftest." + filebase.post_rev('.');
			if ( sam.volume() > 0 )  pf->sampling(sam);
			pf->change_type(nudatatype);
			write_img(outfile, pf, 0);
			if ( signif )  {
				outfile = filebase.pre_rev('.') + "_ftest_sig." + filebase.post_rev('.');
				write_img(outfile, pf->next, 0);
			}
		}

		delete pf;
	}

	// Student's t-test image, equal or unequal variance
	Bimage*  pt = NULL;
	for ( i=0; i<2; i++ ) if ( sttflag[i] ) {
		if ( i == 0 )        // equal variance
			pt = img_students_t_test_equal(set1, set2, p1, p2, weight);
		else if ( i == 1 )   // unequal variance
			pt = img_students_t_test_unequal(set1, set2, p1, p2, weight);
		if ( filebase.length() )  {
			if ( i == 0 )  outfile = filebase.pre_rev('.') + "_stte." + filebase.post_rev('.');
			if ( i == 1 )  outfile = filebase.pre_rev('.') + "_sttu." + filebase.post_rev('.');
			if ( sam.volume() > 0 )  pt->sampling(sam);
			pt->change_type(nudatatype);
			write_img(outfile, pt, 0);
			if ( signif )  {
				if ( i == 0 )  outfile = filebase.pre_rev('.') + "_stte_sig." + filebase.post_rev('.');
				if ( i == 1 )  outfile = filebase.pre_rev('.') + "_sttu_sig." + filebase.post_rev('.');
				write_img(outfile, pt->next, 0);
			}
		}
		delete pt;
	}

	// Compute standard deviation image from variance
	if ( stdev )  {
		img_fom_sqrt(p1);
		if ( p2 )  img_fom_sqrt(p2);
		if ( filebase.length() )  {
			if ( set2 )
				outfile = filebase.pre_rev('.') + "_std1." + filebase.post_rev('.');
			else
				outfile = filebase.pre_rev('.') + "_std." + filebase.post_rev('.');
			write_img(outfile, p1->next, 0);
			if ( p2 )  {
				outfile = filebase.pre_rev('.') + "_std2." + filebase.post_rev('.');
				write_img(outfile, p2->next, 0);
			}
		}
	}

	delete p1;
	delete p2;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/************************************************************************
@brief 	Calculates the square root of the FOM portion of an image.
@param 	*p			image.
@return Bimage* 	square root of the FOM of the input image.

	The new data is calculated as:
		new_datum = sqrt(datum)
	FOM min. must be 0 or greater.

**/
int   img_fom_sqrt(Bimage* p)
{
	if ( !p->next ) {
		cerr << "Error: No FOM block found in " << p->file_name() << endl;
		return -1;
	}

	long   			i, datasize = (long) p->size().volume()*p->images();
	float*          fom = (float *) p->next->data_pointer();

	for ( i=0; i<datasize; i++ )  {
		if ( fom[i] < 0 )  {
			cerr << "Cannot take the square root of the FOM of this image because a value is less than zero." << endl;
			return -1;
		}
		else
			fom[i] = sqrt(fom[i]);
	}

	return 0;
}

/**
@brief 	Applies the Student's t-test to one set of images.
@param 	n			number of files in the set.
@param 	*p			average and variance of the set.
@param 	weight		list of weights.
@return Bimage* 	significant level image (floating point).

	First, average (avg) and variance (var) images are calculated:
			var = [1/(N-1)] sum(x - avgx)^2
	Second, the t value is computed and returned as an image:
			t = avg / sqrt[ var/N ]
	Finally, the significance is calculated and returned as the FOM of the image:
			sig = betai(dof/2, 0.5, dof/(dof + t*t))
	where dof (degrees-of-freedom) is calculated as:
			dof = ws*(1 - 1.0L/N)
	and ws is the weight sum.
	All images are converted to floating point.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.
	Milligan and Flicker (1987) J. Cell Biol. 105:29-39.

**/
Bimage* 	img_students_t_test(int n, Bimage* p, vector<double>& weight)
{
	Bimage*			pt = new Bimage(Float, TSimple, p->size(), 1);
	pt->origin(p->image->origin());
	pt->sampling(p->sampling(0));
	
	float*			var = (float *) p->next->data_pointer();

	pt->next = new Bimage(Float, TSimple, pt->size(), 1);
	float*			siglev = (float *) pt->next->data_pointer();
	
	long  			i, datasize = p->data_size();
	double			ws(0), d;
	
	for ( i=0; i<n; i++ ) ws += weight[i];
	
	double			dof = ws*(1 - 1.0L/n);
	double			hdof = 0.5*dof;
	double			avg_var = 0;
	
	for ( i=0; i<datasize; i++ ) avg_var += var[i];	
	avg_var /= datasize;
	
	if ( verbose ) {
		cout << "Student's t test:" << endl;
		cout << "Number:                         " << n << endl;
		cout << "Average variance:               " << avg_var << endl;
		cout << "Degrees of freedom:             " << dof << endl << endl;
	}

	avg_var *= avg_var;
	
	for ( i=0; i<datasize; i++ ) {
		d = var[i]/n;
		if ( d < avg_var ) d = avg_var;
		d = sqrt(d);
		pt->set(i, (*p)[i]/d);
		siglev[i] = betai(hdof, 0.5, dof/(dof + (*pt)[i]*(*pt)[i]));
	}

	return pt;
}

/**
@author	David Belnap and Bernard Heymann
@brief 	Applies the Student's t-test to two sets of images, assumes the two distributions have approximately the same variance
@param 	set1 		number of files in set 1.
@param 	set2 		number of files in set 2.
@param 	*p1  		average image with variance (FOM), set 1
@param 	*p2 		average image with variance (FOM), set 2
@param 	weight		list of weights.
@return Bimage* 	t-test image and significance image (as FOM).

	Input average images for sets 1 and 2, with variance as FOM portion
	of the images.  The numerator portion of the variance (Vn) is 
	computed from the variance.  The "standard error of the difference 
	of the means" is computed:
	SD = sqrt[ ( (Vn1 + Vn2) / degrees_of_freedom ) * (1/set1 + 1/set2) ]

	and used to compute the t value:
				t = ( avg1 - avg2 ) / SD
	
	Finally, the significance level (probability) is computed.

Reference: 	Press W.H. et al (1992) Numerical Recipes in C.
	Milligan and Flicker (1987) J. Cell Biol. 105:29-39.

**/
Bimage* 	img_students_t_test_equal(int set1, int set2, Bimage* p1, Bimage* p2, vector<double>& weight)
{
	// Set arrays and image parameters
	Bimage*			pt = new Bimage(Float, TSimple, p1->size(), 1);
	pt->origin(p1->image->origin());
	pt->sampling(p1->sampling(0));
	pt->next = new Bimage(Float, TSimple, pt->size(), 1);
	
	long			i, j, datasize = p1->sizeX()*p1->sizeY()*p1->sizeZ()*p1->images();

	float*			var1 = (float *) p1->next->data_pointer();   // pointers to variance data
	float*			var2 = (float *) p2->next->data_pointer();
	float*			sig  = (float *) pt->next->data_pointer();

	// Compute Student's t-test
	if ( verbose )  {
		cout << endl << "Student's t-test, assuming same variance in sets 1 and 2" << endl;
		cout << "Set 1 contains " << set1 << " images" << endl;
		cout << "Set 2 contains " << set2 << " images" << endl;
	}

	double          ws1, ws2;
	for ( i=0, j=0, ws1=0; i<set1; i++, j++ )  ws1 += weight[j];
	for ( i=0, ws2=0; i<set2; i++, j++ )  ws2 += weight[j];
	
	double			inv_weight_sum = (1.0L / ws1) + (1.0L / ws2);

	double          dof1 = ws1*(1.0L - 1.0L/set1);   // degrees of freedom
	double          dof2 = ws2*(1.0L - 1.0L/set2);
	double          dof  = dof1 + dof2;
	double          hdof = 0.5*dof;
	double          SD;
	double			Vn1, Vn2;

	long            zero_count1 = 0, zero_count2 = 0;	// Counts numbers of zero values in t-test
	for ( i=0; i<datasize; i++ )  {
		Vn1 = var1[i] * dof1;				// Compute numerator of variance from the variance
		Vn2 = var2[i] * dof2;
		SD = sqrt(  ( (Vn1 + Vn2) / dof ) * inv_weight_sum  );    // Standard error of the difference of the means
		if ( SD > 0 )  {                               // test for zero in denominator
			pt->set(i, ( (*p1)[i] - (*p2)[i] ) / SD);         // t-test value
			if ( (*pt)[i] == 0 )  zero_count2++;
			sig[i] = betai(hdof, 0.5L, dof/(dof + (*pt)[i]*(*pt)[i]) );   // significance level
		} else  {
			pt->set(i, 0);
			sig[i] = 1;
			zero_count1++;
		}
	}
	
	if ( verbose & VERB_FULL )  {
		cout << zero_count1 << " pixels had a zero value for the standard error of the differences of the means, and t was set to zero." << endl;
		cout << zero_count2 << " pixels had standard deviation > 0 and t = 0." << endl;
		cout << datasize << " = total pixels in image" << endl;
	}

	return pt; 
}



/**
@author	David Belnap
@brief 	Applies the Student's t-test to two sets of images, assumes the two distributions have unequal variance
@param 	set1		number of files in set 1.
@param 	set2		number of files in set 2.
@param 	*p1			average image with variance (FOM), set 1
@param 	*p2			average image with variance (FOM), set 2
@param 	weight		list of weights.
@return Bimage* 	t-test image and significance image (as FOM).

	Input average images for sets 1 and 2, with variance as FOM portion
	of the images.  The t value is computed:
			t = ( avg1 - avg2 ) / sqrt[ var1/ws1 + var2/ws2 ]

	(ws = weighted sum, default = N).  The significance level 
	(probability) is computed.  Output is in floating point.

Reference: 	Press W.H. et al (1992) Numerical Recipes in C.
	Milligan and Flicker (1987) J. Cell Biol. 105:29-39.

**/
Bimage* 	img_students_t_test_unequal(int set1, int set2, Bimage* p1, Bimage* p2, vector<double>& weight)
{
	// Set arrays and image parameters
	Bimage*			pt = new Bimage(Float, TSimple, p1->sizeX(), p1->sizeY(), p1->sizeZ(), 1);
	pt->origin(p1->image->origin());
	pt->sampling(p1->sampling(0));
	pt->next = new Bimage(Float, TSimple, pt->size(), 1);
	
	long			i, j, datasize = p1->sizeX()*p1->sizeY()*p1->sizeZ()*p1->images();

	float*			var1 = (float *) p1->next->data_pointer();   // pointers to variance data
	float*			var2 = (float *) p2->next->data_pointer();
	float*			sig  = (float *) pt->next->data_pointer();	// initialize figure-of-merit subimage structure, for significance level

	// Compute Student's t-test
	if ( verbose )  {
		cout << endl << "Student's t-test, assuming unequal variance in sets 1 and 2" << endl;
		cout << "Set 1 contains " << set1 << " images" << endl;
		cout << "Set 2 contains " << set2 << " images" << endl;
	}

	double          ws1, ws2;
	for ( i=0, j=0, ws1=0; i<set1; i++, j++ )  ws1 += weight[j];
	for ( i=0, ws2=0; i<set2; i++, j++ )  ws2 += weight[j];

	double			dof1 = ws1*(1.0L - 1.0L/set1);   // degrees of freedom
	double			dof2 = ws2*(1.0L - 1.0L/set2);
	double          dof;
	double			Vn1, Vn2;

	long            zero_count1 = 0, zero_count2 = 0;	// Counts numbers of zero values in t-test
	for ( i=0; i<datasize; i++ )  {
		if ( var1[i] > 0 || var2[i] > 0 )  {  // test for zero in denominator
			Vn1 = var1[i]/ws1;
			Vn2 = var2[i]/ws2;
			pt->set(i, ( (*p1)[i] - (*p2)[i] ) / sqrt(Vn1 + Vn2));
			if ( (*pt)[i] == 0 )  zero_count2++;
			dof = (Vn1 + Vn2)*(Vn1 + Vn2) / (Vn1*Vn1/dof1 + Vn2*Vn2/dof2);
			sig[i] = betai(0.5L*dof, 0.5L, dof/(dof + (*pt)[i]*(*pt)[i]) );   // significance level
		} else  {
			pt->set(i, 0);
			sig[i] = 1;
			zero_count1++;
		}
	}

	if ( verbose & VERB_FULL )  {
		cout << zero_count1 << " pixels had zero variance and t was set to zero." << endl;
		cout << zero_count2 << " pixels had variance not equal to 0 and t = 0." << endl;
		cout << datasize << " = total pixels in image" << endl;
	}

	return pt; 
}

/**
@author	David Belnap
@brief 	Applies the F-test to two sets of images
@param 	set1		number of files in set 1.
@param 	set2		number of files in set 2.
@param 	*p1			average image with variance (FOM), set 1
@param 	*p2			average image with variance (FOM), set 2
@param 	weight		list of weights.
@return Bimage*		F-test image and significance image (as FOM).

	The F value is computed for each pixel:
			F(i) = var1(i) / var2(i)

	Variance (as FOM of an average image) is input.  The significance 
	level is computed and stored as FOM of F-test image.  Output is
	in floating point.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
Bimage* 	img_f_test(int set1, int set2, Bimage* p1, Bimage* p2, vector<double>& weight)
{
	// Set image parameters and arrays
	Bimage*			pf = new Bimage(Float, TSimple, p1->size(), 1);
	pf->origin(p1->image->origin());
	pf->sampling(p1->sampling(0));
	pf->next = new Bimage(Float, TSimple, pf->size(), 1);
	
	long			i, j, datasize = p1->sizeX()*p1->sizeY()*p1->sizeZ()*p1->images();

	float*			var1 = (float *) p1->next->data_pointer();   // pointers to variance data
	float*			var2 = (float *) p2->next->data_pointer();
	float*			sig  = (float *) pf->next->data_pointer();

	// Compute F-test
	if ( verbose )  {
		cout << endl << "F-test" << endl;
		cout << "Set 1 contains " << set1 << " images" << endl;
		cout << "Set 2 contains " << set2 << " images" << endl;
	}

	double          ws1, ws2;
	for ( i=0, j=0, ws1=0; i<set1; i++, j++ )  ws1 += weight[j];
	for ( i=0, ws2=0; i<set2; i++, j++ )  ws2 += weight[j];

	double			dof1 = ws1*(1.0L - 1.0L/set1);
	double			dof2 = ws2*(1.0L - 1.0L/set2);
	double          df1, df2;		// degrees of freedom, sets 1 and 2

	for ( i=0; i<datasize; i++ )  {
		if ( var1[i] < SMALLFLOAT ||  var2[i] < SMALLFLOAT )  {
			pf->set(i, 0);
			sig[i] = 1;
			j++;
		} else  {
			if ( var1[i] > var2[i] )  {  // make f the ratio of the larger variance to the smaller one
				pf->set(i, var1[i] / var2[i]);
				df1 = dof1;
				df2 = dof2;
			} else  {
				pf->set(i, var2[i] / var1[i]);
				df1 = dof2;
				df2 = dof1;
			}
			sig[i] = 2*betai(0.5L*df2, 0.5L*df1, df2/(df2 + df1*(*pf)[i]) );   // significance level
			if ( sig[i] > 1 ) sig[i] = 2 - sig[i];
		}
	}
	if ( verbose & VERB_FULL )  {
		cout << j << " pixels in set 1 or set 2 had zero variance, and the F-values were set to zero." << endl;
		cout << datasize << " = total pixels in image" << endl;
	}

	return pf;
}


