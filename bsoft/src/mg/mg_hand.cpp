/**
@file	mg_hand.cpp
@brief	Functions for comparing projections of two oppositely handed reconstructions to images of tilted specimens
@author	David Belnap
@date	Created: 20011003
@date	Modified: 20151008 (BH)
**/


#include "mg_hand.h"
#include "mg_select.h"
#include "Matrix3.h"
#include "qsort_functions.h"
#include "utilities.h"


// Global variables
extern int   verbose;   // Level of output to the screen


/*
@Object: struct Hand_Stats
@brief	Statistical parameters for the hand functions (handedness tilt 
	experiment), and sums needed to compute them.  Also, counts keep 
	score of hand A vs. hand B.
@Features:
**/
struct Hand_Stats {
	int		CountA;     // Count of all particles that matched best with hand A, B, or neither
	int		CountB;
	int		Tie;
	double	Avg1A;      // Average and Standard deviations:
	double	Avg2A;      // control, view 1 hand A (1A); and hands A (2A) and B (2B)
	double	Avg2B;
	double	Std1A;
	double	Std2A;
	double	Std2B;
	double	AvgAbsDiff_2Aphi;   // Absolute differences:  measured angles - predicted angles
	double	AvgAbsDiff_2Athe;
	double	AvgAbsDiff_2Apsi;
	double	AvgAbsDiff_2Bphi;
	double	AvgAbsDiff_2Bthe;
	double	AvgAbsDiff_2Bpsi;
	double	StdAbsDiff_2Aphi;
	double	StdAbsDiff_2Athe;
	double	StdAbsDiff_2Apsi;
	double	StdAbsDiff_2Bphi;
	double	StdAbsDiff_2Bthe;
	double	StdAbsDiff_2Bpsi;
	double	Sum_1A;            // Sums for computing averages and their standard deviations
	double	Sum_2A;            // of correlation coefficients
	double	Sum_2B;
	double	Sum_1Asq;
	double	Sum_2Asq;
	double	Sum_2Bsq;
	double	Sum_AbsDiffPhi2A;    // Sums for computing averages and their standard deviations
	double	Sum_AbsDiffThe2A;    // of the absoluted differences between measured and predicted
	double	Sum_AbsDiffPsi2A;    // angles
	double	Sum_AbsDiffPhi2B;
	double	Sum_AbsDiffThe2B;
	double	Sum_AbsDiffPsi2B;
	double	Sum_AbsDiffPhi2Asq;
	double	Sum_AbsDiffThe2Asq;
	double	Sum_AbsDiffPsi2Asq;
	double	Sum_AbsDiffPhi2Bsq;
	double	Sum_AbsDiffThe2Bsq;
	double	Sum_AbsDiffPsi2Bsq;
};


// Internal function prototypes
int   field_get_handedness(Bimage* mapA, Bimage* mapB, Bmicrograph* mg1, Bmicrograph* mg2, Hand_Stats* Global, 
		double rad_min, double rad_max, double res_min, double res_max, double AmB_min, double AB_min, int diff_out, 
		int origins2, Bstring outimg, int fieldnum);
int   get_handedness_from_tilt_pair(Bimage* mapA, Bimage* mapB, Bmicrograph* mg1, Bmicrograph* mg2,
		Bparticle* particle1, Bparticle* particle2, 
		Hand_Stats* Pair, double rad_min, double rad_max, double res_min, double res_max, double AmB_min, double AB_min, 
		int diff_out, int origins2, Bstring outimg, int fieldnum);

int   hand_compute_stats(Hand_Stats* stats, int diff_out);
int   hand_prepare_image(Bimage* p, double resmin, double resmax);
int   hand_print_outheader(int which_one, int diff_out);
int   hand_print_stats(Hand_Stats* stats, int which_one, int diff_out);
int   hand_write_img(Bimage* img, Bstring imgname, Bstring outimg, int fieldnum, int num);



/**
@brief 	Determines handedness for all selected particles (sel > 0) in a 
	project.
@param	*mapA			3D map (hand A)
@param	*project		Project structure
@param	*mg_ang			angles for micrograph selection, views 1 & 2
@param	*mg_index		indices for micrograph selection, views 1 & 2
@param	*mg_select		micrograph selection criteria, views 1 & 2
@param	rad_min			minimum radius for FOM calculation (pixels)
@param	rad_max			maximum radius for FOM calculation (pixels)
@param	res_min			minimum resolution for FOM calculation (angs.)
@param	res_max			maximum resolution for FOM calculation (angs.)
@param	AmB_min			|FOMA - FOMB| must be this value or greater
@param	AB_min			min. acceptable value for FOM of correct hand
@param	diff_out		output difference of measured & predicted orientations
@param	origins2		flag to determine origins for second view
@param	outimg			prefix & suffix for output projection files
@return int				0.

	Gets tilt-axis direction and rotation angle for a micrograph pair.
	Then, loops through the selected particles in pair.  Sends data to
	function get_handedness_from_tilt_pair, where handedness is
	determined.  Statistics are calculated for each pair and sums for
	global statistics are tabulated.

**/
int  project_get_handedness(Bimage* mapA, Bproject* project, double* mg_ang, 
			int* mg_index, int* mg_select, double rad_min, double rad_max,
			double res_min, double res_max, double AmB_min, double AB_min, 
			int diff_out, int origins2, Bstring outimg)
{
	project->fom_tag[1] = FOM_HAND_A;
	project->fom_tag[2] = FOM_HAND_B;

	Bfield*      field=NULL;			// Placeholder for loop over the fields-of-view in project
	int          handedness;			// Output to calling function:  1=hand A, 2=hand B, 0=undetermined
	int          ifield;				// index for field block
	Bmicrograph  *mg1=NULL, *mg2=NULL;	// Structures to hold input data for first and second views of a micrograph tilt pair
	int          n_Global;				// Number of particle-image pairs analyzed, for statistical calculations
	Bimage       *mapB=NULL;			// Mirror image of input map (mapA)

	Hand_Stats*  Global = new Hand_Stats;	// Global statistics (i.e. for all micrographs)
	memset(Global, 0, sizeof(Hand_Stats));	// zero sums for global statistics

	if ( verbose )  {
		cout << "Minimum Radius (pixels):        " << rad_min << endl;
		cout << "Maximum Radius (pixels):        " << rad_max << endl;
		cout << "High-Pass filter (angs.):       " << res_min << endl;
		cout << "Low-Pass filter (angs.):        " << res_max << endl << endl;
	}


	// Copy enantiomer and invert it through origin where equation is
	//                  mirror(x,y,z) = enantiomer(-x,-y,-z)
	// if origin is given as (0,0,0).  mapB is hand B map.
	mapA->calculate_background();  // first be sure background is set correctly
	Vector3<double>		scale(-1,-1,-1), translate;
	Matrix3				mat(1);
//	mapB = img_transform_copy(mapA, mapA->size(), scale, mapA->image->origin(), translate, mat, FILL_BACKGROUND, 0);
	mapB = mapA->transform(mapA->size(), scale, mapA->image->origin(), translate, mat, FILL_BACKGROUND, 0);

	if (verbose)  {
		cout << "Hand A map origin:              " << mapA->image->origin() << endl;
		cout << "Hand B map origin:              " << mapB->image->origin() << endl;
		cout << "Map Sampling (angs./pixel):     " << mapA->sampling(0) << endl << endl;
	}


	// Determine handedness
	for (  field = project->field, ifield=1;  field;  field = field->next, ifield++  )  {
		// select first view
		mg1 = field_find_micrograph(field, mg_select[0], mg_index[0], mg_ang[0]);
		// select second view
		mg2 = field_find_micrograph(field, mg_select[1], mg_index[1], mg_ang[1]);
		// find handedness for a field
		field_get_handedness(mapA, mapB, mg1, mg2, Global, rad_min, rad_max,   
			res_min, res_max, AmB_min, AB_min, diff_out, origins2, outimg, ifield);
	}


	// Calculate and print global statistics to screen
	hand_compute_stats(Global, diff_out);
	if (verbose)  {
		cout << endl << endl;
		hand_print_stats(Global, 1, diff_out);
	}


	// If 95% of particles agree with one handedness, handedness=1 if hand A or 2 if hand B.  0 if < 95%.
	n_Global = Global->CountA + Global->CountB + Global->Tie;
	if ( (Global->CountA/n_Global) >= 0.95 )  {
		handedness = 1;
		if (verbose)  cout << "Output map has hand A (same as input map)." << endl;
	} else if ( (Global->CountB/n_Global) >= 0.95 )  {
		handedness = 2;
		if (verbose)  cout << "Output map has hand B (opposite of input map)." << endl;
	} else  {
		handedness = 0;
		if (verbose) {
			cout << "Output map is undetermined because the number of selected" << endl;
			cout << "particles favoring hand A or B was not 95% or greater." << endl;
		}
	}

	
	// Memory clean-up
	delete mapB;
	delete Global;


	return handedness;
}



/*
@brief 	Determines handedness for all selected particles (sel > 0) for a 
	field--specifically a selected tilt pair of micrographs from the 
	field.
@param 	*mapA			3D map (hand A)
@param 	*mapB			3D map (hand B)
@param 	*mg1           	Micrograph for view 1 of field
@param 	*mg2           	Micrograph for view 2 of field
@param 	*Global        	Statistics for all fields-of-view
@param	rad_min       	minimum radius for FOM calculation (pixels)
@param	rad_max       	maximum radius for FOM calculation (pixels)
@param	res_min       	minimum resolution for FOM calculation (angs.)
@param	res_max       	maximum resolution for FOM calculation (angs.)
@param	AmB_min       	|FOMA - FOMB| must be this value or greater
@param	AB_min        	min. acceptable value for FOM of correct hand
@param	diff_out      	output difference of measured & predicted orientations
@param	origins2      	flag to determine origins for second view
@param	outimg        	prefix & suffix for output projection files
@param	fieldnum      	index of current field (for image output)
@return int				error code.

	Gets tilt-axis direction and rotation angle for a micrograph pair.
	Then, loops through the selected particles in pair.  Sends data to
	function get_handedness_from_tilt_pair, where handedness is
	determined.  Statistics are calculated for each pair and sums for
	global statistics are tabulated.

**/
int   field_get_handedness(Bimage* mapA, Bimage* mapB, Bmicrograph* mg1, Bmicrograph* mg2, Hand_Stats* Global, 
		double rad_min, double rad_max, double res_min, double res_max, double AmB_min, double AB_min, int diff_out, 
		int origins2, Bstring outimg, int fieldnum)
{

	Bparticle*    particle1 = NULL;   // particles being processed
	Bparticle*    particle2 = NULL;

	Hand_Stats*   Pair = new Hand_Stats;   // Statistics for current micrograph pair
	memset(Pair, 0, sizeof(Hand_Stats));   // zero sums for pair statistics

	// Get the pair's tilt-axis direction and the rotation angle from the mg1 and mg2 structures
	if ( mg1->tilt_axis != mg2->tilt_axis )  {
		error_show("Error in field_get_handedness", __FILE__, __LINE__);
		cerr << "The tilt-axis direction is not the same or is invalid in this pair:" << endl;
		cerr << mg1->fpart << "    axis dir = " << mg1->tilt_axis*180.0/M_PI << endl;
		cerr << mg2->fpart << "    axis dir = " << mg2->tilt_axis*180.0/M_PI << endl;
		cerr << "You may have a pair mismatch or may need to enter a value for both micrographs." << endl;
	    return -1;
	}
	
	if ( verbose ) {
		cout << "Tilt pair:" << endl;
		cout << "Micrograph 1:                   " << mg1->fpart << endl;
		cout << "Micrograph 2:                   " << mg2->fpart << endl;
		cout << "Tilt-axis direction:            " << mg1->tilt_axis*(180/M_PI) << endl;
		cout << "Rotation angle:                 " << (mg2->tilt_angle - mg1->tilt_angle)*(180/M_PI) << endl << endl;
	}

	// Loop through all particles in first micrograph
	//   -It is assumed that the ID numbers in the first and second micrographs refer to the same imaged object.
	//   -To be analyzed, a particle must
	//      a) have the identical number in second micrograph
	//		b) sel > 0
	if (verbose)  hand_print_outheader(0, diff_out);   // Header for screen output
	for ( particle1=mg1->part;  particle1;  particle1=particle1->next ) if ( particle1->sel )  {
		// loop through particles in second mg to find identical particle ID
		for ( particle2=mg2->part;  particle2 && particle1->id != particle2->id;  particle2=particle2->next ) ;  
		if ( particle2 )
			get_handedness_from_tilt_pair(mapA, mapB, mg1, mg2, particle1, particle2, Pair, 
					rad_min, rad_max, res_min, res_max, AmB_min, AB_min, diff_out, origins2, outimg, fieldnum);
	}


	// Calculate and output pair statistics
	hand_compute_stats(Pair, diff_out);
	if ( verbose )
		hand_print_stats(Pair, 0, diff_out);

	// Tabulate sums for global statistics
	Global->Sum_1A   += Pair->Sum_1A;
	Global->Sum_2A   += Pair->Sum_2A;
	Global->Sum_2B   += Pair->Sum_2B;
	Global->Sum_1Asq += Pair->Sum_1Asq;
	Global->Sum_2Asq += Pair->Sum_2Asq;
	Global->Sum_2Bsq += Pair->Sum_2Bsq;
	Global->CountA   += Pair->CountA;
	Global->CountB   += Pair->CountB;
	Global->Tie      += Pair->Tie;
	if (diff_out > 0)  {
		Global->Sum_AbsDiffPhi2A   += Pair->Sum_AbsDiffPhi2A;
		Global->Sum_AbsDiffThe2A   += Pair->Sum_AbsDiffThe2A;
		Global->Sum_AbsDiffPsi2A   += Pair->Sum_AbsDiffPsi2A;
		Global->Sum_AbsDiffPhi2B   += Pair->Sum_AbsDiffPhi2B;
		Global->Sum_AbsDiffThe2B   += Pair->Sum_AbsDiffThe2B;
		Global->Sum_AbsDiffPsi2B   += Pair->Sum_AbsDiffPsi2B;
		Global->Sum_AbsDiffPhi2Asq += Pair->Sum_AbsDiffPhi2Asq;
		Global->Sum_AbsDiffThe2Asq += Pair->Sum_AbsDiffThe2Asq;
		Global->Sum_AbsDiffPsi2Asq += Pair->Sum_AbsDiffPsi2Asq;
		Global->Sum_AbsDiffPhi2Bsq += Pair->Sum_AbsDiffPhi2Bsq;
		Global->Sum_AbsDiffThe2Bsq += Pair->Sum_AbsDiffThe2Bsq;
		Global->Sum_AbsDiffPsi2Bsq += Pair->Sum_AbsDiffPsi2Bsq;
	}


	delete Pair;

	return 0;

}  // End of function field_get_handedness



/*
@brief	Determines handedness from tilt pair of images, (two views of the
	same specimen in which the specimen was rotated about a known
	axis and by a known angle).

@param 	*mapA			3D map (hand A)
@param 	*mapB			3D map (hand B)
@param	*mg1			Micrograph for view 1
@param	*mg2			Micrograph for view 2
@param 	*particle1		Particle from view 1
@param 	*particle2		Particle from view 2
@param	*Pair			Statistics for current micrograph pair
@param 	rad_min			minimum radius for FOM calculation (pixels)
@param 	rad_max			maximum radius for FOM calculation (pixels)
@param 	res_min			minimum resolution for FOM calculation (angs.)
@param 	res_max			maximum resolution for FOM calculation (angs.)
@param 	AmB_min			|FOMA - FOMB| must be this value or greater
@param 	AB_min			min. acceptable value for FOM of correct hand
@param	diff_out		output difference of measured & predicted orientations
@param	origins2		flag to determine origins for second view
@param	outimg			prefix & suffix for output projection files
@param	fieldnum		index of current field (for image output)
@return	int				0.

	-Calculate two predicted orientations for second view from first 
	 orientation, the tilt-axis direction, and the rotation angle.  One 
	 orientation is the predicted orientation for hand A (handedness 
	 of input 3D map and first orientation) in the second view.  The 
	 other is the predicted orientation for hand B (mirror of hand A 
	 map).
	-Projections made of input map at first orientation and predicted 
	 second orientation of hand A.
	-Projection made of mirror map at predicted second orientation of
	 hand B.
	-Projections of predicted second orientations are compared to image
	 of second view.
	-As a control, projection of first orientation is compared to image
	 of first view.
**/
int  get_handedness_from_tilt_pair(Bimage* mapA, Bimage* mapB, Bmicrograph* mg1, Bmicrograph* mg2,
		Bparticle* particle1, Bparticle* particle2, 
		Hand_Stats* Pair, double rad_min, double rad_max, double res_min, double res_max, double AmB_min, double AB_min, 
		int diff_out, int origins2, Bstring outimg, int fieldnum)
{
	double  axis_dir(mg1->tilt_axis);	// Angle of the tilt axis on the micrograph
	double  rotation_angle(mg2->tilt_angle - mg1->tilt_angle);	// Angle of rotation about the tilt axis = tilt_angle (view2) - tilt_angle (view1)

	double  max_radius;       	 // radius of largest sphere that will fit inside the 3D map
	double  cc1A;             	 // Correlation coefficient for first view and predicted second views
	double  cc2A, cc2B;       	 // model projections vs. particle image
	double  diff2A_phi(0);       	 // measured angle - predicted angle
	double  diff2A_the(0);
	double  diff2A_psi(0);
	double  diff2B_phi(0);
	double  diff2B_the(0);
	double  diff2B_psi(0);
	int    ThresholdOK(1);   	 // flag for |A-B| threshold

	Vector3<double>  translate1, translate2A, translate2B;  // Translational-shift arrays for particle images
	Vector3<double>  translate0;                            // for cross-correlation option

	Euler  angles1A;			// Euler angles for first view
	Euler  angles2A, angles2B;	// Euler angles for second view, hands A & B
	Euler  euler2;				// Euler angles for second view (input data)
	View	view1A(particle1->view);	// View vector and angle for first view, hand A
	View	view1B(view1A);
	view1B.angle(angle_set_negPI_to_PI(view1B.angle()+M_PI));
	View	view2;				// View vector and angle for second view (input data)

	Matrix3 rot_matrix1A(1);    // 3x3 matrix for first and predicted second orientations,
	Matrix3 rot_matrix1B(1);    // for hands A and B
	Matrix3 rot_matrix2A(1);
	Matrix3 rot_matrix2B(1);
	Matrix3 rot_matrixtilt(1);  // 3x3 matrix for the tilt axis

	Bimage*   img1  = NULL;  // a particle image from first micrograph
	Bimage*   img2A = NULL;  // a particle image from second micrograph
	Bimage*   img2B = NULL;  // a copy of img2A, but it may have a different origin
	Bimage*   prj1A = NULL;  // a projection image corresponding to the particle image from first micrograph
	Bimage*   prj2A = NULL;  // a projection image from hand A map, corresponding to the particle image from second micrograph
	Bimage*   prj2B = NULL;  // a projection image from hand B map, corresponding to the particle image from second micrograph


	// General parameters
	max_radius = mapA->maximum_included_radius();  // reduces unnecessary computations in rotate_project and img_find_shift


	// Compute predicted orientations for second view (for both hand A and hand B)
	//   -The addition of Pi for the 1B orientation must be added to the psi 1A value, it cannot be added to the
	//    tilt-axis direction.  This is because the tilt axis does not change--and neither does the first view.
	//    Only the psi value for the first view changes.
	rot_matrix1A = view1A.matrix();  // rotation matrices for views 1A & 1B
	rot_matrix1B = view1B.matrix();

	rot_matrixtilt = Matrix3(rotation_angle,angle_set_negPI_to_PI(axis_dir));  // rotation matrix for tilting operation

	rot_matrix2A = rot_matrix1A * rot_matrixtilt;  // rotation matrices for predicted views 2A & 2B
	rot_matrix2B = rot_matrix1B * rot_matrixtilt;

	// Get orientation angles
	angles1A = Euler(view1A);
	angles2A = Euler(rot_matrix2A);  // get predicted angles from matrices
	angles2B = Euler(rot_matrix2B);

	if ( verbose & VERB_PROCESS )  {
		cout << "Input orientation (phi, theta, psi)        = " << 
			angles1A.phi()*(180/M_PI) << " " << angles1A.theta()*(180/M_PI) << " " << angles1A.psi()*(180/M_PI) << endl;
		cout << "Predicted orientation 2A (phi, theta, psi) = " << 
			angles2A.phi()*(180/M_PI) << " " << angles2A.theta()*(180/M_PI) << " " << angles2A.psi()*(180/M_PI) << endl;
		cout << "Predicted orientation 2B (phi, theta, psi) = " << 
			angles2B.phi()*(180/M_PI) << " " << angles2B.theta()*(180/M_PI) << " " << angles2B.psi()*(180/M_PI) << endl;
	}
	if ( verbose & VERB_FULL )  {
		cout << "Matrix for first view, hand A" << endl;
		cout << rot_matrix1A << endl;
		cout << "Matrix for tilt" << endl;
		cout << rot_matrixtilt << endl;
		cout << "Matrix for second view, hand A" << endl;
		cout << rot_matrix2A << endl;
		cout << "Matrix for second view, hand B" << endl;
		cout << rot_matrix2B << endl;
	}


	// Read images and convert to double point, set pixel size to value from micrograph structure (if different)
	img1  = read_img(mg1->fpart,1,(particle1->id - 1));
	img1->change_type(Float);
	img1->sampling(mg1->pixel_size);
	img1->origin(particle1->ori);
	
	img2A = read_img(mg2->fpart,1,(particle2->id - 1));
	img2A->change_type(Float);
	img2A->sampling(mg2->pixel_size);
	img2A->origin(particle2->ori);

	img2B = img2A->copy();  // copy made for hand 2B test because origin may differ

	Vector3<double>		oriA(mapA->image->origin()[0], mapA->image->origin()[1], 0);
	Vector3<double>		oriB(mapB->image->origin()[0], mapB->image->origin()[1], 0);

	// Get translational offset (x & y only, z offset = 0) and set origins of data images
	//   -The translational offset is the origin of the particle image minus the center of the map array
	translate1 = particle1->ori - oriA;

	if (origins2 == 0)  {   // Get image origins (second view) from input parameter file
		translate2A = particle2->ori - oriA;
		translate2B = particle2->ori - oriB;
	} else if (origins2 == 1)  {  // Determine origins of second view via cross-correlation with projections
		prj2A = mapA->rotate_project(rot_matrix2A, translate0, max_radius);  // project with no shifts
		prj2B = mapB->rotate_project(rot_matrix2B, translate0, max_radius);

//		img_find_shift(img2A, prj2A, NULL, res_max, res_min, max_radius, 0, 1);    // hand A
		img2A->find_shift(prj2A, NULL, res_max, res_min, max_radius, 0, 1);    // hand A
		translate2A = prj2A->image->origin() - img2A->image->origin();
		img2A->image->origin(translate2A + oriA);

//		img_find_shift(img2B, prj2B, NULL, res_max, res_min, max_radius, 0, 1);    // hand B
		img2B->find_shift(prj2B, NULL, res_max, res_min, max_radius, 0, 1);    // hand B
		translate2B = prj2B->image->origin() - img2B->image->origin();
//		img2B->image->origin(translate2B + oriB);
		img2B->origin(translate2B + oriB);

		delete prj2A;
		delete prj2B;
	} else if (origins2 == 2)  {   // Get image origins (second view) from Bsubimage structure
		translate2A = img2A->image->origin() - oriA;
		translate2B = img2B->image->origin() - oriB;
	}


	// Compute projections for 1A, 2A, and 2B orientations, origins put at respective image origin
	prj1A = mapA->rotate_project(rot_matrix1A, translate1,  max_radius);
	prj2A = mapA->rotate_project(rot_matrix2A, translate2A, max_radius);
	prj2B = mapB->rotate_project(rot_matrix2B, translate2B, max_radius);


	// Write projections to disk, if desired
	if ( outimg.length() )  {
		hand_write_img(prj1A, "prj1a", outimg, fieldnum, particle1->id);
		hand_write_img(prj2A, "prj2a", outimg, fieldnum, particle1->id);
		hand_write_img(prj2B, "prj2b", outimg, fieldnum, particle1->id);
	}


	//  Projections have same origin and pixel size as images
//	prj1A->image->origin(img1->image->origin());
//	prj2A->image->origin(img2A->image->origin());
//	prj2B->image->origin(img2B->image->origin());
	prj1A->origin(0, img1->image->origin());
	prj2A->origin(0, img2A->image->origin());
	prj2B->origin(0, img2B->image->origin());

	prj1A->sampling(img1->sampling(0));
	prj2A->sampling(img2A->sampling(0));
	prj2B->sampling(img2B->sampling(0));


	//  Prepare images for comparison test by setting statistics and applying a bandpass filter
	hand_prepare_image(img1, res_min, res_max);
	hand_prepare_image(img2A,res_min, res_max);
	hand_prepare_image(img2B,res_min, res_max);
	hand_prepare_image(prj1A,res_min, res_max);
	hand_prepare_image(prj2A,res_min, res_max);
	hand_prepare_image(prj2B,res_min, res_max);


	// Debugging output
	if ( outimg.length() && (verbose & VERB_DEBUG) )  {
		hand_write_img(prj1A, "prj1a", outimg, fieldnum, particle1->id);
		hand_write_img(prj2A, "prj2a", outimg, fieldnum, particle1->id);
		hand_write_img(prj2B, "prj2b", outimg, fieldnum, particle1->id);
		hand_write_img(img1,  "img1",  outimg, fieldnum, particle1->id);
		hand_write_img(img2A, "img2a", outimg, fieldnum, particle1->id);   // img2a and img2b will appear identical
		hand_write_img(img2B, "img2b", outimg, fieldnum, particle1->id);   // but their origins could be different
	}
	if ( verbose & VERB_DEBUG )  {
		cout << "General parameters" << endl;
		cout << "  radius min,max:      " <<  rad_min << " " << rad_max << endl;
		cout << "  resolution min,max:  " <<  res_min << " " << res_max << endl;
		cout << "  micrograph 1,2:      " << mg1->id << " " << mg2->id << endl;
		cout << "  particle   1,2:      " << particle1->id << " " << particle2->id << endl;
		cout << "Image 1 parameters" << endl;
		cout << "  background:          " << img1->background(long(0)) << endl;
		cout << "  origin (x,y,z):      " << img1->image->origin() << endl;
		cout << "  pixel size (x,y,z):  " << img1->sampling(0) << endl;
		cout << "Projection 1A parameters" << endl;
		cout << "  background:          " << img2A->background(long(0)) << endl;
		cout << "  origin (x,y,z):      " << prj1A->image->origin() << endl;
		cout << "  pixel size (x,y,z):  " << prj1A->sampling(0) << endl;
		cout << "Image 2A parameters" << endl;
		cout << "  background:          " << img2B->background(long(0)) << endl;
		cout << "  origin (x,y,z):      " << img2A->image->origin() << endl;
		cout << "  pixel size (x,y,z):  " << img2A->sampling(0) << endl;
		cout << "Projection 2A parameters" << endl;
		cout << "  background:          " << prj1A->background(long(0)) << endl;
		cout << "  origin (x,y,z):      " << prj2A->image->origin() << endl;
		cout << "  pixel size (x,y,z):  " << prj2A->sampling(0) << endl;
		cout << "Image 2B parameters" << endl;
		cout << "  background:          " << prj2A->background(long(0)) << endl;
		cout << "  origin (x,y,z):      " << img2B->image->origin() << endl;
		cout << "  pixel size (x,y,z):  " << img2B->sampling(0) << endl;
		cout << "Projection 2B parameters" << endl;
		cout << "  background:          " << prj2B->background(long(0)) << endl;
		cout << "  origin (x,y,z):      " << prj2B->image->origin() << endl;
		cout << "  pixel size (x,y,z):  " << prj2B->sampling(0) << endl;
	}


	// Compare images to projections
//	cc1A = img_correlation_coefficient_radial(img1,  prj1A, rad_min, rad_max);
//	cc2A = img_correlation_coefficient_radial(img2A, prj2A, rad_min, rad_max);
//	cc2B = img_correlation_coefficient_radial(img2B, prj2B, rad_min, rad_max);
	cc1A = img1->correlate(prj1A, rad_min, rad_max, NULL, 0);
	cc2A = img2A->correlate(prj2A, rad_min, rad_max, NULL, 0);
	cc2B = img2B->correlate(prj2B, rad_min, rad_max, NULL, 0);
	
	particle1->fom[1] = cc1A;
	particle2->fom[1] = cc2A;
	particle2->fom[2] = cc2B;


	// Print CCs to screen, check thresholds, selection option:
	//     -If 1) |difference| of CC2A-CC2B or 2) largest CC (CC2A or CC2B)
	//      is not equal to or greater than their respective thresholds,
	//      discard this particle from statistics and deselected it
	//     -User may also set selection numbers for hand A (1) or B (2)
	if (verbose) cout << setprecision(4) << particle1->id << tab << cc1A << tab << cc2A << tab << cc2B << tab << cc2A-cc2B;

	if ( fabs(cc2A-cc2B) >= AmB_min )  {                         // CC difference > difference threshold
		if ( (cc2A < AB_min) && (cc2B < AB_min) )  {  // both CCs less than CC threshold
			ThresholdOK = 0;
		} else  {                                                               // at least one CC is >= CC threshold
			if ( cc2A > cc2B )  {   // Hand A
				particle1->sel = 1;
				particle2->sel = 1;
			} else  {  // Hand B
				particle1->sel = 2;
				particle2->sel = 2;
			}
			if ( verbose && (diff_out > 0) )  {   // No Deselections:  print predicted orientations to screen
				cout << setprecision(2) << tab << (180/M_PI)*angles2A.phi() << tab << (180/M_PI)*angles2A.theta() << tab << (180/M_PI)*angles2A.psi();
				cout << tab << (180/M_PI)*angles2B.phi() << tab << (180/M_PI)*angles2B.theta() << tab << (180/M_PI)*angles2B.psi();
			}
		}
	} else  {
		ThresholdOK = 0;
	}
	
	if ( !ThresholdOK ) {
		particle1->sel = 0;
		particle2->sel = 0;
//		if (verbose)  cout << "\t|A - B| < threshold, deselected (0)";
	}
	if ( verbose && !diff_out )  cout << tab << ThresholdOK;
	if (verbose) cout << endl;


	// Compute difference between measured and predicted angles if desired and if |A-B| threshold is okay
	//   -absolute difference will be tallied for pair and global statistics
	if ( (ThresholdOK == 1) && (diff_out > 0) )  {

		view2  = particle2->view;
		euler2 = Euler(view2);
//		euler2.check();
//		angles2A.check();
//		angles2B.check();

		diff2A_phi = angle_set_negPI_to_PI(euler2.phi() - angles2A.phi());  // use angle_set_negPI_to_PI to get min. difference
		diff2A_the = euler2.theta() - angles2A.theta();
		diff2A_psi = angle_set_negPI_to_PI(euler2.psi() - angles2A.psi());
		diff2B_phi = angle_set_negPI_to_PI(euler2.phi() - angles2B.phi());
		diff2B_the = euler2.theta() - angles2B.theta();
		diff2B_psi = angle_set_negPI_to_PI(euler2.psi() - angles2B.psi());

		if ( verbose && (diff_out >= 2) )  {
			cout << "              Measured - Predicted  =";
			cout << " " << (180/M_PI)*diff2A_phi;
			cout << " " << (180/M_PI)*diff2A_the;
			cout << " " << (180/M_PI)*diff2A_psi;
			cout << " " << (180/M_PI)*diff2B_phi;
			cout << " " << (180/M_PI)*diff2B_the;
			cout << " " << (180/M_PI)*diff2B_psi << endl << endl;
		}
	}


	// Tabulate sums for pair statistics if |A-B| threshold is okay
	if (ThresholdOK == 1)  {
		Pair->Sum_1A   += cc1A;
		Pair->Sum_2A   += cc2A;
		Pair->Sum_2B   += cc2B;
		Pair->Sum_1Asq += cc1A*cc1A;
		Pair->Sum_2Asq += cc2A*cc2A;
		Pair->Sum_2Bsq += cc2B*cc2B;
		if (cc2A > cc2B)			  // Count to see which handedness wins the contest
			Pair->CountA++;
		else if (cc2B > cc2A)
			Pair->CountB++;
		else
			Pair->Tie++;
		if (diff_out > 0)  {
			Pair->Sum_AbsDiffPhi2A   += fabs(diff2A_phi);
			Pair->Sum_AbsDiffThe2A   += fabs(diff2A_the);
			Pair->Sum_AbsDiffPsi2A   += fabs(diff2A_psi);
			Pair->Sum_AbsDiffPhi2B   += fabs(diff2B_phi);
			Pair->Sum_AbsDiffThe2B   += fabs(diff2B_the);
			Pair->Sum_AbsDiffPsi2B   += fabs(diff2B_psi);
			Pair->Sum_AbsDiffPhi2Asq += diff2A_phi*diff2A_phi;
			Pair->Sum_AbsDiffThe2Asq += diff2A_the*diff2A_the;
			Pair->Sum_AbsDiffPsi2Asq += diff2A_psi*diff2A_psi;
			Pair->Sum_AbsDiffPhi2Bsq += diff2B_phi*diff2B_phi;
			Pair->Sum_AbsDiffThe2Bsq += diff2B_the*diff2B_the;
			Pair->Sum_AbsDiffPsi2Bsq += diff2B_psi*diff2B_psi;
		}
	}


	// Memory cleanup
	delete img1;
	delete img2A;
	delete img2B;
	delete prj1A;
	delete prj2A;
	delete prj2B;
	
	return 0;


}  // End of function get_handedness_from_tilt_pair



/*
@brief	Computes statistics from data collected in Hand_Stats structure.
@param	*stats		Handedness statistical parameters
@param	diff_out	Switch for measured-predicted angle output
@return	int 		0.

	Computes average and standard deviation for view1 projection vs.
	image1, handA-view2 projection vs. image2, and handB-view2 projection
	vs. image2 correlation coefficients.  Also computes these parameters
	for differences between measured and predicted angles, if user
	desires.
**/
int		hand_compute_stats(Hand_Stats* stats, int diff_out)
{
	int			N;   // Number of data points
	double	d;

	N = stats->CountA + stats->CountB + stats->Tie;   // total number of selected particles
	if (N == 0)  return -1;                          // no data = no statistics

	// Averages
	stats->Avg1A = stats->Sum_1A/N;
	stats->Avg2A = stats->Sum_2A/N;
	stats->Avg2B = stats->Sum_2B/N;

	if (diff_out > 0)  {            // Absolute differences:  measured angles - predicted angles
		stats->AvgAbsDiff_2Aphi = stats->Sum_AbsDiffPhi2A/N;
		stats->AvgAbsDiff_2Athe = stats->Sum_AbsDiffThe2A/N;
		stats->AvgAbsDiff_2Apsi = stats->Sum_AbsDiffPsi2A/N;
		stats->AvgAbsDiff_2Bphi = stats->Sum_AbsDiffPhi2B/N;
		stats->AvgAbsDiff_2Bthe = stats->Sum_AbsDiffThe2B/N;
		stats->AvgAbsDiff_2Bpsi = stats->Sum_AbsDiffPsi2B/N;
	}

	// Standard deviations
	if (N < 2)  {      // prevents divide-by-zero errors if only one data point
		stats->Std1A = 0;    stats->Std2A = 0;    stats->Std2B = 0;
		stats->StdAbsDiff_2Aphi = 0;    stats->StdAbsDiff_2Athe = 0;    stats->StdAbsDiff_2Apsi = 0;
		stats->StdAbsDiff_2Bphi = 0;    stats->StdAbsDiff_2Bthe = 0;    stats->StdAbsDiff_2Bpsi = 0;
	} else  {
		stats->Std1A = sqrt( (stats->Sum_1Asq - (stats->Sum_1A*stats->Sum_1A)/N) / (N - 1) );
		stats->Std2A = sqrt( (stats->Sum_2Asq - (stats->Sum_2A*stats->Sum_2A)/N) / (N - 1) );
		stats->Std2B = sqrt( (stats->Sum_2Bsq - (stats->Sum_2B*stats->Sum_2B)/N) / (N - 1) );
		if (diff_out > 0)  {
			d = stats->Sum_AbsDiffPhi2Asq - (stats->Sum_AbsDiffPhi2A*stats->Sum_AbsDiffPhi2A)/N;
			if ( d > 0 ) stats->StdAbsDiff_2Aphi = sqrt( d / (N-1) );
			else stats->StdAbsDiff_2Aphi = 0;
			d = stats->Sum_AbsDiffPhi2Asq - (stats->Sum_AbsDiffThe2Asq*stats->Sum_AbsDiffThe2A)/N;
			if ( d > 0 ) stats->StdAbsDiff_2Athe = sqrt( d / (N-1) );
			else stats->StdAbsDiff_2Athe = 0;
			d = stats->Sum_AbsDiffPsi2Asq - (stats->Sum_AbsDiffPsi2A*stats->Sum_AbsDiffPsi2A)/N;
			if ( d > 0 ) stats->StdAbsDiff_2Apsi = sqrt( d / (N-1) );
			else stats->StdAbsDiff_2Apsi = 0;
			d = stats->Sum_AbsDiffPhi2Bsq - (stats->Sum_AbsDiffPhi2B*stats->Sum_AbsDiffPhi2B)/N;
			if ( d > 0 ) stats->StdAbsDiff_2Bphi = sqrt( d / (N-1) );
			else stats->StdAbsDiff_2Bphi = 0;
			d = stats->Sum_AbsDiffThe2Bsq - (stats->Sum_AbsDiffThe2B*stats->Sum_AbsDiffThe2B)/N;
			if ( d > 0 ) stats->StdAbsDiff_2Bthe = sqrt( d / (N-1) );
			else stats->StdAbsDiff_2Bthe = 0;
			d = stats->Sum_AbsDiffPsi2Bsq - (stats->Sum_AbsDiffPsi2B*stats->Sum_AbsDiffPsi2B)/N;
			if ( d > 0 ) stats->StdAbsDiff_2Bpsi = sqrt( d / (N-1) );
			else stats->StdAbsDiff_2Bpsi = 0;
		}
	}

	return 0;
}	



/*
@brief	Prepares image for comparison with another image
@param 	*p			image.
@param 	resmin		Lower resolution cut-off value.
@param 	resmax		Higher resolution cut-off value.
@return	int			0.

	Updates image background and statistics, and applies a bandpass 
	filter to the image.
**/
int		hand_prepare_image(Bimage* p, double resmin, double resmax)
{

	p->statistics();
	p->fspace_bandpass(resmax, resmin, 0);
	p->correct_background();
	p->subtract_background();

	return 0;

}



/*
@brief	Prints the output header to the screen.
@param	which_one	switch for Pair or Global header
@param	diff_out	switch for measured-predicted angle output
@return	int     	0.

	Two headers:  pair and global
**/
int 	hand_print_outheader(int which_one, int diff_out)
{
	cout << "Particle images vs. projections at predicted orientations" << endl << endl;
	cout << "\tCorrelation Coefficients";
	if (which_one == 0)  {   /* for Pair statistics */
		if (diff_out>0)  cout << "\tPredicted Orientations";
		cout << endl << endl;
		cout << "\tprj.\tprj.\tprj." << endl;
		cout << "\tview1\tview2\tview2\tprjs" << endl;
		cout << "\thandA\thandA\thandB\tvs." << endl;
		cout << "Part.\tvs.\tvs.\tvs.\timage2";
		if (diff_out>0)  cout << "\tview2, handA\tview2, handB";
		cout << endl;
		cout << "Id\timage1\timage2\timage2\tA - B";
		if (diff_out>0)  cout << "\tPhi\tTheta\tPsi\tPhi\tTheta\tPsi";
		else cout << "\tselected";
		cout << endl;
		cout << "================================================";
		if (diff_out>0)  cout << "======================================";
		cout << endl;
	}
	else  {                  /* for Global statistics */
		if (diff_out>0)  cout << "\t|Measured Orientation -";
		if (diff_out>0)  cout << endl << "\tPredicted Orientation|" << endl;
		else cout << endl << endl;
		cout << "\tprj.\t prj.\tprj." << endl;
		cout << "\tview1\tview2\tview2" << endl;
		cout << "\thandA\thandA\thandB" << endl;
		cout << "\tvs.\tvs.\tvs.";
		if (diff_out>0)  cout << "\tview2, handA\tview2, handB";
		cout << endl;
		cout << "Global\timage1\timage2\timage2";
		if (diff_out>0)  cout << "\tPhi\tTheta\tPsi\tPhi\tTheta\tPsi";
		cout << endl;
		cout << "====================================";
		if (diff_out>0)  cout << "============================================";
		cout << endl;
	}
	return 0;
}	



/*
@brief	Prints statistics to the screen.
@param	*stats       	Handedness statistical parameters
@param	which_one     	Switch for Pair or Global statistics
@param	diff_out      	Switch for measured-predicted angle output
@return	int           	0.
**/
int 	hand_print_stats(Hand_Stats* stats, int which_one, int diff_out)
{
	double  ToDeg = 180/M_PI;

	cout << setprecision(4);
	
	if (which_one == 0)  {   /* for Pair statistics */
		if (diff_out > 0)  {
			cout << "*Average\t" << stats->Avg1A << tab << stats->Avg2A << " " <<stats->Avg2B << tab << stats->Avg2A-stats->Avg2B;
			cout << tab << ToDeg*stats->AvgAbsDiff_2Aphi << tab << ToDeg*stats->AvgAbsDiff_2Athe << tab << ToDeg*stats->AvgAbsDiff_2Apsi;
			cout << tab << ToDeg*stats->AvgAbsDiff_2Bphi << tab << ToDeg*stats->AvgAbsDiff_2Bthe << tab << ToDeg*stats->AvgAbsDiff_2Bpsi << endl;
			cout << "Std dev\t" << stats->Std1A << tab << stats->Std2A << tab << stats->Std2B;
			cout << tab << ToDeg*stats->StdAbsDiff_2Aphi << tab << ToDeg*stats->StdAbsDiff_2Athe << tab << ToDeg*stats->StdAbsDiff_2Apsi;
			cout << tab << ToDeg*stats->StdAbsDiff_2Bphi << tab << ToDeg*stats->StdAbsDiff_2Bthe << tab << ToDeg*stats->StdAbsDiff_2Bpsi << endl;
		} else  {
			cout << "Average\t" << stats->Avg1A << tab << stats->Avg2A << tab << stats->Avg2B << tab << stats->Avg2A-stats->Avg2B << endl;
			cout << "Std dev\t" << stats->Std1A << tab << stats->Std2A << tab << stats->Std2B << endl;
		}
		cout << "Count\t\t" << stats->CountA << tab << stats->CountB << endl;
		if (stats->Tie > 0)  cout << "NOTE:  " << stats->Tie << " particles in this pair did not favor either hand A or hand B" << endl;
		if (diff_out > 0)  cout << "*Average for pred. orientations is average of |measured angle - predicted angle|" << endl;
	} else  {                  /* for Global statistics */
		hand_print_outheader(which_one, diff_out);     // Header for screen output
		if (diff_out > 0)  {
			cout << "Average  " << stats->Avg1A << tab << stats->Avg2A << tab << stats->Avg2B;
			cout << tab << ToDeg*stats->AvgAbsDiff_2Aphi << tab << ToDeg*stats->AvgAbsDiff_2Athe << tab << ToDeg*stats->AvgAbsDiff_2Apsi;
			cout << tab << ToDeg*stats->AvgAbsDiff_2Bphi << tab << ToDeg*stats->AvgAbsDiff_2Bthe << tab << ToDeg*stats->AvgAbsDiff_2Bpsi << endl;
			cout << "Std dev  " << stats->Std1A << tab << stats->Std2A << tab << stats->Std2B;
			cout << tab << ToDeg*stats->StdAbsDiff_2Aphi << tab << ToDeg*stats->StdAbsDiff_2Athe << tab << ToDeg*stats->StdAbsDiff_2Apsi;
			cout << tab << ToDeg*stats->StdAbsDiff_2Bphi << tab << ToDeg*stats->StdAbsDiff_2Bthe << tab << ToDeg*stats->StdAbsDiff_2Bpsi << endl;
		} else  {
			cout << "Average\t" << stats->Avg1A << tab << stats->Avg2A << tab << stats->Avg2B << endl;
			cout << "Std dev\t" << stats->Std1A << tab << stats->Std2A << tab << stats->Std2B << endl;
		}
		cout << "Count\t\t" << stats->CountA << tab << stats->CountB << endl;
		if (stats->Tie > 0)  cout << "NOTE:  " << stats->Tie << " particles in this analysis did not favor either hand A or hand B" << endl;
		cout << endl << "NOTE:  Hand A is handedness of input map.  Mirror of input map is hand B." << endl;
	}
	
	return 0;
}	

/**
@brief 	Sets consistent selection values in all (including unused) micrographs.
@param 	*project       	Project structure
@param 	*mg_ang        	angles for micrograph selection, views 1 & 2
@param 	*mg_index      	indices for micrograph selection, views 1 & 2
@param 	*mg_select     	micrograph selection criteria, views 1 & 2
@param	sel_consist 	1 or 2, set other selection values to those of this view
@return int        		0.

	If a field-of-view contains more than two micrographs, the user can
	only use two of them to do the handedness determination.  The unused
	micrographs may need their selection values set to the same as those
	in the used micrographs (views 1 and 2).  In addition, there may be
	differences in the initial selection values for views 1 and 2.  This 
	routine loops through the micrographs in a field-of-view and sets all
	selection values to those of view 1 or 2.

**/
int   hand_select_consist(Bproject* project, double* mg_ang, int* mg_index, int* mg_select, int sel_consist)
{
	Bfield*       field=NULL;        // Placeholder for loop over the fields-of-view in project
	int           pcount1, pcount2;  // count of particles in micrographs
	int           i;                 // counter
	int           index;             // index for standard view
	Bmicrograph*  mg=NULL;           // Micrograph to be reset
	Bmicrograph*  mgstd=NULL;        // Standard view (1 or 2) of a micrograph tilt pair
	Bparticle*    particle=NULL;     // Particles from micrograph to be reset
	Bparticle*    particlestd=NULL;  // Particles from "standard" view (1 or 2), selection values copied from here

	if ( sel_consist < 1  ||  sel_consist > 2 )  {
		cerr << "Error: " << sel_consist << " is a wrong value for setting consistent selection values." << endl << "Set to 1 or 2." << endl;
		bexit(-1);
	}

	index = sel_consist - 1;

	if ( verbose )
		cout << endl << "Selection values in all micrographs set to those of micrograph view "
			<< sel_consist << " (standard)." << endl << "standard micrograph  reset micrograph  (field-of-view):" << endl;
	for ( field=project->field; field; field=field->next )  {  // loop through fields-of-view
		mgstd = field_find_micrograph(field, mg_select[index], mg_index[index], mg_ang[index]);   // select sel_consist view
		pcount1 = micrograph_count_particles(mgstd);
		for ( mg=field->mg; mg; mg=mg->next )  {          // if micrograph not sel_consist view, then reset selection values
			if ( mg != mgstd )  {
				if ( verbose )  cout << mgstd->id << "  " << mg->id << "  (" << field->id << ")" << endl;
				pcount2 = micrograph_count_particles(mg);
				if ( pcount2 == pcount1 )   // number of particles in each micrograph (mg and mgstd) is equivalent
					for ( particlestd=mgstd->part, particle=mg->part;  particlestd && particle;  particlestd=particlestd->next, particle=particle->next )
						particle->sel = particlestd->sel;
				else  {                     // number of particles in mg and mgstd differ
					i = 0;
					for ( particlestd=mgstd->part;  particlestd;  particlestd=particlestd->next )
						for ( particle=mg->part; particle; particle=particle->next )
							if ( particle->id == particlestd->id )  {
								particle->sel = particlestd->sel;
								i++;
							}
					if ( verbose )  {
						cout << "Micrograph " << mgstd->id << " has " << pcount1 << " particles.  Micrograph " << mg->id << " has " << pcount2 << " particles." << endl;
						cout << "Only " << i << " particle select values in micrograph " << mg->id << " were changed." << endl;
					}
				}
			}
		}
	}

	return 0;
}


/**
@brief	Appends micrograph number, image number, and a string tag to a file name.  Then it writes out an image file.
@param 	*img    	Image to write out
@param	imgname		String tag for image, appended to file name
@param	outimg		Prefix and suffix for file name
@param 	fieldnum  	Index number for current field-of-view
@param 	num    		Index number for image
@return	int             0.
**/
int		hand_write_img(Bimage* img, Bstring imgname, Bstring outimg, int fieldnum, int num)
{
	outimg = outimg.pre_rev('.') + Bstring(fieldnum, "_%04d") + Bstring(num, "_%04d_") + imgname + "." + outimg.post_rev('.');
	write_img(outimg, img, 0);

	return 0;
}

