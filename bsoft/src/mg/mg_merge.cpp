/**
@file	mg_merge.cpp
@brief	Library functions to merge images
@author David Belnap and Bernard Heymann
@date	Created: 20030410
@date	Modified: 20150805(BH)
**/
 
#include "mg_merge.h"
#include "mg_select.h"
#include "utilities.h"
 
// Definition of the global variables 
extern int    verbose;     // Level of output to the screen


/**
@author	Bernard Heymann and David Belnap
@brief 	Shifts and rotates each image as defined in individual shift vectors.
@param 	*p				image(s) to be rotated and shifted (converted to floating point).
@param 	*origin	array of 3-value origin vectors, one for each image.
@param 	*shift	array of 3-value shift vectors, one for each image.
@param 	angle				global rotation angle to apply to all images.
@return int						number of images.

	Each image in a Bimage structure is shifted by an unique amount but 
	rotated by the same angle, no scaling or resizing is done.
	Intended for use in merging single particle images from a focal, or
	other, series of particles from micrographs.

**/
int			img_unique_shift_global_rotate(Bimage* p, Vector3<float>* origin, Vector3<float>* shift, float angle)
{
	Vector3<double>	start;
	Vector3<double>	scale(1,1,1), axis(0,0,1);
	Vector3<long>	size(p->size());
	Matrix3			mat = Matrix3(axis, angle);
	Bimage*			p1;
	
	p->change_type(Float);
	
	long   n;
	
	for ( n=0; n<p->images(); n++ ) {
//		p1 = p->extract(n);
//		img_transform(p1, size, scale, origin[n], shift[n], mat, FILL_BACKGROUND, 0);
		p1 = p->transform(n, size, scale, origin[n], shift[n], mat, FILL_BACKGROUND, 0);
		p->replace(n, p1);
		delete p1;
	}
	
	return n;
}


/**
@author David Belnap and Bernard Heymann
@brief 	Aligns and merges (sums) 2D particle images from a focal, 
	or other, series and writes to new image file(s).

	Find the reference micrograph based on the input criterion.  Read 
	the reference particle images first, followed by the non-reference
	images.  Find shifts and origins, shift each non-reference image to 
	align the two origins and rotate by the difference in in-plane
	rotational angle (if any) of the two micrographs.  Return the 
	normalized sum of images.  Images are converted to floating point.

@param 	*project       Project parameter structure
@param 	mg_ref_select		0, closest-to-focus; 1, furthest-from-focus; 2, by index; 3, by rotation angle
@param 	mg_index         	Reference by its index in field (for mg_ref_select=2)
@param 	mg_rot_ang        Reference by micrograph rotation angle (for mg_ref_select=3)
@param 	mg_ori_select     	0, from parameter file; 1, from cross-correlation; 2, from images; -1, no alignment
@param 	outimg			prefix and extension for output files
@return int        				error code.
**/
int			mg_particle_merge_series(Bproject* project, int mg_ref_select, int mg_index, 
				float mg_rot_ang, int mg_ori_select, Bstring outimg)
{

	long			nmgs;					// total number of micrographs in field
	long			i, j;					// indices
	long			ifield;					// index for field block
	long			long_ids(0);			// flag to test for long micrograph ids
	double			inplane_angle(0);		// rotation for an entire micrograph
	long			datasize;				// number of data points in images
	Vector3<float>*	shifts = NULL;			// array of translational shifts
	Vector3<float>*	origins = NULL;			// array of origins
	Bstring			filename;				// outimg + micrograph ids
	Bstring			insert("_merged.");		// Insert for default file name

	Bimage*			pimg = NULL;			// Particle image being processed
	Bimage*			pref = NULL;			// Reference particle image(s)
	Bimage*			psum = NULL;			// Summed image(s)
	
	Bfield*			field;					// Field currently being processed
	Bmicrograph*	mg = NULL;				// Micrograph currently being processed
	Bmicrograph*	mg_ref = NULL;			// Reference micrograph, the one to which other images are aligned
	Bmicrograph*	mg2 = NULL;				// Micrograph place holder
	Bparticle*		particle = NULL;		// Particle currently being processed
	Bparticle*		particle_ref = NULL;	// Reference particle


	if ( verbose )  {
		cout << "Merging particle images on micrograph series:" << endl;
		if ( mg_ori_select ==  0 )       cout << "Origins obtained from parameter files." << endl;
		else if ( mg_ori_select ==  1 )  cout << "Origins obtained from cross-correlation with reference micrograph." << endl;
		else if ( mg_ori_select ==  2 )  cout << "Origins obtained from particle image files." << endl;
		else if ( mg_ori_select == -1 )  cout << "Images are merged without translational or rotational alignment." << endl;
	}
  
	// Loop over all fields in project
	for (  field = project->field, ifield=0;  field;  field = field->next, ifield++  )  {

		// Determine reference micrograph for current field
		mg_ref = field_find_micrograph(field, mg_ref_select, mg_index, mg_rot_ang);

		// Read in reference image file as the image to which the other(s) will be added
		//   -copy to pref if user wishes to find origins
		if (verbose & VERB_RESULT)
			cout << endl << "Reference image file for field " << field->id << ":  " << mg_ref->fpart << endl;
		
		psum = read_img(mg_ref->fpart, 1, -1);
		if ( !psum )  {
			error_show("Error in mg_particle_merge_series", __FILE__, __LINE__);
			cerr << "File " << mg_ref->fpart << " was not found or read in properly." << endl;
			return -1;
		}
		psum->change_type(Float);  // convert to floating point
		if ( mg_ori_select == 1 )  pref = psum->copy();   // as reference for finding origins

		datasize = psum->data_size();
		nmgs = field_count_micrographs(field);   // number of micrographs in current field, for normalizing summed images

		// Loop through micrographs in field,
		// if not reference micrograph, align (if desired) each
		// image to corresponding reference image, sum images
		for ( mg = field->mg;  mg;  mg = mg->next )  if ( mg != mg_ref )  {
			if (verbose & VERB_RESULT)
				cout << "Aligning and adding images in " << mg->fpart << " to reference in " << mg_ref->fpart << endl;

			// Calculate rotation for all particles in micrograph, if alignment is desired
			if ( mg_ori_select > -1 )  {
				inplane_angle =  mg->rot_angle - mg_ref->rot_angle;
				if (verbose & VERB_RESULT)
					cout << "in-plane rotation angle = " << (180/M_PI)*inplane_angle << endl;
			}

			// Get images from "non-reference" micrograph, convert to float
			pimg = read_img(mg->fpart, 1, -1);
			pimg->change_type(Float);

			// Check that psum and pimg have the same number of images
			if ( pimg->images() != psum->images() )  {
				cerr << "WARNING:  The number of images is not consistent (" << psum->images() << " and " << pimg->images() << ") and may cause program to crash." << endl;
			}

			// Get origins and translational shift of non-std images, then align to reference images
			if ( mg_ori_select > -1 )  {
				origins = new Vector3<float>[pimg->images()];

				// get origins from input parameter file
				if ( mg_ori_select == 0 )  {
					shifts  = new Vector3<float>[psum->images()];
					for ( i=0, particle=mg->part, particle_ref=mg_ref->part; particle; i++, particle=particle->next, particle_ref=particle_ref->next)  {
						shifts[i] = particle_ref->ori - particle->ori;
						origins[i] = particle->ori; 
					}

				}

				// find new origins via cross-correlation
				else if ( mg_ori_select == 1 )  {
//					img_find_shift(pimg, pref, NULL, pref->sampling(0)[0]*2, 1e10, pref->sizeX()/4.0, 0, 1);
					pimg->find_shift(pref, NULL, pref->sampling(0)[0]*2, 1e10, pref->sizeX()/4.0, 0, 1);
					for ( i=0, particle=mg->part;  particle;  i++, particle=particle->next) {
						shifts[i] = pref->image[i].origin() - pimg->image[i].origin();
						origins[i] = pimg->image[i].origin();
					}
				}

				// get origins from Bsubimage structure
				else if ( mg_ori_select == 2 )  {
					shifts  = new Vector3<float>[psum->images()];
					for ( i=0, particle=mg->part, particle_ref=mg_ref->part;  particle;  i++, particle=particle->next, particle_ref=particle_ref->next)  {
						origins[i] = pimg->image[i].origin();
						shifts[i] = psum->image[i].origin() - origins[i];
					}
				}

				// align non-reference images to reference images
				img_unique_shift_global_rotate(pimg, origins, shifts, inplane_angle);
				delete[] shifts;
				delete[] origins;

			}

			// Add non-reference image to corresponding reference image
			//   -sum background values
			for ( i=0, j=0; i<datasize; i++, j++ ) psum->set(i, (*psum)[i] + (*pimg)[j]);
			for ( i=0; i<psum->images(); i++ )   psum->background(i, psum->background(i) + pimg->background(i));

			delete pimg;                         // free allocated memory for non-reference image

		}   // End of for loop (micrographs in field) and if statement (analyze if not a reference micrograph)
		if ( mg_ori_select == 1 )  delete pref;   // used as reference for finding origins

		// Normalize summed images
		for ( i=0; i<datasize; i++ )  psum->set(i, (*psum)[i] / nmgs);
		for ( i=0; i<psum->images(); i++ )   psum->background(i, psum->background(i) / nmgs);

		// Write out summed, normalized images
		long_ids = 0;
		if ( outimg.length() )  {  // If adding micrograph ids to filename, make sure all ids are <= 6 characters
			for ( mg = field->mg;  mg;  mg = mg->next )  if ( (mg->id.length() > 6) || (mg->id.length() == 0) )  long_ids++;
			if ( long_ids )
				cerr << "Error: Micrograph ids are too long to add together, revert to default file name." << endl;
		}  else  long_ids = 1;

		if ( !long_ids )  {  // Add micrograph ids or "_merged" to an input file name
			filename = outimg.pre_rev('.') + "_" + mg_ref->id + "." + outimg.post_rev('.');
			for ( mg = field->mg;  mg;  mg = mg->next )  if ( mg != mg_ref )  {  // Test if micrograph is reference, if not add id to file name
				filename = filename.pre_rev('.') + "_" + mg->id + "." + filename.post_rev('.');
			}
		} else  filename = mg_ref->fpart.pre_rev('.') + insert + mg_ref->fpart.post_rev('.');

		if ( verbose & VERB_RESULT )
			cout << "Writing summed images to " << filename << endl;  // Write images to file
		
		write_img(filename, psum, 0);
		delete psum;

		// Delete non-reference micrographs from project structure
		for ( mg = field->mg; mg; ) {
			mg2 = mg->next;
			if ( mg != mg_ref ) micrograph_kill(mg);
			mg = mg2;
		}
		field->mg = mg_ref;
		mg_ref->next = NULL;
		mg_ref->block = ifield;
		mg_ref->fpart = filename;

	}   // end of loop over all fields in project

	return 0;
}

/**
@author David Belnap
@brief 	Gives corresponding particles in a field-of-view the same 
	orientation, figure-of-merit (FOM), and selection--allows selection
	of best orientation based on best FOM.
@param 	*project		Project structure for micrographs
@param 	*orientations	Project structure containing orientations to use
@param	fom_diff		Threshold for difference between FOMs
@return int        		error code.

	Intended for use in un-merging data from a micrograph series, where 
	the particle images had been merged previously and there is only one
	set of parameters per field.
	  Set orientations, FOMs, and selection for corresponding particles 
	in a field-of-view to that found in the project structure named 
	orientations.  If user wishes to select the orientation and
	selection value with the highest FOM, the highest FOM (within the
	specified threshold) is selected.
	  Differences in micrograph rotation angles (with respect to the
	"reference" micrograph that contains the orientations to be applied)
	are applied to the output orientation.

**/
int			mg_particle_unmerge(Bproject* project, Bproject* orientations, float fom_diff)
{

	Bfield*         field = NULL;             // Field-of-view currently being processed
	Bfield*         fieldo = NULL;            // Corresponding field-of-view with orientations to be used
	Bmicrograph*    mg = NULL;                // Micrograph currently being processed
	Bmicrograph*    mgo = NULL;               // Corresponding micrograph with orientations to be used
	Bparticle*      particle = NULL;          // Particle currently being processed
	Bparticle*      particleo = NULL;         // Corresponding micrograph with orientation to be used
	int             change;                   // Flag to substitute values from orientations structure
	int             ccount(0);               // Counts number of values that are changed
	int             num, numo;                // For checking that the number of particles is consistent
	int             total;                    // Total number of particles, changed or unchanged
	int             ucount(0);               // Counts number of values that are not changed
	float           rotang_diff;              // Difference in rotation angle between mgo and mg

	project->euler_flag = orientations->euler_flag;  // project gets euler and omega flags of the orientations, since
	project->omega_flag = orientations->omega_flag;  // these general parameters relate to the orientation

	for (  field=project->field, fieldo=orientations->field;  field && fieldo;  field=field->next, fieldo=fieldo->next  )  {
		if ( field->id != fieldo->id )  {
			error_show("Error in mg_particle_unmerge", __FILE__, __LINE__);
			cerr << "Field-of-view IDs do not match.  Please check your data." << endl;
			cerr << "ID of reference field:  " << fieldo->id << endl;
			cerr << "ID of other field:      " << field->id << endl;
			return -1;
		}
		if ( verbose )  cout << "\nUnmerging field:  " << fieldo->id << endl;
		for ( mgo=fieldo->mg;  mgo;  mgo=mgo->next )  {
			if ( verbose & VERB_PROCESS )  cout << "   reference micrograph:  " << mgo->id << endl;
			for ( mg=field->mg;  mg;  mg=mg->next )  {
				if ( verbose & VERB_PROCESS )  cout << "  unmerging " << mg->id << endl;
				numo = micrograph_count_particles(mgo);   num = micrograph_count_particles(mg);
				if ( num != numo )  {
					error_show("Error in mg_particle_unmerge", __FILE__, __LINE__);
					cerr << "Error:  The number of particles in micrographs " << mg->id 
						<< " (" << num << ") and " << mgo->id << " (" << numo 
						<< ") do not match." << endl;
					cerr << "Please check your data" << endl;
					return -1;
				}
				rotang_diff = mg->rot_angle - mgo->rot_angle;
				for ( particle=mg->part, particleo=mgo->part;  particle && particleo;  particle=particle->next, particleo=particleo->next )  {
					if ( particle->id != particleo->id )  {
						error_show("Error in mg_particle_unmerge", __FILE__, __LINE__);
						cerr << "Error:  Particle IDs do not match.  Please check your data" << endl;
						return -1;
					}
					change = 1;
					if ( fom_diff > 0 ) {
						if ( particle->fom[0] - particleo->fom[0] > fom_diff ) {  // Don't change if particle FOM is better by fom_diff
							change = 0;
							ucount++;
						}
					}
					if ( change )  {
						particle->view = particleo->view;
						particle->fom[0]  = particleo->fom[0];
						particle->sel  = particleo->sel;
						particle->view[3] += rotang_diff;   // applies difference in micrograph rotation angle
						ccount++;
					}
				}
			}
		}
	}

	total = ccount + ucount;
	if ( verbose )  {
		cout << "\nUn-merge option:" << endl;
		if ( fom_diff > 0 )  {
			cout << "Total number of particles:      " << total << endl;
			cout << "Number of unchanged particles:  " << ucount << endl;
			cout << "FOM-difference threshold:       " << fom_diff << endl << endl;
		} else
			cout << "Particle orientations changed:  " << ccount << endl << endl;
	}

	return 0;
}


