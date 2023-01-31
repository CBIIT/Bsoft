/**
@file	mg_tomography.cpp
@brief	Functions to manage tomographic series of images
@author	Bernard Heymann
@date	Created: 20020416
@date	Modified: 20200306
**/

#include "mg_processing.h"
#include "mg_tomography.h"
#include "mg_select.h"
#include "mg_img_proc.h"
#include "rwimg.h"
#include "matrix_linear.h"
#include "simplex.h"
#include "linked_list.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Initializes markers from the reference.
@param 	*mg			micrograph to initialize.
@param 	*model		marker model.
@return int			0.
**/
int			mg_marker_init(Bmicrograph* mg, Bmarker* model)
{
	if ( !model )
		return error_show("No markers found!", __FILE__, __LINE__);
	
	Bmarker*		mark = NULL;
	Bmarker*		mark2 = NULL;

	for ( mark = model; mark; mark=mark->next ) {
		mark2 = (Bmarker *) add_item((char **) &mark2, sizeof(Bmarker));
		if ( !mg->mark ) mg->mark = mark2;
		mark2->id = mark->id;
		mark2->fom = 1;
		mark2->sel = mark->sel;
	}
	
	return 0;
}

/**
@brief 	Calculates the predicted positions of markers from the reference.
@param 	*mg				micrograph to update.
@param 	*model			marker model.
@param 	oriref			reference origin.
@param 	update_location	flag to update location based on marker model.
@return int				0.

	If the marker does not exist in the micrograph, a new one is generated.

**/
int			mg_marker_update(Bmicrograph* mg, Bmarker* model, Vector3<double> oriref, int update_location)
{
	if ( !model )
		return error_show("No model markers found!", __FILE__, __LINE__);
	
	int				upd;
	Bmarker*		mark = NULL;
	Bmarker*		modmark = NULL;
	
//	if ( !mg->mark ) mg_marker_init(mg, model);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mg_marker_update: mg=" << mg->id << endl;
		
	for ( modmark = model; modmark; modmark = modmark->next ) {
		upd = update_location;
		for ( mark = mg->mark; mark && mark->id != modmark->id; mark = mark->next ) ;
		if ( !mark ) {
			mark = (Bmarker *) add_item((char **) &mg->mark, sizeof(Bmarker));
			mark->id = modmark->id;
			mark->fom = 1;
			mark->sel = modmark->sel;
			upd = 1;
		}
		if ( upd )
			mark->loc = mg_location_from_3D_model(modmark->loc, oriref, mg->matrix, mg->origin);
	}
	
	return 0;
}

Vector3<double>	mg_3D_location_from_marker(Vector3<double> loc, Vector3<double> ori3D, Matrix3 mat, Vector3<double> ori)
{
	Vector3<double>	loc3D = loc - ori;
	loc3D = mat * loc3D;
	loc3D += ori3D;
//	loc3D[2] = 0;
	
	return loc3D;
}

/**
@brief 	Fills in missing markers in the 3D model from micrograph markers.
@param 	*mg			micrograph to update.
@param 	*model		marker model.
@param 	oriref		reference origin.
@return int			0.

	If the marker does not exist in the model, a new one is generated.

**/
int			mg_update_model(Bmicrograph* mg, Bmarker* model, Vector3<double> oriref)
{
	if ( !model )
		return error_show("No model markers found!", __FILE__, __LINE__);
	
	Bmarker*		mark = NULL;
	Bmarker*		modmark = NULL;
	
	for ( mark = mg->mark; mark; mark = mark->next ) {
		for ( modmark = model; modmark && mark->id != modmark->id; modmark = modmark->next ) ;
		if ( !modmark ) {
			modmark = (Bmarker *) add_item((char **) &model, sizeof(Bmarker));
			modmark->id = mark->id;
			modmark->fom = 1;
			modmark->sel = mark->sel;
			modmark->loc = mg_3D_location_from_marker(mark->loc, oriref, mg->matrix, mg->origin);
		}
	}
	
	return 0;
}

/**
@brief 	Calculates the 2D location of a marker from a 3D model.
@param 	loc3D			3D marker location.
@param 	mat				micrograph view matrix.
@param 	ori				micrograph origin.
@return Vector3<double>	transformed 2D marker location.
**/
Vector3<double>	mg_location_from_3D_model(Vector3<double> loc3D, Matrix3 mat,
				Vector3<double> ori)
{
	Matrix3			tmat = mat.transpose();
	Vector3<double>	mloc = tmat * loc3D;
	mloc += ori;
	mloc[2] = 0;
	
	return mloc;
}

Vector3<double>	mg_location_from_3D_model(Vector3<double> loc3D, Matrix3 mat,
				Vector3<double> ori, Vector3<double> scale)
{
	Matrix3			tmat = mat.transpose();
	Vector3<double>	mloc = tmat * loc3D;
	mloc /= scale;
	mloc += ori;
	mloc[2] = 0;
	
	return mloc;
}

Vector3<double>	mg_location_from_3D_model(Vector3<double> loc3D, Vector3<double> ori3D,
				Matrix3 mat, Vector3<double> ori)
{
	Vector3<double>	mloc = loc3D - ori3D;
	return mg_location_from_3D_model(mloc, mat, ori);
}

Vector3<double>	mg_location_from_3D_model(Vector3<double> loc3D, Vector3<double> ori3D,
				Matrix3 mat, Vector3<double> ori, Vector3<double> scale)
{
	Vector3<double>	mloc = loc3D - ori3D;
//	cout << loc3D << endl << ori3D << endl << ori << endl << scale << endl;
	return mg_location_from_3D_model(mloc, mat, ori, scale);
}

/**
@brief 	Generates an image representing a gold particle.
@param 	size	size of image to generate.
@param 	radius		gold particle radius in voxel units.
@return Bimage* 			gold particle image.

	An image of the size of an input sub-image is generated with a gold
	particle (black or negative density) at its center.

**/
Bimage*		img_gold_particle(Vector3<long> size, double radius)
{
	Bimage*			p = new Bimage(Float, TSimple, size, 1);
	Vector3<double> origin(p->size()/2);
	double			outrad = sqrt(2.0)*radius;
	
	p->sphere(origin, outrad, 2, FILL_USER, 1);
	p->sphere(origin, radius, 2, FILL_USER, -1);
	
	p->origin(origin);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG img_gold_particle: rad=" << radius << endl;
		write_img("gold.map", p, 0);
	}
	
	p->statistics();
	
	return p;
}

/**
@brief 	Finds gold particles in micrograph images.
@param 	*p			image.
@param	img_select	selected subimage, -1 if all.
@param 	radius		gold particle radius in voxel units.
@param 	edge		edge width to eliminate particles.
@param 	cutoff		FOM cutoff, if 0 then half of maximum FOM.
@return Bmarker* 	list of markers.

	An image of the size of an input sub-image is generated with a gold
	particle (black or negative density) at its center.
	Each sub-image is cross-correlated with the gold particle image
	and the position within the image reported.
	The FOM cutoff is used to select cross-correlation peaks, except it
	is reset for the following:
		cutoff == 0 ==> cutoff = FOMavg + FOMstd
		cutoff <  0 ==> cutoff = FOMmax/2

**/
Bmarker*	img_find_gold_particles(Bimage* p, int img_select, double radius, long edge, double cutoff)
{
	if ( radius < 1 ) return NULL;
	
	int			img_start(0);
	int			img_end(p->images());
	if ( img_select > -1 ) {
		if ( img_select >= p->images() ) img_select = p->images() - 1;
		img_start = img_select;
		img_end = img_start + 1;
	}

	Bimage*			pgold = img_gold_particle(p->size(), radius);
	
    long			i, j, x, y, z, n, nn, np(0);
	long			xx, yy, zz;
	double			fom_avg(0), fom_std(0), fom_max(0);
	Vector3<long>	ve(p->size()/2);
	ve = ve.min(edge);
	
	if ( verbose & VERB_LABEL ) {
		cout << "Finding particles:" << endl;
		cout << "Particle radius:                " << radius*p->sampling(0)[0] 
			<< " A " << radius << " pixels" << endl;
		cout << "Edge:                           " << edge << " pixels" << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_find_gold_particles: radius=" << radius << endl;
	
	Bimage*			pone = NULL;
	Bimage*			pcc = NULL;
	
	long			imgsize(p->sizeX()*p->sizeY()*p->sizeZ());
	long			fomsize(imgsize);
	if ( img_select < 0 ) fomsize *= p->images();
	float*			fom = new float[fomsize];
	for ( i=0; i<fomsize; i++ ) fom[i] = 0;
	
//	if ( verbose & VERB_PROCESS )
//		cout << "Image\tx\ty\tz" << endl;
	for ( i=0, n=img_start; n<img_end; n++ ) {	// i is the external counter
		pone = p->extract(n);
		pcc = pone->cross_correlate(pgold, 2*p->sampling(0)[0], 8*radius*p->sampling(0)[0]);
		delete pone;
		pone = pcc->find_peaks(3);
		delete pcc;
		for ( j=0; j<imgsize; j++, i++ ) if ( (*pone)[j] ) {
			fom[i] = (*pone)[j];
			fom_avg += fom[i];
			fom_std += fom[i]*fom[i];
			if ( fom_max < fom[i] ) fom_max = fom[i];
			np++;
		}
		delete pone;
		if ( verbose & VERB_PROCESS )
			cout << "Image " << n << ": " << np << " gold fiducials" << endl;
	}

	fom_avg /= np;
	fom_std = fom_std/np - fom_avg*fom_avg;
	if ( fom_std > 0 ) fom_std = sqrt(fom_std);
	else fom_std = 0;
	if ( cutoff == 0 ) cutoff = fom_avg + fom_std;
	else if ( cutoff < 0 ) cutoff = fom_max/2;
	if ( verbose )
		cout << "Cutoff = " << cutoff << endl;
		
	Bmarker*		mark = NULL;
	Bmarker*		m = NULL;

	if ( np ) {
		if ( verbose )
			cout << "Marker\tx\ty\tz\tn\tFOM" << endl;
		for ( j=nn=0, n=img_start; n<img_end; n++, nn++ ) {
			for ( z=0; z<p->sizeZ(); z++ ) {
				zz = (long)z - pgold->image->origin()[2];
				if ( zz < 0 ) zz += p->sizeZ();
				if ( zz >= ve[2] && zz < p->sizeZ()-ve[2] )
						for ( y=0; y<p->sizeY(); y++ ) {
					yy = (long)y - pgold->image->origin()[1];
					if ( yy < 0 ) yy += p->sizeY();
					if ( yy >= ve[1] && yy < p->sizeY()-ve[1]  )
							for ( x=0; x<p->sizeX(); x++ ) {
						xx = (long)x - pgold->image->origin()[0];
						if ( xx <  0 ) xx += p->sizeX();
						if ( xx >= ve[0] && xx < p->sizeX()-ve[0] ) {
							i = p->index(x, y, z, nn);
							if ( fom[i] > cutoff ) {
								m = (Bmarker *) add_item((char **) &m, sizeof(Bmarker));
								if ( !mark ) mark = m;
								j++;
								m->id = j;
								m->loc[0] = x - pgold->image->origin()[0];
								m->loc[1] = y - pgold->image->origin()[1];
								m->loc[2] = z - pgold->image->origin()[2];
								m->img_num = n;
								m->sel = 1;
								m->fom = fom[i];
								if ( m->loc[0] < 0 ) m->loc[0] += p->sizeX();
								if ( m->loc[1] < 0 ) m->loc[1] += p->sizeY();
								if ( m->loc[2] < 0 ) m->loc[2] += p->sizeZ();
								if ( verbose )
									cout << m->id << tab << m->loc[0] << tab << 
										m->loc[1] << tab << m->loc[2] << tab << 
										m->img_num << tab << m->fom << endl;
							}
						}
					}
				}
			}
		}
	}
	
	delete pgold;

	delete[] fom;
	
	if ( verbose & VERB_LABEL )
		cout << "Number of particles found:      " << j << endl << endl;
	
	return mark;
}

/**
@brief 	Normalizes a set of images in a project to a desired average and standard deviation.
@param 	*project			project structure.
@param 	avg					desired average.
@param 	std					desired standard deviation (if 0, use defaults).
@param 	norm_type			type of determining the effective average and standard deviation:
								0=simple, 1=Gaussian, 2=Poisson.
@param 	datatype			data type for normalized images.
@param 	setinputZslices		convert z-slices in the input to 2D images.
@param 	setoutputZslices	convert output image back to z-slices.
@param 	cutmin				truncate to minimum.
@param 	cutmax				truncate to maximum.
@param 	replace_threshold	replace maxima above this threshold.
@return int					0, <0 on error.

	The effective average and standard deviation for each image is obtained
	in one of three ways:
		0.		The simple avergae and standard devaition for the image.
		1.		Gaussian fit of the histogram.
		2.		Poisson fit of the histogram.
	A histogram of an image is calculated with a given number of bins.
	The histogram is fit to a Gaussian or Poisson function with exclusion of a
	small number of bins in the histogram (defined as outliers).
	The effective average and standard deviation are used to 
	rescale the data for each image.

**/
int			project_mass_normalize(Bproject* project, double avg, double std, int norm_type,
				DataType datatype, int setinputZslices, int setoutputZslices,
				double cutmin, double cutmax, double replace_threshold)
{
	if ( !project->field ) return -1;
	
	int				nimg(0);
	Bfield* 		field;
	Bmicrograph*	mg;
	Bimage*			p;
	Bstring			oldname;
	Bstring			normname;
	DataType		nudatatype;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( oldname != mg->fmg ) {
				nimg++;
				p = read_img(mg->fmg, 1, -1);
				if ( !p ) {
					error_show("project_mass_normalize", __FILE__, __LINE__);
					return -1;
				}
				nudatatype = datatype;
				if ( nudatatype == Unknown_Type )
					nudatatype = p->data_type();
				else if ( nudatatype > p->data_type() )
					p->change_type(nudatatype);
	
				if ( setinputZslices ) 
					if ( p->slices_to_images() < 0 ) {
						error_show("project_mass_normalize", __FILE__, __LINE__);
						return -1;
					}

				if ( cutmin < cutmax ) p->truncate_to_min_max(cutmin, cutmax);
	
				if ( replace_threshold != 0 )
					p->replace_maxima(replace_threshold);
				
				if ( norm_type >= 0 )
					p->normalize(avg, std, norm_type);
				
				if ( setoutputZslices ) 
					if ( p->images_to_slices() < 0 ) {
						error_show("project_mass_normalize", __FILE__, __LINE__);
						return -1;
					}
				
				if ( p->data_type() != nudatatype ) {
					if ( nudatatype == UCharacter ) p->truncate_to_min_max(0, 255);
					p->change_type(nudatatype);
				}
				
				normname = mg->fmg.pre_rev('.') + "_norm." + mg->fmg.post_rev('.');
				write_img(normname, p, 0);
				delete p;
			}
			oldname = mg->fmg;
			mg->fmg = normname;
		}
	}
	
	return 0;
}

double		lambda_R(Bsimplex& simp)
{
	long			i;
	double			R(0), df, fs(0), fs2(0);
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); i++ ) {
		df = cos(x[i] + simp.parameter(2));
		if ( df )
			df = f[i] - simp.parameter(0)*exp(-simp.parameter(1)/df);
		else
			df = f[i] - simp.parameter(0);
		fs += f[i];
		fs2 += f[i] * f[i];
		R += df*df;
	}
	
	df = i*fs2-fs*fs;
	if ( df > 0 ) R = sqrt(R/df);
	else R = sqrt(R/i);
			
	return R;
}

vector<double>	simplex_fit_tilt_intensities(vector<double>& a, vector<double>& I)
{
	long			i;
	double			Imax(0);
	
	for ( i=0; i<I.size(); ++i ) if ( Imax < I[i] ) Imax = I[i];
	
	Bsimplex		simp(1, 3, 0, a.size(), a, I);
	
	simp.parameter(0, 2*Imax);
	simp.parameter(1, 0.6);
	simp.parameter(2, -0.01);
	simp.limits(0, 0.9*Imax, 5*Imax);
	simp.limits(1, 0.01, 10);
	simp.limits(2, -0.5, 0.5);
	
//	simp.show();
	
	long			maxiter(10);
	double			R = simp.run(10000, 0.0001, lambda_R);
//	simp.show();
	for ( i=0; i<maxiter && ( !isfinite(R) || R > 1 ); ++i ) {
		simp.parameter(0, 1.5*Imax);
		simp.parameter(1, 0.7);
		simp.parameter(2, 0.01);
		R = simp.run(10000, 0.0001, lambda_R);
//		simp.show();
	}
	
	vector<double>	result;
	for ( i=0; i<3; ++i )
		result.push_back(simp.parameter(i));
	result.push_back(R);
	
	return result;
}

/**
@brief 	Transfers the micrograph intensities for a field into a plot.
@param 	*project			project structure.
@return Bplot*				plot structure.

	The micrographs should be unmodified gain-corrected images obtained from the detector.

**/
Bplot*		project_intensity_plot(Bproject* project)
{
	if ( !project->field ) return NULL;

	long			i, j, k;
	Bfield* 		field = project->field;
	Bmicrograph*	mg = NULL;
	Bimage*			p = NULL;
	
	long			ncol(3);
	long			nmg = field_count_mg_selected(field);
	double			amin(0), amax(0), Imin(1e30), Imax(-1e30);
	Bstring			title("Intensity-tilt curve: " + project->filename);
//	cout << "Plot title: " << title << endl;
	Bplot*			plot = new Bplot(1, nmg, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; ++i ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Tilt angle");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(0);
	plot->page(0).column(1).element_size(5);
	plot->page(0).column(1).label("Intensity");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(2).type(2);
	plot->page(0).column(2).label("Int(calc)");
	plot->page(0).column(2).axis(3);
	plot->page(0).column(1).color(1,0,0);
	plot->page(0).column(2).color(0,0,1);
	
//	cout << "nmg=" << nmg << endl;
	
//	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg->intensity < 1e-10 ) {
				if ( mg->fmg.length() ) {
					p = read_img(mg->fmg, 0, mg->img_num);
//					cout << mg->img_num << tab << p->image->average() << endl;
					if ( fabs(p->image->average()) < 1e-20 ) {
						delete p;
						p = read_img(mg->fmg, 1, mg->img_num);
						p->statistics();
//						cout << mg->img_num << tab << p->image->average() << endl;
					}
				} else if ( mg->fframe.length() ) {
					p = read_img(mg->fframe, 0, mg->img_num);
					if ( fabs(p->image->average()) < 1e-20 ) {
						delete p;
						p = read_img(mg->fframe, 1, mg->img_num);
						p->statistics();
					}
				}
				if ( !p ) {
					error_show("project_thickness_to_lambda_ratio", __FILE__, __LINE__);
					return NULL;
				} else {
					mg->intensity = p->image->average();
					delete p;
				}
			}
			if ( mg->intensity < 0.1 ) mg->select = 0;
		}

		for ( i=0, j=plot->rows(), k=2*j, mg = field->mg; mg; mg = mg->next )
				if ( mg->select ) {
			(*plot)[i] = mg->tilt_angle;
			(*plot)[j] = mg->intensity;
			(*plot)[k] = 0;
			if ( amin > mg->tilt_angle ) amin = mg->tilt_angle;
			if ( amax < mg->tilt_angle ) amax = mg->tilt_angle;
			if ( Imin > mg->intensity ) Imin = mg->intensity;
			if ( Imax < mg->intensity ) Imax = mg->intensity;
			i++; j++; k++;
		}
//	}

	double		f(Imax-Imin);
	if ( -amin < amax ) amin = -amax;
	else amax = -amin;
	Imin = f*long(Imin/f);
	Imax = f*long(Imax/f+1);
	plot->page(0).column(0).min(amin*180.0/M_PI);
	plot->page(0).column(0).max(amax*180.0/M_PI);
	plot->page(0).column(1).min(Imin);
	plot->page(0).column(1).max(Imax);
	plot->page(0).column(2).min(Imin);
	plot->page(0).column(2).max(Imax);

//	cout << "setup done " << endl;
	
	return plot;
}

/**
@brief 	Determines the thickness to mean free path ratio from a tilt series.
@param 	*project			project structure.
@param	*plot				plot structure with intensity-tilt data.
@param	flag				flag to apply tilt angle adjustment.
@return double				ratio.

	The micrographs should be unmodified gain-corrected images obtained from the detector.

**/
double		project_fit_intensities(Bproject* project, Bplot* plot, int flag)
{
	long			i, j, k;
	Bfield* 		field = project->field;
	Bmicrograph*	mg = NULL;
	long			nmg = field_count_mg_selected(field);

//	cout << "nmg = " << nmg << endl;

	if ( nmg < 3 ) {
		if ( verbose )
			cerr << "Too few micrographs to determine the thickness!" << endl;
		return 0;
	}
	
	if ( !plot )
		plot = project_intensity_plot(project);

	vector<double>	a(plot->rows()), I(plot->rows());
	for ( i=0, j=plot->rows(); i<plot->rows(); ++i, ++j ) {
		a[i] = (*plot)[i];
		I[i] = (*plot)[j];
	}
	
//	cout << "starting fit" << endl;
	vector<double>	result = simplex_fit_tilt_intensities(a, I);
	
	double			I0 = result[0];
	double			tlr = result[1];
	double			da = result[2];
	double			R = result[3];
	
	double			Ip;
	if ( verbose )
		cout << endl << "Tilt\tInt\tInt(calc)" << endl;
	for ( i=0, j=plot->rows(), k=2*j; i<plot->rows(); ++i, ++j, ++k ) {
		Ip = I0*exp(-tlr/cos(a[i]+da));
		(*plot)[i] = a[i]*180.0/M_PI;
		(*plot)[k] = Ip;
		if ( verbose )
			cout << a[i]*180.0/M_PI << tab << I[i] << tab << Ip << endl;
	}

	Bstring			txt;
	txt = Bstring(I0, "Direct beam intensity: %lg");
	plot->page(0).add_text(txt);
	txt = Bstring(tlr, "Tilt-lambda ratio: %lg");
	plot->page(0).add_text(txt);
	txt = Bstring(da*180.0/M_PI, "Specimen tilt angle: %lg");
	plot->page(0).add_text(txt);
	txt = Bstring(R, "Rfit: %lg");
	plot->page(0).add_text(txt);
	
	if ( flag )
		for ( mg = field->mg; mg; mg = mg->next )
			mg->tilt_angle += da;
	
	if ( verbose ) {
		cout << endl << "Calculating tilt intensity parameters:" << endl;
		cout << "Residual:                       " << R << endl;
		cout << "Direct beam intensity:          " << I0 << endl;
		cout << "Specimen tilt:                  " << da*180.0/M_PI << " degrees" << endl;
		cout << "Thickness to lambda ratio:      " << tlr << endl << endl;
	}
	
	return tlr;
}


/**
@brief 	Determines the thickness to mean free path ratio from a tilt series.
@param 	*project			project structure.
@return double				ratio.

	The micrographs should be unmodified gain-corrected images obtained from the detector.

**/
double		project_thickness_to_lambda_ratio(Bproject* project)
{
	if ( !project->field ) return -1;

	long			i, j, k;
	Bfield* 		field = project->field;
	Bmicrograph*	mg = NULL;
	Bimage*			p = NULL;
	
	long			nmg = field_count_mg_selected(field);
	vector<double>	a(nmg), I(nmg);
	
//	cout << "nmg=" << nmg << endl;
	
//	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg->intensity < 1e-10 ) {
				if ( mg->fmg.length() ) {
					p = read_img(mg->fmg, 0, mg->img_num);
//					cout << mg->img_num << tab << p->image->average() << endl;
					if ( fabs(p->image->average()) < 1e-20 ) {
						delete p;
						p = read_img(mg->fmg, 1, mg->img_num);
						p->statistics();
//						cout << mg->img_num << tab << p->image->average() << endl;
					}
				} else if ( mg->fframe.length() ) {
					p = read_img(mg->fframe, 0, mg->img_num);
					if ( fabs(p->image->average()) < 1e-20 ) {
						delete p;
						p = read_img(mg->fframe, 1, mg->img_num);
						p->statistics();
					}
				}
				if ( !p ) {
					error_show("project_thickness_to_lambda_ratio", __FILE__, __LINE__);
					return -1;
				} else {
					mg->intensity = p->image->average();
					delete p;
				}
			}
			if ( mg->intensity < 0.1 ) mg->select = 0;
		}

		for ( i=0, mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			a[i] = mg->tilt_angle;
			I[i] = mg->intensity;
			i++;
		}
//	}

	nmg = i;
	
	if ( nmg < 3 ) {
		if ( verbose )
			cerr << "Too few micrographs to determine the thickness!" << endl;
		return 0;
	}

	vector<double>	result = simplex_fit_tilt_intensities(a, I);
	
	double			I0 = result[0];
	double			tlr = result[1];
	double			da = result[2];
	double			R = result[3];
	
	double			Ip;
	if ( verbose )
		cout << endl << "Tilt\tInt\tInt(calc)" << endl;
	for ( i=0, j=a.size(), k=2*j; i<a.size(); ++i, ++j, ++k ) {
		Ip = I0*exp(-tlr/cos(a[i]+da));
		if ( verbose )
			cout << a[i]*180.0/M_PI << tab << I[i] << tab << Ip << endl;
	}
	
	for ( mg = field->mg; mg; mg = mg->next )
		mg->tilt_angle += da;
	
	if ( verbose ) {
		cout << endl << "Calculating tilt intensity parameters:" << endl;
		cout << "Residual:                       " << R << endl;
		cout << "Direct beam intensity:          " << I0 << endl;
		cout << "Specimen tilt:                  " << da*180.0/M_PI << " degrees" << endl;
		cout << "Thickness to lambda ratio:      " << tlr << endl << endl;
	}
	
	return tlr;
}

/**
@brief 	Determines the mean free path from a tilt series and thickness.
@param 	*project			project structure.
@param	thickness			tomogram thickness (in angstrom).
@return double				ratio.

	The micrographs should be unmodified gain-corrected images obtained from the detector.

**/
double		project_lambda(Bproject* project, double thickness)
{
	if ( !project->field ) return -1;

	double			tlr = project_thickness_to_lambda_ratio(project);

	double			lambda = thickness/tlr;
	
	if ( verbose ) {
		cout << "Thickness:                      " << thickness << " A" << endl;
		cout << "Proportionality parameter:      " << lambda << " A" << endl << endl;
	}
	
	return lambda;
}

/**
@brief 	Determines the thickness from a tilt series and mean free path.
@param 	*project			project structure.
@param	lambda				proportionality parameter (in angstrom).
@return double				ratio.

	The micrographs should be unmodified gain-corrected images obtained from the detector.

**/
double		project_thickness(Bproject* project, double lambda)
{
	if ( !project->field ) return -1;

	double			tlr = project_thickness_to_lambda_ratio(project);

	double			thickness = lambda*tlr;
	
	if ( verbose ) {
		cout << "Thickness:                      " << thickness << " A" << endl << endl;
		cout << "Proportionality parameter:      " << lambda << " A" << endl << endl;
	}
	
	return lambda;
}

/**
@brief 	Sorts markers for each micrograph by ID number.
@param 	*project	micrograph project.
@return int			0.

	Requires the matrices in the micrograph structures to be defined.

**/
int			project_sort_markers_by_id(Bproject* project)
{
	Bfield*			field;
	Bmicrograph*	mg;
	
	if ( verbose )
		cout << "Sorting markers" << endl << endl;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			markers_sort_by_id(&mg->mark);
		
	return 0;
}

/**
@brief 	Calculates the residuals for an alignment solution.
@param 	*project	micrograph project.
@param 	show		flag to indicate showing lists.
@return double		average residual.

	Requires the matrices in the micrograph structures to be defined.

**/
double		project_tomo_residuals(Bproject* project, int show)
{
	if ( project_count_markers(project) < 1 ) return 0;
	
	Bmarker*			m, *mr;
	Bfield*				field = NULL;
	Bmicrograph*		mg;
	Bmicrograph*		mg_ref = NULL;
	Breconstruction*	rec = project->rec;
	Vector3<double>		loc, shift, tilt_axis;
	Quaternion			q;
	Transform			t;
	Matrix3				mat, mat2;
	double				d(0), dsum(0), dmgsum, w(0), wmg, mgfom, avg_axis(0), wa, wax(0);
	
	for ( field = project->field; field && !mg_ref; field = field->next )
		if ( field->select ) mg_ref = field_find_zero_tilt_mg(field);
	
//	long				nmg = count_list((char *) field->mg);

	if ( mg_ref->origin[0] < 1 || mg_ref->origin[1] < 1 )
		project_set_nominal_mg_origins(project);

	if ( mg_ref->matrix.determinant() < 0.5 )
		project_mg_tilt_to_matrix(project);

	Bstring			id("1");
	if ( !project->rec ) {
		rec = reconstruction_add(&project->rec, id);
		rec->origin = mg_ref->origin;
	}
	
	if ( !rec->mark )
		project_calculate_model(project);
		
	for ( mr = rec->mark; mr; mr = mr->next ) {
		mr->err = 0;
		mr->res = 0;
		mr->sel = 0;
	}

	if ( verbose && show ) {
		cout << "Tomography alignment results:" << endl;
		cout << "Mg\tShiftX\tShiftY\tScaleX\tScaleY\tAxis\tTilt\tLevel\tn\tResid\tFOM" << endl;
	}
	
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
/*		if ( field != project->field ) {
			mg_ref2 = field_find_zero_tilt_mg(field);
			if ( mg_ref2->mark ) {
				t = marker_find_transform(mg_ref2->mark, mg_ref->mark, rec->origin);
//				cout << "axis=" << t.axis*180/M_PI << tab << "angle=" << t.angle*180/M_PI << endl;
				field->matrix = Matrix3(t.axis, t.angle);
			}
		} else field->matrix = Matrix3(1);*/
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			shift = mg->origin - mg_ref->origin;
			dmgsum = wmg = mgfom = 0;
//			mat = mg->matrix.transpose();	// Handedness
//			mat = mg->matrix * field->matrix;
			for ( mr = rec->mark; mr; mr = mr->next ) {
				for ( m = mg->mark; m && m->id != mr->id; m = m->next ) ;
				if ( m ) {
//					loc = mg_location_from_3D_model(mr->loc, rec->origin, mat, mg->origin, mg->scale*rec->scale);
					loc = mg_location_from_3D_model(mr->loc, rec->origin, mg->matrix, mg->origin, mg->scale*rec->scale);
					m->err = loc - m->loc;
					d = m->err.length2();
					if ( m->sel > 0 ) {
						dmgsum += d;
						mgfom += m->fom;
						mr->err += m->err;
						mr->res += d;
						mr->sel += 1;
						wmg += 1;
					}
					d = sqrt(d);
					m->res = d;
					if ( ( verbose & VERB_FULL ) && show )
						cout << mg->id << tab << m->id << tab << loc[0] << tab << loc[1] << tab << loc[2] << tab << d << endl;
				}
			}
			dsum += dmgsum;
			w += wmg;
			if ( wmg ) {
				mg->fom = sqrt(dmgsum/wmg);
				mgfom /= wmg;
			}
			wa = 1 - exp(-mg->tilt_angle*mg->tilt_angle/0.5);
			if ( fabs(mg->tilt_angle) > M_PI/45.0 ) {
				avg_axis += wa*mg->tilt_axis;
				wax += wa;
			}
			if ( verbose && show )
				cout << mg->id << tab << fixed << setprecision(2) << shift[0] << tab << shift[1] << tab
					<< setprecision(4) << mg->scale[0] << tab << mg->scale[1] << tab
					<< mg->tilt_axis*180.0/M_PI << tab
					<< mg->tilt_angle*180.0/M_PI << tab
					<< mg->level_angle*180.0/M_PI << tab
					<< wmg << tab << mg->fom << tab << mgfom << endl;
		}
	}
	
	if ( wax ) avg_axis /= wax;
	if ( w ) d = sqrt(dsum/w);

	// Fix tilt axis angles - particularly for the zero-tilt micrograph
//	for ( field = project->field; field; field = field->next ) if ( field->select )
//		for ( mg = field->mg; mg; mg = mg->next )
//			if ( fabs(mg->tilt_axis - avg_axis) > M_PI/90.0 ) mg->tilt_axis = avg_axis;

	for ( mr = rec->mark; mr; mr = mr->next ) {
		if ( mr->sel ) mr->err /= mr->sel;
		mr->res = sqrt(mr->res/mr->sel);
	}

	if ( verbose && show ) {
		cout << "Average tilt axis angle:        " << avg_axis*180.0/M_PI << endl << endl;
		cout << "Marker errors:\nMarker\tn\tx\ty\tResid\tFOM" << endl;
		for ( mr = rec->mark; mr; mr = mr->next )
			cout << mr->id << tab << mr->sel << tab << 
				setprecision(3) << mr->err[0] << tab << mr->err[1] << tab << 
				mr->res << tab << setprecision(4) << mr->fom << endl;
		cout << "Average residual = " << d << " pixels" << endl << endl;
	}

	return d;
}

/**
@brief 	Calculates the residuals for an alignment solution of one micrograph.
@param 	*mg			micrograph.
@param 	*model		model marker list.
@param 	oriref		reference origin.
@return double		average residual.

	Requires the matrices in the micrograph structures to be defined.

**/
double		mg_tomo_residuals(Bmicrograph* mg, Bmarker* model, Vector3<double> oriref)
{
	Bmarker*			m, *mr;
	Vector3<double>		loc;
	double				d(0), dsum(0), w(0);
	
	for ( mr = model; mr; mr = mr->next ) {
		for ( m = mg->mark; m && m->id != mr->id; m = m->next ) ;
		if ( m ) {
			loc = mg_location_from_3D_model(mr->loc, oriref, mg->matrix, mg->origin, mg->scale);
			m->err = loc - m->loc;
			d = m->err.length2();
			if ( m->sel > 0 ) {
				dsum += d;
				w += 1;
			}
			d = sqrt(d);
			m->res = d;
		}
	}
	
	if ( w ) d = sqrt(dsum/w);

	return d;
}

/**
@brief 	Analyzes marker errors.
@param 	*project	micrograph project.
@return double		0.

	Requires the marker errors to be calculated.

**/
double		project_tomo_errors(Bproject* project)
{
	Bfield*				field = project->field;
	Bmicrograph*		mg;
	Bmarker				*m1, *m2;
	
	int					i, n(0);
	double				scale = 10;
	int					h[1000];
	for ( i=0; i<1000; i++ ) h[i] = 0;

	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg->next; mg = mg->next ) {
			for ( m1=mg->mark, m2=mg->next->mark; m1 && m2; m1=m1->next, m2=m2->next ) {
				i = (int) (scale*(m1->err - m2->err).length() + 0.5);
				if ( i < 1000 ) {
					h[i]++;
					if ( n < i ) n = i;
				}
			}
		}
	}
	
	if ( verbose ) {
		cout << "Histogram of marker refinement shifts:" << endl;
		cout << "Shift\tCount" << endl;
		for ( i=0; i<=n; i++ )
			cout << setprecision(2) << i/scale << tab << h[i] << endl;
		cout << endl;
	}
	
	return 0;
}

/**
@brief 	Calculates the z-coordinates of the marker model.
@param 	*project	project parameter structure.
@return int			0.

	If the project does not have markers, it returns without any change.
	If the reconstruction model is not defined, it is created.

**/
int			project_calculate_model(Bproject* project)
{
	if ( project_count_markers(project) < 1 ) return 0;
	
	int					nm;
	Bfield*				field;
	Bmicrograph*		mg;
	Bmarker				*m, *mr;
	Vector3<double>		t;
	Bmicrograph*		mg_ref = NULL;
	Breconstruction*	rec = project->rec;

	for ( field = project->field; field && !mg_ref; field = field->next )
		if ( field->select ) mg_ref = field_find_zero_tilt_mg(field);

	if ( !mg_ref->mark )
		for ( mg_ref = project->field->mg; mg_ref && !mg_ref->mark; mg_ref = mg_ref->next ) ;
	
	if ( !mg_ref )
		return error_show("project_calculate_model: No markers defined!", __FILE__, __LINE__);
 
	Bstring				id("1");
	if ( !project->rec ) {
		rec = reconstruction_add(&project->rec, id);
		rec->origin = mg_ref->origin;
		rec->mark_radius = mg_ref->mark_radius;
	}
	
	if ( verbose & VERB_FULL )
		cout << "Calculating a 3D marker model" << endl << endl;

	// Set 3D model if not already present
	if ( !rec->mark )
		rec->mark = markers_copy(mg_ref->mark);

	// Fill in extra markers
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( m = mg->mark; m; m = m->next ) {
				for ( mr = rec->mark; mr && mr->id != m->id; mr = mr->next ) ;
				if ( !mr ) {
					mr = (Bmarker *) add_item((char **) &rec->mark, sizeof(Bmarker));
					if ( !rec->mark ) rec->mark = mr;
					mr->id = m->id;
					mr->loc = m->loc;
					mr->res = m->res;
					mr->fom = m->fom;
					mr->sel = m->sel;
				}
			}
		}
	}
	
	// Calculate the micrograph origins
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			t = 0;
			nm = 0;
			for ( mr = rec->mark; mr; mr = mr->next ) {
				for ( m = mg->mark; m && m->id != mr->id; m = m->next ) ;
				if ( m ) {
					t += m->loc - mg_location_from_3D_model(mr->loc, rec->origin, mg->matrix, mg->origin);
					nm++;
				}
			}
			if ( nm ) {
				t /= nm;
				mg->origin += t;
			}
		}
	}
	
	// Calculate the z-coordinates
	for ( mr = rec->mark; mr; mr = mr->next ) {
		mr->loc[2] = 0;
		mr->sel = 0;
	}
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( mr = rec->mark; mr; mr = mr->next ) {
				for ( m = mg->mark; m && m->id != mr->id; m = m->next ) ;
				if ( m ) {
					t = mr->loc - rec->origin;
/*					if ( mg->matrix[6] ) {
						mr->loc[2] += (m->loc[0] - mg->origin[0] - mg->matrix[0]*t[0] - mg->matrix[3]*t[1])/mg->matrix[6];
						mr->sel++;
					}
					if ( mg->matrix[7] ) {
						mr->loc[2] += (m->loc[1] - mg->origin[1] - mg->matrix[1]*t[0] - mg->matrix[4]*t[1])/mg->matrix[7];
						mr->sel++;
					}*/
					if ( mg->matrix[2][0] ) {
						mr->loc[2] += (m->loc[0] - mg->origin[0] - mg->matrix[0][0]*t[0] - mg->matrix[1][0]*t[1])/mg->matrix[2][0];
						mr->sel++;
					}
					if ( mg->matrix[2][1] ) {
						mr->loc[2] += (m->loc[1] - mg->origin[1] - mg->matrix[0][1]*t[0] - mg->matrix[1][1]*t[1])/mg->matrix[2][1];
						mr->sel++;
					}
				}
			}
		}
	}
	for ( mr = rec->mark; mr; mr = mr->next ) {
		if ( mr->sel ) mr->loc[2] /= mr->sel;
		mr->sel = 1;
	}
	
	return 0;
}

/**
@brief 	Generates missing micrograph and model markers.
@param 	*project	project parameter structure.
@return int			0, <0 on error.

	First the 3D model is updated to represent all markers from all micrographs.
	Then the missing markers in the micrographs are generated from the
	updated model.

**/
int			project_generate_markers(Bproject* project)
{
	if ( !project->rec ) return -1;
	
	Bfield*				field = project->field;
	Bmicrograph*		mg;
	
	if ( verbose )
		cout << "Generating missing markers" << endl << endl;

	for ( field = project->field; field; field = field->next ) if ( field->select )
		for ( mg = field->mg; mg; mg = mg->next )
			mg_update_model(mg, project->rec->mark, project->rec->origin);

	for ( field = project->field; field; field = field->next ) if ( field->select )
		for ( mg = field->mg; mg; mg = mg->next )
			mg_marker_update(mg, project->rec->mark, project->rec->origin, 0);
		
	return 0;
}

/**
@brief 	Erases or paints markers in an image.
@param 	*p				image parameter structure.
@param 	*mark			linked list of markers.
@param 	marker_radius	radius to mask out markers.
@return long			number of markers erased.

	Markers can be either erased to a background value or painted in with
	a set value.

**/
long		img_erase_markers(Bimage* p, Bmarker* mark, double marker_radius)
{
	if ( marker_radius < 1 ) return 0;
	
	long			nm(0);
	vector<double>	stats;
	for ( ; mark; mark = mark->next, nm++ ) {
		stats = p->stats_within_radii(0, mark->loc, marker_radius, 2*marker_radius);
		p->sphere(mark->loc, marker_radius, 2, FILL_USER, stats[3]);
	}
	
	return nm;
}

/**
@brief 	Erases or paints markers in micrographs.
@param 	*mg				micrograph parameter structure.
@param 	marker_radius	radius to mask out markers.
@return Bimage*			new micrograph image with erased markers.

	Markers can be either erased to a background value or painted in with
	a set value.

**/
Bimage*		mg_erase_markers(Bmicrograph* mg, double marker_radius)
{
	if ( verbose & VERB_FULL )
		cout << "Reading image " << mg->img_num << " (micrograph " << mg->id << ")" << endl;
	
	Bimage*			p = read_img(mg->fmg, 1, mg->img_num);
	if ( !p ) {
		error_show("mg_erase_markers", __FILE__, __LINE__);
		return NULL;
	}
	
	p->change_type(Float);
	
	p->sampling(mg->pixel_size);
	
	p->origin(mg->origin[0], mg->origin[1], 0);
	
	if ( marker_radius < 1 ) return p;

	img_erase_markers(p, mg->mark, marker_radius);
	
	return p;
}

/**
@brief 	Erases or paints markers in micrographs.
@param 	*project		project parameter structure.
@param 	marker_radius	radius to mask out markers.
@return int				0, <0 on error.

	Markers can be either erased to a background value or painted in with
	a set value.

**/
int			project_erase_markers(Bproject* project, double marker_radius)
{
	long 			n, nmg = project_count_micrographs(project);

	if ( nmg < 1 ) return -1;
	
	Bstring			filename;
	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bimage* 		p;
	Bimage* 		pall = NULL;
	
	if ( verbose )
		cout << "Removing markers with radius " << marker_radius << endl << endl;

	for ( field = project->field; field; field = field->next ) {
		nmg = field_count_micrographs(field);
		p = read_img(mg->fmg, 0, 0);
		if ( !p ) {
			error_show("project_erase_markers", __FILE__, __LINE__);
			return -1;
		}
		pall = new Bimage(p->data_type(), p->compound_type(), p->size(), nmg);
		delete p;
		if ( field->mg->pixel_size[0] ) pall->sampling(mg->pixel_size);
		filename = field->mg->fmg.pre_rev('.') + "_me." + field->mg->fmg.post_rev('.');
		for ( n=0, mg = field->mg; mg; mg = mg->next, n++ ) {
			if ( verbose )
				cout << "Reading image " << mg->img_num << " (micrograph " << mg->id << ")" << endl;
			p = mg_erase_markers(mg, marker_radius);
			if ( !p ) {
				error_show("project_remove_markers", __FILE__, __LINE__);
				return -1;
			}
			pall->replace(n, p);
			pall->image[n] = p->image[0]; 
			delete p;
			mg->fmg = filename;
		}
		write_img(filename, pall, 0);
		delete pall;
	}
	
	return 0;
}

/**
@brief 	Resets the model x and y coordinates from the zero-tilt micrograph.
@param 	*mg			reference micrograph parameter structure.
@param 	*model		model marker list.
@return int			markers selected in model.

	Missing model markers are added in and excessive model markers deleted.

**/
int			mg_reset_model(Bmicrograph* mg, Bmarker* model)
{
	if ( !mg->mark ) return 0;
	
	int			n(0), ns(0);
	Bmarker		*m, *mr;
	
	if ( verbose ) {
		cout << "Resetting 3D model x and y coordinates:" << endl;
		cout << "Micrograph:                     " << mg->id << endl << endl;
	}
	
	for ( mr = model; mr; mr = mr->next ) mr->sel = 0;
	
	for ( m = mg->mark; m; m = m->next ) {
		for ( mr = model; mr && m->id != mr->id; mr = mr->next ) ;
		if ( !mr ) {
			mr = (Bmarker *) add_item((char **) &model, sizeof(Bmarker));
			mr->id = m->id;
		}
		mr->loc[0] = m->loc[0];
		mr->loc[1] = m->loc[1];
		mr->sel = 1;
	}

	for ( mr = model; mr; mr = mr->next, n++ ) ns += mr->sel;

	if ( verbose )
		cout << "Markers in model:               " << n << " (" << ns << ")" << endl << endl;

//	markers_delete_non_selected(&model);
	
	return ns;
}

/**
@brief 	Finds all markers in a project.
@param 	*project	project parameter structure.
@param 	edge		edge to exclude.
@param 	add			flag to add rather than replace the markers.
@return long		number of markers.

	Markers are located in each micrograph by cross-correlation.

**/
long		project_find_markers(Bproject* project, long edge, int add)
{
	Bfield*			field = NULL;
	Bmicrograph*	mg = NULL;
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			mg_find_markers(mg, edge, add);
			
	return project_count_markers(project);
}

/**
@brief 	Finds markers in a micrograph.
@param 	*mg			project parameter structure.
@param 	edge		edge to exclude.
@param 	add			flag to add rather than replace the markers.
@return long		number of markers.

	The markers are located by cross-correlation with a synthetic reference
	whose size is defined by the marker radius. Markers close to the edge
	are eliminated using the given edge size. The main intent is to find 
	the seed for a 3D marker model in a zero degree tilt micrograph. 

**/
long		mg_find_markers(Bmicrograph* mg, long edge, int add)
{
	if ( mg->mark_radius < 1 ) {
		cerr << "Error: The marker radius is too low! (" << mg->mark_radius << ")" << endl;
		bexit(-1);
	}
	
	Bimage*			p = read_img(mg->fmg, 1, mg->img_num);
	if ( !p )
		return error_show("mg_find_markers", __FILE__, __LINE__);
	
	Bmarker*		mark = img_find_gold_particles(p, 0, mg->mark_radius, edge, -1);
	
	delete p;
	
	long			nmark = count_list((char *) mark);
	
	if ( nmark ) {
		if ( add ) {
			markers_add(&mg->mark, mark, 2*mg->mark_radius, 2);
		} else {
			kill_list((char *) mg->mark, sizeof(Bmarker));
			mg->mark = mark;
		}
	}

	nmark = count_list((char *) mg->mark);

	if ( verbose )
		cout << "Markers in " << mg->id << ": " << nmark << endl;
	
	return nmark;
}

/**
@brief 	Generates an average particle image from a micrograph.
@param 	*mg			micrograph to extract marker images from.
@param 	size		size of image.
@return Bimage*		composite marker projection image.
**/
Bimage*		mg_composite_particle(Bmicrograph* mg, Vector3<long> size)
{
	if ( !mg->mark ) {
		error_show("No markers found!", __FILE__, __LINE__);
		return NULL;
	}
	
	int				n;

	Bmarker*		mark;
	Bimage*			pmark;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mg_composite_particle: size = " << size << endl;

	Bimage*			p = read_img(mg->fmg, 1, mg->img_num);
	if ( !p ) {
		error_show("mg_composite_particle", __FILE__, __LINE__);
		return NULL;
	}
	
	Bimage*			pavg = new Bimage(Float, TSimple, size, 1);
	pavg->origin(pavg->size()/2);

	Vector3<double>	ghalf(pavg->image->origin());
	Vector3<double>	gorigin;

	for ( n=0, mark = mg->mark; mark; mark = mark->next, n++ ) if ( mark->sel > 0 ) {
		gorigin = mark->loc - ghalf;
		pmark = p->extract(0, gorigin, size);
		pmark->change_type(Float);
		pavg->add(pmark);
		delete pmark;
	}
	
	delete p;
		
	pavg->calculate_background();
	
	pavg->symmetrize_cylinder();
	
	pavg->statistics();
	
	pavg->calculate_background();
	
	return pavg;
}

/**
@brief 	Counts the markers in a project hierarchy.
@param 	*project	project parameter structure.
@return long 		number of markers.
**/
long		project_count_markers(Bproject* project)
{
	if ( !project ) return 0;
	if ( !project->field ) return 0;
	if ( !project->field->mg ) return 0;
	
	long				nmark(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				nmark += count_list((char *) mg->mark);
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			nmark += count_list((char *) rec->mark);
	}
	
	return nmark;
}

/**
@brief 	Shows the markers in a project hierarchy.
@param 	*project	project parameter structure.
@return long 		number of markers.
**/
long		project_show_markers(Bproject* project)
{
	if ( !project ) return 0;
	if ( !project->field ) return 0;
	if ( !project->field->mg ) return 0;
	
	long			i, nmark(0);
	Bfield* 		field;
	Bmicrograph*	mg;

	cout << "Project markers:" << endl;
	for ( field = project->field; field; field = field->next ) {
		cout << "Field: " << field->id << endl;
		for ( mg = field->mg; mg; mg = mg->next ) {
			i = count_list((char *) mg->mark);
			nmark += i;
			cout << "\tMicrograph: " << mg->id << tab << i << " markers\tselection = " << mg->select << endl;
		}
	}
	cout << "Total:\t" << nmark << " markers" << endl << endl;
	
	return nmark;
}

/**
@brief 	Shows the errors in marker coordinates in a project hierarchy.
@param 	*project		project parameter structure.
@param 	error_cutoff	threshold to show errors (distance in pixels).
@return long 			number of errors above threshold.
**/
long		project_show_errors(Bproject* project, double error_cutoff)
{
	if ( !verbose ) return 0;
	if ( !project ) return 0;
	
	long				nmark(0), nerr(0);
	Bfield*				field = project->field;
	Bmicrograph*		mg;
	Bmarker*			m;

	cout << "Project marker errors:" << endl;
	cout << "Micrograph\tMarker\tx\ty\tex\tey\ted" << endl;
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			for ( m = mg->mark; m; m = m->next ) if ( m->sel > 0 ) {
				if ( m->res > error_cutoff ) {
					cout << mg->id << tab << m->id << tab << m->loc[0] << tab << 
						m->loc[1] << tab << m->err[0] << tab << m->err[1] << tab << m->res << endl;
					nerr++;
				}
				nmark++;
			}
		}
	}
	cout << "Total above " << error_cutoff << ":\t" << nerr << " markers (" << nerr*100.0/nmark << "%)" << endl << endl;
	
	return nmark;
}

/**
@brief 	Deletes selected markers in a project hierarchy.
@param 	*project		project parameter structure.
@param 	&deselect_list	list of markers to deselect.
@return long			number of markers.

	All the occurrences of selected markers in the tilt series are deleted.

**/
long		project_deselect_markers(Bproject* project, Bstring& deselect_list)
{
	if ( !project ) return 0;
	
	long				maxid;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bmarker*			mark;
	
	for ( maxid=0, field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( mark = mg->mark; mark; mark = mark->next )
				if ( maxid < mark->id ) maxid = mark->id;

	for ( rec = project->rec; rec; rec = rec->next )
		for ( mark = rec->mark; mark; mark = mark->next )
			if ( maxid < mark->id ) maxid = mark->id;
	
	maxid++;
	
//	int*				sel = new int[maxid];

//	select_numbers(deselect_list, maxid, sel);
	vector<int>		sel = select_numbers(deselect_list, maxid);

	if ( verbose )
		cout << "Deselecting markers" << endl << endl;
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( mark = mg->mark; mark; mark = mark->next )
				if ( sel[mark->id] ) mark->sel = 0;
		}
	}

	for ( rec = project->rec; rec; rec = rec->next ) {
		for ( mark = rec->mark; mark; mark = mark->next ) {
			if ( sel[mark->id] ) mark->sel = 0;
		}
	}
	
//	delete[] sel;
	
	long			n = count_list((char *) project->field->mg->mark);
	
	return n;
}

/**
@brief 	Deletes selected markers in a project hierarchy.
@param 	*project		project parameter structure.
@param 	&delete_list	list of markers to deselect.
@return long			number of markers.

	All the occurrences of selected markers in the tilt series are deleted.

**/
long		project_delete_markers(Bproject* project, Bstring& delete_list)
{
	if ( !project ) return 0;
	
	long				maxid;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bmarker*			mark;
	
	for ( maxid=0, field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( mark = mg->mark; mark; mark = mark->next )
				if ( maxid < mark->id ) maxid = mark->id;
	
	maxid++;
	
//	int*				sel = new int[maxid];

//	select_numbers(delete_list, maxid, sel);
	vector<int>		sel = select_numbers(delete_list, maxid);

	if ( verbose )
		cout << "Deleting markers:" << endl;
	for ( field = project->field; field; field = field->next ) {
		if ( verbose )
			cout << "Field: " << field->id << ":" << endl;
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( verbose )
				cout << "\tMicrograph: " << mg->id << endl;
			for ( mark = mg->mark; mark; ) {
				if ( sel[mark->id] ) {
					if ( verbose )
						cout << "\t\tDeleting marker: " << mark->id << endl;
					mark = (Bmarker *) remove_item((char **)&mg->mark, (char *)mark, sizeof(Bmarker));
				} else {
					mark = mark->next;
				}
			}
		}
	}

	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( verbose )
			cout << "Reconstruction: " << rec->id << endl;
		for ( mark = rec->mark; mark; ) {
			if ( sel[mark->id] ) {
				if ( verbose )
					cout << "\t\tDeleting marker: " << mark->id << endl;
				mark = (Bmarker *) remove_item((char **)&rec->mark, (char *)mark, sizeof(Bmarker));
			} else {
				mark = mark->next;
			}
		}
	}
	
//	delete[] sel;
	
	long			n = count_list((char *) project->field->mg->mark);
	
	return n;
}

/**
@brief 	Renumbers the markers in a project hierarchy.
@param 	*project		project parameter structure.
@return long			number of markers.

	The markers are assumed to correspond across micrographs and reconstructions.
	The existing marker ids are mapped to an array, and new marker ids generated.

**/
long		project_renumber_markers(Bproject* project)
{
	if ( !project ) return 0;
	
	long				maxid, i, j;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bmarker*			mark;
	
	// Find maximum marker id
	for ( maxid=0, field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( mark = mg->mark; mark; mark = mark->next )
				if ( maxid < mark->id ) maxid = mark->id;

	for ( rec = project->rec; rec; rec = rec->next )
		for ( mark = rec->mark; mark; mark = mark->next )
			if ( maxid < mark->id ) maxid = mark->id;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_renumber_markers: maxid=" << maxid << endl;
		
	maxid++;
	
	// Establish all existing marker ids
	int*				id = new int[maxid];
	for ( i=0; i<maxid; i++ ) id[i] = 0;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( mark = mg->mark; mark; mark = mark->next )
				id[mark->id]++;

	for ( rec = project->rec; rec; rec = rec->next )
		for ( mark = rec->mark; mark; mark = mark->next )
			id[mark->id]++;
	
	// Create a mapping from old to new ids
	for ( i=j=0; i<maxid; i++ )
		if ( id[i] ) id[i] = ++j;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_renumber_markers: new maxid=" << j << endl;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( mark = mg->mark; mark; mark = mark->next )
				mark->id = id[mark->id];

	for ( rec = project->rec; rec; rec = rec->next )
		for ( mark = rec->mark; mark; mark = mark->next )
			mark->id = id[mark->id];	
	
	delete[] id;
	
	return j;
}

/**
@brief 	Selects micrographs in a project hierarchy.
@param 	*project	project parameter structure.
@param 	&mg_select	list of micrograph numbers.
@return long 		number of micrographs selected.

	A micrograph is only selected when it has its selection already set
	and it is in the list of numbers.

**/
long		project_mg_select(Bproject* project, Bstring& mg_select)
{
	if ( !project ) return 0;
	if ( mg_select.length() < 1 ) return 0;
	
	long				n(0);
	Bfield*				field = project->field;
	Bmicrograph*		mg;
	
	vector<long>		en = mg_select.split_into_integers(",");

	if ( verbose & VERB_PROCESS )
		cout << "Selecting micrographs:" << endl;
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			mg->select = 0;
			for ( size_t i=0; i<en.size(); i++ ) {
				if ( mg->img_num == en[i] ) {
					mg->select = 1;
					n++;
					if ( verbose & VERB_PROCESS )
						cout << "Micrograph " << mg->id << " (" << mg->img_num << ") selected" << endl;
				}
			}
		}
	}
	if ( verbose & VERB_PROCESS )
		cout << "Total Selected:                 " << n << endl << endl;

	return n;
}

/**
@brief 	Deselects micrographs in a project hierarchy.
@param 	*project	project parameter structure.
@param 	&mg_exclude	list of micrograph numbers.
@return long 		number of micrographs excluded.
**/
long		project_mg_exclude(Bproject* project, Bstring& mg_exclude)
{
	if ( !project ) return 0;
	
	long				i, nmg(0), n(0);
	Bfield*				field = project->field;
	Bmicrograph*		mg;
	
	if ( mg_exclude.contains("none") ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next ) mg->select = 1;
		return 0;
	}

	for ( field = project->field; field; field = field->next )
		for ( nmg=0, mg = field->mg; mg; mg = mg->next )
			if ( nmg <= mg->img_num ) nmg++;
	
//	int*			en = new int[nmg];

//	select_numbers(mg_exclude, nmg, en);
	vector<int>		en = select_numbers(mg_exclude, nmg);
	
	if ( verbose )
		cout << "Deselecting micrographs:" << endl;
	for ( field = project->field; field; field = field->next ) {
		for ( i=0, mg = field->mg; mg; mg = mg->next, i++ ) if ( en[i] ) {
			mg->select = 0;
			n++;
			if ( verbose & VERB_PROCESS )
				cout << "Micrograph " << mg->id << " (" << mg->img_num << ") deselected" << endl;
		}
	}
	if ( verbose )
		cout << "Total deselected:               " << n << endl << endl;

//	delete[] en;
	
	return n;
}

/**
@brief 	Sets the tilt axis angle in a project parameter structure.
@param 	*project	project parameter structure.
@param 	tilt_axis	tilt axis angle (wrt x-axis).
@return int			0.

	The tilt axis angle for all the micrographs are set to the same value.

**/
int			project_set_tilt_axis(Bproject* project, double tilt_axis)
{
	Bfield*			field = project->field;
	Bmicrograph*	mg;
	
	for ( field = project->field; field; field = field->next ) if ( field->select )
		for ( mg = field->mg; mg; mg = mg->next )
			mg->tilt_axis = tilt_axis;
	
	project_mg_tilt_to_matrix(project);
	
	return 0;
}

/**
@brief 	Inverts the tilt axis in a project parameter structure.
@param 	*project	project parameter structure.
@return int			0.

	The tilt axis angle for all the micrographs are set to the same value.

**/
int			project_invert_tilt_axis(Bproject* project)
{
	Bfield*			field = project->field;
	Bmicrograph*	mg;
	
	for ( field = project->field; field; field = field->next ) if ( field->select )
		for ( mg = field->mg; mg; mg = mg->next )
			mg->tilt_axis = angle_set_negPI_to_PI(mg->tilt_axis + M_PI);
	
	project_mg_tilt_to_matrix(project);
	
	return 0;
}

/**
@brief 	Sets the tilt and tilt axis angles in a project parameter structure.
@param 	*project	project parameter structure.
@param 	tilt_start	starting tilt angle (usually negative).
@param 	tilt_step	tilt increment angle.
@return int			0.

	The single tilt series is defined by a starting tilt angle and a
	tilt increment angle, as well as the tilt axis angle.

**/
int			project_set_tilt_angles(Bproject* project, double tilt_start, double tilt_step)
{
	double			tilt;
	Bfield*			field = project->field;
	Bmicrograph*	mg;
	
	if ( fabs(tilt_step) < 0.005 ) {	// Must be greater than 0.5 degrees
		cerr << "Error: The tilt step size is too small (" << tilt_step << ")!" << endl;
		return -1;
	}
	
	for ( field = project->field; field; field = field->next ) if ( field->select )
		for ( tilt=tilt_start, mg = field->mg; mg; mg = mg->next, tilt+=tilt_step )
			mg->tilt_angle = tilt;

	project_mg_tilt_to_matrix(project);
	
	return 0;
}

/**
@brief 	Sets the matrix from the tilt and axis angles for each micrograph.
@param 	*project	project parameter structure.
@return int			0.

	This is in preparation for reconstruction.

**/
int			project_mg_tilt_to_matrix(Bproject* project)
{
	Bfield* 		field;
	Bmicrograph*	mg;
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			mg->matrix = Matrix3(mg->tilt_angle, mg->tilt_axis).transpose();
//			mg->matrix = mg->matrix.transpose();
		}
	
	return 0;
}

/**
@brief 	Calculates tilt and tilt axis angles from matrices in a project parameter structure.
@param 	*project	project parameter structure.
@return int			0.

	From each matrix, a quaternion is calculated, giving the tilt axis
	and rotation angle.
	The level angle is calculated as the arcsin of the z-coordinate of the axis.

**/
int			project_calculate_angles(Bproject* project)
{
	double			first_axis_angle;
//	Vector3<double>	tilt_axis;
//	Quaternion		q;
	Transform		t;
	Matrix3			mat;
	Bfield*			field = project->field;
	Bmicrograph*	mg;

	Bmicrograph*	mg_ref = field_find_zero_tilt_mg(field);
	
	if ( field->next ) {
		mg = field_find_zero_tilt_mg(field->next);
		t = marker_find_transform(mg->mark, mg_ref->mark, mg_ref->origin);
		field->next->matrix = Matrix3(t.axis, t.angle);
	}
	
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		if ( verbose & VERB_FULL ) {
			cout << "Field " << field->id << " matrix:" << endl << field->matrix << endl;
			cout << "Mg\tAxis\tTilt" << endl;
		}
		first_axis_angle = 10;
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
//			if ( verbose & VERB_FULL )
//				cout << "Micrograph " << mg->id << " matrix:" << endl << mg->matrix << endl;
			if ( mg->matrix.determinant() > 0.5 ) {
				mat = field->matrix.transpose() * mg->matrix;
				vector<double>	a = mat.tilt_angles();
				if ( first_axis_angle > 9 ) first_axis_angle = a[0];
				if ( fabs(a[0] - first_axis_angle) > M_PI_2 ) {
					mg->tilt_axis = a[0] + M_PI;
					if ( mg->tilt_axis > M_PI ) mg->tilt_axis -= TWOPI;
					mg->tilt_angle = a[1];
				} else {
					mg->tilt_axis = a[0];
					mg->tilt_angle = -a[1];
				}
				mg->level_angle = a[2];
/*				q = mat.quaternion();
				cout << "quaternion = " << q << endl;
				tilt_axis = q.axis();
				axis_angle = atan2(tilt_axis[1], tilt_axis[0]);
				if ( !isfinite(axis_angle) ) axis_angle = 0;
				if ( first_axis_angle > 9 ) first_axis_angle = axis_angle;
				if ( fabs(axis_angle - first_axis_angle) > M_PI_2 ) {
					mg->tilt_axis = -axis_angle;
					mg->tilt_angle = 2*acos(q[0]);
				} else {
					mg->tilt_axis = axis_angle;
					mg->tilt_angle = -2*acos(q[0]);
				}
				// Calculate level angle
				mg->level_angle = asin(tilt_axis[2]);
				if ( !isfinite(mg->level_angle) ) mg->level_angle = 0;
				if ( mg->level_angle > 0.2 ) mg->level_angle = 0;*/
			}
			if ( verbose & VERB_FULL )
				cout << mg->id << tab << mg->tilt_axis*180/M_PI
					<< tab << mg->tilt_angle*180/M_PI << endl;
		}
		mg = field_find_zero_tilt_mg(field);
		// Fix the zero-tilt axis
		if ( fabs(mg->tilt_axis - first_axis_angle) > M_PI/90.0 ) mg->tilt_axis = first_axis_angle;
	}

	return 0;
}

/**
@brief 	Inverts the micrograph matrices in a project parameter structure.
@param 	*project	project parameter structure.
@return int			0.

	The tilt axis angle for all the micrographs are set to the same value.

**/
int			project_invert_matrices(Bproject* project)
{
	Bfield*			field = project->field;
	Bmicrograph*	mg;
	
	for ( field = project->field; field; field = field->next ) if ( field->select )
		for ( mg = field->mg; mg; mg = mg->next )
			mg->matrix = mg->matrix.transpose();
	
	return 0;
}


/**
@brief 	Selects markers.
@param 	*project	project parameter structure.
@param 	fom			fom cutoff.
@return int			0.

	Selects all micrograph markers above the given fom cutoff.

**/
int			project_mg_marker_select(Bproject* project, double fom)
{
	Bfield*				field;
	Bmicrograph*		mg;
	Bmarker*			mark;
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( mark = mg->mark; mark; mark = mark->next )
				if ( mark->fom >= fom ) mark->sel = 1;
				else mark->sel = 0;
	
	return 0;
}

/**
@brief 	Sets the marker radius in a project parameter structure.
@param 	*project	project parameter structure.
@param 	mark_radius	gold fiducial marker radius (voxels).
@return int			0.
**/
int			project_set_marker_radius(Bproject* project, double mark_radius)
{
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			mg->mark_radius = mark_radius;
	
	for ( rec = project->rec; rec; rec = rec->next )
		rec->mark_radius = mark_radius;
	
	return 0;
}

/**
@brief 	Checks if the markers in a project parameter structure fall within micrographs.
@param 	*project	project parameter structure.
@param 	flags		flags to modify actions.
@return int			0.

	Markers outside the boundaries are dealt with based on the flags variable:
	1	show markers outside the image boundaries.
	2	set marker selections to zero.
	4	set marker fom's to zero.
	8	set marker errors to zero.

**/
int			project_check_markers(Bproject* project, int flags)
{
	Bfield*				field = project->field;
	Bmicrograph*		mg = field_find_micrograph_by_tiltang(field, 0);
	Bmarker*			m;
	Bimage*				p = read_img(mg->fmg, 0, mg->img_num);
	
	if ( flags & 1 ) {
		cout << "Markers outside micrograph boundaries:" << endl;
		cout << "Micrograph\tMarker\tx\ty" << endl;
	}
	
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( m = mg->mark; m; m = m->next ) {
				if ( m->loc[0] < 0 || m->loc[0] > p->sizeX() || m->loc[1] < 0 || m->loc[1] > p->sizeY() ) {
					if ( flags & 1 ) cout << mg->id << tab << m->id << tab << m->loc[0] << tab << m->loc[1] << endl;
					if ( flags & 2 ) m->sel = 0;
					if ( flags & 4 ) m->fom = 0;
					if ( flags & 8 ) m->err = m->res = 0;
				}
			}
		}
	}

	delete p;
	
	return 0;
}

/**
@brief 	Checks if the markers in a project parameter structure fall within micrographs.
@param 	*project	project parameter structure.
@return int			0.
**/
int			project_fix_markers(Bproject* project)
{
	Bfield*				field = project->field;
	Bmicrograph*		mg = field_find_micrograph_by_tiltang(field, 0);
	Breconstruction*	rec = project->rec;
	Bmarker*			m, *mm;
	Bimage*				p = read_img(mg->fmg, 0, mg->img_num);
	
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( m = mg->mark, mm = rec->mark; m; m = m->next, mm = mm->next ) {
				if ( m->loc[0] < 0 || m->loc[0] > p->sizeX() || m->loc[1] < 0 || m->loc[1] > p->sizeY() ) {
					m->loc = mg_location_from_3D_model(mm->loc, rec->origin, mg->matrix, mg->origin);
					m->err = 0;
					m->sel = 0;
					m->fom = m->res = 0;
				}
			}
		}
	}
	
	delete p;
	
	return 0;
}

/**
@brief 	Clears extraneous areas from a micrograph.
@param 	*p			micrograph image.
@param 	tilt_axis	tilt axis angle (radians).
@param 	tilt_angle	tilt angle (radians).
@param 	thickness	intended thickness of reconstruction.
@param 	width		edge smoothing width.
@return int			0.
**/
int			img_clear_extraneous_areas(Bimage* p, double tilt_axis, double tilt_angle,
				long thickness, double width)
{
	if ( p->sizeZ() > 1 ) return -1;
	if ( p->images() > 1 ) return -1;
	
	long				i, xx, yy, cc;
	double				de(0), dr(0), edge, a(GOLDEN/width);
    Vector3<double>		d, r, hr(p->size()/2-3*width);

	p->calculate_background();
	
	double			fill(p->background(long(0)));

	Vector3<double>	vec = Vector3<double>(-sin(tilt_axis), cos(tilt_axis), 0);
	double			dist(1e30);
	if ( vec[0] && dist > fabs(p->sizeX()*0.5/vec[0]) ) dist = fabs(p->sizeX()*0.5/vec[0]);
	if ( vec[1] && dist > fabs(p->sizeY()*0.5/vec[1]) ) dist = fabs(p->sizeY()*0.5/vec[1]);
	dist = dist*cos(tilt_angle) + fabs(0.5*thickness*sin(tilt_angle));

	for ( i=yy=0; yy<p->sizeY(); yy++ ) {
		d[1] = (double)yy - p->image->origin()[1];
		for ( xx=0; xx<p->sizeX(); xx++ ) {
			d[0] = (double)xx - p->image->origin()[0];
			r = d.abs() - hr;
			if ( r[0] < r[1] ) dr = a*r[1];
			else dr = a*r[0];
			de = a*(fabs(vec[0]*d[0] + vec[1]*d[1]) - dist);
			if ( dr > de ) de = dr;
			if ( de > 50 ) edge = 1e30;
			else edge = exp(de);
			if ( !isfinite(edge) ) {
				cerr << xx << " " << yy << ": Value too large or not finite!: " << edge <<endl;
				edge = 1e30;
			}
			for ( cc=0; cc<p->channels(); cc++, i++ )
				p->set(i, ((*p)[i] + fill*edge)/(1+edge));
		}
	}
	
	return 0;
}

/**
@brief 	Clears extraneous areas from a micrograph.
@param 	*mg			micrograph parameter structure.
@param 	*p			micrograph image.
@param 	thickness	intended thickness of reconstruction.
@param 	width		edge smoothing width.
@return int			0.
**/
int			micrograph_clear_extraneous_areas(Bmicrograph* mg, Bimage* p, long thickness, double width)
{
	if ( !p )
		return error_show("micrograph_clear_extraneous_areas", __FILE__, __LINE__);
		
	View			view = View(mg->matrix);
	p->origin(mg->origin);
	p->image->view(view);
	
	img_clear_extraneous_areas(p, mg->tilt_axis, mg->tilt_angle, thickness, width);
	
	return 0;
}

/**
@brief 	Clears areas from a tilt series not used in a reconstruction.
@param 	*project	project parameter structure.
@param 	thickness	intended thickness of reconstruction.
@param 	width		edge smoothing width.
@return int			0.
**/
int			project_clear_extraneous_areas(Bproject* project, long thickness, double width)
{
	long				nmg(0);
	Bimage*				p = NULL;
	Bimage*				pnu = NULL;
	Bfield*				field = project->field;
	Bmicrograph*		mg = field->mg;

	Bstring				newfile = mg->fmg;
	newfile = newfile.insert(newfile.rfind("."), "_clr");
	
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( nmg=0, mg = field->mg; mg; mg = mg->next, nmg++ ) {
			if ( verbose ) 
				cout << "Micrograph " << mg->img_num << ": " << mg->fmg << endl;
			p = read_img(mg->fmg, 0, mg->img_num);
			if ( !p )
				return error_show("project_clear_extraneous_areas", __FILE__, __LINE__);
			if ( !pnu ) {
				pnu = p;
			} else {
				delete p;
			}
		}
	}
	
	if ( verbose ) {
		cout << "Clearing extraneous areas from micrographs:" << endl;
		cout << "Thickness:                      " << thickness << endl;
		cout << "Edge width:                     " << width << endl;
		cout << "New micrograph file name:       " << newfile << endl << endl;
	}

	pnu->images(nmg);	
	pnu->data_alloc();

//	mg = project->field->mg;

/*	View			view;
	Vector3<double>	v(-sin(mg->tilt_axis), cos(mg->tilt_axis), 0);
	double			a, d(1e30);
	if ( v[0] && d > fabs(pnu->sizeX()*0.5/v[0]) ) d = fabs(pnu->sizeX()*0.5/v[0]);
	if ( v[1] && d > fabs(pnu->sizeY()*0.5/v[1]) ) d = fabs(pnu->sizeY()*0.5/v[1]);*/
	
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( nmg=0, mg = field->mg; mg; mg = mg->next, nmg++ ) {
			p = read_img(mg->fmg, 1, mg->img_num);
			mg->fmg = newfile;
//			micrograph_clear_extraneous_areas(mg, p, thickness, width);
			p->origin(mg->origin);
			p->image->view(View(mg->matrix));
			img_clear_extraneous_areas(p, mg->tilt_axis, mg->tilt_angle, thickness, width);
/*			view = View(mg->matrix);
			if ( verbose ) 
				cout << "Micrograph " << nmg+1 << ": " << mg->fmg << endl;
			if ( !p )
				return error_show("project_clear_extraneous_areas", __FILE__, __LINE__);
			p->calculate_background();
			p->origin(mg->origin);
			v = Vector3<double>(-sin(mg->tilt_axis), cos(mg->tilt_axis), 0);
			a = d*cos(mg->tilt_angle) + fabs(0.5*thickness*sin(mg->tilt_angle));
			img_clear_extraneous_areas(p, v, a, width);*/
			pnu->replace(nmg, p);
			pnu->image[nmg] = p->image[0];
			pnu->image[nmg].origin(mg->origin);
			pnu->image[nmg].view(View(mg->matrix));
			delete p;
		}
		write_img(field->mg->fmg, pnu, 0);
	}
		
	delete pnu;

	return 0;
}

/**
@brief 	Determines the most likely rotation axis location for a tilt series.
@param 	*project	project parameter structure.
@return double		average z offset from marker z average.

	All micrographs must have associated markers.

**/
double		project_marker_rotation_axis(Bproject* project)
{
	long	n;
	double			avg_z_offset(0), zs(0), d;
	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bmicrograph*	mg_ref = field_find_low_tilt_mg_with_markers(field);

	for ( n=0, mg = field->mg; mg; mg = mg->next ) {
		if ( mg->tilt_angle > M_PI/180.0 ) {
			zs += fabs(mg->origin.distance(mg_ref->origin)/sin(mg->tilt_angle));
			n++;
		}
	}

	if ( n ) avg_z_offset = zs/n;
	
	cout << "Average z offset:               " << avg_z_offset << endl << endl;

	cout << "Micrograph\tOffset\tAdjusted" << endl;
	for ( mg = field->mg; mg; mg = mg->next ) {
		d = mg->origin.distance(mg_ref->origin);
		cout << mg->id << tab << d << tab << d - fabs(avg_z_offset*sin(mg->tilt_angle)) << endl;
	}
	cout << endl;

	return avg_z_offset;
}

/**
@brief 	Merges the markers for two reconstructions.
@param 	*project	project parameter structure.
@return int			number of markers merged.

	The reconstructions must have corresponding markers.
	The marker locations are set to the average of the two.
	The second reconstruction is removed.

**/
int			project_merge_rec_markers(Bproject* project)
{
	Breconstruction*	rec1 = project->rec;
	Breconstruction*	rec2 = project->rec->next;
	Bmarker*			mark1 = rec1->mark;
	Bmarker*			mark2 = rec2->mark;
		
	long				n;
	double				R(0);

	if ( verbose & VERB_FULL )
		cout << "Marker\tLocation\tResidual" << endl;
	for ( n=0, mark1 = rec1->mark; mark1; mark1 = mark1->next ) if ( mark1->sel ) {
		for ( mark2 = rec2->mark; mark2 && mark2->id != mark1->id; mark2 = mark2->next ) ;
		if ( mark2 && mark2->sel ) {
			mark1->err = mark1->loc - mark2->loc;
			mark1->res = mark1->err.length();
			R += mark1->res*mark1->res;
			n++;
			mark1->loc = (mark1->loc + mark2->loc)/2;	// Set the marker model to the average
			if ( verbose & VERB_FULL )
				cout << mark1->id << tab << mark1->loc << tab << mark1->res << endl;
		}
	}
	
	if ( n ) R = sqrt(R/n);
	
	if ( verbose & VERB_FULL ) {
		cout << "Markers selected:               " << n << endl;
		cout << "Residual:                       " << R << endl << endl;
	}
	
	reconstruction_kill(rec2);
	rec1->next = NULL;
	
	return n;
}


