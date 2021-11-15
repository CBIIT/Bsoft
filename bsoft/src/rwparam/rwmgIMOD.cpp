/**
@file	rwmgIMOD.cpp
@brief	Converts between IMOD files and a micrograph parameter file
@author	Bernard Heymann
@date	Created: 20070501
@date	Modified: 20170120
**/

#include "rwmgIMOD.h"
#include "mg_img_proc.h"
#include "mg_tomography.h"
#include "rwmg.h"
#include "rwimg.h"
#include "linked_list.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


int			read_imod_xf(Bstring& filename, Bproject* project)
{
	Bmicrograph*		mg = NULL;
	double				a[4], dx, dy;

	ifstream			fprm(filename.c_str());
	if ( fprm.fail() ) {
		cerr << "Error: A text file specifying 2D transformations must be given!" << endl;
		return -1;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Reading 2D transformations from " << filename << endl;
	
	for ( mg = project->field->mg; mg; mg = mg->next ) {
		fprm >> a[0] >> a[1] >> a[2] >> a[3] >> dx >> dy;
		mg->tilt_axis = (atan2(-a[0], a[1]) + atan2(-a[3], -a[2]))/2.0;
		mg->origin[0] -= a[0]*dx + a[2]*dy;
		mg->origin[1] -= a[1]*dx + a[3]*dy;
		if ( verbose )
			cout << mg->origin << tab << mg->tilt_axis*180.0/M_PI << endl;
	}
	
	fprm.close();
	
	return 0;
}

int			read_imod_tlt(Bstring& filename, Bproject* project)
{
	Bmicrograph*		mg = NULL;
	double				ta;

	ifstream			fprm(filename.c_str());
	if ( fprm.fail() ) {
		cerr << "Error: A text file specifying tilt angles must be given!" << endl;
		return -1;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Reading tilt angles from " << filename << endl;
	
	for ( mg = project->field->mg; mg; mg = mg->next ) {
		fprm >> ta;
		mg->tilt_angle = ta*M_PI/180.0;
		if ( verbose )
			cout << mg->tilt_angle*180.0/M_PI << endl;
	}
	
	fprm.close();
	
	return 0;
}

int			read_imod_xyz(Bstring& filename, Breconstruction* rec)
{
	int					mark_id(0), first(1);
	double				x, y, z, dum;
	string				word;
	Bmarker*			mark = NULL;

	ifstream			fprm(filename.c_str());
	if ( fprm.fail() ) {
		cerr << "Error: A text file specifying the 3D marker model must be given!" << endl;
		return -1;
	}
	
	fprm >> mark_id >> x >> y >> z >> dum >> dum;
	while ( fprm >> word ) {
		if ( word == "Pix:" ) fprm >> rec->voxel_size[0];
		if ( word == "Dim:" ) break;
	}
	fprm >> rec->origin[0] >> rec->origin[1];
	rec->origin /= 2;
	
	if ( access(rec->frec.c_str(), R_OK) ) {
		if ( verbose )
			cout << "Getting reconstruction parameters from " << filename << endl;
	} else {
		if ( verbose )
			cout << "Getting reconstruction parameters from " << rec->frec << endl;
		Bimage*		p = read_img(rec->frec, 0, -1);
		rec->voxel_size = p->sampling(0);
		rec->origin = p->size()/2;
		delete p;
	}
	if ( verbose )
			cout << setprecision(3) << "Voxel size = " << rec->voxel_size << " Origin = " << rec->origin << endl;
	
	while ( !fprm.eof() ) {
		if ( first ) first = 0;
		else fprm >> mark_id >> x >> y >> z >> dum >> dum;
		if ( fprm.fail() ) break;
		mark = (Bmarker *) add_item((char **) &mark, sizeof(Bmarker));
		if ( !rec->mark ) rec->mark = mark;
		mark->id = mark_id;
//		mark->loc = Vector3<double>(x - rec->origin[0], y - rec->origin[1], z);
		mark->loc = Vector3<double>(x, y, z);
		mark->fom = 1;
		mark->sel = 1;
		if ( verbose )
			cout << mark_id << tab << x << tab << y << tab << z << endl;
	}
	
	fprm.close();

	return 0;
}

int			read_imod_fid(Bstring& filename, Bproject* project)
{
	Bmicrograph*		mg = NULL;
	Bmarker*			mark = NULL;
	
	ifstream			fprm(filename.c_str());
	if ( fprm.fail() ) {
		cerr << "Error: A binary IMOD file specifying the fiducial markers must be given!" << endl;
		return -1;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Reading markers from " << filename << endl;
	
	long				i, j, psize;
	char				tag[8];
	IMOD*				imod = new IMOD;
	OBJT*				objt = new OBJT;
	CONT*				cont = new CONT;
	float*				point;

	tag[4] = 0;
	fprm.read((char *)tag, 4);
	if ( fprm.fail() ) return -2;
	cout << tag << tab;
	fprm.read((char *)tag, 4);
	if ( fprm.fail() ) return -2;
	cout << tag << endl;

	fprm.read((char *)imod, IMODSIZE);
	if ( fprm.fail() ) return -2;

	swapbytes((unsigned char *) &imod->objsize, 4);
	
	if ( verbose & VERB_PROCESS )
		cout << "number of objects = " << imod->objsize << endl;
	
	fprm.read((char *)tag, 4);
	if ( fprm.fail() ) return -2;

	fprm.read((char *)objt, OBJTSIZE);
	if ( fprm.fail() ) return -2;

	swapbytes((unsigned char *) &objt->contsize, 4);
	
	if ( verbose & VERB_PROCESS )
		cout << "number of contours = " << objt->contsize << endl;
	
	for ( i=0; i<objt->contsize; i++ ) {
		fprm.read((char *)tag, 4);
		if ( fprm.fail() ) return -2;

		fprm.read((char *)cont, CONTSIZE);
		if ( fprm.fail() ) return -2;
		
		swapbytes((unsigned char *) &cont->psize, 4);
	
		if ( verbose & VERB_PROCESS )
			cout << "number of points = " << cont->psize << endl;
	
		psize = 3*cont->psize*sizeof(float);
		point = new float[3*cont->psize];

		fprm.read((char *)point, psize);
		if ( fprm.fail() ) return -2;
		
		swapbytes(psize, (unsigned char *) point, 4);
		
		for ( j=0, mg = project->field->mg; mg; mg = mg->next, j+=3 ) {
			mark = (Bmarker *) add_item((char **) &mg->mark, sizeof(Bmarker));
			mark->id = i+1;
			mark->loc = Vector3<double>(point[j], point[j+1], 0);
			mark->fom = 1;
			mark->sel = 1;
			if ( verbose )
				cout << mg->id << tab << mark->id << tab << mark->loc << endl;
		}
		
		delete[] point;
	}
	
	delete imod;
	delete objt;
	delete cont;

	fprm.close();

	return 0;
}

/**
@brief 	Creates a project structure using IMOD parameters.
@param 	*file_list		list of IMOD parameter file names.
@param 	*project		initialized project structure.
@param 	flag			flag to indicate conversion of the image file.
@return int				error code (<0 means failure).

	Requirements:
		Tilt series micrograph image (if 3D then converted to multi-2D)
		2D transform file (.xf)
		Tilt angle file (.tlt)
	Calculations:
		Tilt angle: from tilt angle file
		Tilt axis: ta = (atan2(-A11, A12) + atan2(-A22, -A21))/2
		Mg origin: o = on - Ad
	where
		A:  2x2 transformation matrix
		ta: tilt axis angle
		on: nominal micrograph origin (center of mg)
		o:  aligned micrograph origin
		d:  shift for micrograph

**/
int			read_project_imod(Bstring* file_list, Bproject* project, int flag)
{
	if ( file_list->empty() ) {
		cerr << "No filenames given to create a project!" << endl;
		exit(-1);
	}

	Bimage* 			p = NULL;
	Bstring				base, file, imgfile;
	Bstring*			thisfile = NULL;
	Bstring				ext;
	int					n;
	Bstring				mg_id, rec_id;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	

	// Set up the overall project based on the provided images
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		file = thisfile->pre_rev(':');
		ext = file.extension();
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_project_imod: reading " << file << endl;
		if ( ext.contains("st") || ext.contains("ali") || ext.contains("rec") ) {
			imgfile = file;
			base = imgfile.base();
			p = read_img(imgfile, 0, -1);
			if ( !p ) {
				cerr << "Error: Image file " << imgfile << " not read!" << endl;
				error_show("read_project_imod", __FILE__, __LINE__);
				return -1;
			}
			if ( ext.contains("rec") ) {
				rec = reconstruction_add(&project->rec, base);
				rec->frec = imgfile;
				rec->voxel_size = p->sampling(0);
				rec->origin = p->size()/2;
			} else {
				if ( p->sizeZ() > 1 ) p->slices_to_images();
				field = field_add(&project->field, base);
				for ( n=0; n<p->images(); n++ ) {
					mg_id = base + Bstring(n+1, "_%03d");
					mg = micrograph_add(&mg, mg_id);
					if ( !field->mg ) field->mg = mg;
					mg->block = n;
					mg->img_num = n;
					mg->fmg = imgfile;
					mg->pixel_size = p->sampling(0);
					mg->intensity = p->image[n].average();
					mg->origin = p->size()/2;
					mg->mark_radius = 10;
				}
				if ( verbose & VERB_FULL ) {
					cout << "Creating field " << field->id << ", with " << n << " micrographs" << endl;
					cout << "Pixel size:                             " << field->mg->pixel_size << endl;
					cout << "Origin:                                 " << field->mg->origin << endl;
				}
			}
			delete p;
		}
	}

	if ( imgfile.length() < 1 ) {
		cerr << "Error: No image file provided!" << endl;
		error_show("read_project_imod", __FILE__, __LINE__);
		return -2;
	}
		
	int				err(0);
	
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		ext = thisfile->extension();
		if ( ext.contains("xf") ) {
			err += read_imod_xf(*thisfile, project);
		} else if ( ext.contains("tlt") ) {
			err += read_imod_tlt(*thisfile, project);
		} else if ( ext.contains("xyz") ) {
			if ( !rec ) {
				rec = reconstruction_add(&project->rec, base);
				rec->frec = base + ".rec";
			}
			err += read_imod_xyz(*thisfile, rec);
		} else if ( ext.contains("fid") ) {
			err += read_imod_fid(*thisfile, project);
		}
	}

	project_mg_tilt_to_matrix(project);
	
	return err;
}

int			write_imod_xf(Bstring& imodfile, Bproject* project)
{
	Bstring				paramfile = imodfile.base() + ".xf";
	Bmicrograph*		mg = project->field->mg;
	Vector3<double>		d, origin = micrograph_get_nominal_origin(mg);

	ofstream			fprm(paramfile.c_str());
	if ( fprm.fail() ) {
		error_show("write_imod_xf", __FILE__, __LINE__);
		return -1;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Writing 2D transformations to " << paramfile << endl;
	
	double				a[4];
	
	for ( mg = project->field->mg; mg; mg = mg->next ) {
		a[0] = a[3] = -sin(mg->tilt_axis);
		a[1] = cos(mg->tilt_axis);
		a[2] = -cos(mg->tilt_axis);
		d = origin - mg->origin;
		fprm << setprecision(7);
		for ( int i=0; i<4; ++i ) fprm << " " << setw(12) << a[i];
		fprm << setprecision(3) << " " << setw(12) << a[0]*d[0] + a[1]*d[1]
			<< " " << setw(12) << a[2]*d[0] + a[3]*d[1] << endl;
/*		fprm << setprecision(7) << setw(12) << a[0]
			<< setw(12) << a[1]
			<< setw(12) << a[2]
			<< setw(12) << a[3]
			<< setprecision(3) << setw(12) << a[0]*d[0] + a[1]*d[1]
			<< setw(12) << a[2]*d[0] + a[3]*d[1] << endl;*/
	}
	
	fprm.close();
	
	return 0;
}

int			write_imod_tlt(Bstring& imodfile, Bproject* project)
{
	Bstring				paramfile = imodfile.base() + ".tlt";
	Bmicrograph*		mg = project->field->mg;
	
	ofstream			fprm(paramfile.c_str());
	if ( fprm.fail() ) {
		error_show("write_imod_tlt", __FILE__, __LINE__);
		return -1;
	}

	if ( verbose & VERB_PROCESS )
		cout << "Writing tilt angles to " << paramfile << endl;
	
	for ( mg = project->field->mg; mg; mg = mg->next )
		fprm << mg->tilt_angle*180.0/M_PI << endl;
	
	fprm.close();

	return 0;
}

int			write_imod_xyz(Bstring& imodfile, Bproject* project)
{
	Bstring				paramfile = imodfile.base() + "fid.xyz";
//	Bmicrograph*		mg = project->field->mg;
	Breconstruction*	rec = project->rec;
	Bmarker*			mark = NULL;
	
	ofstream			fprm(paramfile.c_str());
	if ( fprm.fail() ) {
		error_show("write_project_imod", __FILE__, __LINE__);
		return -1;
	}

	if ( verbose & VERB_PROCESS )
		cout << "Writing the 3D marker model to " << paramfile << endl;
	
	for ( mark = rec->mark; mark; mark = mark->next ) {
		fprm << setw(4) << mark->id
			<< setprecision(2) << setw(10) << mark->loc[0] + rec->origin[0]
			<< setw(10) << mark->loc[1] + rec->origin[1]
			<< setw(10) << mark->loc[2]
			<< setw(7) << 1 << setw(5) << mark->id;
		if ( mark == rec->mark )
			fprm << " Pix: " << setprecision(5) << setw(11) << rec->voxel_size[0]
				<< " Dim: " << setprecision(0) << setw(6) << 2*rec->origin[0]
				<< setw(6) << 2*rec->origin[1];
		fprm << endl;
	}
	
	fprm.close();

	return 0;
}
	

	

/**
@brief 	Creates IMOD parameter files from a project structure.
@param 	imodfile		IMOD parameter file name.
@param 	project			project structure.
@return int				0.

	Writes both .xf, .tlt and .xyz files.
	Calculations:
		Tilt angle: from micrograph tilt angle
		Matrix: A = {cos(ta), -sin(ta), sin(ta), cos(ta)}
		Shift:  d = o - on
	where
		A:  2x2 transformation matrix
		ta: tilt axis angle
		on: nominal micrograph origin (center of mg)
		o:  aligned micrograph origin
		d:  shift for micrograph

**/
int			write_project_imod(Bstring& imodfile, Bproject* project)
{
	int				err(0);
	
	err += write_imod_xf(imodfile, project);
	
	err += write_imod_tlt(imodfile, project);
	
	if ( project->rec && project->rec->mark )
		err += write_imod_xyz(imodfile, project);
	
//	err += write_imod_fid(imodfile, project);
	
	return err;
}
