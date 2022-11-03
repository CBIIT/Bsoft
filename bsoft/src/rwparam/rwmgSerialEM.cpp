/**
@file	rwmgSerialEM.cpp
@brief	Converts between a SerialEM MDOC file and a micrograph parameter file
@author	Bernard Heymann
@date	Created: 20190109
@date	Modified: 20210728
**/

#include "rwmgSerialEM.h"
#include "mg_img_proc.h"
#include "mg_tomography.h"
#include "rwmg.h"
#include "rwimg.h"
#include "mdoc.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Creates a project structure using SerialEM MDOC parameters.
@param 	filename		SerialEM MDOC file name.
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
int			read_project_serialem(Bstring& filename, Bproject* project, int flag)
{
	if ( filename.length() < 1 ) {
		cerr << "No filenames given to create a project!" << endl;
		exit(-1);
	}
	
	int					err(0);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_project_serialem: reading " << filename << endl;
	
	MDOCparser			parser(filename.str());
	JSvalue				root = parser.parse();

	if ( verbose ) {
		cout << "Reading " << filename << endl;
		cout << "Number of sections: " << root["Sections"].size() << endl;
	}

//	cout << root << endl;
//	cout << root["ImageFile"].value() << endl;
	
	Bstring				path, base, file, imgfile, text;
	long				img_num;
	double				axis_angle(0);
	
	if ( filename.contains("/") ) path = filename.pre_rev('/');
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_project_serialem: path = " << path << endl;

	imgfile = root["ImageFile"].value();
	if ( imgfile.contains("/") ) imgfile = imgfile.post_rev('/');
	if ( imgfile.contains("\\") ) imgfile = imgfile.post_rev('\\');
	
	Bstring				ext = imgfile.extension();
	if ( ext.length() < 1 ) imgfile += ".st";
	base = imgfile.base();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_project_serialem: file = " << imgfile << endl;

	if ( access(imgfile.c_str(), R_OK) )
		imgfile = path + "/" + imgfile;

	double				volt(0), Cs(2.0), amp(0.07), sc(1);

	if ( root.exists("Voltage") ) volt = root["Voltage"].real();

	JSvalue&			sections = root["Sections"];

	Bfield*				field = field_add(&project->field, base);
	Bstring				mg_id;
	Bmicrograph*		mg = NULL;
	Bimage* 			p = NULL;
	
	for ( auto it = sections.begin(); it != sections.end(); ++it ) {
		JSvalue&		obj = *it;
//		cout << obj << endl;
		if ( obj.exists("ZValue") ) {
			img_num = obj["ZValue"].integer();
			if ( obj.exists("SubFramePath") ) {
				imgfile = obj["SubFramePath"].value();
				if ( imgfile.contains("/") ) imgfile = imgfile.post_rev('/');
				if ( imgfile.contains("\\") ) imgfile = imgfile.post_rev('\\');
				mg_id = imgfile.base();
				if ( access(imgfile.c_str(), R_OK) )
					imgfile = path + "/" + imgfile;
				if ( access(imgfile.c_str(), R_OK) == 0 )
					p = read_img(imgfile, 0, -1);
			} else {
				mg_id = base + Bstring(img_num+1, "_%03d");
				if ( access(imgfile.c_str(), R_OK) == 0 )
					p = read_img(imgfile, 0, img_num);
			}
			if ( !p ) {
				cerr << "Error: Image file " << imgfile << " not read!" << endl;
				error_show("read_project_serialem", __FILE__, __LINE__);
				return -1;
			}
			if ( p->sizeZ() > 1 ) p->slices_to_images();
			mg = micrograph_add(&mg, mg_id);
			if ( !field->mg ) field->mg = mg;
//			cout << "//" << obj["ZValue"] << "//" << endl;
			mg->block = img_num;
			if ( obj.exists("SubFramePath") ) {
				mg->img_num = 0;
				mg->fframe = imgfile;
			} else {
				mg->img_num = img_num;
				mg->fmg = imgfile;
			}
			mg->pixel_size = p->sampling(0);
			if ( obj.exists("PixelSpacing") ) {
				mg->pixel_size[0] = mg->pixel_size[1] = obj["PixelSpacing"].real();
				mg->frame_pixel_size[0] = mg->frame_pixel_size[1] = obj["PixelSpacing"].real();
			}
			mg->origin = p->size()/2;
			if ( obj.exists("TiltAngle") ) mg->tilt_angle = obj["TiltAngle"].real()*M_PI/180.0;
			mg->tilt_axis = axis_angle;
			mg->mark_radius = 10;
			mg->intensity = p->image->average();
			if ( obj.exists("ExposureTime") )
				mg->exposure = obj["ExposureTime"].real();
			if ( obj.exists("ExposureDose") )
				mg->dose = obj["ExposureDose"].real();
			if ( obj.exists("MinMaxMean") ) {
				JSvalue&	mmm = obj["MinMaxMean"];
				if ( obj.exists("CountsPerElectron") )
					sc = obj["CountsPerElectron"].real();
				mg->intensity = mmm[2].real()/sc;
			}
			if ( obj.exists("Magnification") )
				mg->magnification = obj["Magnification"].real();
			mg->ctf = new CTFparam(volt, Cs, amp);
			if ( obj.exists("TargetDefocus") )
				mg->ctf->defocus_average(-obj["TargetDefocus"].real());
			if ( verbose & VERB_FULL ) {
				cout << "Creating micrograph " << mg_id << ", with " << p->images() << " frames" << endl;
				cout << "Pixel size:                     " << mg->pixel_size << endl;
				cout << "Origin:                         " << mg->origin << endl;
				cout << "Tilt angle:                     " << mg->tilt_angle*180.0/M_PI << endl;
				cout << "Intensity:                      " << mg->intensity << endl;
//				mg->ctf->show();
			}
			delete p;
		} if ( obj.exists("T") ) {
			text = obj["T"].value();
			if ( text.contains("SerialEM") ) {
				project->comment = "# " + text + "\n";
			} else if ( text.contains("Tilt axis angle") ) {
				axis_angle = text.post('=').real();
				axis_angle = (axis_angle - 90)*M_PI/180.0;
			}
		}
	}

	project_sort_by_tilt(project);
	project_mg_tilt_to_matrix(project);
	
	return err;
}


/**
@brief 	Creates SerialEM MDOC files from a project structure.
@param 	filename		SerialEM MDOC file name.
@param 	project			project structure.
@return int				0.

**/
int			write_project_serialem(Bstring& filename, Bproject* project)
{
	int				err(0);
	
//	err += write_imod_xf(semfile, project);
	
	cerr << "Error: Not implemented yet!" << endl;
		
	return err;
}
