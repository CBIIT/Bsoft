/**
@file	rwmgRELION.cpp
@brief	Library routines to read and write micrograph parameters in RELION STAR format
@author Bernard Heymann
@date	Created: 20061101
@date	Modified: 20210330
**/

#include "mg_processing.h"
#include "star.h"
#include "mg_tags.h"
#include "rwimg.h"
#include "linked_list.h"
#include "utilities.h"

#include <sys/stat.h>
#include <fcntl.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			relion_to_project(Bstar2& star, Bproject* project)
{
	int				err(0);
	long			i, j, pid(0);
	Bstring			pfile;
	double			px(1), mag(1), defU, defV, prev_defU(0), volt(12e4), Cs(2e7), amp(0.07);

	Bstring			field_id("1"), mg_id("1");
	Bfield*			field = field_add(&project->field, field_id);
	Bmicrograph*	mg = NULL;
	Bparticle*		partlist = NULL;
	Bparticle*		part = NULL;
	Bparticle*		prev_part = NULL;

	Bstring*		mg_name_list = NULL;
	Bstring*		mg_name = NULL;

	for ( auto ib: star.blocks() ) {
		for ( auto il: ib.loops() ) {
			if ( ( i = il.find("rlnMicrographName") ) >= 0 ) {
				for ( auto ir: il.data() ) {
					mg_name = string_add(&mg_name, ir[i].c_str());
					if ( !mg_name_list ) mg_name_list = mg_name;
					part = particle_add(&part, ++pid);
					if ( !partlist ) partlist = part;
					if ( ( j = il.find("rlnImageName") ) >= 0 ) {
						pfile = ir[j];
						part->id = pfile.pre('@').integer();
						pfile = pfile.post('@');
					}
					if ( ( j = il.find("rlnCoordinateX") ) >= 0 )
						part->loc[0] = to_real(ir[j]);
					if ( ( j = il.find("rlnCoordinateY") ) >= 0 )
						part->loc[1] = to_real(ir[j]);
					if ( ( j = il.find("rlnDetectorPixelSize") ) >= 0 )
						px = to_real(ir[j]);
					if ( ( j = il.find("rlnMagnification") ) >= 0 )
						mag = to_real(ir[j]);
					if ( ( j = il.find("rlnDefocusU") ) >= 0 )
						part->def = to_real(ir[j]);
					if ( ( j = il.find("rlnDefocusV") ) >= 0 )
						part->dev = to_real(ir[j]);
					if ( ( j = il.find("rlnDefocusAngle") ) >= 0 )
						part->ast = to_real(ir[j]);
					if ( ( j = il.find("rlnVoltage") ) >= 0 )
						volt = 1e3 * to_real(ir[j]);
					if ( ( j = il.find("rlnSphericalAberration") ) >= 0 )
						Cs = 1e7 * to_real(ir[j]);
					if ( ( j = il.find("rlnAmplitudeContrast") ) >= 0 )
						amp = to_real(ir[j]);
				}
			}
		}
	}

	if ( verbose )
		cout << "Particles:                      " << pid-1 << endl;
	
	// Set up the micrographs
	for ( i=1, part = partlist, mg_name=mg_name_list; part; part = part->next, ++i ) {
		defU = part->def;
		defV = part->dev;
		part->def = (defU + defV)/2;
		part->dev = fabs(defU - defV)/2;
		if ( defU != prev_defU ) {
			mg_id = Bstring(i, "%04d");
			mg = micrograph_add(&mg, mg_id);
			if ( !field->mg ) field->mg = mg;
			if ( mg_name ) mg->fmg = *mg_name;
			if ( pfile.length() ) mg->fpart = pfile;
			mg->sampling = px*1e4;		// Detector pxel size in micron
			if ( mag ) {
				mg->magnification = mag;
				mg->pixel_size[0] = mg->pixel_size[1] = 1e4*px/mag;
			}
			mg->ctf = new CTFparam;
			mg->ctf->volt(volt);
			mg->ctf->Cs(Cs);
			mg->ctf->amp_shift(asin(amp));
			mg->ctf->defocus_average(part->def);
			mg->ctf->defocus_deviation(part->dev);
			mg->ctf->astigmatism_angle(part->ast);
			mg->ctf->zero(1);
			mg->part = part;
			if ( prev_part ) prev_part->next = NULL;
			pid = 1;
		}
		part->pixel_size = mg->pixel_size;
		prev_defU = defU;
		prev_part = part;
		if ( mg_name ) mg_name = mg_name->next;
	}
	
	if ( verbose )
		cout << "Micrographs:                    " << mg_id.integer() << endl;
	
	// Reconfigure the particle files
	Bimage*			p = NULL;

	if ( pfile.length() ) p = read_img(pfile, 0, 0);
	
	for ( mg = field->mg; mg; mg = mg->next ) {
		for ( pid=1, part = mg->part; part; part = part->next, pid++ ) {
			part->id = pid;
			if ( p ) {
				part->ori = p->size()/2;
				part->pixel_size = p->sampling(0);
			}
			if ( part->mag < 0.5 ) part->mag = 1;
		}
		if ( p ) {
			mg->fpart = pfile;
			mg->pixel_size = p->sampling(0);
			mg->box_size = p->size();
		}
	}
	
	if ( p ) delete p;
	
	return err;
}

/*
loop_ 
_rlnMicrographNameNoDW #1 
_rlnMicrographName #2 
_rlnCtfImage #3 
_rlnDefocusU #4 
_rlnDefocusV #5 
_rlnDefocusAngle #6 
_rlnVoltage #7 
_rlnSphericalAberration #8 
_rlnAmplitudeContrast #9 
_rlnMagnification #10 
_rlnDetectorPixelSize #11 
_rlnCtfFigureOfMerit #12 
_rlnCtfMaxResolution #13 
MotionCorr/job161/Micrographs/20171002_1000_A012_G000_H1000_D001_noDW.mrc MotionCorr/job161/Micrographs/20171002_1000_A012_G000_H1000_D001.mrc CtfFind/job163/Micrographs/20171002_1000_A012_G000_H1000_D001_noDW.ctf:mrc 23997.070312 23977.537109    59.369221   300.000000     2.700000     0.070000 10000.000000     0.660000    -0.004115     7.167600 

loop_ 
_rlnCoordinateX #1 
_rlnCoordinateY #2 
_rlnClassNumber #3 
_rlnAutopickFigureOfMerit #4 
_rlnAnglePsi #5 
_rlnImageName #6 
_rlnMicrographName #7 
_rlnVoltage #8 
_rlnDefocusU #9 
_rlnDefocusV #10 
_rlnDefocusAngle #11 
_rlnSphericalAberration #12 
_rlnCtfBfactor #13 
_rlnCtfScalefactor #14 
_rlnPhaseShift #15 
_rlnAmplitudeContrast #16 
_rlnMagnification #17 
_rlnDetectorPixelSize #18 
_rlnCtfMaxResolution #19 
_rlnCtfFigureOfMerit #20 
 6975.412944   852.888919            4     1.220089   165.000000 000001@Extract/job174/Micrographs/20171002_1000_A012_G000_H1000_D001.mrcs MotionCorr/job161/Micrographs/20171002_1000_A012_G000_H1000_D001.mrc   300.000000 23997.070312 23977.537109    59.369221     2.700000     0.000000     1.000000     0.000000     0.070000 10000.000000     1.980000     7.167600    -0.004115 

*/
int 		project_to_relion(Bproject* project, Bstar2& star, int mg_select, int rec_select)
{
	int				err(0);
	long 			nfield(0), nmg(0);
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	BstarBlock&		block = star.add_block("1");
	BstarLoop&		loop = block.add_loop();
	loop.tags()["rlnCoordinateX"] = 0;
	loop.tags()["rlnCoordinateY"] = 1;

	for ( field=project->field; field; field=field->next ) if ( mg_select < 1 || field->select > 0 ) {
		for ( mg=field->mg; mg; mg=mg->next ) if ( mg_select < 1 || mg->select > 0 ) {
			if ( verbose & VERB_FULL )
				cout << "Writing field \"" << field->id << "\", micrograph \"" << mg->id << "\"" << endl;
			for ( part = mg->part; part; part = part->next ) {
				vector<string>&	vs = loop.add_row(2);
				vs[0] = to_string(part->loc[0]);
				vs[1] = to_string(part->loc[1]);
			}
			nmg++;
		}
		nfield++;
	}

	return err;
}


/**
@brief 	Reading micrograph parameters from a Relion file.
@param 	&filename		file name (or comma-delimited list).
@param 	*project		initialized project structure.
@return int				error code (<0 means failure).
**/
int			read_project_relion(Bstring& filename, Bproject* project)
{
 	Bstar2		star;
	star.line_length(200);                // Set the output line length
	
 	if ( star.read(filename.str()) < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return -1;
	}
	
	if ( verbose ) cout << "Reading a Relion file:          " << filename << endl;

	if ( project->comment.length() < 1 )
		project->comment = star.comment();

	return relion_to_project(star, project);
}

/**
@brief 	Writing micrograph parameters to a Relion file.
@param 	*filename		file name.
@param 	*project		project structure.
@param 	mg_select		flag to select micrograph.
@param 	rec_select		flag to select reconstruction.
@return int				error code (<0 means failure).
**/
int			write_project_relion(Bstring& filename, Bproject* project, int mg_select, int rec_select)
{
 	Bstar2		star;
	star.line_length(200);                // Set the output line length

	if ( verbose ) cout << "Writing a Relion file:          " << filename << endl;

	int			err = project_to_relion(project, star, mg_select, rec_select);
	
	if ( err ) return err;
	
	return star.write(filename.str());
}

/**
@brief 	Modifies some micrograph parameters.
@param 	*project		project structure.
@param 	partfile		Relion stack of particles file.
@param 	path			path to write new particle files.
@param 	partext			extension to set the new particle image format.
@return int				error code (<0 means failure).
**/
int			project_split_particles(Bproject* project, Bstring partfile, Bstring path, Bstring partext)
{
	if ( path.length() && path[-1] != '/' ) path += "/";
	if ( partext.length() < 1 ) partext = "pif";

	if ( partfile.length() < 3 ) partfile = project->field->mg->fpart;
	if ( partfile.length() < 3 ) partfile = project->field->mg->part->fpart;
	if ( partfile.length() < 3 ) {
		cerr << "Error: No particle image file specified" << endl;
		bexit(-1);
	}
	
	// Directory for individual particle files
	if ( path.length() ) {
		mkdir(path.c_str(), O_CREAT );
		chmod(path.c_str(), 0755);
	}
	
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bparticle*			part = NULL;
	
	if ( verbose )
		cout << "Splitting " << partfile << " into micrograph-associated particle files" << endl << endl;
	
	Bimage*				img_stack = read_img(partfile, 1, -1);
	Bimage*				img_part = NULL;

	size_t				n(0), m;
	Bstring				partname;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			partname = mg->fmg.post_rev('/').pre_rev('.') + "." + partext;
			if ( path.length() ) partname = path + partname;
			for ( m=0, part = mg->part; part; part = part->next, m++ )
				part->ori = img_stack->size()/2;
			if ( verbose )
				cout << mg->id << tab << partname << tab << m << endl;
			if ( m ) {
//				cout << n << tab << m << endl;
				img_part = img_stack->extract(n, n+m-1);
				img_part->sampling(mg->pixel_size);
				img_part->origin(img_part->size()/2);
				write_img(partname, img_part, 0);
				delete img_part;
				mg->fpart = partname;
				mg->box_size = img_stack->size();
				n += m;
			}
		}
	}
	
	delete img_stack;

	
	return 0;
}




