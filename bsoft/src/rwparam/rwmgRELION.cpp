/**
@file	rwmgRELION.cpp
@brief	Library routines to read and write micrograph parameters in RELION STAR format
@author Bernard Heymann
@date	Created: 20061101
@date	Modified: 20220323
**/

#include "mg_processing.h"
#include "star.h"
#include "mg_tags.h"
#include "rwimg.h"
#include "string_util.h"
#include "utilities.h"

#include <sys/stat.h>
#include <fcntl.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			relion_to_project(Bstar& star, Bproject* project)
{
	int					err(0);
	long				i, j, pid(0), box_size(0);
	long				ogrp(0), nogrp(0);
	long				nmg(0), npart(0);
	Bstring				pfile;
	double				mpx(1), ppx(1), defU(0), defV(0);
	Vector3<double>		shift, shift_ang;
	Euler				euler;
	vector<double>		magx, magy, v;
	CTFparam			cp;
	vector<CTFparam>	cpa;

	Bstring				field_id("1"), mg_id("1");
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bparticle*			part = NULL;
	string				mg_name, pmg_name;

	for ( auto ib: star.blocks() ) {
		if ( verbose )
			cout << "Reading block: " << ib.tag() << endl;
		for ( auto il: ib.loops() ) {
			if ( verbose )
				cout << "Loop data size: " << il.data().size() << endl;
			if ( ( i = il.find("rlnOpticsGroupName") ) >= 0 ) {
				for ( auto ir: il.data() ) {
					nogrp++;
					if ( ( j = il.find("rlnOpticsGroupName") ) >= 0 )
						cp.identifier(ir[j]);
					else
						cp.identifier(to_string(nogrp));
					if ( verbose & VERB_DEBUG )
						cout << "DEBUG relion_to_project: optics_group: " << cp.identifier() << endl;
					if ( ( j = il.find("rlnOpticsGroup") ) >= 0 )
						cp.select(to_integer(ir[j]));
					if ( ( j = il.find("rlnMicrographOriginalPixelSize") ) >= 0 )
						mpx = to_real(ir[j]);
					if ( ( j = il.find("rlnImagePixelSize") ) >= 0 )
						ppx = to_real(ir[j]);
					if ( ( j = il.find("rlnImageSize") ) >= 0 )
						box_size = to_integer(ir[j]);
					if ( ( j = il.find("rlnVoltage") ) >= 0 )
						cp.volt(1e3 * to_real(ir[j]));
					if ( ( j = il.find("rlnSphericalAberration") ) >= 0 )
						cp.Cs(1e7 * to_real(ir[j]));
					if ( ( j = il.find("rlnAmplitudeContrast") ) >= 0 )
						cp.amp_shift(asin(to_real(ir[j])));	// Not sure if this is correct
					if ( ( j = il.find("rlnBeamTiltX") ) >= 0 )
						cp.beam_tiltX(to_real(ir[j]));
					if ( ( j = il.find("rlnBeamTiltY") ) >= 0 )
						cp.beam_tiltY(to_real(ir[j]));
					if ( ( j = il.find("rlnMagMat00") ) >= 0 )
						magx.push_back(to_real(ir[j]));
					if ( ( j = il.find("rlnMagMat11") ) >= 0 )
						magy.push_back(to_real(ir[j]));
					if ( ( j = il.find("rlnOddZernike") ) >= 0 ) {
						v = parse_real_vector(ir[j].substr(1));
						cp.add_zernike_odd(v);
					}
					if ( ( j = il.find("rlnEvenZernike") ) >= 0 ) {
						v = parse_real_vector(ir[j].substr(1));
						cp.add_zernike_even(v);
//						for ( auto d: cp.aberration_even() ) cout << tab << d;
//						cout << endl;
					}
					cpa.push_back(cp);
					if ( verbose & VERB_PROCESS ) {
						cout << "Optics group:           " << nogrp << endl;
						cp.show();
					}
				}
			}
			if ( ( i = il.find("rlnMicrographName") ) >= 0 ) {
				for ( auto ir: il.data() ) {
					mg_name = ir[i];
					if ( verbose & VERB_DEBUG )
						cout << "DEBUG relion_to_project: micrograph: " << mg_name << endl;
					if ( ( j = il.find("rlnOpticsGroup") ) >= 0 )
						ogrp = to_integer(ir[j]);
					if ( mg_name.compare(pmg_name) ) {
						if ( verbose & VERB_FULL )
							cout << "Micrograph: " << mg_name << endl;
						field_id = mg_id = Bstring(++nmg, "%d");
						pmg_name = mg_name;
						field = field_add(&field, field_id);
						if ( !project->field ) project->field = field;
						mg = micrograph_add(&field->mg, mg_id);
						mg->fmg = mg_name;
						mg->pixel_size[0] = mg->pixel_size[1] = mpx;
						mg->box_size[0] = mg->box_size[1] = box_size;
						mg->box_size[2] = 1;
						mg->ctf = new CTFparam;
						if ( nogrp && ogrp <= nogrp ) {
							mg->ctf->update(cpa[ogrp-1]);
							if ( verbose & VERB_FULL )
								cout << nmg << tab << ogrp << endl;
						}
						part = NULL;
						pid = 0;
					}
					part = particle_add(&part, ++pid);
					npart++;
					if ( !mg->part ) mg->part = part;
					if ( ( j = il.find("rlnImageName") ) >= 0 ) {
						pfile = ir[j];
						part->id = pfile.pre('@').integer();
						mg->fpart = pfile.post('@');
					}
					if ( ( j = il.find("rlnDetectorPixelSize") ) >= 0 )
						ppx = to_real(ir[j]);
					if ( ( j = il.find("rlnCoordinateX") ) >= 0 )
						part->loc[0] = to_real(ir[j]);
					if ( ( j = il.find("rlnCoordinateY") ) >= 0 )
						part->loc[1] = to_real(ir[j]);
					if ( ( j = il.find("rlnOriginX") ) >= 0 )
						shift[0] = to_real(ir[j]);
					if ( ( j = il.find("rlnOriginY") ) >= 0 )
						shift[1] = to_real(ir[j]);
					if ( ( j = il.find("rlnOriginXAngst") ) >= 0 )
						shift_ang[0] = to_real(ir[j]);
					if ( ( j = il.find("rlnOriginYAngst") ) >= 0 )
						shift_ang[1] = to_real(ir[j]);
					if ( ( j = il.find("rlnAnglePsi") ) >= 0 )
						euler[0] = to_real(ir[j])*M_PI/180.0;
					if ( ( j = il.find("rlnAngleTilt") ) >= 0 )
						euler[1] = to_real(ir[j])*M_PI/180.0;
					if ( ( j = il.find("rlnAngleRot") ) >= 0 )
						euler[2] = to_real(ir[j])*M_PI/180.0;	// Correct euler conversion: tested 20211224
					if ( ( j = il.find("rlnMagnification") ) >= 0 )
						part->mag = to_real(ir[j]);		// Deprecated in favor of anisotropic pixel size
					if ( ( j = il.find("rlnDefocusU") ) >= 0 )
						defU = to_real(ir[j]);			// Minimum defocus in angstrom
					if ( ( j = il.find("rlnDefocusV") ) >= 0 )
						defV = to_real(ir[j]);			// Maximum defocus in angstrom
					if ( ( j = il.find("rlnDefocusAngle") ) >= 0 )
						part->ast = to_real(ir[j])*M_PI/180.0;	// Angle could be off by π/2
					if ( ( j = il.find("rlnVoltage") ) >= 0 )
						mg->ctf->volt(1e3 * to_real(ir[j]));
					if ( ( j = il.find("rlnSphericalAberration") ) >= 0 )
						mg->ctf->Cs(1e7 * to_real(ir[j]));
					if ( ( j = il.find("rlnAmplitudeContrast") ) >= 0 )
						mg->ctf->amp_shift(asin(to_real(ir[j])));	// Not sure if this is correct
					if ( ( j = il.find("rlnMaxValueProbDistribution") ) >= 0 )
						part->fom[0] = to_real(ir[j]);		// I assume this is a good FOM
					if ( ( j = il.find("rlnClassNumber") ) >= 0 )
						part->sel = to_integer(ir[j]);
					if ( ( j = il.find("rlnGroupNumber") ) >= 0 )
						part->group = to_integer(ir[j]);
					if ( ( j = il.find("rlnCtfFigureOfMerit") ) >= 0 )
						mg->ctf->fom(to_real(ir[j]));
					if ( nogrp && ogrp <= nogrp ) {
						if ( magx.size() == nogrp ) {
							part->pixel_size[0] = ppx / magx[ogrp-1];	// Correct direction: tested 20211229
							part->pixel_size[1] = ppx / magy[ogrp-1];
						} else {
							part->pixel_size[0] = part->pixel_size[1] = ppx;
						}
					} else {
						part->pixel_size[0] = part->pixel_size[1] = ppx;
					}
					if ( ( j = il.find("rlnOriginX") ) >= 0 )
						part->ori = -shift + mg->box_size/2;
					else
						part->ori = -shift_ang/ppx + mg->box_size/2;
					part->view = euler.view();
					part->def = (defU + defV)/2;
					part->dev = fabs(defU - defV)/2;
					part->ast = angle_set_negPI_to_PI(part->ast);
					mg->ctf->defocus_average(part->def);
					mg->ctf->astigmatism(part->dev, part->ast);
					mg->ctf->zero(1);
				}
			}
		}
	}

	if ( verbose ) {
		cout << "Optics groups:                  " << nogrp << endl;
		cout << "Micrographs:                    " << nmg << endl;
		cout << "Particles:                      " << npart << endl;
	}

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
int 		project_to_relion(Bproject* project, Bstar& star, int mg_select, int rec_select)
{
	int				err(0);
	bool			found(0);
	long 			nt(0), i(0), nfield(0), nmg(0), npart(0), ogrp(0), nogrp(0);
	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bparticle*		part = mg->part;
	CTFparam		cp = *mg->ctf;
	vector<string>	cpid;

	BstarBlock&		block = star.add_block("optics");
	if ( verbose )
		cout << "Writing block: " << block.tag() << endl;
	BstarLoop&		loop = block.add_loop();
	loop.tags()["rlnOpticsGroupName"] = nt++;
	loop.tags()["rlnOpticsGroup"] = nt++;
//	loop.tags()["rlnMtfFileName"] = nt++;
	loop.tags()["rlnMicrographOriginalPixelSize"] = nt++;
	loop.tags()["rlnVoltage"] = nt++;
	loop.tags()["rlnSphericalAberration"] = nt++;
	loop.tags()["rlnAmplitudeContrast"] = nt++;
	loop.tags()["rlnImagePixelSize"] = nt++;
	loop.tags()["rlnImageSize"] = nt++;
	loop.tags()["rlnImageDimensionality"] = nt++;
	loop.tags()["rlnBeamTiltX"] = nt++;
	loop.tags()["rlnBeamTiltY"] = nt++;
	loop.tags()["rlnOddZernike"] = nt++;
	loop.tags()["rlnEvenZernike"] = nt++;
//	loop.tags()["rlnCtfDataAreCtfPremultiplied"] = nt++;
//	loop.tags()["rlnMagMat00"] = nt++;
//	loop.tags()["rlnMagMat01"] = nt++;
//	loop.tags()["rlnMagMat10"] = nt++;
//	loop.tags()["rlnMagMat11"] = nt++;
//	cout << "Tags: " << nt << endl;

	for ( field=project->field; field; field=field->next ) if ( mg_select < 1 || field->select > 0 ) {
		for ( mg=field->mg; mg; mg=mg->next ) if ( mg_select < 1 || mg->select > 0 ) {
			cp = *mg->ctf;
//			cp.show();
			found = 0;
			for ( auto it: cpid ) if ( it == cp.identifier() ) found = 1;
			if ( !found ) {
				cpid.push_back(cp.identifier());
				nogrp++;
				i=0;
				vector<string>&	vs = loop.add_row(nt);
				vs[i++] = cp.identifier();
				vs[i++] = to_string(++ogrp);
				vs[i++] = to_string(mg->pixel_size[0]);
				vs[i++] = to_string(1e-3*cp.volt());
				vs[i++] = to_string(1e-7*cp.Cs());
				vs[i++] = to_string(sin(cp.amp_shift()));
				vs[i++] = to_string(part->pixel_size[0]);
				vs[i++] = to_string(mg->box_size[0]);
				vs[i++] = "2";
				vs[i++] = to_string(cp.beam_tiltX());
				vs[i++] = to_string(cp.beam_tiltY());
				vs[i++] = "[" + concatenate(cp.zernike_odd()) + "]";
				vs[i++] = "[" + concatenate(cp.aberration_even_difference()) + "]";
				//	cout << "Values: " << i << endl;
				if ( i != nt )
					cerr << "Error: The number of values (" << i << " must equal the number of tags (" << nt << ")!" << endl;
			}
		}
	}
	
//	cout << "Optics groups:" << endl;
//	for ( auto it: cpid ) cout << it << endl;
	
	BstarBlock&		block2 = star.add_block("particles");
	if ( verbose )
		cout << "Writing block: " << block2.tag() << endl;
	BstarLoop&		loop2 = block2.add_loop();
	nt = 0;
	loop2.tags()["rlnCoordinateX"] = nt++;
	loop2.tags()["rlnCoordinateY"] = nt++;
//	loop2.tags()["rlnAutopickFigureOfMerit"] = nt++;
	loop2.tags()["rlnClassNumber"] = nt++;
	loop2.tags()["rlnAnglePsi"] = nt++;
	loop2.tags()["rlnImageName"] = nt++;
	loop2.tags()["rlnMicrographName"] = nt++;
	loop2.tags()["rlnOpticsGroup"] = nt++;
//	loop2.tags()["rlnCtfMaxResolution"] = nt++;
	loop2.tags()["rlnCtfFigureOfMerit"] = nt++;
	loop2.tags()["rlnDefocusU"] = nt++;
	loop2.tags()["rlnDefocusV"] = nt++;
	loop2.tags()["rlnDefocusAngle"] = nt++;
//	loop2.tags()["rlnCtfBfactor"] = nt++;
//	loop2.tags()["rlnCtfScalefactor"] = nt++;
//	loop2.tags()["rlnPhaseShift"] = nt++;
	loop2.tags()["rlnAngleRot"] = nt++;
	loop2.tags()["rlnAngleTilt"] = nt++;
	loop2.tags()["rlnOriginXAngst"] = nt++;
	loop2.tags()["rlnOriginYAngst"] = nt++;
//	loop2.tags()["rlnNormCorrection"] = nt++;
//	loop2.tags()["rlnLogLikeliContribution"] = nt++;
	loop2.tags()["rlnMaxValueProbDistribution"] = nt++;
//	loop2.tags()["rlnNrOfSignificantSamples"] = nt++;
	loop2.tags()["rlnGroupNumber"] = nt++;
//	loop2.tags()["rlnRandomSubset"] = nt++;
//	cout << "Tags: " << nt << endl;

	for ( field=project->field; field; field=field->next ) if ( mg_select < 1 || field->select > 0 ) {
		for ( mg=field->mg; mg; mg=mg->next ) if ( mg_select < 1 || mg->select > 0 ) {
			if ( verbose & VERB_FULL )
				cout << "Writing field \"" << field->id << "\", micrograph \"" << mg->id << "\"" << endl;
			cp = *mg->ctf;
			for ( i=0; i<cpid.size(); ++i ) if ( cpid[i] == cp.identifier() ) ogrp = i+1;
			for ( part = mg->part; part; part = part->next ) {
				Euler	euler(part->view);
				vector<string>&	vs = loop2.add_row(nt);
				i = 0;
				vs[i++] = to_string(part->loc[0]);
				vs[i++] = to_string(part->loc[1]);
				vs[i++] = to_string(part->sel);
				vs[i++] = to_string(euler.psi()*180.0/M_PI);
				vs[i++] = mg->fpart.str();
				vs[i++] = mg->fmg.str();
				vs[i++] = to_string(cp.select());
				vs[i++] = to_string(cp.fom());
				vs[i++] = to_string(part->def - part->dev);
				vs[i++] = to_string(part->def + part->dev);
				vs[i++] = to_string(part->ast*180.0/M_PI);
				vs[i++] = to_string(euler.phi()*180.0/M_PI);
				vs[i++] = to_string(euler.theta()*180.0/M_PI);
				vs[i++] = to_string(part->ori[0]*part->pixel_size[0]);
				vs[i++] = to_string(part->ori[1]*part->pixel_size[1]);
				vs[i++] = to_string(part->fom[0]);
				vs[i++] = to_string(part->group);
				npart++;
//				cout << "Values: " << i << endl;
				if ( i != nt )
					cerr << "Error: The number of values (" << i << " must equal the number of tags (" << nt << ")!" << endl;
			}
			nmg++;
		}
		nfield++;
	}
	
	if ( verbose ) {
		cout << "Optics groups:                  " << nogrp << endl;
		cout << "Micrographs:                    " << nmg << endl;
		cout << "Particles:                      " << npart << endl;
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
 	Bstar		star;
	star.line_length(200);                // Set the output line length
	
 	if ( star.read(filename.str()) < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return -1;
	}
	
	if ( verbose ) cout << "Reading a Relion file:          " << filename << endl;

	project->comment = "# Relion file: " + filename + "\n";

	if ( star.comment().size() )
		project->comment += star.comment();

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
 	Bstar		star;
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




