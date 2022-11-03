/**
@file	rwmgSTAR.cpp
@brief	Library routines to read and write micrograph parameters in STAR format
@author Bernard Heymann
@date	Created: 20010206
@date	Modified: 20220113
**/

#include "mg_processing.h"
#include "rwmgSTAR.h"
#include "star.h"
#include "mg_tags.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
int			star_to_project(Bstar& star, Bproject* project, int flag);
int 		project_to_star(Bproject* project, Bstar& star, int mg_select, int rec_select);

/**
@brief 	Reading micrograph parameters from STAR files.
@param 	&filename		file name (or comma-delimited list).
@param 	*project		initialized project structure.
@param	flag			update tags.
@return int				error code (<0 means failure).
**/
int			read_project_star(Bstring& filename, Bproject* project, int flag)
{
 	Bstar		star;
	star.line_length(200);                // Set the output line length
	
 	if ( star.read(filename.str()) < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return -1;
	}
	
	if ( project->comment.length() < 1 )
		project->comment = star.comment();
	
	int			err = star_to_project(star, project, flag);
	
	return err;
}

/**
@brief 	Writing micrograph parameters to a STAR file.
@param 	&filename		file name.
@param 	*project		project structure.
@param 	mg_select		flag to only write selected micrographs.
@param 	rec_select		flag to only convert selected reconstructions.
@return int				error code (<0 means failure).
**/
int			write_project_star(const char* filename, Bproject* project, int mg_select, int rec_select)
{
	Bstring		thefile(filename);
	return write_project_star(thefile, project, mg_select, rec_select);
}

int			write_project_star(Bstring& filename, Bproject* project, int mg_select, int rec_select)
{
 	Bstar		star;
	star.line_length(200);                // Set the output line length
	
	int			err = project_to_star(project, star, mg_select, rec_select);
	
	if ( err ) return err;
	
	return star.write(filename.str());
}

/**
@brief 	Replacing old STAR tags with new ones.
@param	&star			STAR database.
@return int				0.
**/
int			mg_star_update_tags(Bstar& star)
{
	star.replace_tag(MICROGRAPH_PARTICLE_FILE, PARTICLE_FILE);
	star.replace_tag(MICROGRAPH_FILAMENT_FILE, FILAMENT_FILE);
	
	star.replace_tag(PARTICLE_NUMBER, PARTICLE_ID);
	star.replace_tag(PARTICLE_MG_X, PARTICLE_X);
	star.replace_tag(PARTICLE_MG_Y, PARTICLE_Y);
	star.replace_tag(PARTICLE_MG_Z, PARTICLE_Z);
	star.replace_tag(PARTICLE_X_ORIGIN, PARTICLE_ORIGIN_X);
	star.replace_tag(PARTICLE_Y_ORIGIN, PARTICLE_ORIGIN_Y);
	star.replace_tag(PARTICLE_Z_ORIGIN, PARTICLE_ORIGIN_Z);
	star.replace_tag(MICROGRAPH_BOX_RADIUS, PARTICLE_BOX_RADIUS);
	star.replace_tag(MICROGRAPH_BOX_RADIUS_X, PARTICLE_BOX_RADIUS_X);
	star.replace_tag(MICROGRAPH_BOX_RADIUS_Y, PARTICLE_BOX_RADIUS_Y);
	star.replace_tag(MICROGRAPH_BOX_RADIUS_Z, PARTICLE_BOX_RADIUS_Z);
	star.replace_tag(MICROGRAPH_BAD_RADIUS, PARTICLE_BAD_RADIUS);
	star.replace_tag(MICROGRAPH_BAD_X, PARTICLE_BAD_X);
	star.replace_tag(MICROGRAPH_BAD_Y, PARTICLE_BAD_Y);
	star.replace_tag(MICROGRAPH_BAD_Z, PARTICLE_BAD_Z);

	star.replace_tag(MICROGRAPH_MARKER_RADIUS, MARKER_RADIUS);
	star.replace_tag(MICROGRAPH_MARKER_ID, MARKER_ID);
	star.replace_tag(MICROGRAPH_MARKER_X, MARKER_X);
	star.replace_tag(MICROGRAPH_MARKER_Y, MARKER_Y);
	star.replace_tag(MICROGRAPH_MARKER_Z, MARKER_Z);
	star.replace_tag(MICROGRAPH_MARKER_ERROR_X, MARKER_ERROR_X);
	star.replace_tag(MICROGRAPH_MARKER_ERROR_Y, MARKER_ERROR_Y);
	star.replace_tag(MICROGRAPH_MARKER_ERROR_Z, MARKER_ERROR_Z);
	star.replace_tag(MICROGRAPH_MARKER_FOM, MARKER_FOM);

	star.replace_tag(MICROGRAPH_VOLTAGE, CTF_VOLTAGE);
	star.replace_tag(MICROGRAPH_CTF_CS, CTF_CS);
	star.replace_tag(MICROGRAPH_CTF_CC, CTF_CC);
	star.replace_tag(MICROGRAPH_CTF_ALPHA, CTF_ALPHA);
	star.replace_tag(MICROGRAPH_CTF_DE, CTF_DE);
	star.replace_tag(MICROGRAPH_CTF_AMP_CONT, CTF_AMP);
	star.replace_tag(MICROGRAPH_CTF_ZERO, CTF_ZERO);
	star.replace_tag(MICROGRAPH_CTF_DEF_AVG, CTF_DEF_AVG);
	star.replace_tag(MICROGRAPH_CTF_DEF_DEV, CTF_DEF_DEV);
	star.replace_tag(MICROGRAPH_CTF_DEF_MIN, CTF_DEF_MIN);
	star.replace_tag(MICROGRAPH_CTF_DEF_MAX, CTF_DEF_MAX);
	star.replace_tag(MICROGRAPH_CTF_AST_ANG, CTF_AST_ANG);
	star.replace_tag(MICROGRAPH_CTF_BASELINE, CTF_BASELINE);
	star.replace_tag(MICROGRAPH_CTF_ENVELOPE, CTF_ENVELOPE);

	return 0;
}

Bframe*		frame_from_starblock(BstarBlock& block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG frame_from_starblock:" << endl;

	long			i, j;
	Bframe*			framelist = NULL;
	Bframe*			frame = NULL;

	for ( auto il: block.loops() ) {
		if ( ( i = il.find(MICROGRAPH_FRAME) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				frame = frame_add(&frame, to_integer(ir[i]));
				if ( !framelist ) framelist = frame;
				if ( ( j = il.find(MICROGRAPH_FRAME_SHIFT_X) ) >= 0 )
					frame->shift[0] = to_real(ir[j]);
				if ( ( j = il.find(MICROGRAPH_FRAME_SHIFT_Y) ) >= 0 )
					frame->shift[1] = to_real(ir[j]);
				if ( ( j = il.find(MICROGRAPH_FRAME_SELECT) ) >= 0 )
					frame->sel = to_integer(ir[j]);
				if ( ( j = il.find(MICROGRAPH_FRAME_FOM) ) >= 0 )
					frame->fom = to_real(ir[j]);
			}
		}
	}

	return framelist;
}

int			project_get_fom_tags(Bstar star, FOMType* fom_tag)
{
	int				f, n(0);
	Bstring			tag;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_get_fom_tags:" << endl;

	for ( f=0; f<NFOM; ++f ) fom_tag[f] = NoFOM;
	
	for ( auto ib: star.blocks() ) {
		for ( auto il: ib.loops() ) {
			if ( il.find(PARTICLE_ID) >= 0 ) {
				for ( f = FOM; f != FOMlast; ++f ) {
					tag = get_fom_tag(FOMType(f));
					if ( il.find(tag.str()) >= 0 )
						fom_tag[n++] = FOMType(f);
					if ( verbose & VERB_DEBUG )
						cout << "DEBUG project_get_fom_tags: tag="  << tag << tab << n << endl;
					if ( n >= NFOM ) return n;
				}
				return n;
			}
		}
	}
	
	return n;
}

Bparticle*	particle_from_starblock(BstarBlock& block, FOMType fom_tag[NFOM])
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG particle_from_starblock:" << endl;

	long			i, j;
	Bparticle*		partlist = NULL;
	Bparticle*		part = NULL;
	Bstring			tag;
	int				f, omega_flag(0);
	Euler			euler;

	for ( auto il: block.loops() ) {
		if ( ( i = il.find(PARTICLE_ID) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				part = particle_add(&part, to_integer(ir[i]));
				if ( !partlist ) partlist = part;
//				cout << part->id << endl;
				if ( ( j = il.find(PARTICLE_FILE) ) >= 0 )
					part->fpart = ir[j];
				if ( ( j = il.find(PARTICLE_GROUP) ) >= 0 )
					part->group = to_integer(ir[j]);
				if ( ( j = il.find(PARTICLE_DEFOCUS) ) >= 0 )
					part->def = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_DEF_DEV) ) >= 0 )
					part->dev = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_AST_ANG) ) >= 0 )
					part->ast = to_real(ir[j])*M_PI/180.0;
				if ( ( j = il.find(PARTICLE_MAGNIF) ) >= 0 )
					part->mag = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_X) ) >= 0 )
					part->loc[0] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_Y) ) >= 0 )
					part->loc[1] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_Z) ) >= 0 )
					part->loc[2] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_PIXEL) ) >= 0 )
					part->pixel_size[0] = part->pixel_size[1] = part->pixel_size[2] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_PIXEL_X) ) >= 0 )
					part->pixel_size[0] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_PIXEL_Y) ) >= 0 )
					part->pixel_size[1] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_PIXEL_Z) ) >= 0 )
					part->pixel_size[2] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_ORIGIN_X) ) >= 0 )
					part->ori[0] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_ORIGIN_Y) ) >= 0 )
					part->ori[1] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_ORIGIN_Z) ) >= 0 )
					part->ori[2] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_VIEW_X) ) >= 0 )
					part->view[0] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_VIEW_Y) ) >= 0 )
					part->view[1] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_VIEW_Z) ) >= 0 )
					part->view[2] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_VIEW_ANGLE) ) >= 0 )
					part->view[3] = to_real(ir[j])*M_PI/180.0;
				else {
					if ( ( j = il.find(PARTICLE_PSI) ) >= 0 ) omega_flag = 0;
					else if ( ( j = il.find(PARTICLE_OMEGA) ) >= 0 ) omega_flag = 1;
					euler[0] = to_real(ir[j])*M_PI/180.0;
					if ( omega_flag ) euler[0] = -euler[0];
					if ( ( j = il.find(PARTICLE_THETA) ) >= 0 )
						euler[1] = to_real(ir[j])*M_PI/180.0;
					if ( ( j = il.find(PARTICLE_PHI) ) >= 0 )
						euler[2] = to_real(ir[j])*M_PI/180.0;
					part->view = euler.view();
				}
				if ( ( j = il.find(PARTICLE_SELECT) ) >= 0 )
					part->sel = to_integer(ir[j]);
//				cout << tab << part->sel << endl;
				for ( f=0; f<NFOM; f++ ) if ( fom_tag[f] ) {
					tag = get_fom_tag(fom_tag[f]);
//					cout << tab << tag << endl;
					if ( ( j = il.find(tag.str()) ) >= 0 )
						part->fom[f] = to_real(ir[j]);
				}
			}
		}
	}

	return partlist;
}

Bfilament*	filament_from_starblock(BstarBlock& block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG filament_from_starblock:" << endl;

	long			i, j;
	Bfilament*		fillist = NULL;
	Bfilament*		fil = NULL;
	Bfilnode*		fnode = NULL;
	int 			f, fp(-1);

	for ( auto il: block.loops() ) {
		if ( ( i = il.find(FILAMENT_ID) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				f = to_integer(ir[i]);
				if ( f != fp ) fil = filament_add(&fil, f);
				if ( !fillist ) fillist = fil;
				fnode = filament_node_add(&fil->node, 0);
				fp = f;
				if ( ( j = il.find(FILAMENT_FILE) ) >= 0 )
					fil->ffil = ir[j];
				if ( ( j = il.find(FILAMENT_NODE_ID) ) >= 0 )
					fnode->id = to_integer(ir[j]);
				if ( ( j = il.find(FILAMENT_NODE_X) ) >= 0 )
					fnode->loc[0] = to_real(ir[j]);
				if ( ( j = il.find(FILAMENT_NODE_Y) ) >= 0 )
					fnode->loc[1] = to_real(ir[j]);
				if ( ( j = il.find(FILAMENT_NODE_Z) ) >= 0 )
					fnode->loc[2] = to_real(ir[j]);
			}
		}
	}

	return fillist;
}

Bbadarea*	badarea_from_starblock(BstarBlock& block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG badarea_from_starblock:" << endl;
	
	long			i, j;
	Bbadarea*		bad = NULL;
	Bbadarea*		badlist = NULL;

	for ( auto il: block.loops() ) {
		if ( ( i = il.find(PARTICLE_BAD_X) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				bad = (Bbadarea *) add_item((char **) &bad, sizeof(Bbadarea));
				if ( !badlist ) badlist = bad;
				bad->loc[0] = to_real(ir[i]);
				if ( ( j = il.find(PARTICLE_BAD_Y) ) >= 0 )
					bad->loc[1] = to_real(ir[j]);
				if ( ( j = il.find(PARTICLE_BAD_Z) ) >= 0 )
					bad->loc[2] = to_real(ir[j]);
			}
		}
	}

	return badlist;
}

Bmarker*	marker_from_starblock(BstarBlock& block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG marker_from_starblock:" << endl;
	
	long			i, j;
	Bmarker*		marklist = NULL;
	Bmarker*		mark = NULL;

	for ( auto il: block.loops() ) {
		if ( ( i = il.find(MARKER_ID) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				mark = (Bmarker *) add_item((char **) &mark, sizeof(Bmarker));
				if ( !marklist ) marklist = mark;
				mark->id = to_integer(ir[i]);
				mark->fom = 1;
				mark->sel = 1;
				if ( ( j = il.find(MARKER_X) ) >= 0 )
					mark->loc[0] = to_real(ir[j]);
				if ( ( j = il.find(MARKER_Y) ) >= 0 )
					mark->loc[1] = to_real(ir[j]);
				if ( ( j = il.find(MARKER_Z) ) >= 0 )
					mark->loc[2] = to_real(ir[j]);
				if ( ( j = il.find(MARKER_ERROR_X) ) >= 0 )
					mark->err[0] = to_real(ir[j]);
				if ( ( j = il.find(MARKER_ERROR_Y) ) >= 0 )
					mark->err[1] = to_real(ir[j]);
				if ( ( j = il.find(MARKER_ERROR_Z) ) >= 0 )
					mark->err[2] = to_real(ir[j]);
				if ( ( j = il.find(MARKER_RESIDUAL) ) >= 0 )
					mark->res = to_real(ir[j]);
				if ( ( j = il.find(MARKER_FOM) ) >= 0 )
					mark->fom = to_real(ir[j]);
				if ( ( j = il.find(MARKER_SELECT) ) >= 0 )
					mark->sel = to_integer(ir[j]);
			}
		}
	}

	return marklist;
}

Bstrucfac*	strucfac_from_starblock(BstarBlock& block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG strucfac_from_starblock:" << endl;
	
	long			i, j;
	Bstrucfac*		sflist = NULL;
	Bstrucfac*		sf = NULL;
	
	for ( auto il: block.loops() ) {
		if ( ( i = il.find(REFLEX_H) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				sf = (Bstrucfac *) add_item((char **) &sf, sizeof(Bstrucfac));
				if ( !sflist ) sflist = sf;
				sf->fom = sf->sel = 1;
				sf->index[0] = to_integer(ir[i]);
				if ( ( j = il.find(REFLEX_K) ) >= 0 )
					sf->index[1] = to_integer(ir[j]);
				if ( ( j = il.find(REFLEX_L) ) >= 0 )
					sf->index[2] = to_integer(ir[j]);
				if ( ( j = il.find(REFLEX_X) ) >= 0 )
					sf->loc[0] = to_real(ir[j]);
				if ( ( j = il.find(REFLEX_Y) ) >= 0 )
					sf->loc[1] = to_real(ir[j]);
				if ( ( j = il.find(REFLEX_Z) ) >= 0 )
					sf->loc[2] = to_real(ir[j]);
				if ( ( j = il.find(REFLEX_AMP) ) >= 0 )
					sf->amp = to_real(ir[j]);
				if ( ( j = il.find(REFLEX_SIGAMP) ) >= 0 )
					sf->sigamp = to_real(ir[j]);
				if ( ( j = il.find(REFLEX_PHI) ) >= 0 )
					sf->phi = to_real(ir[j]);
				if ( ( j = il.find(REFLEX_SIGPHI) ) >= 0 )
					sf->sigphi = to_real(ir[j]);
				if ( ( j = il.find(REFLEX_FOM) ) >= 0 )
					sf->fom = to_real(ir[j]);
				if ( ( j = il.find(REFLEX_STATUS) ) >= 0 )
					sf->sel = to_integer(ir[j]);
			}
		}
	}

	return sflist;
}

Blayerline*	layerline_from_starblock(BstarBlock& block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG layerline_from_starblock:" << endl;
	
	long			i, j;
	Blayerline*		linelist = NULL;
	Blayerline*		line = NULL;

	for ( auto il: block.loops() ) {
		if ( ( i = il.find(LAYERLINE_NUMBER) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				line = (Blayerline *) add_item((char **) &line, sizeof(Blayerline));
				if ( !linelist ) linelist = line;
				line->fom = line->sel = 1;
				line->number = to_integer(ir[i]);
				if ( ( j = il.find(LAYERLINE_ORDER) ) >= 0 )
					line->order = to_integer(ir[j]);
				if ( ( j = il.find(LAYERLINE_DISTANCE) ) >= 0 )
					line->distance = to_real(ir[j]);
				if ( ( j = il.find(LAYERLINE_FREQ) ) >= 0 )
					line->freq = to_real(ir[j]);
				if ( ( j = il.find(LAYERLINE_AMP) ) >= 0 )
					line->amp = to_real(ir[j]);
				if ( ( j = il.find(LAYERLINE_FOM) ) >= 0 )
					line->fom = to_real(ir[j]);
				if ( ( j = il.find(LAYERLINE_SELECT) ) >= 0 )
					line->sel = to_integer(ir[j]);
			}
		}
	}

	return linelist;
}

CTFparam*		ctf_from_starblock(BstarBlock& block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ctf_from_starblock:" << endl;
	
	if ( !block.exists(CTF_VOLTAGE) ) return NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ctf_from_starblock: volt=" << block.real(CTF_VOLTAGE) << endl;

	CTFparam*		ctf = new CTFparam;
	vector<double>	v;
	
	if ( block.exists(CTF_ID) )
		ctf->identifier(block[CTF_ID]);
	if ( block.exists(CTF_SELECT) )
		ctf->select(block.integer(CTF_SELECT));
	if ( block.exists(CTF_FOM) )
		ctf->fom(block.real(CTF_FOM));
	ctf->volt(block.real(CTF_VOLTAGE));
	ctf->focal_length(block.real(CTF_FOCAL));
	if ( block.exists(CTF_APERTURE) )
		ctf->objective_aperture(block.real(CTF_APERTURE));
	if ( block.exists(CTF_AMP_SHIFT) )
		ctf->amp_shift(block.real(CTF_AMP_SHIFT));
	else
		ctf->amp_shift(asin(block.real(CTF_AMP)));
	if ( block.exists(CTF_DEF_AVG) ) {
		ctf->defocus_average(block.real(CTF_DEF_AVG));
//		ctf->defocus_deviation(block.real(CTF_DEF_DEV));
		ctf->astigmatism(block.real(CTF_DEF_DEV), block.real(CTF_AST_ANG) * M_PI/180.0);
	} else {
		double	min = block.real(CTF_DEF_MIN);
		double	max = block.real(CTF_DEF_MAX);
		ctf->defocus_average((max + min)/2);
//		ctf->defocus_deviation((max - min)/2);
		ctf->astigmatism(fabs(max - min)/2, block.real(CTF_AST_ANG) * M_PI/180.0);
	}
//	ctf->astigmatism_angle(block.real(CTF_AST_ANG) * M_PI/180.0);
	ctf->Cs(block.real(CTF_CS));
	ctf->Cc(block.real(CTF_CC));
	ctf->beam_tiltX(block.real(CTF_TILT_X));
	ctf->beam_tiltY(block.real(CTF_TILT_Y));
	if ( block.exists(CTF_ABERRATION_ODD) ) {
		v = parse_real_vector(block[CTF_ABERRATION_ODD].substr(1));
		ctf->aberration_odd_update(v);
	}
	if ( block.exists(CTF_ABERRATION_EVEN) ) {
		v = parse_real_vector(block[CTF_ABERRATION_EVEN].substr(1));
		ctf->aberration_even_update(v);
	}
	ctf->alpha(block.real(CTF_ALPHA));
	ctf->dE(block.real(CTF_DE));

	if ( block.exists(CTF_BASELINE) )
		ctf->parse_baseline_equation(block[CTF_BASELINE]);
	
	if ( block.exists(CTF_ENVELOPE) )
		ctf->parse_envelope_equation(block[CTF_ENVELOPE]);

	return ctf;
}

Bmicrograph*	micrograph_from_starblock(BstarBlock& block, Bstring& mgid, FOMType fom_tag[NFOM])
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG micrograph_from_starblock:" << endl;
	
	Bmicrograph*	mg = NULL;
	mg = micrograph_add(&mg, mgid);
	if ( block.exists(MICROGRAPH_FILE) ) mg->fmg = block[MICROGRAPH_FILE];
	if ( block.exists(MICROGRAPH_FRAMES_FILE) ) mg->fframe = block[MICROGRAPH_FRAMES_FILE];
	if ( block.exists(PARTICLE_FILE) )
		if ( !block.exists_loop(PARTICLE_FILE) )
			mg->fpart = block[PARTICLE_FILE];
	if ( block.exists(FILAMENT_FILE) ) mg->ffil = block[FILAMENT_FILE];
	if ( block.exists(MICROGRAPH_TRANSFORM_FILE) ) mg->fft = block[MICROGRAPH_TRANSFORM_FILE];
	if ( block.exists(MICROGRAPH_POWERSPEC_FILE) ) mg->fps = block[MICROGRAPH_POWERSPEC_FILE];
	mg->img_num = block.integer(MICROGRAPH_NUMBER);
	mg->select = block.integer(MICROGRAPH_SELECT);
	mg->fom = block.real(MICROGRAPH_FOM);
	mg->magnification = block.real(MICROGRAPH_MAGNIFICATION);
	mg->sampling = block.real(MICROGRAPH_SAMPLING);
	if ( block.exists(MICROGRAPH_PIXEL) ) {
		mg->pixel_size[0] = mg->pixel_size[1] = block.real(MICROGRAPH_PIXEL);
	} else {
		mg->pixel_size[0] = block.real(MICROGRAPH_PIXEL_X);
		mg->pixel_size[1] = block.real(MICROGRAPH_PIXEL_Y);
	}
	if ( block.exists(MICROGRAPH_FRAME_PIXEL_X) ) {
		mg->frame_pixel_size[0] = block.real(MICROGRAPH_FRAME_PIXEL_X);
		mg->frame_pixel_size[1] = block.real(MICROGRAPH_FRAME_PIXEL_Y);
	}
	mg->dose = block.real(MICROGRAPH_DOSE);
	mg->exposure = block.real(MICROGRAPH_EXPOSURE);
	mg->intensity = block.real(MICROGRAPH_INTENSITY);
	mg->wri = block.real(MICROGRAPH_WATER_RING);
	mg->origin[0] = block.real(MICROGRAPH_ORIGIN_X);
	mg->origin[1] = block.real(MICROGRAPH_ORIGIN_Y);
	mg->origin[2] = block.real(MICROGRAPH_ORIGIN_Z);
	mg->scale[0] = block.real(MICROGRAPH_SCALE_X);
	mg->scale[1] = block.real(MICROGRAPH_SCALE_Y);
	mg->scale[2] = block.real(MICROGRAPH_SCALE_Z);
	mg->tilt_axis = block.real(MICROGRAPH_TILT_AXIS);
	mg->tilt_angle = block.real(MICROGRAPH_TILT_ANGLE);
	mg->level_angle = block.real(MICROGRAPH_LEVEL_ANGLE);
	mg->rot_angle = block.real(MICROGRAPH_ROT_ANGLE);
	mg->matrix[0][0] = block.real(MICROGRAPH_MATRIX_1_1);
	mg->matrix[0][1] = block.real(MICROGRAPH_MATRIX_1_2);
	mg->matrix[0][2] = block.real(MICROGRAPH_MATRIX_1_3);
	mg->matrix[1][0] = block.real(MICROGRAPH_MATRIX_2_1);
	mg->matrix[1][1] = block.real(MICROGRAPH_MATRIX_2_2);
	mg->matrix[1][2] = block.real(MICROGRAPH_MATRIX_2_3);
	mg->matrix[2][0] = block.real(MICROGRAPH_MATRIX_3_1);
	mg->matrix[2][1] = block.real(MICROGRAPH_MATRIX_3_2);
	mg->matrix[2][2] = block.real(MICROGRAPH_MATRIX_3_3);
	mg->hvec[0] = block.real(MICROGRAPH_HVEC_X);
	mg->hvec[1] = block.real(MICROGRAPH_HVEC_Y);
	mg->hvec[2] = block.real(MICROGRAPH_HVEC_Z);
	mg->kvec[0] = block.real(MICROGRAPH_KVEC_X);
	mg->kvec[1] = block.real(MICROGRAPH_KVEC_Y);
	mg->kvec[2] = block.real(MICROGRAPH_KVEC_Z);
	mg->lvec[0] = block.real(MICROGRAPH_LVEC_X);
	mg->lvec[1] = block.real(MICROGRAPH_LVEC_Y);
	mg->lvec[2] = block.real(MICROGRAPH_LVEC_Z);
	mg->helix_axis = block.real(MICROGRAPH_HELIX_AXIS);
	mg->helix_rise = block.real(MICROGRAPH_HELIX_RISE);
	mg->helix_angle = block.real(MICROGRAPH_HELIX_ANGLE);
	mg->helix_radius = block.real(MICROGRAPH_HELIX_RADIUS);
	mg->tilt_axis *= M_PI/180.0;
	mg->tilt_angle *= M_PI/180.0;
	mg->level_angle *= M_PI/180.0;
	mg->rot_angle  *= M_PI/180.0;
	mg->helix_axis  *= M_PI/180.0;
	mg->helix_angle  *= M_PI/180.0;

	if ( block.exists(PARTICLE_BOX_RADIUS) )
		mg->box_size[0] = mg->box_size[1] = mg->box_size[2] = 
			2 * block.integer(PARTICLE_BOX_RADIUS);
	if ( block.exists(PARTICLE_BOX_RADIUS_X) ) {
		mg->box_size[0] = 2 * block.integer(PARTICLE_BOX_RADIUS_X);
		mg->box_size[1] = 2 * block.integer(PARTICLE_BOX_RADIUS_Y);
		mg->box_size[2] = 2 * block.integer(PARTICLE_BOX_RADIUS_Z);
	}
	if ( block.exists(PARTICLE_BOX_SIZE) )
		mg->box_size[0] = mg->box_size[1] = mg->box_size[2] =
			block.integer(PARTICLE_BOX_SIZE);
	if ( block.exists(PARTICLE_BOX_SIZE_X) ) {
		mg->box_size[0] = block.integer(PARTICLE_BOX_SIZE_X);
		mg->box_size[1] = block.integer(PARTICLE_BOX_SIZE_Y);
		mg->box_size[2] = block.integer(PARTICLE_BOX_SIZE_Z);
	}

	mg->filament_width = block.real(FILAMENT_WIDTH);
	mg->fil_node_radius = block.real(FILNODE_RADIUS);
	mg->bad_radius = block.real(PARTICLE_BAD_RADIUS);
	mg->sf_radius = block.real(REFLEX_RADIUS);
	mg->mark_radius = block.real(MARKER_RADIUS);
	
	mg->ctf = ctf_from_starblock(block);
	
	mg->frame = frame_from_starblock(block);
	
	if ( block.exists_loop(PARTICLE_ID) )
		mg->part = particle_from_starblock(block, fom_tag);

	if ( block.exists_loop(FILAMENT_ID) )
		mg->fil = filament_from_starblock(block);
			
	if ( block.exists_loop(PARTICLE_BAD_X) )
		mg->bad = badarea_from_starblock(block);

	if ( block.exists_loop(MARKER_X) )
		mg->mark = marker_from_starblock(block);

	if ( block.exists_loop(REFLEX_H) )
		mg->sf = strucfac_from_starblock(block);

	if ( block.exists_loop(LAYERLINE_NUMBER) )
		mg->layer = layerline_from_starblock(block);

	return mg;
}

Breconstruction*	reconstruction_from_starblock(BstarBlock& block, FOMType fom_tag[NFOM])
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock:" << endl;
	
	if ( verbose & VERB_FULL )
		cout << "Converting reconstruction" << endl;
	
	Bstring				id("1");
	Breconstruction*	rec = NULL;
	reconstruction_add(&rec, id);
	rec->select = 1;
	
	if ( block.exists(MAP_ID) ) rec->id = block[MAP_ID];
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: ID = " << rec->id << endl;

	if ( block.exists(MAP_RECONSTRUCTION) ) rec->frec = block[MAP_RECONSTRUCTION];
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: frec = " << rec->frec << endl;
	
	if ( block.exists(MAP_TRANSFORM_FILE) ) rec->fft = block[MAP_TRANSFORM_FILE];
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: fft = " << rec->fft << endl;
	
	if ( block.exists(MAP_POWERSPEC_FILE) ) rec->fps = block[MAP_POWERSPEC_FILE];
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: fps = " << rec->fps << endl;
	
	if ( block.exists(PARTICLE_FILE) )
		if ( !block.exists_loop(PARTICLE_FILE) )
			rec->fpart = block[PARTICLE_FILE];
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: fpart = " << rec->fpart << endl;
	
	if ( block.exists(FILAMENT_FILE) ) rec->ffil = block[FILAMENT_FILE];
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: fpart = " << rec->ffil << endl;
	
	if ( block.exists(MAP_MODEL) ) string_add(&rec->model, block[MAP_MODEL].c_str());
	
	rec->select = block.integer(MAP_SELECT);
	rec->fom = block.real(MAP_FOM);

	rec->origin[0] = block.real(MAP_ORIGIN_X);
	rec->origin[1] = block.real(MAP_ORIGIN_Y);
	rec->origin[2] = block.real(MAP_ORIGIN_Z);
	rec->scale[0] = block.real(MAP_SCALE_X);
	rec->scale[1] = block.real(MAP_SCALE_Y);
	rec->scale[2] = block.real(MAP_SCALE_Z);
	if ( block.exists(MAP_VOXEL_SIZE) ) {
		rec->voxel_size[0] = block.real(MAP_VOXEL_SIZE);
		rec->voxel_size[1] = rec->voxel_size[2] = rec->voxel_size[0];
	} else {
		rec->voxel_size[0] = block.real(MAP_VOXEL_SIZE_X);
		rec->voxel_size[1] = block.real(MAP_VOXEL_SIZE_Y);
		rec->voxel_size[2] = block.real(MAP_VOXEL_SIZE_Z);
	}

	if ( block.exists(MAP_SYMMETRY) ) rec->symmetry = block[MAP_SYMMETRY];
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: symmetry = " << rec->symmetry << endl;

	rec->ctf = ctf_from_starblock(block);
	
	if ( block.exists(PARTICLE_BOX_RADIUS) )
		rec->box_size[0] = rec->box_size[1] = rec->box_size[2] = 
			2 * block.integer(PARTICLE_BOX_RADIUS);
	if ( block.exists(PARTICLE_BOX_RADIUS_X) ) {
		rec->box_size[0] = 2 * block.integer(PARTICLE_BOX_RADIUS_X);
		rec->box_size[1] = 2 * block.integer(PARTICLE_BOX_RADIUS_Y);
		rec->box_size[2] = 2 * block.integer(PARTICLE_BOX_RADIUS_Z);
	}
	if ( block.exists(PARTICLE_BOX_SIZE) )
		rec->box_size[0] = rec->box_size[1] = rec->box_size[2] =
			block.integer(PARTICLE_BOX_SIZE);
	if ( block.exists(PARTICLE_BOX_SIZE_X) ) {
		rec->box_size[0] = block.integer(PARTICLE_BOX_SIZE_X);
		rec->box_size[1] = block.integer(PARTICLE_BOX_SIZE_Y);
		rec->box_size[2] = block.integer(PARTICLE_BOX_SIZE_Z);
	}

	rec->filament_width = block.real(FILAMENT_WIDTH);
	rec->fil_node_radius = block.real(FILNODE_RADIUS);
	rec->bad_radius = block.real(PARTICLE_BAD_RADIUS);
	rec->sf_radius = block.real(REFLEX_RADIUS);
	rec->mark_radius = block.real(MARKER_RADIUS);

	rec->view[0] = block.real(MAP_VIEW_X);
	rec->view[1] = block.real(MAP_VIEW_Y);
	rec->view[2] = block.real(MAP_VIEW_Z);
	rec->view[3] = block.real(MAP_VIEW_ANGLE)*M_PI/180.0;

	if ( block.exists(PARTICLE_ID) )
		rec->part = particle_from_starblock(block, fom_tag);
	
	if ( block.exists(FILAMENT_ID) )
		rec->fil = filament_from_starblock(block);
			
	if ( block.exists(PARTICLE_BAD_X) )
		rec->bad = badarea_from_starblock(block);

	if ( block.exists(MARKER_X) )
		rec->mark = marker_from_starblock(block);
	
	if ( block.exists(REFLEX_H) )
		rec->sf = strucfac_from_starblock(block);
	
	return rec;
}

/*
@brief 	Converts micrograph data from a STAR data base to a project structure.
@param	&star		STAR data base.
@param 	*project	project structure.
@param	flag		update tags.
@return int			error code (<0 means failure).

	The function sets up the project hierarchy from a STAR database file.
	Micrograph data blocks must contain the "micrograph.id" tag.
	If the "micrograph.field_of_view_id" is present, the micrographs are
	distributed into multi-micrograph fields-of-view. Otherwise, each
	micrograph is assumed to be a unique field-of-view.
	All angles given in the STAR data base are assumed to be in degrees, 
	and are converted here only to radians.

**/
int			star_to_project(Bstar& star, Bproject* project, int flag)
{
	int				err(0);
	int 			nmg(0), nrec(0);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG star_to_project:" << endl;
	
	if ( flag ) mg_star_update_tags(star);
		
	nmg = star.number_of_blocks(MICROGRAPH_ID);
	if ( nmg < 1 ) nmg = star.number_of_blocks(MICROGRAPH_FILE);
	if ( nmg < 1 ) nmg = star.number_of_blocks(MICROGRAPH_FRAMES_FILE);
	if ( nmg < 1 ) nmg = star.number_of_blocks(MICROGRAPH_TRANSFORM_FILE);
	if ( nmg < 1 ) nmg = star.number_of_blocks(MICROGRAPH_POWERSPEC_FILE);

	nrec = star.number_of_blocks(MAP_ID);
	if ( nrec < 1 ) nrec = star.number_of_blocks(MAP_RECONSTRUCTION);
	if ( nrec < 1 ) nrec = star.number_of_blocks(MAP_TRANSFORM_FILE);
	if ( nrec < 1 ) nrec = star.number_of_blocks(MAP_POWERSPEC_FILE);
	
	if ( nmg < 1 && nrec < 1 ) {
		cerr << "Warning: No micrograph or reconstruction data blocks found" << endl;
		return 0;
	}

	BstarBlock		block;
	Bfield* 		field = NULL;
	Bmicrograph		*mg, *mg2;
	Breconstruction	*rec, *rec2;
	
	for ( rec2 = project->rec; rec2 && rec2->next; rec2 = rec2->next ) ;

	if ( star.exists(MAP_REFERENCE) ) {
		block = star.find(MAP_REFERENCE);
		for ( auto il: block.loops() )
			for ( auto ir: il.data() )
				string_add(&project->reference, ir[0].c_str());
	}
	
	if ( ( verbose & VERB_DEBUG ) && project->reference )
		cout << "DEBUG star_to_project: reference map = " << *(project->reference) << endl;
	
	int 			i, j, f, nfield(0);
	Bstring			field_id, mgid;
	Bstring			base;
	
	// Make sure FOM tags are properly set up
	if ( star.number_of_blocks(PARTICLE_FOM) ) project->fom_tag[0] = FOM;
	if ( star.number_of_blocks(PARTICLE_FOM_CV) ) project->fom_tag[1] = FOM_CV;
	if ( star.number_of_blocks(PARTICLE_HANDA_FOM) ) project->fom_tag[1] = FOM_HAND_A;
	if ( star.number_of_blocks(PARTICLE_HANDB_FOM) ) project->fom_tag[2] = FOM_HAND_B;
	
	project_get_fom_tags(star, project->fom_tag);
	for ( f=i=0; f<NFOM; f++ ) i += project->fom_tag[f];
	if ( i < 2 ) project->fom_tag[0] = FOM;
	if ( verbose & VERB_DEBUG )
		for ( f=0; f<NFOM; f++ ) if ( project->fom_tag[f] )
			cout << "DEBUG star_to_project: fom_tag=" << project->fom_tag[f] << endl;

	// The orientation convention used is set based on whether the view angle
	// array or the phi euler angle array is present in the file.
	// If both are present, the one with the larger number of values wins,
	// otherwise only the view vector and rotation angle
	// representation will be used.
	i = star.number_of_blocks(PARTICLE_PHI);
	j = star.number_of_blocks(PARTICLE_VIEW_ANGLE);
	project->euler_flag = 0;
	project->omega_flag = 0;
	if ( i > 0 && j > 0 )
		cerr << "Error: The STAR file contains both views and Euler angles!" << endl;
	else if ( i > j ) {
		project->euler_flag = 1;
		if ( star.number_of_blocks(PARTICLE_OMEGA) > 0 )
			project->omega_flag = 1;
	}
	
	nmg = i = 0;
	for ( auto ib: star.blocks() ) {
		if ( ib.exists(MICROGRAPH_ID) || ib.exists(MICROGRAPH_FILE) ) {
			mgid = 0;
			field_id = 0;
			/* This deals with old style access to micrograph data blocks */
			if ( ib.exists(MICROGRAPH_ID) ) {
				mgid = ib[MICROGRAPH_ID];
			} else {
				if ( ib.exists(MICROGRAPH_FILE) ) mgid = ib[MICROGRAPH_FILE];
				else if ( ib.exists(MICROGRAPH_FRAMES_FILE) ) mgid = ib[MICROGRAPH_FRAMES_FILE];
				else if ( ib.exists(PARTICLE_FILE) ) mgid = ib[PARTICLE_FILE];
				else if ( ib.exists(FILAMENT_FILE) ) mgid = ib[FILAMENT_FILE];
				else if ( ib.exists(MICROGRAPH_TRANSFORM_FILE) ) mgid = ib[MICROGRAPH_TRANSFORM_FILE];
				else if ( ib.exists(MICROGRAPH_POWERSPEC_FILE) ) mgid = ib[MICROGRAPH_POWERSPEC_FILE];
				if ( mgid.length() > 0 ) mgid = mgid.base();
			}
			if ( mgid.length() > 0 ) {
				if ( ib.exists(MICROGRAPH_FIELD_ID) ) field_id = ib[MICROGRAPH_FIELD_ID];
				if ( field_id.length() < 1 ) field_id = mgid;
				if ( !project->field ) {
					field_add(&field, field_id);
					project->field = field;	// First field
					nfield++;
				} else {
					field = NULL;
					if ( field_id.length() > 0 )
						field = field_find_id(project->field, field_id);
					if ( !field ) {
						field = project->field;
						field = field_add(&field, field_id);
						nfield++;
					}
				}
				mg = micrograph_from_starblock(ib, mgid, project->fom_tag);
				mg->block = i++;
				if ( verbose & VERB_FULL )
					cout << "Reading field \"" << field->id << "\", micrograph \"" << mg->id <<"\"" << endl;
				if ( field->mg ) {
					for ( mg2 = field->mg; mg2->next; mg2 = mg2->next ) ;
					mg2->next = mg;
				} else field->mg = mg;
				mgid = 0;
				field_id = 0;
				nmg++;
			}
		} else if ( ib.exists(MAP_ID) ||
				 ib.exists(MAP_RECONSTRUCTION) ||
				 ib.exists(MAP_TRANSFORM_FILE) ||
				 ib.exists(MAP_POWERSPEC_FILE) ) {
			rec = reconstruction_from_starblock(ib, project->fom_tag);
			if ( !project->rec ) project->rec = rec;
			else rec2->next = rec;
			rec2 = rec;
		} else if ( ib.exists_loop(PARTICLE_ID) ) {
			project->class_avg = particle_from_starblock(ib, project->fom_tag);
		}
	}
	
	star.replace_tag(CTF_DEF_MIN, CTF_DEF_AVG);
	star.replace_tag(CTF_DEF_MAX, CTF_DEF_DEV);
			
	return err;
}

int		frame_to_starblock(Bframe* frame, BstarBlock& block)
{
	BstarLoop&		loop = block.add_loop();
	loop.tags()[MICROGRAPH_FRAME] = 0;
	loop.tags()[MICROGRAPH_FRAME_SHIFT_X] = 1;
	loop.tags()[MICROGRAPH_FRAME_SHIFT_Y] = 2;
	loop.tags()[MICROGRAPH_FRAME_SELECT] = 3;
	loop.tags()[MICROGRAPH_FRAME_FOM] = 4;
	for ( Bframe* f = frame; f; f = f->next ) {
		vector<string>&	vs = loop.add_row(5);
		vs[0] = to_string(f->id);
		vs[1] = to_string(f->shift[0]);
		vs[2] = to_string(f->shift[1]);
		vs[3] = to_string(f->sel);
		vs[4] = to_string(f->fom);
	}
	
	return 0;
}

int		particle_to_starblock(Bparticle* part, BstarBlock& block,
				FOMType fom_tag[NFOM], int euler_flag, int omega_flag)
{
	Bparticle*		p;
	long			f, i, nt(0);
	Bstring			tag;
	Euler			euler;

	if (verbose & VERB_DEBUG )
		cout << "DEBUG particle_to_starblock: fom tags:" << endl;
	
	BstarLoop&		loop = block.add_loop();
	loop.tags()[PARTICLE_ID] = 0;
	loop.tags()[PARTICLE_GROUP] = 1;
	loop.tags()[PARTICLE_DEFOCUS] = 2;
	loop.tags()[PARTICLE_DEF_DEV] = 3;
	loop.tags()[PARTICLE_AST_ANG] = 4;
	loop.tags()[PARTICLE_MAGNIF] = 5;
	loop.tags()[PARTICLE_X] = 6;
	loop.tags()[PARTICLE_Y] = 7;
	loop.tags()[PARTICLE_Z] = 8;
	loop.tags()[PARTICLE_PIXEL_X] = 9;
	loop.tags()[PARTICLE_PIXEL_Y] = 10;
	loop.tags()[PARTICLE_PIXEL_Z] = 11;
	loop.tags()[PARTICLE_ORIGIN_X] = 12;
	loop.tags()[PARTICLE_ORIGIN_Y] = 13;
	loop.tags()[PARTICLE_ORIGIN_Z] = 14;
	if ( euler_flag < 1 ) {
		loop.tags()[PARTICLE_VIEW_X] = 15;
		loop.tags()[PARTICLE_VIEW_Y] = 16;
		loop.tags()[PARTICLE_VIEW_Z] = 17;
		loop.tags()[PARTICLE_VIEW_ANGLE] = 18;
		nt = 19;
	} else {
		loop.tags()[PARTICLE_PHI] = 15;
		loop.tags()[PARTICLE_THETA] = 16;
		loop.tags()[PARTICLE_PSI] = 17;
		nt = 18;
	}
	loop.tags()[PARTICLE_SELECT] = nt++;
	for ( f=0; f<NFOM && fom_tag[f]; ++f ) {
		tag = get_fom_tag(fom_tag[f]);
		loop.tags()[tag.str()] = nt++;
		if ( verbose & VERB_DEBUG )
			cout << nt << tab << tag << endl;
	}
	if ( part->fpart.length() )
		loop.tags()[PARTICLE_FILE] = nt++;
	
	if (verbose & VERB_DEBUG )
		cout << "DEBUG particle_to_starblock: columns=" << nt << endl;
	
	for ( p = part; p; p = p->next ) {
//		cout << p->id << tab << "-" << p->fpart << "-" << endl;
		vector<string>&	vs = loop.add_row(nt);
		vs[0] = to_string(p->id);
		vs[1] = to_string(p->group);
		vs[2] = to_string(p->def);
		vs[3] = to_string(p->dev);
		vs[4] = to_string(p->ast*180.0/M_PI);
		vs[5] = to_string(p->mag);
		vs[6] = to_string(p->loc[0]);
		vs[7] = to_string(p->loc[1]);
		vs[8] = to_string(p->loc[2]);
		vs[9] = to_string(p->pixel_size[0]);
		vs[10] = to_string(p->pixel_size[1]);
		vs[11] = to_string(p->pixel_size[2]);
		vs[12] = to_string(p->ori[0]);
		vs[13] = to_string(p->ori[1]);
		vs[14] = to_string(p->ori[2]);
		if ( euler_flag < 1 ) {
			vs[15] = to_string(p->view[0]);
			vs[16] = to_string(p->view[1]);
			vs[17] = to_string(p->view[2]);
			vs[18] = to_string(p->view[3]*180.0/M_PI);
			i = 19;
		} else {
			euler = Euler(p->view);
			vs[15] = to_string(euler.phi()*180.0/M_PI);
			vs[16] = to_string(euler.theta()*180.0/M_PI);
			vs[17] = to_string(euler.psi()*180.0/M_PI);
			i = 18;
		}
		vs[i++] = to_string(p->sel);
		for ( f=0; i<nt && f<NFOM && fom_tag[f]; ++f, ++i )
			vs[i] = to_string(p->fom[f]);
		if ( p->fpart.length() )
			vs[i] = p->fpart.str();
	}
		
	if (verbose & VERB_DEBUG )
		cout << "DEBUG particle_to_starblock: columns converted=" << i << endl;
	
	return 0;
}

int		badarea_to_starblock(Bbadarea* bad, BstarBlock& block)
{
	BstarLoop&		loop = block.add_loop();
	loop.tags()[PARTICLE_BAD_X] = 0;
	loop.tags()[PARTICLE_BAD_Y] = 1;
	loop.tags()[PARTICLE_BAD_Z] = 2;

	for ( Bbadarea* b = bad; b; b = b->next ) {
		vector<string>&	vs = loop.add_row(3);
		vs[0] = to_string(b->loc[0]);
		vs[1] = to_string(b->loc[1]);
		vs[2] = to_string(b->loc[2]);
	}
	
	return 0;
}

int		filament_to_starblock(Bfilament* fil, BstarBlock& block)
{
	BstarLoop&		loop = block.add_loop();
	loop.tags()[FILAMENT_ID] = 0;
	loop.tags()[FILAMENT_NODE_ID] = 1;
	loop.tags()[FILAMENT_NODE_X] = 2;
	loop.tags()[FILAMENT_NODE_Y] = 3;
	loop.tags()[FILAMENT_NODE_Z] = 4;

	for ( Bfilament* f = fil; f; f=f->next ) {
		for ( Bfilnode* fn = f->node; fn; fn = fn->next ) {
			vector<string>&	vs = loop.add_row(5);
			vs[0] = to_string(f->id);
			vs[1] = to_string(fn->id);
			vs[2] = to_string(fn->loc[0]);
			vs[3] = to_string(fn->loc[1]);
			vs[4] = to_string(fn->loc[2]);
		}
	}

	return 0;
}

int		marker_to_starblock(Bmarker* mark, BstarBlock& block)
{
	BstarLoop&		loop = block.add_loop();
	loop.tags()[MARKER_ID] = 0;
	loop.tags()[MARKER_X] = 1;
	loop.tags()[MARKER_Y] = 2;
	loop.tags()[MARKER_Z] = 3;
	loop.tags()[MARKER_ERROR_X] = 4;
	loop.tags()[MARKER_ERROR_Y] = 5;
	loop.tags()[MARKER_ERROR_Z] = 6;
	loop.tags()[MARKER_RESIDUAL] = 7;
	loop.tags()[MARKER_FOM] = 8;
	loop.tags()[MARKER_SELECT] = 9;

	for ( Bmarker* m = mark; m; m = m->next ) {
		vector<string>&	vs = loop.add_row(10);
		vs[0] = to_string(m->id);
		vs[1] = to_string(m->loc[0]);
		vs[2] = to_string(m->loc[1]);
		vs[3] = to_string(m->loc[2]);
		vs[4] = to_string(m->err[0]);
		vs[5] = to_string(m->err[1]);
		vs[6] = to_string(m->err[2]);
		vs[7] = to_string(m->res);
		vs[8] = to_string(m->fom);
		vs[9] = to_string(m->sel);
	}
	
	return 0;
}

int		strucfac_to_starblock(Bstrucfac* sf, BstarBlock& block)
{
	BstarLoop&		loop = block.add_loop();
	loop.tags()[REFLEX_H] = 0;
	loop.tags()[REFLEX_K] = 1;
	loop.tags()[REFLEX_L] = 2;
	loop.tags()[REFLEX_X] = 3;
	loop.tags()[REFLEX_Y] = 4;
	loop.tags()[REFLEX_Z] = 5;
	loop.tags()[REFLEX_AMP] = 6;
	loop.tags()[REFLEX_SIGAMP] = 7;
	loop.tags()[REFLEX_PHI] = 8;
	loop.tags()[REFLEX_SIGPHI] = 9;
	loop.tags()[REFLEX_FOM] = 10;
	loop.tags()[REFLEX_STATUS] = 11;

	for ( Bstrucfac* f = sf; f; f = f->next ) {
		vector<string>&	vs = loop.add_row(12);
		vs[0] = to_string(f->index[0]);
		vs[1] = to_string(f->index[1]);
		vs[2] = to_string(f->index[2]);
		vs[3] = to_string(f->loc[0]);
		vs[4] = to_string(f->loc[1]);
		vs[5] = to_string(f->loc[2]);
		vs[6] = to_string(f->amp);
		vs[7] = to_string(f->sigamp);
		vs[8] = to_string(f->phi);
		vs[9] = to_string(f->sigphi);
		vs[10] = to_string(f->fom);
		vs[11] = to_string(f->sel);
	}
	
	return 0;
}

int		layerline_to_starblock(Blayerline* line, BstarBlock& block)
{
	BstarLoop&		loop = block.add_loop();
	loop.tags()[LAYERLINE_NUMBER] = 0;
	loop.tags()[LAYERLINE_ORDER] = 1;
	loop.tags()[LAYERLINE_DISTANCE] = 2;
	loop.tags()[LAYERLINE_FREQ] = 3;
	loop.tags()[LAYERLINE_AMP] = 4;
	loop.tags()[LAYERLINE_FOM] = 5;
	loop.tags()[LAYERLINE_SELECT] = 6;
	
	for ( Blayerline* ll = line; ll; ll = ll->next ) {
		vector<string>&	vs = loop.add_row(7);
		vs[0] = to_string(ll->number);
		vs[1] = to_string(ll->order);
		vs[2] = to_string(ll->distance);
		vs[3] = to_string(ll->freq);
		vs[4] = to_string(ll->amp);
		vs[5] = to_string(ll->fom);
		vs[6] = to_string(ll->sel);
	}
	
	return 0;
}

int		ctf_to_starblock(CTFparam* ctf, BstarBlock& block)
{
	float		angle;
	Bstring		s;
	
	block[CTF_ID] = ctf->identifier();
	block[CTF_SELECT] = to_string(ctf->select());
	block[CTF_FOM] = to_string(ctf->fom());
	block[CTF_VOLTAGE] = to_string(ctf->volt());
	block[CTF_FOCAL] = to_string(ctf->focal_length());
	block[CTF_APERTURE] = to_string(ctf->objective_aperture());
	block[CTF_CS] = to_string(ctf->Cs());
	block[CTF_CC] = to_string(ctf->Cc());
	block[CTF_TILT_X] = to_string(ctf->beam_tiltX());
	block[CTF_TILT_Y] = to_string(ctf->beam_tiltY());
	if ( ctf->aberration_odd().size() )
		block[CTF_ABERRATION_ODD] = "[" + concatenate(ctf->aberration_odd()) + "]";
	if ( ctf->aberration_even().size() )
		block[CTF_ABERRATION_EVEN] = "[" + concatenate(ctf->aberration_even_difference()) + "]";
	block[CTF_ALPHA] = to_string(ctf->alpha());
	block[CTF_DE] = to_string(ctf->dE());
	block[CTF_AMP_SHIFT] = to_string(ctf->amp_shift());
	block[CTF_DEF_AVG] = to_string(ctf->defocus_average());
	block[CTF_DEF_DEV] = to_string(ctf->defocus_deviation());
	angle = ctf->astigmatism_angle()*180.0/M_PI;
	block[CTF_AST_ANG] = to_string(angle);
	block[CTF_ZERO] = to_string(ctf->zero(1));
	if ( ctf->baseline_type() ) {
		s = ctf->baseline_equation();
		block[CTF_BASELINE] = s.str();
	}
	if ( ctf->envelope(0) ) {
		s = ctf->envelope_equation();
		block[CTF_ENVELOPE] = s.str();
	}

	return 0;
}

int		micrograph_to_starblock(Bmicrograph* mg, BstarBlock& block,
				FOMType fom_tag[NFOM], int euler_flag, int omega_flag)
{
	float			angle;

	if (verbose & VERB_DEBUG )
		cout << "DEBUG micrograph_to_starblock: " << mg->id << endl;
	
	block[MICROGRAPH_ID] = mg->id.str();
	block[MICROGRAPH_NUMBER] = to_string(mg->img_num);
	block[MICROGRAPH_SELECT] = to_string(mg->select);
	block[MICROGRAPH_FOM] = to_string(mg->fom);
	if ( mg->fmg.length() ) block[MICROGRAPH_FILE] = mg->fmg.str();
	if ( mg->fframe.length() ) block[MICROGRAPH_FRAMES_FILE] = mg->fframe.str();
	if ( mg->fft.length() ) block[MICROGRAPH_TRANSFORM_FILE] = mg->fft.str();
	if ( mg->fps.length() ) block[MICROGRAPH_POWERSPEC_FILE] = mg->fps.str();
	if ( mg->fpart.length() && mg->part && !mg->part->fpart.length() )
		block[PARTICLE_FILE] = mg->fpart.str();
	if ( mg->ffil.length() ) block[FILAMENT_FILE] = mg->ffil.str();
	block[MICROGRAPH_MAGNIFICATION] = to_string(mg->magnification);
	block[MICROGRAPH_SAMPLING] = to_string(mg->sampling);
	block[MICROGRAPH_PIXEL_X] = to_string(mg->pixel_size[0]);
	block[MICROGRAPH_PIXEL_Y] = to_string(mg->pixel_size[1]);
	block[MICROGRAPH_FRAME_PIXEL_X] = to_string(mg->frame_pixel_size[0]);
	block[MICROGRAPH_FRAME_PIXEL_Y] = to_string(mg->frame_pixel_size[1]);
	block[MICROGRAPH_DOSE] = to_string(mg->dose);
	block[MICROGRAPH_EXPOSURE] = to_string(mg->exposure);
	block[MICROGRAPH_INTENSITY] = to_string(mg->intensity);
	block[MICROGRAPH_WATER_RING] = to_string(mg->wri);
	angle = mg->tilt_axis*180.0/M_PI;
	block[MICROGRAPH_TILT_AXIS] = to_string(angle);
	angle = mg->tilt_angle*180.0/M_PI;
	block[MICROGRAPH_TILT_ANGLE] = to_string(angle);
	angle = mg->level_angle*180.0/M_PI;
	block[MICROGRAPH_LEVEL_ANGLE] = to_string(angle);
	angle = mg->rot_angle*180.0/M_PI;
	block[MICROGRAPH_ROT_ANGLE] = to_string(angle);
	block[MICROGRAPH_ORIGIN_X] = to_string(mg->origin[0]);
	block[MICROGRAPH_ORIGIN_Y] = to_string(mg->origin[1]);
	block[MICROGRAPH_ORIGIN_Z] = to_string(mg->origin[2]);
	block[MICROGRAPH_SCALE_X] = to_string(mg->scale[0]);
	block[MICROGRAPH_SCALE_Y] = to_string(mg->scale[1]);
	block[MICROGRAPH_SCALE_Z] = to_string(mg->scale[2]);
	block[MICROGRAPH_MATRIX_1_1] = to_string(mg->matrix[0][0]);
	block[MICROGRAPH_MATRIX_1_2] = to_string(mg->matrix[0][1]);
	block[MICROGRAPH_MATRIX_1_3] = to_string(mg->matrix[0][2]);
	block[MICROGRAPH_MATRIX_2_1] = to_string(mg->matrix[1][0]);
	block[MICROGRAPH_MATRIX_2_2] = to_string(mg->matrix[1][1]);
	block[MICROGRAPH_MATRIX_2_3] = to_string(mg->matrix[1][2]);
	block[MICROGRAPH_MATRIX_3_1] = to_string(mg->matrix[2][0]);
	block[MICROGRAPH_MATRIX_3_2] = to_string(mg->matrix[2][1]);
	block[MICROGRAPH_MATRIX_3_3] = to_string(mg->matrix[2][2]);
	block[MICROGRAPH_HVEC_X] = to_string(mg->hvec[0]);
	block[MICROGRAPH_HVEC_Y] = to_string(mg->hvec[1]);
	block[MICROGRAPH_HVEC_Z] = to_string(mg->hvec[2]);
	block[MICROGRAPH_KVEC_X] = to_string(mg->kvec[0]);
	block[MICROGRAPH_KVEC_Y] = to_string(mg->kvec[1]);
	block[MICROGRAPH_KVEC_Z] = to_string(mg->kvec[2]);
	block[MICROGRAPH_LVEC_X] = to_string(mg->lvec[0]);
	block[MICROGRAPH_LVEC_Y] = to_string(mg->lvec[1]);
	block[MICROGRAPH_LVEC_Z] = to_string(mg->lvec[2]);
	angle = mg->helix_axis*180.0/M_PI;
	block[MICROGRAPH_HELIX_AXIS] = to_string(angle);
	block[MICROGRAPH_HELIX_RISE] = to_string(mg->helix_rise);
	angle = mg->helix_angle*180.0/M_PI;
	block[MICROGRAPH_HELIX_ANGLE] = to_string(angle);
	block[MICROGRAPH_HELIX_RADIUS] = to_string(mg->helix_radius);
	block[PARTICLE_BOX_SIZE_X] = to_string(mg->box_size[0]);
	block[PARTICLE_BOX_SIZE_Y] = to_string(mg->box_size[1]);
	block[PARTICLE_BOX_SIZE_Z] = to_string(mg->box_size[2]);
	block[PARTICLE_BAD_RADIUS] = to_string(mg->bad_radius);
	block[FILAMENT_WIDTH] = to_string(mg->filament_width);
	block[FILNODE_RADIUS] = to_string(mg->fil_node_radius);
	block[REFLEX_RADIUS] = to_string(mg->sf_radius);
	block[MARKER_RADIUS] = to_string(mg->mark_radius);

	if ( mg->ctf ) ctf_to_starblock(mg->ctf, block);

	if ( mg->frame ) frame_to_starblock(mg->frame, block);
			
	long		npart = particle_count(mg->part);
	long		nfil = filament_node_count(mg->fil);
	long		nbad = count_list((char *) mg->bad);
	long		nmark = count_list((char *) mg->mark);
	long		nsf = count_list((char *) mg->sf);
	long		nline = count_list((char *) mg->layer);
		
	if ( npart > 0 )
		particle_to_starblock(mg->part, block, fom_tag, 
			euler_flag, omega_flag);

	if ( nbad > 0)
		badarea_to_starblock(mg->bad, block);

	if ( nmark > 0 )
		marker_to_starblock(mg->mark, block);

	if ( nfil > 0 )
		filament_to_starblock(mg->fil, block);

	if ( nsf > 0 )
		strucfac_to_starblock(mg->sf, block);

	if ( nline > 0 )
		layerline_to_starblock(mg->layer, block);

	return 0;
}

int		reconstruction_to_starblock(Breconstruction* rec, BstarBlock& block,
				FOMType fom_tag[NFOM], int euler_flag, int omega_flag)
{
	float			angle;

	block[MAP_ID] = rec->id.str();

	if ( rec->frec.length() ) block[MAP_RECONSTRUCTION] = rec->frec.str();

	if ( rec->fft.length() ) block[MAP_TRANSFORM_FILE] = rec->fft.str();

	if ( rec->fps.length() ) block[MAP_POWERSPEC_FILE] = rec->fps.str();

	if ( rec->fpart.length() && !rec->part->fpart.length() )
		block[PARTICLE_FILE] = rec->fpart.str();

	if ( rec->ffil.length() ) block[FILAMENT_FILE] = rec->ffil.str();

	if ( rec->model ) {
		BstarLoop&		loop = block.add_loop();
		loop.tags()[MAP_MODEL] = 0;
		for ( Bstring *s = rec->model; s; s = s->next ) {
			vector<string>&	vs = loop.add_row(1);
			vs[0] = s->str();
		}
	}

	block[MAP_SELECT] = to_string(rec->select);
	block[MAP_FOM] = to_string(rec->fom);

	block[MAP_ORIGIN_X] = to_string(rec->origin[0]);
	block[MAP_ORIGIN_Y] = to_string(rec->origin[1]);
	block[MAP_ORIGIN_Z] = to_string(rec->origin[2]);
	block[MAP_SCALE_X] = to_string(rec->scale[0]);
	block[MAP_SCALE_Y] = to_string(rec->scale[1]);
	block[MAP_SCALE_Z] = to_string(rec->scale[2]);
	block[MAP_VOXEL_SIZE_X] = to_string(rec->voxel_size[0]);
	block[MAP_VOXEL_SIZE_Y] = to_string(rec->voxel_size[1]);
	block[MAP_VOXEL_SIZE_Z] = to_string(rec->voxel_size[2]);

	if ( rec->symmetry.length() ) block[MAP_SYMMETRY] = rec->symmetry.str();

	if ( rec->ctf ) ctf_to_starblock(rec->ctf, block);

	block[PARTICLE_BOX_SIZE_X] = to_string(rec->box_size[0]);
	block[PARTICLE_BOX_SIZE_Y] = to_string(rec->box_size[1]);
	block[PARTICLE_BOX_SIZE_Z] = to_string(rec->box_size[2]);
	block[PARTICLE_BAD_RADIUS] = to_string(rec->bad_radius);
	block[FILAMENT_WIDTH] = to_string(rec->filament_width);
	block[FILNODE_RADIUS] = to_string(rec->fil_node_radius);
	block[REFLEX_RADIUS] = to_string(rec->sf_radius);
	block[MARKER_RADIUS] = to_string(rec->mark_radius);

	block[MAP_VIEW_X] = to_string(rec->view[0]);
	block[MAP_VIEW_Y] = to_string(rec->view[1]);
	block[MAP_VIEW_Z] = to_string(rec->view[2]);
	angle = rec->view.angle()*180.0/M_PI;
	block[MAP_VIEW_ANGLE] = to_string(angle);

	long		npart = particle_count(rec->part);
	long		nfil = filament_node_count(rec->fil);
	long		nbad = count_list((char *) rec->bad);
	long		nmark = count_list((char *) rec->mark);
	long		nsf = count_list((char *) rec->sf);
		
	if ( npart > 0 )
		particle_to_starblock(rec->part, block, fom_tag, euler_flag, omega_flag);
		
	if ( nbad > 0)
		badarea_to_starblock(rec->bad, block);
		
	if ( nfil > 0 )
		filament_to_starblock(rec->fil, block);
	
	if ( nmark > 0 )
		marker_to_starblock(rec->mark, block);

	if ( nsf > 0 )
		strucfac_to_starblock(rec->sf, block);

	return 0;
}

/*
@brief	Converts micrograph data from an image project structure to a
	STAR data base.
@param	*project		project parameter structure.
@param	&star			STAR data base.
@param 	mg_select		flag to only convert selected micrographs.
@param 	rec_select		flag to only convert selected reconstructions.
@return	int				error code (<0 means failure).

	The function packes new micrograph data back into an existing STAR
	data base from which the old project parameters were obtained.
	All angles used internally are in radians, and are converted here
	to degrees for output. The angles within the micrograph structure
	are left in radians.
**/
int 		project_to_star(Bproject* project, Bstar& star, int mg_select, int rec_select)
{
	int				err(0);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_to_star:" << endl;

	star.comment(project->comment.str());
	
//	star->split = project->split;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bstring*			file_list = NULL;

	if ( project->reference ) {
		BstarBlock&		block = star.add_block(MAP_REFERENCE);
		BstarLoop&		loop = block.add_loop();
		loop.tags()[MAP_REFERENCE] = 0;
		for ( file_list=project->reference; file_list; file_list=file_list->next ) {
			vector<string>&	vs = loop.add_row(1);
			vs[0] = file_list->str();
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG project_to_star: reference map = " << *file_list << endl;
		}
	}
	
	long	 			i, f, nfield(0), nmg(0), nrec(0), npart(0);

	for ( f=i=0; f<NFOM; f++ ) i += project->fom_tag[f];
	if ( i < 2 ) project->fom_tag[0] = FOM;
	
	if ( verbose & VERB_PROCESS )
		cout << "Converting a project into a STAR database" << endl;
	
	for ( field=project->field; field; field=field->next ) if ( mg_select < 1 || field->select > 0 ) {
		for ( mg=field->mg; mg; mg=mg->next ) if ( mg_select < 1 || mg->select > 0 ) {
			if ( verbose & VERB_FULL )
				cout << "Writing field \"" << field->id << "\", micrograph \"" << mg->id << "\"" << endl;
			BstarBlock&		block = star.add_block(mg->id.str());
			block[MICROGRAPH_FIELD_ID] = field->id.str();
			micrograph_to_starblock(mg, block, project->fom_tag, project->euler_flag, project->omega_flag);
			nmg++;
		}
		nfield++;
	}
	
	if ( project->class_avg ) {
		npart = particle_count(project->class_avg);
		if ( verbose & VERB_FULL )
			cout << "Writing " << npart << " class averages" << endl;
		BstarBlock&		block = star.add_block(CLASS_AVERAGE);
		particle_to_starblock(project->class_avg, block, project->fom_tag,
			project->euler_flag, project->omega_flag);
	}
	
	for ( rec = project->rec; rec; rec = rec->next ) if ( rec_select < 1 || rec->select > 0 ) {
		if ( verbose & VERB_FULL )
			cout << "Writing reconstruction \"" << rec->id << "\"" << endl;
		BstarBlock&		block = star.add_block(rec->id.str());
		reconstruction_to_starblock(rec, block, project->fom_tag,
			project->euler_flag, project->omega_flag);
		nrec++;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << nfield << " fields, " << nmg << " micrographs, " << nrec << " reconstructions written" << endl;
	
	return err;
}

