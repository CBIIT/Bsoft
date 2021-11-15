/**
@file	rwmgSTAR.cpp
@brief	Library routines to read and write micrograph parameters in STAR format
@author Bernard Heymann
@date	Created: 20010206
@date	Modified: 20200202
**/

#include "mg_processing.h"
#include "rwmgSTAR.h"
#include "rwstar.h"
#include "mg_tags.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
int			star_to_project(Bstar* star, Bproject* project);
int 		project_to_star(Bproject* project, Bstar* star, int mg_select, int rec_select);

/**
@brief 	Reading micrograph parameters from STAR files.
@param 	&filename		file name (or comma-delimited list).
@param 	*project		initialized project structure.
@return int						error code (<0 means failure).
**/
int			read_project_star(Bstring& filename, Bproject* project)
{
	int				err(0);
	Bstar*          star = init_star();
	star->line_length = 200;                // Set the output line length

	err = read_star(filename, star);
	
	if ( err == 0 ) {
		if ( project->comment.length() < 1 )
			project->comment = star->comment;
		err = star_to_project(star, project);
	}
	
	kill_star(star);
	
	return err;
}

/**
@brief 	Writing micrograph parameters to a STAR file.
@param 	&filename		file name.
@param 	*project		project structure.
@param 	mg_select			flag to only write selected micrographs.
@param 	rec_select			flag to only convert selected reconstructions.
@return int						error code (<0 means failure).
**/
int			write_project_star(const char* filename, Bproject* project, int mg_select, int rec_select)
{
	Bstring		thefile(filename);
	return write_project_star(thefile, project, mg_select, rec_select);
}

int			write_project_star(Bstring& filename, Bproject* project, int mg_select, int rec_select)
{
	int				err(0);
	Bstar*          star = init_star();
	star->line_length = 200;                // Set the output line length

	err = project_to_star(project, star, mg_select, rec_select);
	
	if ( err == 0 )
		err = write_star(filename, star);
		
	kill_star(star);
	
	return err;
}

/**
@brief 	Replacing old STAR tags with new ones.
@param	*star				STAR database.
@return int						0.
**/
int			mg_star_update_tags(Bstar* star)
{
	item_change_tag(star, MICROGRAPH_PARTICLE_FILE, PARTICLE_FILE);
	item_change_tag(star, MICROGRAPH_FILAMENT_FILE, FILAMENT_FILE);
	
	item_change_tag(star, PARTICLE_NUMBER, PARTICLE_ID);
	item_change_tag(star, PARTICLE_MG_X, PARTICLE_X);
	item_change_tag(star, PARTICLE_MG_Y, PARTICLE_Y);
	item_change_tag(star, PARTICLE_MG_Z, PARTICLE_Z);
	item_change_tag(star, PARTICLE_X_ORIGIN, PARTICLE_ORIGIN_X);
	item_change_tag(star, PARTICLE_Y_ORIGIN, PARTICLE_ORIGIN_Y);
	item_change_tag(star, PARTICLE_Z_ORIGIN, PARTICLE_ORIGIN_Z);
	item_change_tag(star, MICROGRAPH_BOX_RADIUS, PARTICLE_BOX_RADIUS);
	item_change_tag(star, MICROGRAPH_BOX_RADIUS_X, PARTICLE_BOX_RADIUS_X);
	item_change_tag(star, MICROGRAPH_BOX_RADIUS_Y, PARTICLE_BOX_RADIUS_Y);
	item_change_tag(star, MICROGRAPH_BOX_RADIUS_Z, PARTICLE_BOX_RADIUS_Z);
	item_change_tag(star, MICROGRAPH_BAD_RADIUS, PARTICLE_BAD_RADIUS);
	item_change_tag(star, MICROGRAPH_BAD_X, PARTICLE_BAD_X);
	item_change_tag(star, MICROGRAPH_BAD_Y, PARTICLE_BAD_Y);
	item_change_tag(star, MICROGRAPH_BAD_Z, PARTICLE_BAD_Z);

	item_change_tag(star, MICROGRAPH_MARKER_RADIUS, MARKER_RADIUS);
	item_change_tag(star, MICROGRAPH_MARKER_ID, MARKER_ID);
	item_change_tag(star, MICROGRAPH_MARKER_X, MARKER_X);
	item_change_tag(star, MICROGRAPH_MARKER_Y, MARKER_Y);
	item_change_tag(star, MICROGRAPH_MARKER_Z, MARKER_Z);
	item_change_tag(star, MICROGRAPH_MARKER_ERROR_X, MARKER_ERROR_X);
	item_change_tag(star, MICROGRAPH_MARKER_ERROR_Y, MARKER_ERROR_Y);
	item_change_tag(star, MICROGRAPH_MARKER_ERROR_Z, MARKER_ERROR_Z);
	item_change_tag(star, MICROGRAPH_MARKER_FOM, MARKER_FOM);

	item_change_tag(star, MICROGRAPH_VOLTAGE, CTF_VOLTAGE);
	item_change_tag(star, MICROGRAPH_CTF_CS, CTF_CS);
	item_change_tag(star, MICROGRAPH_CTF_CC, CTF_CC);
	item_change_tag(star, MICROGRAPH_CTF_ALPHA, CTF_ALPHA);
	item_change_tag(star, MICROGRAPH_CTF_DE, CTF_DE);
	item_change_tag(star, MICROGRAPH_CTF_AMP_CONT, CTF_AMP);
	item_change_tag(star, MICROGRAPH_CTF_ZERO, CTF_ZERO);
	item_change_tag(star, MICROGRAPH_CTF_DEF_AVG, CTF_DEF_AVG);
	item_change_tag(star, MICROGRAPH_CTF_DEF_DEV, CTF_DEF_DEV);
	item_change_tag(star, MICROGRAPH_CTF_DEF_MIN, CTF_DEF_MIN);
	item_change_tag(star, MICROGRAPH_CTF_DEF_MAX, CTF_DEF_MAX);
	item_change_tag(star, MICROGRAPH_CTF_AST_ANG, CTF_AST_ANG);
	item_change_tag(star, MICROGRAPH_CTF_BASELINE, CTF_BASELINE);
	item_change_tag(star, MICROGRAPH_CTF_ENVELOPE, CTF_ENVELOPE);

	return 0;
}

Bframe*		frame_from_starblock(Bstar_block* block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG frame_from_starblock:" << endl;
	
	Bstar_item*		item;
	Bstring*		data;
	Bframe*			framelist = NULL;
	Bframe*			frame = NULL;
	
	item = item_find(block, MICROGRAPH_FRAME);
	if ( item )
		for ( data=item->data; data; data=data->next ) {
			frame = frame_add(&frame, data->integer());
			if ( !framelist ) framelist = frame;
		}
	if ( ( item = item_find(block, MICROGRAPH_FRAME_SHIFT_X) ) )
		for ( data=item->data, frame=framelist; data && frame; data=data->next, frame=frame->next )
			if ( data->length() > 1 ) frame->shift[0] = data->real();
	if ( ( item = item_find(block, MICROGRAPH_FRAME_SHIFT_Y) ) )
		for ( data=item->data, frame=framelist; data && frame; data=data->next, frame=frame->next )
			if ( data->length() > 1 ) frame->shift[1] = data->real();
	if ( ( item = item_find(block, MICROGRAPH_FRAME_SELECT) ) )
		for ( data=item->data, frame=framelist; data && frame; data=data->next, frame=frame->next )
			if ( data->length() > 1 ) frame->sel = data->integer();
	if ( ( item = item_find(block, MICROGRAPH_FRAME_FOM) ) )
		for ( data=item->data, frame=framelist; data && frame; data=data->next, frame=frame->next )
			if ( data->length() > 1 ) frame->fom = data->real();

	return framelist;
}

Bparticle*	particle_from_starblock(Bstar_block* block, FOMType fom_tag[NFOM])
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG particle_from_starblock:" << endl;
	
	Bstar_item*		item;
	Bstring*		data;
	Bparticle*		partlist = NULL;
	Bparticle*		part = NULL;
	
	Bstring			tag;
	int				f, omega_flag(0);
	Euler			euler;
	
	item = item_find(block, PARTICLE_ID);
	if ( item )
		for ( data=item->data; data; data=data->next ) {
			part = particle_add(&part, data->integer());
			if ( !partlist ) partlist = part;
		}
	if ( ( item = item_find(block, PARTICLE_FILE) ) ) {
		if ( item->loop >= 0 )
			for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
				if ( data->length() > 1 ) part->fpart = *data;
	}
	item = item_find(block, PARTICLE_GROUP);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->group = data->integer();
	item = item_find(block, PARTICLE_DEFOCUS);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->def = data->real();
	item = item_find(block, PARTICLE_DEF_DEV);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->dev = data->real();
	item = item_find(block, PARTICLE_AST_ANG);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->ast = data->real() * M_PI/180.0;
	item = item_find(block, PARTICLE_MAGNIF);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->mag = data->real();
	item = item_find(block, PARTICLE_X);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->loc[0] = data->real();
	item = item_find(block, PARTICLE_Y);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->loc[1] = data->real();
	item = item_find(block, PARTICLE_Z);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->loc[2] = data->real();
	item = item_find(block, PARTICLE_PIXEL);
	if ( item ) {
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->pixel_size[0] = part->pixel_size[1] = part->pixel_size[2] = data->real();
	} else {
		item = item_find(block, PARTICLE_PIXEL_X);
		if ( item )
			for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
				part->pixel_size[0] = data->real();
		item = item_find(block, PARTICLE_PIXEL_Y);
		if ( item )
			for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
				part->pixel_size[1] = data->real();
		item = item_find(block, PARTICLE_PIXEL_Z);
		if ( item )
			for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
				part->pixel_size[2] = data->real();
	}
	item = item_find(block, PARTICLE_ORIGIN_X);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->ori[0] = data->real();
	item = item_find(block, PARTICLE_ORIGIN_Y);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->ori[1] = data->real();
	item = item_find(block, PARTICLE_ORIGIN_Z);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->ori[2] = data->real();
	item = item_find(block, PARTICLE_VIEW_X);
	if ( item ) {
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->view[0] = data->real();
		item = item_find(block, PARTICLE_VIEW_Y);
		if ( item )
			for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
				part->view[1] = data->real();
		item = item_find(block, PARTICLE_VIEW_Z);
		if ( item )
			for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
				part->view[2] = data->real();
		item = item_find(block, PARTICLE_VIEW_ANGLE);
		if ( item )
			for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
				part->view[3] = data->real()*M_PI/180.0;
	} else {
		omega_flag = 0;
		item = item_find(block, PARTICLE_PSI);
		if ( !item ) {
			item = item_find(block, PARTICLE_OMEGA);
			omega_flag = 1;
		}
		if ( item )
			for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next ) {
				part->view[0] = data->real()*M_PI/180.0;
				if ( omega_flag ) part->view[0] = -part->view[0];
			}
		item = item_find(block, PARTICLE_THETA);
		if ( item )
			for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
				part->view[1] = data->real()*M_PI/180.0;
		item = item_find(block, PARTICLE_PHI);
		if ( item )
			for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
				part->view[2] = data->real()*M_PI/180.0;
		for ( part=partlist; part; part=part->next ) {
			euler = Euler(part->view[0], part->view[1], part->view[2]);
			part->view = euler.view();
		}
	}
	for ( f=0; f<NFOM; f++ ) if ( fom_tag[f] ) {
		tag = get_fom_tag(fom_tag[f]);
		if ( tag.length() > 1 ) {
			item = item_find(block, tag);
			if ( item )
				for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
					part->fom[f] = data->real();
		}
	}
	item = item_find(block, PARTICLE_SELECT);
	if ( item )
		for ( data=item->data, part=partlist; data && part; data=data->next, part=part->next )
			part->sel = data->integer();

	return partlist;
}

Bfilament*	filament_from_starblock(Bstar_block* block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG filament_from_starblock:" << endl;
	
	Bstar_item*		item;
	Bstring*		data;
	Bfilament*		fillist = NULL;
	Bfilament*		fil = NULL;
	Bfilnode*		fnode = NULL;

	int 			f, fp;

	item = item_find(block, FILAMENT_ID);
	for ( fp=-1, data=item->data; data; data=data->next ) {
		f = data->integer();
		if ( f != fp ) fil = filament_add(&fillist, f);
		filament_node_add(&fil->node, 0);
		fp = f;
	}
	if ( ( item = item_find(block, FILAMENT_FILE) ) ) {
		if ( item->loop >= 0 )
			for ( data=item->data, fil=fillist; data && fil; fil=fil->next )
				fil->ffil = *data;
	}

	item = item_find(block, FILAMENT_NODE_ID);
	if ( item )
		for ( data=item->data, fil=fillist; data && fil; fil=fil->next )
			for ( fnode=fil->node; data && fnode; data=data->next, fnode=fnode->next )
				fnode->id = data->integer();
	item = item_find(block, FILAMENT_NODE_X);
	if ( item )
		for ( data=item->data, fil=fillist; data && fil; fil=fil->next )
			for ( fnode=fil->node; data && fnode; data=data->next, fnode=fnode->next )
				fnode->loc[0] = data->real();
	item = item_find(block, FILAMENT_NODE_Y);
	if ( item )
		for ( data=item->data, fil=fillist; data && fil; fil=fil->next )
			for ( fnode=fil->node; data && fnode; data=data->next, fnode=fnode->next )
				fnode->loc[1] = data->real();
	item = item_find(block, FILAMENT_NODE_Z);
	if ( item )
		for ( data=item->data, fil=fillist; data && fil; fil=fil->next )
			for ( fnode=fil->node; data && fnode; data=data->next, fnode=fnode->next )
				fnode->loc[2] = data->real();

	return fillist;
}

Bbadarea*	badarea_from_starblock(Bstar_block* block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG badarea_from_starblock:" << endl;
	
	Bstar_item*		item;
	Bstring*		data;
	Bbadarea*		bad = NULL;
	Bbadarea*		badlist = NULL;
	
	long			j, nbad = item_get_number_for_block(block, PARTICLE_BAD_X);
	
	for ( j=0; j<nbad; j++ ) {
		bad = (Bbadarea *) add_item((char **) &bad, sizeof(Bbadarea));
		if ( !badlist ) badlist = bad;
	}
	item = item_find(block, PARTICLE_BAD_X);
	if ( item )
		for ( data=item->data, bad=badlist; data && bad; data=data->next, bad=bad->next )
			bad->loc[0] = data->real();
	item = item_find(block, PARTICLE_BAD_Y);
	if ( item )
		for ( data=item->data, bad=badlist; data && bad; data=data->next, bad=bad->next )
			bad->loc[1] = data->real();
	item = item_find(block, PARTICLE_BAD_Z);
	if ( item )
		for ( data=item->data, bad=badlist; data && bad; data=data->next, bad=bad->next )
			bad->loc[2] = data->real();

	return badlist;
}

Bmarker*	marker_from_starblock(Bstar_block* block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG marker_from_starblock:" << endl;
	
	Bstar_item*		item;
	Bstring*		data;
	Bmarker*		marklist = NULL;
	Bmarker*		mark = NULL;
	
	long			j, nmark = item_get_number_for_block(block, MARKER_X);
	
	for ( j=0; j<nmark; j++ ) {
		mark = (Bmarker *) add_item((char **) &mark, sizeof(Bmarker));
		if ( !marklist ) marklist = mark;
		mark->id = j+1;
		mark->fom = 1;
		mark->sel = 1;
	}
	item = item_find(block, MARKER_ID);
	if ( item )
		for ( data=item->data, mark=marklist; data && mark; data=data->next, mark=mark->next )
			mark->id = data->integer();
	item = item_find(block, MARKER_X);
	if ( item )
		for ( data=item->data, mark=marklist; data && mark; data=data->next, mark=mark->next )
			mark->loc[0] = data->real();
	item = item_find(block, MARKER_Y);
	if ( item )
		for ( data=item->data, mark=marklist; data && mark; data=data->next, mark=mark->next )
			mark->loc[1] = data->real();
	item = item_find(block, MARKER_Z);
	if ( item )
		for ( data=item->data, mark=marklist; data && mark; data=data->next, mark=mark->next )
			mark->loc[2] = data->real();
	item = item_find(block, MARKER_ERROR_X);
	if ( item )
		for ( data=item->data, mark=marklist; data && mark; data=data->next, mark=mark->next )
			mark->err[0] = data->real();
	item = item_find(block, MARKER_ERROR_Y);
	if ( item )
		for ( data=item->data, mark=marklist; data && mark; data=data->next, mark=mark->next )
			mark->err[1] = data->real();
	item = item_find(block, MARKER_ERROR_Z);
	if ( item )
		for ( data=item->data, mark=marklist; data && mark; data=data->next, mark=mark->next )
			mark->err[2] = data->real();
	item = item_find(block, MARKER_RESIDUAL);
	if ( item )
		for ( data=item->data, mark=marklist; data && mark; data=data->next, mark=mark->next )
			mark->res = data->real();
	item = item_find(block, MARKER_FOM);
	if ( item )
		for ( data=item->data, mark=marklist; data && mark; data=data->next, mark=mark->next )
			mark->fom = data->real();
	item = item_find(block, MARKER_SELECT);
	if ( item )
		for ( data=item->data, mark=marklist; data && mark; data=data->next, mark=mark->next )
			mark->sel = data->integer();

	return marklist;
}

Bstrucfac*	strucfac_from_starblock(Bstar_block* block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG strucfac_from_starblock:" << endl;
	
	Bstar_item*		item;
	Bstring*		data;
	Bstrucfac*		sflist = NULL;
	Bstrucfac*		sf = NULL;
	
	long			j, nsf = item_get_number_for_block(block, REFLEX_H);
	
	for ( j=0; j<nsf; j++ ) {
		sf = (Bstrucfac *) add_item((char **) &sf, sizeof(Bstrucfac));
		if ( !sflist ) sflist = sf;
		sf->fom = sf->sel = 1;
	}
	item = item_find(block, REFLEX_H);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->index[0] = data->integer();
	item = item_find(block, REFLEX_K);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->index[1] = data->integer();
	item = item_find(block, REFLEX_L);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->index[2] = data->integer();
	item = item_find(block, REFLEX_X);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->loc[0] = data->real();
	item = item_find(block, REFLEX_Y);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->loc[1] = data->real();
	item = item_find(block, REFLEX_Z);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->loc[2] = data->real();
	item = item_find(block, REFLEX_AMP);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->amp = data->real();
	item = item_find(block, REFLEX_SIGAMP);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->sigamp = data->real();
	item = item_find(block, REFLEX_PHI);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->phi = data->real();
	item = item_find(block, REFLEX_SIGPHI);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->sigphi = data->real();
	item = item_find(block, REFLEX_FOM);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->fom = data->real();
	item = item_find(block, REFLEX_STATUS);
	if ( item )
		for ( data=item->data, sf=sflist; data && sf; data=data->next, sf=sf->next )
			sf->sel = data->integer();

	return sflist;
}

Blayerline*	layerline_from_starblock(Bstar_block* block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG layerline_from_starblock:" << endl;
	
	Bstar_item*		item;
	Bstring*		data;
	Blayerline*		linelist = NULL;
	Blayerline*		line = NULL;
	
	long			j, nline = item_get_number_for_block(block, LAYERLINE_NUMBER);
	
	for ( j=0; j<nline; j++ ) {
		line = (Blayerline *) add_item((char **) &line, sizeof(Blayerline));
		if ( !linelist ) linelist = line;
		line->fom = line->sel = 1;
	}
	item = item_find(block, LAYERLINE_NUMBER);
	if ( item )
		for ( data=item->data, line=linelist; data && line; data=data->next, line=line->next )
			line->number = data->integer();
	item = item_find(block, LAYERLINE_ORDER);
	if ( item )
		for ( data=item->data, line=linelist; data && line; data=data->next, line=line->next )
			line->order = data->integer();
	item = item_find(block, LAYERLINE_DISTANCE);
	if ( item )
		for ( data=item->data, line=linelist; data && line; data=data->next, line=line->next )
			line->distance = data->real();
	item = item_find(block, LAYERLINE_FREQ);
	if ( item )
		for ( data=item->data, line=linelist; data && line; data=data->next, line=line->next )
			line->freq = data->real();
	item = item_find(block, LAYERLINE_AMP);
	if ( item )
		for ( data=item->data, line=linelist; data && line; data=data->next, line=line->next )
			line->amp = data->real();
	item = item_find(block, LAYERLINE_FOM);
	if ( item )
		for ( data=item->data, line=linelist; data && line; data=data->next, line=line->next )
			line->fom = data->real();
	item = item_find(block, LAYERLINE_SELECT);
	if ( item )
		for ( data=item->data, line=linelist; data && line; data=data->next, line=line->next )
			line->sel = data->integer();

	return linelist;
}

CTFparam*		ctf_from_starblock(Bstar_block* block)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ctf_from_starblock:" << endl;
	
	if ( item_get_number_for_block(block, CTF_DEF_AVG) < 1 && 
		item_get_number_for_block(block, CTF_DEF_MIN) < 1 )
			return NULL;
	
	CTFparam*		ctf = new CTFparam;
	
	ctf->volt(item_get_float(block, CTF_VOLTAGE));
	if ( item_get_number_for_block(block, CTF_FOCAL) )
		ctf->focal_length(item_get_float(block, CTF_FOCAL));
	if ( item_get_number_for_block(block, CTF_APERTURE) )
		ctf->objective_aperture(item_get_float(block, CTF_APERTURE));
	if ( item_get_number_for_block(block, CTF_AMP_SHIFT) )
		ctf->amp_shift(item_get_float(block, CTF_AMP_SHIFT));
	else
		ctf->amp_shift(asin(item_get_float(block, CTF_AMP)));
	if ( item_get_number_for_block(block, CTF_DEF_AVG) > 0 ) {
		ctf->defocus_average(item_get_float(block, CTF_DEF_AVG));
		ctf->defocus_deviation(item_get_float(block, CTF_DEF_DEV));
	} else {
		double	min = item_get_float(block, CTF_DEF_MIN);
		double	max = item_get_float(block, CTF_DEF_MAX);
		ctf->defocus_average((max + min)/2);
		ctf->defocus_deviation((max - min)/2);
	}
	ctf->astigmatism_angle(item_get_float(block, CTF_AST_ANG) * M_PI/180.0);
	ctf->Cs(item_get_float(block, CTF_CS));
	ctf->Cc(item_get_float(block, CTF_CC));
	ctf->alpha(item_get_float(block, CTF_ALPHA));
	ctf->dE(item_get_float(block, CTF_DE));
	ctf->zero(1);

	Bstar_item*	item;
	
	if ( ( item = item_find(block, CTF_BASELINE) ) )
		ctf->parse_baseline_equation(*(item->data));
	
	if ( ( item = item_find(block, CTF_ENVELOPE) ) )
		ctf->parse_envelope_equation(*(item->data));

	return ctf;
}

Bmicrograph*	micrograph_from_starblock(Bstar_block* block, Bstring& mg_id, FOMType fom_tag[NFOM])
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG micrograph_from_starblock:" << endl;
	
	Bstar_item*		item;
	
	int 			npart(0), nfil(0), nbad, nmark, nsf, nline;
	
	Bmicrograph*	mg = NULL;
	mg = micrograph_add(&mg, mg_id);
	if ( ( item = item_find(block, MICROGRAPH_FILE) ) ) mg->fmg = *(item->data);
	if ( ( item = item_find(block, MICROGRAPH_FRAMES_FILE) ) ) mg->fframe = *(item->data);
	if ( ( item = item_find(block, PARTICLE_FILE) ) ) {
		if ( item->loop < 0 )
			mg->fpart = *(item->data);
	}
	if ( ( item = item_find(block, FILAMENT_FILE) ) ) mg->ffil = *(item->data);
	if ( ( item = item_find(block, MICROGRAPH_TRANSFORM_FILE) ) ) mg->fft = *(item->data);
	if ( ( item = item_find(block, MICROGRAPH_POWERSPEC_FILE) ) ) mg->fps = *(item->data);
	mg->img_num = item_get_integer(block, MICROGRAPH_NUMBER);
	mg->select = item_get_integer(block, MICROGRAPH_SELECT);
	mg->fom = item_get_float(block, MICROGRAPH_FOM);
	mg->magnification = item_get_float(block, MICROGRAPH_MAGNIFICATION);
	mg->sampling = item_get_float(block, MICROGRAPH_SAMPLING);
	if ( ( item = item_find(block, MICROGRAPH_PIXEL) ) ) {
		mg->pixel_size[0] = mg->pixel_size[1] = item_get_float(block, MICROGRAPH_PIXEL);
	} else {
		mg->pixel_size[0] = item_get_float(block, MICROGRAPH_PIXEL_X);
		mg->pixel_size[1] = item_get_float(block, MICROGRAPH_PIXEL_Y);
	}
	mg->dose = item_get_float(block, MICROGRAPH_DOSE);
	mg->intensity = item_get_float(block, MICROGRAPH_INTENSITY);
	mg->wri = item_get_float(block, MICROGRAPH_WATER_RING);
	mg->origin[0] = item_get_float(block, MICROGRAPH_ORIGIN_X);
	mg->origin[1] = item_get_float(block, MICROGRAPH_ORIGIN_Y);
	mg->origin[2] = item_get_float(block, MICROGRAPH_ORIGIN_Z);
	mg->scale[0] = item_get_float(block, MICROGRAPH_SCALE_X);
	mg->scale[1] = item_get_float(block, MICROGRAPH_SCALE_Y);
	mg->scale[2] = item_get_float(block, MICROGRAPH_SCALE_Z);
	mg->tilt_axis = item_get_float(block, MICROGRAPH_TILT_AXIS);
	mg->tilt_angle = item_get_float(block, MICROGRAPH_TILT_ANGLE);
	mg->level_angle = item_get_float(block, MICROGRAPH_LEVEL_ANGLE);
	mg->rot_angle = item_get_float(block, MICROGRAPH_ROT_ANGLE);
	mg->matrix[0][0] = item_get_float(block, MICROGRAPH_MATRIX_1_1);
	mg->matrix[0][1] = item_get_float(block, MICROGRAPH_MATRIX_1_2);
	mg->matrix[0][2] = item_get_float(block, MICROGRAPH_MATRIX_1_3);
	mg->matrix[1][0] = item_get_float(block, MICROGRAPH_MATRIX_2_1);
	mg->matrix[1][1] = item_get_float(block, MICROGRAPH_MATRIX_2_2);
	mg->matrix[1][2] = item_get_float(block, MICROGRAPH_MATRIX_2_3);
	mg->matrix[2][0] = item_get_float(block, MICROGRAPH_MATRIX_3_1);
	mg->matrix[2][1] = item_get_float(block, MICROGRAPH_MATRIX_3_2);
	mg->matrix[2][2] = item_get_float(block, MICROGRAPH_MATRIX_3_3);
	mg->hvec[0] = item_get_float(block, MICROGRAPH_HVEC_X);
	mg->hvec[1] = item_get_float(block, MICROGRAPH_HVEC_Y);
	mg->hvec[2] = item_get_float(block, MICROGRAPH_HVEC_Z);
	mg->kvec[0] = item_get_float(block, MICROGRAPH_KVEC_X);
	mg->kvec[1] = item_get_float(block, MICROGRAPH_KVEC_Y);
	mg->kvec[2] = item_get_float(block, MICROGRAPH_KVEC_Z);
	mg->lvec[0] = item_get_float(block, MICROGRAPH_LVEC_X);
	mg->lvec[1] = item_get_float(block, MICROGRAPH_LVEC_Y);
	mg->lvec[2] = item_get_float(block, MICROGRAPH_LVEC_Z);
	mg->helix_axis = item_get_float(block, MICROGRAPH_HELIX_AXIS);
	mg->helix_rise = item_get_float(block, MICROGRAPH_HELIX_RISE);
	mg->helix_angle = item_get_float(block, MICROGRAPH_HELIX_ANGLE);
	mg->helix_radius = item_get_float(block, MICROGRAPH_HELIX_RADIUS);
	mg->tilt_axis *= M_PI/180.0;
	mg->tilt_angle *= M_PI/180.0;
	mg->level_angle *= M_PI/180.0;
	mg->rot_angle  *= M_PI/180.0;
	mg->helix_axis  *= M_PI/180.0;
	mg->helix_angle  *= M_PI/180.0;

	if ( item_get_number_for_block(block, PARTICLE_BOX_RADIUS) )
		mg->box_size[0] = mg->box_size[1] = mg->box_size[2] = 
			2 * item_get_integer(block, PARTICLE_BOX_RADIUS);
	if ( item_get_number_for_block(block, PARTICLE_BOX_RADIUS_X) ) {
		mg->box_size[0] = 2 * item_get_integer(block, PARTICLE_BOX_RADIUS_X);
		mg->box_size[1] = 2 * item_get_integer(block, PARTICLE_BOX_RADIUS_Y);
		mg->box_size[2] = 2 * item_get_integer(block, PARTICLE_BOX_RADIUS_Z);
	}
	if ( item_get_number_for_block(block, PARTICLE_BOX_SIZE) )
		mg->box_size[0] = mg->box_size[1] = mg->box_size[2] =
			item_get_integer(block, PARTICLE_BOX_SIZE);
	if ( item_get_number_for_block(block, PARTICLE_BOX_SIZE_X) ) {
		mg->box_size[0] = item_get_integer(block, PARTICLE_BOX_SIZE_X);
		mg->box_size[1] = item_get_integer(block, PARTICLE_BOX_SIZE_Y);
		mg->box_size[2] = item_get_integer(block, PARTICLE_BOX_SIZE_Z);
	}

	mg->filament_width = item_get_float(block, FILAMENT_WIDTH);
	mg->fil_node_radius = item_get_float(block, FILNODE_RADIUS);
	mg->bad_radius = item_get_float(block, PARTICLE_BAD_RADIUS);
	mg->sf_radius = item_get_float(block, REFLEX_RADIUS);
	mg->mark_radius = item_get_float(block, MARKER_RADIUS);
	
	mg->ctf = ctf_from_starblock(block);
	
	mg->frame = frame_from_starblock(block);
	
	npart = item_get_number_for_block(block, PARTICLE_ID);
	if ( verbose & VERB_FULL )
		cout << npart << " particles" << endl;
	if ( npart )
		mg->part = particle_from_starblock(block, fom_tag);

	nfil = item_get_number_for_block(block, FILAMENT_ID);
	if ( nfil ) {
		if ( verbose & VERB_FULL )
			cout << nfil << " filament nodes" << endl;
		mg->fil = filament_from_starblock(block);
	}
			
	nbad = item_get_number_for_block(block, PARTICLE_BAD_X);
	if ( nbad ) 
		mg->bad = badarea_from_starblock(block);

	nmark = item_get_number_for_block(block, MARKER_X);
	if ( nmark )
		mg->mark = marker_from_starblock(block);

	nsf = item_get_number_for_block(block, REFLEX_H);
	if ( nsf ) 
		mg->sf = strucfac_from_starblock(block);

	nline = item_get_number_for_block(block, LAYERLINE_NUMBER);
	if ( nline ) 
		mg->layer = layerline_from_starblock(block);

	return mg;
}

Breconstruction*	reconstruction_from_starblock(Bstar_block* block, FOMType fom_tag[NFOM])
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock:" << endl;
	
	if ( verbose & VERB_FULL )
		cout << "Converting reconstruction" << endl;
	
	Bstar_item*			item;
	Bstring				id("1");
	Breconstruction*	rec = NULL;
	reconstruction_add(&rec, id);
	
	rec->block = block->number;
	rec->select = 1;
	
	if ( ( item = item_find(block, MAP_ID) ) ) rec->id = *(item->data);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: ID = " << rec->id << endl;

	if ( ( item = item_find(block, MAP_RECONSTRUCTION) ) ) rec->frec = *(item->data);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: frec = " << rec->frec << endl;
	
	if ( ( item = item_find(block, MAP_TRANSFORM_FILE) ) ) rec->fft = *(item->data);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: fft = " << rec->fft << endl;
	
	if ( ( item = item_find(block, MAP_POWERSPEC_FILE) ) ) rec->fps = *(item->data);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: fps = " << rec->fps << endl;
	
	if ( ( item = item_find(block, PARTICLE_FILE) ) ) {
		if ( item->loop < 0 )
			rec->fpart = *(item->data);
	}
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: fpart = " << rec->fpart << endl;
	
	if ( ( item = item_find(block, FILAMENT_FILE) ) ) rec->ffil = *(item->data);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: fpart = " << rec->ffil << endl;
	
	if ( ( item = item_find(block, MAP_MODEL) ) ) string_add(&rec->model, *(item->data));
	
	rec->select = item_get_integer(block, MAP_SELECT);
	rec->fom = item_get_float(block, MAP_FOM);

	rec->origin[0] = item_get_float(block, MAP_ORIGIN_X);
	rec->origin[1] = item_get_float(block, MAP_ORIGIN_Y);
	rec->origin[2] = item_get_float(block, MAP_ORIGIN_Z);
	rec->scale[0] = item_get_float(block, MAP_SCALE_X);
	rec->scale[1] = item_get_float(block, MAP_SCALE_Y);
	rec->scale[2] = item_get_float(block, MAP_SCALE_Z);
	if ( ( item = item_find(block, MAP_VOXEL_SIZE) ) ) {
		rec->voxel_size[0] = item_get_float(block, MAP_VOXEL_SIZE);
		rec->voxel_size[1] = rec->voxel_size[2] = rec->voxel_size[0];
	} else {
		rec->voxel_size[0] = item_get_float(block, MAP_VOXEL_SIZE_X);
		rec->voxel_size[1] = item_get_float(block, MAP_VOXEL_SIZE_Y);
		rec->voxel_size[2] = item_get_float(block, MAP_VOXEL_SIZE_Z);
	}

	if ( ( item = item_find(block, MAP_SYMMETRY) ) ) rec->symmetry = *(item->data);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_from_starblock: symmetry = " << rec->symmetry << endl;

	if ( item_find(block, CTF_VOLTAGE) ) rec->ctf = ctf_from_starblock(block);
	
	if ( item_get_number_for_block(block, PARTICLE_BOX_RADIUS) )
		rec->box_size[0] = rec->box_size[1] = rec->box_size[2] = 
			2 * item_get_integer(block, PARTICLE_BOX_RADIUS);
	if ( item_get_number_for_block(block, PARTICLE_BOX_RADIUS_X) ) {
		rec->box_size[0] = 2 * item_get_integer(block, PARTICLE_BOX_RADIUS_X);
		rec->box_size[1] = 2 * item_get_integer(block, PARTICLE_BOX_RADIUS_Y);
		rec->box_size[2] = 2 * item_get_integer(block, PARTICLE_BOX_RADIUS_Z);
	}
	if ( item_get_number_for_block(block, PARTICLE_BOX_SIZE) )
		rec->box_size[0] = rec->box_size[1] = rec->box_size[2] =
			item_get_integer(block, PARTICLE_BOX_SIZE);
	if ( item_get_number_for_block(block, PARTICLE_BOX_SIZE_X) ) {
		rec->box_size[0] = item_get_integer(block, PARTICLE_BOX_SIZE_X);
		rec->box_size[1] = item_get_integer(block, PARTICLE_BOX_SIZE_Y);
		rec->box_size[2] = item_get_integer(block, PARTICLE_BOX_SIZE_Z);
	}

	rec->filament_width = item_get_float(block, FILAMENT_WIDTH);
	rec->fil_node_radius = item_get_float(block, FILNODE_RADIUS);
	rec->bad_radius = item_get_float(block, PARTICLE_BAD_RADIUS);
	rec->sf_radius = item_get_float(block, REFLEX_RADIUS);
	rec->mark_radius = item_get_float(block, MARKER_RADIUS);

	rec->view[0] = item_get_float(block, MAP_VIEW_X);
	rec->view[1] = item_get_float(block, MAP_VIEW_Y);
	rec->view[2] = item_get_float(block, MAP_VIEW_Z);
	rec->view[3] = item_get_float(block, MAP_VIEW_ANGLE)*M_PI/180.0;

	long		npart = item_get_number_for_block(block, PARTICLE_ID);
	if ( npart ) {
		if ( verbose & VERB_FULL )
			cout << "Reading reconstruction \"" << rec->id << "\" with " << npart << " particles" << endl;
		rec->part = particle_from_starblock(block, fom_tag);
	}
	
	long		nfil = item_get_number_for_block(block, FILAMENT_ID);
	if ( nfil ) {
		if ( verbose & VERB_FULL )
			cout << "Reading reconstruction \"" << rec->id << "\", with " << nfil << " filament nodes" << endl;
		rec->fil = filament_from_starblock(block);
	}
			
	long		nbad = item_get_number_for_block(block, PARTICLE_BAD_X);
	if ( nbad ) 
		rec->bad = badarea_from_starblock(block);

	long		nmark = item_get_number_for_block(block, MARKER_X);
	if ( nmark ) {
		if ( verbose & VERB_FULL )
			cout << "Reading reconstruction \"" << rec->id << "\" with " << nmark << " markers" << endl;
		rec->mark = marker_from_starblock(block);
	}
	
	long		nsf = item_get_number_for_block(block, REFLEX_H);
	if ( nsf ) {
		if ( verbose & VERB_FULL )
			cout << "Reading reconstruction \"" << rec->id << "\" with " << nsf << " structure factors" << endl;
		rec->sf = strucfac_from_starblock(block);
	}
	
	return rec;
}

/*
@brief 	Converts micrograph data from a STAR data base to a project structure.
@param	*star		STAR data base.
@param 	*project	project structure.
@return int			error code (<0 means failure).

	The function sets up the project hierarchy from a STAR database file.
	Micrograph data blocks must contain the "micrograph.id" tag.
	If the "micrograph.field_of_view_id" is present, the micrographs are
	distributed into multi-micrograph fields-of-view. Otherwise, each
	micrograph is assumed to be a unique field-of-view.
	All angles given in the STAR data base are assumed to be in degrees, 
	and are converted here only to radians.

**/
int			star_to_project(Bstar* star, Bproject* project)
{
	int				err(0);
	int 			nmg(0), nrec(0);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG star_to_project:" << endl;
	
	mg_star_update_tags(star);
		
	nmg = item_get_number(star, MICROGRAPH_ID);
	if ( nmg < 1 ) nmg = item_get_number(star, MICROGRAPH_FILE);
	if ( nmg < 1 ) nmg = item_get_number(star, MICROGRAPH_FRAMES_FILE);
	if ( nmg < 1 ) nmg = item_get_number(star, MICROGRAPH_TRANSFORM_FILE);
	if ( nmg < 1 ) nmg = item_get_number(star, MICROGRAPH_POWERSPEC_FILE);

	nrec = item_get_number(star, MAP_ID);
	if ( nrec < 1 ) nrec = item_get_number(star, MAP_RECONSTRUCTION);
	if ( nrec < 1 ) nrec = item_get_number(star, MAP_TRANSFORM_FILE);
	if ( nrec < 1 ) nrec = item_get_number(star, MAP_POWERSPEC_FILE);
	
	if ( nmg < 1 && nrec < 1 ) {
		cerr << "Warning: No micrograph or reconstruction data blocks found" << endl;
		return 0;
	}

	Bstar_block*	block;
	Bstar_item*		item;
	Bstring*		data;
	
	Bfield* 		field = NULL;
	Bmicrograph		*mg, *mg2;
	Breconstruction	*rec, *rec2;
	
	for ( rec2 = project->rec; rec2 && rec2->next; rec2 = rec2->next ) ;

	block = block_find_with_tag(star, MAP_REFERENCE);
	if ( block ) {
		item = item_find(block, MAP_REFERENCE);
		if ( item )
			for ( data=item->data; data; data=data->next )
				string_add(&project->reference, data->c_str());
	}
	
	if ( ( verbose & VERB_DEBUG ) && project->reference )
		cout << "DEBUG star_to_project: reference map = " << *(project->reference) << endl;
	
	int 			i, j, f, nfield(0);
	Bstring			field_id, mg_id;
	Bstring			base;
	
	// Make sure FOM tags are properly set up
	if ( item_get_number(star, PARTICLE_FOM) ) project->fom_tag[0] = FOM;
	if ( item_get_number(star, PARTICLE_FOM_CV) ) project->fom_tag[1] = FOM_CV;
	if ( item_get_number(star, PARTICLE_HANDA_FOM) ) project->fom_tag[1] = FOM_HAND_A;
	if ( item_get_number(star, PARTICLE_HANDB_FOM) ) project->fom_tag[2] = FOM_HAND_B;
	
	for ( f=i=0; f<NFOM; f++ ) i += project->fom_tag[f];
	if ( i < 2 ) project->fom_tag[0] = FOM;

	// The orientation convention used is set based on whether the view angle
	// array or the phi euler angle array is present in the file.
	// If both are present, the one with the larger number of values wins,
	// otherwise only the view vector and rotation angle
	// representation will be used.
	i = item_get_number(star, PARTICLE_PHI);
	j = item_get_number(star, PARTICLE_VIEW_ANGLE);
	project->euler_flag = 0;
	project->omega_flag = 0;
	if ( i > 0 && j > 0 )
		cerr << "Error: The STAR file contains both views and Euler angles!" << endl;
	else if ( i > j ) {
		project->euler_flag = 1;
		if ( item_get_number(star, PARTICLE_OMEGA) > 0 )
			project->omega_flag = 1;
	}
	
	nmg = 0;
	for ( i=0, block = star->block; block; block = block->next, i++ ) {
		if ( item_get_number_for_block(block, MICROGRAPH_ID) ) {
			mg_id = 0;
			field_id = 0;
			/* This deals with old style access to micrograph data blocks */
			if ( ( item = item_find(block, MICROGRAPH_ID) ) ) {
				mg_id = *(item->data);
			} else {
				if ( ( item = item_find(block, MICROGRAPH_FILE) ) ) mg_id = *(item->data);
				else if ( ( item = item_find(block, MICROGRAPH_FRAMES_FILE) ) ) mg_id = *(item->data);
				else if ( ( item = item_find(block, PARTICLE_FILE) ) ) mg_id = *(item->data);
				else if ( ( item = item_find(block, FILAMENT_FILE) ) ) mg_id = *(item->data);
				else if ( ( item = item_find(block, MICROGRAPH_TRANSFORM_FILE) ) ) mg_id = *(item->data);
				else if ( ( item = item_find(block, MICROGRAPH_POWERSPEC_FILE) ) ) mg_id = *(item->data);
				if ( mg_id.length() > 0 ) mg_id = mg_id.base();
			}
			if ( mg_id.length() > 0 ) {
				if ( ( item = item_find(block, MICROGRAPH_FIELD_ID) ) ) field_id = *(item->data);
				if ( field_id.length() < 1 ) field_id = mg_id;
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
				mg = micrograph_from_starblock(block, mg_id, project->fom_tag);
				mg->block = i;
				if ( verbose & VERB_FULL )
					cout << "Reading field \"" << field->id << "\", micrograph \"" << mg->id <<"\"" << endl;
				if ( field->mg ) {
					for ( mg2 = field->mg; mg2->next; mg2 = mg2->next ) ;
					mg2->next = mg;
				} else field->mg = mg;
				mg_id = 0;
				field_id = 0;
				nmg++;
			}
		} else if ( item_get_number_for_block(block, MAP_ID) ||
				 item_get_number_for_block(block, MAP_RECONSTRUCTION) ||
				 item_get_number_for_block(block, MAP_TRANSFORM_FILE) ||
				 item_get_number_for_block(block, MAP_POWERSPEC_FILE) ) {
			rec = reconstruction_from_starblock(block, project->fom_tag);
			if ( !project->rec ) project->rec = rec;
			else rec2->next = rec;
			rec2 = rec;
		} else if ( item_get_number_for_block(block, PARTICLE_ID) ) {
			project->class_avg = particle_from_starblock(block, project->fom_tag);
		}
	}
	
	item_change_tag(star, CTF_DEF_MIN, CTF_DEF_AVG);
	item_change_tag(star, CTF_DEF_MAX, CTF_DEF_DEV);
			
	return err;
}

int 		item_put_filament_list(Bstar_block* block, const char* tag, Bfilament* fil, const char* format)
{
	if ( !block ) return -1;

	Bstring			s;

	Bstar_item* 	item = item_find_or_make(block, tag);
	Bstring*		data = NULL;
	
	if ( item->data ) {
		string_kill(item->data);
		item->data = NULL;
	}	

	item->type = NumberItem;
	item->loop = -1;

	if ( strcmp(tag, FILAMENT_ID) == 0 )
		item->type = StringItem;

	Bfilnode*		fnode = NULL;

	for ( ; fil; fil=fil->next )
		for ( fnode=fil->node; fnode; fnode=fnode->next ) {
			if ( strcmp(tag, FILAMENT_ID) == 0 ) s = Bstring(fil->id, format);
			else if ( strcmp(tag, FILAMENT_NODE_ID) == 0 ) s = Bstring(fnode->id, format);
			else if ( strcmp(tag, FILAMENT_NODE_X) == 0 ) s = Bstring(fnode->loc[0], format);
			else if ( strcmp(tag, FILAMENT_NODE_Y) == 0 ) s = Bstring(fnode->loc[1], format);
			else if ( strcmp(tag, FILAMENT_NODE_Z) == 0 ) s = Bstring(fnode->loc[2], format);
			data = string_add(&data, s);
			if ( !item->data ) item->data = data;
			if ( item->maxlen < s.length() ) item->maxlen = s.length();
		}
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_put_filament_list: " << tag << "=" << *(item->data) << endl;
	
	return 0;
}

int		frame_to_starblock(Bframe* frame, Bstar_block* block)
{
	long			loop;
	
	item_put_list(block, MICROGRAPH_FRAME, (char *) frame, 
		(char *)&frame->id - (char *)frame, "%4d");
	item_put_list(block, MICROGRAPH_FRAME_SHIFT_X, (char *) frame,
		(char *)&frame->shift[0] - (char *)frame, "%7.3lf");
	item_put_list(block, MICROGRAPH_FRAME_SHIFT_Y, (char *) frame,
		(char *)&frame->shift[1] - (char *)frame, "%7.3lf");
	item_put_list(block, MICROGRAPH_FRAME_FOM, (char *) frame,
		(char *)&frame->fom - (char *)frame, "%7.4lf");
	item_put_list(block, MICROGRAPH_FRAME_SELECT, (char *) frame,
		(char *)&frame->sel - (char *)frame, "%4d");
	loop = item_index(block, MICROGRAPH_FRAME);
	loop_set_identifier(block, loop, 5, MICROGRAPH_FRAME, MICROGRAPH_FRAME_SHIFT_X,
		MICROGRAPH_FRAME_SHIFT_Y, MICROGRAPH_FRAME_FOM, MICROGRAPH_FRAME_SELECT);
	
	return 0;
}

int		particle_to_starblock(Bparticle* part, Bstar_block* block, 
				FOMType fom_tag[NFOM], int euler_flag, int omega_flag)
{
	Bparticle*		p;
	Bstar_item*		item;
	Bstring*		data = NULL;
	
	long			j, f, npart, loop;
	Euler			euler;
	Bstring			tag;
	float*			psi = NULL;
	float*			theta = NULL;
	float*			phi = NULL;
	
	item_put_list(block, PARTICLE_ID, (char *) part, 
		(char *)&part->id - (char *)part, "%4d");
	item_put_list(block, PARTICLE_GROUP, (char *) part, 
		(char *)&part->group - (char *)part, "%4d");
	if ( part->def ) {
		item_put_list(block, PARTICLE_DEFOCUS, (char *) part,
			(char *)&part->def - (char *)part, "%7.0lf");
		item_put_list(block, PARTICLE_DEF_DEV, (char *) part,
			(char *)&part->dev - (char *)part, "%7.0lf");
		item_put_angle_list(block, PARTICLE_AST_ANG, (char *) part,
			(char *)&part->ast - (char *)part, "%7.2lf");
	}
	item_put_list(block, PARTICLE_MAGNIF, (char *) part,
		(char *)&part->mag - (char *)part, "%7.4lf");
	item_put_list(block, PARTICLE_X, (char *) part, 
		(char *)&part->loc[0] - (char *)part, "%7.2lf");
	item_put_list(block, PARTICLE_Y, (char *) part, 
		(char *)&part->loc[1] - (char *)part, "%7.2lf");
	item_put_list(block, PARTICLE_Z, (char *) part, 
		(char *)&part->loc[2] - (char *)part, "%7.2lf");
	item_put_list(block, PARTICLE_PIXEL_X, (char *) part,
		(char *)&part->pixel_size[0] - (char *)part, "%7.4lf");
	item_put_list(block, PARTICLE_PIXEL_Y, (char *) part,
		(char *)&part->pixel_size[1] - (char *)part, "%7.4lf");
	item_put_list(block, PARTICLE_PIXEL_Z, (char *) part,
		(char *)&part->pixel_size[2] - (char *)part, "%7.4lf");
	item_put_list(block, PARTICLE_ORIGIN_X, (char *) part,
		(char *)&part->ori[0] - (char *)part, "%7.3lf");
	item_put_list(block, PARTICLE_ORIGIN_Y, (char *) part, 
		(char *)&part->ori[1] - (char *)part, "%7.3lf");
	item_put_list(block, PARTICLE_ORIGIN_Z, (char *) part, 
		(char *)&part->ori[2] - (char *)part, "%7.3lf");
	if ( euler_flag < 1 ) {
		item_put_list(block, PARTICLE_VIEW_X, (char *) part, 
			(char *)&part->view[0] - (char *)part, "%7.4lf");
		item_put_list(block, PARTICLE_VIEW_Y, (char *) part, 
			(char *)&part->view[1] - (char *)part, "%7.4lf");
		item_put_list(block, PARTICLE_VIEW_Z, (char *) part, 
			(char *)&part->view[2] - (char *)part, "%7.4lf");
		item_put_angle_list(block, PARTICLE_VIEW_ANGLE, (char *) part, 
			(char *)&part->view[3] - (char *)part, "%7.2lf");
	} else {
		npart = particle_count(part);
		psi = new float[npart];
		theta = new float[npart];
		phi = new float[npart];
		for ( j=0, p=part; p; p=p->next, j++ ) {
			euler = Euler(p->view);
			phi[j] = euler.phi()*180.0/M_PI;
			theta[j] = euler.theta()*180.0/M_PI;
			psi[j] = euler.psi()*180.0/M_PI;
			if ( omega_flag ) psi[j] = -euler.psi()*180.0/M_PI;
		}
		item_put_float_array(block, PARTICLE_PHI, npart, phi, "%7.2f");
		item_put_float_array(block, PARTICLE_THETA, npart, theta, "%7.2f");
		if ( omega_flag ) {
			item_put_float_array(block, PARTICLE_OMEGA, npart, psi, "%7.2f");
			item_delete_from_block(block, PARTICLE_PSI);
		} else {
			item_put_float_array(block, PARTICLE_PSI, npart, psi, "%7.2f");
			item_delete_from_block(block, PARTICLE_OMEGA);
		}
		if ( psi ) delete[] psi;
		if ( theta ) delete[] theta;
		if ( phi ) delete[] phi;
	}
	for ( f=0; f<NFOM; f++ ) if ( fom_tag[f] ) {
		tag = get_fom_tag(fom_tag[f]);
		if ( tag.length() > 1 )
			item_put_list(block, tag.c_str(), (char *) part, 
				(char *)&part->fom[f] - (char *)part, "%7.4lf");
	}
	item_put_list(block, PARTICLE_SELECT, (char *) part, 
		(char *)&part->sel - (char *)part, "%4d");
	if ( part->fpart.length() ) {
		item = item_find_or_make(block, PARTICLE_FILE);
		if ( item->data ) {
			string_kill(item->data);
			item->data = NULL;
		}
		for ( p=part; p; p=p->next ) {
			if ( p->fpart.length() < 1 ) data = string_add(&data, "?");
			else data = string_add(&data, p->fpart);
			if ( !item->data ) item->data = data;
			if ( item->maxlen < p->fpart.length() ) item->maxlen = p->fpart.length();
		}
	}
	loop = item_index(block, PARTICLE_ID);
	loop_set_identifier(block, loop, 33, PARTICLE_ID, PARTICLE_GROUP,
		PARTICLE_DEFOCUS, PARTICLE_DEF_DEV, PARTICLE_AST_ANG, PARTICLE_MAGNIF,
		PARTICLE_X, PARTICLE_Y, PARTICLE_Z,
		PARTICLE_PIXEL_X, PARTICLE_PIXEL_Y, PARTICLE_PIXEL_Z,
		PARTICLE_ORIGIN_X, PARTICLE_ORIGIN_Y, PARTICLE_ORIGIN_Z,
		PARTICLE_VIEW_X, PARTICLE_VIEW_Y, PARTICLE_VIEW_Z, PARTICLE_VIEW_ANGLE,
		PARTICLE_PHI, PARTICLE_THETA, PARTICLE_PSI, PARTICLE_OMEGA,
		PARTICLE_SELECT, PARTICLE_FOM, PARTICLE_CC, PARTICLE_FOM_CV, PARTICLE_FOM_SNR,
		PARTICLE_HANDA_FOM, PARTICLE_HANDB_FOM, PARTICLE_PFT_CC,
		PARTICLE_PRJ_CC, PARTICLE_CMP_CC, PARTICLE_RFACTORAB);
	if ( part->fpart.length() ) loop_set_identifier(block, loop, 1, PARTICLE_FILE);
	
	return 0;
}

int		badarea_to_starblock(Bbadarea* bad, Bstar_block* block)
{
	long			loop;
	
	item_put_list(block, PARTICLE_BAD_X, (char *) bad, 
		(char *)&bad->loc[0] - (char *)bad, "%7.2lf");
	item_put_list(block, PARTICLE_BAD_Y, (char *) bad, 
		(char *)&bad->loc[1] - (char *)bad, "%7.2lf");
	item_put_list(block, PARTICLE_BAD_Z, (char *) bad, 
		(char *)&bad->loc[2] - (char *)bad, "%7.2lf");
	loop = item_index(block, PARTICLE_BAD_X);
	loop_set_identifier(block, loop, 3, PARTICLE_BAD_X, PARTICLE_BAD_Y, PARTICLE_BAD_Z);
	
	return 0;
}

int		filament_to_starblock(Bfilament* fil, Bstar_block* block)
{
	long			loop;
	
	item_put_filament_list(block, FILAMENT_ID, fil, "%6d");
	item_put_filament_list(block, FILAMENT_NODE_ID, fil, "%6d");
	item_put_filament_list(block, FILAMENT_NODE_X, fil, "%7.2lf");
	item_put_filament_list(block, FILAMENT_NODE_Y, fil, "%7.2lf");
	item_put_filament_list(block, FILAMENT_NODE_Z, fil, "%7.2lf");
	loop = item_index(block, FILAMENT_ID);
	loop_set_identifier(block, loop, 5, FILAMENT_ID, 
		FILAMENT_NODE_ID, FILAMENT_NODE_X, FILAMENT_NODE_Y, FILAMENT_NODE_Z);
	
	return 0;
}

int		marker_to_starblock(Bmarker* mark, Bstar_block* block)
{
	long			loop;
	
	item_put_list(block, MARKER_ID, (char *) mark, 
		(char *)&mark->id - (char *)mark, "%5d");
	item_put_list(block, MARKER_X, (char *) mark, 
		(char *)&mark->loc[0] - (char *)mark, "%7.2f");
	item_put_list(block, MARKER_Y, (char *) mark, 
		(char *)&mark->loc[1] - (char *)mark, "%7.2f");
	item_put_list(block, MARKER_Z, (char *) mark, 
		(char *)&mark->loc[2] - (char *)mark, "%7.2f");
	item_put_list(block, MARKER_ERROR_X, (char *) mark, 
		(char *)&mark->err[0] - (char *)mark, "%7.2f");
	item_put_list(block, MARKER_ERROR_Y, (char *) mark, 
		(char *)&mark->err[1] - (char *)mark, "%7.2f");
	item_put_list(block, MARKER_ERROR_Z, (char *) mark, 
		(char *)&mark->err[2] - (char *)mark, "%7.2f");
	item_put_list(block, MARKER_RESIDUAL, (char *) mark, 
		(char *)&mark->res - (char *)mark, "%7.4f");
	item_put_list(block, MARKER_FOM, (char *) mark, 
		(char *)&mark->fom - (char *)mark, "%7.4f");
	item_put_list(block, MARKER_SELECT, (char *) mark, 
		(char *)&mark->sel - (char *)mark, "%5d");
	loop = item_index(block, MARKER_X);
	loop_set_identifier(block, loop, 10, MARKER_ID,
		MARKER_X, MARKER_Y, MARKER_Z,
		MARKER_ERROR_X, MARKER_ERROR_Y, MARKER_ERROR_Z,
		MARKER_RESIDUAL, MARKER_FOM, MARKER_SELECT);
	
	return 0;
}

int		strucfac_to_starblock(Bstrucfac* sf, Bstar_block* block)
{
	long			loop;
	
	item_put_list(block, REFLEX_H, (char *) sf, 
		(char *)&sf->index[0] - (char *)sf, "%4d");
	item_put_list(block, REFLEX_K, (char *) sf, 
		(char *)&sf->index[1] - (char *)sf, "%4d");
	item_put_list(block, REFLEX_L, (char *) sf, 
		(char *)&sf->index[2] - (char *)sf, "%4d");
	item_put_list(block, REFLEX_X, (char *) sf, 
		(char *)&sf->loc[0] - (char *)sf, "%7.2lf");
	item_put_list(block, REFLEX_Y, (char *) sf, 
		(char *)&sf->loc[1] - (char *)sf, "%7.2lf");
	item_put_list(block, REFLEX_Z, (char *) sf, 
		(char *)&sf->loc[2] - (char *)sf, "%7.2lf");
	item_put_list(block, REFLEX_AMP, (char *) sf, 
		(char *)&sf->amp - (char *)sf, "%7.2lf");
	item_put_list(block, REFLEX_SIGAMP, (char *) sf, 
		(char *)&sf->sigamp - (char *)sf, "%7.2lf");
	item_put_list(block, REFLEX_PHI, (char *) sf, 
		(char *)&sf->phi - (char *)sf, "%7.2lf");
	item_put_list(block, REFLEX_SIGPHI, (char *) sf, 
		(char *)&sf->sigphi - (char *)sf, "%7.2lf");
	item_put_list(block, REFLEX_FOM, (char *) sf, 
		(char *)&sf->fom - (char *)sf, "%7.4lf");
	item_put_list(block, REFLEX_STATUS, (char *) sf, 
		(char *)&sf->sel - (char *)sf, "%4d");
	loop = item_index(block, REFLEX_H);
	loop_set_identifier(block, loop, 12, REFLEX_H, REFLEX_K, REFLEX_L,
		REFLEX_X, REFLEX_Y, REFLEX_Z,
		REFLEX_AMP, REFLEX_SIGAMP, REFLEX_PHI, REFLEX_SIGPHI,
		REFLEX_FOM, REFLEX_STATUS);
	
	return 0;
}

int		layerline_to_starblock(Blayerline* line, Bstar_block* block)
{
	long			loop;
	
	item_put_list(block, LAYERLINE_NUMBER, (char *) line, 
		(char *)&line->number - (char *)line, "%4d");
	item_put_list(block, LAYERLINE_ORDER, (char *) line, 
		(char *)&line->order - (char *)line, "%4d");
	item_put_list(block, LAYERLINE_DISTANCE, (char *) line, 
		(char *)&line->distance - (char *)line, "%7.2lf");
	item_put_list(block, LAYERLINE_FREQ, (char *) line, 
		(char *)&line->freq - (char *)line, "%7.4lf");
	item_put_list(block, LAYERLINE_AMP, (char *) line, 
		(char *)&line->amp - (char *)line, "%7.2lf");
	item_put_list(block, LAYERLINE_FOM, (char *) line, 
		(char *)&line->fom - (char *)line, "%7.4lf");
	item_put_list(block, LAYERLINE_SELECT, (char *) line, 
		(char *)&line->sel - (char *)line, "%4d");
	loop = item_index(block, LAYERLINE_NUMBER);
	loop_set_identifier(block, loop, 7, LAYERLINE_NUMBER, LAYERLINE_ORDER, 
		LAYERLINE_DISTANCE, LAYERLINE_FREQ, LAYERLINE_AMP, LAYERLINE_FOM, LAYERLINE_SELECT);
	
	return 0;
}

int		ctf_to_starblock(CTFparam* ctf, Bstar_block* block)
{
	float		angle;
	Bstring		string;
	
	item_put_float(block, CTF_VOLTAGE, ctf->volt(), "%f");
	item_put_float(block, CTF_FOCAL, ctf->focal_length(), "%f");
	item_put_float(block, CTF_APERTURE, ctf->objective_aperture(), "%f");
	item_put_float(block, CTF_CS, ctf->Cs(), "%f");
	item_put_float(block, CTF_CC, ctf->Cc(), "%f");
	item_put_float(block, CTF_ALPHA, ctf->alpha(), "%f");
	item_put_float(block, CTF_DE, ctf->dE(), "%f");
	item_put_float(block, CTF_AMP_SHIFT, ctf->amp_shift(), "%f");
	item_put_float(block, CTF_DEF_AVG, ctf->defocus_average(), "%f");
	item_put_float(block, CTF_DEF_DEV, ctf->defocus_deviation(), "%f");
	angle = ctf->astigmatism_angle()*180.0/M_PI;
	item_put_float(block, CTF_AST_ANG, angle, "%f");
	item_put_float(block, CTF_ZERO, ctf->zero(1), "%f");
	if ( ctf->baseline_type() ) {
		string = ctf->baseline_equation();
		item_put_string(block, CTF_BASELINE, string);
	}
	if ( ctf->envelope(0) ) {
		string = ctf->envelope_equation();
		item_put_string(block, CTF_ENVELOPE, string);
	}

	return 0;
}

int		micrograph_to_starblock(Bmicrograph* mg, Bstar_block* block,
				FOMType fom_tag[NFOM], int euler_flag, int omega_flag)
{
	float			angle;

	item_put_string(block, MICROGRAPH_ID, mg->id);
	block->tag = mg->id;
	item_put_integer(block, MICROGRAPH_NUMBER, mg->img_num, "%d");
	item_put_integer(block, MICROGRAPH_SELECT, mg->select, "%d");
	item_put_float(block, MICROGRAPH_FOM, mg->fom, "%f");
	if ( mg->fmg.length() ) item_put_string(block, MICROGRAPH_FILE, mg->fmg);
	if ( mg->fframe.length() ) item_put_string(block, MICROGRAPH_FRAMES_FILE, mg->fframe);
	if ( mg->fft.length() ) item_put_string(block, MICROGRAPH_TRANSFORM_FILE, mg->fft);
	if ( mg->fps.length() ) item_put_string(block, MICROGRAPH_POWERSPEC_FILE, mg->fps);
	if ( mg->fpart.length() && mg->part && !mg->part->fpart.length() )
		item_put_string(block, PARTICLE_FILE, mg->fpart);
	if ( mg->ffil.length() ) item_put_string(block, FILAMENT_FILE, mg->ffil);
	item_put_float(block, MICROGRAPH_MAGNIFICATION, mg->magnification, "%f");
	item_put_float(block, MICROGRAPH_SAMPLING, mg->sampling, "%f");
	item_put_float(block, MICROGRAPH_PIXEL_X, mg->pixel_size[0], "%f");
	item_put_float(block, MICROGRAPH_PIXEL_Y, mg->pixel_size[1], "%f");
	item_put_float(block, MICROGRAPH_DOSE, mg->dose, "%f");
	item_put_float(block, MICROGRAPH_INTENSITY, mg->intensity, "%f");
	item_put_float(block, MICROGRAPH_WATER_RING, mg->wri, "%f");
	angle = mg->tilt_axis*180.0/M_PI;
	item_put_float(block, MICROGRAPH_TILT_AXIS, angle, "%f");
	angle = mg->tilt_angle*180.0/M_PI;
	item_put_float(block, MICROGRAPH_TILT_ANGLE, angle, "%f");
	angle = mg->level_angle*180.0/M_PI;
	item_put_float(block, MICROGRAPH_LEVEL_ANGLE, angle, "%f");
	angle = mg->rot_angle*180.0/M_PI;
	item_put_float(block, MICROGRAPH_ROT_ANGLE, angle, "%f");
	item_put_float(block, MICROGRAPH_ORIGIN_X, mg->origin[0], "%f");
	item_put_float(block, MICROGRAPH_ORIGIN_Y, mg->origin[1], "%f");
	item_put_float(block, MICROGRAPH_ORIGIN_Z, mg->origin[2], "%f");
	item_put_float(block, MICROGRAPH_SCALE_X, mg->scale[0], "%f");
	item_put_float(block, MICROGRAPH_SCALE_Y, mg->scale[1], "%f");
	item_put_float(block, MICROGRAPH_SCALE_Z, mg->scale[2], "%f");
	item_put_float(block, MICROGRAPH_MATRIX_1_1, mg->matrix[0][0], "%f");
	item_put_float(block, MICROGRAPH_MATRIX_1_2, mg->matrix[0][1], "%f");
	item_put_float(block, MICROGRAPH_MATRIX_1_3, mg->matrix[0][2], "%f");
	item_put_float(block, MICROGRAPH_MATRIX_2_1, mg->matrix[1][0], "%f");
	item_put_float(block, MICROGRAPH_MATRIX_2_2, mg->matrix[1][1], "%f");
	item_put_float(block, MICROGRAPH_MATRIX_2_3, mg->matrix[1][2], "%f");
	item_put_float(block, MICROGRAPH_MATRIX_3_1, mg->matrix[2][0], "%f");
	item_put_float(block, MICROGRAPH_MATRIX_3_2, mg->matrix[2][1], "%f");
	item_put_float(block, MICROGRAPH_MATRIX_3_3, mg->matrix[2][2], "%f");
	item_put_float(block, MICROGRAPH_HVEC_X, mg->hvec[0], "%f");
	item_put_float(block, MICROGRAPH_HVEC_Y, mg->hvec[1], "%f");
	item_put_float(block, MICROGRAPH_HVEC_Z, mg->hvec[2], "%f");
	item_put_float(block, MICROGRAPH_KVEC_X, mg->kvec[0], "%f");
	item_put_float(block, MICROGRAPH_KVEC_Y, mg->kvec[1], "%f");
	item_put_float(block, MICROGRAPH_KVEC_Z, mg->kvec[2], "%f");
	item_put_float(block, MICROGRAPH_LVEC_X, mg->lvec[0], "%f");
	item_put_float(block, MICROGRAPH_LVEC_Y, mg->lvec[1], "%f");
	item_put_float(block, MICROGRAPH_LVEC_Z, mg->lvec[2], "%f");
	angle = mg->helix_axis*180.0/M_PI;
	item_put_float(block, MICROGRAPH_HELIX_AXIS, angle, "%f");
	item_put_float(block, MICROGRAPH_HELIX_RISE, mg->helix_rise, "%f");
	angle = mg->helix_angle*180.0/M_PI;
	item_put_float(block, MICROGRAPH_HELIX_ANGLE, angle, "%f");
	item_put_float(block, MICROGRAPH_HELIX_RADIUS, mg->helix_radius, "%f");
	item_put_float(block, PARTICLE_BOX_SIZE_X, mg->box_size[0], "%f");
	item_put_float(block, PARTICLE_BOX_SIZE_Y, mg->box_size[1], "%f");
	item_put_float(block, PARTICLE_BOX_SIZE_Z, mg->box_size[2], "%f");
	item_put_float(block, PARTICLE_BAD_RADIUS, mg->bad_radius, "%f");
	item_put_float(block, FILAMENT_WIDTH, mg->filament_width, "%f");
	item_put_float(block, FILNODE_RADIUS, mg->fil_node_radius, "%f");
	item_put_float(block, REFLEX_RADIUS, mg->sf_radius, "%f");
	item_put_float(block, MARKER_RADIUS, mg->mark_radius, "%f");

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

int		reconstruction_to_starblock(Breconstruction* rec, Bstar_block* block,
				FOMType fom_tag[NFOM], int euler_flag, int omega_flag)
{
	float			angle;

	item_put_string(block, MAP_ID, rec->id);
	block->tag = rec->id;

	if ( rec->frec.length() ) item_put_string(block, MAP_RECONSTRUCTION, rec->frec);

	if ( rec->fft.length() ) item_put_string(block, MAP_TRANSFORM_FILE, rec->fft);

	if ( rec->fps.length() ) item_put_string(block, MAP_POWERSPEC_FILE, rec->fps);

	if ( rec->fpart.length() && !rec->part->fpart.length() )
		item_put_string(block, PARTICLE_FILE, rec->fpart);

	if ( rec->ffil.length() ) item_put_string(block, FILAMENT_FILE, rec->ffil);

	if ( rec->model )
		item_put_list(block, MAP_MODEL, (char *)rec->model, 
				(char *)rec->model->c_str() - (char *)rec->model, "%s");

	item_put_integer(block, MAP_SELECT, rec->select, "%d");
	item_put_float(block, MAP_FOM, rec->fom, "%f");

	item_put_float(block, MAP_ORIGIN_X, rec->origin[0], "%f");
	item_put_float(block, MAP_ORIGIN_Y, rec->origin[1], "%f");
	item_put_float(block, MAP_ORIGIN_Z, rec->origin[2], "%f");
	item_put_float(block, MAP_SCALE_X, rec->scale[0], "%f");
	item_put_float(block, MAP_SCALE_Y, rec->scale[1], "%f");
	item_put_float(block, MAP_SCALE_Z, rec->scale[2], "%f");
	item_put_float(block, MAP_VOXEL_SIZE_X, rec->voxel_size[0], "%f");
	item_put_float(block, MAP_VOXEL_SIZE_Y, rec->voxel_size[1], "%f");
	item_put_float(block, MAP_VOXEL_SIZE_Z, rec->voxel_size[2], "%f");

	if ( rec->symmetry.length() ) item_put_string(block, MAP_SYMMETRY, rec->symmetry);

	if ( rec->ctf ) ctf_to_starblock(rec->ctf, block);

	item_put_float(block, PARTICLE_BOX_SIZE_X, rec->box_size[0], "%f");
	item_put_float(block, PARTICLE_BOX_SIZE_Y, rec->box_size[1], "%f");
	item_put_float(block, PARTICLE_BOX_SIZE_Z, rec->box_size[2], "%f");
	item_put_float(block, PARTICLE_BAD_RADIUS, rec->bad_radius, "%f");
	item_put_float(block, FILAMENT_WIDTH, rec->filament_width, "%f");
	item_put_float(block, FILNODE_RADIUS, rec->fil_node_radius, "%f");
	item_put_float(block, REFLEX_RADIUS, rec->sf_radius, "%f");
	item_put_float(block, MARKER_RADIUS, rec->mark_radius, "%f");

	item_put_float(block, MAP_VIEW_X, rec->view[0], "%f");
	item_put_float(block, MAP_VIEW_Y, rec->view[1], "%f");
	item_put_float(block, MAP_VIEW_Z, rec->view[2], "%f");
	angle = rec->view.angle()*180.0/M_PI;
	item_put_float(block, MAP_VIEW_ANGLE, angle, "%f");

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
@param	*star			STAR data base.
@param 	mg_select		flag to only convert selected micrographs.
@param 	rec_select		flag to only convert selected reconstructions.
@return	int				error code (<0 means failure).

	The function packes new micrograph data back into an existing STAR
	data base from which the old project parameters were obtained.
	All angles used internally are in radians, and are converted here
	to degrees for output. The angles within the micrograph structure
	are left in radians.
**/
int 		project_to_star(Bproject* project, Bstar* star, int mg_select, int rec_select)
{
	int				err(0);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_to_star:" << endl;

	star->comment = project->comment;
	
	star->split = project->split;
	
	Bstar_block*		block = NULL;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bstring*			file_list = NULL;

	if ( project->reference ) {
		block = (Bstar_block *) add_item((char **)&star->block, sizeof(Bstar_block));
		item_put_string_list(block, MAP_REFERENCE, project->reference);
		if ( verbose & VERB_DEBUG )
			for ( file_list=project->reference; file_list; file_list=file_list->next )
				cout << "DEBUG project_to_star: reference map = " << *file_list << endl;
	}
	
	long	 			i, f, nfield(0), nmg(0), nrec(0), npart(0);

	for ( f=i=0; f<NFOM; f++ ) i += project->fom_tag[f];
	if ( i < 2 ) project->fom_tag[0] = FOM;
	
	if ( verbose & VERB_PROCESS )
		cout << "Converting a project into a STAR database" << endl;
	
	if ( project->euler_flag < 1 ) {
		item_delete_all(star, PARTICLE_OMEGA);
		item_delete_all(star, PARTICLE_PSI);
		item_delete_all(star, PARTICLE_THETA);
		item_delete_all(star, PARTICLE_PHI);
	} else {
		item_delete_all(star, PARTICLE_VIEW_X);
		item_delete_all(star, PARTICLE_VIEW_Y);
		item_delete_all(star, PARTICLE_VIEW_Z);
		item_delete_all(star, PARTICLE_VIEW_ANGLE);
	}
	
	for ( field=project->field; field; field=field->next ) if ( mg_select < 1 || field->select > 0 ) {
		for ( mg=field->mg; mg; mg=mg->next ) if ( mg_select < 1 || mg->select > 0 ) {
			if ( verbose & VERB_FULL )
				cout << "Writing field \"" << field->id << "\", micrograph \"" << mg->id << "\"" << endl;
			block = (Bstar_block *) add_item((char **)&star->block, sizeof(Bstar_block));
			item_put_string(block, MICROGRAPH_FIELD_ID, field->id);
			micrograph_to_starblock(mg, block, project->fom_tag, project->euler_flag, project->omega_flag);
			nmg++;
		}
		nfield++;
	}
	
	if ( project->class_avg ) {
		npart = particle_count(project->class_avg);
		if ( verbose & VERB_FULL )
			cout << "Writing " << npart << " class averages" << endl;
		block = (Bstar_block *) add_item((char **)&star->block, sizeof(Bstar_block));
		block->tag = CLASS_AVERAGE;
		particle_to_starblock(project->class_avg, block, project->fom_tag, 
			project->euler_flag, project->omega_flag);
	}
	
	for ( rec = project->rec; rec; rec = rec->next ) if ( rec_select < 1 || rec->select > 0 ) {
		if ( verbose & VERB_FULL )
			cout << "Writing reconstruction \"" << rec->id << "\"" << endl;
		block = (Bstar_block *) add_item((char **)&star->block, sizeof(Bstar_block));
		reconstruction_to_starblock(rec, block, project->fom_tag, 
			project->euler_flag, project->omega_flag);
		nrec++;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << nfield << " fields, " << nmg << " micrographs, " << nrec << " reconstructions written" << endl;
	
	return err;
}

