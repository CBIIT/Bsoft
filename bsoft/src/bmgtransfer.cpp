/**
@file	bmgtransfer.cpp
@brief	Transfer pieces of information between micrograph paramater files.
@author Bernard Heymann
@date	Created: 20140701
@date	Modified: 20140703
**/

#include "mg_processing.h"
#include "rwmg.h"
#include "mg_ctf.h"
#include "file_util.h"
#include "utilities.h"
#include "options.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

long		project_replace_parameters(Bproject* project, Bproject* proj_source, Bstring& which);

// Usage assistance
const char* use[] = {
" ",
"Usage: bmgtransfer [options] input.star/input.img [input.star/input.img]",
"------------------------------------------------------------------------",
"Transfer pieces of information from input micrograph parameter files",
"to a master destination parameter file.",
"Transfers are linked to field ID's and can occur in one of two ways:",
"1. From one source micrograph to all destination micrographs in the field.",
"2. From each source micrograph to the corresponding destination micrograph in the field.",
" ",
"Actions:",
"-ctf                     Replace CTF parameters.",
"-particles               Replace particle parameters.",
"-origins                 Replace particle origins.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
" ",
"Input:",
"-master destination.star Master parameter file into which to merge information from other input files.",
" ",
"Output:",
"-output output.star      Output parameter file (required any time \"-split\" is used to indicate the type).",
"-split 3                 Split micrograph or reconstruction data blocks into individual files:",
"                         Argument: 1-6: number of digits inserted before extension",
"                         Argument: \"id\": micrograph ID's are used as file names (no -output required).",
"                         Argument: \"field\": field ID's are used as file names (no -output required).",
"-dump file.txt           Dump all the particle origins and views.",
" ",
NULL
};

int main (int argc, char **argv)
{
	// Initialize optional variables
	Bstring			which;						// Specification of parameters to transfer
	Bstring			masterfile;					// Input master parameter file name
	Bstring			outfile;					// Output parameter file name
	Bstring			dumpfile;					// File to dump particle info to
	int				split(0);					// Output one big STAR file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "ctf" ) {
			if ( which.length() ) which += ",";
			which += "ctf";
		}
		if ( curropt->tag == "particles" ) {
			if ( which.length() ) which += ",";
			which += "particles";
		}
		if ( curropt->tag == "origins" ) {
			if ( which.length() ) which += ",";
			which += "origins";
		}
		if ( curropt->tag == "master" )
			masterfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "split" ) {
			if ( curropt->value.contains("id") || curropt->value.contains("ID") ) split = 9;
			else if ( curropt->value.contains("field") || curropt->value.contains("FIELD") ) split = -9;
			else if ( ( split = curropt->value.integer() ) < 1 )
				cerr << "-split: An integer must be specified!" << endl;
			else
				if ( split > 6 ) split = 6;
		}
		if ( curropt->tag == "dump" )
			dumpfile = curropt->filename();
    }
	option_kill(option);
	
	Bproject*		project = NULL;
	
	if ( masterfile.length() )
		project = read_project(masterfile);
	else {
		cerr << "Error: A master parameter file must be specified!" << endl;
		bexit(-1);
	}
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project2 = read_project(file_list);
		
	string_kill(file_list);

	if ( ! project2 ) {
		cerr << "Error: No input parameter files read!" << endl;
		bexit(-1);
	}

	if ( which.length() )
		project_replace_parameters(project, project2, which);

	project_kill(project2);
	
	if ( dumpfile.length() ) project_dump(project, dumpfile);

	
	if ( project ) {
		if ( split < 0 )
			project_split_field_write(project);
		else if ( outfile.length() || split == 9 ) {
			project->split = split;
			write_project(outfile, project, 0, 0);
		}
	}
	
	project_kill(project);

	bexit(0);
}


long		mg_replace_parameters(Bmicrograph* mg, Bmicrograph* mg_source, Bstring& which, Vector3<double> origin)
{
	Bparticle		*part, *part_source;
	Vector3<double>	shift = mg->origin - origin;
	
	if ( which.contains("ctf") ) {
		if ( mg_source->fps.length() > 0 ) mg->fps = mg_source->fps;
		if ( mg_source->magnification ) mg->magnification = mg_source->magnification;
		if ( mg_source->sampling ) mg->sampling = mg_source->sampling;
		if ( mg_source->ctf ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->update(mg_source->ctf);
		}
	}

	if ( which.contains("part") ) {
		mg->box_size = mg_source->box_size;
		mg->bad_radius = mg_source->bad_radius;
		particle_kill(mg->part);
		mg->part = particle_copy(mg_source->part);
		for ( part = mg->part; part; part = part->next )
			part->loc += shift;
	}

	if ( which.contains("ori") ) {
		for ( part = mg->part, part_source = mg_source->part; part && part_source; part = part->next, part_source = part_source->next )
			part->ori = part_source->ori;
	}
	
	return 1;
}

long		field_replace_parameters(Bfield* field, Bfield* field_source, Bstring& which)
{
	Bmicrograph		*mg, *mg_source;
	Vector3<double>	origin = field->mg->origin;	// Assumes the first micrograph is the reference

	long			nmg = field_count_micrographs(field);
	long			nmg_source = field_count_micrographs(field_source);
	
	if ( nmg_source == 1 ) {	// Replaces parameters from one micrograph into many
		mg_source = field_source->mg;
		for ( nmg=0, mg = field->mg; mg; mg = mg->next )
			nmg += mg_replace_parameters(mg, mg_source, which, origin);
	} else if ( nmg_source == nmg ) {	// Replaces parameters from into corresponding micrographs
		for ( nmg=0, mg = field->mg, mg_source = field_source->mg; mg && mg_source; mg = mg->next, mg_source = mg_source->next )
			nmg += mg_replace_parameters(mg, mg_source, which, origin);
	} else {
		cerr << "Error in updating parameters: The numbers of micrographs are unsuitable!" << endl;
		cerr << "	Master field " << field->id << " has " << nmg << " micrographs" << endl;
		cerr << "	Source field " << field_source->id << " has " << nmg_source << " micrographs" << endl;
	}
	
	return nmg;
}

long		rec_replace_parameters(Breconstruction* rec, Breconstruction* rec_source, Bstring& which)
{
	Bparticle		*part, *part_source;

	if ( which.contains("ctf") ) {
		if ( rec_source->fps.length() > 0 ) rec->fps = rec_source->fps;
		if ( rec_source->ctf ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			rec->ctf->update(rec_source->ctf);
		}
	}

	if ( which.contains("part") ) {
		rec->box_size = rec_source->box_size;
		rec->bad_radius = rec_source->bad_radius;
		particle_kill(rec->part);
		rec->part = particle_copy(rec_source->part);
	}
	
	if ( which.contains("ori") ) {
		for ( part = rec->part, part_source = rec_source->part; part && part_source; part = part->next, part_source = part_source->next )
			part->ori = part_source->ori;
	}
	
	return 1;
}

long		project_replace_parameters(Bproject* project, Bproject* proj_source, Bstring& which)
{
	long			nmg(0), nrec(0);
	Bfield 			*field, *field_source;
	Breconstruction	*rec, *rec_source;
	
	if ( verbose )
		cout << "Replacing parameters:           " << which << endl;

	if ( project->select < 1 ) {
		for ( field_source = proj_source->field; field_source; field_source = field_source->next ) {
			for ( field = project->field; field && field->id != field_source->id; field = field->next ) ;
			if ( field )
				nmg += field_replace_parameters(field, field_source, which);
		}
		if ( verbose )
			cout << "Micrographs updated:            " << nmg << endl << endl;
	} else {
		for ( rec_source = proj_source->rec; rec_source; rec_source = rec_source->next ) {
			for ( rec = project->rec; rec && rec->id != rec_source->id; rec = rec->next ) ;
			if ( rec )
				nrec += rec_replace_parameters(rec, rec_source, which);
		}
		if ( verbose )
			cout << "Reconstructions updated:        " << nmg << endl << endl;
	}
	

	return nmg;
}

