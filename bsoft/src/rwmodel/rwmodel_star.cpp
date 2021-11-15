/**
@file	rwmodel_star.cpp
@brief	Library routines to read and write STAR model parameters
@author Bernard Heymann
@date	Created: 20060919
@date	Modified: 20210819
**/

#include "rwmodel.h"
#include "star.h"
#include "model_links.h"
#include "model_util.h"
#include "model_tags.h"
#include "file_util.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

map<string, int>	comptype_tags()
{
	map<string, int>	tag;
	tag[COMPTYPE_ID] = 0;
	tag[COMPTYPE_FILENAME] = 1;
	tag[COMPTYPE_NUMBER] = 2;
	tag[COMPTYPE_MASS] = 3;
	tag[COMPTYPE_FOM] = 4;
	tag[COMPTYPE_SELECT] = 5;

	return tag;
}

map<string, int>	linktype_tags()
{
	map<string, int>	tag;
	tag[LINKTYPE_ID1] = 0;
	tag[LINKTYPE_ID2] = 1;
	tag[LINKTYPE_LENGTH] = 2;
	tag[LINKTYPE_DISTANCE] = 3;
	tag[LINKTYPE_KLENGTH] = 4;
	tag[LINKTYPE_KDISTANCE] = 5;
	tag[LINKTYPE_FOM] = 6;
	tag[LINKTYPE_SELECT] = 7;

	return tag;
}

map<string, int>	angletype_tags()
{
	map<string, int>	tag;
	tag[ANGLETYPE_ID1] = 0;
	tag[ANGLETYPE_ID2] = 1;
	tag[ANGLETYPE_ID3] = 2;
	tag[ANGLETYPE_ANGLE] = 3;
	tag[ANGLETYPE_KANGLE] = 4;
	tag[ANGLETYPE_FOM] = 5;
	tag[ANGLETYPE_SELECT] = 6;

	return tag;
}

map<string, int>	component_tags()
{
	map<string, int>	tag;
	tag[COMPONENT_ID] = 0;
	tag[COMPONENT_TYPE_ID] = 1;
	tag[COMPONENT_X] = 2;
	tag[COMPONENT_Y] = 3;
	tag[COMPONENT_Z] = 4;
	tag[COMPONENT_VIEW_X] = 5;
	tag[COMPONENT_VIEW_Y] = 6;
	tag[COMPONENT_VIEW_Z] = 7;
	tag[COMPONENT_VIEW_ANGLE] = 8;
	tag[COMPONENT_RADIUS] = 9;
	tag[COMPONENT_RED] = 10;
	tag[COMPONENT_GREEN] = 11;
	tag[COMPONENT_BLUE] = 12;
	tag[COMPONENT_ALPHA] = 13;
	tag[COMPONENT_DENSITY] = 14;
	tag[COMPONENT_FOM] = 15;
	tag[COMPONENT_SELECT] = 16;

	return tag;
}

map<string, int>	link_tags()
{
	map<string, int>	tag;
	tag[COMPLINK_1] = 0;
	tag[COMPLINK_2] = 1;
	tag[COMPLINK_ANGLE] = 2;
	tag[COMPLINK_RADIUS] = 3;
	tag[COMPLINK_LENGTH] = 4;
	tag[COMPLINK_RED] = 5;
	tag[COMPLINK_GREEN] = 6;
	tag[COMPLINK_BLUE] = 7;
	tag[COMPLINK_ALPHA] = 8;
	tag[COMPLINK_FOM] = 9;
	tag[COMPLINK_SELECT] = 10;

	return tag;
}

/**
@brief 	Reads STAR model parameters.
@param 	*file_list	list of model parameter file names.
@return Bmodel*		model parameters.
**/
Bmodel*		read_model_star(Bstring* file_list)
{
	Bstar2			star;
	star.line_length(200);                // Set the output line length

	while ( file_list ) {
		star.read(file_list->str());
		file_list = file_list->next;
	}
	
	if ( verbose )
		cout << "Converting from STAR to model" << endl;
	
//	item_change_tag(star, COMPTYPE_MAP_NUMBER, COMPTYPE_NUMBER);

	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	long			i(0), j, k, l;
	Bstring			path;
	string			id;
	RGBA<float>		rgba;
	map<string,Bcomponent*>	htab;
	
//	for ( auto ib = star.blocks().begin(); ib != star.blocks().end(); ++ib ) {
	for ( auto ib: star.blocks() ) {
//		cout << ib.tag() << endl;
		path = ib.file_name();
		path = path.substr(0, path.rfind("/"));
//		cout << path << endl;
//		ib.show_tags();
		id = ib.at(MODEL_ID);
//		if ( id.length() < 1 ) id = Bstring(++i, "%d");
		if ( id.size() < 1 ) id = to_string(++i);
		if ( model ) mp = model->add(id);
		else mp = model = new Bmodel(id);
		mp->model_type(ib.at(MODEL_TYPE));
		mp->mapfile(find_file(ib.at(MODEL_MAP_FILENAME), path).str());
		mp->image_number(stol(ib.at(MODEL_MAP_NUMBER)));
		mp->symmetry(ib.at(MODEL_SYM));
		mp->handedness(stol(ib.at(MODEL_HAND)));
		mp->FOM(stod(ib.at(MODEL_FOM)));
		mp->select(stol(ib.at(MODEL_SELECT)));

//		for ( auto il = ib->loops().begin(); il != ib->loops().end(); ++il ) {
		for ( auto il: ib.loops() ) {
//			il.show_tags();
			if ( ( k = il.find(COMPTYPE_ID) ) >= 0 ) {
//				cout << "COMPTYPE_ID found: " << COMPTYPE_ID << endl;
				Bcomptype*			comptype = NULL;
				for ( auto ir: il.data() ) {
					id = ir[k];
					if ( comptype ) comptype = comptype->add(id);
					else mp->type = comptype = new Bcomptype(id);
					if ( ( j = il.find(COMPTYPE_FILENAME) ) >= 0 ) comptype->file_name(ir[j]);
					if ( ( j = il.find(COMPTYPE_NUMBER) ) >= 0 ) comptype->image_number(stol(ir[j]));
					if ( ( j = il.find(COMPTYPE_MASS) ) >= 0 ) comptype->mass(stod(ir[j]));
					if ( ( j = il.find(COMPTYPE_FOM) ) >= 0 ) comptype->FOM(stod(ir[j]));
					if ( ( j = il.find(COMPTYPE_SELECT) ) >= 0 ) comptype->select(stol(ir[j]));
				}
			} else if ( ( k = il.find(COMPONENT_ID) ) >= 0 ) {
//				cout << "COMPONENT_ID found: " << COMPONENT_ID << endl;
				l = il.find(COMPONENT_TYPE_ID);
				Bcomponent*		comp = NULL;
				Vector3<double>	loc;
				View2<float>	view;
				for ( auto ir: il.data() ) {
					id = ir[k];
					if ( comp ) comp = comp->add(id);
					else mp->comp = comp = new Bcomponent(id);
					htab[id] = comp;
					id = ir[l];
					comp->type(mp->type->find(id));
					if ( ( j = il.find(COMPONENT_X) ) >= 0 ) loc[0] = stod(ir[j]);
					if ( ( j = il.find(COMPONENT_Y) ) >= 0 ) loc[1] = stod(ir[j]);
					if ( ( j = il.find(COMPONENT_Z) ) >= 0 ) loc[2] = stod(ir[j]);
					comp->location(loc);
					if ( ( j = il.find(COMPONENT_VIEW_X) ) >= 0 ) view[0] = stod(ir[j]);
					if ( ( j = il.find(COMPONENT_VIEW_Y) ) >= 0 ) view[1] = stod(ir[j]);
					if ( ( j = il.find(COMPONENT_VIEW_Z) ) >= 0 ) view[2] = stod(ir[j]);
					if ( ( j = il.find(COMPONENT_VIEW_ANGLE) ) >= 0 ) view[3] = stod(ir[j]) *M_PI/180.0;
					comp->view(view);
					if ( ( j = il.find(COMPONENT_RADIUS) ) >= 0 ) comp->radius(stod(ir[j]));
					if ( ( j = il.find(COMPONENT_RED) ) >= 0 ) rgba[0] = stod(ir[j]);
					if ( ( j = il.find(COMPONENT_GREEN) ) >= 0 ) rgba[1] = stod(ir[j]);
					if ( ( j = il.find(COMPONENT_BLUE) ) >= 0 ) rgba[2] = stod(ir[j]);
					if ( ( j = il.find(COMPONENT_ALPHA) ) >= 0 ) rgba[3] = stod(ir[j]);
					comp->color(rgba);
					if ( ( j = il.find(COMPONENT_DENSITY) ) >= 0 ) comp->density(stod(ir[j]));
					if ( ( j = il.find(COMPONENT_FOM) ) >= 0 ) comp->FOM(stod(ir[j]));
					if ( ( j = il.find(COMPONENT_SELECT) ) >= 0 ) comp->select(stol(ir[j]));
				}
			} else if ( ( k = il.find(COMPLINK_1) ) >= 0 ) {
//				cout << "COMPLINK_1 found: " << COMPLINK_1 << endl;
				if ( htab.size() < 1 ) {
					cerr << "No components found to link up!" << endl;
					bexit(-1);
				}
				l = il.find(COMPLINK_2);
				Blink*			link = NULL;
				Bcomponent*		comp1 = NULL;
				Bcomponent*		comp2 = NULL;
				for ( auto ir: il.data() ) {
//					comp1 = mp->comp->find(ir[k]);
//					comp2 = mp->comp->find(ir[l]);
					comp1 = htab[ir[k]];
					comp2 = htab[ir[l]];
					if ( link ) link = link->add(comp1, comp2);
					else link = mp->link = new Blink(comp1, comp2);
					if ( ( j = il.find(COMPLINK_ANGLE) ) >= 0 ) link->angle(stod(ir[j]));
					if ( ( j = il.find(COMPLINK_RADIUS) ) >= 0 ) link->radius(stod(ir[j]));
					if ( ( j = il.find(COMPLINK_LENGTH) ) >= 0 ) link->length(stod(ir[j]));
					if ( ( j = il.find(COMPLINK_RED) ) >= 0 ) rgba[0] = stod(ir[j]);
					if ( ( j = il.find(COMPLINK_GREEN) ) >= 0 ) rgba[1] = stod(ir[j]);
					if ( ( j = il.find(COMPLINK_BLUE) ) >= 0 ) rgba[2] = stod(ir[j]);
					if ( ( j = il.find(COMPLINK_ALPHA) ) >= 0 ) rgba[3] = stod(ir[j]);
					link->color(rgba);
					if ( ( j = il.find(COMPLINK_FOM) ) >= 0 ) link->FOM(stod(ir[j]));
					if ( ( j = il.find(COMPLINK_SELECT) ) >= 0 ) link->select(stol(ir[j]));
				}
//				cout << "Setting up link lengths" << endl;
				for ( link=mp->link; link; link=link->next )
					link->length(link->comp[0]->location().distance(link->comp[1]->location()));
			}
		}
	}

	model->comment(star.comment());
	
	models_process(model, model_setup_links);

	return model;
}

/**
@brief 	Writes STAR model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	split		flag to split into data blocks.
@return int			models written.
**/
int			write_model_star(Bstring& filename, Bmodel* model, int split)
{
	Bstring			id;
 	Bstar2			star;
	
	Bmodel*			mp = NULL;
	Bcomptype*		comptype = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;
	
	star.comment(model->comment());
	star.line_length(200);
	
//	item_change_tag(star, COMPTYPE_MAP_NUMBER, COMPTYPE_NUMBER);
	
	int				n, err(0);
	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
//		cout << "writing model " << mp->identifier() << endl;
		BstarBlock&	block = star.add_block(mp->identifier());
		block[MODEL_ID] = mp->identifier();
		block[MODEL_TYPE] = mp->model_type();
		block[MODEL_SYM] = mp->symmetry();
		block[MODEL_HAND] = to_string(mp->handedness());
		block[MODEL_MAP_FILENAME] = mp->mapfile();
		block[MODEL_MAP_NUMBER] = to_string(mp->image_number());
		block[MODEL_FOM] = to_string(mp->FOM());
		block[MODEL_SELECT] = to_string(mp->select());
//		cout << "writing component types" << endl;
		if ( mp->type ) {
			BstarLoop&			loop = block.add_loop();
			loop.tags() = comptype_tags();
			for ( comptype = mp->type; comptype; comptype = comptype->next ) {
				vector<string>&	vs = loop.add_row(6);
				vs[0] = comptype->identifier();
				vs[1] = comptype->file_name();
				vs[2] = to_string(comptype->image_number());
				vs[3] = to_string(comptype->mass());
				vs[4] = to_string(comptype->FOM());
				vs[5] = to_string(comptype->select());
			}
		}
//		cout << "writing components" << endl;
		if ( mp->comp ) {
			BstarLoop&			loop = block.add_loop();
			loop.tags() = component_tags();
			for ( comp = mp->comp; comp; comp = comp->next ) {
				vector<string>&	vs = loop.add_row(17);
				vs[0] = comp->identifier();
				if ( comp->type() ) id = comp->type()->identifier();
				else id = "?";
				vs[1] = id.str();
				vs[2] = to_string(comp->location()[0]);
				vs[3] = to_string(comp->location()[1]);
				vs[4] = to_string(comp->location()[2]);
				vs[5] = to_string(comp->view()[0]);
				vs[6] = to_string(comp->view()[1]);
				vs[7] = to_string(comp->view()[2]);
				vs[8] = to_string(comp->view()[3]*180.0/M_PI);
				vs[9] = to_string(comp->radius());
				vs[10] = to_string(comp->color()[0]);
				vs[11] = to_string(comp->color()[1]);
				vs[12] = to_string(comp->color()[2]);
				vs[13] = to_string(comp->color()[3]);
				vs[14] = to_string(comp->density());
				vs[15] = to_string(comp->FOM());
				vs[16] = to_string(comp->select());
			}
		}
//		cout << "writing links" << endl;
		if ( mp->link ) {
			BstarLoop&			loop = block.add_loop();
			loop.tags() = link_tags();
			for ( link = mp->link; link; link = link->next ) {
				vector<string>&	vs = loop.add_row(11);
				vs[0] = link->comp[0]->identifier();
				vs[1] = link->comp[1]->identifier();
				vs[2] = to_string(link->angle()*180.0/M_PI);
				vs[3] = to_string(link->radius());
				vs[4] = to_string(link->length());
				vs[5] = to_string(link->color()[0]);
				vs[6] = to_string(link->color()[1]);
				vs[7] = to_string(link->color()[2]);
				vs[8] = to_string(link->color()[3]);
				vs[9] = to_string(link->FOM());
				vs[10] = to_string(link->select());
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_star: " << filename << endl;

	err = star.write(filename.str(), split);
	
	if ( err < 0 ) return err;
	
	return  n;
}

