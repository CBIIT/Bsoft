/**
@file	rwmodel_cmm.cpp
@brief	Library routines to read and write Chimera marker model parameters
@author Bernard Heymann
@date	Created: 20060919
@date	Modified: 20161004
**/

#include "rwmodel.h"
#include "file_util.h"
#include "linked_list.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

struct TagValue {
	TagValue*	next;
	Bstring		tag;
	Bstring		value;
} ;

TagValue*	tag_value_parse(char* string)
{
	TagValue*	tvlist = NULL;
	TagValue*	tv = NULL;
	
	Bstring		s(string);
	s = s.within('<', '>');
	
	Bstring*	slist = s.split();
	Bstring*	sp;
	
	for ( sp = slist; sp; sp = sp->next ) {
		tv = (TagValue *) add_item((char **) &tv, sizeof(TagValue));
		if ( !tvlist ) tvlist = tv;
//		cout << "sp = " << *sp << endl;
		if ( sp->contains("=") ) {
			tv->tag = sp->pre('=');
			tv->value = sp->post('\"');
			while ( tv->value.index('\"', 0) < 0 ) { sp = sp->next; tv->value += " " + *sp; }
			tv->value = tv->value.pre('\"');
		} else tv->tag = *sp;
	}
	
	string_kill(slist);
	
	return tvlist;
}

int			tag_value_kill(TagValue* tv)
{
	TagValue*	tv2;
	
	while ( tv ) {
		tv2 = tv->next;
		tv->tag = 0;
		tv->value = 0;
		delete[] (char *)tv;
		tv = tv2;
	}
	
	return 0;
}

/**
@brief 	Reads Chimera marker model parameters.
@param 	*file_list	list of model parameter file names.
@return Bmodel*		model parameters.
**/
Bmodel*		read_model_chimera(Bstring* file_list)
{
	int				i;
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomptype*		ct = NULL;
	Bcomponent*		comp = NULL;
	Bcomponent*		comp1 = NULL;
	Bcomponent*		comp2 = NULL;
	Blink*			link = NULL;
	Bstring			id, path;
	Bstring*		filename;
	ifstream		fmod;
	char			aline[1024];
	RGBA<float>		rgba(1,1,1,1);	// Default white
	TagValue*		tvlist = NULL;
	TagValue*		tv = NULL;

	for ( i=1, filename = file_list; filename; filename = filename->next, i++ ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << *filename << endl;
		fmod.open(filename->c_str());
		if ( fmod.fail() ) return NULL;
		path = filename->pre_rev('/');
		while ( fmod.getline(aline, 1024) ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG read_model_chimera: " << aline << endl;
			tvlist = tag_value_parse(aline);
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG read_model_chimera: tag=" << tvlist->tag << endl;
			if ( tvlist->tag == "marker_set" ) {
				link = NULL;
				comp = NULL;
				mp = model_add(&mp, "1");
				if ( !model ) model = mp;
				for ( tv = tvlist; tv; tv = tv->next ) if ( tv->value.length() ) {
					if ( verbose & VERB_DEBUG )
						cout << "DEBUG read_model_chimera: tag=" << tvlist->tag << " value=" << tv->value << endl;
					if ( tv->tag == "name" ) mp->identifier(tv->value.str());
					if ( tv->tag == "type" ) mp->model_type(tv->value.str());
					if ( tv->tag == "hand" ) mp->handedness(tv->value.integer());
					if ( tv->tag == "symmetry" ) mp->symmetry(tv->value.str());
					if ( tv->tag == "file" ) mp->mapfile(find_file(tv->value, path).str());
					if ( tv->tag == "img_num" ) mp->image_number(tv->value.integer());
					if ( tv->tag == "fom" ) mp->FOM(tv->value.real());
					if ( tv->tag == "select" ) mp->select(tv->value.integer());
				}
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG read_model_chimera: model id = " << mp->identifier() << endl;
			} else if ( tvlist->tag == "comment" ) {
				while ( fmod.getline(aline, 1024) && strncmp(aline, "</comment>", 10) ) {
					mp->comment(mp->comment() + aline);
					mp->comment(mp->comment() + "\n");
				}
			} else if ( tvlist->tag == "type" ) {
				for ( tv = tvlist; tv; tv = tv->next )
					if ( tv->tag == "id" ) id = tv->value;
				ct = mp->add_type(id);
				if ( !mp->type ) mp->type = ct;
				for ( tv = tvlist; tv; tv = tv->next ) {
					if ( tv->tag == "file" ) ct->file_name(tv->value.str());
					if ( tv->tag == "num" ) ct->image_number(tv->value.integer());
					if ( tv->tag == "mass" ) ct->mass(tv->value.real());
					if ( tv->tag == "fom" ) ct->FOM(tv->value.real());
					if ( tv->tag == "select" ) ct->select(tv->value.integer());
				}
			} else if ( tvlist->tag == "marker" ) {
				for ( tv = tvlist; tv; tv = tv->next )
					if ( tv->tag == "id" ) id = tv->value;
//				comp = component_add(&comp, id);
//				if ( !mp->comp ) mp->comp = comp;
				if ( comp ) comp = comp->add(id.str());
				else comp = mp->add_component(id.str());
				for ( tv = tvlist; tv; tv = tv->next ) {
					if ( tv->tag == "type" )
//						comp->type = model_add_type_by_id(mp, tv->value);
						comp->type(mp->add_type(tv->value.str()));
					if ( tv->tag == "x" ) comp->location()[0] = tv->value.real();
					if ( tv->tag == "y" ) comp->location()[1] = tv->value.real();
					if ( tv->tag == "z" ) comp->location()[2] = tv->value.real();
					if ( tv->tag == "vx" ) comp->view()[0] = tv->value.real();
					if ( tv->tag == "vy" ) comp->view()[1] = tv->value.real();
					if ( tv->tag == "vz" ) comp->view()[2] = tv->value.real();
					if ( tv->tag == "va" ) comp->view()[3] = tv->value.real()*M_PI/180.0;
					if ( tv->tag == "radius" ) comp->radius(tv->value.real());
					if ( tv->tag == "r" ) comp->color()[0] = tv->value.real();
					if ( tv->tag == "g" ) comp->color()[1] = tv->value.real();
					if ( tv->tag == "b" ) comp->color()[2] = tv->value.real();
					if ( tv->tag == "a" ) comp->color()[3] = tv->value.real();
					if ( tv->tag == "density" ) comp->density(tv->value.real());
					if ( tv->tag == "fom" ) comp->FOM(tv->value.real());
					if ( tv->tag == "select" ) comp->select(tv->value.integer());
				}
			} else if ( tvlist->tag == "link" ) {
				comp1 = comp2 = NULL;
				for ( tv = tvlist; tv; tv = tv->next ) {
					if ( tv->tag == "id1" )
						for ( comp = mp->comp; comp; comp = comp->next )
							if ( comp->identifier() == tv->value.str() ) comp1 = comp;
					if ( tv->tag == "id2" )
						for ( comp = mp->comp; comp; comp = comp->next )
							if ( comp->identifier() == tv->value.str() ) comp2 = comp;
				}
				link = link_add(&link, comp1, comp2, 0, 1);
				if ( !mp->link ) mp->link = link;
				for ( tv = tvlist; tv; tv = tv->next ) {
					if ( tv->tag == "radius" ) link->radius(tv->value.real());
					if ( tv->tag == "r" ) rgba[0] = tv->value.real();
					if ( tv->tag == "g" ) rgba[1] = tv->value.real();
					if ( tv->tag == "b" ) rgba[2] = tv->value.real();
					if ( tv->tag == "a" ) rgba[3] = tv->value.real();
					if ( tv->tag == "fom" ) link->FOM(tv->value.real());
					if ( tv->tag == "select" ) link->select(tv->value.integer());
				}
				link->color(rgba);
			}
			tag_value_kill(tvlist);
		}
		fmod.close();
	}

	return model;
}

/**
@brief 	Writes Chimera marker model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@return int			models written.
**/
int			write_model_chimera(Bstring& filename, Bmodel* model)	
{
	int				n;
	Bmodel*			mp = NULL;
	Bcomptype*		ct = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;
	Bstring			onename;

	ofstream		fmod;

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		if ( model->next )
			onename = filename.pre_rev('.') + Bstring(n+1, "_%04d.") + filename.post_rev('.');
		else
			onename = filename;
		fmod.open(onename.c_str());
		if ( fmod.fail() ) return  -1;
		fmod << "<?xml version=\"1.0\"?>" << endl;
		fmod << "<marker_set name=\"" << mp->identifier() << "\" type=\"" << mp->model_type()
			<< "\" hand=\"" << mp->handedness() << "\" symmetry=\"" << mp->symmetry()
			<< "\" file=\"" << mp->mapfile() << "\" img_num=\"" << mp->image_number()
			<< "\" fom=\"" << mp->FOM()<< "\" select=\"" << mp->select()<< "\">" << endl;
		fmod << "<comment>" << endl << mp->comment() << endl << "</comment>" << endl;
		for ( ct = mp->type; ct; ct = ct->next ) {
			fmod << "<type id=\"" << ct->identifier() << "\" file=\"" << ct->file_name()
				<< "\" num=\"" << ct->image_number() << "\" mass=\"" << ct->mass()
				<< "\" fom=\"" << ct->FOM()<< "\" select=\"" << ct->select()<< "\"/>" << endl;
		}
		for ( comp = mp->comp; comp; comp = comp->next ) {
			fmod << "<marker id=\"" << comp->identifier() << "\" type=\"" << comp->type()->identifier()
				<< "\" x=\"" << comp->location()[0] << "\" y=\"" << comp->location()[1]
				<< "\" z=\"" << comp->location()[2] << "\" vx=\"" << comp->view()[0]
				<< "\" vy=\"" << comp->view()[1] << "\" vz=\"" << comp->view()[2]
				<< "\" va=\"" << comp->view().angle()*180.0/M_PI
				<< "\" r=\"" << comp->color()[0] << "\" g=\"" << comp->color()[1]
				<< "\" b=\"" << comp->color()[2] << "\" a=\"" << comp->color()[3]
				<< "\" radius=\"" << comp->radius() << "\" density=\"" << comp->density()
				<< "\" fom=\"" << comp->FOM()<< "\" select=\"" << comp->select()<< "\"/>" << endl;
		}
		for ( link = mp->link; link; link = link->next ) {
			fmod << "<link id1=\"" << link->comp[0]->identifier() << "\" id2=\"" << link->comp[1]->identifier()
				<< "\" r=\"" << link->color()[0] << "\" g=\"" << link->color()[1]
				<< "\" b=\"" << link->color()[2] << "\" a=\"" << link->color()[3]
				<< "\" radius=\"" << link->radius() << "\" fom=\"" << link->FOM()
				<< "\" select=\"" << link->select()<< "\"/>" << endl;
		}
		fmod << "</marker_set>" << endl;
		fmod.close();
	}
	
	return  n;
}


