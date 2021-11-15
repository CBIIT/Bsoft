/**
@file	rwmodel_vega.cpp
@brief	Library routines to read and write Vega model parameters
@author Bernard Heymann
@date	Created: 20060919
@date	Modified: 20120211
**/

#include "rwmodel.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads Vega model parameters.
@param 	*file_list	list of model parameter file names.
@return Bmodel*		model parameters.
**/
Bmodel*		read_model_vega(Bstring* file_list)
{
	int				i, n, m, l1, l2, l3, readflag;
	float			x, y, z;
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;
	Bstring			id;
	Bstring			comptype("VER");
	Bstring*		filename;
	ifstream		fmod;
	char			aline[200];
	Bcomponent**	car;

	for ( i=1, filename = file_list; filename; filename = filename->next, i++ ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << *filename << endl;
		fmod.open(filename->c_str());
		if ( fmod.fail() ) return  NULL;
		readflag = m = 0;
//		while ( fmod.getline(aline, 200) != NULL ) {
		while ( fmod.getline(aline, 200) && fmod.good() ) {
			if ( readflag ) {
				sscanf(aline, "%d %f %f %f %d %d %d", &n, &x, &y, &z, &l1, &l2, &l3);
				if ( n > 0 ) {
					id = Bstring(n, "%d");
//					comp = component_add(&comp, id);
//					if ( !mp->comp ) mp->comp = comp;
//					comp = mp->add_component(id);
					if ( comp ) comp = comp->add(n);
					else model->comp = comp = new Bcomponent(n);
					comp->location(Vector3<float>(x, y, z));
					comp->flag[0] = l1;	// The flags carry the indices of the linked components
					comp->flag[1] = l2;
					comp->flag[2] = l3;
					comp->select(n);
//					comp->type = model_add_type_by_id(mp, comptype);
					comp->type(mp->add_type(comptype));
					if ( m < n ) m = n;
				} else readflag = 0;
			}
			if ( strstr(aline, "writegraph3d planar") ) {
				link = NULL;
				comp = NULL;
//				id = Bstring(i, "%d");
//				mp = model_add(&mp, id);
//				if ( !model ) model = mp;
				if ( mp ) mp = mp->add(i);
				else model = mp = new Bmodel(i);
				mp->model_type(mp->identifier());
				readflag = 1;
			}
		}
		fmod.close();
		m++;
		car = new Bcomponent*[m];
		for ( comp = mp->comp; comp; comp = comp->next ) car[comp->select()] = comp;
		for ( comp = mp->comp; comp; comp = comp->next ) {
			for ( n=0; n<comp->link.size() && comp->flag[n]>0; n++ ) {
				comp->link[n] = car[comp->flag[n]];
				if ( comp->select()< comp->flag[n] ) {
					link = link_add(&link, comp, comp->link[n], 0, 1);
					if ( !mp->link ) mp->link = link;
				}
				comp->flag[n] = 1;	// flag=1 indicates a link
			}
			comp->select(1);
		}
		delete[] car;
	}

	return model;
}

/**
@brief 	Writes Vega model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@return int			models written.
**/
int			write_model_vega(Bstring& filename, Bmodel* model)
{
	int				n;
	Bmodel*			mp = NULL;
	Bcomponent*		comp;
	Bstring			onename;
	
	ofstream		fmod;

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		if ( model->next )
			onename = filename.pre_rev('.') + Bstring(n+1, "_%04d.") + filename.post_rev('.');
		else
			onename = filename;
		fmod.open(onename.c_str());
		if ( fmod.fail() ) return  -1;
		fmod << ">>writegraph3d planar <<" << endl;
		for ( comp = mp->comp; comp; comp = comp->next )
			fmod << stol(comp->identifier()) << " " <<
				comp->location() << " " <<
				stol(comp->link[0]->identifier()) << " " <<
				stol(comp->link[1]->identifier()) << " " <<
				stol(comp->link[2]->identifier()) << endl;
		fmod << "0" << endl;
		fmod.close();
	}
	
	return 0;
}
