/**
@file	rwmodel_mol.cpp
@brief	Library routines to read and write atomic model parameters
@author Bernard Heymann
@date	Created: 20060919
@date	Modified: 20100302
**/

#include "rwmodel.h"
#include "rwmolecule.h"
#include "model_links.h"
#include "model_util.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads molecular model parameters.
@param 	*file_list	list of model parameter file names.
@param 	&paramfile	parameter file.
@return Bmodel*		model parameters.
**/
Bmodel*		read_model_molecule(Bstring* file_list, Bstring& paramfile)
{
	Bstring*		filename;
    Bstring    		atom_select("all");
	Bstring			id;
	Bmolgroup*		molgroup;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	Bbond*			bond;
	
	int				i;
	RGBA<float>		rgba(1,1,1,1);
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_model_molecule: " << *file_list << endl;
	
	for ( i=1, filename = file_list; filename; filename = filename->next, i++ ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << *filename << endl;
		molgroup = read_molecule(*filename, atom_select, paramfile);
//		id = Bstring(i, "%d");
//		mp = model_add(&mp, id);
//		if ( !model ) model = mp;
		if ( mp ) mp = mp->add(i);
		else mp = model = new Bmodel(i);
		mp->model_type(mp->identifier());
		if ( molgroup ) {
			if ( molgroup->id.length() && molgroup->id[0] != ' ') mp->identifier(molgroup->id);
			mp->symmetry(molgroup->pointgroup.str());
			comp = NULL;
			link = NULL;
			for ( bond = molgroup->bond; bond; bond = bond->next ) {
				link = (Blink *) add_item((char **) &link, sizeof(Blink));
				if ( !mp->link ) mp->link = link;
				link->radius(1);
				link->select(1);
				link->color(rgba);
			}
			for ( mol = molgroup->mol; mol; mol = mol->next ) {
				for ( res = mol->res; res; res = res->next ) {
					for ( atom = res->atom; atom; atom = atom->next ) {
//						id = Bstring(atom->num, "%d");
//						comp = component_add(&comp, id);
//						if ( !mp->comp ) mp->comp = comp;
//						comp = mp->add_component(id);
						if ( comp ) comp = comp->add(atom->num);
						else comp = mp->comp = new Bcomponent(atom->num);
						comp->location(atom->coord);
						comp->FOM(atom->b);
						comp->select((int) (atom->q + 0.999));
						id = atom->type;
						id = id.no_space();
//						comp->type = model_add_type_by_id(mp, id);
						comp->type(mp->add_type(id));
						for ( bond = molgroup->bond, link = mp->link; bond && link; bond = bond->next, link = link->next ) {
							if ( atom == bond->atom1 ) link->comp[0] = comp;
							else if ( atom == bond->atom2 ) link->comp[1] = comp;
						}
					}
				}
			}
			molgroup_kill(molgroup);
		}
	}

//	model_list_setup_links(model);
	models_process(model, model_setup_links);
	
	return model;
}

/**
@brief 	Writes molecular model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@return int			models written.
**/
int			write_model_molecule(Bstring& filename, Bmodel* model)
{
	int				i, n;
	Bstring			restype("UNK");
	Bmolgroup*		molgroup = NULL;
	Bmolgroup*		mg = NULL;
	Bmolgroup*		mgc = NULL;
	Bmolecule*		mol = NULL;
	Bresidue*		res = NULL;
	Batom*  		atom = NULL, *atom2;
	Bbond*			bond = NULL;
	
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;

	for ( n=i=0, mp = model; mp; mp = mp->next, n++ ) {
		mg = molgroup_init();
		if ( !molgroup ) molgroup = mg;
		else mgc->next = mg;
		mgc = mg;
//		mg->comment = model->comment().pre('\n') + "\nREMARK Model";
		mg->id = mp->identifier();
		mg->pointgroup = mp->symmetry();
		mol = molecule_add(&mg->mol, mg->id);
		res = residue_add(&mol->res, restype);
		res->num = 1;
		atom = NULL;
		for ( comp = mp->comp; comp; comp = comp->next ) {
			atom = atom_add(&atom, comp->type()->identifier().c_str());
			if ( !res->atom ) res->atom = atom;
			res->insert[0] = ' ';
			atom->num = stoi(comp->identifier());
			atom->coord = comp->location();
			atom->sel = comp->select();
			atom->q = comp->select();
			atom->b = comp->FOM();
		}
		bond = NULL;
		for ( link = mp->link; link; link = link->next ) {
			for ( comp = mp->comp, atom = res->atom; comp && comp != link->comp[0] && atom; 
				comp = comp->next, atom = atom->next ) ;
			for ( comp = mp->comp, atom2 = res->atom; comp && comp != link->comp[1] && atom2; 
				comp = comp->next, atom2 = atom2->next ) ;
			if ( atom && atom2 ) {
				bond = (Bbond *) bond_add(&bond, atom, atom2, link->length(), 1);
				if ( !mg->bond ) mg->bond = bond;
			} else {
				cerr << "Error: Problem with link " << link->comp[0]->identifier() << " to " << link->comp[1]->identifier() << endl;
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_molecule: atoms converted: " << i << endl;

	molgroup_list_write(filename, molgroup);

	molgroup_list_kill(molgroup);
	
	return  n;
}

