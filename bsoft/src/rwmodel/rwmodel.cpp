/**
@file	rwmodel.cpp
@brief	Library routines to read and write model parameters
@author Bernard Heymann
@date	Created: 20060919
@date	Modified: 20221115
**/

#include "rwmodel.h"
#include "rwmodel_param.h"
#include "rwmodel_star.h"
#include "rwmodel_xml.h"
#include "rwmodel_cmm.h"
#include "rwmodel_bild.h"
#include "rwmodel_mol.h"
#include "rwmodel_pdb.h"
#include "rwmodel_cif.h"
#include "rwmodel_vega.h"
#include "rwmodel_xyz.h"
#include "model_links.h"
#include "file_util.h"
#include "linked_list.h"
#include "Color.h"
#include "utilities.h"

// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern string	command;		// Command line

// Internal function prototypes


Bmodel*		read_model(Bstring& filename)	
{
	return read_model(&filename);
}

Bmodel*		read_model(Bstring* file_list)	
{
	Bstring		paramfile;
	return read_model(file_list, paramfile);
}

Bmodel*		read_model(Bstring& filename, Bstring& paramfile)
{
	return read_model(&filename, paramfile);
}

/**
@brief 	Reads model parameters.
@param 	*file_list	list of model parameter file names.
@param 	&paramfile	parameter file.
@return Bmodel*		model parameters.
**/
Bmodel*		read_model(Bstring* file_list, Bstring& paramfile)
{
	Bstring			ext = file_list->extension();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_model: extension=" << ext << endl;

	Bmodel*			model = NULL;
	Bstring			path;

    if ( ext.contains("star") )
		model = read_model_star(file_list);
    else if ( ext.contains("xml") )
		model = read_model_xml(file_list);
    else if ( ext.contains("cmm") )
		model = read_model_chimera(file_list);
    else if ( ext.contains("bld") ||  ext.contains("bild") )
		model = read_model_bild(file_list);
	else if ( ext.contains("v3d") )
		model = read_model_vega(file_list);
	else if ( ext.contains("pdb") || ext.contains("ent") )
		model = read_model_pdb(file_list, paramfile);
	else if ( ext.contains("cif") )
		model = read_model_cif(*file_list, paramfile);
	else if ( ext.contains("xyz") )
		model = read_model_xyz(file_list, paramfile);
    else
		model = read_model_molecule(file_list, paramfile);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_model: id=" << model->identifier() << endl;

	if ( !model ) {
		cerr << "Error: File with extension " << ext << " not read!" << endl;
		error_show(file_list->c_str(), __FILE__, __LINE__);
	} else {
//		if ( !model->link ) model_list_setup_links(model);
		path = file_list->pre_rev('/');
		model_check(model, path);
		map<string,Bcomptype>	param = read_atom_properties(paramfile);
		model->update_component_types(param);
	}
	
//	if ( model->type ) model->type->show();
	
	return model;
}

/**
@brief 	Writes model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@return int			number of models.
**/
int			write_model(Bstring& filename, Bmodel* model)
{
	return write_model(filename, model, 0);
}

/**
@brief 	Writes model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	split		number of digits for writing multiple models.
@return int			number of models.
**/
int			write_model(Bstring& filename, Bmodel* model, int split)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model: starting" << endl;

	if ( !model ) return 0;
	
	int				n(0);

	model->comment() += command_line_time2();

	if ( filename.empty() ) {
		cerr << "Error: No model filename!" << endl;
		return  -1;
	}
	
	Bstring			path(0);
	model_check(model, path);
	
	Bstring			ext = filename.extension();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model: extension=" << ext << endl;

	if ( verbose & VERB_LABEL )
	    cout << "Writing file:                   " << filename << endl;

//	Bstring			path = filename.pre_rev('/');
//	model_check(model, path);
	
    if ( ext.contains("star") )
		n = write_model_star(filename, model, split);
    else if ( ext.contains("xml") )
		n = write_model_xml(filename, model);
    else if ( ext.contains("cmm") )
		n = write_model_chimera(filename, model, split);
    else if ( ext.contains("bld") ||  ext.contains("bild") )
		n = write_model_bild(filename, model, split);
    else if ( ext.contains("v3d") )
		n = write_model_vega(filename, model, split);
	else if ( ext.contains("pdb") )
		n = write_model_pdb(filename, model, split);
	else if ( ext.contains("cif") )
		n = write_model_cif(filename, model);
	else if ( ext.contains("xyz") )
		n = write_model_xyz(filename, model, split);
    else
		n = write_model_molecule(filename, model);
	
	if ( n < 0 ) {
		cerr << "Error: File with extension " << ext << " not written!" << endl;
		error_show(filename.c_str(), __FILE__, __LINE__);
	}
	
	return  n;
}

/**
@brief 	Adds a model to a linked list.
@param 	**model		model list.
@param 	id			model identifier.
@return Bmodel*		new model.

	The function allocates memory for a new model structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bmodel*		model_add(Bmodel** model, string id)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_add: id=" << id << endl;
	Bmodel* 		this_model = *model;
	if ( !this_model ) return this_model = new Bmodel(id);
	return this_model->add(id);
}

/**
@brief 	Adds a component to a linked list.
@param	**comp		component list.
@param 	&id			component identifier.
@return Bcomponent	new component.


	The function allocates memory for a new component structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bcomponent*	component_add(Bcomponent** comp, string& id)
{
	Bcomponent* 		comp1 = *comp;

	if ( comp1 ) comp1 = comp1->add(id);
	else *comp = comp1 = new Bcomponent(id);
	
	return comp1;
}

Bcomponent*	component_add(Bcomponent** comp, long number)
{
	string			id(to_string(number));
	return component_add(comp, id);
}

/**
@brief 	Adds a component link to a linked list.
@param	**link		component link list.
@param	*comp1		first component.
@param	*comp2		second component.
@param 	length		length of link.
@param 	radius		display radius.
@return Bcomponent	new component.

	The function allocates memory for a new link structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Blink*		link_add(Blink** link, Bcomponent* comp1, Bcomponent* comp2, double length, double radius)
{
	Blink*			link1 = *link;

	// Find if the link exists
	if ( comp1->find_link_exists(comp2) && comp2->find_link_exists(comp1) ) {
		link1 = link1->find(comp1, comp2);
		if ( verbose & VERB_FULL )
			cerr << "Error: " << comp1->identifier() << " already linked to " << comp2->identifier() << endl;
		return link1;
	}
	
	comp1->add_link(comp2);
	comp2->add_link(comp1);

	if ( link1 ) link1 = link1->add(comp1, comp2, length, radius);
	else link1 = new Blink(comp1, comp2, length, radius);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG link_add: " << comp1->identifier() << " - " << comp2->identifier() << endl;

	return link1;
}

Blink*		link_add(Blink** link, Bcomponent* comp1, Bcomponent* comp2)
{
	return link_add(link, comp1, comp2, 0, 0);
}

/**
@brief 	Sets all the map file names of selected models.
@param 	*model		model parameters.
@param 	&mapfile	map file name.
@return int			0.
**/
int			model_set_map_filenames(Bmodel* model, Bstring& mapfile)
{
	Bmodel*		mp;
	
	for ( mp = model; mp; mp = mp->next )
		if ( mp->select()) mp->mapfile(mapfile.str());
	
	return 0;
}

/**
@brief 	Reset the component types.
@param 	*model		model.
@param 	&set_type	component type.
@return int			number of models.

	Sets all the component types to the given string.

**/
int			model_set_type(Bmodel* model, Bstring& set_type)
{
	if ( !model ) return 0;
	
	int				n;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomptype*		ct;

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		comp_type_list_kill(mp->type);
		mp->type = NULL;
		ct = mp->add_type(set_type);
		for ( comp = mp->comp; comp; comp = comp->next )
			comp->type(ct);
	}
	
	return  n;
}

/**
@brief 	Change a component type name.
@param 	*model			model.
@param 	&change_type	component type.
@return int				number of models.

	Sets all the component types to the given string.

**/
int			model_change_type(Bmodel* model, Bstring& change_type)
{
	if ( !model ) return 0;
	
	int				n;
	Bmodel*			mp;
	Bcomptype*		ct;
	
	Bstring			ot(change_type.pre(',')), nt(change_type.post(','));

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		for ( ct = mp->type; ct; ct = ct->next )
			if ( ct->identifier() == ot.str() ) ct->identifier(nt.str());
	}
	
	return  n;
}

/**
@brief 	Checks model properties.
@param 	*model		model.
@param	path		search path to find map files.
@return int			0.
**/
int			model_check(Bmodel* model, Bstring path)
{
	if ( !model ) return 0;
	
//	long			nmod(0), ncomp(0), nlink(0);
	Bmodel*			mp = NULL;
//	Bcomptype*		ct = NULL;
//	Bcomponent*		comp = NULL;
//	Blink*			link = NULL;
//	Bstring			nutype("UNK");

 	
	for ( mp = model; mp; mp = mp->next ) {
		if ( path.length() ) {
			if ( mp->mapfile().empty() ) mp->mapfile("?");
			else mp->mapfile(find_file(mp->mapfile().c_str(), path).str());
		}
		mp->check();
/*		for ( comp = mp->comp; comp; comp = comp->next, ncomp++ ) {
			comp->view().normalize();
			if ( comp->radius() < 1e-6 ) comp->radius(1);
			if ( !comp->type() ) {
				cerr << "Component " << comp->identifier() << " has no type assigned" << endl;
				comp->type(mp->add_type(nutype));
			}
		}
		for ( ct = mp->type; ct; ct = ct->next ) {
			if ( ct->file_name().empty() ) ct->file_name("?");
			if ( ct->mass() < 1 ) ct->mass(1);
		}
		for ( link = mp->link; link; link = link->next, nlink++ ) {
			if ( link->length() < 1e-6 )
				link->length(link->comp[0]->location().distance(link->comp[1]->location()));
			if ( link->radius() < 1e-6 )
				link->radius(0.1*link->length());
		}
		mp->update_type_counts();*/
	}


	
	return 0;
}

/**
@brief 	Copies a model.
@param 	*model		model parameters.
@return Bmodel*		model copy.

	Copies all the models in a list.

**/
Bmodel*		model_list_copy(Bmodel* model)
{
	if ( !model ) return  NULL;
	
	Bmodel*			model_nu = NULL;
	Bmodel*			mp;
	Bmodel*			mp_new;
	
	for ( mp = model; mp; mp = mp->next ) {
		if ( model_nu ) {
			mp_new->next = mp->copy();
			mp_new = mp_new->next;
		} else {
			model_nu = mp_new = mp->copy();
		}
	}
	
	return model_nu;
}

/**
@brief 	Copies a model.
@param 	*model		model parameters.
@return Bmodel*		model copy.

	Copies only the first model in the list.

**/
Bmodel*		model_copy(Bmodel* model)
{
	if ( !model ) return  NULL;
	
	int				i;
	
	Bcomptype*		type;
	Bcomponent*		comp;
	Blink*			link;
	Bpolygon*		poly;

	Bcomptype*		type_nu = NULL;
	Bcomponent*		comp_new = NULL;
	Blink*			link_nu = NULL;
	Bpolygon*		poly_nu = NULL;
	
	Bmodel*			model_nu = new Bmodel(model->identifier());
	
	model_nu->model_type(model->model_type());
	model_nu->comment(model->comment());
	model_nu->symmetry(model->symmetry());
	model_nu->mapfile(model->mapfile());
	model_nu->image_number(model->image_number());
	model_nu->select(model->select());
	model_nu->FOM(model->FOM());
	model_nu->handedness(model->handedness());
	
	for ( type = model->type; type; type = type->next ) {
		type_nu = model_nu->add_type(type->identifier(), type->file_name(), type->image_number());
		type_nu->component_count(type->component_count());
		type_nu->mass(type->mass());
		type_nu->FOM(type->FOM());
		type_nu->select(type->select());
	}
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		comp_new = model_nu->add_component(comp);
		comp_new->type(model_nu->add_type(comp->type()->identifier()));
	}
	
	for ( link = model->link; link; link = link->next ) {
		link_nu = (Blink *) add_item((char **) &link_nu, sizeof(Blink));
		if ( !model_nu->link ) model_nu->link = link_nu;
		for ( comp = model->comp, comp_new = model_nu->comp; comp; comp = comp->next, comp_new = comp_new->next ) {
			if ( link->comp[0] == comp ) link_nu->comp[0] = comp_new;
			if ( link->comp[1] == comp ) link_nu->comp[1] = comp_new;
		}
		link_nu->length(link->length());
		link_nu->radius(link->radius());
		link_nu->color(link->color());
	}
	
	for ( poly = model->poly; poly; poly = poly->next ) {
		poly_nu = (Bpolygon *) add_item((char **) &poly_nu, sizeof(Bpolygon));
		if ( !model_nu->poly ) model_nu->poly = poly_nu;
		poly_nu->normal(poly->normal());
		poly_nu->closed(poly->closed());
		for ( comp = model->comp, comp_new = model_nu->comp; comp; comp = comp->next, comp_new = comp_new->next ) {
			for ( i=0; poly->comp[i]; i++ )
				if ( poly->comp[i] == comp ) poly_nu->comp[i] = comp_new;
		}
	}

	model_setup_links(model_nu);
	
	return model_nu;
}

/**
@brief 	Deallocates memory for a list of components.
@param	*comp	component list.
@return int					total number of components.
**/
int			component_list_kill(Bcomponent* comp)
{
	int				n(0);
	Bcomponent*		c = NULL;
	Bcomponent*		c2 = NULL;
	
	for ( c = comp; c; ) {
		c2 = c->next;
		delete c;
		c = c2;
		n++;
	}
	
	return  n;
}

/**
@brief 	Deallocates memory for a list of component types.
@param	*type		component type list.
@return int					total number of component types.
**/
int			comp_type_list_kill(Bcomptype* type)
{
	int				n(0);
	Bcomptype*		ct = NULL;
	Bcomptype*		ct2 = NULL;
	
	for ( ct = type; ct; ) {
		ct2 = ct->next;
		delete ct;
		ct = ct2;
		n++;
	}
	
	return  n;
}

/**
@brief 	Deallocates memory for a list of component links.

	Only the first model in the list is processed.

@param 	*model		model.
@return int					total number of component links.
**/
int			model_link_list_kill(Bmodel* model)
{
	if ( !model ) return 0;
	
	int				i, n(0);
	Bcomponent*		comp;
	Blink*			link = NULL;
	Blink*			link2 = NULL;
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		for ( i=0; i<comp->link.size(); i++ ) {
			comp->link[i] = NULL;
			comp->flag[i] = 0;
		}
	}
	
	for ( link = model->link; link; ) {
		link2 = link->next;
		delete link;
		link = link2;
		n++;
	}
	
	return  n;
}

/**
@brief 	Deletes a link.

	The link in the model link list is removed.
	The associated references to the link in the component link arrays
	are removed and the link arrays reorganized.

@param	**link_list	pointer to list of links.
@param	*comp	one component in the link.
@param 	i				index for second component in link array of first component.
@return Bdistmat*			new distance matrix structure.
**/
int			link_kill(Blink** link_list, Bcomponent* comp, int i)
{
	int			j, k;
	Blink*		link, *plink;
	Bcomponent*	comp2 = comp->link[i];
	
	j = (*link_list)->index(comp2, comp);
	
	if ( j < 0 ) {
		cerr << "Error: Component " << comp->identifier() << " is not linked to component " << comp2->identifier() << endl << endl;
		return  -1;
	}

	if ( verbose & VERB_FULL )
		cout << "Deleting: " << comp->identifier() << " (" << comp->select() << ") " << comp2->identifier() << " (" << comp2->select() << ")" << endl;
	
	for ( k=i; k+1<comp->link.size() && comp->link[k+1]; k++ ) {
		comp->link[k] = comp->link[k+1];
		comp->flag[k] = comp->flag[k+1];
	}
	comp->link[k] = NULL;
	comp->flag[k] = 0;
	
	for ( k=j; k+1<comp2->link.size() && comp2->link[k+1]; k++ ) {
		comp2->link[k] = comp2->link[k+1];
		comp2->flag[k] = comp2->flag[k+1];
	}
	comp2->link[k] = NULL;
	comp2->flag[k] = 0;
	
	for ( link = plink = *link_list; link;  ) {
		if ( ( link->comp[0] == comp && link->comp[1] == comp2 ) ||
				( link->comp[1] == comp && link->comp[0] == comp2 ) ) {
			if ( link == plink ) {
				*link_list = plink = link->next;
				delete link;
				link = plink;
			} else {
				plink->next = link->next;
				delete link;
				link = plink->next;
			}
		} else {
			plink = link;
			link = link->next;
		}
	}
	
	return 0;
}

/**
@brief 	Deletes a link.

	The link in the model link list is removed.
	The associated references to the link in the component link arrays
	are removed and the link arrays reorganized.

@param	**link_list	pointer to list of links.
@param	*comp	one component in the link.
@param	*comp2	second component in the link.
@return Bdistmat*			new distance matrix structure.
**/
int			link_kill(Blink** link_list, Bcomponent* comp, Bcomponent* comp2)
{
	int			i = (*link_list)->index(comp, comp2);
	
	if ( !comp->link[i] ) {
		cerr << "Error: Component " << comp->identifier() << " is not linked to component " << comp2->identifier() << endl << endl;
		return  -1;
	}

	link_kill(link_list, comp, i);
	
	return 0;
}

/**
@brief 	Deallocates memory for a list of polygons.
@param 	*poly		polygon list.
@return int					total number of polygons.
**/
int			poly_list_kill(Bpolygon* poly)
{
	int				n(0);
	Bpolygon*		p = NULL;
	Bpolygon*		p2 = NULL;
	
	for ( p = poly; p; ) {
		p2 = p->next;
		delete p;
		p = p2;
		n++;
	}
	
	return  n;
}

/**
@brief 	Deallocates memory for links associated with a component.
@param	*comp	associated componet.
@param	**link		pointer to component link list.
@return int					number of component links deallocated.
**/
int			comp_associated_links_kill(Bcomponent* comp, Blink** link)
{
	int				n(0);
	Blink*			l = *link;
	Blink*			lp = NULL;

	while ( l ) {
		if ( l->comp[0] == comp || l->comp[1] == comp ) {
			if ( lp ) {
				lp->next = l->next;
				delete l;
				l = lp->next;
			} else {
				*link = l->next;
				delete l;
				l = *link;
			}
			n++;
		} else {
			lp = l;
			l = l->next;
		}
	}
	
	return  n;
}

/**
@brief 	Deallocates all memory in the list.
@param 	*model		model parameters.
@return int					0.
**/
int 		model_kill(Bmodel* model)
{
	Bmodel*			mp = NULL;
	Bmodel*			mp2 = NULL;
	
	for ( mp = model; mp; ) {
		mp2 = mp->next;
		model_link_list_kill(mp);
		component_list_kill(mp->comp);
		comp_type_list_kill(mp->type);
		poly_list_kill(mp->poly);
		delete mp;
		mp = mp2;
	}
	
	return 0;
}

/**
@brief 	Associates a model file with a component type.
@param 	*model			the model.
@param 	&associate_type	component type.
@param 	&associate_file	component file name.
@return int				number of types associated.

	Model files can be coordinates or maps.

**/
int			model_associate(Bmodel* model, Bstring& associate_type, Bstring& associate_file)
{
	if ( !model ) return 0;
	
	int				n(0);
	Bmodel*			mp = NULL;
	Bcomptype*		ct = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Associating component " << associate_type << " with file " << associate_file << endl << endl;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( ct = mp->type; ct; ct = ct->next ) if ( ct->select() ) {
			if ( ct->identifier() == associate_type.str() ) {
				ct->file_name(associate_file.str());
				n++;
			}
		}
	}
	
	return  n;
}

/**
@brief 	Associates a mass with a component type.
@param 	*model			model list.
@param 	&associate_type	component type.
@param 	mass			component type mass.
@return int				number of types associated.
**/
int			model_associate_mass(Bmodel* model, Bstring& associate_type, double mass)
{
	if ( !model ) return 0;
	
	int				n(0);
	Bmodel*			mp = NULL;
	Bcomptype*		ct = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Associating component " << associate_type << " with mass = " << mass << endl << endl;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( ct = mp->type; ct; ct = ct->next ) if ( ct->select() ) {
			if ( ct->identifier() == associate_type.str() ) {
				ct->mass(mass);
				n++;
			}
		}
	}
	
	return  n;
}

/**
@brief 	Sets the filenames of all selected component types to the given string.
@param 	*model		model parameters.
@param 	&filename	component file name.
@return int			number of component types set.

	The image numbers are sequentially set as well.

**/
int			model_set_comptype_filenames(Bmodel* model, Bstring& filename)
{
	if ( !model ) return 0;
	
	int				n(0);
	Bmodel*			mp;
	string			fn(filename.str());
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() )
		n += mp->set_type_filenames(fn);
	
	return  n;
}


/**
@brief 	Set the display radius for all components to a specific value.
@param 	*model		model parameters.
@param 	comprad		component display radius.
@return long			number of components selected.
**/
long		model_set_component_radius(Bmodel* model, double comprad)
{
	if ( !model ) return 0;

	int				n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() )
		n += mp->set_component_radius(comprad);

	return n;
}

