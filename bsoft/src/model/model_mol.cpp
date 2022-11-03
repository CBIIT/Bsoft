/**
@file	model_mol.cpp
@brief	Library routines for processing molecular models
@author Bernard Heymann
@date	Created: 20220215
@date	Modified: 20220514
**/

#include "model_mol.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Determines the element from the component description.
@param 	*comp			component.
@return string			element.
**/
string		component_element(Bcomponent* comp)
{
		string		cel = comp->description().substr(0,2);
		
		cel[1] = tolower(cel[1]);
		
		if ( cel[1] == ' ' ) cel.resize(1);
		
		return cel;
}


/**
@brief 	Determines the element from the type identifier.
@param 	*comp			component.
@param	atompar			atom type parameters.
@return string			element, empty if not found.
**/
string		component_element(Bcomponent* comp, map<string,Bcomptype>& atompar)
{
/*	string		cel = comp->type()->identifier().substr(0,1);
	
	if ( atompar.find(cel) == atompar.end() ) {
		cel = comp->type()->identifier().substr(0,2);
		cel[1] = tolower(cel[1]);
	}
*/
	string		cel = component_element(comp);

	if ( atompar.find(cel) == atompar.end() ) {
		cerr << "Warning: Element " << cel << " not found!" << endl;
		cel.clear();
	}

	return cel;
}

/**
@brief 	Calculates the elemental composition.
@param 	*model			model.
@param	atompar			atom type parameters.
@return JSvalue			composition.
**/
JSvalue		model_elements(Bmodel* model, map<string,Bcomptype>& atompar)
{
	JSvalue			el(JSobject);

	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomptype*		ct;
	string			cel;

    for ( mp = model; mp; mp = mp->next ) {
		for( comp = mp->comp; comp; comp = comp->next ) {
			cel = component_element(comp, atompar);
			if ( cel.length() ) {
				if ( el.exists(cel) )
					el[cel] += 1;
				else
					el[cel] = 1;
			}
			if ( atompar.find(cel) != atompar.end() ) {
//				cout << "Before conversion:" << endl;
//				atompar[cel].show();
				ct = model->add_type(&atompar[cel]);
				comp->type(ct);
//				cout << "After conversion:" << endl;
//				ct->show();
//				cout << cel << tab << "Z = " << comp->type()->index() << endl;
			}
		}
	}

	if ( verbose )
		cout << el << endl;
	
	return el;
}
