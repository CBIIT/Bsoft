/**
@file	rwmodel_xyz.cpp
@brief	Library routines to read and write Kirkland's xyz atomic model parameters
@author Bernard Heymann
@date	Created: 20220426
@date	Modified: 20220426
**/

#include "rwmodel.h"
#include "rwmodel_param.h"
#include "model_mol.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads Kirkland's xyz atomic model parameters.
@param 	*file_list		list of model parameter file names.
@param 	&paramfile		parameter file for atomic Z numbers.
@return Bmodel*			model parameters.
**/
Bmodel*		read_model_xyz(Bstring* file_list, Bstring& paramfile)
{
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Bcomptype*		ct = NULL;
	string			id("1"), type("XYZ");
	Bstring*		filename;
	ifstream		fmod;
	string			s;
	long			Z, natom(0);
	double			occ, wobble;
	Vector3<double>	size, loc;

	map<string,Bcomptype> 	atompar = read_atom_properties(paramfile);
	map<long,string>		Zsymbol;
	
	for ( auto a: atompar ) Zsymbol[a.second.index()] = a.first;

	for ( filename = file_list; filename; filename = filename->next ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << *filename << endl;
		fmod.open(filename->c_str());
		if ( fmod.fail() ) return  NULL;
		if ( model ) mp = model->add(id);
		else mp = model = new Bmodel(id);
		mp->model_type(type);
		mp->select(1);
		getline(fmod, s);
		fmod >> size[0] >> size[1] >> size[2];
		while ( fmod.good() ) {
			fmod >> Z >> loc[0] >> loc[1] >> loc[2] >> occ >> wobble;
			if ( Z > 0 ) {
				natom++;
				if ( comp ) comp = comp->add(natom);
				else comp = model->comp = new Bcomponent(natom);
				comp->location(loc);
				comp->density(occ);
				comp->FOM(wobble);
				comp->select(1);
				ct = model->add_type(Zsymbol[Z]);
				comp->type(ct);
				ct->index(Z);
				comp->description(Zsymbol[Z]);
			}
		}
		fmod.close();
	}

	return model;
}

/**
@brief 	Writes Kirkland's xyz atomic model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@return int			models written.
**/
int			write_model_xyz(Bstring& filename, Bmodel* model)
{
	int					n, Z(1);
	Bmodel*				mp = NULL;
	Bcomponent*			comp;
	Bstring				onename;
	double				occ(1), wobble(0);
	Vector3<double>		loc, size, off;
	Bstring				paramfile;

	map<string,Bcomptype> 	atompar = read_atom_properties(paramfile);

	model->calculate_bounds();
	size = model->maximum() - model->minimum();
	if ( size[1] > size[0] ) size[0] = size[1];
	size[0] = size[1] = 1.5*size[0];
	off = (size - model->maximum() - model->minimum())*0.5;
	
	ofstream		fmod;

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		if ( model->next )
			onename = filename.pre_rev('.') + Bstring(n+1, "_%04d.") + filename.post_rev('.');
		else
			onename = filename;
		fmod.open(onename.c_str());
		if ( fmod.fail() ) return  -1;
		fmod << command_line().c_str() << endl;
		fmod << setw(16) << size[0] << setw(16) << size[1] << setw(16) << size[2] << endl;
 		for ( comp = mp->comp; comp; comp = comp->next ) {
 			Z = atompar[component_element(comp)].index();
			loc = comp->location() + off;
//			occ = comp->density();
//			wobble = comp->FOM();
			fmod << setw(5) << Z << setw(14) << loc[0] << setw(14)
				<< loc[1] << setw(14) << loc[2] << setw(14) << occ
				<< setw(14) << wobble << endl;
		}
		fmod << setw(5) << -1 << endl;
		fmod.close();
	}
	
	return 0;
}
