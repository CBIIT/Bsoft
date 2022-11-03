/**
@file	model_select.cpp
@brief	Library routines used for model and component selection
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20220323
**/

#include "rwmodel.h"
#include "model_poly.h"
#include "model_views.h"
#include "model_compare.h"
#include "model_util.h"
#include "rwimg.h"
#include "qsort_functions.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates selection statistics.
@param 	*model		model parameters.
@return long		number of components selected.

	The FOM is assumed to be a value from 0 to 1.

**/
long		model_selection_stats(Bmodel* model)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_selection_stats: " << model->identifier() << endl;
	
	long			i, h;
	long			nmod(0), ncomp(0), nlink(0), npoly(0);
	long			nmods(0), ncomps(0), nlinks(0), npolys(0);
	long			n, nmt(0), nct(0), nv[1000];
	int				val[MAXLINK];
	double			fomavg(0), fomstd(0), fmin(1e30), fmax(-1e30);
	Vector3<double>	cmin(1e10, 1e10, 1e10), cmax(-1e10, -1e10, -1e10), com;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Bcomptype*		type = NULL;
	Blink*			link = NULL;
	Bpolygon*		poly = NULL;
	Bstring*		modtypelist = NULL;
	Bstring*		modtype = NULL;
	Bstring*		comptypelist = NULL;
	Bstring*		comptype = NULL;
	Bstring*		temptype = NULL;
	Bstring*		typefilelist = NULL;
	Bstring*		typefile = NULL;
	
	for ( i=0; i<MAXLINK; i++ ) val[i] = 0;
	for ( i=0; i<1000; i++ ) nv[i] = 0;

	for ( mp = model; mp; mp = mp->next, nmod++ ) {
		for ( comp = mp->comp; comp; comp = comp->next ) ncomp++;
		for ( link = mp->link; link; link = link->next ) nlink++;
		for ( poly = mp->poly; poly; poly = poly->next ) npoly++;
	}
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_selection_stats: nmod=" << nmod << " ncomp=" << ncomp
			<< " nlink=" << nlink << " npoly=" << npoly << endl;
	
	// Compile temporary lists of model and component types
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		nmods++;
		if ( mp->model_type().length() < 1 ) mp->model_type("Unknown");
		for ( temptype = modtypelist; temptype && *temptype != mp->model_type(); temptype = temptype->next ) ;
		if ( !temptype ) {
			modtype = string_add(&modtype, mp->model_type().c_str());
			if ( !modtypelist ) modtypelist = modtype;
			nmt += 3;
		}
		for ( type = mp->type; type; type = type->next ) {
			for ( temptype = comptypelist; temptype && *temptype != type->identifier(); temptype = temptype->next ) ;
			if ( !temptype ) {
				comptype = string_add(&comptype, type->identifier().c_str());
				if ( !comptypelist ) comptypelist = comptype;
				typefile = string_add(&typefile, type->file_name().c_str());
				if ( !typefilelist ) typefilelist = typefile;
				nct++;
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_selection_stats: nmods=" << nmods << " nmt=" << nmt << " nct=" << nct << endl;
	
	if ( nmods < 1 ) {
		cerr << "No models selected!" << endl << endl;
		return 0;
	}
	
	if ( nmt < 3 ) nmt = 3;
	if ( nct < 1 ) nct = 1;
	
	// Sort the temporary lists
	if ( modtypelist ) string_sort(modtypelist, 0, 1);
	if ( comptypelist ) string_sort(comptypelist, 0, 1);
	
	// Rearrange associated file names
	for ( temptype = comptypelist, typefile = typefilelist; temptype && typefile; temptype = temptype->next, typefile = typefile->next ) {
		for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
			for ( type = mp->type; type && *temptype != type->identifier(); type = type->next ) ;
			if ( type ) *typefile = type->file_name() ;
		}
	}
	
	int*			nvert = new int[nmt];
	int*			nmodtype = new int[nmt];
	double*			modfom = new double[nmt];
	int*			ncomptype = new int[nct];
	double*			compfom = new double[nct];
	
	for ( i=0; i<nmt; i++ ) modfom[i] = nvert[i] = nmodtype[i] = 0;
	for ( i=0; i<nct; i++ ) compfom[i] = ncomptype[i] = 0;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG model_selection_stats: model id=" << mp->identifier() << endl;
		for ( n=0, comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) n++;
		if ( n < 1000 ) nv[n]++;
		for ( i=0, temptype = modtypelist; temptype && *temptype != mp->model_type(); temptype = temptype->next, i+=3 ) ;
		i += mp->handedness() + 1;
		if ( i < nmt ) {
			nvert[i] = n;
			nmodtype[i]++;
			modfom[i] += mp->FOM();
		}
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG model_selection_stats: component id=" << comp->identifier() << " ncomps=" << ncomps << endl;
			for ( i=0, temptype = comptypelist; temptype && *temptype != comp->type()->identifier(); temptype = temptype->next, i++ ) ;
			if ( i < nct ) {
				ncomptype[i]++;
				compfom[i] += comp->FOM();
			}
			com += comp->location();
			cmin = cmin.min(comp->location());
			cmax = cmax.max(comp->location());
			if ( fmin > comp->FOM() ) fmin = comp->FOM();
			if ( fmax < comp->FOM() ) fmax = comp->FOM();
			fomavg += comp->FOM();
			fomstd += comp->FOM()*comp->FOM();
			ncomps++;
			for ( i=0; i<comp->link.size() && comp->link[i]; i++ ) ;
			if ( i < comp->link.size() ) val[i]++;
		}
		for ( link = mp->link; link; link = link->next )
			if ( link->select() ) nlinks++;
		for ( poly = mp->poly; poly; poly = poly->next )
			if ( poly->select() ) npolys++;
	}
	
	if ( ncomps ) {
		com /= ncomps;
		fomavg /= ncomps;
		fomstd = fomstd/ncomps - fomavg*fomavg;
		if ( fomstd > 0 ) fomstd = sqrt(fomstd);
		else fomstd = 0;
	}
	
	if ( verbose ) {
		cout << "Models:                         " << nmod << " (" << nmods << ")" << endl;
		cout << "Components:                     " << ncomp << " (" << ncomps << ")" << endl;
		cout << "Links:                          " << nlink << " (" << nlinks << ")" << endl;
		cout << "Polygons:                       " << npoly << " (" << npolys << ")" << endl << endl;
		if ( ncomps ) {
			cout << "Center-of-mass:                 " << com << endl;
			cout << "Coordinate minima:              " << cmin << endl;
			cout << "Coordinate maxima:              " << cmax << endl;
			cout << "Selected FOM minimum & maximum: " << fmin << " " << fmax << endl;
			cout << "Selected FOM average & stdev:   " << fomavg << " " << fomstd << endl << endl;
		}
		cout << "Model sizes:\nVert\tCount" << endl;
		for ( i=n=0; i<1000; i++ ) if ( nv[i] ) {
			cout << i << tab << nv[i] << endl;
			n += nv[i];
		}
		cout << "Total\t" << n << endl << endl;
		cout << "Model types:\nType\tVert\tHand\tCount\tFOM" << endl;
		for ( i=0, modtype = modtypelist; i<nmt && modtype; modtype = modtype->next, i+=3 )
			for ( h=i; h<i+3; h++ )
				if ( nmodtype[h] ) cout << *modtype << tab << nvert[h] << tab
					<< h%3 - 1 << tab << nmodtype[h] << tab << modfom[h]/nmodtype[h] << endl;
		cout << endl;
		cout << "Component types:\nType\tHand\tCount\tFOM\tFile" << endl;
		for ( i=0, comptype = comptypelist, typefile = typefilelist;
				i<nct && comptype && typefile; comptype = comptype->next, typefile = typefile->next, i++ ) {
			if ( ncomptype[i] ) compfom[i] /= ncomptype[i];
			cout << *comptype << tab << component_hand(*comptype) << tab <<
				ncomptype[i] << tab << compfom[i] << tab << *typefile << endl;
		}
		cout << endl;
		cout << "Valence\tCount" << endl;
		for ( i=0; i<MAXLINK; i++ ) if ( val[i] ) cout << i << tab << val[i] << endl;
		cout << endl;
	}
	
	string_kill(modtypelist);
	string_kill(comptypelist);
	string_kill(typefilelist);
	delete[] nvert;
	delete[] nmodtype;
	delete[] modfom;
	delete[] ncomptype;
	delete[] compfom;
	
	return ncomps;
}

/**
@brief 	Shows the distrubution of selections.
@param 	*model	model parameters.
@return long		number of selection levels.


**/
long 		model_show_selection(Bmodel* model)
{
	long			i;
	vector<long>	lev;
	Bmodel*			mp;
	Bcomponent*		comp;
	
	for ( mp = model; mp; mp = mp->next )
		for ( comp = mp->comp; comp; comp = comp->next ) {
			while ( lev.size() <= comp->select() ) lev.push_back(0);
			lev[comp->select()]++;
		}
	
	cout << "Component selection distribution:" << endl;
	cout << "Level\tCount" << endl;
	for ( i = 0; i<lev.size(); ++i )
		cout << i << tab << lev[i] << endl;
	cout << endl;
	
	return lev.size();
}

/**
@brief 	Selects models, components and component types.
@param 	*model		model parameters.
@param 	&selstr		selection string.
@return long		number of selections made.

	The selection string can have one of the following formats:
		#model@component
		#model%comp_type
		^model_type@component
		^model_type%comp_type
	Only elements originally selected is considered, except where the "." is used.

**/
long 		model_select(Bmodel* model, Bstring& selstr)
{
	long			n(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomptype*		type;
	Blink*			link;
	
	Bstring			mod_id;
	Bstring			mod_type_id;
	Bstring			comp_id;
	Bstring			comp_type_id;
	Bstring*		list, *one;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_select: selstr=" << selstr << endl;

	if ( selstr[0] == '#' ) {
		if ( selstr.contains("@") ) {
			mod_id = selstr.substr(1, selstr.index('@')-1);
			comp_id = selstr.post('@');
		} else if ( selstr.contains("%") ) {
			mod_id = selstr.substr(1, selstr.index('%')-1);
			comp_type_id = selstr.post('%');
		} else {
			mod_id = selstr.substr(1, selstr.length());
		}
	} else if ( selstr[0] == '^' ) {
		if ( selstr.contains("@") ) {
			mod_type_id = selstr.substr(1, selstr.index('@')-1);
			comp_id = selstr.post('@');
		} else if ( selstr.contains("%") ) {
			mod_type_id = selstr.substr(1, selstr.index('%')-1);
			comp_type_id = selstr.post('%');
		} else {
			mod_type_id = selstr.substr(1, selstr.length());
		}
	} else if ( selstr[0] == '@' ) {
		comp_id = selstr.substr(1, selstr.length());
	} else if ( selstr[0] == '%' ) {
		comp_type_id = selstr.substr(1, selstr.length());
	}
	
	if ( verbose ) {
		cout << "Selection identifiers:" << endl;
		if ( mod_id.length() )
			cout << "Model:                          " << mod_id << endl;
		if ( mod_type_id.length() )
			cout << "Model type:                     " << mod_type_id << endl;
		if ( comp_id.length() )
			cout << "Component:                      " << comp_id << endl;
		if ( comp_type_id.length() )
			cout << "Component type:                 " << comp_type_id << endl;
		cout << endl;
	}

	if ( mod_id == "." || mod_type_id == "." ) {
		for ( mp = model; mp; mp = mp->next ) mp->select(1);
	} else {
		if ( mod_id.length() ) {
			list = mod_id.split(",");
			for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
				for ( one = list; one && mp->identifier() != one->str(); one = one->next ) ;
				if ( !one ) mp->select(0);
			}
			string_kill(list);
		}
		if ( mod_type_id.length() ) {
			list = mod_type_id.split(",");
			for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
				for ( one = list; one && mp->model_type() != one->str(); one = one->next ) ;
				if ( !one ) mp->select(0);
			}
			string_kill(list);
		}
	}
	
	if ( comp_type_id == "." ) {
		for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
			for ( type = mp->type; type; type = type->next ) type->select(1);
			for ( comp = mp->comp; comp; comp = comp->next ) comp->select(1);
		}
	} else if ( comp_type_id.length() ) {
		list = comp_type_id.split(",");
		for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
			for ( type = mp->type; type; type = type->next ) if ( type->select() ) {
				for ( one = list; one && type->identifier() != one->str(); one = one->next ) ;
				if ( !one ) type->select(0);
			}
			for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
				comp->select(comp->type()->select());
//				for ( one = list; one && comp->type()->identifier() != one->str(); one = one->next ) ;
//				if ( !one ) comp->select(0);
			}
		}
		string_kill(list);
	}

	if ( comp_id == "." ) {
		for ( mp = model; mp; mp = mp->next ) if ( mp->select() )
			for ( comp = mp->comp; comp; comp = comp->next ) comp->select(1);
	} else if ( comp_id.length() ) {
		list = comp_id.split(",");
		for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
			for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
				for ( one = list; one && comp->identifier() != one->str(); one = one->next ) ;
				if ( !one ) comp->select(0);
			}
		}
		string_kill(list);
	}

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) n++;
		for ( link = mp->link; link; link = link->next )
			if ( link->comp[0]->select() && link->comp[1]->select() ) link->select(1);
	}

	return n;
}

/**
@brief 	Resets the selection to all models.
@param 	*model		model parameters.
@param	number		selection number to select.
@return long			number of components selected.
**/
long		model_select(Bmodel* model, long number)
{
	long			nsel(0), ntot(0);
	Bmodel*			mp = NULL;
	Bcomponent*		comp;

	if ( verbose )
		cout << "Selecting components with selection number " << number << endl;

	for ( mp = model; mp; mp = mp->next )
		for ( comp = mp->comp; comp; comp = comp->next ) ntot++;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			if ( comp->select() != number ) comp->select(0);
			nsel += ( comp->select() > 0 );
		}
	}
	
	if ( verbose )
		cout << "Number of components selected:  " << nsel << " (" << nsel*100.0/ntot << " %)" << endl << endl;

	return nsel;
}

/**
@brief 	Resets the selection to all models.
@param 	*model		model parameters.
@return long		number of models selected.
**/
long		model_select_all(Bmodel* model)
{
	long			nsel(0);
	Bmodel*			mp = NULL;

	for ( nsel=0, mp = model; mp; mp = mp->next, nsel++ ) mp->select(1);
	
	return nsel;
}

/**
@brief 	Resets the selection to all unknown models.
@param 	*model		model parameters.
@return long		number of models selected.
**/
long		model_select_unknowns(Bmodel* model)
{
	long			nsel(0);
	Bmodel*			mp;
	
	for ( nsel=0, mp = model; mp; mp = mp->next ) {
		mp->select(0);
		if ( mp->model_type().length() < 1 ) {
			mp->select(1);
			nsel++;
		}
	}
	
	return nsel;
}

/**
@brief 	Resets the selection to all components.
@param 	*model		model parameters.
@return long		number of components selected.

	Only one model is modified.

**/
long		model_reset_selection(Bmodel* model)
{
	long			nsel(0);
	Bcomptype*		type;
	Bcomponent*		comp;
	Blink*			link;
	Bpolygon*		poly;

	if ( verbose & VERB_FULL )
		cout << "Model " << model->identifier() << ": reset selection" << endl;
	
	model->select(1);
	for ( type = model->type; type; type = type->next )
		type->select(1);
	for ( comp = model->comp; comp; comp = comp->next, nsel++ )
		comp->select(1);
	for ( link = model->link; link; link = link->next )
		link->select(1);
	for ( poly = model->poly; poly; poly = poly->next )
		poly->select(1);
	
	return nsel;
}

/**
@brief 	Unsets the selection to all components but not models.
@param 	*model		model parameters.
@return long		number of components deselected.
**/
long		model_unset_selection(Bmodel* model)
{
	int				nsel(0);
	Bcomptype*		type;
	Bcomponent*		comp;
	Blink*			link;
	Bpolygon*		poly;

	for ( type = model->type; type; type = type->next )
		type->select(0);
	for ( comp = model->comp; comp; comp = comp->next, nsel++ )
		comp->select(0);
	for ( link = model->link; link; link = link->next )
		link->select(0);
	for ( poly = model->poly; poly; poly = poly->next )
		poly->select(0);
	
	return nsel;
}

/**
@brief 	Inverts the selection of all components.
@param 	*model		model parameters.
@return long		number of components selected.
**/
long		model_invert_selection(Bmodel* model)
{
	int				nsel(0);
	Bcomponent*		comp;

	for ( comp = model->comp; comp; comp = comp->next, nsel++ )
		if ( comp->select() ) comp->select(0);
		else {
			comp->select(1);
			nsel++;
		}
	
	return nsel;
}

/**
@brief 	Selects sets of components, each set with the same size.
@param 	*model		parameter structure with all parameters.
@param 	size		number of components in each set.
@param 	flag		flag to not count across model boundaries.
@return long		number of components selected.

	Sets up sets of components, each set identified as a number in the
	selection array.

**/
long		model_select_sets(Bmodel* model, int size, int flag)
{
	if ( !model || size <= 0 ) return 0;
	
	long			nsel(0), count(0), number(0), ntot(0);
	
	Bmodel*			mp;
	Bcomponent*		comp;
	
	if ( verbose & VERB_LABEL )
		cout << "Generating sets of " << size << " components each" << endl;
	
	for ( nsel=0, mp = model; mp; mp = mp->next ) {
		if ( flag && count ) {
			count = 0;
			number++;
		}
		for ( comp = mp->comp; comp; comp = comp->next, ntot++ ) {
			if ( comp->select() > 0 ) {
				if ( count >= size ) {
					number++;
					count = 0;
				}
				if ( number < 1 ) number = 1;
				comp->select(number);
				count++;
				nsel++;
			}
		}
	}
	
	if ( verbose & VERB_LABEL ) {
		cout << "Number of sets generated:       " << number << endl;
		cout << "Number of components selected:  " << nsel << " (" << nsel*100.0/ntot << " %)" << endl << endl;
	}
	
	return nsel;
}

/**
@brief 	Selects components within bounds for a list of models.
@param 	*model			list of models.
@param 	&start			minimum coordinates.
@param 	&end			maximum coordinates.
@return long				number of components selected.

**/
long		models_select_within_bounds(Bmodel* model, Vector3<double>& start, Vector3<double>& end)
{
	long		nsel(0);
	
	for ( Bmodel* m = model; m; m = m->next )
		nsel += m->select_within_bounds(start, end);
	
	return nsel;
}


/**
@brief 	Selects models within a range of the number of components.
@param 	*model		model parameters.
@param 	ncomp_min	minimum number of components.
@param 	ncomp_max	maximum number of components.
@return long		number of models selected.
**/
long		model_select_number_of_components(Bmodel* model, int ncomp_min, int ncomp_max)
{
	if ( ncomp_min < 1 ) return 0;
	if ( ncomp_min > ncomp_max ) ncomp_max = ncomp_min;
	
	long			n, nsel(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( n=0, comp = mp->comp; comp; comp = comp->next ) n++;
		if ( n >= ncomp_min && n <= ncomp_max ) {
			mp->select(1);
			nsel++;
		} else mp->select(0);
	}

	if ( verbose & VERB_PROCESS )
		cout << "Models with " << ncomp_min << " to " << ncomp_max << " components selected: " << nsel << endl << endl;
	
	return nsel;
}

/**
@brief 	Selects a random number of components.
@param 	*model		model parameters.
@param 	number		number of components.
@return int			0.
**/
int			model_select_random(Bmodel* model, long number)
{
	model->deselect_all();
	
	random_seed();

	vector<Bcomponent*>	carr = model->component_array();
	long		ncomp(carr.size()), nsel(0);
	double		irm = ncomp*1.0L/get_rand_max();

	if ( verbose & VERB_PROCESS )
		cout << "Randomly selecting " << number << " components" << endl;
	
	long		j;
	while ( nsel < number ) {
		j = irm*random();
		if ( j < ncomp && carr[j]->select() == 0 ) {
			carr[j]->select(1);
			nsel++;
		}
	}
		
	return 0;
}

/**
@brief 	Determines if a polyhedron is closed given a rule.
@param 	*model			model parameters.
@param 	closure_rule	1=valency, 2=order.
@param 	val_order		magnitude of valency or order.
@return long			number of models selected.

	Polyhedron closure is arbitrarily decided by a rule:
	1.	Fixed valency.
	3.	Fixed polygon order.
	Models that fail the given rule are deselected.
	If the given valency or order is zero, it is set to the value for
	the first component or polygon.

**/
long		model_select_closed(Bmodel* model, int closure_rule, int val_order)
{
	long			n(0), i;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bpolygon*		poly;
	
	if ( closure_rule == 2 && !model->poly ) model_poly_generate(model);
	
	if ( verbose ) {
		cout << "Selecting models based on fixed ";
		if ( closure_rule == 1 ) cout << "valency of " << val_order << endl << endl;
		else cout << "polygon order of " << val_order << endl << endl;
	}
	
	for ( mp = model; mp; mp = mp->next ) {
		if ( !mp->comp ) mp->select(0);
		if ( closure_rule == 1 ) {
			for ( comp = mp->comp; comp && mp->select(); comp = comp->next ) {
				for ( i=0; i<comp->link.size() && comp->link[i]; i++ ) ;
				if ( val_order ) {
					if ( i != val_order ) mp->select(0);
				} else {
					val_order = i;
				}
			}
		} else if ( closure_rule == 2 ) {
			for ( poly = mp->poly; poly && mp->select(); poly = poly->next ) {
				for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) ;
				if ( val_order ) {
					if ( i != val_order ) mp->select(0);
				} else {
					val_order = i;
				}
			}
		}
		if ( mp->select() ) n++;
	}

	return n;
}

/**
@brief 	Selects fullerene type models.
@param 	*model		model structure.
@return long		number of models selected.

	Each model is tested for the presence of polygons that have orders
	only of five or six.

**/
long		model_select_fullerene(Bmodel* model)
{
	long			n(0);
	Bmodel*			mp;
	Bpolygon*		poly;
	
	if ( !model->poly ) model_poly_generate(model);
	
	if ( verbose )
		cout << "Selecting fullerene models" << endl << endl;
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( poly = mp->poly; poly && mp->select(); poly = poly->next )
			if ( poly->size() < 5 || poly->size() > 6 ) mp->select(0);
		if ( mp->select() ) n++;
	}
	
	return n;
}

/**
@brief 	Selects non-fullerene type models.
@param 	*model		model structure.
@return long		number of models selected.

	Each model is tested for the presence of polygons that have orders
	different from five or six.

**/
long		model_select_non_fullerene(Bmodel* model)
{
	long			n(0);
	Bmodel*			mp;
	Bpolygon*		poly;
	
	if ( !model->poly ) model_poly_generate(model);
	
	if ( verbose )
		cout << "Selecting non-fullerene models" << endl << endl;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		mp->select(0);
		for ( poly = mp->poly; poly && mp->select() == 0; poly = poly->next )
			if ( poly->size() < 5 || poly->size() > 6 ) mp->select(1);
		if ( mp->select() ) n++;
	}
	
	return n;
}

/**
@brief 	Selects fullerene type models.
@param 	*model		model structure.
@param 	valence		component valence
@return long		number of models selected.

	Each model is tested for the presence of polygons that have orders
	only of five or six.

**/
long		model_select_valence(Bmodel* model, int valence)
{
	long			n(0), v;
	Bmodel*			mp;
	Bcomponent*		comp;
	
	if ( verbose )
		cout << "Selecting components with valency " << valence << endl << endl;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			comp->select(0);
			for ( v=0; v<comp->link.size() && comp->link[v]; v++ ) ;
			if ( v == valence ) comp->select(1);
		}
		n++;
	}
	
	return n;
}

/**
@brief 	Selects model polygons with a given order.
@param 	*model		model to color.
@param 	order		polygon order.
@return long		0.
**/
long		model_select_polygons(Bmodel* model, int order)
{
	if ( order < 1 ) return 0;
	
	if ( !model->poly ) model_poly_generate(model);
	
	long			i;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bpolygon*		poly;
	
	if ( verbose )
		cout << "Selecting polygons with order " << order << endl << endl;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) comp->select(0);
		for ( poly = mp->poly; poly; poly = poly->next ) if ( poly->size() == order ) {
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ )
				poly->comp[i]->select(1);
		}
	}
	
	return 0;
}

/**
@brief 	Selects the first number of components.
@param 	*model		model parameters.
@param 	first		first number of components to select.
@return long			number of components selected.
**/
long 		model_select_first(Bmodel* model, int first)
{
	long			n, nsel(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	
	for ( n=0, mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next, n++ ) {
			if ( comp->select() > 0 && n < first ) {
				comp->select(1);
				nsel++;
			} else comp->select(0);
		}
	}

	return nsel;
}

/**
@brief 	Selects components within a shell.
@param 	*model		model parameters.
@param 	center		center of the shell.
@param 	minrad		minimum radius of the sphere.
@param 	maxrad		radius of the sphere.
@return long			number of components selected.
**/
long 		model_select_within_shell(Bmodel* model, Vector3<double> center, double minrad, double maxrad)
{
	if ( maxrad < 1e-3 ) maxrad = 1e30;
	
	long			nsel(0);
	double			d;
	Bmodel*			mp;
	Bcomponent*		comp;
	
	for ( mp = model; mp; mp = mp->next ) if ( model->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			d = comp->location().distance(center);
			if ( d < minrad || d > maxrad )
				comp->select(-1);
			else
				nsel++;
		}
	}

	if ( verbose )
		cout << "Number of components selected:  " << nsel << endl << endl;

	return nsel;
}

/**
@brief 	Selects components within a mask.
@param 	*model		model parameters.
@param 	pmask		first number of components to select.
@return long			number of components selected.
**/
long 		model_select_in_mask(Bmodel* model, Bimage* pmask)
{
	long			nsel(0), ntot(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	Vector3<double>	loc;
	
	if ( verbose )
		cout << "Deselecting components outside mask " << pmask->file_name() << endl;

	for ( ntot=0, mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next, ntot++ ) {
			loc = comp->location()/pmask->sampling(0) + pmask->image->origin();
//			cout << comp->identifier() << tab << loc << tab << pmask->get(0, loc) << endl;
			if ( pmask->get(0, loc) < 0.5 ) comp->select(0);
			else nsel += (comp->select() > 0);
		}
	}

	if ( verbose )
		cout << "Number of components selected:  " << nsel << " (" << nsel*100.0/ntot << " %)" << endl << endl;

	return nsel;
}

/**
@brief 	Deletes tagged models, components and links.
@param 	**model		pointer to model parameters.
@return long			remaining number of components.

	A model, component or link is tagged for deletion by a negative selection flag.

**/
long 		model_delete(Bmodel** model)
{
	long			nmod(0), ncomp(0), nlink(0), nmdel(0), ncdel(0), nldel(0);
	Bmodel*			mp_prev;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomponent*		comp_prev;
	Blink*			link;
	Blink*			link_prev;
	
	if ( verbose )
		cout << "Deleting models and model elements" << endl << endl;

	for ( mp = mp_prev = *model; mp; ) {
		if ( mp->select() < 0 ) {
//			if ( verbose & VERB_PROCESS )
//				cout << "Model: " << mp->identifier() << endl;
			if ( mp == *model ) {
				mp_prev = *model = mp->next;
				mp->next = NULL;
				model_kill(mp);
				mp = mp_prev;
			} else {
				mp_prev->next = mp->next;
				mp->next = NULL;
				model_kill(mp);
				mp = mp_prev->next;
			}
			nmdel++;
		} else {
			mp_prev = mp;
			mp = mp->next;
			nmod++;
		}
	}
	
	for ( mp = *model; mp; mp = mp->next ) {
//		if ( verbose & VERB_PROCESS )
//			cout << "Model: " << mp->identifier() << endl;
		for ( link = link_prev = mp->link; link;  ) {
			if ( link->select() < 0 || link->comp[0]->select() < 0 || link->comp[1]->select() < 0 ) {
				vector<Bcomponent*>		nulinks;
				for ( auto comp: link->comp[0]->link )
					if ( comp != link->comp[1] ) nulinks.push_back(comp);
				link->comp[0]->link = nulinks;
				nulinks.clear();
				for ( auto comp: link->comp[1]->link )
					if ( comp != link->comp[0] ) nulinks.push_back(comp);
				link->comp[1]->link = nulinks;
				if ( link == mp->link ) {
					link_prev = mp->link = link->next;
					delete link;
					link = link_prev->next;
				} else {
					link_prev->next = link->next;
					delete link;
					link = link_prev->next;
				}
				nldel++;
			} else {
				link_prev = link;
				link = link->next;
				nlink++;
			}
		}
		for ( comp = comp_prev = mp->comp; comp; ) {
			if ( comp->select() < 0 ) {
				if ( comp == mp->comp ) {
					comp_prev = mp->comp = comp->next;
					delete comp;
					comp = comp_prev;
				} else {
					comp_prev->next = comp->next;
					delete comp;
					comp = comp_prev->next;
				}
				ncdel++;
			} else {
				comp_prev = comp;
				comp = comp->next;
				ncomp++;
			}
		}
	}
	
	nmod += nmdel;
	ncomp += ncdel;
	nlink += nldel;

	if ( verbose & VERB_PROCESS ) {
		if ( nmod )  cout << "Number of models deleted:       " << nmdel << " (" << nmdel*100.0/nmod << " %)" << endl;
		if ( ncomp ) cout << "Number of components deleted:   " << ncdel << " (" << ncdel*100.0/ncomp << " %)" << endl;
		if ( nlink ) cout << "Number of links deleted:        " << nldel << " (" << nldel*100.0/nlink << " %)" << endl;
		cout << endl;
	}
	
	return ncomp - ncdel;
}

/**
@brief 	Deletes components based on their types.
@param 	*model		model parameters.
@param 	&comptype	component type.
@return long			number of components remaining.
**/
long 		model_delete_comp_type(Bmodel* model, Bstring& comptype)
{
	long			n(0);
	Bmodel*			mp;
	Bcomponent*		comp, *comp_prev;
	Bcomptype*		ct, *ctp;
	Blink*			link, *link_prev;
	
	for ( mp = model; mp; mp = mp->next ) {
		ct = mp->find_type(comptype.str());
		for ( link = link_prev = mp->link; link;  ) {
			if ( link->comp[0]->type() == ct || link->comp[1]->type() == ct ) {
				if ( link == mp->link ) {
					link_prev = mp->link = link->next;
					delete link;
					link = link_prev;
				} else {
					link_prev = link->next;
					delete link;
					link = link_prev->next;
				}
			} else {
				link_prev = link;
				link = link->next;
			}
		}
		for ( comp = comp_prev = mp->comp; comp; ) {
			if ( comp->type() == ct ) {
				if ( comp == mp->comp ) {
					comp_prev = mp->comp = comp->next;
					delete comp;
					comp = comp_prev;
				} else {
					comp_prev->next = comp->next;
					delete comp;
					comp = comp_prev->next;
				}
			} else {
				comp_prev = comp;
				comp = comp->next;
				n++;
			}
		}
		for ( ct = ctp = mp->type; ct;  ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG model_delete_comp_type: type=" << ct->identifier() << endl;
			if ( ct->identifier() == comptype.str() ) {
				if ( ct == mp->type ) {
					ctp = mp->type = ct->next;
					delete ct;
					ct = ctp;
				} else {
					ctp->next = ct->next;
					delete ct;
					ct = ctp->next;
				}
			} else {
				ctp = ct;
				ct = ct->next;
			}
		}
	}
	
	return n;
}

/**
@brief 	Deletes selected components and associated links from the model.
@param 	**model		pointer to model parameters.
@return long		remaining number of components.
**/
long 		model_delete_non_selected(Bmodel** model)
{
	long			nmod(0), ncomp(0), nlink(0), nmdel(0), ncdel(0), nldel(0);
	Bmodel*			mp_prev;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomponent*		comp_prev;
	Blink*			link;
	Blink*			link_prev;
	
	if ( verbose & VERB_PROCESS )
		cout << "Deleting non-selected components" << endl;
	
	for ( mp = mp_prev = *model; mp; ) {
		if ( verbose & VERB_PROCESS )
			cout << "Model: " << mp->identifier() << endl;
		if ( mp->select() < 1 ) {
			if ( mp == *model ) {
				mp_prev = *model = mp->next;
				mp->next = NULL;
				model_kill(mp);
				mp = mp_prev;
			} else {
				mp_prev->next = mp->next;
				mp->next = NULL;
				model_kill(mp);
				mp = mp_prev->next;
			}
			nmdel++;
		} else {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG model_delete_non_selected: deleting links" << endl;
			for ( link = link_prev = mp->link; link;  ) {
				if ( link->select() < 1 ||
						( link->comp[0] && link->comp[0]->select() < 1 ) ||
						( link->comp[1] && link->comp[1]->select() < 1 ) ) {
					if ( link == mp->link ) {
						link_prev = mp->link = link->next;
						delete link;
						link = link_prev;
					} else {
						link_prev = link->next;
						delete link;
						link = link_prev->next;
					}
					nldel++;
				} else {
					link_prev = link;
					link = link->next;
					nlink++;
				}
			}
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG model_delete_non_selected: deleting components" << endl;
			for ( comp = comp_prev = mp->comp; comp; ) {
				if ( comp->select() < 1 ) {
					if ( comp == mp->comp ) {
						comp_prev = mp->comp = comp->next;
						delete comp;
						comp = comp_prev;
					} else {
						comp_prev->next = comp->next;
						delete comp;
						comp = comp_prev->next;
					}
					ncdel++;
				} else {
					comp_prev = comp;
					comp = comp->next;
					ncomp++;
				}
			}
			mp_prev = mp;
			mp = mp->next;
			nmod++;
		}
	}
	
//	model_count_component_types(*model);
	(*model)->update_type_counts();
	
	nmod += nmdel;
	ncomp += ncdel;
	nlink += nldel;

	if ( verbose & VERB_PROCESS ) {
		if ( nmod )  cout << "Number of models deleted:       " << nmdel << " (" << nmdel*100.0/nmod << " %)" << endl;
		if ( ncomp ) cout << "Number of components deleted:   " << ncdel << " (" << ncdel*100.0/ncomp << " %)" << endl;
		if ( nlink ) cout << "Number of links deleted:        " << nldel << " (" << nldel*100.0/nlink << " %)" << endl;
		cout << endl;
	}
	
	return ncomp - ncdel;
}

/**
@brief 	Converts selection numbers to different types.
@param 	*model			model parameters.
@param 	*comp_type		linked list of strings.
@param 	&filename		component type file name.
@return long			number of component types assigned.
**/
long		model_type_from_selection(Bmodel* model, Bstring* comp_type, Bstring& filename)
{
	long			i, nname, nselnum, ntype;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Bcomptype*		ct = NULL;
	Bstring*		ct_name;
	Bstring			base = model->type->identifier();
	Bstring			id;
	
	for ( nname = 0, ct_name = comp_type; ct_name; ct_name = ct_name->next )
		if ( ct_name->length() > 0 ) nname++;

	for ( ntype = 0, mp = model; mp; mp = mp->next ) {
		for ( i=0, ct = mp->type; ct; ct = ct->next ) i++;
		if ( ntype < i ) ntype = i;
	}
				
	for ( nselnum = 0, mp = model; mp; mp = mp->next )
		for ( comp = mp->comp; comp; comp = comp->next )
			if ( nselnum < comp->select() ) nselnum = comp->select();
	
	if ( verbose )
		cout << "Creating %d component types" << nselnum << endl << endl;
	
	if ( ntype != nselnum ) {
		if ( nname == nselnum ) {
			for ( mp = model; mp; mp = mp->next ) {
				comp_type_list_kill(mp->type);
				mp->type = NULL;
				if ( verbose & VERB_FULL )
					cout << "Model: " << mp->identifier() << endl;
				for ( i=0, ct_name = comp_type; i<nselnum && ct_name; i++, ct_name = ct_name->next ) {
					if ( filename.length() )
						ct = mp->add_type(ct_name->str(), filename.str(), 0);
					else
						ct = mp->add_type(ct_name->str());
					if ( verbose & VERB_FULL )
						cout << "\tType: " << ct->identifier() << endl;
				}
			}
		} else  {
			if ( comp_type && comp_type->length() ) base = *comp_type;
			for ( mp = model; mp; mp = mp->next ) {
				comp_type_list_kill(mp->type);
				mp->type = NULL;
				if ( verbose & VERB_FULL )
					cout << "Model: " << mp->identifier() << endl;
				for ( i=0; i<nselnum; i++ ) {
					id = base + Bstring(i+1, "%03d");
					if ( filename.length() )
						ct = mp->add_type(id.str(), filename.str(), 0);
					else
						ct = mp->add_type(id.str());
					if ( verbose & VERB_FULL )
						cout << "\tType: " << ct->identifier() << endl;
				}
			}
		}
	} else if ( nname == ntype ) {
		for ( mp = model; mp; mp = mp->next )
			for ( ct = mp->type, ct_name = comp_type; ct && ct_name; ct = ct->next, ct_name = ct_name->next )
				ct->identifier() = ct_name->str();
	}
	
	for ( mp = model; mp; mp = mp->next )
		for ( i=0, ct = mp->type; ct; ct = ct->next, i++ ) {
			ct->file_name(filename.str());
			ct->image_number(i);
		}
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) {
			for ( i=1, ct = mp->type; i<nselnum && i<comp->select() && ct; i++, ct = ct->next ) ;
			if ( ct ) comp->type(ct);
		}
	}
	
	return nselnum;
}

/**
@brief 	Deselects components based on the FOM.
@param 	*model			model parameters.
@param 	fom_cutoff		FOM threshold.
@return long			number of components selected.

	The FOM is assumed to be a value from 0 to 1 and components with
	FOM's below the given cutoff are deselected.

**/
long		model_fom_deselect(Bmodel* model, double fom_cutoff)
{
	long			nsel(0);
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;

	for ( nsel=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) {
			if ( comp->FOM() < fom_cutoff ) comp->select(0);
			if ( comp->select() > 0 ) nsel++;
		}
	}
	
	return nsel;
}

/**
@brief 	Deselects components based on the FOM relative to the maximum FOM.
@param 	*model			model parameters.
@param 	fom_fraction	FOM threshold.
@return long			number of components selected.

	The FOM/FOMmax ratio is calculated and components with ratios below
	the given fraction are deselected.

**/
long		model_fom_max_fraction_deselect(Bmodel* model, double fom_fraction)
{
	if ( fom_fraction > 1 ) fom_fraction = 0.5;
	
	long			nsel(0);
	double			fom_max = 1e-6, fom_cut;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;

	for ( nsel=0, mp = model; mp; mp = mp->next ) if ( mp->select() )
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 )
			if ( fom_max < comp->FOM() ) fom_max = comp->FOM();

	fom_cut = fom_fraction*fom_max;
	
	for ( nsel=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			if ( comp->FOM() < fom_cut ) comp->select(0);
			if ( comp->select() > 0 ) nsel++;
		}
	}
	
	return nsel;
}

/**
@brief 	Deselects components based on the distance from the origin.
@param 	*model		model parameters.
@param 	minrad		minimum distance from origin.
@param 	maxrad		maximum distance from the origin.
@return long		number of components selected.

	All components inside a minmum radius and beyond a maximum radius are deselected.

**/
long		models_radius_deselect(Bmodel* model, double minrad, double maxrad)
{
	if ( maxrad < 1e-3 ) maxrad = 1e30;
	
	long			nsel(0), ntot(0);
	double			len;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	
	if ( verbose )
		cout << "Deselecting components outside " << minrad << " - " << maxrad << " A" << endl;

	for ( ntot=0, mp = model; mp; mp = mp->next )
		for ( comp = mp->comp; comp; comp = comp->next ) ntot++;

	for ( nsel=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			len = comp->location().length();
			if ( len < minrad || len > maxrad ) comp->select(-1);
			if ( comp->select() > 0 ) nsel++;
		}
	}

	if ( verbose )
		cout << "Number of components selected:  " << nsel << " (" << nsel*100.0/ntot << " %)" << endl << endl;

	return nsel;
}

/**
@brief 	Selects components based on the FOM.
@param 	*model		model parameters.
@param 	fom_step	FOM step size.
@return long		number of components.

	The FOM is assumed to represent energy and therefore values below
	the given cutoff are selected.

**/
long		model_fom_histogram(Bmodel* model, double fom_step)
{
	if ( fom_step < 0.001 ) fom_step = 0.001;
	
	long			n;
	double			fmin = 1e30, fmax = -1e30, pc;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;

	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			if ( comp->FOM() < fmin ) fmin = comp->FOM();
			if ( comp->FOM() > fmax ) fmax = comp->FOM();
			n++;
		}
	}
	
	int				i, sum, nhc(0), nh = (int) ((fmax - fmin)/fom_step + 1);
	int*			hist = new int[nh];
	for ( i=0; i<nh; i++ ) hist[i] = 0;

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			i = (int) ((comp->FOM() - fmin)/fom_step);
			if ( i >= 0 && i < nh ) {
				hist[i]++;
				nhc++;
			}
		}
	}
	
	cout << "FOM histogram:\nBin\tValue\tCount\tSum\t%" << endl;
	for ( i=0, sum=0; i<nh; i++ ) {
		sum += hist[i];
		pc = sum*100.0/n;
		cout << i << tab << (i+0.5)*fom_step + fmin << tab << hist[i] << tab << sum << tab << pc << endl;
	}
	cout << endl;
	
	delete[] hist;
	
	return n;
}

/**
@brief 	Ranks components based on the FOM.
@param 	*model			model parameters.
@param 	nrank			number of levels.
@return long			number of selected components.

	Each selected component is ranked into one of a number of levels.

**/
long		model_fom_ranking(Bmodel* model, int nrank)
{
	if ( nrank < 2 ) nrank = 2;
	
	long			i, n(0), nsel;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Bcomptype*		ct = NULL, *ctp;
	Bstring			ct_id;
	
	if ( verbose )
		cout << "Ranking selection into %d levels" << nrank << endl << endl;

//	model_count_component_types(model);
	model->update_type_counts();

	for ( nsel=0, mp = model; mp; mp = mp->next ) if ( mp->select() )
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) nsel++;
	
	int_float*		rank = new int_float[nsel];

	for ( nsel=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			rank[nsel].i = nsel;
			rank[nsel].f = comp->FOM();
			nsel++;
		}
	}

	qsort((void *) rank, nsel, sizeof(int_float), 
			(int (*)(const void *, const void *)) QsortLargeToSmallIntFloat);
	
	int*				sel = new int[nsel];
	for ( i=0; i<nsel; i++ ) sel[i] = 0;

	for ( i=0; i<nsel; i++ ) {
		n = nrank - (i*nrank)/nsel;
		sel[rank[i].i] = n;
	}
	
	for ( nsel=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG model_fom_ranking: model=" << mp->identifier() << " map=" << mp->mapfile() << endl;
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			comp->select(sel[nsel]);
			ct = comp->type();
			ct->component_count_decrement();
			ct_id = ct->identifier().c_str() + Bstring(comp->select(), "%03d");
//			comp->type = model_add_type_by_id_and_filename(mp, ct_id, ct->file_name() , ct->image_number());
			comp->type(mp->add_type(ct_id.str(), ct->file_name() , ct->image_number()));
			comp->type()->select(comp->select());
			comp->type()->component_count_increment();
			comp->type()->mass(ct->mass());
			nsel++;
		}
		for ( ct = ctp = mp->type; ct;  ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG model_fom_ranking: type=" << ct->identifier() << " num=" << ct->component_count() << endl;
			if ( ct->component_count() < 1 ) {
				if ( ct == mp->type ) {
					ctp = mp->type = ct->next;
					delete ct;
					ct = ctp;
				} else {
					ctp->next = ct->next;
					delete ct;
					ct = ctp->next;
				}
			} else {
				ctp = ct;
				ct = ct->next;
			}
		}
	}

	delete[] sel;
	delete[] rank;

	return n;
}

/**
@brief 	Deletes components that are too close to others.
@param 	**model			model handle.
@param 	distance		distance to test for.
@return long			0.

	If a component is too close to another, the component further down
	in the list is deleted.

**/
long		model_delete_overlapped_components(Bmodel** model, double distance)
{
	Bmodel*			mp;
	Bcomponent		*comp1, *comp2;
	
	for ( mp = *model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp1 = mp->comp; comp1->next; comp1 = comp1->next ) if ( comp1->select() ) {
			for ( comp2 = comp1->next; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
				if ( comp1->location().distance(comp2->location()) < distance )
					comp2->select(0);
			}
		}
	}
	
	model_delete_non_selected(model);
	
	return 0;
}

/**
@brief 	Averages components within a distance.
@param 	*model			model.
@param 	distance		overlap distance.
@return long			0.

	Components within the given distance are averaged and the first
	instance is retained while the others are deleted.

**/
long		model_average_overlapped_components(Bmodel* model, double distance)
{
	Matrix			dm = model_distance_matrix(model, 0);

	if ( verbose ) {
		cout << "Consolidating components:" << endl;
		cout << "Distance criterion:             " << distance << endl;
		cout << "Number of components:           " << dm.rows() << endl << endl;
	}
	
	long			n, i, j;
	Bcomponent*		comp;
	Bcomponent*		comp2;

	for ( i=0, comp = model->comp; comp; comp = comp->next, i++ ) if ( comp->select() ) {
		n = 1;
		for ( j=0, comp2 = model->comp; comp2; comp2 = comp2->next, j++ ) if ( i != j ) {
			if ( dm[i][j] < distance ) {
				comp->shift(comp2->location());
				comp->FOM(comp->FOM() + comp2->FOM());
				n++;
				comp2->select(0);
			}
		}
		if ( n > 1 ) {
			comp->scale(1/n);
			comp->FOM(comp->FOM() / n);
		}
	}
	
	return model_delete_non_selected(&model);
}


/**
@brief	Prunes overlapping components based on a distance criterion.
@param 	*model			model structure to be modified.
@param 	mindist			distance criterion.
@return long				number of remaining components.

	The first component in any pair of overlapping components is kept.

**/
long		model_prune_simple(Bmodel* model, double mindist)
{
	long			n, ns, ndel(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomponent*		comp2;

	for ( n=ns=0, comp = model->comp; comp; comp = comp->next, n++ ) if ( comp->select() > 0 ) ns++;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Pruning a set of solutions based on distance:" << endl;
		cout << "Number of components:           " << n << endl;
		cout << "Number selected:                " << ns << endl;
		cout << "Overlap distance:               " << mindist << endl << endl;
	}

	for ( ndel=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp->next; comp = comp->next ) if ( comp->select() > 0 ) {
			for ( comp2 = comp->next; comp2; comp2 = comp2->next ) if ( comp2->select() > 0 ) {
				if ( comp->location().distance(comp2->location()) < mindist ) {
					comp2->select(-1);
					ndel++;
				}
			}
		}
	}
	
	if ( verbose & VERB_FULL )
		cout << "Number of components deleted:   " << ndel << " (" << ndel*100.0/n << " %)" << endl;
	
	model_delete(&model);

	return n - ndel;
}

/**
@brief 	Deselects components based on overlap distance and the FOM.
@param 	*model			model parameters.
@param 	distance		overlap distance.
@return long				number of components selected.

	The selected component with the highest FOM is chosen as the first 
	component in the solution set and all components within the distance 
	criterion from it are deselected. The selected component with the next 
	highest FOM is then chosen and again, all components within the distance 
	criterion from it are deselected.
	This is repeated until there are no more selected components to choose.
	Only the first model is processed.

**/
long		model_prune_fom_old(Bmodel* model, double distance)
{
	int				h, i, j, k, n, ns, nsel, imax;
	double			fom_max, fom;
	Bcomponent*		comp = NULL;
	Bcomponent*		comp2 = NULL;

	for ( n=ns=0, comp = model->comp; comp; comp = comp->next, n++ ) if ( comp->select() > 0 ) ns++;
	
	if ( verbose ) {
		cout << "Pruning a set of solutions based on distance and FOM:" << endl;
		cout << "Number of components:           " << n << endl;
		cout << "Number selected:                " << ns << endl;
		cout << "Overlap distance:               " << distance << endl << endl;
	}
	
	int*			sol = new int[ns];
	char*			a = new char[ns*ns];
	double*			f = new double[ns];
	
	for ( i=0; i<ns; i++ ) f[i] = sol[i] = 0;
	for ( i=0; i<ns*ns; i++ ) a[i] = 0;

	// Sets up the adjacency matrix
	for ( i=k=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		for ( comp2 = model->comp; comp2; comp2 = comp2->next ) if ( comp2->select() > 0 ) {
			if ( comp != comp2 )
				if ( comp->location().distance(comp2->location()) < distance ) a[k] = 1;
			k++;
		}
		f[i] = comp->FOM();
		i++;
	}
	
	// Successively find the highest FOM and set deletion tag for all others too close
	for ( h=0, imax = 0; h<ns && imax > -1; h++ ) {
		for ( i=0, imax = -1, fom_max = -1; i<ns; i++ ) {
			if ( fom_max < f[i] ) {
				fom_max = f[i];
				imax = i;
			}
		}
		if ( imax > -1 ) {
			sol[imax] = 1;		// Tag for keeping
			f[imax] = -1;
			i = imax*ns;
			for ( j=0; j<ns; j++ ) if ( a[i+j] ) {
				f[j] = -1;
				sol[j] = -1;	// Tag for deletion
			}
		}
	}
	
	for ( i=nsel=0, fom=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		comp->select(sol[i++]);
		if ( comp->select() > 0 ) {
			nsel++;
			fom += comp->FOM();
		}
	}
	
	if ( nsel ) fom /= nsel;

	if ( verbose ) {
		cout << "Components selected:" << endl;
		cout << "ID\tType\tFOM\tDist" << endl;
		for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 )
			cout << comp->identifier() << tab << comp->type()->identifier() << tab << comp->FOM() << tab << comp->location().length() << endl;
		cout << "Number selected:                " << nsel << " (" << nsel*100.0/n << " )" << endl;
		cout << "Average FOM:                    " << fom << endl << endl;
	}
	
	delete[] sol;
	delete[] a;
	delete[] f;

	model_delete(&model);

	return nsel;
}

long		model_prune_fom(Bmodel* model, double distance)
{
	int				n, ns;
	Bcomponent*		comp = NULL;

	for ( n=ns=0, comp = model->comp; comp; comp = comp->next, n++ ) if ( comp->select() > 0 ) ns++;
	
	if ( verbose ) {
		cout << "Pruning a set of solutions based on distance and FOM:" << endl;
		cout << "Number of components:           " << n << endl;
		cout << "Number selected:                " << ns << endl;
		cout << "Overlap distance:               " << distance << endl << endl;
	}
	
	double							d;
	multimap<double,Bcomponent*>	cmm;
	
	for ( comp = model->comp; comp; comp = comp->next )
		if ( comp->select() > 0 ) cmm.insert(pair<double,Bcomponent*>(comp->FOM(), comp));
	
	for ( auto it1 = cmm.rbegin(); it1 != cmm.rend(); ++it1 )
			if ( it1->second->select() > 0 ) {
//		cout << it1->second->FOM() << endl;
		for ( auto it2 = it1; it2 != cmm.rend(); ++it2 )
				if ( it2 != it1 && it2->second->select() > 0 ) {
			d = it1->second->location().distance(it2->second->location());
			if ( d < distance ) it2->second->select(-1);
		}
	}
	
	return model_delete(&model);
}

/**
@brief 	Averages similar components and reduces redundancy.
@param 	*model		model parameters.
@return lon			number of clusters.

	A view distance matrix is calculated for all the components.
	All components within a small distance (~10% of model radius) are averaged
	into the first component and the rest are deselected.
	Only the first model is processed.

**/
long		model_prune_similar(Bmodel* model)
{
	long			i, j, n, ns, c(0);
	double			a;
	Bcomponent*		comp = NULL;
	Bcomponent*		comp1 = NULL;

	Matrix			dmat = model_distance_matrix(model, 1);

	for ( n=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		comp->view()[3] = angle_set_negPI_to_PI(comp->view().angle());
		comp->select(-1);
		n++;
	}
	
	for ( i=0, comp = model->comp; comp->next; comp = comp->next ) if ( comp->select() > 0 ) {
		for ( j=i+1, comp1 = comp->next; comp1; comp1 = comp1->next ) if ( comp1->select() > 0 ) {
			if ( dmat[i][j] < 0.1 ) {
				if ( comp->select() < 1 ) comp->select(++c);
				comp1->select(comp->select());
			}
			j++;
		}
		i++;
	}

	if ( verbose ) {
		cout << "Number of clusters:             " << c << endl;
		cout << "Averaging similar components:" << endl << endl;
	}
	
	for ( ns=0, comp = model->comp; comp; comp = comp->next )
		if ( ns < comp->select() ) ns = comp->select();
		
	for ( i=1; i<=ns; i++ ) {
		for ( n=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() == i ) {
			if ( n ) {
				comp->select(-1);
				comp1->location(comp1->location() + comp->location());
				comp1->view()[0] += comp->view()[0];
				comp1->view()[1] += comp->view()[1];
				comp1->view()[2] += comp->view()[2];
				a = comp1->view().angle()/n;
				comp1->view()[3] += a + angle_set_negPI_to_PI(comp->view().angle() - a);
			} else {
				comp1 = comp;		// First component
			}
			n++;
		}
		if ( n ) {
			comp1->scale(1.0/n);
			comp1->view()[0] /= n;
			comp1->view()[1] /= n;
			comp1->view()[2] /= n;
			comp1->view()[3] /= n;
			comp1->view().normalize();
		}
	}

	model_delete(&model);

	return ns;
}


Bimage*		model_composite_map(Bmodel* model, Vector3<int> size, 
				Vector3<double> origin, Vector3<double> sampling)
{
	if ( sampling[0] <= 0) sampling = sampling.max(1);

	long	i;
	double			v, t;
	Vector3<double>	loc, start, end;
	Bmodel*			mp = model;
	Bcomponent*		comp;
	
	Bimage*			p = new Bimage(Float, TSimple, size[0], size[1], size[2], 1);
	if ( origin[0] == 0 && origin[1] == 0 && origin[2] == 0 )
		origin = p->size()/2;
	p->origin(origin);
	p->sampling(sampling);
	
	long	datasize = p->data_size();
	
	if ( verbose ) {
		cout << "Assembling a model map:" << endl;
		cout << "Size:                           " << p->size() << endl;
		cout << "Origin:                         " << p->image->origin() << endl;
		cout << "Sampling:                       " << p->sampling(0) << " A/voxel" << endl;
	}

//	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			loc[0] = (long) (comp->location()[0]/p->sampling(0)[0] + p->image->origin()[0] + 0.5);
			loc[1] = (long) (comp->location()[1]/p->sampling(0)[1] + p->image->origin()[1] + 0.5);
			loc[2] = (long) (comp->location()[2]/p->sampling(0)[2] + p->image->origin()[2] + 0.5);
			p->sphere(loc, comp->radius()/p->sampling(0)[0], 0, FILL_USER, 1);
		}
//	}

	// Calculate the mask volume
	for ( i=0, v=t=0; i<datasize; i++ ) {
		if ( (*p)[i] > 0.5 ) v += 1;
		t += (*p)[i];
	}
	
	if ( verbose ) {
		cout << "Assembled map volume:           " << v * sampling.volume() << " A3" << endl;
		cout << "Assembled map integral:         " << t * sampling.volume() << endl;
		cout << "Integral volume ratio:          " << t/v << endl << endl;
	}
	
	return p;
}

long		comp_select_one(Bcomponent* comp, long i, long n, double* ovm)
{
	long			j, k;
	Bcomponent*	comp2 = comp;
	
	cout << i << ":";
	
	for ( j=i+1, k=i*n+j, comp = comp->next; j<n && comp; j++, k++, comp = comp->next ) if ( comp->select() > 0 ) {
		if ( ovm[k] > 0 ) comp->select(0);
		else comp_select_one(comp, j, n, ovm);
	}
	
	for ( j=0; comp2; j++, comp2 = comp2->next ) cout << " " << j+1 << ":" << comp2->select();
	cout << endl;
	
	return 0;
}


double*		model_overlap_matrix(Bmodel* model, double distance)
{
	long			n, i, j;
	double			d;
	Bmodel*			mp = model;
	Bcomponent*		comp;
	Bcomponent*		comp2;	
	
//	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->select() )
//		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) n++;
	
	for ( n=0, comp = mp->comp; comp; comp = comp->next ) n++;
	
	if ( verbose ) {
		cout << "Calculating an overlap matrix:" << endl;
		cout << "Distance criterion:             " << distance << endl;
		cout << "Number of components:           " << n << endl << endl;
	}

	double*			ovm = new double[n*n];
	for ( i=0; i<n*n; i++ ) ovm[i] = 0;
	

//	for ( i=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( i=0, comp = mp->comp; comp->next; comp = comp->next, i++ ) {
			for ( j=i+1, comp2 = comp->next; comp2; comp2 = comp2->next, j++ ) {
				d = comp->location().distance(comp2->location())/distance;
				if ( d < 2 ) {
					ovm[i*n+j] = ovm[j*n+i] = 1 - (d/4)*(3-d*d/4);
				}
			}
		}
//	}
	
	for ( i=0; i<n; i++ ) {
		for ( j=0; j<n; j++ )
			cout << tab << ovm[i*n+j];
		cout << endl;
	}
	
	comp_select_one(model->comp, 0, n, ovm);

	return ovm;
}


/**
@brief 	Deselects components based on overlap distance and the FOM.
@param 	*model		model parameters.
@param 	distance	overlap distance.
@return long		number of components selected.

	The selected component with the highest FOM is chosen as the first 
	component in the solution set and all components within the distance 
	criterion from it are deselected. The selected component with the next 
	highest FOM is then chosen and again, all components within the distance 
	criterion from it are deselected.
	This is repeated until there are no more selected components to choose.
	Only the first model is processed.

**/
long		model_prune_fit(Bmodel* model, double distance)
{
	long			n, ns, nsel(0);
//	double			fom(0);
	Vector3<double>	min, max;
	Vector3<double>	sampling(5,5,5);
	Bcomponent*		comp = NULL;
//	Bcomponent*		comp1 = NULL;

	for ( n=ns=0, comp = model->comp; comp; comp = comp->next, n++ ) if ( comp->select() > 0 ) {
		comp->radius(distance);
		min = min.min(comp->location());
		max = max.max(comp->location());
		comp->view()[3] = angle_set_negPI_to_PI(comp->view().angle());
		ns++;
	}
	max = max.max(-min);
	max += 2*distance;
	max = max/sampling;
	
	Vector3<int>	size((int) (2*max[0]), (int) (2*max[1]), (int) (2*max[2]));
	Vector3<double>	origin = max;
	Bimage*			p = model_composite_map(model, size, origin, sampling);
	
//	fom = 1/p->FOM_maximum();
	
	write_img("ass.map", p, 0);
	
	delete p;
	

//	double*			ovm = model_overlap_matrix(model, distance);
//	delete[] ovm;
/*	
	for ( i=0; i<10 && fom < 0.9; i++ ) {
		// Change selection
		// Calculate new FOM
		p = model_composite_map(model, size, origin, sampling);
		fom = 1/p->FOM_maximum();
		delete p;
	}
*/
/*
	int				h, i, j, k, n, ns, nsel, imax;
	double			fom_max, fom;
	Bcomponent*		comp = NULL;
	Bcomponent*		comp2 = NULL;

	for ( n=ns=0, comp = model->comp; comp; comp = comp->next, n++ ) if ( comp->select() > 0 ) ns++;
	
	if ( verbose ) {
		cout << "Pruning a set of solutions based on distance and FOM:" << endl;
		cout << "Number of components:           %d\n", n);
		cout << "Number selected:                %d\n", ns);
		cout << "Overlap distance:               %g\n\n", distance);
	}
	
	int*			sol = new int[ns];
	char*			a = new char[ns*ns];
	double*			f = new double[ns];
	
	for ( i=0; i<ns; i++ ) f[i] = sol[i] = 0;
	for ( i=0; i<ns*ns; i++ ) a[i] = 0;

	// Sets up the adjacency matrix
	for ( i=k=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		for ( comp2 = model->comp; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			if ( comp->location().distance(comp2->location()) < distance ) a[k] = 1;
//			cout << " %d", a[k]);
			k++;
		}
//		cout << endl;
		f[i] = comp->FOM();
		i++;
	}
	
	for ( h=0; h<ns && imax > -1; h++ ) {
		for ( i=0, imax = -1, fom_max = -1; i<ns; i++ ) {
			if ( fom_max < f[i] ) {
				fom_max = f[i];
				imax = i;
			}
		}
		if ( imax > -1 ) {
			sol[imax] = 1;
			f[imax] = -1;
			i = imax*ns;
			for ( j=0; j<ns; j++ ) if ( a[i+j] ) f[j] = -1;
		}
	}
	
	for ( i=nsel=0, fom=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		comp->select() = sol[i++];
		if ( comp->select() > 0 ) {
			nsel++;
			fom += comp->FOM();
		}
	}
	
	if ( nsel ) fom /= nsel;

	if ( verbose ) {
		cout << "Components selected:" << endl;
		cout << "Number\tID\tType\tFOM" << endl;
		for ( i=0, comp = model->comp; comp; comp = comp->next, i++ ) if ( comp->select() > 0 )
			cout << "%d\t%s\t%s\t%g\n", i+1, comp->identifier().c_str(), comp->type.c_str(), comp->fom);
		cout << "Number selected:                %d (%g %%)\n", nsel, nsel*100.0/n);
		cout << "Average FOM:                    %g\n\n", fom);
	}
	
	delete[] sol;
	delete[] a;
	delete[] f;
*/
	return nsel;
}

/**
@brief 	Deselects components based on overlap distance and the FOM.
@param 	*model		model parameters.
@param 	sampling	grid sampling.
@return long		number of components selected.

	Each component location is mapped to a grid point and the component FOM
	is assigned to the grid point if it is larger.
	Peaks are determined within each 3x3x3 kernel: a peak is defined as
	having a higher value than any of its 26 neigbors.
	The component associated with the value at a peak is selected.
	Only the first model is processed.

**/
long		model_prune_large(Bmodel* model, double sampling)
{
	if ( !model->comp ) return -1;
	
	long	i, ii, x, y, z, xx, yy, zz, imax, n, ng, np, nsel, vol;
	double			dmax;
	Vector3<double>	min = model->comp->location();
	Vector3<double>	max = model->comp->location();
	Vector3<int>	size;
	Bcomponent*		comp = NULL;
	
	for ( n=nsel=0, comp = model->comp; comp; comp = comp->next, n++ ) if ( comp->select() > 0 ) nsel++;
	
	if ( verbose ) {
		cout << "Pruning a set of solutions based on grid sampling and FOM:" << endl;
		cout << "Number of components:           " << n << endl;
		cout << "Number selected:                " << nsel << endl;
		cout << "Grid sampling:                  " << sampling << endl;
	}
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		min = min.min(comp->location());
		max = max.max(comp->location());
	}
	min -= 2*sampling;
	max += 2*sampling;
	for ( i=0; i<3; i++ ) size[i] = (int) ((max[i] - min[i])/sampling);
	vol = (long) size.volume();
	
	double*			grid = new double[vol];
	double*			peak = new double[vol];
	for ( i=0; i<vol; i++ ) grid[i] = peak[i] = 0;
	
	if ( verbose )
		cout << "Grid volume:                    " << vol << endl << endl;
	
	// Set up grid with maxima
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		x = (long) ((comp->location()[0] - min[0])/sampling);
		y = (long) ((comp->location()[1] - min[1])/sampling);
		z = (long) ((comp->location()[2] - min[2])/sampling);
		i = (z*size[1] + y)*size[0] + x;
		if ( i >= vol ) cerr << "Volume exceeded: " << i << endl;
		if ( grid[i] < comp->FOM() ) grid[i] = comp->FOM();
	}
	
	if ( verbose & VERB_FULL )
		cout << "Find peaks in grid" << endl;

	// Find peaks within a 3x3x3 kernel
	for ( ng=np=0, z=1; z<size[2]-1; z++ ) {
		for ( y=1; y<size[1]-1; y++ ) {
			for ( x=1; x<size[0]-1; x++ ) {
				i = (z*size[1] + y)*size[0] + x;
				if ( grid[i] > 0 ) {
					imax = i;
					dmax = 0;
					for ( zz=z-1; zz<=z+1; zz++ ) {
						for ( yy=y-1; yy<=y+1; yy++ ) {
							for ( xx=x-1; xx<=x+1; xx++ ) {
								ii = (zz*size[1] + yy)*size[0] + xx;
								if ( dmax < grid[ii] ) {
									dmax = grid[ii];
									imax = ii;
								}
							}
						}
					}
					if ( imax == i ) {
						peak[i] = grid[i];
						np++;
					}
					ng++;
				}
			}
		}
	}
	
	if ( verbose & VERB_FULL )
		cout << "Find components associated with peaks" << endl;
	
	// Find models associated with peaks and select those
	for ( nsel=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		x = (long) ((comp->location()[0] - min[0])/sampling);
		y = (long) ((comp->location()[1] - min[1])/sampling);
		z = (long) ((comp->location()[2] - min[2])/sampling);
		i = (z*size[1] + y)*size[0] + x;
		if ( i >= vol ) cerr << "Volume exceeded: " << i << endl;
		comp->select(-1);
		if ( peak[i] > 0 ) if ( fabs(peak[i] - comp->FOM()) < 1e-10 ) {
			comp->select(1);
			nsel++;
		}
	}
	
	if ( verbose ) {
		cout << "Points placed on grid:          " << ng << endl;
		cout << "Peaks:                          " << np << endl;
		cout << "Components selected:           " << nsel << endl << endl;
	}
	
	delete[] grid;
	delete[] peak;

	model_delete(&model);

	return nsel;
}

/**
@brief 	Finds the components that overlap with the reference.
@param 	*model		model parameters.
@param 	&reffile	reference molecule file name.
@param 	distance	distance for overlap.
@return long		number of components selected.

	All components that are within the distance criterion of the components
	in the reference model are deselected.

**/
long		model_find_overlap(Bmodel* model, Bstring& reffile, double distance)
{
	Bmodel*			modref = read_model(reffile);

	int				nsel(0);
	Bmodel*			mp;
	Bmodel*			mpr;
	Bcomponent*		comp;
	Bcomponent*		compr;	
	
	if ( verbose ) {
		cout << "Deselecting components:" << endl;
		cout << "Reference model:                " << reffile << endl;
		cout << "Distance criterion:             " << distance << endl;
	}
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			for ( mpr = modref; mpr && comp->select(); mpr = mpr->next ) {
				for ( compr = mpr->comp; compr && comp->select(); compr = compr->next ) {
					if ( comp->location().distance(compr->location()) < distance ) comp->select(0);
				}
			}
			if ( comp->select() > 0 ) nsel++;
		}
	}

	model_kill(modref);
	
	if ( verbose )
		cout << "Number of components retained:  " << nsel << endl << endl;
	
	return nsel;
}

