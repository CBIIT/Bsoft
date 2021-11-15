/**
@file	model_poly.cpp
@brief	Functions to manipulate polyhedral coordinate files
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20150208
**/

#include "model_poly.h"
#include "model_transform.h"
#include "model_views.h"
#include "model_util.h"
#include "matrix_linear.h"
#include "model_links.h"
#include "model_compare.h"
#include "math_util.h"
#include "random_numbers.h"
#include "symmetry.h"
#include "Vector3.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Generates links between the vertices to define the polyhedron faces.
@param 	*model		model structure.
@return	int					0.

	Two vertices are linked only when they are on the surface of the 
	polyhedron, i.e., all other vertices are on one side of the pair 
	of vertices. This is only true if there are no other vertices v
	with the following property:
		v•v0 > v0•v0
	where v0 is the average of the two vertices considered for linkage.
	Only the first model is processed.
**/
int			model_poly_faces(Bmodel* model)
{
	int				found, nlink(0);
	double			d0;
	Vector3<double>	v0;
	Bcomponent*		comp = NULL;
	Bcomponent*		comp1 = NULL;
	Bcomponent*		comp2 = NULL;
	Blink*			link = NULL;
	
	for ( comp1 = model->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
		for ( comp2 = comp1->next; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			v0 = (comp1->location() + comp2->location())/2;
			d0 = v0.length2();
			for ( found = 0, comp = model->comp; comp && !found; comp = comp->next ) 
				if ( comp != comp1 && comp != comp2 )
					if ( v0.scalar(comp->location()) > d0 ) found = 1;
			if ( !found ) {
				link = link_add(&link, comp1, comp2, 0, 0);
				if ( !model->link ) model->link = link;
				nlink++;
			}
		}
	}

	if ( verbose )
		cout << "Number of links generated:      " << nlink << endl << endl;
	
	return 0;
}

/*
@brief 	Generates polygons based on a vertex network.
@param	*comp			current end vertex used to select next link.
@param 	ilink			previously selected link index.
@param 	nlink			number of links so far in polygon.
@param 	*poly			polygon structure.
@return int 			number of vertices in polygon.

	A choice is made between the two non-originating connections by calculating 
	the angle between the vertex normal and the cross product of the originating
	and non-originating bond. If this angle is smaller than 90 degrees, it is 
	selected to continue following. Each pair of bonds are flagged to not traverse
	the same polygon more than once.

**/
int			poly_get_connectivity(Bcomponent* comp, int ilink, int nlink, Bpolygon* poly)
{
	if ( poly->closed() ) return nlink;
	if ( !comp->link[ilink] ) return nlink;
	
	int				i, j;
	double			a, mina = TWOPI;
	Vector3<double>	e1, e2, cp;
	Bcomponent*		compsel = NULL;
	
	poly->comp.push_back(comp);
	nlink++;
	
	// Exit when the polygon is too large
	if ( nlink >= MAXLINK ) return nlink;
	
	// Find the closest link in one direction
	for ( i=0, j=-1; i<comp->link.size() && comp->link[i]; i++ ) {
		if ( i != ilink ) {
			e1 = comp->link[ilink]->location() - comp->location();
			e2 = comp->link[i]->location() - comp->location();
			cp = comp->force().cross(e1);
			a = e1.angle(e2);
			if ( e2.angle(cp) > M_PI_2 ) a = TWOPI - a;
			if ( mina > a ) {
				mina = a;
				j = i;
			}
		}
	}
	if ( j < 0 ) return 0;			// No new links found
	if ( comp->flag[j] ) return 0;	// Check if already encountered
	comp->flag[j] = 1;				// Set flag value
	
	compsel = comp->link[j];			// One more comp in polygon
	if ( compsel == poly->comp[0] ) {	// End of polygon
		poly->closed(1);
	} else {
		for ( i=0; i<compsel->link.size() && compsel->link[i] != comp; i++ ) ;
		if ( i<compsel->link.size() ) nlink = poly_get_connectivity(compsel, i, nlink, poly);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << tab << compsel->identifier();
	
	return nlink;
}

/**
@brief 	Generates polygons based on a vertex network.
@param 	*model	model structure.
@return int				0.

	The search startegy is to start at a vertex and then search always turning into
	the same direction. First the outward pointing normal for each vertex is
	calculated. Then the connectivity is followed always turning in the same 
	direction at each vertex.

**/
int			model_poly_generate(Bmodel* model)
{
	int				i, j, m, nmod, nmods, n, nc, npm(20);
	double			a, maxang;
	Vector3<double>	vec1, vec2, com, rloc;
	vector<int>		np(npm,0);
	vector<int>		npc(npm,0);
	vector<int>		nv(1000000,0);
	
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomponent*		comp2;
	Bpolygon*		poly = NULL;
	
//	model->show_links();

	if ( verbose & VERB_PROCESS )
		cout << "Generating a polygon list" << endl << endl;
		
	if ( verbose & VERB_FULL ) {
		cout << "Vertex normals:" << endl;
		cout << "ID\tLinks\tx\ty\tz" << endl;
	}
	
	// Calculate the vertex normals
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		com = model_center_of_mass(mp);
		for ( comp = mp->comp; comp; comp = comp->next ) {
			rloc = comp->location() - com;
			for ( i=0; i<comp->link.size() && comp->link[i]; i++ ) ;
//			cout << comp->identifier() << tab << i << tab << rloc << endl;
			comp->force(Vector3<float>(0,0,0));
			if ( i > 1 ) {
				maxang = TWOPI/(i - 1.0L);
				for ( i=1; i<comp->link.size() && comp->link[i]; i++ ) {
					vec1 = comp->link[i]->location() - comp->location();
					for ( j=0; j<i; j++ ) {
						vec2 = comp->link[j]->location() - comp->location();
						a = vec1.angle(vec2);
						if ( a > 0.001 && a < maxang ) {
							vec2 = comp->location().normal(comp->link[i]->location(), comp->link[j]->location());
							if ( vec2.angle(rloc) > M_PI_2 ) vec2 = -vec2;
							comp->force(comp->force() + vec2);
						}
					}
				}
			}
			if ( comp->force() == 0 ) comp->force(rloc);
			comp->force().normalize();
			if ( verbose & VERB_FULL )
				cout << comp->identifier() << tab << i << tab << comp->force() << endl;
		}
	}
	
	for ( mp = model; mp; mp = mp->next )
		for ( comp = mp->comp; comp; comp = comp->next )
			for ( i=0; i<comp->flag.size(); i++ ) comp->flag[i] = 0;

	// Generate the polygon list
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		poly_list_kill(mp->poly);
		mp->poly = NULL;
		poly = (Bpolygon *) add_item((char **) &mp->poly, sizeof(Bpolygon));
		for ( n=0, comp = mp->comp; comp; comp = comp->next ) {
			for ( i=n=0; i<comp->link.size() && comp->link[i]; i++ ) {
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG model_poly_generate: " << comp->identifier() << tab << i << endl;
				n = poly_get_connectivity(comp, i, 0, poly);
				if ( verbose & VERB_DEBUG )
					cout << "\tn=" << n << endl;
				if ( n > 2 && n < MAXLINK )
					poly = (Bpolygon *) add_item((char **) &poly, sizeof(Bpolygon));
				else
//					for ( j=0; j<poly->comp.size(); j++ ) poly->comp[j] = NULL;
					poly->comp.clear();
			}
		}	
		// Prune the last item
		if ( n < 3 ) {
			if ( mp->poly->next ) {
				for ( poly = mp->poly; poly->next->next; poly = poly->next ) ;
				delete poly->next;
				poly->next = NULL;
			} else {
				delete mp->poly;
				mp->poly = NULL;
			}
		}
	}

	for ( mp = model; mp; mp = mp->next )
		for ( comp = mp->comp; comp; comp = comp->next )
			for ( i=0; i<comp->flag.size(); i++ ) comp->flag[i] = 0;

	if ( verbose & VERB_FULL )
		cout << "Polygon normals:\nID\tVert\tx\ty\tz" << endl;
	for ( n=nc=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( poly = mp->poly; poly; poly = poly->next, n++ ) {
			poly->select(1);
			poly->normal(Vector3<float>(0,0,0));
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				if ( i+1 < poly->size() ) comp = poly->comp[i+1];
				else comp = poly->comp[0];
				if ( i+2 < poly->size() ) comp2 = poly->comp[i+2];
				else if ( i+1 < poly->size() ) comp2 = poly->comp[0];
				else comp2 = poly->comp[1];
				poly->normal(poly->normal() + comp->location().normal(poly->comp[i]->location(), comp2->location()));
				if ( poly->closed() ) {
					for ( j=0; j<poly->comp.size() && poly->comp[i]->link[j]; j++ )
						if ( poly->comp[i]->link[j] == comp )
							poly->comp[i]->flag[j] = poly->size();
				}
			}
			poly->normal().normalize();
			if ( poly->size() < npm ) {
				np[poly->size()]++;
				npc[poly->size()] += poly->closed();
				nc += poly->closed();
			}
			if ( verbose & VERB_FULL )
				cout << n+1 << tab << poly->size() << tab << poly->normal() << endl;
		}
	}

	for ( nmod=nmods=0, mp = model; mp; mp = mp->next, nmod++ ) if ( mp->select() ) {
		nmods++;
		for ( comp = mp->comp; comp; comp = comp->next ) {
			for ( i=m=1; i<comp->link.size() && comp->link[i]; i++ ) {
				m *= 10;
				for ( j=0; j<i; j++ )
					if ( comp->flag[j] > comp->flag[i] )
						swap(comp->flag[i], comp->flag[j]);
			}
			for ( i=j=0; i<comp->link.size() && comp->link[i]; i++ ) {
				j += m*comp->flag[i];
				m /= 10;
			}
			nv[j]++;
		}
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Models:                         " << nmod << " (" << nmods << ")" << endl;
		if ( np[3] ) cout << "Triangles:                      " << np[3] << " (" << npc[3] << ")" << endl;
		if ( np[4] ) cout << "Tetragons:                      " << np[4] << " (" << npc[4] << ")" << endl;
		if ( np[5] ) cout << "Pentagons:                      " << np[5] << " (" << npc[5] << ")" << endl;
		if ( np[6] ) cout << "Hexagons:                       " << np[6] << " (" << npc[6] << ")" << endl;
		if ( np[7] ) cout << "Heptagons:                      " << np[7] << " (" << npc[7] << ")" << endl;
		if ( np[8] ) cout << "Octagons:                       " << np[8] << " (" << npc[8] << ")" << endl;
		if ( np[8] ) cout << "Nonagons:                       " << np[9] << " (" << npc[9] << ")" << endl;
		if ( np[10] ) cout << "Decagons:                       " << np[10] << " (" << npc[10] << ")" << endl;
		cout << "Polygons:                       " << n << " (" << nc << ")" << endl << endl;
		cout << "Vertex types:\nIndex\tCount" << endl;
		for ( i=0; i<1000000; i++ ) if ( nv[i] > 0 ) cout << i << tab << nv[i] << endl;
		cout << endl;
	}
	
	return 0;
}

/**
@brief 	Determines the vertex type based on adjacent polygons.
@param 	*model		model structure.
@return int					0.

	The polygon order is written into flags for vertex links in one
	direction for each polygon.
	New component types are generated.

**/
int			model_vertex_types(Bmodel* model)
{
	int				i, j;
	Bstring			id;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bpolygon*		poly;
	
	if ( !model->poly ) model_poly_generate(model);
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( poly = mp->poly; poly; poly = poly->next ) {
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				if ( i+1 < poly->size() ) comp = poly->comp[i+1];
				else comp = poly->comp[0];
				for ( j=0; j<poly->comp[i]->link.size() && poly->comp[i]->link[j] != comp; j++ ) ;
				poly->comp[i]->flag[j] = poly->size();
			}
		}
		comp_type_list_kill(mp->type);
		mp->type = NULL;
		for ( comp = mp->comp; comp; comp = comp->next ) {
//			comp->type = 0;
//			for ( i=0; i<MAXLINK && comp->link[i]; i++ )
//				comp->type += Bstring(comp->flag[i], "%d");
//			comp->type = comp->type.canonical(1);
//			comp->type = model_add_type_by_id(mp, comp->type);
			id = 0;
			for ( i=0; i<comp->link.size() && comp->link[i]; i++ )
				id += Bstring(comp->flag[i], "%d");
			id = id.canonical(1);
//			comp->type = model_add_type_by_id(mp, id);
			comp->type(mp->add_type(id));
		}
	}
	
	return 0;
}

Bstring		component_6digit_type(Bcomponent* comp)
{
	int				i, j, t;
	double			a;
	Vector3<double>	v;
	Bstring			id;
	
	if ( !comp->link[2] ) return 0;
	
	v = comp->link[0]->location().normal(comp->link[1]->location(), comp->link[2]->location());
	a = comp->location().angle(v);
	
	if ( a < M_PI_2 ) {
		j = 0;
		t = 1;
	} else {
		j = 2;
		t = -1;
	}
	
//	comp->type = 0;
//	for ( i=0; comp->link[i]; i++, j+=t )
//		comp->type += Bstring(comp->flag[j]/100, "%d") +  Bstring(comp->flag[j]%10, "%d");
//	comp->type = comp->type.canonical(2);
	id = 0;
	for ( i=0; i<comp->link.size() && comp->link[i]; i++, j+=t )
		id += Bstring(comp->flag[j]/100, "%d") +  Bstring(comp->flag[j]%10, "%d");
	id = id.canonical(2);
	
	return id;
}


/**
@brief 	Determines the vertex type based on adjacent and opposed polygons.
@param 	*model		model structure.
@return int					0.

	The link flag of each component is asigned such that the order of the
	rigth adjacent polygon is in the first digit, that of the left adjacent
	polygon in the second digit, and that of the opposing polygon in the 
	third digit.
	New component types are generated.

**/
int			model_extended_vertex_types(Bmodel* model)
{
	int				i, j, k;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomponent*		comp2;
	Bpolygon*		poly;
	Bstring			id;
	
	if ( !model->poly ) model_poly_generate(model);

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next )
			for ( j=0; comp->link[j]; j++ ) comp->flag[j] = 0;
		for ( poly = mp->poly; poly; poly = poly->next ) {
			for ( i=0; poly->comp[i]; i++ ) {
				if ( i+1 < poly->size() ) comp = poly->comp[i+1];
				else comp = poly->comp[0];
				for ( j=0; poly->comp[i]->link[j] != comp; j++ ) ;
				poly->comp[i]->flag[j] += 100*poly->size();		// First digit = right polygon 
				for ( j=0; poly->comp[i] != comp->link[j]; j++ ) ;
				comp->flag[j] += 10*poly->size();				// Second digit = left polygon
			}
		}
		for ( poly = mp->poly; poly; poly = poly->next ) {
			for ( i=0; poly->comp[i]; i++ ) {
				if ( i+1 < poly->size() ) comp = poly->comp[i+1];
				else comp = poly->comp[0];
				if ( i+2 < poly->size() ) comp2 = poly->comp[i+2];
				else if ( i+1 < poly->size() ) comp2 = poly->comp[0];
				else comp2 = poly->comp[1];
				for ( j=0; poly->comp[i]->link[j] != comp; j++ ) ;
				for ( k=0; comp->link[k] == poly->comp[i] || comp->link[k] == comp2; k++ ) ;
				poly->comp[i]->flag[j] += comp->flag[k]/100;	// Third digit = opposed polygon
			}
		}
		comp_type_list_kill(mp->type);
		mp->type = NULL;
		for ( comp = mp->comp; comp; comp = comp->next ) {
			id = component_6digit_type(comp);
//			comp->type = model_add_type_by_id(mp, id);
			comp->type(mp->add_type(id));
			component_hand(id);
		}
	}
	
	return 0;
}

/**
@brief 	Calculates the dual of a polyhedral network.
@param 	*model		model structure.
@param 	order			order of polygons to convert to vertices, < 3 = all.
@return Bmodel*				new model structure with the dual.

	The polygons are first defined to calculate the dual network that
	has vertices at the polygon centers.

**/
Bmodel*		model_poly_dual(Bmodel* model, int order)
{
	long			i, n(0);
	double			side, rad, length(0);
	Bstring			id;
	Bmodel*			mp;
	Bpolygon*		poly = NULL;
	Bcomponent*		comp = NULL;
	Bmodel*			model_new = NULL;
	Bmodel*			mp_new = NULL;

	Bstring			ctstr[10];
	ctstr[0] = "VER";
	ctstr[1] = "MON";
	ctstr[2] = "DI";
	ctstr[3] = "TRI";
	ctstr[4] = "TET";
	ctstr[5] = "PEN";
	ctstr[6] = "HEX";
	ctstr[7] = "HEP";
	ctstr[8] = "OCT";
	ctstr[9] = "NON";
		
	if ( verbose )
		cout << "Generating a polyhedron dual" << endl << endl;

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
//		mp_new = (Bmodel *) add_item((char **) &mp_new, sizeof(Bmodel));
//		if ( !model_new ) model_new = mp_new;
		if ( mp_new ) mp_new = mp_new->add(mp->identifier());
		else model_new = mp_new = new Bmodel(mp->identifier());
//		mp_new->identifier() = mp->identifier();
		mp_new->model_type(mp->model_type());
		mp_new->select(1);
		comp = NULL;
		for ( n=0, poly = mp->poly; poly; poly = poly->next )
				if ( poly->closed() && ( order < 3 || poly->size() == order ) ) {
//			id = Bstring(++n, "%d");
//			comp = component_add(&comp, id);
//			if ( !mp_new->comp ) mp_new->comp = comp;
			if ( comp ) comp = comp->add(++n);
			else mp_new->comp = comp = new Bcomponent(++n);
			if ( poly->size() > 2 && poly->size() < 10 ) id = ctstr[poly->size()];
			else id = ctstr[0];
//			comp->type = model_add_type_by_id(mp_new, id);
			comp->type(mp_new->add_type(id));
/*			switch ( poly->size() ) {
				case 3: comp->type = "TRI"; break;
				case 4: comp->type = "TET"; break;
				case 5: comp->type = "PEN"; break;
				case 6: comp->type = "HEX"; break;
				case 7: comp->type = "HEP"; break;
				case 8: comp->type = "OCT"; break;
				case 9: comp->type = "NON"; break;
				default: comp->type = "VER";
			}*/
			comp->location(Vector3<float>(0,0,0));
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ )
				comp->shift(poly->comp[i]->location());
			comp->scale(1.0/i);
			comp->radius(poly->comp[0]->radius());
//			comp->type = model_add_type_by_id(mp_new, comp->type);
			side = poly->comp[0]->location().distance(poly->comp[poly->size()-1]->location());
			rad = poly->comp[0]->location().distance(comp->location());
			for ( i=1; poly->comp[i]; i++ ) {
				side += poly->comp[i]->location().distance(poly->comp[i-1]->location());
				rad += poly->comp[i]->location().distance(comp->location());
			}
			side /= i;
			rad /= i;
			length += sqrt(4*rad*rad - side*side);
		}
	}

	length /= n;
	
//	model_stats(model_new);
	
	model_link_list_generate(model_new, 1.2*length);

	model_poly_generate(model_new);
	
	return model_new;
}

/**
@brief 	Analyzes a model for polygon regularity and planarity.
@param 	*model	model structure.
@return int						0.
**/
int			model_poly_analyze(Bmodel* model)
{
	if ( !model->poly ) model_poly_generate(model);
	
	cout << "Analyzing a model:" << endl;

	Bmodel*			mp;
	Bpolygon*		poly;
	
	int				i, j, n;
	int				nn[10], np[10], nt[10];

	cout << "Polygon neighbours:\nID\tVert\tTet\tPent\tHex\tHept" << endl;
	for ( j=0; j<10; j++ ) np[j] = nt[j] = 0;
	for ( n=0, mp = model; mp; mp = mp->next ) {
		for ( poly = mp->poly; poly; poly = poly->next, n++ ) {
			for ( j=0; j<10; j++ ) nn[j] = 0;
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				for ( j=0; j<poly->comp[i]->link.size() && poly->comp[i]->link[j]; j++ )
					nn[poly->comp[i]->flag[j]]++;
			}
			nn[i] -= i;								// Delete self-references
			if ( nn[i] < 0 ) nn[i] = 0;
			for ( j=0; j<10; j++ ) nn[j] /= 2;		// Compensate for double counting
			if ( i == 4 ) nt[nn[5]]++;
			if ( i == 5 ) np[nn[5]]++;
			cout << n << tab << i << tab << nn[4] << tab << nn[5] << tab << nn[6] << tab << nn[7] << endl;
		}
	}

	model_poly_links(model);
	
	model_poly_angles(model);
	
	model_poly_regularity(model);
	
	model_poly_planarity(model);

//	Vector3<double>	pax = model_principal_axes(model, NULL);
	model_principal_axes(model);
/*
	double			Ap, As, Vs, g=GOLDEN;
	Vs = (M_PI*4.0/3.0)*pax.volume();
	Ap = 4.0*M_PI*pow((pow((double)pax[0]*pax[1],g)+pow((double)pax[0]*pax[2],g)+pow((double)pax[1]*pax[2],g))/3.0, 1/g);
	As = pow(M_PI, 1.0/3.0) * pow(6*Vs, 2.0/3.0);
	cout << "Ellipsoid volume:               " << Vs << " A3" << endl;
//	cout << "Polyhedral volume:              " << Vp << " A3" << endl;
	cout << "Sphericity:                     " << As/Ap << endl;
	cout << "Axes ratio:                     " << pax[0]/pax[2] << endl;
	cout << "Oblateness:                     " << (pax[1] - pax[2])/(pax[0] - pax[2]) << endl << endl;
*/
	return 0;
}

/**
@brief 	Calculates all the model links.
@param 	*model		model structure.
@return int					number of links.
**/
int			model_poly_links(Bmodel* model)
{
	int					n, nlink;
	double				d, da, ds, dat(0), dst(0);
	Bmodel*				mp;
	Blink*				link;
	
	if ( verbose )
		cout << "\nLink lengths:\nID\tAvg\tStd" << endl;
	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		da = ds = 0;
		for ( nlink=0, link = mp->link; link; link = link->next, n++, nlink++ ) {
			d = link->comp[0]->location().distance(link->comp[1]->location());
			link->FOM(d*d);
			da += d;
			ds += link->FOM();
		}
		if ( nlink ) {
			dat += da;
			dst += ds;
			da /= nlink;
			ds = sqrt(ds/nlink - da*da);
		}
		if ( verbose )
			cout << mp->identifier() << tab << da << tab << ds << endl;
	}
	
	if ( n ) {
		dat /= n;
		dst = dst/n - dat*dat;
		if ( dst > 0 ) dst = sqrt(dst);
		else dst = 0;
	}

	if ( verbose )
		cout << endl << "Overall link length:            " << dat << " (" << dst << ")" << endl << endl;
	
	return n;
}

/**
@brief 	Calculates all the polygon angles.
@param 	*model		model structure.
@return int					number of angles.

	The angles for each polygon is calculated and averaged. The overall 
	statistics for every polygon order is shown.

**/
int			model_poly_angles(Bmodel* model)
{
	int					i, n, nc, na[10], mna;
	double				a, a2, pavg, pstd, avg[10], std[10], mavg, mstd, aavg(0), astd(0);
	Vector3<double>		dcen, polycen;
	Bmodel*				mp;
	Bpolygon*			poly;
	Bcomponent*			comp;
	Bcomponent*			comp2;
	Vector3<double>		v1, v2;

	for ( i=0; i<10; i++ ) {
		na[i] = 0;
		avg[i] = std[i] = 0;
	}
	
	if ( verbose )
		cout << endl << "Angles:\nID\tVert\tAvg\tStd" << endl;
	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) comp->FOM(0);
		mavg = mstd = 0;
		mna = 0;
		for ( poly = mp->poly; poly; poly = poly->next, n++ ) {
			for ( i=0, pavg=pstd=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				if ( i+1 < poly->size() ) comp = poly->comp[i+1];
				else comp = poly->comp[0];
				if ( i+2 < poly->size() ) comp2 = poly->comp[i+2];
				else if ( i+1 < poly->size() ) comp2 = poly->comp[0];
				else comp2 = poly->comp[1];
				v1 = poly->comp[i]->location() - comp->location();
				v2 = comp2->location() - comp->location();
				a = v1.angle(v2);
				a2 = a*a;
				pavg += a;
				pstd += a2;
				avg[poly->size()] += a;
				std[poly->size()] += a2;
				na[poly->size()]++;
				mavg += a;
				mstd += a2;
				mna++;
				a = M_PI*(1 - 2.0/poly->size()) - a;
				comp->FOM(comp->FOM() + a*a);
			}
			pavg /= poly->size();
			if ( pstd ) {
				pstd = pstd/poly->size() - pavg*pavg;
				if ( pstd > 0 ) pstd = sqrt(pstd);
				else pstd = 0;
			}
		}
		aavg += mavg;
		astd += mstd;
		if ( mna ) {
			mavg /= mna;
			mstd = mstd/mna - mavg*mavg;
			if ( mstd > 0 ) mstd = sqrt(mstd);
			else mstd = 0;
		}
		for ( nc=0, comp = mp->comp; comp; comp = comp->next, nc++ )
			comp->FOM(sqrt(comp->FOM()/3));
		mp->FOM(mavg);
		if ( verbose )
			cout << mp->identifier() << tab << nc << tab << mavg*180.0/M_PI << tab << mstd*180.0/M_PI << endl;
	}

	if ( verbose )
		cout << endl << "Overall polygon angle average:\nOrder\tNumber\tAvg\tStd" << endl;
	for ( i=n=0, pavg=pstd=0; i<10; i++ ) if ( na[i] ) {
		pavg += avg[i];
		pstd += std[i];
		n += na[i];
		avg[i] /= na[i];
		if ( std[i] ) {
			std[i] = std[i]/na[i] - avg[i]*avg[i];
			if ( std[i] > 0 ) std[i] = sqrt(std[i]);
			else std[i] = 0;
		}
		if ( verbose )
			cout << i << tab << na[i] << tab << avg[i]*180.0/M_PI << tab << std[i]*180.0/M_PI << endl;
	}
	pavg /= n;
	if ( pstd ) pstd = sqrt(pstd/n - pavg*pavg);
	if ( verbose )
		cout << "All\t" << n << tab << pavg*180.0/M_PI << tab << pstd*180.0/M_PI << endl << endl;
	
	return n;
}

/**
@brief 	Analyzes a model for polygon regularity.
@param 	*model		model structure.
@return double				standard deviation from regularity.

	Regularity is defined as the adherence to a constant distance of each
	vertex from the polygon center. The polygon area is:
	     n * s^2       1 + cos(2*PI/n)
	A = ------- sqrt(-----------------)
	       4           1 - cos(2*PI/n)
	where n is the number of vertices in the polygon.
	The contribution of each polygon to the polyhedral volume is:
	V = A * dc / 3
	where dc is the distance of the polygon center to the polyhedral center.

**/
double		model_poly_regularity(Bmodel* model)
{
	int					i, n, mn;
	double				d, davg, mavg, ds, dc, dr, dr_all(0);
	double				dstd, mstd, dstd_all(0), A, Ap(0), Vp(0);
	Vector3<double>		dcen, polycen;
	Bmodel*				mp;
	Bpolygon*			poly;

	Vector3<double>		com = model_center_of_mass(model);

	if ( verbose )
		cout << "Polygon regularity:" << endl;
	if ( verbose & VERB_LABEL )
		cout << "Model\tDist\tStd" << endl;
	if ( verbose & VERB_PROCESS )
		cout << "ID\tVert\tDist\tStd\tRdist\tSpar\tRad" << endl;
	for ( n=0, mp = model; mp; mp = mp->next ) {
		mn = 0;
		mavg = mstd = 0;
		for ( poly = mp->poly; poly; poly = poly->next, n++ ) {
			polycen = 0;
			dcen = 0;
			davg = dstd = 0;
			ds = dc = 0;
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				polycen += poly->comp[i]->location();				// Center of polygon
				dcen += poly->comp[i]->location() - com;
				if ( poly->comp[i+1] ) ds += poly->comp[i]->location().distance(poly->comp[i+1]->location());
				else if ( poly->closed() ) ds += poly->comp[i]->location().distance(poly->comp[0]->location());
			}
			polycen /= i;
			dcen /= i;
			dc = dcen.length();
			ds /= i + poly->closed() - 1;
			A = i*(ds*ds/4)*sqrt((1+cos(TWOPI/i))/(1-cos(TWOPI/i)));	// Polygon area
			Ap += A;
			Vp += A*dc/3;	// Contribution to polyhedral volume
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				d = polycen.distance(poly->comp[i]->location());
				davg += d;
				dstd += d*d;
			}
			mavg += davg;
			mstd += dstd;
			mn += i;
			davg /= i;
			dstd = dstd/i - davg*davg;
//			dr = davg*sqrt(2*(1 - cos(TWOPI/i)));	// Distance relative to the spar
			dr = davg*2*sin(M_PI/i);				// Distance relative to the spar
			dr_all += dr;
			if ( dstd > 0 ) {
				dstd_all += dstd;
				dstd = sqrt(dstd);
			} else dstd = 0;
			if ( verbose & VERB_PROCESS )
				cout << n << tab << i << tab << davg << tab << dstd << tab << dr << tab << ds << tab << dc << endl;
		}
		if ( mn ) mavg /= mn;
		if ( mstd ) mstd = sqrt(mstd/mn - mavg*mavg);
		if ( verbose & VERB_LABEL )
			cout << mp->identifier() << tab << mavg << tab << mstd << endl;
	}
	
	dstd_all = sqrt(dstd_all/n);
	dr_all /= n;
	cout << "Regularity deviation:           " << dstd_all << " A (" << n << ")" << endl;
	cout << "Overall spar length:            " << dr_all << " A" << endl;
	cout << "Polyhedral volume:              " << Vp << " A3" << endl << endl;
	
	return dstd_all;
}

/**
@brief 	Analyzes a model for polygon planarity.
@param 	*model		model structure.
@return double				standard deviation from planarity.

	A plane is fit through the polygon vertices and the normal calculated from:
		n•p = d
	where n is the normal vector, p is a point in the plane, and d is the offset.
	The polygon planarity is defined as the root-mean-square-deviation from 
	the fitted plane.

**/
double		model_poly_planarity(Bmodel* model)
{
	int					i, j, m, n;
	double				offset, d, dstd, dstd_all(0);
	vector<double>		b(3);
	Matrix				a(3,3);
	Vector3<double>		polycen, loc, normal;
	Bmodel*				mp;
	Bpolygon*			poly;

	cout << "Polygon planarity:\nID\tVert\tStd" << endl;
	dstd_all = 0;
	for ( n=0, mp = model; mp; mp = mp->next ) {
		for ( poly = mp->poly; poly; poly = poly->next, n++ ) {
			for ( i=0; i<3; i++ ) b[i] = 0;
//			for ( i=0; i<9; i++ ) a[i] = 0;
			a.fill(0);
			polycen = 0;
			dstd = 0;
			for ( m=0; m<poly->comp.size() && poly->comp[m]; m++ ) {
				loc = poly->comp[m]->location();
				polycen += loc;
				for ( i=0; i<3; i++ ) {
					for ( j=0; j<=i; j++ )
						a[j][i] += loc[i]*loc[j];
					b[i] += loc[i];
				}
			}
			polycen /= m;
			for ( i=1; i<3; i++ )
				for ( j=0; j<i; j++ )
					a[i][j] = a[j][i];
			offset = 0;
			if ( a[0][0] == 0 ) {
				normal = Vector3<double>(1, 0, 0);
			} else if ( a[1][1] == 0 ) {
				normal = Vector3<double>(0, 1, 0);
			} else if ( a[2][2] == 0 ) {
				normal = Vector3<double>(0, 0, 1);
			} else {
				a.LU_decomposition(b);
				normal = Vector3<double>(b[0], b[1], b[2]);
				normal.normalize();
				offset = polycen.scalar(normal);
			}
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				d = poly->comp[i]->location().scalar(normal) - offset;
				dstd += d*d;
			}
			dstd /= i;
			dstd_all += dstd;
			dstd = sqrt(dstd);
			cout << n << tab << m << tab << dstd << endl;
		}
	}
	dstd_all = sqrt(dstd_all/n);
	cout << "Planarity deviation:            " << dstd_all << " A (" << n << ")" << endl << endl;	

	return dstd_all;
}

/**
@brief 	Calculates the different energy terms for all models.
@param 	*model		model structure.
@param 	angle_ref		reference angle (<=0 to use the polygon angle).
@return double				0.

	The angular energy is calculated either with a given reference angle,
	or with the nominal polygon inner angle as reference.
	Regularity is defined as the adherence to a constant distance of each
	vertex from the polygon center. The polygon area is:
	     n * s^2       1 + cos(2*PI/n)
	A = ------- sqrt(-----------------)
	       4           1 - cos(2*PI/n)
	where n is the number of vertices in the polygon.
	The contribution of each polygon to the polyhedral volume is:
	V = A * dc / 3
	where dc is the distance of the polygon center to the polyhedral center.

**/
double		model_poly_energy(Bmodel* model, double angle_ref)
{
	int					i, j, m, n, nmod, ncomp;
	double				d, ds, dr;
	double				El, Ea, Er, Ep;
	double				Elt(0), Eat(0), Ert(0), Ept(0);
	double				offset;
	vector<double>		b(3);
	Matrix				a(3,3);
	Vector3<double>		polycen, loc, normal, v1, v2;
	Bmodel*				mp;
	Bcomponent*			comp;
	Bcomponent*			comp2;
	Blink*				link;
	Bpolygon*			poly;

	if ( verbose )
		cout << "Model ID\t#vert\tElink   \tEangle   \tEpolyreg\tEpolyplan" << endl;
	for ( nmod = 0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		El = Ea = Er = Ep = 0;
		for ( ncomp=0, comp = mp->comp; comp; comp = comp->next ) ncomp++;
		for ( n=0, link = mp->link; link; link = link->next ) {
			d = link->length() - link->comp[0]->location().distance(link->comp[1]->location());
			El += d*d;
			n++;
		}
		if ( El ) El = sqrt(El/n);
		for ( n=0, poly = mp->poly; poly; poly = poly->next ) {
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				if ( i+1 < poly->size() ) comp = poly->comp[i+1];
				else comp = poly->comp[0];
				if ( i+2 < poly->size() ) comp2 = poly->comp[i+2];
				else if ( i+1 < poly->size() ) comp2 = poly->comp[0];
				else comp2 = poly->comp[1];
				v1 = poly->comp[i]->location() - comp->location();
				v2 = comp2->location() - comp->location();
				if ( angle_ref > 0 ) d = angle_ref - v1.angle(v2);
				else d = M_PI*(1 - 2.0/poly->size()) - v1.angle(v2);
				Ea += d*d;
				n++;
			}
		}
		if ( Ea ) Ea = sqrt(Ea/n);
		for ( n=0, poly = mp->poly; poly; poly = poly->next ) {
			polycen = 0;
			ds = 0;
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				polycen += poly->comp[i]->location();				// Center of polygon
				if ( poly->comp[i+1] ) ds += poly->comp[i]->location().distance(poly->comp[i+1]->location());
				else if ( poly->closed() ) ds += poly->comp[i]->location().distance(poly->comp[0]->location());
			}
			polycen /= i;
			ds /= i + poly->closed() - 1;
			dr = ds/(2*sin(M_PI/i));
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				d = dr - polycen.distance(poly->comp[i]->location());
				Er += d*d;
				n++;
			}
		}
		if ( Er ) Er = sqrt(Er/n);
		for ( n=0, poly = mp->poly; poly; poly = poly->next ) {
			for ( i=0; i<3; i++ ) b[i] = 0;
//			for ( i=0; i<9; i++ ) a[i] = 0;
			a.fill(0);
			polycen = 0;
			for ( m=0; m<poly->comp.size() && poly->comp[m]; m++ ) {
				loc = poly->comp[m]->location();
				polycen += loc;
				for ( i=0; i<3; i++ ) {
					for ( j=0; j<=i; j++ )
						a[j][i] += loc[i]*loc[j];
					b[i] += loc[i];
				}
			}
			polycen /= m;
			for ( i=1; i<3; i++ )
				for ( j=0; j<i; j++ )
					a[i][j] = a[j][i];
			offset = 0;
			if ( a[0][0] == 0 ) {
				normal = Vector3<double>(1, 0, 0);
			} else if ( a[1][1] == 0 ) {
				normal = Vector3<double>(0, 1, 0);
			} else if ( a[2][2] == 0 ) {
				normal = Vector3<double>(0, 0, 1);
			} else {
				a.LU_decomposition(b);
				normal = Vector3<double>(b[0], b[1], b[2]);
				normal.normalize();
				offset = polycen.scalar(normal);
			}
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
				d = poly->comp[i]->location().scalar(normal) - offset;
				Ep += d*d;
				n++;
			}
		}
		if ( Ep ) Ep = sqrt(Ep/n);
		if ( verbose )
			cout << mp->identifier() << tab << ncomp << tab << El << tab << Ea << tab << Er << tab << Ep << endl;
		Elt += El;
		Eat += Ea;
		Ert += Er;
		Ept += Ep;
		nmod++;
	}
	
	if ( nmod ) {
		Elt /= nmod;
		Eat /= nmod;
		Ert /= nmod;
		Ept /= nmod;
	}
	
	if ( verbose )
		cout << "Average  \t \t" << Elt << tab << Eat << tab << Ert << tab << Ept << endl << endl;
	
	return 0;
}

/**
@brief 	Calculates the number of edges shared by pentagons.
@param 	*model		model structure.
@return int					number of edges shared by pentagons for last model.
**/
int			model_poly_pentagon_adjacency(Bmodel* model)
{
	if ( !model->poly ) model_poly_generate(model);

	int				i, n(0), n3, nn;
	Bmodel*			mp;
	Bcomponent*		comp;
	Blink*			link;
	Bpolygon*		poly;
	
	if ( verbose )
		cout << "Pentagon adjacency:\nModel\tNp\tN3\tNN" << endl;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) comp->select(0);
		for ( poly = mp->poly; poly; poly = poly->next ) if ( poly->size() == 5 ) {
			for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ )
				poly->comp[i]->select_increment();
		}
		for ( nn=0, link = mp->link; link; link = link->next )
			if ( link->comp[0]->select() > 2 && link->comp[1]->select() > 2 ) nn++;
		for ( n=n3=0, comp = mp->comp; comp; comp = comp->next ) {
			if ( comp->select() > 2 ) {
				n += comp->select();
				n3++;
			} else if ( comp->select() == 2 ) n++;
			comp->select(1);
		}
		n /= 2;
		if ( verbose )
			cout << mp->identifier() << tab << n << tab << n3 << tab << nn << endl;
	}

	if ( verbose )
		cout << endl;
	
	return n;
}

#define	NSMAX		100

/**
@brief 	Finds the symmetry axes associated with polyhedral elements.
@param 	*model		model structure.
@param 	threshold	cutoff to flag a symmetry axis.
@return Bstring		symmetry string.

	Symmetry axes are associated with specific elements:
		link		reflection and 2-fold axis.
		vertex		n-fold axis based on vertex order.
		polygon		n-fold axis based on polygon order.
	Only the first model in the list is processed.

**/
Bstring		model_poly_find_symmetry(Bmodel* model, double threshold)
{
	Bstring				sym_label("C1");
	if ( !model ) return sym_label;
	
	if ( !model->poly ) model_poly_generate(model);
	
	long			i, j, n, ns(0), npf, order;
	double			R;
	long*			pf;
	Bcomponent*		comp;
	Blink*			link;
	Bpolygon*		poly;
	Transform		t;
	Vector3<double>	polycen;

	Bsymop*			op = new Bsymop[NSMAX];
	
	t.origin = model_center_of_mass(model);
	t.angle = M_PI;
	
	double			Gr = model_gyration_radius(model);
	
	// Reflections
	if ( verbose & VERB_FULL )
		cout << "Link\tOrder\tSym\tR" << endl;
	for ( n=1, link = model->link; link; link = link->next, n++ ) {
		t.axis = link->comp[1]->location() - link->comp[0]->location();
		R = model_reflect_and_compare(model, t.axis, t.origin);
		R /= Gr;
		if ( R < threshold ) {
			if ( verbose & VERB_FULL )
				cout << n << "\t1\t1\t" << R << endl;
			op[ns].axis(t.axis[0], t.axis[1], t.axis[2]);
			op[ns].order(1);
			ns++;
		}
		t.axis = t.axis.cross(link->comp[0]->location() - t.origin);
		R = model_reflect_and_compare(model, t.axis, t.origin);
		R /= Gr;
		if ( R < threshold ) {
			if ( verbose & VERB_FULL )
				cout << n << "\t1\t1\t" << R << endl;
			op[ns].axis(t.axis[0], t.axis[1], t.axis[2]);
			op[ns].order(1);
			ns++;
		}
	}

	// Rotations
	for ( n=1, link = model->link; link; link = link->next, n++ ) {
		t.axis = (link->comp[0]->location() + link->comp[1]->location())/2.0 - t.origin;
		R = model_rotate_and_compare(model, t);
		R /= Gr;
		if ( R < threshold ) {
			if ( verbose & VERB_FULL )
				cout << n << "\t2\t2\t" << R << endl;
			op[ns].axis(t.axis[0], t.axis[1], t.axis[2]);
			op[ns].order(2);
			op[ns].angle(t.angle);
			ns++;
		}
	}
		
	if ( verbose & VERB_FULL )
		cout << "Vertex\tOrder\tSym\tR" << endl;
	for ( n=1, comp = model->comp; comp; comp = comp->next, n++ ) {
		for ( i=0; i<comp->link.size() && comp->link[i]; i++ ) ;
		if ( i < 1 ) {
			if ( verbose & VERB_FULL )
				cerr << "No links for component " << comp->identifier() << endl;
		} else {
			t.axis = comp->location() - t.origin;
			t.angle = M_PI*2.0/i;
			pf = prime_factors(i, npf);
			order = 0;
			for ( j=0; j<npf; j++ ) {
				t.angle = M_PI*2.0/pf[j];
				R = model_rotate_and_compare(model, t);
				R /= Gr;
				if ( R < threshold ) {
					if ( verbose & VERB_FULL )
						cout << n << tab << i << tab << pf[j] << tab << R << endl;
					if ( order < pf[j] ) order = pf[j];
				}
			}
			if ( npf && ( pf[0] != i ) ) {
				t.angle = M_PI*2.0/i;
				R = model_rotate_and_compare(model, t);
				R /= Gr;
				if ( R < threshold ) {
					if ( verbose & VERB_FULL )
						cout << n << tab << i << tab << i << tab << R << endl;
					order = i;
				}
			}
			delete[] pf;
			if ( order ) {
				op[ns].axis(t.axis[0], t.axis[1], t.axis[2]);
				op[ns].order(order);
				op[ns].angle(M_PI*2.0/order);
				ns++;
			}
		}
	}

	if ( verbose & VERB_FULL )
		cout << "Polygon\tOrder\tSym\tR" << endl;
	for ( n=1, poly = model->poly; poly; poly = poly->next, n++ ) {
		polycen = 0;
		for ( i=0; poly->comp[i]; i++ )
			polycen += poly->comp[i]->location();
		polycen /= i;
		t.axis = polycen - t.origin;
		t.axis.normalize();
		pf = prime_factors(i, npf);
		order = 0;
		for ( j=0; j<npf; j++ ) {
			t.angle = M_PI*2.0/pf[j];
			R = model_rotate_and_compare(model, t);
			R /= Gr;
			if ( R < threshold ) {
				if ( verbose & VERB_FULL )
					cout << n << tab << i << tab << pf[j] << tab << R << endl;
				if ( order < pf[j] ) order = pf[j];
			}
		}
		if ( pf[0] != i ) {
			t.angle = M_PI*2.0/i;
			R = model_rotate_and_compare(model, t);
			R /= Gr;
			if ( R < threshold ) {
				if ( verbose & VERB_FULL )
					cout << n << tab << i << tab << i << tab << R << endl;
				order = i;
			}
		}
		delete[] pf;
		if ( order ) {
			op[ns].axis(t.axis[0], t.axis[1], t.axis[2]);
			op[ns].order(order);
			op[ns].angle(M_PI*2.0/order);
			ns++;
		}
	}

	if ( verbose & VERB_FULL )
		cout << endl;

	for ( i=0; i<ns; i++ ) op[i].normalize();
	
	double			a, maxang = M_PI/18.0;
	int				na(0);
	vector<int>		nord(NSMAX,0);
	
	// Eliminate redundant axes
	for ( i=1; i<ns; i++ ) {
		for ( j=0; j<i; j++ ) {
			if ( op[i].order() == op[j].order() ) {
				a = fabs(op[i].axis().angle(op[j].axis()));
				if ( a < maxang ) op[i].order(0);
				else if ( fabs(M_PI - a) < maxang ) op[i].order(0);
			}
		}
	}

	for ( i=0; i<ns; i++ ) if ( op[i].order() ) nord[op[i].order()]++;
	
	if ( verbose & VERB_FULL ) {
		cout << "Order\tx\ty\tz" << endl;
		for ( i=0; i<ns; i++ ) if ( op[i].order() )
			cout << op[i].order() << tab << op[i].axis() << endl;
		cout << endl;
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Order\tAxes" << endl;
		for ( i=1; i<NSMAX; i++ )
			if ( nord[i] ) cout << i << tab << nord[i] << endl;
		cout << endl;
	}
	
	for ( i=1; i<NSMAX; i++ )
		if ( nord[i] ) na += nord[i];

	Vector3<double>		vecz(0,0,1), vecx(1,0,0);

	if ( nord[5] > 3 && nord[3] > 5 && nord[2] > 7 ) {
		sym_label = "I";
		if ( nord[1] ) sym_label = "Ih";
		for ( i=0; i<ns && op[i].order() != 2; i++ ) ;
		if ( op[i].order() == 2 ) {
			vecz = op[i].axis();
			for ( j=0; j<ns; j++ ) if ( op[j].order() == 5 )
				if ( op[i].axis().angle(op[j].axis()) < M_PI/5.0 )
					vecx = vecz.cross(op[j].axis());
		}
	} else if ( nord[4] > 1 && nord[3] > 2 && nord[2] > 3 ) {
		sym_label = "O";
		if ( nord[1] ) sym_label = "Oh";
		for ( i=0; i<ns && op[i].order() != 4; i++ ) ;
		if ( op[i].order() == 4 ) {
			vecz = op[i].axis();
			for ( j=0; j<ns; j++ ) if ( op[j].order() == 4 )
				if ( fabs(M_PI_2 - op[i].axis().angle(op[j].axis())) < maxang )
					vecx = op[j].axis();
		}
	} else if ( nord[3] > 2 && nord[2] > 1 ) {
		sym_label = "T";
		for ( i=0; i<ns && op[i].order() != 2; i++ ) ;
		if ( op[i].order() == 2 ) {
			vecz = op[i].axis();
			for ( j=0; j<ns; j++ ) if ( op[j].order() == 2 )
				if ( fabs(M_PI_2 - op[i].axis().angle(op[j].axis())) < maxang )
					vecx = op[j].axis();
		}
		if ( nord[1] ) {
			for ( j=0; j<ns; j++ ) if ( op[j].order() == 1 )
				if ( fabs(M_PI_2 - vecz.angle(op[j].axis())) < maxang )
					sym_label = "Td";
			for ( j=0; j<ns; j++ ) if ( op[j].order() == 1 ) {
				a = fabs(vecz.angle(op[j].axis()));
				if ( a < maxang || fabs(M_PI - a) < maxang )
					sym_label = "Th";
			}
		}
	} else if ( nord[2] > 1 ) {
		for ( order=0, i=2; i<NSMAX; i++ ) if ( nord[i] ) order = i;
		sym_label = Bstring(order, "D%d");
		// D2 is special because it can be oriented 3 different ways
		if ( order == 2 && nord[1] ) {
			for ( j=0; j<ns && op[j].order() != 1; j++ ) ;
			for ( i=0; i<ns; i++ ) if ( op[i].order() == 2 )
				if ( fabs(op[i].axis().angle(op[j].axis())) < maxang ||
					fabs(M_PI_2 - op[i].axis().angle(op[j].axis())) < maxang )
						break;
		} else {
			for ( i=0; i<ns && op[i].order() != order; i++ ) ;
		}
		if ( op[i].order() == order ) {
			vecz = op[i].axis();
			for ( j=0; j<ns; j++ ) if ( op[j].order() == 2 )
				if ( fabs(M_PI_2 - op[i].axis().angle(op[j].axis())) < maxang )
					vecx = op[j].axis();
		}
		if ( nord[1] ) {
			for ( j=0; j<ns; j++ ) if ( op[j].order() == 1 )
				if ( fabs(M_PI_2 - vecz.angle(op[j].axis())) < maxang )
					sym_label = Bstring(order, "D%dd");
			for ( j=0; j<ns; j++ ) if ( op[j].order() == 1 ) {
				a = fabs(vecz.angle(op[j].axis()));
				if ( a < maxang || fabs(M_PI - a) < maxang )
					sym_label = Bstring(order, "D%dh");
			}
		}
	} else {
		sym_label = "C1";
		for ( order=0, i=1; i<NSMAX; i++ ) if ( nord[i] ) order = i;
		if ( order ) {
			if ( order == 1 ) sym_label = "Cs";
			else sym_label = Bstring(order, "C%d");
			for ( i=0; i<ns && op[i].order() != order; i++ ) ;
			if ( op[i].order() ) {
				vecx = op[i].axis().cross(vecz);
				vecz = op[i].axis();
				vecx.normalize();
				// Use any reflection planes to further define orientation
				if ( order>1 && nord[1] ) {
					for ( j=0; j<ns; j++ ) if ( op[j].order() == 1 ) {
						if ( fabs(M_PI_2 - vecz.angle(op[j].axis())) < maxang ) {
							vecx = op[j].axis();
							sym_label = Bstring(order, "C%dv");
						}
					}
					for ( j=0; j<ns; j++ ) if ( op[j].order() == 1 ) {
						a = fabs(vecz.angle(op[j].axis()));
						if ( a < maxang || fabs(M_PI - a) < maxang )
							sym_label = Bstring(order, "C%dh");
					}
				}
			}
		}
	}
	
	t.axis = Vector3<double>(vecz[1], -vecz[0], 0);
	t.axis.normalize();
	t.angle = acos(vecz[2]);
	Matrix3			mat = Matrix3(t.axis, t.angle);
//	cout << mat << endl;
	t = Transform(mat);
	vecx = mat * vecx;
	t.axis = Vector3<double>(0, 0, 1);
	t.angle = atan2(-vecx[1], vecx[0]);
	mat = Matrix3(t.axis, t.angle) * mat;
	t = Transform(mat);
	
	delete[] op;
	
	string			symstr(sym_label.str());
	regex			rmstr("[svdh]");
	model->symmetry(symstr);
	
	model->handedness(1);
//	if ( model->symmetry().contains("s") || model->symmetry().contains("v") ||
//			model->symmetry().contains("d") || model->symmetry().contains("h") )
	if ( regex_match(symstr, rmstr) ) model->handedness(0);
	else model_poly_hand(model);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Point group:                    " << sym_label << endl;
		cout << "Hand:                           " << model->handedness() << endl;
		cout << "Z-axis:                         " << vecz << endl;
		cout << "X-axis:                         " << vecx << endl << endl;
	}
	
	model_rotate(model, t);
	
	return sym_label;
}

/**
@brief 	Calculates the polyhedron hand.
@param 	*model		model structure.
@return int					hand sign.

	The hand of a polyhedron is based on the sign of the rotational strength 
	element associated with the first two eigenvectors of the adjacency matrix.
	It is assumed that the symmetry has been determined.
	The hand sign is:
		0	no handedness.
		1	one enantiomorph.
		-1	other enantiomorph.

**/
int			model_poly_hand(Bmodel* model)
{
	if ( model->symmetry().length() < 1 )
		model_poly_find_symmetry(model, 0.1);

	model->handedness(0);

	regex			rmstr("[svdh]");
	if ( regex_match(model->symmetry(), rmstr) )
		return model->handedness();

	model->handedness(1);

	int				i, j, k, ncomp;
	Bcomponent*		comp;
	
	for ( ncomp = 0, comp = model->comp; comp; comp = comp->next ) ncomp++;
	
	if ( ncomp < 2 ) return 0;
	
	Matrix			a = model_adjacency_matrix(model);

	vector<double>	d = a.jacobi_rotation();
	a.eigen_sort(d);
	
	Matrix			a2 = model_adjacency_matrix(model);

	double				R;
	Vector3<double>		v1, v2;
	Vector3<double>*	coor = new Vector3<double>[ncomp];

	for ( i=0, comp = model->comp; comp; comp = comp->next, i++ ) {
		coor[i] = comp->location();
		coor[i].normalize();
	}
	
	for ( i=j=0, v1=v2=0; i<ncomp; i++, j+=ncomp ) {
		v1 += coor[i]*(a[i][0]*a[i][1]);
		for ( k=0; k<ncomp; k++ )
			if ( a2[i][k] )
				v2 += coor[i].cross(coor[k])*(a[i][0]*a[k][1]);
	}
	R = v1.scalar(v2);
	
	if ( R > 0 ) model->handedness(1);
	else model->handedness(-1);
	
	delete[] coor;
	
	return model->handedness();
}

/**
@brief 	Compares a model with reference models based on the eigenvalues of the adjacency matrix.
@param 	*model			model structure.
@param 	*refmodel		reference model(s) to compare.
@return int						0.

	The comparison is based on the eigenvectors of the adjacency matrix,
	which are related to spherical harmonics.
	There are 3 P(sigma) eigenvectors giving the vertex coordinates.
	These are usually (but not always) vectors 2, 3, and 4 ordered by eigenvalue.
	The eigenvalues are characteristic for a polyhedron, although they may not be unique.
	For every model identified, the model type is set from the reference ID.
	The reference model selection is incremented to indicate the count. 

	Dover Publications, Inc., Mineola, New York, pages 101 - 104.
Reference: 	Fowler, P.W. and Manolopoulos, D.E. (2006) An Atlas of Fullerenes. 

**/
int			model_poly_compare(Bmodel* model, Bmodel* refmodel)
{
	int			i, j, nref, n, d, nmod(0), ntot(0);
	Bmodel*		mp;
	Bmodel*		mp_ref;
	Bmodel*		mp_pick;

	// Get all the reference eigenvalue arrays
	for ( nref=0, mp_ref = refmodel; mp_ref; mp_ref = mp_ref->next )
		if ( mp_ref->select() ) nref++;
	
	if ( verbose ) 
		cout << "Comparing models with " << nref << " references" << nref << endl << endl;
	
//	int*		nr = new int[nref];
//	double**	vr = new double*[nref];
	vector<int>	nr(nref);
	vector<vector<double>>	vr(nref);
	
	for ( i = 0, mp_ref = refmodel; mp_ref; mp_ref = mp_ref->next ) if ( mp_ref->select() ) {
		nr[i] = count_list((char *)mp_ref->comp);
		vr[i] = model_poly_eigenvalues(mp_ref, 0);
		mp_ref->FOM(0);
		i++;
	}
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		mp_pick = NULL;
		n = count_list((char *)mp->comp);
		if ( n ) {
			model_poly_find_symmetry(mp, 0.1);
			vector<double>	v = model_poly_eigenvalues(mp, 0);
			for ( i = 0, mp_ref = refmodel; mp_ref && mp_pick == NULL; mp_ref = mp_ref->next ) if ( mp_ref->select() ) {
				if ( n == nr[i] ) {
					for ( j = 0, d = 0; j<n; j++ )
						if ( fabs(v[j] - vr[i][j]) > 1e-6 ) d++;
					if ( d < 1 ) {
						mp_pick = mp_ref;
						if ( verbose & VERB_FULL )
							cout << "Test: " << mp->identifier() << tab << mp_ref->identifier() << tab << 
								mp->symmetry() << tab << mp_ref->symmetry() << endl;
					}
				}
				i++;
			}
		}
		if ( mp_pick ) {
			cout << mp->identifier() << tab << n << tab << mp->symmetry() << tab << mp->handedness() << tab <<
				mp_pick->identifier() << tab << mp_pick->symmetry() << endl;
			mp->model_type(mp_pick->identifier());
			mp_pick->FOM(mp_pick->FOM() + 1);
			mp->FOM(1);
			nmod++;
		} else {
			cout << mp->identifier() << tab << n << tab << mp->symmetry() << "\t-\t-" << endl;
			mp->FOM(0);
		}
		ntot++;
	}

	for ( mp = model; mp; mp = mp->next )
		if ( mp->select() ) mp->select((int) mp->FOM());

	cout << "Models identified:              " << nmod << "/" << ntot << " (" << nmod*100.0/ntot << ")" << endl << endl;
	
//	for ( i = 0; i<nref; i++ ) delete[] vr[i];
//	delete[] nr;
//	delete[] vr;
		
	return 0;
}

/**
@brief 	Generates the eigenvalues of the adjacency matrix for a model.
@param 	*model			model structure (modified with the topological coordinates).
@param 	show			flag to show eigenvalues.
@return vector<double>	eigenvalues.

	The eigenvectors of the adjacency matrix are related to spherical harmonics.
	The eigenvalues are characteristic for a polyhedron, although they may not be unique.
	Only the first model in the list is processed.

	Dover Publications, Inc., Mineola, New York, pages 101 - 104.
Reference: 	Fowler, P.W. and Manolopoulos, D.E. (2006) An Atlas of Fullerenes. 

**/
vector<double>	model_poly_eigenvalues(Bmodel* model, int show)
{
	int				i, j, k, ncomp;
	vector<double>	d;
	Bcomponent*		comp;
	
	for ( ncomp = 0, comp = model->comp; comp; comp = comp->next ) ncomp++;
	
	if ( ncomp < 1 ) return d;
	
	Matrix			a = model_adjacency_matrix(model);

	d = a.jacobi_rotation();
	a.eigen_sort(d);
	
	if ( verbose & VERB_FULL ) {
		for ( k=i=0; i<ncomp; i++ ) {
			cout << i << ": ";
			for ( j=0; j<ncomp; j++, k++ )
				cout << " " << a[i][j];
			cout << ": " << d[i] << endl;
		}
	}
	
	if ( ( verbose & (VERB_PROCESS|VERB_FULL) ) && show ) {
		cout << "Eigenvalues:";
		for ( i=0; i<ncomp; i++ ) cout << tab << d[i];
		cout << endl;
	}
	
	return d;
}

// Recursively looks for all linked components with the same selection setting
int			comp_count_connected(Bcomponent* comp)
{
	int			i, n = 1;
	int			s = comp->select();
	comp->select(0);
	
	for ( i=0; i<comp->link.size() && comp->link[i]; i++ ) if ( comp->link[i]->select() == s )
		n += comp_count_connected(comp->link[i]);

	return n;
}

int			model_number_connected_clusters(Bmodel* model)
{
	int			i, c(0);
	Bcomponent*	comp;
	
//	if ( verbose & VERB_FULL )
		cout << "Connected clusters:" << endl;
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		i = comp_count_connected(comp);
		c++;
//		if ( verbose & VERB_FULL )
			cout << c << tab << i << endl;
	}
	
	return c;
}

/**
@brief 	Generates coordinates for components based on the eigenvectors of the adjacency matrix.
@param 	*model			model structure.
@return dvector<double>	eigenvalues.

	The eigenvectors of the adjacency matrix are related to spherical harmonics.
	There are 3 P(sigma) eigenvectors giving the vertex coordinates.
	These are usually (but not always) vectors 2, 3, and 4 ordered by eigenvalue.
	The actual 3 P(sigma) eigenvectors are identified as those having exactly
	a single node, i.e., closely connected vertices cluster together in every dimension.
	The eigenvalues are characteristic for a polyhedron, although they may not be unique.

	Dover Publications, Inc., Mineola, New York, pages 101 - 104.
Reference: 	Fowler, P.W. and Manolopoulos, D.E. (2006) An Atlas of Fullerenes. 

**/
vector<double>	model_poly_sphere_coor(Bmodel* model)
{
	int				i, j, k, m, ncomp;
	double			len, s[3] = {1,1,1};
	Bcomponent*		comp;
	Blink*			link;
	
	for ( ncomp = 0, comp = model->comp; comp; comp = comp->next ) ncomp++;
	
	Matrix			a = model_adjacency_matrix(model);

	vector<double>	d = a.jacobi_rotation();
	a.eigen_sort(d);
	
	int				p[3] = {1,2,3};
	
	if ( ncomp > 50 ) {
		for ( i=1, m=0; i<ncomp && m<3; i++ ) {
			for ( j=0, comp = model->comp; j<ncomp && comp; j++, comp = comp->next ) {
				if ( a[j][i] > 0 ) comp->select(1);
				else comp->select(-1);
			}
			k = model_number_connected_clusters(model);
			if ( k == 2 ) p[m++] = i;
			if ( verbose & VERB_FULL )
				cout << "Eigenvector " << i << ": clusters =" << k << endl;
		}
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Eigenvectors used:";
		for ( i=0; i<3; i++ ) cout << tab << p[i];
		cout << endl;
	}
	
	for ( i=0; i<3; i++ ) s[i] = 1/sqrt(d[0] - d[p[i]]);
	
	for ( i=0, comp = model->comp; comp; comp = comp->next, i++ )
		comp->location(Vector3<double>(s[0]*a[i][p[0]], s[1]*a[i][p[1]], s[2]*a[i][p[2]]));

	for ( i=0, len=0, link = model->link; link; link = link->next, i++ )
		len += link->comp[0]->location().distance(link->comp[1]->location());
	
	len = i/len;
	
	for ( comp = model->comp; comp; comp = comp->next )
		comp->scale(len);

	return d;
}

