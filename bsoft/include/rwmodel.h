/**
@file	rwmodel.h
@brief	Header file for reading and writing atomic model files
@author Bernard Heymann
@date	Created: 20060919
@date	Modified: 20220929
**/

#include "View2.h"
#include "Bstring.h"
#include "Color.h"

#include <map>

#ifndef _Bmodel_
#define	MAXLINK		10				// Maximum number of links per component

/************************************************************************
@Object: class Bcomptype
@Description:
	Model component type parameter structure.
@Features:
	Reference to a model file, either atomic coordinates or a 3D map.
	Number of components of this type.
	Figure-of-merit.
*************************************************************************/
class Bcomptype {
private:
	void	initialize() {
		next = NULL;
		id = "?";
		ind = 0;
		fmod = "?";
		num = 0;
		cnt = 0;
		mas = 1;
		chrg = 0;
		lrad = 1;
		rad = 1;
		fom = 0;
		sel = 1;
		for ( int i=0; i<11; ++i ) coef[i] = 0;
	}
public:
	Bcomptype*		next;			// Next component type in list
private:
	string			id;				// Component type identifier
	int				ind;			// Type index for reference in a distance matrix
	string			fmod;			// Model file name (molecule or map)
	int				num;			// Image number in model or map file
	float			cnt;			// Number of components
	float			mas;			// Component mass
	float			chrg;			// Component charge
	float			lrad;			// Link radius (bond radius; half link length)
	float			rad;			// Characteristic radius (VdW radius)
	float			coef[11];		// Cross section coefficients
	float			fom;			// Figure-of-merit
	int				sel;			// Selection flag
public:
	Bcomptype() { initialize(); }
	Bcomptype(string s) { initialize(); id = s; }
	Bcomptype(long i) { initialize(); id = to_string(i); }
	Bcomptype(string fn, long img_num) { initialize(); id = "1"; fmod = fn; num = img_num; }
	Bcomptype(string s, string fn, long img_num) { initialize(); id = s; fmod = fn; num = img_num; }
	Bcomptype(Bcomptype* ct) {
		next = NULL;
		id = ct->id;
		ind = ct->ind;
		fmod = ct->fmod;
		num = ct->num;
		cnt = ct->cnt;
		mas = ct->mas;
		chrg = ct->chrg;
		lrad = ct->lrad;
		rad = ct->rad;
		fom = ct->fom;
		sel = ct->sel;
		coefficients(ct->coefficients());
//		cout << "In Bcomptype:" << endl;
//		show();
	}
	void			identifier(string s) { id = s; }
	void			identifier(Bstring s) { id = s.str(); }
	string&			identifier() { return id; }
	void			index(long i) { ind = i; }
	long			index() { return ind; }
	long			index(string s) {
		Bcomptype*		c = find(s);
		if ( c ) return c->ind;
		else return -1;
	}
	void			file_name(string s) { fmod = s; }
	string&			file_name() { return fmod; }
	void			image_number(long i) { num = i; }
	long			image_number() { return num; }
	void			mass(double m) { mas = m; }
	double			mass() { return mas; }
	void			charge(double c) { chrg = c; }
	double			charge() { return chrg; }
	void			link_radius(double r) { lrad = r; }
	double			link_radius() { return lrad; }
	void			radius(double r) { rad = r; }
	double			radius() { return rad; }
	template <typename T>
	void			coefficients(vector<T> c) {
		for ( int i=0; i<11 && i<c.size(); ++i ) coef[i] = c[i];
	}
	vector<double>	coefficients() {
		vector<double>	c(11);
		for ( int i=0; i<11; ++i ) c[i] = coef[i];
		return c;
	}
	void			FOM(double d) { fom = d; }
	double			FOM() { return fom; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
	long			select_increment() { return sel++; }
	void			component_count(float n) { cnt = n; }
	float			component_count() { return cnt; }
	void			component_count_increment() { cnt++; }
	void			component_count_decrement() { cnt--; }
	void			component_count_add(float n) { cnt += n; }
	Bcomptype*		add(string s) {
		Bcomptype* 		c(this);
		while ( c->id != s && c->next ) c = c->next;
		if ( c->next ) return c;
		return c->next = new Bcomptype(s);
	}
	Bcomptype*		add(long i) {
		Bcomptype* 		c(this);
		string			s(to_string(i));
		while ( c->id != s && c->next ) c = c->next;
		if ( c->next ) return c;
		return c->next = new Bcomptype(s);
	}
	Bcomptype*		add(string fn, long img_num) {
		Bcomptype* 		c(this);
		while ( !( c->fmod == fn && c->num == img_num ) && c->next ) c = c->next;
		if ( c->fmod == fn && c->num == img_num ) return c;
		return c->next = new Bcomptype(fn, img_num);
	}
	Bcomptype*		add(string s, string fn, long img_num) {
		Bcomptype* 		c(this);
		while ( !( c->id == s && c->fmod == fn && c->num == img_num ) && c->next ) c = c->next;
		if ( c->fmod == fn && c->num == img_num ) return c;
		return c->next = new Bcomptype(s, fn, img_num);
	}
	Bcomptype*		add(Bcomptype* ct) {
		Bcomptype* 		c(this);
		while ( c->next ) c = c->next;
		return c->next = new Bcomptype(ct);
	}
	Bcomptype*		find(string s) {
		Bcomptype* 		c(this);
		while ( c && c->id != s ) c = c->next;
		return c;
	}
	long			count() {
		long			n(0);
		for ( Bcomptype* c=this; c; c=c->next ) n++;
		return n;
	}
	long			count_selected() {
		long			n(0);
		for ( Bcomptype* c=this; c; c=c->next ) if ( sel ) n++;
		return n;
	}
	void		update(Bcomptype& ct) {
		id = ct.id;
		ind = ct.ind;
		fmod = ct.fmod;
		num = ct.num;
		cnt = ct.cnt;
		mas = ct.mas;
		chrg = ct.chrg;
		lrad = ct.lrad;
		rad = ct.rad;
		fom = ct.fom;
		sel = ct.sel;
		coefficients(ct.coefficients());
	}
	bool			check() {
		if ( id.empty() ) id = "?";
		if ( fmod.empty() ) fmod = "?";
		if ( mas < 1 ) mas = 1;
		return 1;
	}
	void			show() {
		cout << "Identity:                       " << id << endl;
		cout << "Index:                          " << ind << endl;
		cout << "File name:                      " << fmod << endl;
		cout << "Number:                         " << num << endl;
		cout << "Count:                          " << cnt << endl;
		cout << "Mass:                           " << mas << endl;
		cout << "Charge:                         " << chrg << endl;
		cout << "Link radius:                    " << lrad << endl;
		cout << "Radius:                         " << rad << endl;
		cout << "FOM:                            " << fom << endl;
		cout << "Select:                         " << sel << endl;
		cout << "Coefficients:" << endl;
		for ( int i=0; i<11; ++i ) cout << i+1 << tab << coef[i] << endl;
	}
} ;

/************************************************************************
@Object: class Blinktype
@Description:
	A structure for a link type.
@Features:
	This defines all the properties of a link type.
*************************************************************************/
class Blinktype {
private:
	void	initialize() {
		len = dist = 1;
		Klen = Kdist = 1;
		fom = 0;
		sel = 1;
	}
	string		id[2];		// Component identifiers
	float		len;		// Covalent bond length
	float		dist;		// Van der Waals interaction distance
	float		Klen;		// Covalent bond strenght
	float		Kdist;		// Van der Waals interaction strength
	float		fom;		// Figure-of-merit
	int			sel;		// Selection flag
public:
	Blinktype() { initialize(); }
	Blinktype(string& id1, string& id2) {
		initialize();
		id[0] = id1;
		id[1] = id2;
	}
	Blinktype(string& id1, string& id2, double d) {
		initialize();
		id[0] = id1;
		id[1] = id2;
		len = d;
	}
	void			identifier(string s, int i) { id[i] = s; }
	string&			identifier(int i) { return id[i]; }
	void			length(double d) { len = d; }
	double			length() { return len; }
	void			distance(double d) { dist = d; }
	double			distance() { return dist; }
	void			Klength(double d) { Klen = d; }
	double			Klength() { return Klen; }
	void			Kdistance(double d) { Kdist = d; }
	double			Kdistance() { return Kdist; }
	void			FOM(double d) { fom = d; }
	double			FOM() { return fom; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
} ;

/************************************************************************
@Object: class Bangletype
@Description:
	Model angle parameter structure.
@Features:
	List of components.
*************************************************************************/
class Bangletype {
private:
	void	initialize() {
		ang = 0;
		Kang = 1;
		fom = 0;
		sel = 1;
	}
	string			id[3];			// List of component identifiers
	float			ang;			// Reference angle
	float			Kang;			// Angle interaction strength
	float			fom;			// Figure-of-merit
	int				sel;			// Selection flag
public:
	Bangletype() { initialize(); }
	Bangletype(string& id1, string& id2, string& id3) {
		initialize();
		id[0] = id1;
		id[1] = id2;
		id[2] = id3;
	}
	void			identifier(string s, int i) { id[i] = s; }
	string&			identifier(int i) { return id[i]; }
	void			angle(double d) { ang = d; }
	double			angle() { return ang; }
	void			Kangle(double d) { Kang = d; }
	double			Kangle() { return Kang; }
	void			FOM(double d) { fom = d; }
	double			FOM() { return fom; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
} ;


/************************************************************************
@Object: class Bcomponent
@Description:
	Model component parameter structure.
@Features:
	Type which must be a valid component type. 
	Location.
	Orientation.
	Vector for use as force or displacement.
	Radius.
	Color.
	Figure-of-merit.
	Selection.
*************************************************************************/
class Bcomponent {
private:
	void	initialize() {
		next = NULL;
		typ = NULL;
		viw = View2<float>(0,0,1,0);
		rad = 1;
		rgba = RGBA<float>(1,1,1,1);
		den = 1;
		chrg = 0;
		fom = 0;
		sel = 1;
	}
public:
	Bcomponent*			next;			// Next component in list
private:
	string				id;				// Component identifier
	string				dsc;			// Additional description
	Bcomptype*			typ;			// Component type pointer
	Vector3<float>		loc;			// Location coordinates (angstroms)
	View2<float>		viw; 			// View: 3-value unit vector and angle (radians)
	Vector3<float>		vec;			// Vector for use as force or displacement
	Vector3<float>		vel;			// Vector for use as velocity
	float				rad;			// Display radius
	RGBA<float>			rgba;			// RGBA color
	float				den;			// Component density or intensity
	float				chrg;			// Component charge
	float				fom;			// Figure-of-merit
	int					sel;			// Selection flag
public:
	vector<Bcomponent*>	link;			// Links to other components
	vector<Bcomponent*>	neighbor;		// Neighboring components
	vector<int>			flag;			// Flags for connectors
	Bcomponent() { initialize(); }
	Bcomponent(string s) { initialize(); id = s; }
	Bcomponent(long i) { initialize(); id = to_string(i); }
	Bcomponent(Bcomponent* c) {
		next = NULL;
		id = c->id;
		typ = NULL;
		loc = c->loc;
		viw = c->viw;
		vec = c->vec;
		vel = c->vel;
		rad = c->rad;
		rgba = c->rgba;
		den = c->den;
		fom = c->fom;
		sel = c->sel;
	}
	void			identifier(string s) { id = s; }
	void			identifier(Bstring s) { id = s.str(); }
	string&			identifier() { return id; }
	void			description(string s) { dsc = s; }
	void			description(Bstring s) { dsc = s.str(); }
	string&			description() { return dsc; }
	void			type(Bcomptype* t) { typ = t; }
	Bcomptype*		type() { return typ; }
	void			location(Vector3<double> v) { loc = v; }
	Vector3<float>&	location() { return loc; }
	void			shift(Vector3<double> v) { loc += v; }
	void			scale(Vector3<double> v) { loc *= v; }
	void			scale(double d) { loc *= d; }
	void			view(View2<float> v) { viw = v; }
	void			view(Vector3<float> v) { viw = View2<float>(v); }
	View2<float>&	view() { return viw; }
	void			force(Vector3<double> v) { vec = v; }
	Vector3<float>&	force() { return vec; }
	void			velocity(Vector3<double> v) { vel = v; }
	Vector3<float>&	velocity() { return vel; }
	void			radius(double r) { rad = r; }
	double			radius() { return rad; }
	void			color(RGBA<float> c) { rgba = c; }
	RGBA<float>&	color() { return rgba; }
	void			density(double d) { den = d; }
	double			density() { return den; }
	void			charge(double d) { chrg = d; }
	double			charge() { return chrg; }
	void			FOM(double d) { fom = d; }
	double			FOM() { return fom; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
	long			select_increment() { return sel++; }
	Bcomponent*		add(string s) {
		Bcomponent* 	c(this);
		while ( c->next ) c = c->next;
		return c->next = new Bcomponent(s);
	}
	Bcomponent*		add(long i) {
		Bcomponent* 	c(this);
		while ( c->next ) c = c->next;
		return c->next = new Bcomponent(i);
	}
	Bcomponent*		add(Bcomponent* c) {
		Bcomponent* 	cn(this);
		while ( cn->next ) cn = cn->next;
		return cn->next = new Bcomponent(c);
	}
	void			add_link(Bcomponent* comp) {
		link.push_back(comp);
		flag.push_back(1);
	}
	Bcomponent*		find(string s) {
		Bcomponent* 	c(this);
		while ( c && c->id != s ) c = c->next;
		return c;
	}
	Bcomponent*		find_or_add(string s) {
		Bcomponent* 	c = NULL;
		if ( ( c = find(s) ) ) return c;
		else return add(s);
	}
	Bcomponent*		find_closest(Vector3<double> coor) {
		Bcomponent* 	c(this);
		Bcomponent* 	cc = NULL;
		double			d, dmin(1e30);
		for ( ; c; c = c->next ) if ( c->sel ) {
			d = c->loc.distance(coor);
			if ( dmin > d ) {
				dmin = d;
				cc = c;
			}
		}
		return cc;
	}
	bool			find_link_exists(Bcomponent* c) {
		for ( auto it = link.begin(); it != link.end(); ++it )
			if ( *it == c ) return 1;
		return 0;
	}
	void			find_and_add_links(string s1, string s2) {
		Bcomponent* 	c1 = find(s1);
		Bcomponent* 	c2 = find(s2);
		if ( c1 && c2 ) {
			c1->add_link(c2);
			c2->add_link(c1);
		}
	}
	long			count() {
		long			n(0);
		for ( Bcomponent* c=this; c; c=c->next ) n++;
		return n;
	}
	long			count_selected() {
		long			n(0);
		for ( Bcomponent* c=this; c; c=c->next ) if ( sel ) n++;
		return n;
	}
	Vector3<double>	normal() {
		long			i, j;
		Vector3<double>	v1, v2, n;
		for ( i=1; i<link.size(); ++i ) {
			v1 = link[i]->loc - loc;
			for ( j=0; j<i; ++j ) {
				v2 = link[j]->loc - loc;
				v2 = v2.cross(v1);
				if ( n.angle(v2) > M_PI_2 ) v2 = -v2;
				n += v2;
			}
		}
		if ( n.length2() ) n.normalize();
		return n;
	}
	void			calculate_normals() {
		for ( Bcomponent* c=this; c; c=c->next ) c->view(c->normal());
	}
	bool			check() {
		viw.normalize();
		if ( rad < 1e-6 ) rad = 1;
		if ( !typ ) {
			cerr << "Component " << id << " has no type assigned" << endl;
			return 0;
		}
		return 1;
	}
} ;

/************************************************************************
@Object: class Blink
@Description:
	Model link parameter structure.
@Features:
	Link between components.
	Link rotation angle for orientation.
	Length.
	Radius.
	Color.
	Figure-of-merit.
	Selection.
*************************************************************************/
class Blink {
private:
	void	initialize() {
		next = NULL;
		comp[0] = comp[1] = NULL;
		ang = 0;
		len = 0;
		rad = 1;
		rgba = RGBA<float>(1,1,1,1);
		fom = 0;
		sel = 1;
	}
	void	insert_components(Bcomponent* comp1, Bcomponent* comp2) {
		comp[0] = comp1;
		comp[1] = comp2;
		rgba = (comp1->color() + comp2->color())/2;
		if ( len < 1e-6 ) len = comp1->location().distance(comp2->location());
		if ( rad < 1e-6 ) rad = (comp1->radius() + comp2->radius())/4;
	}
public:
	Blink*			next;			// Next link in list
	Bcomponent*		comp[2];		// Linked components
private:
	float			ang;			// Rotation angle around the link
	float			len;			// Link reference length
	float			rad;			// Display radius
	RGBA<float>		rgba;			// RGBA color
	float			fom;			// Figure-of-merit
	int				sel;			// Selection flag
public:
	Blink() { initialize(); }
	Blink(Bcomponent* comp1, Bcomponent* comp2) {
		initialize();
		insert_components(comp1, comp2);
	}
	Blink(Bcomponent* comp1, Bcomponent* comp2, double d, double r) {
		initialize();
		len = d;
		rad = r;
		insert_components(comp1, comp2);
	}
	Blink(Blink* l) {
		initialize();
		ang = l->ang;
		len = l->len;
		rad = l->rad;
		rgba = l->rgba;
		fom = l->fom;
		sel = l->sel;
	}
	void			angle(double a) { ang = a; }
	double			angle() { return ang; }
	void			length(double d) { len = d; }
	double			length() { return len; }
	void			radius(double r) { rad = r; }
	double			radius() { return rad; }
	void			color(RGBA<float> c) { rgba = c; }
	RGBA<float>&	color() { return rgba; }
	void			adopt_component_color() {
		rgba = (comp[0]->color() + comp[1]->color())/2;
	}
	void			FOM(double d) { fom = d; }
	double			FOM() { return fom; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
	long			select_increment() { return sel++; }
	Blink*			add(Bcomponent* comp1, Bcomponent* comp2) {
		Blink* 			l(this);
		while ( l->next ) l = l->next;
		return l->next = new Blink(comp1, comp2);
	}
	Blink*			add(Bcomponent* comp1, Bcomponent* comp2, double d, double r) {
		Blink* 			l(this);
		while ( l->next ) l = l->next;
		return l->next = new Blink(comp1, comp2, d, r);
	}
	Blink*			add(Blink* l) {
		Blink* 			ln(this);
		while ( ln->next ) ln = ln->next;
		return ln->next = new Blink(l);
	}
	Blink*			find(Bcomponent* comp1, Bcomponent* comp2) {
		Blink* 			l(this);
		for ( l = this; l; l = l->next ) {
			if ( l->comp[0] == comp1 && l->comp[1] == comp2 ) return l;
			if ( l->comp[0] == comp2 && l->comp[1] == comp1 ) return l;
		}
//		cerr << "Error: Link not found!" << endl;
		return NULL;
	}
	Blink*			find(string s1, string s2) {
		Blink* 			l(this);
		for ( l = this; l; l = l->next ) {
			if ( l->comp[0]->identifier() == s1 && l->comp[1]->identifier() == s2 ) return l;
			if ( l->comp[0]->identifier() == s2 && l->comp[1]->identifier() == s1 ) return l;
		}
		return NULL;
	}
	long			index(Bcomponent* comp1, Bcomponent* comp2) {
		long			i(0);
		Blink* 			l(this);
		for ( l = this; l; l = l->next, ++i ) {
			if ( l->comp[0] == comp1 && l->comp[1] == comp2 ) return i;
			if ( l->comp[0] == comp2 && l->comp[1] == comp1 ) return i;
		}
		cerr << "Error: Link not found!" << endl;
		return -1;
	}
	long			count() {
		long			n(0);
		for ( Blink* l=this; l; l=l->next ) n++;
		return n;
	}
	bool			check() {
		if ( len < 1e-6 )
			len = comp[0]->location().distance(comp[1]->location());
		if ( rad < 1e-6 )
			rad = 0.1*len;
		return 1;
	}
} ;

/************************************************************************
@Object: class Bpolygon
@Description:
	Model polygon parameter structure.
@Features:
	List of components.
	Polygon normal.
	A flag to indicate if it is a closed polygon.
*************************************************************************/
class Bpolygon {
private:
	void	initialize() {
		next = NULL;
		clos = 0;
		fom = 0;
		sel = 1;
	}
public:
	Bpolygon*		next;
	vector<Bcomponent*>	comp;		// List of components
private:
	Vector3<float>	norm;			// Vector normal to the polygon plane
	int				clos;			// Flag to indicate a closed polygon
	float			fom;			// Figure-of-merit
	int				sel;			// Selection flag
public:
	Bpolygon() { initialize(); }
	Bpolygon(Bpolygon* p) {
		next = NULL;
		clos = p->clos;
		norm = p->norm;
		fom = p->fom;
		sel = p->sel;
	}
	long			size() { return comp.size(); }
	void			normal(Vector3<float> n) { norm = n; }
	Vector3<float>&	normal() { return norm; }
	void			closed(int i) { clos = i; }
	int				closed() { return clos; }
	void			FOM(double d) { fom = d; }
	double			FOM() { return fom; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
	long			select_increment() { return sel++; }
	Bpolygon*		add(Bpolygon* p) {
		Bpolygon* 			pn(this);
		while ( pn->next ) pn = pn->next;
		return pn->next = new Bpolygon(p);
	}
	long			count() {
		long			n(0);
		for ( Bpolygon* p=this; p; p=p->next ) n++;
		return n;
	}
} ;

/************************************************************************
@Object: class Bgroup
@Description:
	Componet group parameter structure.
@Features:
	Grouping of components into substructures for each model.
*************************************************************************/
class Bgroup {
public:
	Bgroup*			next;			// Next group in list
private:
	string			id;				// Group identifier
	string			gtp;			// Group type identifier
	string			tp;				// Type identifier
	string			dsc;			// Additional description
	string			sq;				// Sequence
public:
	Bgroup(string s) { next = NULL; id = s; }
	void			identifier(string s) { id = s; }
	string&			identifier() { return id; }
	void			group_type(string s) { gtp = s; }
	string&			group_type() { return gtp; }
	void			type(string s) { tp = s; }
	string&			type() { return tp; }
	void			description(string s) { dsc = s; }
	string&			description() { return dsc; }
	void			sequence(string s) { sq = s; }
	string&			sequence() { return sq; }
	Bgroup*			add(string s) {
		Bgroup* 		g(this);
		while ( g->next ) g = g->next;
		return g->next = new Bgroup(s);
	}
} ;

/************************************************************************
@Object: class Bmodel
@Description:
	Model parameter structure.
@Features:
	Orientation parameters for each model.
*************************************************************************/
class Bmodel {
private:
	void	initialize() {
		next = NULL;
		hand = 0;
		sym = "C1";
		fmap = "?";
		img_num = 0;
		fom = 0;
		sel = 1;
		type = NULL;
		comp = NULL;
		link = NULL;
		poly = NULL;
		group = NULL;
	}
public:
	Bmodel*			next;			// Next model in list
private:
	string			com;			// Model comment string
	string			id;				// Model identifier
	string			type_id;		// Type identifier
	string			fmap;			// Density map file name
	string			fmask;			// Multi-level mask file name
	int				img_num;		// Image number in map file
	string			sym;			// Symmetry label
	int				hand;			// Hand or enantiomorph
	float			fom;			// Figure-of-merit
	int				sel;			// Selection flag
	Vector3<double>	min;			// Coordinate minima
	Vector3<double>	max;			// Coordinate maxima
public:
	Bcomptype*		type;			// Component type list
	Bcomponent*		comp;			// Component list
	Blink*			link;			// Link list
	Bpolygon*		poly;			// Polygon list
	Bgroup*			group;			// Component group list
	Bmodel() { initialize(); }
	Bmodel(string s) { initialize(); id = s; }
	Bmodel(Bstring s) { initialize(); id = s.str(); }
	Bmodel(long i) { initialize(); id = to_string(i); }
	void			comment(string s) { com = s; }
	string&			comment() { return com; }
	void			identifier(string s) { id = s; }
	void			identifier(Bstring s) { id = s.str(); }
	string&			identifier() { return id; }
	void			model_type(string s) { type_id = s; }
	string&			model_type() { return type_id; }
	void			mapfile(string s) { fmap = s; }
	string&			mapfile() { return fmap; }
	void			maskfile(string s) { fmask = s; }
	string&			maskfile() { return fmask; }
	void			image_number(long i) { img_num = i; }
	long			image_number() { return img_num; }
	void			symmetry(string s) { sym = s; }
	string&			symmetry() { return sym; }
	void			handedness(long i) { hand = i; }
	long			handedness() { return hand; }
	void			FOM(double d) { fom = d; }
	double			FOM() { return fom; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
	long			select_increment() { return sel++; }
	long			select_all() {
		long		n(0);
		for ( Bcomponent* c = comp; c; c = c->next ) {
			c->select(1);
			n++;
		}
		return n;
	}
	void			deselect_all() {
		for ( Bcomponent* c = comp; c; c = c->next )
			c->select(0);
	}
	// Deselects outside bounds
	long			select_within_bounds(Vector3<double>& start, Vector3<double>& end) {
		long			nsel(0);
		for ( Bcomponent* c = comp; c; c = c->next ) if ( c->select() ) {
			if ( ( c->location() >= start ) && ( c->location() <= end ) ) nsel++;
			else c->select(0);
		}
		return nsel;
	}
	Vector3<double>	minimum() { return min; }
	Vector3<double>	maximum() { return max; }
	Bmodel*			add(string s) {
		Bmodel* 		m(this);
		while ( m->next ) m = m->next;
		return m->next = new Bmodel(s);
	}
	Bmodel*			add(long i) {
		Bmodel* 		m(this);
		while ( m->next ) m = m->next;
		return m->next = new Bmodel(i);
	}
	Bmodel*			add(Bmodel* model) {
		Bmodel* 		m(this);
		while ( m->next ) m = m->next;
		return m->next = model;
	}
	Bmodel*			find(string s) {
		Bmodel* 		m(this);
		while ( m && m->id != s ) m = m->next;
		return m;
	}
	Bmodel*			copy() {
		Bmodel*			m = new Bmodel(id);
		m->next = NULL;
		m->type_id = type_id;
		m->com = com;
		m->sym = sym;
		m->fmap = fmap;
		m->img_num = img_num;
		m->sel = sel;
		m->fom = fom;
		m->hand = hand;
		m->type = copy_types();
		m->comp = copy_components();
		m->link = copy_links(m->comp);
		m->poly = copy_polygons(m->comp);
		return m;
	}
	Bcomptype*		copy_types() {
		Bcomptype*		ct;
		Bcomptype*		ctn = NULL;
		Bcomptype*		ct1 = NULL;
		for ( ct = type; ct; ct = ct->next ) {
			if ( ct1 ) ct1 = ct1->add(ct);
			else ctn = ct1 = new Bcomptype(ct);
		}
		return ctn;
	}
	Bcomponent*		copy_components() {
		Bcomponent*		c;
		Bcomponent*		cn = NULL;
		Bcomponent*		c1 = NULL;
		for ( c = comp; c; c = c->next ) {
			if ( c1 ) c1 = c1->add(c);
			else cn = c1 = new Bcomponent(c);
			c1->type(add_type(c->type()->identifier()));
			for ( auto it = c->link.begin(); it != c->link.end(); ++it )
				cn->find_and_add_links(c->identifier(), (*it)->identifier());
		}
		return cn;
	}
	Blink*			copy_links(Bcomponent* c) {
		Blink*			l;
		Blink*			ln = NULL;
		Blink*			l1 = NULL;
		for ( l = link; l; l = l->next ) {
			if ( l1 ) l1 = l1->add(l);
			else ln = l1 = new Blink(l);
			l1->comp[0] = c->find(l->comp[0]->identifier());
			l1->comp[1] = c->find(l->comp[1]->identifier());
		}
		return ln;
	}
	Bpolygon*		copy_polygons(Bcomponent* c) {
		Bpolygon*		p;
		Bpolygon*		pn = NULL;
		Bpolygon*		p1 = NULL;
		for ( p = poly; p; p = p->next ) {
			if ( p1 ) p1 = p1->add(p);
			else pn = p1 = new Bpolygon(p);
			for ( auto it=p->comp.begin(); it!=p->comp.end(); ++it )
				pn->comp.push_back(c->find((*it)->identifier()));
		}
		return pn;
	}
	long			count() {
		long			n(0);
		for ( Bmodel* m=this; m; m=m->next ) n++;
		return n;
	}
	Bcomponent*		add_component(string s) {
		if ( comp ) return comp->add(s);
		return comp = new Bcomponent(s);
	}
	Bcomponent*		add_component(long i) {
		if ( comp ) return comp->add(i);
		return comp = new Bcomponent(i);
	}
	Bcomponent*		add_component(Bcomponent* c) {
		if ( comp ) return comp->add(c);
		return comp = new Bcomponent(c);
	}
	Bcomponent*		find_component(string s) {
		if ( comp ) return comp->find(s);
		return NULL;
	}
	Blink*			find_link(Bcomponent* c1, Bcomponent* c2) {
		if ( link ) return link->find(c1, c2);
		return NULL;
	}
	Blink*			find_link(string s1, string s2) {
		if ( link ) return link->find(s1, s2);
		return NULL;
	}
	float			component_count() { return comp->count(); }
	float			component_count_selected() { return comp->count_selected(); }
	long			set_component_radius(double r) {
		long			n(0);
		Bcomponent*		c;
		for ( c = comp; c; c = c->next ) if ( c->select() ) {
			c->radius(r);
			n++;
		}
		return n;
	}
	long			set_component_density(double d) {
		long			n(0);
		Bcomponent*		c;
		for ( c = comp; c; c = c->next ) if ( c->select() ) {
			c->density(d);
			n++;
		}
		return n;
	}
	long			set_component_fom(double f) {
		long			n(0);
		Bcomponent*		c;
		for ( c = comp; c; c = c->next ) if ( c->select() ) {
			c->FOM(f);
			n++;
		}
		return n;
	}
	Bcomponent*		closest_component(Vector3<double> loc) {
		return comp->find_closest(loc);
	}
	void			calculate_normals() { comp->calculate_normals(); }
	long			set_type_filenames(string& fn) {
		long			n(0);
		Bcomptype*		ct;
		for ( ct = type; ct; ct = ct->next ) if ( ct->select() ) {
			ct->file_name(fn);
			n++;
		}
		return n;
	}
	Bcomptype*		add_type(string s) {
		if ( find_type(s) ) return find_type(s);
		if ( type ) return type->add(s);
		return type = new Bcomptype(s);
	}
	Bcomptype*		add_type(Bstring s) {
		if ( find_type(s.str()) ) return find_type(s.str());
		if ( type ) return type->add(s.str());
		return type = new Bcomptype(s.str());
	}
	Bcomptype*		add_type(long i) {
		if ( type ) return type->add(i);
		return type = new Bcomptype(i);
	}
	Bcomptype*		add_type(string fn, long img_num) {
		if ( type ) return type->add(fn, img_num);
		return type = new Bcomptype(fn, img_num);
	}
	Bcomptype*		add_type(string s, string fn, long img_num) {
		if ( find_type(s) ) return find_type(s);
		if ( type ) return type->add(s, fn, img_num);
		return type = new Bcomptype(s, fn, img_num);
	}
	Bcomptype*		add_type(Bcomptype* ct) {
//		cout << "In add_type:" << endl;
//		ct->show();
		if ( find_type(ct->identifier()) ) return find_type(ct->identifier());
		if ( type ) return type->add(ct);
		return type = new Bcomptype(ct);
	}
	Bcomptype*		find_type(string s) {
		if ( type ) return type->find(s);
		return NULL;
	}
	long			component_type_count() { return type->count(); }
	long			component_type_count_selected() { return type->count_selected(); }
	void			update_type_counts() {
		Bcomptype*		ct;
		Bcomponent*		c;
		for ( ct = type; ct; ct = ct->next ) ct->component_count(0);
		for ( c = comp; c; c = c->next )
			if ( c->type() )
				c->type()->component_count_increment();
	}
	void		update_component_types(map<string,Bcomptype>& types) {
		for ( Bcomptype* ct = type; ct; ct = ct->next )
			if ( types.find(ct->identifier()) != types.end() )
				ct->update(types[ct->identifier()]);
	}
	Bgroup*		add_group(string s) {
		if ( group ) return group->add(s);
		return group = new Bgroup(s);
	}
	long			link_count() { return link->count(); }
	long			polygon_count() { return poly->count(); }
	bool			check() {
//		cout << "checking model " << id << endl;
		if ( com.empty() ) com = "?";
		if ( id.empty() ) id = "?";
		if ( type_id.empty() ) type_id = "?";
		if ( fmap.empty() ) fmap = "?";
		if ( !type ) add_type(1);
		for ( Bcomptype* ct = type; ct; ct = ct->next ) ct->check();
		for ( Bcomponent* c = comp; c; c = c->next ) {
			if ( !c->type() ) c->type(type);
			c->check();
		}
		for ( Blink* l = link; l; l = l->next ) l->check();
		update_type_counts();
		calculate_bounds();
//		cout << "finished checking" << endl;
		return 1;
	}
	void			calculate_bounds() {
		min = Vector3<double>(1e30,1e30,1e30);
		max = -min;
		Bcomponent*		c;
		for ( c = comp; c; c = c->next ) if ( c->select() ) {
			min = min.min(c->location());
			max = max.max(c->location());
		}
	}
	Vector3<double>	center_of_coordinates() {
		Vector3<double>	coc;
		if ( !sel ) return coc;
		long			n(0);
		Bcomponent*		c;
		for ( c = comp; c; c = c->next ) if ( c->select() ) {
			coc += comp->location();
			n++;
		}
		if ( n ) coc /= n;
		return coc;
	}
	Vector3<double>	center_of_mass() {
		Vector3<double>	com;
		if ( !sel ) return com;
		double			m(1), w(0);
		Bcomponent*		c;
		for ( c = comp; c; c = c->next ) if ( c->select() ) {
			if ( c->type() ) {
				m = c->type()->mass();
				if ( m <= 0 ) m = 1;
			}
			com += comp->location() * m;
			w += m;
		}
		if ( w ) com /= w;
		return com;
	}
	long			shift(Vector3<double> t) {
		for ( Bcomponent* c = comp; c; c = c->next ) c->shift(t);
		return component_count();
	}
	long			trim(Vector3<double> b) {
		Vector3<double>		start;
		Bcomponent* 		c;
		Bcomponent*			nu_list = NULL;
		Bcomponent*			c2 = NULL;
		for ( c = comp; c; c = c->next ) {
			if ( c->location().within(start, b) ) {
				if ( c2 ) c2 = c2->next = new Bcomponent(c);
				else nu_list = c2 = new Bcomponent(c);
				c2->type(add_type(c->type()->identifier()));
				for ( auto it = c->link.begin(); it != c->link.end(); ++it )
					nu_list->find_and_add_links(c->identifier(), (*it)->identifier());
			}
		}
		for ( c = comp; c; c = c2 ) {
			c2 = c->next;
			delete c;
		}
		comp = nu_list;
		return component_count();
	}
	void			show_components() {
		cout << "Components:" << endl;
		for ( Bcomponent* c = comp; c; c = c->next )
			cout << c->identifier() << tab << c->location() << endl;
	}
	void			show_links() {
		cout << "Links:" << endl;
		for ( Blink* l = link; l; l = l->next )
			cout << l->comp[0]->identifier() << tab << l->comp[1]->identifier() << endl;
	}
	vector<Bcomponent*>	component_array() {
		vector<Bcomponent*>	carr;
		for( Bcomponent* c = comp; c; c = c->next ) if ( c->select() )
			carr.push_back(comp);
		return carr;
	}
} ;

#define _Bmodel_
#endif

/* Function prototypes */
Bmodel*		read_model(Bstring& filename);
Bmodel*		read_model(Bstring* file_list);
Bmodel*		read_model(Bstring& filename, Bstring& paramfile);
Bmodel*		read_model(Bstring* file_list, Bstring& paramfile);
int			write_model(Bstring& filename, Bmodel* model);
int			write_model(Bstring& filename, Bmodel* model, int split);
Bmodel*		model_add(Bmodel** model, string id);
Bcomponent*	component_add(Bcomponent** comp, string& id);
Bcomponent*	component_add(Bcomponent** comp, unsigned long number);
Blink*		link_add(Blink** link, Bcomponent* comp1, Bcomponent* comp2, double length, double radius);
Blink*		link_add(Blink** link, Bcomponent* comp1, Bcomponent* comp2);
int			model_set_map_filenames(Bmodel* model, Bstring& mapfile);
int			model_set_type(Bmodel* model, Bstring& set_type);
int			model_change_type(Bmodel* model, Bstring& change_type);
int			model_check(Bmodel* model, Bstring path);
Bmodel*		model_list_copy(Bmodel* model);
int			component_list_kill(Bcomponent* comp);
int			comp_type_list_kill(Bcomptype* type);
int			model_link_list_kill(Bmodel* model);
int			link_kill(Blink** link_list, Bcomponent* comp, int i);
int			link_kill(Blink** link_list, Bcomponent* comp, Bcomponent* comp2);
int			poly_list_kill(Bpolygon* poly);
int			comp_associated_links_kill(Bcomponent* comp, Blink** link);
int 		model_kill(Bmodel* model);
int			model_associate(Bmodel* model, Bstring& associate_type, Bstring& associate_file);
int			model_associate_mass(Bmodel* model, Bstring& associate_type, double mass);
int			model_set_comptype_filenames(Bmodel* model, Bstring& filename);
long		model_set_component_radius(Bmodel* model, double comprad);


