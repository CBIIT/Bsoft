/**
@file	options.h 
@brief	Header file for option handlers 
@author Bernard Heymann 
@date	Created: 20010613
@date	Modified: 20220718
**/

#include "Bstring.h"
#include "Vector3.h"
#include "View.h"
#include "View2.h"
#include "Euler.h"
#include "symmetry.h"
#include "UnitCell.h"
#include "Complex.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	A class for a command line input option.

	Each option is defined by a tag (i.e., the string following a hyphen)
	and a value (the next string after the tag).
**/
class Boption {
public:
	Boption*	next;
	Bstring		tag;
	Bstring		value;
private:
	int			err;
public:
	Boption(char* arg, char* arg1, const char* use[]) : next(NULL), err(0) {
//		cout << arg << tab << arg1 << endl;
		Bstring		usetag;
		for ( int n=0; use[n] != NULL; n++ ) {
			usetag = use[n];
//			cout << usetag << endl;
			if ( usetag[0] == '-' ) {
				if ( usetag.find(arg, 0) == 0 ) {
					if ( tag.length() ) {
						cerr << "Error: Ambiguous option " << arg << endl;
						err++;
					} else {
						tag = usetag.pre(' ').post('-');
						if ( usetag.post(' ').pre(' ').length() ) value = arg1;
						if ( verbose & VERB_DEBUG )
							cout << "DEBUG Boption: usetag=" << tag << " value=" << value << endl;
					}
				}
			}
		}
		if ( tag.length() < 1 ) {
			cerr << "Error: Option " << arg << " is not defined!" << endl;
			err++;
		}
	}
/*	int 		verbosity(char* optarg) {
		Bstring		arg(optarg);
		return verbosity(arg);
	}
	int 		verbosity(Bstring& optarg) {
		verbose = optarg.integer();
		Bstring			arg = optarg.lower();
		if ( arg.contains("non") ) verbose |= VERB_NONE;
		if ( arg.contains("res") ) verbose |= VERB_RESULT;
		if ( arg.contains("lab") ) verbose |= VERB_LABEL;
		if ( arg.contains("pro") ) verbose |= VERB_PROCESS;
		if ( arg.contains("sta") ) verbose |= VERB_STATS;
		if ( arg.contains("ful") ) verbose |= VERB_FULL;
		if ( arg.contains("tim") ) verbose |= VERB_TIME;
		if ( arg.contains("mem") ) verbose |= VERB_MEMORY;
		if ( arg.contains("deb") ) verbose |= VERB_DEBUG;
		if ( arg.contains("star") ) verbose |= VERB_DEBUG_STAR;
		if ( arg.contains("dm") ) verbose |= VERB_DEBUG_DM;
 		return verbose;
	}*/
	int			errors() {
		int				n(0);
		Boption*		opt;
		for ( opt = this; opt; opt = opt->next ) n += opt->err;
		return n;
	}
	Bstring		filename() {
		if ( value.length() < 1 ) {
			cerr << "-" << tag << ": A valid file name must be specified!" << endl;
			err++;
		}
		if ( value.contains(",") ) return value.pre(',');
		else return value;
	}
	DataType 	datatype() {
		if ( value.length() < 1 ) {
      	 	cerr << "-datatype: A data type must be specified!" << endl;
			err++;
		}
		return getdatatype(value[0]);
	}
	long		integer() {
		if ( value.length() < 1 ) {
			cerr << "-" << tag << ": An integer must be specified!" << endl;
			err++;
		}
		return value.integer();
	}
	double		real() {
		if ( value.length() < 1 ) {
			cerr << "-" << tag << ": A real value must be specified!" << endl;
			err++;
		}
		return value.real();
	}
	double		angle() {	// Returns an angle in radians
		if ( value.length() < 1 ) {
			cerr << "-" << tag << ": An angle must be specified!" << endl;
			err++;
		}
		if ( value[0] == 'p' || value[0] == 'P' ) {
			return M_PI;
		} else if ( value.contains("°") ) {
			return value.real() * M_PI/180.0;
		} else if ( value[0] == 'p' ) {
			double			d(1);
			if ( value.contains("/") ) d = 1.0/value.post('/').real();
			else if ( value.contains("*") ) d = value.post('*').real();
			return d*M_PI;
		}
		return value.real();
	}
	double		angle_degrees() {	// Converts from degrees and returns an angle in radians
		if ( value.length() < 1 ) {
			cerr << "-" << tag << ": An angle must be specified!" << endl;
			err++;
		}
		return value.real() * M_PI/180.0;
	}
	template <typename T1, typename T2>
	long	values(T1& v1, T2& v2) {
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() > 0 ) v1 = (T1)d[0];
		if ( d.size() > 1 ) v2 = (T2)d[1];
		return d.size();
	}
	template <typename T1, typename T2, typename T3>
	long	values(T1& v1, T2& v2, T3& v3) {
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() > 0 ) v1 = (T1)d[0];
		if ( d.size() > 1 ) v2 = (T2)d[1];
		if ( d.size() > 2 ) v3 = (T3)d[2];
		return d.size();
	}
	template <typename T1, typename T2, typename T3, typename T4>
	long	values(T1& v1, T2& v2, T3& v3, T4& v4) {
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() > 0 ) v1 = (T1)d[0];
		if ( d.size() > 1 ) v2 = (T2)d[1];
		if ( d.size() > 2 ) v3 = (T3)d[2];
		if ( d.size() > 3 ) v4 = (T4)d[3];
		return d.size();
	}
	template <typename T1, typename T2, typename T3, typename T4, typename T5>
	long	values(T1& v1, T2& v2, T3& v3, T4& v4, T5& v5) {
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() > 0 ) v1 = (T1)d[0];
		if ( d.size() > 1 ) v2 = (T2)d[1];
		if ( d.size() > 2 ) v3 = (T3)d[2];
		if ( d.size() > 3 ) v4 = (T4)d[3];
		if ( d.size() > 4 ) v5 = (T5)d[4];
		return d.size();
	}
	double			fill(int& fill_type) {
		double			f = value.real();
		fill_type = FILL_USER;
		Bstring			arg = value.lower();
		if ( arg[0] == 'a' ) fill_type = FILL_AVERAGE;
		if ( arg[0] == 'b' ) fill_type = FILL_BACKGROUND;
		if ( arg.contains("mi") ) fill_type = FILL_MIN;
		if ( arg.contains("ma") ) fill_type = FILL_MAX;
		return f;
	}
	Vector3<double>	vector3() {
		Vector3<double>	vec;
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() < 2 ) {
			cerr << "-" << tag << ": At least 2 values for the vector must be specified!" << endl;
			err++;
		}
		for ( unsigned long i=0; i<d.size(); i++ ) vec[i] = d[i];
		return vec;
	}
	Vector3<double>	normalized_vector3() {
		Vector3<double>	vec = vector3();
		vec.normalize();
		return vec;
	}
	Vector3<long>	size() {
		Vector3<long>	v(1,1,1);
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() < 1 ) {
			cerr << "-" << tag << ": At least one size value must be specified!" << endl;
			err++;
		}
		for ( unsigned long i=0; i<d.size(); i++ ) v[i] = (long) (d[i] + 0.5);
		if ( d.size() < 2 ) v[2] = v[1] = v[0];
		else if ( d.size() < 3 ) v[2] = 1;
		return v;
	}
	Vector3<long>	size2() {
		Vector3<long>	v = size();
		v[2] = 1;
		return v;
	}
	Vector3<double>	origin() {
		Vector3<double>	ori;
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() < 1 ) {
			cerr << "-" << tag << ": At least one origin value must be specified!" << endl;
			err++;
		}
		for ( unsigned long i=0; i<d.size(); i++ ) ori[i] = d[i];
		return ori;
	}
	Vector3<double>	scale() {
		Vector3<double>	sc(1,1,1);
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() < 1 ) {
			cerr << "-" << tag << ": At least one scale value must be specified!" << endl;
			err++;
		}
		for ( unsigned long i=0; i<d.size() && i<3; i++ ) sc[i] = d[i];
		if ( d.size() < 2 ) sc[2] = sc[1] = sc[0];
		return sc;
	}
	template <typename T1, typename T2>
	long		line(Vector3<T1>& st, Vector3<T2>& en) {
		long 			i, n(0);
		vector<double>	d = value.split_into_doubles(",");
		n = d.size();
		if ( n < 6 ) {
			cerr << "-" << tag << ": Both start and end coordinates must be specified!" << endl;
			err++;
		}
		for ( i=0; i<n && i<3; i++ ) st[i] = (T1) d[i];
		for ( i=3; i<n && i<6; i++ ) en[i-3] = (T2) d[i];
		return n;
	}
	template <typename T1, typename T2>
	long		box(Vector3<T1>& st, Vector3<T2>& sz) {
		long 			i, n(0);
		vector<double>	d = value.split_into_doubles(",");
		n = d.size();
		if ( n < 4 ) {
			cerr << "-" << tag << ": At least one size value must be specified!" << endl;
			err++;
		}
		for ( i=0; i<n && i<3; i++ ) st[i] = (T1) d[i];
		for ( i=3; i<n && i<6; i++ ) sz[i-3] = (T2) d[i];
		if ( n < 4 ) sz[2] = sz[1] = sz[0];
		else if ( n < 5 ) sz[2] = 1;
		if ( sz.volume() < 1 ) sz = sz.max(1);
		return n;
	}
	template <typename T1, typename T2>
	long		box(vector<T1>& st, vector<T2>& sz) {
		long 			i, n(0);
		vector<double>	d = value.split_into_doubles(",");
		n = d.size();
		if ( n < 4 ) {
			cerr << "-" << tag << ": At least one size value must be specified!" << endl;
			err++;
		}
		for ( i=0; i<n && i<3; i++ ) st[i] = (T1) d[i];
		for ( i=3; i<n && i<6; i++ ) sz[i-3] = (T2) d[i];
		if ( n < 4 ) sz[2] = sz[1] = sz[0];
		else if ( n < 5 ) sz[2] = 1;
		if ( volume(sz) < 1 ) sz = vecmax(sz, 1);
		return n;
	}
	Bstring		symmetry_string() {
		if ( value.length() < 1 ) {
			cerr << "-symmetry: A symmetry must be specified!" << endl;
			err++;
		}
		Bsymmetry	sym(value);
		return sym.label();
	}
	Bsymmetry	symmetry() {
		if ( value.length() < 1 ) {
			cerr << "-symmetry: A symmetry must be specified!" << endl;
			err++;
		}
		return Bsymmetry(value);
	}
	View		view() {
		View			v;
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() < 3 ) {
			cerr << "-View: At least 3 values for the vector must be specified!" << endl;
			err++;
		}
		for ( unsigned long i=0; i<d.size(); i++ ) v[i] = d[i];
		v[3] *= M_PI/180;
		v.normalize();
		return v;
	}
	View2<double>	view2() {
		View2<double>	v;
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() < 3 ) {
			cerr << "-View: At least 3 values for the vector must be specified!" << endl;
			err++;
		}
		for ( unsigned long i=0; i<d.size(); i++ ) v[i] = d[i];
		v[3] *= M_PI/180;
		v.normalize();
		return v;
	}
	Euler		euler() {
		Euler			e;
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() < 3 ) {
			cerr << "-Euler: All three angles must be specified!" << endl;
			err++;
		}
		for ( long i=0; i<3; i++ ) e[i] = d[i]*M_PI/180.0;
		return e;
	}
 	double			real_unit(Bstring& s) {
		double			m = s.real();
		// Convert distances to angstrom
		if ( s.contains("m") )	// millimeter
			m *= 1e7;
		else if ( s.contains("u") )	// micrometer
			m *= 1e4;
		else if ( s.contains("n") )	// nanometer
			m *= 10;
		// Convert masses to Dalton
		else if ( s.contains("k") || value.contains("K") )	// kilo
			m *= 1e3;
		else if ( s.contains("M") )	// Mega
			m *= 1e6;
		else if ( s.contains("g") || value.contains("G") )	// Giga
			m *= 1e9;
		// Convert degrees to radians
		else if ( s.contains("d") )
			m *= M_PI/180.0;
		return m;
	}
	double			real_units() {
		return real_unit(value);
	}
/* 	double			real_units() {
		double			m = value.real();
		// Convert distances to angstrom
		if ( value.contains("m") )	// millimeter
			m *= 1e7;
		else if ( value.contains("u") )	// micrometer
			m *= 1e4;
		else if ( value.contains("n") )	// nanometer
			m *= 10;
		// Convert masses to Dalton
		else if ( value.contains("k") || value.contains("K") )	// kilo
			m *= 1e3;
		else if ( value.contains("M") )	// Mega
			m *= 1e6;
		else if ( value.contains("g") || value.contains("G") )	// Giga
			m *= 1e9;
		// Convert degrees to radians
		else if ( value.contains("d") )
			m *= M_PI/180.0;
		return m;
	}*/
	template <typename T1, typename T2>
 	long			real_units(T1& v1, T2& v2) {
 		Bstring*	s = value.split(",");
 		v1 = real_unit(*s);
 		s = s->next;
 		if ( s ) v2 = real_unit(*s);
 		else return 1;
 		return 2;
 	}
	template <typename T1, typename T2, typename T3>
 	long			real_units(T1& v1, T2& v2, T3& v3) {
 		Bstring*	s = value.split(",");
 		v1 = real_unit(*s);
 		s = s->next;
 		if ( s ) v2 = real_unit(*s);
 		else return 1;
 		s = s->next;
 		if ( s ) v3 = real_unit(*s);
 		else return 2;
 		return 3;
 	}
	UnitCell		unit_cell() {
		UnitCell		uc;
		vector<double>	d = value.split_into_doubles(",");
		if ( d.size() < 6 ) {
			cerr << "-unitcell: All 6 values for the unit cell must be specified!" << endl;
			err++;
		}
		for ( unsigned long i=0; i<d.size(); i++ ) uc[i] = d[i];
		uc.degrees_to_radians();
		uc.set_angle_range();
		return uc;
	}
	ComplexConversion complex_conversion() {
		ComplexConversion	conv(NoConversion);
		if ( value[0] == 'r' ) conv = Real;
		else if ( value[0] == 'i' ) {
			conv = Imaginary;
			if ( value[1] == 'n' ) conv = Intensity;
		} else if ( value[0] == 'a' || value[0] == 'A' ) conv = Amplitude;
		else if ( value[0] == 'I' ) conv = Intensity;
		return conv;
	}
	int				ctf_action() {
		int			action(0);
		if ( value.length() < 1 )
			cerr << "-action: An action must be specified!" << endl;
		else {
			action = value.integer();
			value = value.lower();
			if ( value.contains("fli") ) action = 1;
			if ( value.contains("app") ) action = 2;
			if ( value.contains("cor") ) action = 3;
			if ( value.contains("wie") ) action = 4;
			if ( value.contains("bas") ) action = 5;
			if ( value.contains("baseline2") ) action = 6;
			if ( value.contains("basef") ) action = 7;
			if ( value.contains("basec") ) action = 8;
			if ( value.contains("env") ) action = 9;
			if ( value.contains("app") &&
				value.contains("env") ) action = 10;
			if ( value.contains("prep") ) action = 11;
			if ( value.contains("fit") ) action = 12;
			if ( value.contains("prep") &&
				value.contains("fit") ) action = 13;
		}
		return action;
	}
} ;
	 
// Function prototypes 
Boption*	get_option_list(const char* use[], int argc, char* argv[], int& optind);
int			option_kill(Boption* option);
int 		get_option_verbose(char* optarg);
int 		get_option_verbose(Bstring& optarg);
//double		get_option_mass(Bstring& optarg);

