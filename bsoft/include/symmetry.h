/**
@file	symmetry.h
@brief	Header file for general symmetry functions
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20210116
**/

#include "View.h"
#include "Bstring.h"

#ifndef _Bsymmetry_
/************************************************************************
@Object: class Bsymop
@Description:
	General symmetry operator.
@Features:
	Both point groups and space groups are covered.
	A symmetry operator describes the complete series of rotations around 
		an axis, with the number of operations given by the order field
		and the rotation angle given by the angle field. The shift field
		is used in gliding operations, such as needed for helical symmetry.
*************************************************************************/
class Bsymop {
private:
	Vector3<double>	ax;		// Symmetry axis
	int				n;		// Number of times the operation is applied, 1=identity
	double			ang;	// Rotation angle, derived from order when 0
	double			sh;		// Shift up symmetry axis for helix
	void	initialize() {
		n = 0;
		ang = 0;
		sh = 0;
	}
public:
	Bsymop() { initialize(); }
	Bsymop(double x, double y, double z, int nn, double a) :
		ax(Vector3<double>(x,y,z)), n(nn), ang(a), sh(0) {}
	void			order(int i) { n = i; }
	int				order() { return n; }
	void			angle(double a) { ang = a; }
	double			angle() { return ang; }
	void			shift(double s) { sh = s; }
	double			shift() { return sh; }
	void			axis(Vector3<double> a) { ax = a; }
	void			axis(double x, double y, double z) { ax = Vector3<double>(x,y,z); }
	Vector3<double>&	axis() { return ax; }
	double			normalize() { return ax.normalize(); }
	Matrix3			matrix() {
		return Matrix3(ax, M_PI*2.0/n);
	}
} ;

/************************************************************************
@Object: class Bsymmetry
@Description:
	General symmetry structure.
@Features:
	Both point groups and space groups are covered.
*************************************************************************/
struct Bsymmetry {
private:
	int 			pnt;		// Point group (< 1 if a crystal)
	int 			sp;			// Space group (< 1 if not a crystal)
	Bstring			lbl;		// Symmetry label string
	vector<Bsymop>	op; 		// Symmetry operators
private:
	void	initialize() {
		pnt = 0;
		sp = 0;
		lbl = "C1";
	}
	Bstring		clean_symstring(Bstring& sym) {
		int					j(0);
		lbl = "C1";
		if ( sym.length() < 1 ) return lbl;
		// Remove leading blanks and convert to upper case
		lbl = sym.no_space();
		// Alternate nomenclature for point groups
		if ( isdigit(lbl[0]) ) {
			j = lbl.integer();
			if ( lbl.contains("532") && lbl.contains("90") )
				lbl = "I90";
			else if ( lbl.contains("532") )
				lbl = "I";
			else if ( lbl.contains("432") )
				lbl = "O";
			else if ( lbl.contains("23") )
				lbl = "T";
			else if ( j > 1 && lbl[1] == 2 )
				lbl = Bstring(j, "D%d");
			else if ( j > 0 )
				lbl = Bstring(j, "C%d");
			else
				lbl = "C1";
		}
		lbl[0] = toupper(lbl[0]);
		return lbl;
	}
public:
	Bsymmetry() { initialize(); }
	Bsymmetry(Bstring sym);
//	~Bsymmetry() { op.clear(); }
	void			point(int i) { pnt = i; }
	int				point() { return pnt; }
	int				point() const { return pnt; }
	void			space(int i) { sp = i; }
	int				space() { return sp; }
	void			label(Bstring& s) { lbl = s; }
	Bstring&		label() { return lbl; }
	int		order() {
		int		ns(1);
		for ( auto it = op.begin(); it != op.end(); ++it ) ns *= it->order();
		return ns;
	}
	int			operations() { return op.size(); }
	Bsymop&		operator[](int i) { return op[i]; }
	vector<Matrix3>	matrices() {
		vector<Matrix3>		mat;
		mat.push_back(Matrix3(1));
		for ( auto it = op.begin(); it != op.end(); ++it )
			for ( int i = 1; i < it->order(); ++i )
				mat.push_back(Matrix3(it->axis(), i*M_PI*2.0L/it->order()));
		return mat;
	}
	void		transform(Matrix3& mat) {
		for ( auto it = op.begin(); it != op.end(); ++it )
			it->axis(mat * it->axis());
	}
} ;
#define _Bsymmetry_
#endif

// Function prototypes
string		symmetry_helical_label(double helix_rise, double helix_angle,
				int dyad_axis, int cyclic, double seam_shift);
View		view_symmetry_reference(Bsymmetry& sym);
Matrix3		symmetry_rotate_to_axis(Bsymmetry& sym, long axis, long axis_flag);
vector<Vector3<double>>	symmetry_get_axes(Bsymmetry& sym);
View*		symmetry_get_all_views(Bsymmetry& sym, View asu_view);
View*		symmetry_get_all_views(Bsymmetry& sym, View* views);
//Matrix3*	symmetry_get_all_matrices(Bsymmetry& sym, long& nsym);
vector<Matrix3>	symmetry_get_all_matrices(Bsymmetry& sym);
View* 		asymmetric_unit_views(Bsymmetry& sym, double theta_step, double phi_step, int flag);
View* 		asymmetric_unit_views(Bsymmetry& sym, double theta_step, double phi_step, double alpha_step, int flag);
View*		side_views(Bsymmetry& sym, double side_ang, double theta_step, double phi_step);
View*		side_views(Bsymmetry& sym, double side_ang, double theta_step, double phi_step, double alpha_step);
int 		change_views_to_asymmetric_unit(Bsymmetry& sym, View* view);
View 		find_asymmetric_unit_view(Bsymmetry& sym, View theview);
View 		find_closest_symmetric_view(Bsymmetry& sym, View view_ref, View view);
View* 		reference_asymmetric_unit_views(Bsymmetry& sym);
View 		random_symmetric_view(View& asu_view, Bsymmetry& sym);
int 		test_asymmetric_unit_view(View theview, Bsymmetry& sym);
int			sym_show_matrices(Bsymmetry& sym);
int			sym_show_operational_matrices(Bsymmetry& sym);
int			sym_show_pdb_matrices(Bsymmetry& sym);

