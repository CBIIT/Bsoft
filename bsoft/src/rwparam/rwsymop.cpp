/**
@file	rwsymop.cpp
@brief	Library routines to read and write symmetry operators
@author Bernard Heymann
@date	Created: 19991225
@date	Modified: 20210328
**/

#include "rwsymop.h"
#include "star.h"
#include "sym_tags.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
char* 		read_symop_star(Bstring& filename, int spacegroup, int& nsym);
char* 		read_symop_lib(Bstring& filename, int spacegroup, int& nsym);
int 		write_symop_star(Bstring& filename, int spacegroup, int nsym, char* symop, int line_len);
int 		write_pointgroup_star(Bstring& filename, Bsymmetry& sym, View ref_view);
float*		sym_matrices_from_text_list(int nsym, char* symop, int line_len);

/**
@brief 	Reading crystallographic symmetry operators.
@param 	&filename		file name.
@param 	spacegroup		crystal space group number.
@param 	&nsym			number of symmetry operators.
@return float* 			set of 12-value symmetry matrices.

	The symmetry operators are encoded as a set of matrices.

**/
float* 		read_symat(Bstring& filename, int spacegroup, int& nsym)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symat: Getting symmetry matrices from: " << filename << endl;
	
	char* 		symop = read_symop(filename, spacegroup, nsym);
	
	// Set up the symmetry operator matrix
	float*		mat = sym_matrices_from_text_list(nsym, symop, 80);
	
	delete[] symop;
	
	return mat;
}

/**
@brief 	Reading crystallographic symmetry operators.
@param 	&symopfile		file name.
@param 	spacegroup		crystal space group number.
@param 	&nsym			number of symmetry operators.
@return char* 			set of 12-value symmetry matrices.

	The symmetry operators are encoded as 80 character lines.

**/
char* 		read_symop(Bstring& symopfile, int spacegroup, int& nsym)
{
	// No space group has more than 192 operators
	if ( spacegroup < 1 || spacegroup > 100000 ) return 0;
	
	if ( symopfile.length() < 1 ) {
		symopfile = "symop.star";
		symopfile = parameter_file_path(symopfile);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symop: Reading symmetry operator file " << symopfile << endl;

	Bstring		ext = symopfile.extension();
	
	char*		symop = NULL;
	
	if ( spacegroup > 1 ) {
		if ( ext.contains("star") || ext.contains("cif") ) {
			symop = read_symop_star(symopfile, spacegroup, nsym);
		} else {
			symop = read_symop_lib(symopfile, spacegroup, nsym);
		}
		if ( !symop )
			cerr << "Error: No symmetry operator file read!" << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symop: Symmetry operator file read: " << symopfile << endl;
	
	return symop;
}

/**
@brief 	Writing crystallographic symmetry operators.
@param 	&filename		file name.
@param 	spacegroup		crystal space group number.
@return int				error code (<0 means failure).
**/
int 		write_symat(Bstring& filename, int spacegroup)
{
	int				nsym(0), err(0);
	Bstring			temp;
	char*			symop = read_symop(temp, spacegroup, nsym);
	if ( !symop || !nsym ) return -1;
	
	Bstring		symatfile = filename;
	
	err = write_symop_star(symatfile, spacegroup, nsym, symop, 80);
	
	delete[] symop;
	
	return err;
}

/**
@brief 	Writing point group symmetry operators.
@param 	&filename			file name.
@param 	&symmetry_string	symmetry string.
@param 	ref_view			reference view.
@return int					error code (<0 means failure).
**/
int 		write_pointgroup(Bstring& filename, Bstring& symmetry_string, View ref_view)
{
	Bsymmetry	sym(symmetry_string);
	
	int			err = write_pointgroup(filename, sym, ref_view);

	return err;
}

int 		write_pointgroup(Bstring& filename, Bsymmetry& sym, View ref_view)
{
	int			err(0);

	if ( filename.contains(".star") )
		err = write_pointgroup_star(filename, sym, ref_view);
	else {
		cerr << "Error: File type for " << filename << " not supported!" << endl;
		err = -1;
	}

	return err;
}

// Find space group label and operators in a STAR format file
char* 		read_symop_star(Bstring& filename, int spacegroup, int& nsym)
{
 	Bstar2					star;
	
 	if ( star.read(filename.str()) < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return NULL;
	}

	int				i(0), j;
	char*			symop = NULL;

	for ( auto ib: star.blocks() ) {
		if ( spacegroup == ib.integer(SYMMETRY_NUMBER) ) {
			for ( auto il: ib.loops() ) {
				if ( ( j = il.find(SYMMETRY_EQUIVXYZ) ) >= 0 ) {
					nsym = il.data().size();
					symop = new char[nsym*80];
					for ( auto ir: il.data() ) {
						string	symstr(ir[j]);
						strcpy(symop+80*i, symstr.c_str());
						i++;
					}
				}
			}
			break;
		}
	}
	
	return symop;
}

// Find space group label and operators in "symop.lib" or similar file
char* 		read_symop_lib(Bstring& filename, int spacegroup, int& nsym)
{
	int				i, j, k(0), notfound(1), number(0), nlines;
	char			aline[80], symall[4000];
	for ( i=0; i<4000; i++ ) symall[i] = 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symop_lib: Symmetry operator library: " << filename << endl << endl;

	Bstring			atfile;
	Bstring			symopfile;
	if ( filename.empty() ) symopfile = filename;
	else symopfile = "symop.lib";
	
	if ( access(symopfile.c_str(), R_OK) != 0 )
		symopfile = parameter_file_path(symopfile);

	ifstream		fsym;
	fsym.open(symopfile.c_str());
	if ( fsym.fail() ) return NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symop_lib: Symmetry operator library: " << symopfile << endl << endl;
	
	while ( notfound && fsym.getline(aline, 80) ) {
		sscanf( aline, "%d %d", &number, &nlines );
		if ( number == spacegroup ) notfound = 0;
	}
	
	if ( notfound ) {
		fsym.close();
		return NULL;
	}
	
	for ( i=nsym=0; i<nlines; i++ ) {
		fsym.getline(aline, 80);
		for ( j=0; j<(int)strlen(aline); j++ ) {
			if ( aline[j] != ' ' && aline[j] != '\n' ) {
				symall[k] = tolower(aline[j]);
				k++;
			}
			if ( aline[j] == '*' ) nsym++;
		}
		symall[k] = '*';
		k++;
		nsym++;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symop_lib: Spacegroup=" << number << " Noperators=" << nsym << endl;

	char*		symop = new char[nsym*80];
	for ( i=0; i<nsym*80; i++ ) symop[i] = 0;

	j = 0;
	for ( i=0; i<nsym; i++ ) {
		k = 0;
		while( symall[j] != '*' && k<80 ) {
			symop[i*80+k] = symall[j];
			j++;
			k++;
		}
		j++;
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_symop_lib: " << &symop[80*i] << endl;
	}
	
	fsym.close();
	
	return symop;
}

int 		write_symop_star(Bstring& filename, int spacegroup, int nsym, char* symop, int line_len)
{
	int				i;
	Bstring			s;
	Bstar2			star;

	star.comment("# Symmetry operators\n\n");

	BstarBlock&		block = star.add_block(to_string(spacegroup));

	block[SYMMETRY_NUMBER] = to_string(spacegroup);

	BstarLoop&		loop = block.add_loop();
	loop.tags()[SYMMETRY_EQUIVID] = 0;
	loop.tags()[SYMMETRY_EQUIVXYZ] = 1;

	for ( i=0; i<nsym; ++i ) {
		s = Bstring(symop+line_len*i);
		vector<string>&	vs = loop.add_row(2);
		vs[0] = to_string(i+1);
		vs[1] = s.str();
	}
	
	return star.write(filename.str());
}

int 		write_pointgroup_star(Bstring& filename, Bsymmetry& sym, View ref_view)
{
 	Bstar2			star;

	star.comment("# Symmetry from bsym\n\n");

	int				i;
	Matrix3			mat = ref_view.matrix();
	Vector3<double>	v;

	BstarBlock&		block = star.add_block(sym.label().str());

	block[SYMMETRY_POINT_GROUP] = sym.label().str();
	block[SYMMETRY_PG_NUMBER] = to_string(sym.point());

	BstarLoop&		loop = block.add_loop();
	loop.tags()[SYMMETRY_AXIS_ORDER] = 0;
	loop.tags()[SYMMETRY_AXIS_X] = 1;
	loop.tags()[SYMMETRY_AXIS_Y] = 2;
	loop.tags()[SYMMETRY_AXIS_Z] = 3;

	for ( i=0; i<sym.operations(); i++ ) {
		v = mat * sym[i].axis();
		vector<string>&	vs = loop.add_row(4);
		vs[0] = to_string(sym[i].order());
		vs[1] = to_string(v[0]);
		vs[2] = to_string(v[1]);
		vs[3] = to_string(v[2]);
	}
	
	return star.write(filename.str());
}

/**
@brief 	Calculates symmetry matrices from a list of strings.
@param 	nsym			number of symmetry operators.
@param 	*symop			array of symmetry operator lines.
@param 	line_len		length of text line in the array.
@return float* 			a set of 12-value symmetry matrices.

	The list of strings is expected to be packed into a single character
	array with a fixed length for each string. Each string encodes a
	symmetry operation in terms of x, y and z operations in reciprocal
	space.

**/
float*		sym_matrices_from_text_list(int nsym, char* symop, int line_len)
{
	// Set up the symmetry operator matrix
	int 		i, j, k, l;
	float*		mat = new float[nsym*12];
	for ( i=0; i<nsym*12; i++ ) mat[i] = 0;

	char		op[200];
	for ( i=0; i<nsym; i++ ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG sym_matrices_from_text_list: Symmetry operator " << i+1 << ":" << endl;
		k = 0;
		for ( j=0; j<3; j++ ) {
			l = 0;
			memset(op, 0, line_len);
			while ( k<line_len && symop[i*line_len+k] != ',' ) {
				op[l] = tolower(symop[i*line_len+k]);
				k++;
				l++;
			}
			op[l] = 0;
			if ( strstr(op,"-x") ) mat[i*12+4*j] = -1;
			else if ( strstr(op,"x") ) mat[i*12+4*j] = 1;
			if ( strstr(op,"-y") ) mat[i*12+4*j+1] = -1;
			else if ( strstr(op,"y") ) mat[i*12+4*j+1] = 1;
			if ( strstr(op,"-z") ) mat[i*12+4*j+2] = -1;
			else if ( strstr(op,"z") ) mat[i*12+4*j+2] = 1;
			if ( strstr(op,"1/2") ) mat[i*12+4*j+3] = 0.5;
			if ( strstr(op,"1/4") ) mat[i*12+4*j+3] = 0.25;
			if ( strstr(op,"3/4") ) mat[i*12+4*j+3] = 0.75;
			if ( strstr(op,"1/3") ) mat[i*12+4*j+3] = 1.0/3.0;
			if ( strstr(op,"2/3") ) mat[i*12+4*j+3] = 2.0/3.0;
			if ( strstr(op,"1/6") ) mat[i*12+4*j+3] = 1.0/6.0;
			if ( strstr(op,"5/6") ) mat[i*12+4*j+3] = 5.0/6.0;
			k++;
			if ( verbose & VERB_DEBUG )
				cout << "|" << mat[i*12+4*j] << " " << mat[i*12+4*j+1]
					<< " " << mat[i*12+4*j+2] << "|   |" << mat[i*12+4*j+3] << "|" << endl;
		}
	}
	
	return 0;
}
