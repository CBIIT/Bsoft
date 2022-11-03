/**
@file	rwresprop.cpp
@brief	Library routines to read and write residue properties
@author Bernard Heymann
@date	Created: 19991114
@date	Modified: 20210328
**/

#include "rwresprop.h"
#include "star.h"
#include "mol_tags.h"
#include "seq_util.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
Bresidue_type* 	residue_type_add(Bresidue_type** rtype, char c);
Bresidue_type* 	read_residue_prop_star(Bstring& filename);
Bresidue_type* 	read_residue_properties(Bstring& filename);
int 			write_residue_prop_star(Bstring& filename, Bresidue_type* rt);
Bresidue_matrix*	read_residue_mat_star(Bstring& filename);


// Masses from Amos Bairoch (SWISS-PROT)
// Volumes from:
//		Tsai J, Taylor R, Chothia C, Gerstein M.
//		The packing density in proteins: standard radii and volumes.
//		J Mol Biol. 1999 Jul 2;290(1):253-66.
// Side-chain extensions according to Walther et al. 1996, 255, 536-553.
//		only for ACFGILMVW
//		others are guesses with sigma = 0.5
// Hydrophobicity according to the GES scale
struct resprop { char code; float mass, vol, ext, extsd, charge, hphob; }
defaultprops[] = {
	{'-',   0.00,   0.0,  0.0,  0.00,  0.0,  0.0},
	{'*',   0.00,   0.0,  0.0,  0.00,  0.0,  0.0},
	{'A',  89.09,  89.0,  3.4,  0.10,  0.0, -1.6},
	{'B', 132.61, 118.0,  5.5,  0.50, -0.5,  7.0},
	{'C', 121.15, 103.0,  4.4,  0.24,  0.0, -2.0},
	{'D', 133.10, 114.0,  5.5,  0.50, -1.0,  9.2},
	{'E', 147.13, 139.0,  6.1,  0.50, -1.0,  8.2},
	{'F', 165.19, 199.0,  6.7,  0.32,  0.0, -3.7},
	{'G',  75.07,  64.0,  2.3,  0.07,  0.0, -1.0},
	{'H', 155.16, 157.0,  6.0,  0.50,  0.5,  3.0},
	{'I', 131.17, 163.0,  5.1,  0.22,  0.0, -3.1},
	{'K', 146.19, 165.0,  6.7,  0.50,  1.0,  8.8},
	{'L', 131.17, 163.0,  5.5,  0.20,  0.0, -2.8},
	{'M', 149.21, 166.0,  6.1,  0.45,  0.0, -3.4},
	{'N', 132.12, 122.0,  5.5,  0.50,  0.0,  4.8},
	{'P', 115.13, 122.0,  4.4,  0.50,  0.0,  0.2},
	{'Q', 146.15, 147.0,  6.1,  0.50,  0.0,  4.1},
	{'R', 174.23, 191.0,  7.3,  0.50,  1.0, 12.0},
	{'S', 105.09,  94.0,  4.4,  0.50,  0.0, -0.6},
	{'T', 119.12, 120.0,  4.5,  0.50,  0.0, -1.2},
	{'V', 117.15, 138.0,  4.5,  0.12,  0.0, -2.6},
	{'W', 204.23, 226.0,  7.4,  0.61,  0.0, -1.9},
	{'X', 128.16, 144.0,  5.5,  0.50,  0.0,  0.0},
	{'Y', 181.19, 195.0,  7.3,  0.50,  0.0,  0.7},
	{'Z', 146.64, 143.0,  6.1,  0.50, -0.5,  6.2},
} ;

char defaultcode[] = "-*ABCDEFGHIKLMNPQRSTVWXYZ";
float defaultsim[] = {
	1.0,1.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,
	1.0,1.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,
	-4.0,-4.0,4.0,-2.0,0.0,-2.0,-1.0,-2.0,0.0,-2.0,-1.0,-1.0,-1.0,-1.0,-2.0,-1.0,-1.0,-1.0,1.0,0.0,0.0,-3.0,0.0,-2.0,-1.0,
	-4.0,-4.0,-2.0,4.0,-3.0,4.0,1.0,-3.0,-1.0,0.0,-3.0,0.0,-4.0,-3.0,3.0,-2.0,0.0,-1.0,0.0,-1.0,-3.0,-4.0,-1.0,-3.0,1.0,
	-4.0,-4.0,0.0,-3.0,9.0,-3.0,-4.0,-2.0,-3.0,-3.0,-1.0,-3.0,-1.0,-1.0,-3.0,-3.0,-3.0,-3.0,-1.0,-1.0,-1.0,-2.0,-2.0,-2.0,-3.0,
	-4.0,-4.0,-2.0,4.0,-3.0,6.0,2.0,-3.0,-1.0,-1.0,-3.0,-1.0,-4.0,-3.0,1.0,-1.0,0.0,-2.0,0.0,-1.0,-3.0,-4.0,-1.0,-3.0,1.0,
	-4.0,-4.0,-1.0,1.0,-4.0,2.0,5.0,-3.0,-2.0,0.0,-3.0,1.0,-3.0,-2.0,0.0,-1.0,2.0,0.0,0.0,-1.0,-2.0,-3.0,-1.0,-2.0,4.0,
	-4.0,-4.0,-2.0,-3.0,-2.0,-3.0,-3.0,6.0,-3.0,-1.0,0.0,-3.0,0.0,0.0,-3.0,-4.0,-3.0,-3.0,-2.0,-2.0,-1.0,1.0,-1.0,3.0,-3.0,
	-4.0,-4.0,0.0,-1.0,-3.0,-1.0,-2.0,-3.0,6.0,-2.0,-4.0,-2.0,-4.0,-3.0,0.0,-2.0,-2.0,-2.0,0.0,-2.0,-3.0,-2.0,-1.0,-3.0,-2.0,
	-4.0,-4.0,-2.0,0.0,-3.0,-1.0,0.0,-1.0,-2.0,8.0,-3.0,-1.0,-3.0,-2.0,1.0,-2.0,0.0,0.0,-1.0,-2.0,-3.0,-2.0,-1.0,2.0,0.0,
	-4.0,-4.0,-1.0,-3.0,-1.0,-3.0,-3.0,0.0,-4.0,-3.0,4.0,-3.0,2.0,1.0,-3.0,-3.0,-3.0,-3.0,-2.0,-1.0,3.0,-3.0,-1.0,-1.0,-3.0,
	-4.0,-4.0,-1.0,0.0,-3.0,-1.0,1.0,-3.0,-2.0,-1.0,-3.0,5.0,-2.0,-1.0,0.0,-1.0,1.0,2.0,0.0,-1.0,-2.0,-3.0,-1.0,-2.0,1.0,
	-4.0,-4.0,-1.0,-4.0,-1.0,-4.0,-3.0,0.0,-4.0,-3.0,2.0,-2.0,4.0,2.0,-3.0,-3.0,-2.0,-2.0,-2.0,-1.0,1.0,-2.0,-1.0,-1.0,-3.0,
	-4.0,-4.0,-1.0,-3.0,-1.0,-3.0,-2.0,0.0,-3.0,-2.0,1.0,-1.0,2.0,5.0,-2.0,-2.0,0.0,-1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,-1.0,
	-4.0,-4.0,-2.0,3.0,-3.0,1.0,0.0,-3.0,0.0,1.0,-3.0,0.0,-3.0,-2.0,6.0,-2.0,0.0,0.0,1.0,0.0,-3.0,-4.0,-1.0,-2.0,0.0,
	-4.0,-4.0,-1.0,-2.0,-3.0,-1.0,-1.0,-4.0,-2.0,-2.0,-3.0,-1.0,-3.0,-2.0,-2.0,7.0,-1.0,-2.0,-1.0,-1.0,-2.0,-4.0,-2.0,-3.0,-1.0,
	-4.0,-4.0,-1.0,0.0,-3.0,0.0,2.0,-3.0,-2.0,0.0,-3.0,1.0,-2.0,0.0,0.0,-1.0,5.0,1.0,0.0,-1.0,-2.0,-2.0,-1.0,-1.0,3.0,
	-4.0,-4.0,-1.0,-1.0,-3.0,-2.0,0.0,-3.0,-2.0,0.0,-3.0,2.0,-2.0,-1.0,0.0,-2.0,1.0,5.0,-1.0,-1.0,-3.0,-3.0,-1.0,-2.0,0.0,
	-4.0,-4.0,1.0,0.0,-1.0,0.0,0.0,-2.0,0.0,-1.0,-2.0,0.0,-2.0,-1.0,1.0,-1.0,0.0,-1.0,4.0,1.0,-2.0,-3.0,0.0,-2.0,0.0,
	-4.0,-4.0,0.0,-1.0,-1.0,-1.0,-1.0,-2.0,-2.0,-2.0,-1.0,-1.0,-1.0,-1.0,0.0,-1.0,-1.0,-1.0,1.0,5.0,0.0,-2.0,0.0,-2.0,-1.0,
	-4.0,-4.0,0.0,-3.0,-1.0,-3.0,-2.0,-1.0,-3.0,-3.0,3.0,-2.0,1.0,1.0,-3.0,-2.0,-2.0,-3.0,-2.0,0.0,4.0,-3.0,-1.0,-1.0,-2.0,
	-4.0,-4.0,-3.0,-4.0,-2.0,-4.0,-3.0,1.0,-2.0,-2.0,-3.0,-3.0,-2.0,-1.0,-4.0,-4.0,-2.0,-3.0,-3.0,-2.0,-3.0,11.0,-2.0,2.0,-3.0,
	-4.0,-4.0,0.0,-1.0,-2.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-2.0,-1.0,-1.0,0.0,0.0,-1.0,-2.0,-1.0,-1.0,-1.0,
	-4.0,-4.0,-2.0,-3.0,-2.0,-3.0,-2.0,3.0,-3.0,2.0,-1.0,-2.0,-1.0,-1.0,-2.0,-3.0,-1.0,-2.0,-2.0,-2.0,-1.0,2.0,-1.0,7.0,-2.0,
	-4.0,-4.0,-1.0,1.0,-3.0,1.0,4.0,-3.0,-2.0,0.0,-3.0,1.0,-3.0,-1.0,0.0,-1.0,3.0,0.0,0.0,-1.0,-2.0,-3.0,-1.0,-2.0,4.0,
} ;

/**
@brief 	Reading residue properties from parameter files.
@param 	&filename		file name.
@return Bresidue_type*			the residue property structure, NULL on failure.
**/
Bresidue_type*	get_residue_properties(Bstring& filename)
{
	if ( verbose & VERB_DEBUG )	
		cout << "DEBUG get_residue_properties: Initializing atomic parameters" << endl;
	
	Bresidue_type*	rt_curr = NULL;
	Bresidue_type*	rt_first = NULL;
	
	int 			i;
		
	// Residue parameter file
	Bstring			rtfile;
	if ( filename.length() ) rtfile = filename;
	else rtfile = "res_prop.star";
	
	Bstring			propfile = parameter_file_path(rtfile);
    Bstring   		ext = rtfile.extension();
	if ( ext.length() ) {
		if ( ext.contains("star") )
			rt_first = read_residue_prop_star(propfile);
		else
			rt_first = read_residue_properties(propfile);
	}
	
	if ( !rt_first ) {
		if ( verbose )
			cerr << "Warning: Residue property file " << rtfile << " not opened! Using default properties" << endl;
		for ( i=0; i<MAXRES; i++ ) {
			rt_curr = residue_type_add(&rt_curr, defaultprops[i].code);
			if ( i == 1 ) rt_first = rt_curr;
			rt_curr->mass = defaultprops[i].mass;
			rt_curr->vol = defaultprops[i].vol;
			rt_curr->ext = defaultprops[i].ext;
			rt_curr->extsd = defaultprops[i].extsd;
			rt_curr->charge = defaultprops[i].charge;
			rt_curr->hphob = defaultprops[i].hphob;
		}
	}
	
	for ( i=0, rt_curr = rt_first; rt_curr; rt_curr = rt_curr->next, i++ ) ;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_residue_properties: Number of residue types = " << i << endl;
		
    return rt_first;
}

/**
@brief 	Writing residue properties from a parameter file.
@param 	&filename		file name.
@param 	*rt	the residue property structure.
@return int					error code (<0 means failure).
**/
int 			write_residue_properties(Bstring& filename, Bresidue_type* rt)
{
	int		err = write_residue_prop_star(filename, rt);
	
	return err;
}

/**
@brief 	Reading a residue matrix from a file.
@param 	&filename		file name.
@return Bresidue_matrix*		residue matrix property structure, NULL on failure.
**/
Bresidue_matrix*	get_residue_matrix(Bstring& filename)
{
	int					i;
	Bresidue_matrix*	resmat = NULL;
	
	// Matrix file
	Bstring				matfile;
	if ( filename.c_str() ) matfile = filename;
	else matfile = "blosum62.star";
	
	Bstring				propfile = parameter_file_path(matfile);
    Bstring				ext = matfile.extension();
	if ( ext.length() ) {
		if ( ext.contains("star") )
			resmat = read_residue_mat_star(propfile);
	}
		
	if ( !resmat ) {
		resmat = new Bresidue_matrix;
		memset(resmat, 0, sizeof(Bresidue_matrix));
		resmat->n = MAXRES;
		resmat->c = defaultcode;
		resmat->m = new float[resmat->n*resmat->n];
		for ( i=0; i<resmat->n*resmat->n; i++ ) resmat->m[i] = defaultsim[i];
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG get_residue_matrix: Matrix from file " << matfile << endl;
		for ( int i=0; i<resmat->n; i++ ) {
			cout << resmat->c[i];
			for ( int j=0; j<resmat->n; j++ )
				cout << tab << resmat->m[resmat->n*j+i];
			cout << endl;
		}
	}
	
	return resmat;
}

/**
@brief 	Deallocating a residue matrix.
@param 	*resmat	the residue matrix property structure.
@return int						0.
**/
int				residue_matrix_kill(Bresidue_matrix* resmat)
{
	if ( !resmat ) return 0;
	
	resmat->c = 0;
	
	if ( resmat->m ) delete[] resmat->m;
	
	delete resmat;
	
	return 0;
}

Bresidue_type*  residue_type_add(Bresidue_type** rtype, char c)
{
	Bresidue_type*          this_rtype = *rtype;
	Bresidue_type*          new_rtype = new Bresidue_type;
	memset(new_rtype, 0, sizeof(Bresidue_type));

	new_rtype->c = c;
	getcode3(c, new_rtype->cod);

	if ( !this_rtype )
		*rtype = new_rtype;
	else {
		while ( this_rtype->next ) this_rtype = this_rtype->next;
		this_rtype->next = new_rtype;
	}

	return new_rtype;
}

Bresidue_type* 	read_residue_prop_star(Bstring& filename)
{
 	Bstar			star;
	
 	if ( star.read(filename.str()) < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return NULL;
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_residue_prop_star: " << filename << endl;

	int				i, j, n(0);
	Bresidue_type*	rt = NULL;
	Bresidue_type*	rt_first = NULL;

	for ( auto ib: star.blocks() ) {
		for ( auto il: ib.loops() ) {
			if ( ( i = il.find(RESPROP_CODE1) ) >= 0 ) {
				for ( auto ir: il.data() ) {
					rt = residue_type_add(&rt, ir[i][0]);
					if ( !rt_first ) rt_first = rt;
					if ( ( j = il.find(RESPROP_MASS) ) >= 0 )
						rt->mass = to_real(ir[j]);
					if ( ( j = il.find(RESPROP_VOLUME) ) >= 0 )
						rt->vol = to_real(ir[j]);
					if ( ( j = il.find(RESPROP_EXTENSION) ) >= 0 )
						rt->ext = to_real(ir[j]);
					if ( ( j = il.find(RESPROP_EXTSTDEV) ) >= 0 )
						rt->extsd = to_real(ir[j]);
					if ( ( j = il.find(RESPROP_HYDROPHOBICITY) ) >= 0 )
						rt->hphob = to_real(ir[j]);
					if ( ( j = il.find(RESPROP_CHARGE) ) >= 0 )
						rt->charge = to_real(ir[j]);
					if ( ( j = il.find(RESPROP_H) ) >= 0 )
						rt->comp[0] = to_real(ir[j]);
					if ( ( j = il.find(RESPROP_C) ) >= 0 )
						rt->comp[1] = to_real(ir[j]);
					if ( ( j = il.find(RESPROP_N) ) >= 0 )
						rt->comp[2] = to_real(ir[j]);
					if ( ( j = il.find(RESPROP_O) ) >= 0 )
						rt->comp[3] = to_real(ir[j]);
					if ( ( j = il.find(RESPROP_S) ) >= 0 )
						rt->comp[4] = to_real(ir[j]);
					n++;
				}
			}
		}
	}
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_residue_prop_star: n=" << n << endl;
		
	return rt_first;
}

int 	write_residue_prop_star(Bstring& filename, Bresidue_type* rt)
{
	string			id("Residue_properties");
	Bstar			star;
	Bresidue_type*	t = NULL;

	star.comment("# Residue properties\n\n");

	BstarBlock&		block = star.add_block(id);

	BstarLoop&		loop = block.add_loop();
	loop.tags()[RESPROP_CODE1] = 0;
	loop.tags()[RESPROP_MASS] = 1;
	loop.tags()[RESPROP_VOLUME] = 2;
	loop.tags()[RESPROP_EXTENSION] = 3;
	loop.tags()[RESPROP_EXTSTDEV] = 4;
	loop.tags()[RESPROP_HYDROPHOBICITY] = 5;
	loop.tags()[RESPROP_CHARGE] = 6;
	loop.tags()[RESPROP_H] = 7;
	loop.tags()[RESPROP_C] = 8;
	loop.tags()[RESPROP_N] = 9;
	loop.tags()[RESPROP_O] = 10;
	loop.tags()[RESPROP_S] = 11;

	for ( t = rt; t; t = t->next ) {
		vector<string>&	vs = loop.add_row(12);
		vs[0] = t->c;
		vs[1] = to_string(t->mass);
		vs[2] = to_string(t->vol);
		vs[3] = to_string(t->ext);
		vs[4] = to_string(t->extsd);
		vs[5] = to_string(t->hphob);
		vs[6] = to_string(t->charge);
		vs[7] = to_string(t->comp[0]);
		vs[8] = to_string(t->comp[1]);
		vs[9] = to_string(t->comp[2]);
		vs[10] = to_string(t->comp[3]);
		vs[11] = to_string(t->comp[4]);
	}
	
	return star.write(filename.str());
}

Bresidue_type* 	read_residue_properties(Bstring& filename)
{
	Bresidue_type*	rt_curr = NULL;
	Bresidue_type*	rt_first = NULL;
	
    // Open property file read only
    FILE    	    *fprop;
    if ( ( fprop = fopen(filename.c_str(), "r") ) == NULL ) return NULL;
    
    if ( verbose )
		cout << "Reading file:                   " << filename << endl;
	
    // The property list column headers are recognized by one of a set of
	// expected strings (such "Mass" or "MW" or "VOL")
    int         i, n;
	int 		headflag(0), nhead(0), nrow(0);
	int 		imass = -1, ivol = -1, iext = -1, iextsd = -1, icharge = -1, ihphob = -1;
	float		val[32];
	int 		nrechead = 8;
	const char*	rechead[] = {" MW"," MAS"," VOL"," EXT"," EXTSD"," HP"," HYD"," CH"}; 
	
	// Parse the file
    char	   	aline[128], *aptr;
	Bstring		string, *itemlist, *item;
    while ( fgets(aline, 100, fprop) != NULL ) {
		string = aline;
		if ( string.length() < 3 ) 			// Blank lines = end of section
			headflag = 0;
		if ( headflag > 1 ) {				// Read a line of properties
			rt_curr = residue_type_add(&rt_curr, aline[0]);
			if ( !rt_first ) rt_first = rt_curr;
			aptr = aline + 1;
	    	for ( i=0; i<nhead; i++ ) {
				sscanf(aptr, "%f%n", &val[i], &n);
				aptr += n;
			}
			if ( imass > -1 ) rt_curr->mass = val[imass];
			if ( ivol > -1 ) rt_curr->vol = val[ivol];
			if ( iext > -1 ) rt_curr->ext = val[iext];
			if ( iextsd > -1 ) rt_curr->extsd = val[iextsd];
			if ( icharge > -1 ) rt_curr->charge = val[icharge];
			if ( ihphob > -1 ) rt_curr->hphob = val[ihphob];
			nrow++;
	    }
		if ( headflag == 0 ) 	// Check for property headers
			for ( i=0; i<nrechead; i++ )
				if ( string.contains(rechead[i]) ) headflag = 1;
		if ( headflag == 1 ) {
			itemlist = string.split();
			for ( item = itemlist, nhead = 0; item; item = item->next, nhead++ ) {
				if ( item->contains("MW") || item->contains("MAS") ) imass = nhead;
				if ( item->contains("VOL") ) ivol = nhead;
				if ( item->contains("EXT") ) iext = nhead;
				if ( item->contains("EXTSD") ) iextsd = nhead;
				if ( item->contains("HP") || item->contains("HYD") ) ihphob = nhead;
				if ( item->contains("CH") || item->contains("Q") ) icharge = nhead;
			}
			string_kill(itemlist);
			headflag = 100;
			nrow = 0;
		}
	}
	
    fclose(fprop);
    
    return rt_first;
}

Bresidue_matrix*	read_residue_mat_star(Bstring& filename)
{
/*	int 				nres2 = item_get_number(star, RESMATRIX_SUBSTITUTION);
	int 		nres = (int) sqrt(nres2);
	char*		mat_code = item_get_string_array(star, RESMATRIX2_CODE1, -1);
	char*		code1 = new char[nres+1];
	for ( i=0; i<nres; i++ ) code1[i] = mat_code[2*i];
	float*		sim_mat = item_get_float_array(star, RESMATRIX_SUBSTITUTION, -1);*/
	return NULL;
}


