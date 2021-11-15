/**
@file	rwHKL.cpp
@brief	Routines to read and write HKL reflection files
@author Bernard Heymann
@date	Created: 19981229
@date	Modified: 20120211
**/

#include "rwHKL.h"
#include "Complex.h"
#include "linked_list.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a HKL structure factor file format.
This function reads a HKL file typically produced with MTZDUMP from
		X-ray crystallographic structure factor files in the CCP4 package.
	The first line as a space-delimited list of headers for the columns.
	All subsequent lines are rows containing the indices and data.
	Column headers:
		H K L		indices.
		AMP 		amplitudes.
		SIGAMP		amplitude deviations.
		PHI 		phases.
		SIGPHI		phase deviations.
		FOM 		figure-of-merit values.
		FREE		flag for R-free calculations
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int					error code (<0 means failure).
**/
int 	readHKL(Bimage* p, int readdata)
{
    ifstream		fhkl;
    long			i, j, m(1), n(0);

    fhkl.open(p->file_name());
    if ( fhkl.fail() ) return -1;

    char    	    aline[MAXLINELEN];
	
    fhkl.getline(aline, MAXLINELEN);	    // First line with column headers
    for ( i=0; i<strlen(aline); i++ ) {
    	if ( isspace(aline[i]) || ispunct(aline[i]) ) aline[i]=' ';
    	aline[i] = toupper(aline[i]);
    }
	
	Bstring			headerline(aline);
	Bstring*		one;
	Bstring*		header = headerline.split();
	m = count_list((char *)header);

    int		ih(0); 	    	    	// Find the H index header item
    int		ik(0); 	    	    	// Find the K index header item
    int		il(0); 	    	    	// Find the L index header item
    int		iamp(0);	    	    // Find AMP, PHI, FOM and FREE columns
 	int 	isigamp(0);
	int		iphi(0);
    int 	isigphi(0);
    int		ifom(0);
	int		ifre(0);

	for ( one=header, i=0; one && i<m; one=one->next, i++ ) {
		if ( *one == "H" || *one == "IH" ) ih = i;
		if ( *one == "K" || *one == "IK" ) ik = i;
		if ( *one == "L" || *one == "IL" ) il = i;
		if ( one->contains("AMP") ) iamp = i;
		if ( one->contains("SIGAMP") ) isigamp = i;
		if ( one->contains("PH") ) iphi = i;
		if ( one->contains("SIGPH") ) isigphi = i;
		if ( one->contains("FOM") ) ifom = i;
		if ( one->contains("FRE") ) ifre = i;
	}

 	if ( verbose & VERB_PROCESS ) {
	    cout << "Headers:                        ";
		for ( one=header; one; one=one->next ) cout << " " << *one;
    	cout << endl;
	}

    float   	v[64], fommax(0);  	    	    	// Variable input list
	char*		aptr;
	int 		k, hmin(0), hmax(0), kmin(0), kmax(0), lmin(0), lmax(0);
    while ( fhkl.getline(aline, MAXLINELEN) ) {
    	if ( strlen(aline) > 10 ) {
    		for ( i=0; i<64; ++i ) v[i] = 0;
			aptr = aline;
			i = 0;
			while ( strlen(aptr) > 1 ) {
				if ( sscanf(aptr, "%f%n", &v[i], &k) < 1 && i < m ) {
					fhkl.getline(aline, MAXLINELEN);
					aptr = aline;
					sscanf(aptr, "%f%n", &v[i], &k);
				}
				aptr += k;
				i++;
			}
    		if ( i ) {
				if ( hmin > v[ih] ) hmin = (int) v[ih];
				if ( ik < m && kmin > v[ik] ) kmin = (int) v[ik];
				if ( il < m && lmin > v[il] ) lmin = (int) v[il];
				if ( hmax < v[ih] ) hmax = (int) v[ih];
				if ( ik < m && kmax < v[ik] ) kmax = (int) v[ik];
				if ( il < m && lmax < v[il] ) lmax = (int) v[il];
				if ( ifom < m && fommax < v[ifom] ) fommax = v[ifom];
	    		n++;
			}
    	}
	}
	
	fhkl.close();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwHKL: First pass done" << endl;
	
	p->images(1);
	p->size(1,1,1);
	if ( p->sizeX() < (long)2*(hmax + 1) ) p->sizeX(2*(hmax + 1));
	if ( p->sizeX() < (long)2*(-hmin) ) p->sizeX(-2*hmin);
	if ( p->sizeY() < (long)2*(kmax + 1) ) p->sizeY(2*(kmax + 1));
	if ( p->sizeY() < (long)2*(-kmin) ) p->sizeY(-2*kmin);
	if ( p->sizeZ() < (long)2*(lmax + 1) ) p->sizeZ(2*(lmax + 1));
	if ( p->sizeZ() < (long)2*(-lmin) ) p->sizeZ(-2*lmin);
//	if ( p->FOM_maximum() < 1 ) p->FOM_maximum(1); 	// Preserve fractional FOM values
	if ( kmin == kmax ) p->sizeY(1);
	if ( lmin == lmax ) p->sizeZ(1);
	
	int 		x(0), y(0), z(0);
	p->data_type(Float);
	p->compound_type(TComplex);
	p->channels(2);
	p->fourier_type(Standard);
	

	if ( !readdata ) return 0;

    Complex<float>*	data = (Complex<float>*) p->data_alloc();
	float*		fom = NULL;
	if ( ifom ) {
		p->next = new Bimage(Float, TSimple, p->size(), p->images());
		fom = (float *) p->next->data_pointer();
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwHKL: Memory allocated" << endl;
	
	
    fhkl.open(p->file_name());
    if ( fhkl.fail() ) return -1;

	fhkl.getline(aline,128);		// Ignore the header line
	
    while ( fhkl.getline(aline, MAXLINELEN) ) {
		if ( verbose & VERB_DEBUG )
			cout << aline << endl;
    	if ( strlen(aline) > 10 ) {
    		for ( i=0; i<64; ++i ) v[i] = 0;
			aptr = aline;
			i = 0;
			while ( strlen(aptr) > 1 ) {
				if ( sscanf(aptr, "%f%n", &v[i], &k) < 1 && i < m ) {
					fhkl.getline(aline, MAXLINELEN);
					aptr = aline;
					sscanf(aptr, "%f%n", &v[i], &k);
				}
				aptr += k;
				i++;
			}
    		if ( i ) {
				x = (int) v[ih];
				if ( x < 0 ) x += p->sizeX();
				if ( ik ) {
					y = (int) v[ik];
					if ( y < 0 ) y += p->sizeY();
				}
				if ( il ) {
					z = (int) v[il];
					if ( z < 0 ) z += p->sizeZ();
				}
				j = (z*p->sizeY() + y)*p->sizeX() + x;
				if ( iphi ) {
					v[iphi] *= M_PI/180.0;
					data[j] = Complex<float>(v[iamp]*cos(v[iphi]), v[iamp]*sin(v[iphi]));
				}
				if ( p->next && ifom ) fom[j] = v[ifom]/fommax;
			}
    	}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwHKL: File read" << endl;	
	
    fhkl.close();
    
	string_kill(header);
	
    return 0;
}

/**
@brief	Writing a HKL structure factor file format.
This function writes a HKL structure factor file.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeHKL(Bimage* p)
{
	if ( p->compound_type() < TComplex ) p->simple_to_complex();
	p->change_type(Float);
	
    ofstream		fhkl;
    long   			i, x, y, z, nw(0), freeflag(0);
    long			h, k, l;
	float			fom1;
	float*			fom = NULL;
	if ( p->next ) fom = (float *) p->next->data_pointer();
	Complex<float>*	data = (Complex<float> *) p->data_pointer();

    fhkl.open(p->file_name());
    if ( fhkl.fail() ) return -1;
	
    fhkl << "  H   K   L   AMP  PHI  FOM  FREE" << endl;
	
	for ( i=z=0; z<p->sizeZ(); z++ ) {
		l = z;
		if ( l > (p->sizeZ() - 1)/2 ) l -= p->sizeZ();
		for ( y=0; y<p->sizeY(); y++ ) {
			k = y;
			if ( k > (p->sizeY() - 1)/2 ) k -= p->sizeY();
			for ( x=0; x<p->sizeX(); x++, i++ ) {
				h = x;
				if ( h > (p->sizeX() - 1)/2 ) h -= p->sizeX();
				fom1 = 1;
				if ( data[i].amp() < 1e-10 ) fom1 = 0;
				if ( fom ) fom1 = fom[i];
				if ( fom1 > 1e-10 ) {
		    	    fhkl << h << " " << k << " " << l << " " << data[i].amp() << " " 
						<< data[i].phi()*180/M_PI << " " << fom1 << " " << freeflag << endl;
    	    		nw++;
				}
			}
    	}
    }
    
    fhkl.close();
    
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG writeHKL: Number of structure factors written: " << nw << endl;
	}

    return 0;
}

