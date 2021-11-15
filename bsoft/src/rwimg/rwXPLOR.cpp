/**
@file	rwXPLOR.cpp
@brief	Routines to read and write XPLOR reflection files
@author Bernard Heymann
@date	Created: 19981229
@date	Modified: 20130429
**/

#include "rwXPLOR.h"
#include "linked_list.h"
#include "Complex.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Internal function prototypes
int 	readXPLORmap(Bimage* p, int readdata);
int 	readXPLORsf(Bimage* p, int readdata);
int 	writeXPLORmap(Bimage* p);
int 	writeXPLORsf(Bimage* p);

/**
@brief	Reading a XPLOR map or structure factor file format.
A minimal implementation.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int					error code (<0 means failure).
**/
int 	readXPLOR(Bimage* p, int readdata)
{
	int 			err = 0;
	string			ext(extension(p->file_name()));
	
	if ( ext.find("xpl") != string::npos || ext.find("cns") != string::npos ) err = readXPLORmap(p, readdata);
	else if ( ext.find("rfl") != string::npos ) err = readXPLORsf(p, readdata);
	else err = -1;

	return err;
}

/*
@brief	Reading a XPLOR map file format.
A map text format used in X-ray crystallography.
@param	*p		the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int				error code (<0 means failure).
**/
int 	readXPLORmap(Bimage* p, int readdata)
{
	ifstream		fimg(p->file_name());
	
    if ( fimg.fail() ) return -1;

    char    	    aline[MAXLINELEN];
	long			i, j, x, y, z, xs(1), xf(1), ys(1), yf(1), zs(1), zf(1), n(0);
	double			avg, std;
	UnitCell		uc;

	p->images(1);
	
    fimg.getline(aline, MAXLINELEN);	    // First line blank
    fimg.getline(aline, MAXLINELEN);	    // Second line with number of title lines
	sscanf(aline, "%ld", &n);
	for ( i=0; i<n; i++ ) {
		fimg.getline(aline, MAXLINELEN);	// Title line
		p->label(p->label() + aline);
	}
    fimg.getline(aline, MAXLINELEN);	    // Dimensions and origin
	sscanf(aline, "%ld %ld %ld %ld %ld %ld %ld %ld %ld", &x, &xs, &xf,
		&y, &ys, &yf, &z, &zs, &zf);
    fimg.getline(aline, MAXLINELEN);	    // Unit cell
	sscanf(aline, "%lf %lf %lf %lf %lf %lf",
		&uc[0], &uc[1], &uc[2], &uc[3], &uc[4], &uc[5]);
    fimg.getline(aline, MAXLINELEN);	    // Line with "ZYX"
	
	p->size(xf-xs+1, yf-ys+1, zf-zs+1);
	p->channels(1);
//	p->image->origin(-xs, -ys, -zs);
	p->origin(-xs, -ys, -zs);
	p->data_type(Float);
	uc.degrees_to_radians();
	p->unit_cell(uc);
	p->sampling(uc.a()/x, uc.b()/y, uc.c()/z);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readXPLORmap: size=" << p->size() << endl;
	
	if ( !readdata ) {
		fimg.close();
		return 0;
	}
	
	long	slice_size = p->sizeX()*p->sizeY();
	float*			data = (float *) p->data_alloc();

	for ( i=z=0; z<p->sizeZ(); z++ ) {
		fimg.getline(aline, MAXLINELEN);			// Slice number
//		cout << "reading slice:" << z << " " << aline << endl;
		for ( j=0; j<slice_size; j+=6 ) {
			fimg.getline(aline, MAXLINELEN);
//			cout << "reading data:" << i << " " << aline << endl;
			if ( !fimg.fail() ) {
				n = sscanf(aline, "%f %f %f %f %f %f", &data[i], &data[i+1], 
					&data[i+2], &data[i+3], &data[i+4], &data[i+5]);
				i += n;
			}
		}
	}

    fimg.getline(aline, MAXLINELEN);	    // -9999
    fimg.getline(aline, MAXLINELEN);	    // Average and standard deviation
	sscanf(aline, "%lf %lf", &avg, &std);
	p->average(avg);
	p->standard_deviation(std);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readXPLORmap: avg=" << p->average() << " " << p->standard_deviation() << endl;
	
	fimg.close();
	
	return 0;
}

/*
@brief	Reading a XPLOR structure factor file format.
A structure factor text format used in X-ray crystallography.
	A minimal implementation of the XPLOR reflection file format which
	recognizes only the indices, amplitudes and phases on a line.
@param	*p		the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int				error code (<0 means failure).
**/
int 	readXPLORsf(Bimage* p, int readdata)
{
    long			i = 0, j, h, k, l, x, y, z, n = 0;
	long			hmin = 0, hmax = 0, kmin = 0, kmax = 0, lmin = 0, lmax = 0;
	float			amp, phi;

	ifstream		fhkl(p->file_name());
	
    if ( fhkl.fail() ) return -1;

    char    	    aline[MAXLINELEN], tag[20];
    fhkl.getline(aline, MAXLINELEN);	    // First line with number of reflections
	sscanf( aline, "%s %ld", tag, &n);
	
    fhkl.getline(aline, MAXLINELEN);	    // Second line with anomalous tag
	
    while ( fhkl.getline(aline, MAXLINELEN) ) {
		if ( strstr(aline, "INDE") || strstr(aline, "inde") ) {
			i = sscanf( aline, "%s %ld %ld %ld %s %f %f", tag, &h, &k, &l,
				tag, &amp, &phi );		// A reflection line
    		if ( i ) {
				if ( hmin > h ) hmin = h;
				if ( kmin > k ) kmin = k;
				if ( lmin > l ) lmin = l;
				if ( hmax < h ) hmax = h;
				if ( kmax < k ) kmax = k;
				if ( lmax < l ) lmax = l;
	    		n++;
			}
		}
	}
			
	p->size(1,1,1);
	if ( p->sizeX() < (long)2*hmax + 1 ) p->sizeX(2*(hmax + 1));
	if ( p->sizeX() < 1 - (long)2*hmin ) p->sizeX(-2*(hmin - 1));
	if ( p->sizeY() < (long)2*kmax + 1 ) p->sizeY(2*(kmax + 1));
	if ( p->sizeY() < 1 - (long)2*kmin ) p->sizeY(-2*(kmin - 1));
	if ( p->sizeZ() < (long)2*lmax + 1 ) p->sizeZ(2*(lmax + 1));
	if ( p->sizeZ() < 1 - (long)2*lmin ) p->sizeZ(-2*(lmin - 1));
	
	p->data_type(Float);
	p->compound_type(TComplex);
	p->channels(2);

	fhkl.close();
	
	if ( !readdata ) return 0;
	
    Complex<float>*	data = (Complex<float>*) p->data_alloc();
	
	fhkl.open(p->file_name());
	
    if ( fhkl.fail() ) return -1;

	fhkl.getline(aline, MAXLINELEN);		// Ignore the header line

    while ( fhkl.getline(aline, MAXLINELEN) ) {
		if ( strstr(aline, "INDE") || strstr(aline, "inde") ) {
			i = sscanf( aline, "%s %ld %ld %ld %s %f %f", tag, &h, &k, &l,
				tag, &amp, &phi );		// A reflection line
    		if ( i ) {
				x = h;
				if ( x < 0 ) x += p->sizeX();
				y = k;
				if ( y < 0 ) y += p->sizeY();
				z = l;
				if ( z < 0 ) z += p->sizeZ();
				j = (z*p->sizeY() + y)*p->sizeX() + x;
				phi *= M_PI/180;
				data[j] = Complex<float>(amp*cos(phi), amp*sin(phi));
			}
    	}
	}

    fhkl.close();
    
    return 0;
}

/**
@brief	Writing a XPLOR map or structure factor file format.
A text format used in X-ray crystallography.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeXPLOR(Bimage* p)
{
	int 			err(0);
	string			ext(extension(p->file_name()));
	
	if ( ext.find("xpl") != string::npos || ext.find("cns") != string::npos ) err = writeXPLORmap(p);
	else if ( ext.find("rfl") != string::npos ) err = writeXPLORsf(p);
	
	return err;
}

/*
@brief	Writing a XPLOR structure factor file format.
A structure factor text format used in X-ray crystallography.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeXPLORmap(Bimage* p)
{
	ofstream		fimg(p->file_name());
	
	if ( fimg.fail() ) return -1;

	long			i, j, x, y, z, n;
	long			xs, xf, ys, yf, zs, zf;
	long			slice_size = p->sizeX()*p->sizeY();

	UnitCell		uc = p->unit_cell();
	
	x = (long) (uc.a()/p->sampling(0)[0] + 0.5);
	y = (long) (uc.b()/p->sampling(0)[1] + 0.5);
	z = (long) (uc.c()/p->sampling(0)[2] + 0.5);
//	xs = (long) -(0.5 + p->image->origin()[0]);
//	ys = (long) -(0.5 + p->image->origin()[1]);
//	zs = (long) -(0.5 + p->image->origin()[2]);
	xs = (long) -p->image->origin()[0];
	ys = (long) -p->image->origin()[1];
	zs = (long) -p->image->origin()[2];
	xf = (long) (xs + p->sizeX() - 1);
	yf = (long) (ys + p->sizeY() - 1);
	zf = (long) (zs + p->sizeZ() - 1);

	fimg << endl;			// First line blank
	vector<string> 		slist = split(p->label(), '\n');
	n = slist.size();
	fimg << setw(8) << n << endl;	    // Second line with number of title lines
	for ( auto it = slist.begin(); it != slist.end(); ++it )
		fimg << "REMARKS " << *it << endl;	// Title line
	fimg << right << setw(8) << x << setw(8) << xs << setw(8) << xf
		<< setw(8) << y << setw(8) << ys << setw(8) << yf
		<< setw(8) << z << setw(8) << zs << setw(8) << zf << endl;  // Dimensions and origin
	fimg << setw(12) << uc.a() << setw(12) << uc.b() << setw(12) << uc.c()
		<< setw(12) << uc.alpha()*180/M_PI << setw(12) << uc.beta()*180/M_PI
		<< setw(12) << uc.gamma()*180/M_PI << endl;	// Unit cell
	fimg << "ZYX" << endl;				// Line with "ZYX"
	
	for ( i=z=0; z<p->sizeZ(); z++ ) {
		fimg << setw(8) << z << endl;
		for ( j=n=0; j<slice_size; j++, i++ ) {
    		fimg << setw(12) << (*p)[i];
			n++;
			if ( n == 6 ) {
				n = 0;
				fimg << endl;
			}
		}
		if ( n > 0 ) fimg << endl;
	}
 
 	fimg << "-9999" << endl;							// -9999
	fimg << setw(12) << p->average() << setw(12) << p->standard_deviation() << endl;	// Average and standard deviation
	
	fimg.close();
	
    return 0;
}

/*
@brief	Writing a XPLOR structure factor file format.
A structure factor text format used in X-ray crystallography.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeXPLORsf(Bimage* p)
{
    long   			i, x, y, z, n = 0, nw = 0, freeflag = 0;
	long			data_size = p->sizeX()*p->sizeY()*p->sizeZ();
    long			h, k, l;
	float			fom1 = 1;
	float*			fom = NULL;
	if ( p->next ) fom = (float *) p->next->data_pointer();
	Complex<float>*	data = (Complex<float> *) p->data_pointer();

	ofstream		fhkl(p->file_name());
	
	if ( fhkl.fail() ) return -1;

	for ( i=0; i<data_size; i++ ) if ( data[i].amp() >= 1e-10 ) n++;

    fhkl << " NREFlection=" << n << endl;
    fhkl << " ANOMalous=FALSe { equiv. to HERMitian=TRUE}" << endl;
/* Not necessary
    fhkl << " DECLare NAME=FOBS  DOMAin=RECIprocal  TYPE=REAL END" << endl;
    fhkl << " DECLare NAME=PHASE DOMAin=RECIprocal  TYPE=REAL END" << endl;
    fhkl << " DECLare NAME=FOM   DOMAin=RECIprocal  TYPE=REAL END" << endl;
    fhkl << " DECLare NAME=TEST  DOMAin=RECIprocal  TYPE=INTEGER END" << endl;
*/
	for ( i=z=0; z<p->sizeZ(); z++ ) {
		l = z;
		if ( l > (long)(p->sizeZ() - 1)/2 ) l -= p->sizeZ();
		for ( y=0; y<p->sizeY(); y++ ) {
			k = y;
			if ( k > (long)(p->sizeY() - 1)/2 ) k -= p->sizeY();
			for ( x=0; x<p->sizeX(); x++, i++ ) {
				h = x;
				if ( h > (long)(p->sizeX() - 1)/2 ) h -= p->sizeX();
				fom1 = 1;
				if ( data[i].amp() < 1e-10 ) fom1 = 0;
				if ( fom ) fom1 = fom[i];
				if ( fom1 > 1e-10 ) {
    				fhkl << " INDEX " << h << " " << k << " " << l 
						<< "  FOBS=" << data[i].amp() 
						<< "  PHASE=" << data[i].phi()*180/M_PI 
						<< "  FOM=" << fom1 << "  TEST=" << freeflag << endl;
    	    		nw++;
				}
			}
    	}
    }
    
    fhkl.close();
    
    return 0;
}
