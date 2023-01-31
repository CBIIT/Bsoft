/**
@file	rwASCII.cpp
@brief	Functions for reading and writing ASCII files
@author Bernard Heymann
@date	Created: 20000318
@date	Modified: 20120211
**/

#include "rwASCII.h"
#include "Complex.h"
#include "linked_list.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading an ASCII or text image format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int			error code (<0 means failure).

This function reads an ASCII image file with up to five dimensions
	in the order:
			c (channels), x, y, z, n (number of images)
	Default data type is Float
	Data is given as real (R) or complex (R and I or A and P),
			and may include an optional FOM (F)
	Column labels:
			Images: C X Y Z N R I F
			Structure factors: H K L A P R I F
**/
int 	readASCII(Bimage* p, int readdata)
{
    ifstream		fimg;
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;
	
	long 			i, j, k, l, m = 1;
	double			v[16], min[16], max[16];
	long			order[16];
    char			aline[MAXLINELEN];
	Vector3<double>	ori;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readASCII: Start reading text image" << endl;
	
	// Initialize order and extrema arrays
	for ( i=0; i<16; i++ ) {
		order[i] = -1;
		min[i] = 1e30;
		max[i] = -1e30;
		v[i] = 0;
	}
	
	// The first line contains the column labels
	fimg.getline( aline, MAXLINELEN);
    for ( size_t i=0; i<strlen(aline); i++ ) {
    	if ( isspace(aline[i]) || ispunct(aline[i]) ) aline[i]=' ';
    	aline[i] = toupper(aline[i]);
    }
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readASCII: Header line: " << aline << endl;
	
	Bstring			headerline(aline);
	Bstring*		one;
	Bstring*		label = headerline.split();
	m = count_list((char *)label);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readASCII: " << m << " column labels:";
		for ( one=label; one; one=one->next ) cout << *one << " ";
		cout << endl;
	}
	
	if ( verbose & VERB_DEBUG ) {
		if ( p->fourier_type() == NoTransform ) cout << "DEBUG readASCII: An image" << endl;
		else cout << "DEBUG readASCII: A transform" << endl;
	}
	
	// Associate each label with an index and a value variable
	for ( one=label, i=0; one && i<m; one=one->next, i++ ) {
		if ( order[0] < 0 && (*one)[0] == 'X' ) order[0] = i;
		if ( order[0] < 0 && (*one)[0] == 'H' ) {
			order[0] = i;
			p->compound_type(TComplex);
			p->fourier_type(Standard);
			p->channels(2);
		}
		if ( order[1] < 0 && (*one)[0] == 'Y' ) order[1] = i;
		if ( order[1] < 0 && (*one)[0] == 'K' ) order[1] = i;
		if ( order[2] < 0 && (*one)[0] == 'Z' ) order[2] = i;
		if ( order[2] < 0 && (*one)[0] == 'L' ) order[2] = i;
		if ( order[3] < 0 && (*one)[0] == 'N' ) order[3] = i;
		if ( order[4] < 0 && (*one)[0] == 'R' ) order[4] = i;
		if ( order[4] < 0 && (*one)[0] == 'A' ) order[4] = i;
		if ( order[4] < 0 && *one == "VX" ) {
			order[4] = i;
			p->compound_type(TVector3);
			p->channels(3);
		}
		if ( order[4] < 0 && *one == "Red" ) {
			order[4] = i;
			p->compound_type(TRGB);
			p->channels(3);
		}
		if ( order[4] < 0 && *one == "Cyn" ) {
			order[4] = i;
			p->compound_type(TCMYK);
			p->channels(4);
		}
		if ( order[5] < 0 && (*one)[0] == 'I' ) order[5] = i;
		if ( order[5] < 0 && (*one)[0] == 'P' ) order[5] = i;
		if ( order[5] < 0 && *one == "VY" ) order[5] = i;
		if ( order[5] < 0 && *one == "Grn" ) order[5] = i;
		if ( order[5] < 0 && *one == "Mag" ) order[5] = i;
		if ( order[6] < 0 && (*one)[0] == 'P' ) order[6] = i;
		if ( order[6] < 0 && *one == "VZ" ) order[6] = i;
		if ( order[6] < 0 && *one == "Blu" ) order[6] = i;
		if ( order[6] < 0 && *one == "Yel" ) order[6] = i;
		if ( order[7] < 0 && (*one)[0] == 'P' ) order[7] = i;
		if ( order[7] < 0 && *one == "Alf" ) {
			order[7] = i;
			p->compound_type(TRGBA);
			p->channels(4);
		}
		if ( order[7] < 0 && *one == "blK" ) order[7] = i;
		if ( order[8] < 0 && (*one)[0] == 'F' ) order[8] = i;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readASCII: order: " << order[0] << " " << order[1] << " " << order[2] << " " << 
			order[3] << " " << order[4] << " " << order[5] << " " << order[6] << " " << order[7] << " " << order[8] << " (m=" << m << ")" << endl;
	
	// Pass once through the file to get the number of data points and extrema
	i = 0;
 	while ( !fimg.eof() ) {
		fimg.getline(aline, MAXLINELEN);
		if ( sscanf(aline, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &v[0], &v[1], &v[2], 
				&v[3], &v[4], &v[5], &v[6], &v[7], &v[8]) < m ) {
			cout << "Missing data at line " << i << endl;
		} else {
			for ( j=0; j<9; j++ ) {
				if ( min[j] > v[j] ) min[j] = v[j];
				if ( max[j] < v[j] ) max[j] = v[j];
			}
/*			if ( verbose & VERB_DEBUG ) {
				cout << "DEBUG readASCII: " << i;
				for ( j=0; j<m; j++ ) cout << " " << v[j];
				cout << endl;
			}*/
			i++;
		}
	}
	
	fimg.close();
	
	// Set some parameters for the image
	p->data_type(Float);
	p->channels(1);
	if ( order[3] > -1 ) p->images((long) (max[order[3]] + 1));
	else p->images(1);
	
	if ( p->fourier_type() == NoTransform ) {
		if ( order[0] > -1 ) {
			p->sizeX((long) (max[order[0]] - min[order[0]] + 1));
			ori[0] = min[order[0]];
		}
		if ( order[1] > -1 ) {
			p->sizeY((long) (max[order[1]] - min[order[1]] + 1));
			ori[1] = min[order[1]];
		}
		if ( order[2] > -1 ) {
			p->sizeZ((long) (max[order[2]] - min[order[2]] + 1));
			ori[2] = min[order[2]];
		}
		p->origin(ori);
	} else {
		if ( order[0] > -1 ) {
			p->sizeX((int) fabs(min[order[0]]));
			if ( p->sizeX() < fabs(max[order[0]]) ) p->sizeX((long) fabs(max[order[0]]));
			p->sizeX(p->sizeX() + 1);
		}
		if ( order[1] > -1 ) {
			p->sizeY((int) fabs(min[order[1]]));
			if ( p->sizeY() < fabs(max[order[1]]) ) p->sizeY((long) fabs(max[order[1]]));
			p->sizeY(p->sizeY() + 1);
		}
		if ( order[2] > -1 ) {
			p->sizeZ((int) fabs(min[order[2]]));
			if ( p->sizeZ() < fabs(max[order[2]]) ) p->sizeZ((long) fabs(max[order[2]]));
			p->sizeZ(p->sizeZ() + 1);
		}
	}
	
	for ( i=4; i<8; i++ ) {
		if ( order[i] > -1 ) {
			if ( p->minimum() > min[order[i]] ) p->minimum(min[order[i]]);
			if ( p->maximum() < max[order[i]] ) p->maximum(max[order[i]]);
		}
	}
	
	p->size(p->size().max(1));

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readASCII: size=" << p->size() << " channels=" << p->channels() << " images=" << p->images() << endl;

	if ( !readdata ) return 0;
	
	// Allocate memory for the data
	p->data_alloc();
	float*			fom = NULL;
	if ( order[8] > -1 ) {
		p->next = new Bimage(Float, TSimple, p->size(), p->images());
		fom = (float *) p->next->data_pointer();
	}

	// Rewind and read the data into the image
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;

	fimg.getline( aline, MAXLINELEN);	// Header labels

	long			c(0), x(0), y(0), z(0), n(0);
	i = 0;
	while ( !fimg.eof() ) {
		fimg.getline(aline, MAXLINELEN);
		if ( sscanf(aline, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &v[0], &v[1], &v[2], 
				&v[3], &v[4], &v[5], &v[6], &v[7], &v[8]) < m ) {
			cout << "Missing data at line " << i << endl;
		} else {
			if ( order[0] > -1 ) x = (long) v[order[0]];
			if ( order[1] > -1 ) y = (long) v[order[1]];
			if ( order[2] > -1 ) z = (long) v[order[2]];
			if ( order[3] > -1 ) n = (long) v[order[3]];
			if ( x < 0 ) x += p->sizeX(); 	// Wrapping for structure factors
			if ( y < 0 ) y += p->sizeY();
			if ( z < 0 ) z += p->sizeZ();
			j = ((n*p->sizeZ() + z)*p->sizeY() + y)*p->sizeX() + x;
			for ( c=0, k=j*p->channels(), l=4; c<p->channels(); c++, k++, l++ ) p->set(k, v[order[l]]);
			if ( order[8] > -1 ) fom[j] = v[order[8]];
			i++;
			if ( verbose & VERB_DEBUG ) {
				cout << "DEBUG readASCII: " << i;
				for ( j=0; j<m; j++ ) cout << " " << v[j];
				for ( c=0, k=j*p->channels(); c<p->channels(); c++, k++ ) cout << " " << (*p)[k];
				cout << endl;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Data points read:               " << i << endl << endl;
	
	fimg.close();
	
	string_kill(label);
	
	return 0;
}

/**
@brief	Writing an ASCII or text image format.
@param	*p			the image structure.
@return	int			error code (<0 means failure).

This function writes an ascii image file with up to five dimensions
	in the order:
			c (channels), x, y, z, n (number of images)
	Default data type is Float
	Data is given as real (R) or complex (R and I or A and P),
			and may include an optional FOM (F)
	Column labels:
			Images: C X Y Z N R I F
			Structure factors: H K L A P R I F
**/
int 	writeASCII(Bimage* p)
{
    ofstream        fimg;
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;
	
	long 			i(0), nlab(128);
	long 			order[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	char			label[nlab];
	for ( i=0; i<nlab; ++i ) label[i] = ' ';
	
	// Get the data order and labels
	i = 0;
	if ( p->sizeX() > 1 ) {
		order[0] = i;
		label[10*i] = 'X';
		if ( p->fourier_type() ) label[10*i] = 'H';
		i++;
	}
	if ( p->sizeY() > 1 ) {
		order[1] = i;
		label[10*i] = 'Y';
		if ( p->fourier_type() ) label[10*i] = 'K';
		i++;
	}
	if ( p->sizeZ() > 1 ) {
		order[2] = i;
		label[10*i] = 'Z';
		if ( p->fourier_type() ) label[10*i] = 'L';
		i++;
	}
	if ( p->images() > 1 ) {
		order[3] = i;
		strncpy(label + 10*i, "Nimage", 7);
		i++;
	}
	
	switch ( p->compound_type() ) {
		case TSimple: strncpy(label + 10*i, "Real", 5); i++; break;
		case TComplex: strncpy(label + 10*i, "Real", 5); i++;
			strncpy(label + 10*i, "Imag", 5); i++; break;
		case TVector2: strncpy(label + 10*i, "VX", 4); i++;
			strncpy(label + 10*i, "VY", 4); i++; break;
		case TVector3: strncpy(label + 10*i, "VX", 4); i++;
			strncpy(label + 10*i, "VY", 4); i++;
			strncpy(label + 10*i, "VX", 4); i++; break;
		case TView: strncpy(label + 10*i, "VX", 4); i++;
			strncpy(label + 10*i, "VY", 4); i++;
			strncpy(label + 10*i, "VX", 4); i++;
			strncpy(label + 10*i, "Ang", 4); i++; break;
		case TRGB: strncpy(label + 10*i, "Red", 4); i++;
			strncpy(label + 10*i, "Grn", 4); i++;
			strncpy(label + 10*i, "Blu", 4); i++; break;
		case TRGBA: strncpy(label + 10*i, "Red", 4); i++;
			strncpy(label + 10*i, "Grn", 4); i++;
			strncpy(label + 10*i, "Blu", 4); i++;
			strncpy(label + 10*i, "Alf", 4); i++; break;
		case TCMYK: strncpy(label + 10*i, "Cyn", 4); i++;
			strncpy(label + 10*i, "Mag", 4); i++;
			strncpy(label + 10*i, "Yel", 4); i++;
			strncpy(label + 10*i, "blK", 4); i++; break;
		default: strncpy(label + 10*i, "Real", 5); i++; break;
	}

	if ( p->next ) {
		order[7] = i;
		memcpy(label + 10*i, "FOM", 3);
		i++;
	}
//	strncpy(label + 10*i, "\0", 1);
	label[10*i] = 0;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG writeASCII: label: " << label << endl;
		cout << "DEBUG writeASCII: order: " << order[0] << " " << order[1] << " " << order[2] << " " << 
			order[3] << " " << order[4] << " " << order[5] << " " << order[6] << " " << order[7] << " " << order[8] << " (m=" << i << ")" << endl;
	}
	
	fimg << label << endl;
	
	long			c, x, y, z, n;
	float*			fom = (float *) p->next->data_pointer();
	
	for ( i=n=0; n<p->images(); n++ ) {
		for ( z=0; z<p->sizeZ(); z++ ) {
			for ( y=0; y<p->sizeY(); y++ ) {
				for ( x=0; x<p->sizeX(); x++ ) {
					if ( order[0] > -1 ) fimg << " " << x;
					if ( order[1] > -1 ) fimg << " " << y;
					if ( order[2] > -1 ) fimg << " " << z;
					if ( order[3] > -1 ) fimg << " " << n;
					for ( c=0; c<p->channels(); c++, i++ ) fimg << " " << (*p)[i];
					if ( order[4] > -1 ) fimg << " " << fom[i];
					fimg << endl;
				}
			}
		}
	}
	
	fimg.close();
	
	return 0;
}

