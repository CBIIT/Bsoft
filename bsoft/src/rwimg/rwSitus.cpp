/**
@file	rwSitus.cpp
@brief	Functions for reading and writing Situs files
@author Bernard Heymann
@date	Created: 20150316
@date	Modified: 20150316
**/

#include "rwSitus.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a Situs image format.
@param	*p		the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int					error code (<0 means failure).
This function reads an Situs image file:
	The first line contains (in order):
		sampling
		-origin		3 values
		size		3 values
	Data type returned is Float
**/
int 	readSitus(Bimage* p, int readdata)
{
    ifstream		fimg;
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSitus: Start reading text image" << endl;
	
	size_t			i(0), n(0);
	double			v;
	Vector3<double>	ori;
	Vector3<long>	size;
	
	fimg >> v;
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSitus: Pixel size = " << v << endl;
	fimg >> ori[0] >> ori[1] >> ori[2];
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSitus: Origin = " << -ori << endl;
	fimg >> size[0] >> size[1] >> size[2];
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readSitus: Size = " << size << endl;
	
	p->size(size);
	p->images(1);
	p->channels(1);
	p->data_type(Float);
	p->sampling(v,v,v);
	p->origin(-ori);
	p->data_alloc();
	
	if ( readdata) {
		n = (size_t) size[0]*size[1]*size[2];
	
		while ( !fimg.eof() && i < n ) {
			fimg >> v;
			if ( fabs(v) < SMALLFLOAT ) v = 0;
			p->set(i++, v);
		}
	
		if ( i < n )
			cerr << "Error: The data read (" << i << ") is less than the size specified (" << n << ")" << endl;
	}
	
	fimg.close();
	
	return 0;
}

/**
@brief	Writing a Situs image format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeSitus(Bimage* p)
{
    ofstream        fimg;
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;
	
	long			i;
	
	fimg << p->sampling(0)[0];
	for ( i=0; i<3; i++ ) fimg << " " << -p->image->origin()[i];
	for ( i=0; i<3; i++ ) fimg << " " << p->size()[i];
	fimg << endl;
	
	for ( i=0; i<p->data_size(); i++ ) {
		if ( i%10 ) fimg << " ";
		else fimg << endl;
		fimg << (*p)[i];
	}
	
	fimg << endl;
	
	fimg.close();
	
	return 0;
}

