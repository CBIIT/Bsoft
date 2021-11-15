/**
@file	rwkernel.cpp
@brief	Library routines to read and write kernel files
@author Bernard Heymann
@date	Created: 20031102
@date 	Modified: 20120211
**/

#include "rwkernel.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reads a text kernel file.
@param	*p			the image structure.
@return	int			0, <0 on error.
**/
int 		readKernel(Bimage *p)
{
	string		thisname;
	
	if ( access(p->file_name().c_str(), R_OK) == 0 )
		thisname = p->file_name();
	else
		thisname = parameter_file_path(p->file_name());

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readKernel: filename=" << thisname << endl;

    ifstream    	fkernel(thisname);
    if ( fkernel.fail() ) {
		error_show(p->file_name(), __FILE__, __LINE__);
		return -1;
	}
	
	p->file_name(thisname);
	p->data_type(Float);
	
	int				i;
	Vector3<long>	kernelsize;

	for ( i=0; i<3 && !fkernel.eof(); i++ )
		fkernel >> kernelsize[i];
		
	if ( kernelsize.volume() < 1 ) {
		cerr << "Error: The kernel size is wrong!" << endl;
		return -1;
	}
	
	p->size(kernelsize);
	p->images(1);
	p->channels(1);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readKernel: kernel size=" << p->size() << endl;

	p->origin(p->sizeX()/2, p->sizeY()/2, p->sizeZ()/2);
	
	long			nkernel = p->sizeX()*p->sizeY()*p->sizeZ();
	float*			kernel = (float *) p->data_alloc();
	
	for ( i=0; i<nkernel && !fkernel.eof(); i++ )
		fkernel >> kernel[i];
		
	fkernel.close();
	
	if ( i < nkernel ) {
		cerr << "Error: The kernel does not contain all the required values!" << endl;
		return -1;
	}

	if ( verbose & VERB_DEBUG ) {
		for ( i=0; i<nkernel; i++ )
			cout << kernel[i] << " ";
		cout << endl;
	}
		
	return 0;
}

/**
@brief	Writes a text kernel file.
@param	*p			the image structure.
@return	int			0.
**/
int 		writeKernel(Bimage *p)
{
	p->change_type(Float);
	
    ofstream    	fkernel(p->file_name().c_str());
    if ( fkernel.fail() ) {
		error_show(p->file_name().c_str(), __FILE__, __LINE__);
		return -1;
	}
	
	int			i, x, y, z;
	
	fkernel << p->sizeX() << " " << p->sizeY() << " " << p->sizeZ() << endl;
	for ( i=z=0; z<p->sizeZ(); z++ ) {
		for ( y=0; y<p->sizeY(); y++ ) {
			for ( x=0; x<p->sizeX(); x++, i++ )
				fkernel << (*p)[i] << " ";
			fkernel << endl;
		}
	}
	fkernel << endl;
	
	fkernel.close();
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG write_kernel_file: kernel size=" << p->size() << endl; 
		for ( i=0; i<p->sizeX()*p->sizeY()*p->sizeZ(); i++ )
			cout << (*p)[i] << " ";
		cout << endl;
	}
	
	return 0;
}

