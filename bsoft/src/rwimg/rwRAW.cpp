/**
@file	rwRAW.cpp
@brief	Functions for reading and writing RAW files
@author Bernard Heymann
@date	Created: 19990724
@date 	Modified: 20150910

	A RAW file is a generalized image file defined as consisting of an optional header,
		and a contiguous, non-compressed block of data
**/

#include "rwRAW.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a raw block of data in a file.
@param	*p			the image structure.
@param 	img_select		image selection in multi-image file (-1 = all images).
@return	int 				0.
A 5D format for generalized access to data in Bsoft.
	File format extensions: .raw and #
							The special symbol '#' appended to a valid
							file name is taken to indicate that the
							file should be interpreted as a raw data
							block. To further specify how the data block
							should be read, several variable-value pairs
							can be appended, each pair started with a '#'.
	Variable-value pairs:	h=header_size
							d=data_type
							x=size_x,size_y,size_z
							p=page_size_x,page_size_y,page_size_z
							a=padding_bytes_between_pages
							s=sampling_x,sampling_y,sampling_z
							c=number_of_channels
							n=number_of_images
							i=image_selected
							f=transform_type (n=NoTransform, s=Standard, c=Centered, h=Hermitian, q=CentHerm)
							b=1=swap_bytes
							v=1=vax_floating_point
**/
int 	readRAW(Bimage* p, int img_select)
{
	// Get the part of the file name with the specifications	
//	Bstring			specs(p->file_name().post('#'));
	Bstring			specs(p->file_name().substr(p->file_name().find("#")+1));
//	string			specs(p->file_name().substr(p->file_name().find("#")+1));
//	transform(specs.begin(), specs.end(), specs.begin(), ::tolower);
	
	// Get the real filename without specifications
//	p->file_name(p->file_name().pre('#'));
	p->file_name(p->file_name().substr(0, p->file_name().find('#')));
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG rwRAW: Filename = " << p->file_name() << endl;
		cout << "DEBUG rwRAW: Specs = " << specs << endl;
	}
		
	// Set defaults
	p->data_offset(0);
	long			i(0), offset(0), pad(0), c(1), n(1);
	int				swap(0), vax(0);
	char			transform_type('n');
	Vector3<long>	size(1,1,1), pagesize(1,1,1);
	Vector3<double>	sam(1,1,1);
	
	// Extract the parameters from the filename
	Bstring			item, value;
	for ( specs = specs.lower(); specs.contains("="); specs = specs.post('#') ) {
		item = specs.pre('#');
		value = item.post('=');
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG rwRAW: item=" << item << " value=" << value << endl;
		if ( item[0] == 'h' ) sscanf(value.c_str(), "%ld", &offset);
		if ( item[0] == 'd' ) p->data_type(getdatatype(value[0]));
		if ( item[0] == 'x' ) sscanf(value.c_str(), "%ld,%ld,%ld", &size[0], &size[1], &size[2]);
		if ( item[0] == 'p' ) sscanf(value.c_str(), "%ld,%ld,%ld", &pagesize[0], &pagesize[1], &pagesize[2]);
		if ( item[0] == 's' ) i = sscanf(value.c_str(), "%lf,%lf,%lf", &sam[0], &sam[1], &sam[2]);
		if ( item[0] == 'c' ) sscanf(value.c_str(), "%ld", &c);
		if ( item[0] == 'n' ) sscanf(value.c_str(), "%ld", &n);
		if ( item[0] == 'i' ) sscanf(value.c_str(), "%d", &img_select);
		if ( item[0] == 'f' ) sscanf(value.c_str(), "%c", &transform_type);
		if ( item[0] == 'b' ) sscanf(value.c_str(), "%d", &swap);
		if ( item[0] == 'v' ) sscanf(value.c_str(), "%d", &vax);
		if ( item[0] == 'a' ) sscanf(value.c_str(), "%ld", &pad);
	}
	if ( i == 1 ) {
		if ( size[1] > 1 ) sam[1] = sam[0];
		if ( size[2] > 1 ) sam[2] = sam[0];
	}
	
	p->data_offset(offset);
	p->channels(c);
//	p->images(n);
	
	switch ( transform_type ) {
		case 'n': p->fourier_type(NoTransform); break;
		case 's': p->fourier_type(Standard); break;
		case 'c': p->fourier_type(Centered); break;
		case 'h': p->fourier_type(Hermitian); break;
		case 'q': p->fourier_type(CentHerm); break;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG rwRAW: Specs: offset=" << p->data_offset() << " data type=" << p->data_type() << endl;
		cout << "DEBUG rwRAW:        size=" << size << endl;
		cout << "DEBUG rwRAW:        channels=" << p->channels() << "images=" << n << endl; 
		cout << "DEBUG rwRAW:        transform=" << transform_type << " swap=" << swap << " vax=" << vax << " pad=" << pad << endl;
	}
		
	// Check parameters
	p->compound_type(TSimple);
	if ( p->channels() == 2 ) p->compound_type(TComplex);
	if ( p->channels() == 3 ) p->compound_type(TRGB);
	if ( p->channels() == 4 ) p->compound_type(TRGBA);
	if ( pad ) pagesize[2] = 1;
	
	// If only half of a transform is stored, it need to be handled
	long			xstore = size[0];
	if ( p->fourier_type() == Hermitian || p->fourier_type() == CentHerm )
		xstore = size[0]/2 + 1;
	
	if ( p->data_type() == Bit ) {
		xstore = (size[0] - 1)/8 + 1;
		pagesize[0] = 8*xstore;
	}
	
	ifstream*		fimg = new ifstream;
    fimg->open(p->file_name());
    if ( fimg->fail() ) return -1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwRAW: File opened" << endl;
		
	fimg->seekg(0, ios_base::end);
	long			len = fimg->tellg();
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwRAW: File length = " << len << endl;

	long			datasize = len - p->data_offset();
	long			datatypesize = p->data_type_size();
	long			specsize = xstore*size[1]*size[2]*p->channels()*n*datatypesize;
	
	if ( specsize > datasize ) {
		cout << "File size=" << datasize << " < spec size=" << specsize << " (" <<
			xstore << "," << size[1] << "," << size[2] << "," << p->channels() << "," << n << "," << datatypesize << ")" << endl;
		p->channels(1);
		size[2] = 1;
		size[0] = xstore = (long) sqrt(1.0*datasize);
		if ( p->fourier_type() == Hermitian || p->fourier_type() == CentHerm )
			size[0] = 2*(xstore - 1);
		size[1] = datasize/size[0];
		p->data_type(UCharacter);
		cout << "Assume a 2D image with 8 bit data of size " << p->sizeX() << "x" << p->sizeY() << endl;
	}

	if ( pagesize[0] < 2 ) pagesize[0] = size[0];
	if ( pagesize[1] < 2 ) pagesize[1] = size[1];
	
	if ( img_select > -1 ) {
		if ( img_select >= n ) img_select = n - 1;
		n = 1;
	} else img_select = -1;
	
	p->images(n);
	
	p->size(size);
	p->page_size(pagesize);
	p->sampling(sam);
	
	
	p->read_data(fimg, img_select, swap, vax, pad);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwRAW: File read" << endl;

	fimg->close();
	
	delete fimg;
	
	return 0;
}

/**
@brief	Writing a raw block of data in a file.
@param	*p			the image structure.
@return	int 				0.
A 5D format for generalized access to data in Bsoft.
**/
int 	writeRAW(Bimage* p)
{
	long 			datasize = p->alloc_size();
	
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write((char *)p->data_pointer(), datasize);
	
	fimg.close();
	
	return 0;
}

