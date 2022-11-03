/**
@file	Bimage.cpp
@brief	Methods for the image class
@author Bernard Heymann
@date	Created: 20110603
@date 	Modified: 20220105
**/

#include "Bimage.h"
#include "utilities.h"

#include <iostream>
using namespace std;

// Definition of the global variables 
extern int 		verbose;		// Level of output to the screen

/**
@brief Initializes a sub-image.

	All values are set to zero except mag and fom, set to one.

**/
Bsub_image::Bsub_image()
{
	min = max = avg = std = 0;
	bkg = 0;
	ux = uy = uz = 1;
	ox = oy = oz = 0;
	vx = vy = vz = angle = 0;
	mag = 1;
	fom = 1;
	sel = 1;
//	cout << "sub-image initialized" << endl;
}

/**
@brief Initializes an image.

	The internal initialization function is called.

**/
Bimage::Bimage()
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage: Bimage size = " << sizeof(Bimage) << endl;
	
	initialize();
}

/**
@brief The copy constructor.

	The internal copy function is called.

**/
Bimage::Bimage(const Bimage& p)
{
//	cout << "Copy image" << endl;
	internal_copy(p);
}

/**
@brief Reads an image from a file.

	The internal initialization function is called.

**/
Bimage::Bimage(Bstring& fn, int readdata, int img_select)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage: Bimage size = " << sizeof(Bimage) << endl;
	
	initialize();

}

/**
@brief Initializes an image.
@param 	type	data type.
@param 	ctype	compound type.
@param 	nx		x dimension.
@param 	ny		y dimension.
@param 	nz		z dimension.
@param 	nn		number of sub-images.

	The internal initialization function is called.
	The basic properties and size is set up and the data allocated.

**/
Bimage::Bimage(DataType type, CompoundType ctype, long nx, long ny,
				long nz, long nn)
{
	initialize(ctype, nx, ny, nz, nn);

	data_type(type);

	data_alloc_and_clear();
}

/**
@brief Initializes an image.
@param 	type	data type.
@param 	ctype	compound type.
@param 	size	image size.
@param 	nn		number of sub-images.

	The internal initialization function is called.
	The basic properties and size is set up and the data allocated.

**/
Bimage::Bimage(DataType type, CompoundType ctype, Vector3<long> size, long nn)
{
	initialize(ctype, size[0], size[1], size[2], nn);

	data_type(type);
	
	data_alloc_and_clear();
}

/**
 @brief Initializes an image.
 @param 	type	data type.
 @param 	ctype	compound type.
 @param 	size	image size.
 @param 	nn		number of sub-images.
 
 The internal initialization function is called.
 The basic properties and size is set up and the data allocated.
 
 **/
Bimage::Bimage(DataType type, CompoundType ctype, vector<long> size, long nn)
{
	initialize(ctype, size[0], size[1], size[2], nn);
	
	data_type(type);
	
	data_alloc_and_clear();
}

/**
@brief Initializes a generic multi-channel image.
@param 	type	data type.
@param 	nc		number of channels.
@param 	nx		x dimension.
@param 	ny		y dimension.
@param 	nz		z dimension.
@param 	nn		number of sub-images.

	The internal initialization function is called.
	The basic properties and size is set up and the data allocated.

**/
Bimage::Bimage(DataType type, long nc, long nx, long ny,
				long nz, long nn)
{
	initialize(nc, nx, ny, nz, nn);
	
	data_type(type);
	
	data_alloc_and_clear();
}

/**
@brief Initializes a generic multi-channel image.
@param 	type	data type.
@param 	nc		number of channels.
@param 	size	image size.
@param 	nn		number of sub-images.

	The internal initialization function is called.
	The basic properties and size is set up and the data allocated.

**/
Bimage::Bimage(DataType type, long nc, Vector3<long> size, long nn)
{
	initialize(nc, size[0], size[1], size[2], nn);
	
	data_type(type);
	
	data_alloc_and_clear();
}

/**
@brief 	Initializes an image from a 2D matrix.
@param 	mat			matrix.
@param 	scale		integer scale for enlarging the image.

	The internal initialization function is called.
	The basic properties and size is set up and the data allocated.
	The matrix contents are transferred with dimension scaling.

**/
Bimage::Bimage(Matrix& mat, long scale)
{
	if ( scale < 1 ) scale = 1;
	if ( scale > 10 ) scale = 10;

	initialize(1, scale*mat.columns(), scale*mat.rows(), 1, 1);

	data_type(Float);

	data_alloc_and_clear();
		
	long		i, j, k, xx, yy;
		
	for ( k=yy=0; yy<y; ++yy ) {
		i = yy/scale;
		for ( xx=0; xx<x; ++xx, ++k ) {
			j = xx/scale;
			set(k, mat[i][j]);
		}
	}
}

/**
@brief Image destructor.

	All properties and associated data are deallocated.

**/
Bimage::~Bimage()
{
	if ( next ) {
		try { delete next; }
		catch (...) { cerr << "Failed to deallocate next image!" << endl; }
	}
	next = NULL;
	
	// Strings
	id = 0;
	
//	metadata.erase();
	
	if ( image ) {
		try { delete[] image; }
		catch (...) { cerr << "Failed to deallocate sub-images!" << endl; }
	}
	image = NULL;
	
	if ( d.uc ) {
		try { delete[] d.uc; }
		catch (...) { cerr << "Failed to deallocate data block!" << endl; }
	}
	d.uc = NULL;
}

/*
@brief Image initialization function.

	All properties are set to default values.

**/
void		Bimage::initialize()
{
	// Set parameter defaults	
	x = y = z = c = n = 0;
	px = py = pz = 0;
	sn = sz = 0;
	ss = 1;

	data_type(UCharacter);
	compoundtype = TSimple;
	fouriertype = NoTransform;
	
	offset = 0;
	
	min = max = avg = std = 0;
	smin = smax = 0;

	ucell = UnitCell();

	// Strings
	id = 0;

	// Metadata
	metadata = JSvalue(JSobject);
	next = NULL;
	image = NULL;
	d.uc = NULL;
	
	data_size();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::initialize: Image setup done" << endl;
}

void		Bimage::initialize(long nc, long nx,
						long ny, long nz, long nn)
{
	initialize();

	c = nc;
	x = nx;
	y = ny;
	z = nz;
	
	if ( c > 1 ) compoundtype = TMulti;
	
	images(nn);
}

void		Bimage::initialize(CompoundType ctype, long nx,
						long ny, long nz, long nn)
{
	long		nc;
	
	switch ( ctype ) {
		case TSimple: nc = 1; break;
		case TComplex: nc = 2; break;
		case TVector3: nc = 3; break;
		case TView: nc = 4; break;
		case TRGB: nc = 3; break;
		case TRGBA: nc = 4; break;
		case TCMYK: nc = 4; break;
		default: nc = 1;
	}

	initialize(nc, nx, ny, nz, nn);

	compoundtype = ctype;
}

CompoundType	Bimage::guess_compoundtype(long nc)
{
	switch ( nc ) {
		case 1: return TSimple;
		case 2: return TComplex;
		case 3: return TRGB;
		case 4: return TRGBA;
	}
	return TMulti;
}

bool		Bimage::check_compoundtype(long nc, CompoundType ct)
{
	switch ( ct ) {
		case TSimple: return nc == 1;
		case TComplex: return nc == 2;
		case TVector2: return nc == 2;
		case TVector3: return nc == 3;
		case TView: return nc == 4;
		case TRGB: return nc == 3;
		case TRGBA: return nc == 4;
		case TCMYK: return nc == 4;
		case TMulti: return 1;
	}
	return 1;
}

long		Bimage::compound_type(CompoundType ct)
{
	compoundtype = ct;
	
	switch ( ct ) {
		case TSimple: c = 1; break;
		case TComplex: c = 2; break;
		case TVector2: c = 2; break;
		case TVector3: c = 3; break;
		case TView: c = 4; break;
		case TRGB: c = 3; break;
		case TRGBA: c = 4; break;
		case TCMYK: c = 4; break;
		case TMulti: ; break;
	}
	return c;
}

void		Bimage::images(long nn) {
	if ( nn <= 0 )
		cerr << "Error in Bimage::images: The number of images must be greater than 1" << endl;
	if ( nn != n && image ) { delete[] image; image = NULL; }
	n = nn;
	if ( !image ) {
		image = new Bsub_image[nn];
		meta_data_update();
	}
}

void		Bimage::internal_copy(const Bimage& p)
{
	offset = p.offset;
	
	data_type(p.datatype);
	compoundtype = p.compoundtype;
	fouriertype = p.fouriertype;
	c = p.c;
	x = p.x;
	y = p.y;
	z = p.z;
//	n = p.n;
	px = p.px;
	py = p.py;
	pz = p.pz;
	
	min = p.min;
	max = p.max;
	avg = p.avg;
	std = p.std;

	ucell = p.ucell;

	// Strings
	id = p.id;
	
	next = NULL;
	
	images(p.n);
//	cout << "image pointer:" << image << "." << endl;
	if ( p.image )
		for ( long j=0; j<n; j++ ) image[j] = p.image[j];

	// Metadata
	metadata = p.metadata;

	data_size();
	long	allocsize = alloc_size();

	if ( d.uc ) delete[] d.uc;
	d.uc = NULL;
	if ( p.d.uc ) {
		d.uc = new unsigned char[allocsize];
		for ( long j=0; j<datasize; j++ ) set(j, p[j]);
	}
//	cout << "data pointer:" << &d.uc << "." << endl;
	
}

/**
@brief Checks an image for reasonable properties.

**/
void			Bimage::check()
{
	long			j;
	
	check_sampling();
	
	Vector3<double>	u(image->sampling());
	
	if ( fouriertype > Standard ) fouriertype = Standard;
	
	if ( c < 1 ) c = 1;
	if ( px < 1 ) px = x;
	if ( py < 1 ) py = y;
	if ( pz < 1 ) pz = z;
	if ( px > 2*x ) px = x;
	if ( py > 2*y ) py = y;
	if ( pz > 2*z ) pz = z;
	if ( n < 1 ) n = 1;
//	if ( ux <= 1e-6 ) ux = 1;
//	if ( uy <= 1e-6 ) uy = 1;
//	if ( uz <= 1e-6 ) uz = 1;
//	if ( ux > 1e6 ) ux = 1;
//	if ( uy > 1e6 ) uy = 1;
//	if ( uz > 1e6 ) uz = 1;
//	if ( fabs(ux - uy) > 100 ) ux = uy = uz = 1;
	if ( !compoundtype ) compoundtype = TSimple;
	
	switch ( compoundtype ) {
		case TSimple:	j = 1;
		case TComplex:	j = 2;
		case TVector3:	j = 3;
		case TView:		j = 4;
		case TRGB:		j = 3;
		case TRGBA:		j = 4;
		case TCMYK:		j = 4;
		default: j = c;
	}
	
	if ( j != c ) {
		cerr << "Error in Bimage::check: The compound type size (" << j << ") and number of channels (" << c << ") don't agree!" << endl;
		bexit(-1);
	}
	
	data_size();
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::check: Data type=" << datatype << " Compound type=" << compoundtype << endl;
		cout << "DEBUG Bimage::check: Images & channels: " << n << " " << c << endl;
		cout << "DEBUG Bimage::check: Size: " << x << " " << y << " " << z << endl;
		cout << "DEBUG Bimage::check: Page: " << px << " " << py << " " << pz << endl;
		cout << "DEBUG Bimage::check: Units: " << u[0] << " " << u[1] << " " << u[2] << endl;
		cout << "DEBUG Bimage::check: min=" << min << " max=" << max << " avg=" << avg << " std=" << std << endl;
	}

	if ( !image ) {
		cerr << "Error in Bimage::check: The sub-image list has not been allocated!" << endl;
		bexit(-1);
	}

	int			check_stats = 0;
	
	if ( !isfinite(min) || !isfinite(max) || !isfinite(avg) || !isfinite(std) ) check_stats = 1;

	if ( (min>=max) || (std<=0) || (avg<min) || (avg>max) ) check_stats = 1;
	
	if ( min < dtmin ) check_stats = 1;
	if ( max > dtmax ) check_stats = 1;
	
	if ( check_stats ) statistics();
	
	if ( smin == 0 && smax == 0 ) {
		smin = min;
		smax = max;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::check: Min, max, mean, stdev: " << 
				min << " " << max << " " << avg << " " << std << endl;
		cout << "DEBUG Bimage::check: Background:";
//		for ( j=0; j<n; j++ ) cout << " " << image[j].background();
//		cout << endl;
		for ( j=0; j<n; j++ ) cout << " " << background(j);
		cout << endl;
	}
	
	View		view;
	for ( j=0; j<n; j++ ) {
		if ( !isfinite(image[j].background()) ) image[j].background(0);
		else if ( fabs((double)image[j].background()) < 1e-30 ) image[j].background(0);
		if ( !isfinite(image[j].origin()[0]) || !isfinite(image[j].origin()[0])
			| !isfinite(image[j].origin()[0]) ) image[j].origin(0,0,0);
		view = image[j].view();
		view.normalize();
		image[j].view(view);
		if ( image[j].magnification() <= 0 )
			image[j].magnification(1);
	}

	UnitCell		uc = unit_cell();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::check: Unit cell: " << uc.a() << " " << uc.b() << " " << uc.c()
			 << " " << uc.alpha() << " " << uc.beta() << " " << uc.gamma() << endl;
	
	
	if ( uc.a() < 1 ) uc.a(u[0]*x);
	if ( uc.b() < 1 ) uc.b(u[1]*y);
	if ( uc.c() < 1 ) uc.c(u[2]*z);
	if ( uc.a() > 100*x*u[0] ) uc.a(u[0]*x);
	if ( uc.b() > 100*y*u[1] ) uc.b(u[1]*y);
	if ( uc.c() > 100*z*u[2] ) uc.c(u[2]*z);
	if ( uc.alpha() > TWOPI || uc.beta() > TWOPI || uc.gamma() > TWOPI ) uc.degrees_to_radians();
	if ( uc.alpha() < 0.001 || uc.alpha() > M_PI ) uc.alpha(M_PI_2);
	if ( uc.beta() < 0.001 || uc.beta() > M_PI ) uc.beta(M_PI_2);
	if ( uc.gamma() < 0.001 || uc.gamma() > M_PI ) uc.gamma(M_PI_2);
	
	unit_cell(uc);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::check: Unit cell: " << uc.a() << " " << uc.b() << " " << uc.c()
			 << " " << uc.alpha() << " " << uc.beta() << " " << uc.gamma() << endl;
/*	
	Bstring		t;
	if ( lbl.length() ) {
		j = lbl.length() - 1;
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::check: lbl length=" << j << " last=" << lbl[j] << endl;
//		if ( verbose & VERB_DEBUG )
//			cout << "DEBUG Bimage::check: lbl=" << lbl << "/end/" << endl;
		while ( j > 0 && isspace(lbl[j]) ) j--;
		j++;
		if ( j < lbl.length() ) lbl = lbl.left(j);
	}
*/	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::check: Finished checking parameters" << endl;
}

/**
@brief Check if two images are the same size.
@param 	*p			image to check.
@return bool 		1 if yes, 0 if no.

	The function returns the answer and leaves the calling function to 
	deal with the result.

**/
bool 		Bimage::check_if_same_size(Bimage* p)
{
	if ( !p ) {
		error_show("Both images must be defined", __FILE__, __LINE__);
		return 0;
	}
	
	if ( x != p->x || y != p->y || z != p->z || c != p->c || n != p->n ) {
		cerr << "Error: Image dimensions not the same:" << endl;
		cerr << "Image 1: c=" << c << " x=" << x << " y=" << y << " z=" << z << " n=" << n << endl;
		cerr << "Image 2: c=" << p->c << " x=" << p->x << " y=" << p->y << " z=" << p->z << " n=" << p->n << endl;
		return 0;
	}
    
	return 1;
}

/**
@brief Check if two images are the same size.
@param 	*p			image to check.
@return bool 		1 if yes, 0 if no.

	The function returns the answer and leaves the calling function to 
	deal with the result.

**/
bool 		Bimage::check_if_same_image_size(Bimage* p)
{
	if ( !p ) {
		error_show("Both images must be defined", __FILE__, __LINE__);
		return 0;
	}
	
	if ( x != p->x || y != p->y || z != p->z || c != p->c ) {
		cerr << "Error: Image dimensions not the same:" << endl;
		cerr << "Image 1: c=" << c << " x=" << x << " y=" << y << " z=" << z << " n=" << n << endl;
		cerr << "Image 2: c=" << p->c << " x=" << p->x << " y=" << p->y << " z=" << p->z << " n=" << p->n << endl;
		return 0;
	}
    
	return 1;
}

/**
@brief Checks that the sampling is properly specified.

	If the sampling is zero or less it is set to 1.
	If a dimension size is one, the sampling is set to 1.

**/
void		Bimage::check_sampling()
{
	Vector3<double>	u;

	for ( long nn=0; nn<n; ++nn ) {
		u = image[nn].sampling();
		for ( long i=0; i<3; ++i )
			if ( size()[i] < 2 ) u[i] = 1;
		image[nn].sampling(u);
	}	
}

/**
@brief Checks that the resolution falls within reasonable limits.
@param &resolution	resolution (modified).

	The resolution is set to between the real size and the sampling.

**/
void		Bimage::check_resolution(double& resolution)
{
	if ( resolution > real_size()[0] ) {
		resolution = real_size()[0];
		return;
	}

	double			invsum2(0), rm;
	Vector3<double>	u(image->sampling());
	
	if ( x > 1 ) invsum2 += 1/(u[0]*u[0]);
	if ( y > 1 ) invsum2 += 1/(u[1]*u[1]);
	if ( z > 1 ) invsum2 += 1/(u[2]*u[2]);
	rm = 2/sqrt(invsum2);
	if ( resolution < rm ) resolution = rm;
}
/*
void		Bimage::check_resolution(double& resolution)
{
	Vector3<double>	u(image->sampling());
	if ( x > 1 && resolution < u[0] ) resolution = u[0];
	if ( y > 1 && resolution < u[1] ) resolution = u[1];
	if ( z > 1 && resolution < u[2] ) resolution = u[2];
	if ( resolution > real_size()[0] ) resolution = real_size()[0];
}
*/

/**
@brief Check if this image has the same number of channels and data type as another.
@param 	*p			reference image.
@return bool 		1 if yes, 0 if no.

	The function returns the answer and leaves the calling function to 
	deal with the result.

**/
bool 		Bimage::compatible(Bimage* p)
{
	int 	non(0);

	if ( c != p->c ) {
		cout << metadata["filename"] << " has an incompatible number of channels:  " << c << endl;
		non ++;
	}
	
	if ( datatype != p->datatype ) {
		cout << metadata["filename"] << " has an incompatible data type:   " << datatype << endl;
		non ++;
	}
	
	return (non==0);
}

/**
@brief	Update metadata from the sub-image information.

	The sub-image information is encoded in the metadata in JSON format.

**/
void		Bimage::meta_data_update()
{
	if ( n <= 0 ) {
		cerr << "Error in Bimage::meta_data_update: The number of images must be greater than 0" << endl;
		return;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::meta_data_update: Updating metadata" << endl;

	metadata["version"] = BVERSION; 		// Version number
	metadata["time"] = time(NULL);
	metadata["user"] = get_user_name();

	JSvalue					arr(JSarray);

	if ( metadata.exists("image") )
		arr = metadata["image"].array();
	
	for ( long i=0; i<n; ++i ) {
		Bsub_image&		si = image[i];
		if ( i >= arr.size() )
			arr.push_back(JSvalue(JSobject));
		JSvalue&		img = arr[i];
		img["number"] = i+1;
		img["sampling"] = si.sampling().array();
		img["origin"] = si.origin().array();
		img["view"] = si.view().array();
//		img["magnification"] = si.magnification();
		img["background"] = si.background();
		img["average"] = si.average();
		img["standard_deviation"] = si.standard_deviation();
		img["minimum"] = si.minimum();
		img["maximum"] = si.maximum();
		img["fom"] = si.FOM();
		img["select"] = si.select();
	}

	metadata["image"] = arr;
}

/**
@brief	Update sub-image information from the metadata.

	The metadata in JSON format is transferred to the sub-image information.

**/
void		Bimage::update_from_meta_data()
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG update_from_meta_data: converting JSON header" << endl;
		
	vector<double>			sam(3), ori(3), vw(4);
	JSvalue					arr(JSarray);
	
	if ( metadata.exists("image") )
		arr = metadata["image"].array();
	
	for ( auto it = arr.begin(); it != arr.end(); ++it ) {
		JSvalue&		img = *it;
		long			i = img["number"].integer() - 1;
		if ( i >= n ) break;
		Bsub_image&		si = image[i];
		if ( img.exists("sampling") ) {
			sam = img["sampling"].array_real();
			si.sampling(sam);
		}
		if ( img.exists("origin") ) {
			ori = img["origin"].array_real();
			si.origin(ori);
		}
		if ( img.exists("view") ) {
			vw = img["view"].array_real();
			si.view(vw);
		}
		if ( img.exists("magnification") ) si.magnification(img["magnification"].real());
		if ( img.exists("background") ) si.background(img["background"].real());
		if ( img.exists("average") ) si.average(img["average"].real());
		if ( img.exists("standard_deviation") ) si.standard_deviation(img["standard_deviation"].real());
		if ( img.exists("minimum") ) si.minimum(img["minimum"].real());
		if ( img.exists("maximum") ) si.maximum(img["maximum"].real());
		if ( img.exists("fom") ) si.FOM(img["fom"].real());
		if ( img.exists("select") ) si.select(img["select"].integer());
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG update_from_meta_data: converting done" << endl;
}

/**
@brief	Erases all sub-image records from the meta data except one.
@param 	img_num			sub-image number.
**/
void		Bimage::meta_data_retain_one_image(long img_num)
{
	JSvalue					arr(JSarray);
	
	if ( metadata.exists("image") )
		arr = metadata["image"].array();
	else return;
	
	JSvalue					img = arr[img_num];
	JSvalue					arrnu(JSarray);

	img["number"] = 1;

	arrnu.push_back(img);

	metadata["image"] = arrnu;

//	cout << "image number = " << img_num << tab << "new array = " << metadata["image"] << endl;
}

/**
@brief 	Allocate image data.
@param 	nbytes			number of bytes to allocate.
@return unsigned char*	pointer to data.

	Any old data block is deleted.

**/
unsigned char*	Bimage::data_alloc(long nbytes)
{
	if ( !data_size() ) {
		cerr << "Error in Bimage::data_alloc: Data size is zero!" << endl;
		cerr << "(c=" << c << " x=" << x << " y=" << y << " z=" << z << " n=" << n << ")" << endl;
		exit(-1);
	}
	
	long			ds(datasize*data_type_size());
	if ( datatype == Bit ) ds /= 8;
	
	if ( nbytes != ds ) {
		cerr << "Error in Bimage::data_alloc: Data size is not equal to the allocation requested!" << endl;
		cerr << "(" << ds << " != " << nbytes << ")" << endl;
		exit(-1);
	}
	
	if ( d.uc ) delete [] d.uc;
	
	d.uc = new unsigned char[nbytes];
	
	if ( !d.uc ) {
		cerr << "Error: Image data not allocated! (" << nbytes << ")" << endl;
		exit(-1);
	}
	
//	for ( long j=0; j<nbytes; j++ ) d.uc[j] = 0;

	data_size();

	return d.uc;
}

/**
@brief Allocate image data with size parameters.
@param 	type	data type.
@param 	ctype	compound type.
@param 	nx		x dimension.
@param 	ny		y dimension.
@param 	nz		z dimension.
@param 	nn		number of sub-images.
@return unsigned char*	pointer to data.

	The image parameters are set and a new data block allocated.

**/
unsigned char*	Bimage::data_alloc(DataType type, CompoundType ctype, long nx, 
							   long ny, long nz, long nn)
{
	if ( image ) {
		if ( n != nn ) {
			delete[] image;
			image = NULL;
		}
	}

	data_type(type);
	
	c = compound_type_size();
	x = nx; y = ny; z = nz;
	
	if ( datatype == Bit ) px = 8*((x-1)/8 + 1);
	
	images(nn);
	
	return data_alloc();
}

/**
@brief Assign image data.
@param 	*nudata			allocated data pointer.
@return unsigned char*	pointer to data.

	Note: The size of the data block must correspond to the image dimensions!

**/
unsigned char*	Bimage::data_assign(unsigned char* nudata)
{
	if ( !nudata ) {
		cerr << "Error: Image data not assigned!" << endl;
		exit(-1);
	}
	
	if ( d.uc ) delete [] d.uc;
	
	d.uc = nudata;
	
	data_size();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::data_assign: datasize=" << datasize << endl;
	
	return d.uc;
}

/**
@brief Deallocates the image data.
**/
void		Bimage::data_delete()
{
	if ( d.uc ) delete [] d.uc;

	d.uc = NULL;
}


/**
@brief Assigns an image.

	The internal copy function is called.

**/
Bimage&		Bimage::operator=(const Bimage& p)
{
//	cout << "Assign image" << endl;
	internal_copy(p);

//	cout << "Copy image done" << endl;
	return *this;
}

/**
@brief Returns the data value at the given index.
@param 	j		index.
@return double	value.

	The elemental data value is returned in double precision.

**/
double		Bimage::operator[](long j) const
{
	if ( j >= datasize ) return 0;
	
	switch ( datatype ) {
		case Bit:		return ((d.uc[j/8] << (j%8)) & 0x80) == 0x80;
		case UCharacter:		return d.uc[j];
		case SCharacter:		return d.sc[j];
		case UShort:	return d.us[j];
		case Short:		return d.ss[j];
		case UInteger:	return d.ui[j];
		case Integer:	return d.si[j];
		case ULong:		return d.ul[j];
		case Long:		return d.sl[j];
		case Float:		return d.f[j];
		case Double:	return d.d[j];
		default:		return 0;
	}
}

/**
@brief Returns the data value at the given coordinates.
@param 	nn		image index.
@param 	xx		x coordinate.
@param 	yy		y coordinate.
@param 	zz		z coordinate.
@param 	cc		channel index.
@return double	value.

	The elemental data value is returned in double precision.

**/
double		Bimage::get(long nn, long xx, long yy, long zz, long cc)
{
	if ( nn >= n ) return 0;
	if ( xx >= x ) return 0;
	if ( yy >= y ) return 0;
	if ( zz >= z ) return 0;

	long		j = (((nn*z + zz)*y + yy)*x + xx)*c + cc;
	
	return (*this)[j];
}

/**
@brief 	Returns the data value at the given coordinates.
@param 	nn		image index.
@param 	vox		voxel coordinates.
@param 	cc		channel index.
@return double	value.

	The elemental data value is returned in double precision.

**/
double		Bimage::get(long nn, Vector3<double> vox, long cc)
{
	if ( nn >= n ) return 0;
	if ( vox[0] < 0 || vox[0] > x-1 ) return 0;
	if ( vox[1] < 0 || vox[1] > y-1 ) return 0;
	if ( vox[2] < 0 || vox[2] > z-1 ) return 0;

	long		j = (((nn*z + (long)(vox[2]+0.5))*y + (long)(vox[1]+0.5))*x + (long)(vox[0]+0.5))*c + cc;
	
	return (*this)[j];
}

/**
@brief 	Returns an array with all channel data at the given coordinates.
@param 	nn				image index.
@param 	vox				voxel coordinates.
@return vector<double>	value array.

	The values of all channels are returned in double precision.

**/
vector<double>	Bimage::values(long nn, Vector3<double> vox)
{
	long			i;
	vector<double>	values(c,0);
	
	if ( nn >= n ) return values;
	if ( vox[0] < 0 || vox[0] > x-1 ) return values;
	if ( vox[1] < 0 || vox[1] > y-1 ) return values;
	if ( vox[2] < 0 || vox[2] > z-1 ) return values;

	long		j = (((nn*z + (long)(vox[2]+0.5))*y + (long)(vox[1]+0.5))*x + (long)(vox[0]+0.5))*c;
	for ( i=0; i<c; i++, j++ ) values[i] = (*this)[j];
	
	return values;
}

/**
@brief Returns a complex value at the given index.
@param 	j	index.
@return Complex<double>	the complex floating point value.

	The index refers to compound values.

**/
Complex<double>	Bimage::complex(long j)
{
	long		i2 = 2*j;
	return Complex<double>((*this)[i2],(*this)[i2+1]);
}

/**
@brief Returns a 3-value vector at the given index.
@param 	j	index.
@return Vector3<double>	the 3-value vector.

	The index refers to compound values.

**/
Vector3<double>	Bimage::vector3(long j)
{
	long		i3 = 3*j;
	return Vector3<double>((*this)[i3],(*this)[i3+1],(*this)[i3+2]);
}

/**
@brief Returns a color value at the given index.
@param 	j	index.
@return RGB<double>		the color value.

	The index refers to compound values.

**/
RGB<double>		Bimage::rgb(long j)
{
	long		i3 = 3*j;
	return RGB<double>((*this)[i3],(*this)[i3+1],(*this)[i3+2]);
}

/**
@brief Returns a color value at the given index.

	The index refers to compound values.

@param 	j	index.
@return RGBA<double>	the color value.
**/
RGBA<double>	Bimage::rgba(long j)
{
	long		i4 = 4*j;
	return RGBA<double>((*this)[i4],(*this)[i4+1],(*this)[i4+2],(*this)[i4+3]);
}

/**
@brief Returns a color value at the given index.
@param 	j				index.
@return CMYK<double>	the color value.

	The index refers to compound values.

**/
CMYK<double>	Bimage::cmyk(long j)
{
	long		i4 = 4*j;
	return CMYK<double>((*this)[i4],(*this)[i4+1],(*this)[i4+2],(*this)[i4+3]);
}

TypePointer		Bimage::fill_value(double v)
{
	TypePointer			fp;
	fp.uc = new unsigned char[data_type_size()];
	
	switch ( datatype ) {
		case Bit:			fp.uc[0] = (v)? 1: 0; break;
		case UCharacter:	fp.uc[0] = (unsigned char)v; break;
		case SCharacter:	fp.sc[0] = (signed char)v; break;
		case UShort:		fp.us[0] = (unsigned short)v; break;
		case Short:			fp.ss[0] = (short)v; break;
		case UInteger:		fp.ui[0] = (unsigned int)v; break;
		case Integer:		fp.si[0] = (int)v; break;
		case ULong:			fp.ul[0] = (unsigned long)v; break;
		case Long:			fp.sl[0] = (long)v; break;
		case Float:			fp.f[0] = v; break;
		case Double:		fp.d[0] = v; break;
		default: ;
	}
	
	return fp;
}

/**
@brief Sets a single value at the given index.
@param 	j		index.
@param 	v		value.

**/
void		Bimage::set(long j, double v)
{
	if ( j >= datasize ) return;
	
	long		k;

	if ( v < dtmin ) v = dtmin;
	else if ( v > dtmax ) v = dtmax;
	
	switch ( datatype ) {
		case Bit:		k = (j/x)*px + j%x; if ( 2*v > min+max ) d.uc[k/8] |= 0x80 >> k%8; break;
		case UCharacter:	d.uc[j] = (unsigned char)v; break;
		case SCharacter:	d.sc[j] = (signed char)v; break;
		case UShort:	d.us[j] = (unsigned short)v; break;
		case Short:		d.ss[j] = (short)v; break;
		case UInteger:	d.ui[j] = (unsigned int)v; break;
		case Integer:	d.si[j] = (int)v; break;
		case ULong:		d.ul[j] = (unsigned long)v; break;
		case Long:		d.sl[j] = (long)v; break;
		case Float:		d.f[j] = v; break;
		case Double:	d.d[j] = v; break;
		default: ;
	}
}

/**
@brief Sets a complex value at the given index.
@param 	j		index.
@param 	cv		complex value.

	The index refers to compound values.

**/
void		Bimage::set(long j, Complex<double> cv)
{
	long		i2 = 2*j;
	
	set(i2++, cv.real());
	set(i2, cv.imag());
}

/**
@brief Sets a color value at the given index.
@param 	j		index.
@param 	color	color value.

	The index refers to compound values.

**/
void		Bimage::set(long j, RGB<double> color)
{
	long		i3 = 3*j;

	if ( i3 >= datasize ) return;
	
	set(i3++, color.r());
	set(i3++, color.g());
	set(i3, color.b());
}

/**
@brief Sets a color value at the given index.
@param 	j		index.
@param 	color	color value.

	The index refers to compound values.

**/
void		Bimage::set(long j, RGBA<double> color)
{
	long		i4 = 4*j;

	if ( i4 >= datasize ) return;
	
	set(i4++, color.r());
	set(i4++, color.g());
	set(i4++, color.b());
	set(i4, color.a());
}

/**
@brief Sets a color value at the given index.
@param 	j		index.
@param 	color	color value.

	The index refers to compound values.

**/
void		Bimage::set(long j, CMYK<double> color)
{
	long		i4 = 4*j;

	if ( i4 >= datasize ) return;
	
	set(i4++, color.c());
	set(i4++, color.m());
	set(i4++, color[1]);
	set(i4, color.k());
}

/**
@brief Sets a 3-value vector at the given index.
@param 	j		index.
@param 	vec		3-value vector.

	The index refers to compound values.

**/
void		Bimage::set(long j, Vector3<double> vec)
{
	long		i3 = 3*j;

	if ( i3 >= datasize ) return;
	
	set(i3++, vec[0]);
	set(i3++, vec[1]);
	set(i3, vec[2]);
}

/**
@brief Sets a 3-value vector at the given index.
@param 	j			index.
@param 	view		4-value view.

	The index refers to compound values.
**/
void		Bimage::set(long j, View view)
{
	long		i4 = 4*j;

	if ( i4 >= datasize ) return;
	
	set(i4++, view[0]);
	set(i4++, view[1]);
	set(i4++, view[2]);
	set(i4, view.angle());
}

/**
@brief Adds a value at a given location to neigboring data elements.
@param 	xx			x location.
@param 	yy			y location.
@param 	zz			z location.
@param 	nn			image number (4 th dimension).
@param 	v			value to add.

	Inverse of trilinear interpolation.
**/
void		Bimage::add(double xx, double yy, double zz, long nn, double v)
{
	long			slicesize = x*y;
	long			imgsize = slicesize*z;
	long			ix = (long) xx;
	long			iy = (long) yy;
	long			iz = (long) zz;
	long			jx, jy, jz, xk, yk, zk;
	double			fx = xx - ix;
	double			fy = yy - iy;
	double			fz = zz - iz;
	double			w;
	
	for ( zk=0; zk<2; zk++, iz++ ) {
		fz = 1 - fz;
		jz = nn*imgsize + iz*slicesize + iy*x;
		for ( yk=0, jy=jz; yk<2; yk++, jy+=x ) {
			fy = 1 - fy;
			for ( xk=0, jx=jy+ix; xk<2; xk++, jx++ ) {
				fx = 1 - fx;
				w = fx*fy*fz;
				add(jx, v * w);
			}
		}
	}
}

/**
@brief Sets a value at a given location to neigboring data elements if it is larger.
@param 	xx			x location.
@param 	yy			y location.
@param 	zz			z location.
@param 	nn			image number (4 th dimension).
@param 	v			value to assess and set.

	All neighboring voxels are checked and set to the given value if it is larger.
**/
void		Bimage::set_max(double xx, double yy, double zz, long nn, double v)
{
	long	slicesize = x*y;
	long	imgsize = slicesize*z;
	long	ix = (long) xx;
	long	iy = (long) yy;
	long	iz = (long) zz;
	long	jx, jy, jz, xk, yk, zk;
	
	for ( zk=0; zk<2; zk++, iz++ ) {
		jz = nn*imgsize + iz*slicesize + iy*x;
		for ( yk=0, jy=jz; yk<2; yk++, jy+=x ) {
			for ( xk=0, jx=jy+ix; xk<2; xk++, jx++ ) {
				if ( (*this)[jx] < v ) set(jx, v);
			}
		}
	}
}

/**
@brief Averages in the xy plane when scaling is less than 1.
@param 	cc			channel.
@param 	xf			x location.
@param 	yf			y location.
@param 	zf			z location.
@param 	nn			image number.
@param 	iscale		inverse scale.
@return double		average.

	All neighboring voxels in x and y are averaged and result returned.
**/
double		Bimage::average2D(long cc, double xf, double yf, 
							double zf, long nn, double iscale)
{
	long			zz = (long) zf;
	long			j, ix, iy, num(0), ij = (nn*z + zz)*y*x;
	double			xx, yy, avg(0);
	double			xmax = xf + iscale - 0.5, ymax = yf + iscale - 0.5;
	
	for ( yy=yf; yy < ymax; yy += 1 ) {
		iy = (long) yy;
		if ( iy < y ) for ( xx=xf; xx < xmax; xx += 1 ) {
			ix = (long) xx;
			if ( ix < x ) {
				j = (ij + iy*x + ix)*c + cc;
				avg += (*this)[j];
				num++;
			}
		}
	}
	
	if ( num ) avg /= num;
	
	return avg;
}

/**
@brief Averages in 3D when scaling is less than 1.
@param 	cc			channel.
@param 	xf			x location.
@param 	yf			y location.
@param 	zf			z location.
@param 	nn			image number.
@param 	iscale		inverse scale.
@return double		average.

	All neighboring voxels in x, y and z are averaged and result returned.
**/
double		Bimage::average(long cc, double xf, double yf, 
							double zf, long nn, double iscale)
{
	long			j, ix, iy, iz, ij, num(0);
	double			xx, yy, zz, avg(0);
	double			xmax = xf + iscale - 0.5;
	double 			ymax = yf + iscale - 0.5;
	double 			zmax = zf + iscale - 0.5;
	if ( xmax > x ) xmax = x;
	if ( ymax > y ) ymax = y;
	if ( zmax > z ) zmax = z;
	
	for ( zz=zf; zz < zmax; zz += 1 ) {
		iz = (long) zz;
		ij = (nn*z + iz)*y*x;
		for ( yy=yf; yy < ymax; yy += 1 ) {
			iy = (long) yy;
			for ( xx=xf; xx < xmax; xx += 1 ) {
				ix = (long) xx;
				j = (ij + iy*x + ix)*c + cc;
				avg += (*this)[j];
				num++;
			}
		}
	}
	
	if ( num ) avg /= num;
	
	return avg;
}

/**
@brief 	Interpolates using a given location.
@param 	cc			channel.
@param 	xx			x location.
@param 	yy			y location.
@param 	zz			z location.
@param 	nn			image number.
@param 	fill		fill value if outside boundaries.
@return double		interpolated value.

	Trilinear interpolation.
**/
double		Bimage::interpolate(long cc, double xx,
				double yy, double zz, long nn, double fill) const
{
//	cout << cc << tab << xx << tab << yy << tab << zz << tab << nn << endl;
	if ( cc >= c ) return fill;
	if ( xx < 0 || xx >= x ) return fill;
	if ( yy < 0 || yy >= y ) return fill;
	if ( zz < 0 || zz >= z ) return fill;
	if ( nn >= n ) return fill;

	long			ix = (long) xx;
	long			iy = (long) yy;
	long			iz = (long) zz;

	long			nx = (x < ix + 2)? 1: 2;
	long			ny = (y < iy + 2)? 1: 2;
	long			nz = (z < iz + 2)? 1: 2;

	long			i, xk, yk, zk;
	double			fx = xx - ix;
	double			fy = yy - iy;
	double			fz = zz - iz;
	
	double			value(0), w(0), ws(0), wyz;

	iz = (nn*z + iz)*y;
	
	for ( zk=0; zk<nz; zk++, iz+=y ) {
		fz = 1.0L - fz;
		iy = (long) yy;
		iy = (iz + iy)*x;
		for ( yk=0; yk<ny; yk++, iy+=x ) {
			fy = 1.0L - fy;
			ix = (long) xx;
			i = (iy + ix)*c + cc;
			wyz = fy*fz;
			for ( xk=0; xk<nx; xk++, i+=c ) {
				fx = 1.0L - fx;
				w = fx*wyz;
				ws += w;
				value += (*this)[i] * w;
			}
		}
	}
	
	if ( ws ) value /= ws;
	else value = fill;
	
	return value;
}
/*
double		Bimage::interpolate(long cc, double xx,
				double yy, double zz, long nn, double fill) const
{
//	cout << cc << tab << xx << tab << yy << tab << zz << tab << nn << endl;
	if ( cc >= c ) return fill;
	if ( xx < 0 || xx >= x ) return fill;
	if ( yy < 0 || yy >= y ) return fill;
	if ( zz < 0 || zz >= z ) return fill;
	if ( nn >= n ) return fill;
	
	long		ix = (long) xx;
	long		iy = (long) yy;
	long		iz = (long) zz;
	long		xc(x*c), xyc(xc*y), xyzc(xyc*z);
	double		fx = xx - ix;
	double		fy = yy - iy;
	double		fz = zz - iz;
	double		value(0), w(0), ws(0), wyz;
	
	nn *= xyzc;
	iz *= xyc;
	iy *= xc;
	ix *= c;
	long		iz1 = (iz < xyzc - xyc)? iz+xyc: iz;
	long		iy1 = (iy < xyc - xc)? iy+xc: iy;
	long		ix1 = (ix < xc - c)? ix+c: ix;
	
	for ( ; iz<iz1; iz+=xyc ) {
		fz = 1 - fz;
		for ( iy = ((long) yy)*xc; iy<iy1; iy+=xc ) {
			fy = 1 - fy;
			wyz = fy*fz;
			for ( ix = ((long) xx)*c + cc; ix<ix1; ix+=c ) {
				fx = 1 - fx;
				w = fx*wyz;
				ws += w;
				value += (*this)[nn+iz+iy+ix] * w;
			}
		}
	}

	if ( ws ) value /= ws;
	
	return value;
}
*/

/**
@brief 	Interpolates using a given location with wrapping.
@param 	cc			channel.
@param 	xx			x location.
@param 	yy			y location.
@param 	zz			z location.
@param 	nn			image number.
@return double		interpolated value.

	Trilinear interpolation with wrapping around periodic image boundaries.

**/
double		Bimage::interpolate_wrap(long cc, double xx,
				double yy, double zz, long nn) const
{
//	cout << cc << tab << xx << tab << yy << tab << zz << tab << nn << endl;
	if ( cc >= c ) return 0;
	if ( nn >= n ) return 0;
	while ( xx < 0 ) xx += x;
	while ( yy < 0 ) yy += y;
	while ( zz < 0 ) zz += z;
	while ( xx >= x ) xx -= x;
	while ( yy >= y ) yy -= y;
	while ( zz >= z ) zz -= z;

	long			nx(2), ny(2), nz(2);
	if ( x < 2 ) nx = 1;
	if ( y < 2 ) ny = 1;
	if ( z < 2 ) nz = 1;
	
	long			ix = (long) xx;
	long			iy = (long) yy;
	long			iz = (long) zz;
	long			i, in(nn*z), xk, yk, zk, jy, jz;
	double			fx = xx - ix;
	double			fy = yy - iy;
	double			fz = zz - iz;
	
	double			value(0), wyz;

	for ( zk=0; zk<nz; zk++, iz++ ) {
		if ( iz >= z ) iz -= z;
		fz = 1 - fz;
		jz = (in + iz)*y;
		iy = (long) yy;
		for ( yk=0; yk<ny; yk++, iy++ ) {
			if ( iy >= y ) iy -= y;
			fy = 1 - fy;
			jy = (jz + iy)*x;
			ix = (long) xx;
			wyz = fy*fz;
			for ( xk=0; xk<nx; xk++, ix++ ) {
				if ( ix >= x ) ix -= x;
				i = (jy + ix)*c + cc;
				fx = 1 - fx;
				value += (*this)[i] * (fx*wyz);
			}
		}
	}
	
	return value;
}

/**
@brief 	Returns the size of the datatype in bits.
@return long			data type size in bits.

	The Bit type returns 1.

**/
long		Bimage::data_type_bits() const
{
	long	bits(0);
	
	switch ( datatype ) {
		case Bit: 						bits = 1; break;
		case UCharacter: case SCharacter:	bits = 8*sizeof(char); break;
		case UShort: case Short:		bits = 8*sizeof(short); break;
		case UInteger: case Integer: 	bits = 8*sizeof(int); break;
		case ULong: case Long: 			bits = 8*sizeof(long); break;
		case Float:						bits = 8*sizeof(float); break;
		case Double:					bits = 8*sizeof(double); break;
		default: 		bits = 0;
	}
	
	return bits;
}

/**
@brief 	Returns the size of the datatype in bytes.
@return long			data type size.

	The Bit type returns 1.

**/
long		Bimage::data_type_size() const
{
	long	typesize(0);
	
	switch ( datatype ) {
		case Bit:
		case UCharacter: case SCharacter:	typesize = sizeof(char); break;
		case UShort: case Short:typesize = sizeof(short); break;
		case UInteger: case Integer: 	typesize = sizeof(int); break;
		case ULong: case Long: 	typesize = sizeof(long); break;
		case Float:				typesize = sizeof(float); break;
		case Double:			typesize = sizeof(double); break;
		default: typesize = 0;
	}
	
	return typesize;
}

/**
@brief Get the minimum of a datatype.
@return double			minimum value.

	This function is used for returning the minimum value
	a data type can hold.

**/
double		Bimage::data_type_min()
{
	double			min(0);
	
	switch ( datatype ) {
		case Bit: min = 0; break;
		case UCharacter: min = 0; break;
		case SCharacter: min = -128; break;
		case UShort: min = 0; break;
		case Short: min = SHRT_MIN; break;
		case UInteger: min = 0; break;
		case Integer: min = INT_MIN; break;
		case ULong: min = 0; break;
		case Long: min = LONG_MIN; break;
		case Float: min = -FLT_MAX; break;
		case Double: min = -DBL_MAX; break;
		default: break;
	}

	return min;
}

/**
@brief Get the maximum of a datatype.
@return double			maximum value.

	This function is used for returning the maximum value
	a data type can hold.

**/
double		Bimage::data_type_max()
{
	double			max(0);
	
	switch ( datatype ) {
		case Bit: max = 1; break;
		case UCharacter: max = 255; break;
		case SCharacter: max = 127; break;
		case UShort: max = USHRT_MAX; break;
		case Short: max = SHRT_MAX; break;
		case UInteger: max = UINT_MAX; break;
		case Integer: max = INT_MAX; break;
		case ULong: max = ULONG_MAX; break;
		case Long: max = LONG_MAX; break;
		case Float: max = FLT_MAX; break;
		case Double: max = DBL_MAX; break;
		default: break;
	}

	return max;
}

/**
@brief Get the string representation of a datatype.
@return Bstring			string representation of a datatype.
**/
Bstring		Bimage::data_type_string()
{
	Bstring		string;
	
	switch ( datatype ) {
		case Bit:               string = "bit"; break;
		case UCharacter:		string = "unsigned char"; break;
		case SCharacter:		string = "signed char"; break;
		case UShort:            string = "unsigned short"; break;
		case Short:             string = "short"; break;
		case UInteger:			string = "unsigned int"; break;
		case Integer:			string = "int"; break;
		case ULong:				string = "unsigned long"; break;
		case Long:				string = "long"; break;
		case Float: 			string = "float"; break;
		case Double: 			string = "double"; break;
		default: string = "unknown";
	}
	
	return string;
}

/**
@brief Get the string representation of a datatype.
@return Bstring			string representation of a datatype.
**/
Bstring		Bimage::compound_type_string()
{
	Bstring		string;
	
	switch ( compoundtype ) {
		case TSimple: string = "Simple"; break;
		case TComplex: string = "Complex"; break;
		case TVector3: string = "Vector3"; break;
		case TView: string = "View"; break;
		case TRGB: string = "RGB"; break;
		case TRGBA: string = "RGBA"; break;
		case TCMYK: string = "CMYK"; break;
		case TMulti: string = "Multiple"; break;
		default: string = "unknown"; break;
	}
	
	return string;
}

/**
@brief Determines the replacement data type.
@return DataType				replacement data type.

	An integer data type is switched between signed and unsigned.

**/
void		Bimage::fix_type()
{
	switch ( datatype ) {
		case UCharacter: datatype = SCharacter; break;
		case SCharacter: datatype = UCharacter; break;
		case UShort: datatype = Short; break;
		case Short: datatype = UShort; break;
		case UInteger: datatype = Integer; break;
		case Integer: datatype = UInteger; break;
		case ULong: datatype = Long; break;
		case Long: datatype = ULong; break;
		default: break;
	}

	data_type(datatype);
}

/**
@brief Get the data type indicated by a single letter code.
@param 	letter 	letter indicating data type.

	This function is used in optional command-line arguments to indicate 
	a new data type for an image.

**/
void		Bimage::change_type(char letter)
{
	DataType		nutype;
	
	switch ( letter ) {
		case '1': nutype = Bit; break;
		case 'b': nutype = UCharacter; break;
		case 'c': nutype = SCharacter; break;
		case 'u': nutype = UShort; break;
		case 's': nutype = Short; break;
		case 'j': nutype = UInteger; break;
		case 'i': nutype = Integer; break;
		case 'k': nutype = ULong; break;
		case 'l': nutype = Long; break;
		case 'f': nutype = Float; break;
		case 'd': nutype = Double; break;
		case 'S': nutype = Short; break;
		case 'I': nutype = Integer; break;
		case 'F': nutype = Float; break;
		case 'D': nutype = Double; break;
		default: nutype = Unknown_Type; break;
	}
	
	change_type(nutype);
}

/**
@brief Get the data type from a string.
@param 	*string 	string describing the data type.
**/
void		Bimage::change_type(char* string)
{
	if ( strlen(string) == 1 ) {
		change_type(string[0]);
		return;
	}

	size_t   		j;
	DataType		nutype = Unknown_Type;
		
	for ( j=0; j<strlen(string); j++ ) string[j] = tolower(string[j]);
	
	if ( strstr(string, "unsigned") ) {
		if ( strstr(string, "char") ) nutype = UCharacter;
		else if ( strstr(string, "short") ) nutype = UShort;
		else if ( strstr(string, "int") ) nutype = UInteger;
		else if ( strstr(string, "long") ) nutype = ULong;
	} else {
		if ( strstr(string, "bit") ) nutype = Bit;
		else if ( strstr(string, "byte") ) nutype = UCharacter;
		else if ( strstr(string, "char") ) nutype = SCharacter;
		else if ( strstr(string, "short") ) nutype = Short;
		else if ( strstr(string, "int") ) nutype = Integer;
		else if ( strstr(string, "long") ) nutype = Long;
		else if ( strstr(string, "float") ) nutype = Float;
		else if ( strstr(string, "double") ) nutype = Double;
	}
	
	change_type(nutype);
}

/**
@brief Change the data to the new type.
@param 	nutype 	new data type.
**/
void		Bimage::change_type(DataType nutype)
{
	if ( datatype == nutype || nutype == Unknown_Type ) return;
	
	if ( max - min < 1e-37 ) max = min + 1;
	
	long			j, k, l, m;
	long			oldtypesize = data_type_size();
	double			v(0);
//	double			omn = dtmin;
	
	TypePointer		p = d;
	
	DataType		oldtype = datatype;
	
	data_type(nutype);
	double			mn = dtmin;
	double			mx = dtmax;

	if ( datatype == Bit ) px = 8*((x-1)/8+1);

	long			typesize = data_type_size();

	long			allocsize = alloc_size();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::change_type: datatype=" << datatype << " datasize=" << datasize << endl;	
	
	double			scale(1), shift(0);
	
	// Integer switches between signed and unsigned
	if ( datatype >= UCharacter && datatype <= Long ) {
		if ( oldtype >= Float ) {		// floating point to integer
			if ( min >= 0 ) scale = mx/max;
			else {
				scale = (mx - mn)/(max - min);
				shift = mn - min*scale;
			}
			if ( compoundtype == TComplex && datatype == Short ) {
				scale = 0.5*mx/max;
				shift = 0;
			}
		} else {
			if ( min == 0 && max == 1 ) {	// masks
				scale = 1;
				shift = 0;
			} else if ( oldtypesize == typesize ) {
				if ( min >= mn && max <= mx ) {
					shift = 0;
				} else if ( datatype == oldtype + 1 ) {	// unsigned -> signed
					shift = mn;
				} else if ( datatype == oldtype - 1 ) {	// signed -> unsigned
					shift = -min;
				}
			} else if ( oldtype > datatype && ( min < mn || max > mx ) ) {	// larger to smaller integer
				scale = (mx - mn)/(max - min);
				shift = mn - min*scale;
			}
		}
	}
	
	if ( datatype == Float && oldtype == Double && ( min < mn || max > mx ) ) {
			scale = (mx - mn)/(max - min);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::change_type: scale=" << scale << " shift=" << shift << endl;

	if ( oldtypesize != typesize || oldtype == Bit || datatype == Bit ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::change_type: allocsize=" << allocsize << endl;
		d.uc = new unsigned char[allocsize];
		for ( j=0; j<allocsize; j++ ) d.uc[j] = 0;
	}
	
	if ( oldtype == Bit ) {
		for ( j=l=0; l<n*y*z; l++ )
			for ( k=0, m=l*px; k<x; k++, j++, m++ ) set(j, ((p.uc[m/8] << (m%8)) & 0x80) == 0x80);
	} else for ( j=l=0; j<datasize; j++, l++ ) {
		switch ( oldtype ) {
			case UCharacter:	v = p.uc[j]; break;
			case SCharacter:	v = p.sc[j]; break;
			case UShort:		v = p.us[j]; break;
			case Short:			v = p.ss[j]; break;
			case UInteger:		v = p.ui[j]; break;
			case Integer:		v = p.si[j]; break;
			case ULong:			v = p.ul[j]; break;
			case Long:			v = p.sl[j]; break;
			case Float:			v = p.f[j]; break;
			case Double:		v = p.d[j]; break;
			default: ;
		}
		v = scale*v + shift;
		k = j;
		if ( datatype == Bit ) {
			if ( l >= x ) {
				for ( ; l<px; l++, k++ ) ;
				l = 0;
			}
		}
		set(k, v);
	}

	if ( oldtypesize != typesize || oldtype == Bit || datatype == Bit )
		delete [] p.uc;
	
	if ( statistics() )
		cerr << tab << "in Bimage::change_type" << endl;
}

/**
@brief 	Splits the channels into individual sub-images.
@return Bimage*			new image.

	The channels are converted to successive sets of images.

**/
Bimage*		Bimage::split_channels()
{
	Bimage*			pnu = copy_header(n*c);
	pnu->compound_type(TSimple);
	pnu->channels(1);
	pnu->data_alloc();
	
	long			i, j, k, m, nn, cc, imgsize(x*y*z), offset(n*imgsize);

	for ( nn=i=m=0; nn<n; ++nn )
		for ( j=0; j<imgsize; ++j, ++m )
			for ( cc=0, k=m; cc<c; ++cc, ++i, k+=offset )
				pnu->set(k, (*this)[i]);

	pnu->statistics();
	
	return pnu;
}

/**
@brief 	Splits the channels into individual images in a linked list.
@return Bimage*			new image.

	Each channels is converted to a separate image in the list.

**/
Bimage*		Bimage::split_channels_to_images()
{
	Bimage*			pnu = copy_header(n);
	pnu->compound_type(TSimple);
	pnu->channels(1);
	pnu->data_alloc();
	
	Bimage*			p1 = pnu;
	
	long			i, j, k, nn, cc, imgsize(x*y*z);

	for ( cc=0; cc<c; ++cc ) {
		if ( cc ) p1 = p1->next = p1->copy();
		for ( nn=k=0, i=cc; nn<n; ++nn )
			for ( j=0; j<imgsize; ++j, ++k, i+=c )
				p1->set(k, (*this)[i]);
		p1->statistics();
	}
	
	return pnu;
}

/**
@brief 	Combines images as channels.
@param 	nc			number of channels to create.
@param 	ct			compound type.

	The input image is replaced with a new image.
	If the compound type is single value, the compound type
	is inferred from the number of channels.
**/
void		Bimage::combine_channels(long nc, CompoundType ct)
{
	if ( nc < 2 ) return;
	
	if ( ct == TSimple ) ct = guess_compoundtype(nc);
	
	long			nun, nn, cc, j, m, k;
	long			tsize = data_type_size();
	long			imgsize = x*y*z;
	unsigned char*	nudata = new unsigned char[alloc_size()];

	nun = n/nc;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::combine_channels: " << nc << " channels" << endl;
		cout << "DEBUG Bimage::combine_channels: " << ct << " compound type" << endl;
	}
	
	for ( nn=k=0; nn<nun; nn++ ) {
		for ( j=0; j<imgsize; j++ ) {
			for ( cc=0, m=nn*imgsize*nc+j; cc<nc; cc++, k+=tsize, m+=imgsize ) {
				memcpy(nudata+k, d.uc+m*tsize, tsize);
			}
		}
   }

	data_assign(nudata);
	
	Bsub_image*		nusub = new Bsub_image[nun];
	
	if ( image ) {
		for ( nn=j=0; nn<nun; nn++, j+=nc )
			nusub[nn] = image[j];
		delete[] image;
	}

	image = nusub;
	
	n = nun;
	c = nc;
	compoundtype = ct;

	statistics();
}


/**
@brief Zeroes the first voxel in each image.

**/
void		Bimage::zero_fourier_origin()
{
	long			j, nn;
	long			imagesize = x*y*z;
	Complex<double>	cv;
	
	for ( j=nn=0; nn<n; nn++, j+=imagesize ) set(j, cv);
	
}


/**
@brief Returns the radius of the enclosed sphere or circle.
@return double		radius of largest sphere that will fit in an image.

	Returns the radius of the largest sphere or circle that will fit
	inside a three- or two-dimensional image, respectively, or the 
	midpoint of a one-dimensional image.
	Assumes that the x, y, and z dimensions are the number of data points
	in those respective directions, and that the length (in pixels)
	of the respective dimensions is x, y, or z minus one.  The minimum 
	length is two times the radius of the largest sphere that will fit 
	inside the image.

**/
double		Bimage::maximum_included_radius()
{
	double			max_radius;
	long			min_dimension(x);
	if ( (y > 1) && (y < min_dimension) )  min_dimension = y;
	if ( (z > 1) && (z < min_dimension) )  min_dimension = z;

	max_radius = (min_dimension - 1)/2;

	return max_radius;
}

/**
@brief Changes the slices in a 3D image into a set of 2D images.
@return int					error code (<0 means failure).
**/
int 		Bimage::slices_to_images()
{
	if ( z == 1 ) return 0;
	
	if ( n > 1 ) {
		cerr << "Error: " << metadata["filename"] << " is already a multi-image file!" << endl;
		error_show("Bimage::slices_to_images", __FILE__, __LINE__);
		return -1;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::slices_to_images: " << z << " slices to images" << endl;
	
	long			nn;
	Bsub_image		subimg = image[0];

	if ( image ) delete [] image;
	
	n = z;
	z = 1;
	
	image = new Bsub_image[n];
	
	for ( nn=0; nn<n; nn++ )
		image[nn] = subimg;
	
/*	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::slices_to_images: metadata:" << endl << metadata << endl;
	
	for ( nn=1; nn<n; nn++ )
		metadata["image"].array().push_back(metadata["image"][0]);
*/	
	origin(default_origin());
	
	return 0;
}

/**
@brief Changes the 2D images to slices in a 3D image.
@return int					error code (<0 means failure).
**/
int 		Bimage::images_to_slices()
{
	if ( z > 1 ) {
		cerr << "Error: " << metadata["filename"] << " is already a 3D image file!" << endl;
		error_show("Bimage::images_to_slices", __FILE__, __LINE__);
		return -1;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::images_to_slices: " << n << " images to slices" << endl;
	
	Bsub_image		subimg = image[0];

	if ( image ) delete [] image;
	
	z = n;
	n = 1;
	
	image = new Bsub_image[n];
	
	image[0] = subimg;

	origin(default_origin());
	
	return 0;
}

/**
@brief 	Sets the sub-image selections based on a list.
@param 	list		list of sub-images to select.
@return long			number of images selected.

	The data is not altered.

**/
long	Bimage::set_subset_selection(Bstring list)
{
	long			nimg(0);
	
	vector<int>		imgnum = select_numbers(list, n);
	
	for ( long nn = 0; nn < n; ++nn ) {
		if ( imgnum[nn] ) image[nn].select(1);
		else image[nn].select(0);
		nimg += imgnum[nn];
	}
	
	return nimg;
}

/**
@brief 	Retains or deletes sub-images from a multi-image structure.
@param 	list		list of sub-images to retain or delete.
@param	retain		0=delete list, 1=retain list.
@return long	error code (<0 means failure) or number of sub-images remaining.

	Sub-images specified in a list are either retained or deleted.
	The new data replaces the old data.

**/
long	Bimage::delete_images(Bstring list, int retain)
{
	long			ni, no, ns(0);
	
	vector<int>		imgnum = select_numbers(list, n);
	
	for ( ni=0; ni<n; ni++ ) if ( imgnum[ni] ) ns++;
	if ( !retain ) ns = n - ns;

	long			imgsize = x*y*z*c*data_type_size();

	Bsub_image*		subimg = new Bsub_image[ns];
	unsigned char*	nudata = new unsigned char[ns*imgsize];
	
	for ( ni=no=0; ni<n; ni++ ) if ( retain == imgnum[ni] ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG delete_images: retaining " << ni << endl;
		memcpy(nudata+no*imgsize, d.uc+ni*imgsize, imgsize);
		subimg[no] = image[ni];
		no++;
	}

	delete[] image;
	image = subimg;
	n = ns;
	
	data_assign(nudata);

	return ns;
}

long	Bimage::select_images(Bstring list)
{
	return delete_images(list, 1);
}

/**
@brief	Read image data in a generalized style.
@param	*fimg		file descriptor: file opened in calling function.
@param 	img_select	image selection: if -1, all images, if >= 0, one image.
@param 	sb			flag activates byte swapping.
@param 	vax 		indicate vax style floating point - activates conversion.
@param 	pad 		any interspersed control or separation bytes.
@return	unsigned char*	data pointer, NULL if reading failed.
The whole file or a single image from a file may be read.
	The data is read in the largest blocks possible for efficiency.
	Any interspersed padding and page sizes not matching the data size
	contribute to inefficiency in reading.
	Swapping:	sb = 1:	swap bytes the size of the data type
				sb > 1:	swap these number of bytes regardless of the data type
**/
unsigned char*	Bimage::read_data(ifstream* fimg, int img_select, int sb, int vax, long pad)
{
	if ( page_size().volume() < 1 ) page_size(x, y, z);
	
	// If only half of a transform is stored, it need to be handled
	long				xstore = x;
	long				xpage = px;
	if ( fouriertype == Hermitian || fouriertype == CentHerm ) {
		xstore = x/2 + 1;
		if ( px > xstore ) xpage = xstore;
	}
	
	if ( datatype == Bit ) {
		xstore = xpage = (x - 1)/8 + 1;
		page_size(8*xstore, 1, 1);
	}
	
	// Reset select to get the correct offset
	if ( img_select < 0 ) img_select = 0;
	
	long			i, iz, iy, j, jz, jy, nn, xx, yy, zz, pxx, pyy, pzz, npages, off;
	long			datatypesize = data_type_size();
	long			elementsize = datatypesize*c;
	long			datasize = n*xstore*y*z*elementsize;
	long			pagesize = xpage*py*pz*elementsize;
	long			readsize, pagemax = 1073741824;
	unsigned char*	data = data_alloc(datasize);
	unsigned char*	page = NULL;

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::read_data: type size=" << datatypesize << endl;
		cout << "DEBUG Bimage::read_data: data size: " << size() << " c=" << c << " n=" << n << " (" << datasize << ")" << endl;
		cout << "DEBUG Bimage::read_data: page size: " << page_size() << " (" << pagesize << ")" << endl;
		cout << "DEBUG Bimage::read_data: sb=" << sb << " vax=" << vax << " offset=" << offset << endl;
	}
	
	if ( !data ) { 
		error_show("No image data block allocated!", __FILE__, __LINE__);
		return NULL; 
	}

	if ( x == px && y == py && z == pz ) {		// 3D block
		off = offset + img_select*(pagesize + pad);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::read_data: Reading 3D blocks: off=" << off << " pad=" << pad << endl;
		fimg->seekg(off);
		for ( nn=0; nn<n; nn++ ) {
			for ( j=0; j<pagesize; j+=pagemax ) {
				readsize = pagesize - j;
				if ( readsize > pagemax ) readsize = pagemax;
				fimg->read((char *) (data + nn*pagesize + j), readsize);
			}
			if ( pad ) fimg->seekg(pad, ios_base::cur);
		}
	} else if ( x == px && y == py && pz == 1 ) {	// 2D page (MFF)
		off = offset + img_select*z*(pagesize + pad);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::read_data: Reading 2D pages: off=" << off << " pad=" << pad << endl;
		fimg->seekg(off);
		for ( nn=0; nn<n; nn++ ) {
			for ( pzz=0; pzz<z; pzz+=pz ) {
				i = (nn*z + pzz)*pagesize;
				fimg->read((char *) (data + i), pagesize);
				if ( pad ) fimg->seekg(pad, ios_base::cur);
			}
		}
	} else {					// More complicated paging (DSN6 and BRIX)
		off = offset;
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::read_data: Reading single values: pagesize=" << 
					pagesize << " off=" << off << " pad=" << pad << endl;
		npages = 1;
		npages *= (z - 1)/pz + 1;
		npages *= (y - 1)/py + 1;
		npages *= (xstore - 1)/xpage + 1;
		off = offset + img_select*npages*(pagesize + pad);
		fimg->seekg(off);
		page = new unsigned char[pagesize];
		for ( nn=0; nn<n; nn++ ) {
			for ( pzz=0; pzz<z; pzz+=pz ) {
				for ( pyy=0; pyy<y; pyy+=py ) {
					for ( pxx=0; pxx<xstore; pxx+=xpage ) {	//cout << pxx << tab << pyy << tab << pzz << endl;
						fimg->read((char *)page, pagesize);
						if ( sb > 1 )
							 for ( i=0; i<pagesize; i+=sb ) swapbytes(page+i, sb);
						for ( zz=0; zz<pz && zz<z-pzz; zz++ ) {
							iz = (nn*z + zz + pzz)*y;
							jz = zz*py;
							for ( yy=0; yy<py && yy<y-pyy; yy++ ) {
								iy = (iz + yy + pyy)*xstore;
								jy = (jz + yy)*xpage;
								for ( xx=0; xx<xpage && xx<xstore-pxx; xx++ ) {	 //cout << xx << tab << yy << tab << zz << endl;
									i = (iy + xx + pxx)*elementsize;
									j = (jy + xx)*elementsize;
									memcpy(data+i, page+j, elementsize);
								}
							}
						}
						if ( pad ) fimg->seekg(pad, ios_base::cur);
					}
				}
			}
		}
		if ( page ) delete [] page;
	}
	
	// Convert vax format and swap bytes if required
    if ( vax && sb < 2 && ( datatype == Float ) ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::read_data: Converting vax image data with sb = " << sb << endl;
		for ( i=0; i<datasize; i+=4 )
    	    vax2ieee(data+i, 1-sb);
    } else if ( datatypesize > 1 ) {
		if ( sb == 1 ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG Bimage::read_data: Swapping image data with typesize = " << datatypesize << endl;
			for ( i=0; i<datasize; i+=datatypesize ) 
				swapbytes(data+i, datatypesize);
		} else if ( sb > 1 ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG Bimage::read_data: Swapping image data with typesize = " << sb << endl;
			for ( i=0; i<datasize; i+=sb ) 
				swapbytes(data+i, sb);
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::read_data: Finished reading and converting data" << endl;

	return data;
}

int			Bimage::write(Bstring& fn)
{
	if ( verbose )
		cout << "Write image to " << fn << endl;
	
	ofstream		fimg(fn.c_str());
	if ( fimg.fail() )  return -1;

	fimg.write((char *)&datatype, 8);
	fimg.write((char *)&c, 5*8);
	fimg.write((char *)d.uc, alloc_size());
	
	fimg.close();
	
	return 0;
}

int			Bimage::unpack_transform(unsigned char* data, FourierType tf)
{
	return unpack_transform(-1, data, tf);
}

int			Bimage::unpack_transform(int img_select, unsigned char* data, FourierType tf)
{
	if ( tf <= Standard ) {
		cerr << "Error: The transform is not in the correct format to unpack!" << endl;
		exit(-1);
	}

	if ( verbose & VERB_PROCESS )
		cout << "Unpacking Fourier transform:    " << tf << endl << endl;
	
	long			i, h, xx, yy, zz, nn, hx, hy, hz;
	long			xstore = x, ox(0), oy(0), oz(0);
	long			elsize = c*data_type_size();
	Complex<double>	cv;

	long			nmin = (img_select >= 0)? img_select: 0;
	long			nmax = (img_select >= 0)? img_select: n - 1;
	
	switch ( tf ) {
		case Centered:
			ox = x/2;
			oy = y/2;
			oz = z/2;
			break;
		case Hermitian:
			xstore = x/2 + 1;
			break;
		case CentHerm:
			oy = y/2;
			oz = z/2;
			xstore = x/2 + 1;
			break;
		default: break;
	}
		
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::unpack_transform: data size: " << size() << endl;
		cout << "DEBUG Bimage::unpack_transform: element size: " << elsize << endl;
		cout << "DEBUG Bimage::unpack_transform: origin: " << ox << tab << oy << tab << oz << endl;
		cout << "DEBUG Bimage::unpack_transform: xstore: " << xstore << endl;
	}
	
	// Transfer the data
	for ( nn=nmin; nn<=nmax; nn++ ) {
		for ( hz=0; hz<z; hz++ ) {
			zz = (hz>=oz)? hz - oz: z + hz - oz;
			for ( hy=0; hy<y; hy++ ) {
				yy = (hy>=oy)? hy - oy: y + hy - oy;
				for ( hx=0; hx<xstore; hx++, data+=elsize ) {
					xx = (hx>=ox)? hx - ox: x + hx - ox;
					i = index(xx, yy, zz, nn);
					memcpy(data_pointer() + i*elsize, data, elsize);
				}
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::unpack_transform: transfer done" << endl;
	
	// Fill in the hermitian/friedel partners
	if ( xstore < x ) {
		for ( nn=nmin; nn<=nmax; nn++ ) {
			for ( zz=0; zz<z; zz++ ) {
				hz = (zz)? z - zz: 0;
				for ( yy=0; yy<y; yy++ ) {
					hy = (yy)? y - yy: 0;
					for ( xx=xstore; xx<x; xx++ ) {
						hx = x - xx;
						i = index(xx, yy, zz, nn);
						h = index(hx, hy, hz, nn);
						cv = complex(h);
						cv = cv.conj();
						set(i, cv);
					}
				}
			}
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::unpack_transform: friedel half done" << endl;
	}
	
	fourier_type(Standard);
	
	return 0;
}

int			Bimage::pack_transform(unsigned char* data, FourierType tf)
{
	return pack_transform(-1, data, tf);
}

int			Bimage::pack_transform(int img_select, unsigned char* data, FourierType tf)
{
	if ( tf <= Standard ) {
		cerr << "Error: The transform is not in the correct format to unpack!" << endl;
		exit(-1);
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Packing Fourier transform:    " << tf << endl << endl;
	
	long			i, xx, yy, zz, nn, hx, hy, hz;
	long			xstore = x, ox(0), oy(0), oz(0);
	long			elsize = c*data_type_size();
	
	long			nmin = (img_select >= 0)? img_select: 0;
	long			nmax = (img_select >= 0)? img_select: n - 1;
	
	switch ( tf ) {
		case Centered:
			ox = x/2;
			oy = y/2;
			oz = z/2;
			break;
		case Hermitian:
			xstore = x/2 + 1;
			break;
		case CentHerm:
			oy = y/2;
			oz = z/2;
			xstore = x/2 + 1;
			break;
		default: break;
	}
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::pack_transform: Element size = " << elsize << endl;
	
//	cout << "xstore = " << xstore << endl;
	
	// Transfer the data
	for ( nn=nmin; nn<=nmax; nn++ ) {
		for ( hz=0; hz<z; hz++ ) {
			zz = (hz>=oz)? hz - oz: z + hz - oz;
			for ( hy=0; hy<y; hy++ ) {
				yy = (hy>=oy)? hy - oy: y + hy - oy;
				for ( hx=0; hx<xstore; hx++, data+=elsize ) {
					xx = (hx>=ox)? hx - ox: x + hx - ox;
					i = index(xx, yy, zz, nn);
					memcpy(data, data_pointer() + i*elsize, elsize);
				}
			}
		}
	}
	
	return 0;
}


/**
@brief 	Prints out header information for an image.
@return int					error code (<0 means failure).
**/
int		Bimage::information()
{
	tm*			t = get_localtime();
	cout << "Time:                           " << t->tm_year+1900 << "-" << t->tm_mon+1 
		<< "-" << t->tm_mday << " " << t->tm_hour << ":" << t->tm_min << ":" << t->tm_sec << endl;
	if ( metadata.exists("user") ) cout << "User:                           " << metadata["user"] << endl;
	if ( metadata.exists("type") ) cout << "Type:                           " << metadata["type"] << endl;
	cout << "Data offset:                    " << offset << endl;
	cout << "Number of images:               " << n << endl;
	cout << "Dimensions:                     " << x << " " << y << " " << z << " voxels" << endl;
	cout << "Page dimensions:                " << px << " " << py << " " << pz << " voxels" << endl;
	cout << "Channels:                       " << c << endl;
	Bstring		datatype_string = data_type_string();
	cout << "Data type:                      " << datatype_string << 
			" (size = " << data_type_size() << ")" << endl;
	if ( fouriertype ) {
		cout << "Transform:                      ";
		switch ( fouriertype ) {
			case Standard: cout << "Standard" << endl; break;
			case Centered: cout << "Centered" << endl; break;
			case Hermitian: cout << "Hermitian" << endl; break;
			case CentHerm: cout << "Centered hermitian" << endl; break;
			default: break;
		}
	}
	Bstring		 compoundtype_string = compound_type_string();
	cout << "Compound type:                  " << compoundtype_string << endl;
	Vector3<double>	u(image->sampling());
	cout << "Voxel units/sampling:           " << u[0] << " " << u[1] << " " << u[2] << endl;
	cout << "Min, max, ave, std:             " << 
			min << " " << max << " " << avg << " " << std << endl;
	
    cout << "Text label:                     (length = " << label().length() << ")" << endl;
	cout << label() << endl;
	
	UnitCell		uc = unit_cell();
	Matrix3			frac_mat;
	if ( space_group() > 0 ) {
		frac_mat = uc.skew_matrix();
		cout.setf( ios::fixed, ios::floatfield );
		cout << "Unit cell dimensions:           " << setprecision(2) <<
				uc.a() << " " << uc.b() << " " << uc.c() << endl;
		cout << "Unit cell angles:               " << 
				uc.alpha()*180/M_PI << " " << uc.beta()*180/M_PI << " " << uc.gamma()*180/M_PI << endl;
		cout << "Fractionalization matrix:       " << setprecision(5) <<
				 endl << frac_mat;
//		cout << "Fractional translation:         " <<
//				frac_mat[9] << " " << frac_mat[10] << " " << frac_mat[11] << endl;
	    cout << "Space group:                    " << space_group() << endl;
	}
	
	if ( (symmetry()).length() )
		cout << "Symmetry label:                 " << symmetry() << endl;
	
	cout << endl;
	
	return 0;
}

/**
@brief 	Prints out header information for all sub-images.
@return int			error code (<0 means failure).
**/
int			Bimage::subimage_information()
{
	long		nn;
	
	cout << "\nImage\tMin\tMax\tAvg\tStd\tBkg" << endl;
	for ( nn=0; nn<n; nn++ )
		cout << nn+1 << tab << setprecision(4) << 
			image[nn].minimum() << tab <<
			image[nn].maximum() << tab <<
			image[nn].average() << tab <<
			image[nn].standard_deviation() << tab <<
			image[nn].background() << endl;

	cout << endl;

	cout << "Image\tSampling\t\tOrigin\t\t\tView" << endl;
	for ( nn=0; nn<n; nn++ )
		cout << nn+1 << tab << setprecision(3) <<
			image[nn].sampling() << tab << setprecision(2) <<
			image[nn].origin() << tab << setprecision(4) <<
			image[nn].view().vector3() << tab <<
			image[nn].view().angle()*180.0/M_PI << endl;
	
	cout << endl;

//	cout << "metadata type = " << metadata.type() << endl;

//	cout << metadata << endl;
	
	return 0;
}

/**
@brief 	Prints out moments for all sub-images.
@return int			error code (<0 means failure).
**/
int			Bimage::moments(long max_order)
{
	long		nn;
	
	for ( nn=0; nn<n; ++nn ) moments(max_order, nn);

	cout << endl;
	
	return 0;
}

/**
@brief 	Prints out moments for one sub-image.
@param	max_order	maximum order.
@param	nn			sub-image.
@return int			error code (<0 means failure).
**/
int			Bimage::moments(long max_order, long nn)
{
	if ( max_order < 2 ) max_order = 2;
	if ( max_order > 20 ) max_order = 20;
	
	cout << "Moments for image " << nn << ":" << endl;
	
	long			i, j, k, ds(image_size());
	double			v, vc, cm[max_order], m[max_order];
	
	for ( k=0; k<max_order; ++k ) cm[k] = m[k] = 0;
	
	for ( i=0, j=nn*ds; i<ds; ++i, ++j ) {
		v = (*this)[j];
		vc = 1;
		for ( k=0; k<max_order; ++k ) {
			vc *= v;
			m[k] += vc;
		}
	}

	for ( k=0; k<max_order; ++k ) m[k] /= ds;
	
	for ( i=0, j=nn*ds; i<ds; ++i, ++j ) {
		v = (*this)[j] - m[0];
		vc = 1;
		for ( k=0; k<max_order; ++k ) {
			vc *= v;
			cm[k] += vc;
		}
	}
	
	for ( k=0; k<max_order; ++k ) cm[k] /= ds;

	v = sqrt(cm[1]);	// StDev
	vc = 1;
	
	cout << "Order\tMoment\tCentral\tPearson" << endl;
	for ( k=0; k<max_order; ++k ) {
		vc *= v;
		cout << k+1 << tab << m[k] << tab << cm[k] << tab << cm[k]/vc << endl;
	}
	
	if ( max_order > 2) cout << "Skewness = " << cm[2]/(cm[1]*sqrt(cm[1])) << endl;
	if ( max_order > 3) cout << "Kurtosis = " << cm[3]/(cm[1]*cm[1]) - 3 << endl;
	cout << endl;
	
	return 0;
}

/**
@brief 	Prints out header information associated with a tag string.
@param	tag			tag string.

Recognized tags:
	size, ch(annel), x, y, z, im(age), sam(pling), ori(gen), view,
	min(imum), max(imum), av(erage), st(andard deviation), var(iance), stat(istics),
	text, json
	
**/
void		Bimage::get(Bstring tag)
{
	if ( tag == "size" ) cout << size() << endl;
	else if ( tag.contains("ch") ) cout << c << endl;
	else if ( tag == "x" ) cout << x << endl;
	else if ( tag == "y" ) cout << y << endl;
	else if ( tag == "z" ) cout << z << endl;
	else if ( tag.contains("im") ) cout << n << endl;
	else if ( tag.contains("sam") || tag.contains("pix") ) cout << sampling(0) << endl;
	else if ( tag.contains("ori") ) cout << image->origin() << endl;
	else if ( tag == "view" ) cout << image->view() << endl;
	else if ( tag.contains("min") ) cout << min << endl;
	else if ( tag.contains("max") ) cout << max << endl;
	else if ( tag.contains("av") ) cout << avg << endl;
	else if ( tag.contains("st") ) cout << std << endl;
	else if ( tag.contains("var") ) cout << std*std << endl;
	else if ( tag.contains("stat") ) cout << min << tab << max << tab << avg << tab << std << endl;
	else if ( tag == "text" ) cout << label() << endl;
	else if ( tag == "json" ) cout << metadata << endl;
}

/**
@brief Copies the header information and data of an image into a new image structure.
@return Bimage* 			copied image, NULL if copy failed.
**/
Bimage* 	Bimage::copy()
{
	long			j;
	Bimage*			img = copy_header();
	if ( !img ) {
		cerr << "Error: Bimage:copy failed!" << endl;
		return img;
	}
	
//	if ( verbose & VERB_DEBUG )
//		cout << "DEBUG Bimage::copy: alloc_size=" << alloc_size() << " data_size=" << data_size() << endl;
	
	img->data_alloc();
	for ( j=0; j<img->alloc_size(); j++ ) img->d.uc[j] = d.uc[j];
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::copy: done" << endl;
	
	return img;
}

/**
@brief Copies the header information and data of an image into a new image structure.
@param	nu_nimg		new number of images (if < 1, keep the old number)
@return Bimage* 	copied image, NULL if copy failed.
**/
Bimage* 	Bimage::copy(long nu_nimg)
{
	if ( nu_nimg < 1 ) return copy();
	
	long			j;
	Bimage*			img = copy_header(nu_nimg);
	if ( !img ) {
		cerr << "Error: Bimage:copy failed!" << endl;
		return img;
	}
	
//	if ( verbose & VERB_DEBUG )
//		cout << "DEBUG Bimage::copy: alloc_size=" << alloc_size() << " data_size=" << data_size() << endl;
	
	img->data_alloc();
	for ( j=0; j<alloc_size() && j<img->alloc_size(); j++ ) img->d.uc[j] = d.uc[j];
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::copy: done" << endl;
	
	return img;
}


/**
@brief	Copy an image structure into a new one.
@param	nu_nimg		new number of images (if < 1, keep the old number)
@return Bimage*		the new image structure, NULL if copy failed.

	All information from an old image structure is copied into a new one.
	This takes care of the internal arrays in the structure.
	The data pointer is left to point to the original data.

**/
Bimage* 	Bimage::copy_header(long nu_nimg)
{
	long			j;

	// Allocate memory for the image parameter structure
	Bimage*			img = new Bimage;
	if ( !img ) return img;
	
	img->offset = offset;
	
	img->size(x, y, z);
	img->channels(c);
	img->data_type(datatype);
	img->compound_type(compoundtype);
	img->page_size(px, py, pz);
	img->fourier_type(fouriertype);

//	cout << "compoundtype=" << img->compoundtype << " c=" << c << endl;
//	cout << "img compoundtype=" << img->compoundtype << " c=" << img->c << endl;

	// Copy all strings
	img->id = id;
	
	// Statistics
	img->minimum(min);
	img->maximum(max);
	img->average(avg);
	img->standard_deviation(std);
	img->show_minimum(smin);
	img->show_maximum(smax);
	img->show_image(sn);
	img->show_slice(sz);
	img->show_scale(ss);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::copy_header: New image header copied" << endl;
	
	// Copy metadata
	
	//	cout << "copy_header:" << endl << img->metadata << endl;
	
	// Create and copy new sub-images, taking care of metadata
	if ( nu_nimg < 1 ) {
		nu_nimg = n;
		img->metadata = metadata;
	} else {
		img->images(nu_nimg);
	}
	
	img->unit_cell(ucell);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::copy_header: New number of images = " << img->n << endl;
	
	for ( j=0; j<n && j<img->n; j++ ) img->image[j] = image[j];
	for ( ; j<img->n; j++ ) img->image[j] = image[0];

	img->data_size();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::copy_header: Units = " << img->sampling(0) << endl;
	
	return img;
}



/**
@brief Replaces the data with that from the given image.
@param 	*img		source image.
@return int			0.

	The input image must be the same size as the receiving image.

**/
int			Bimage::replace(Bimage* img)
{
	long		j;

	for ( j=0; j<datasize; j++ ) set(j, (*img)[j]);
	
	return 0;
}

/**
@brief Replaces one sub-image in an image structure.
@param 	nn			image number to replace.
@param 	*img		source image.
@param	nr			source sub-image number.
@return int			0.

	The input image must be the same size as the receiving image.

**/
int			Bimage::replace(long nn, Bimage* img, long nr)
{
	long		i, j, k, imgsize(x*y*z*c);

	image[nn] = img->image[nr];
	
	for ( i=nn*imgsize, j=nr*imgsize, k=0; k<imgsize; ++i, ++j, ++k ) set(i, (*img)[j]);
	
	return 0;
}

/**
@brief 	Replaces one sub-image in an image structure.
@param 	nn			image number to replace.
@param 	*img		source image.
@param	nr			source sub-image number.
@param	fill		fill vlue.
@return int			0.

	The input image may have a different size as the receiving image.

**/
int			Bimage::replace(long nn, Bimage* img, long nr, double fill)
{
	long		imgsize(x*y*z*c);
	long		i, j, xx, yy, zz, cc;

	image[nn] = img->image[nr];
	
	for ( j=imgsize*nn, i=imgsize*(nn+1); j<i; j++ ) set(j, fill);
	
	for ( zz=0; zz<z && zz<img->z; zz++ ) {
		for ( yy=0; yy<y && yy<img->y; yy++ ) {
			for ( xx=0; xx<x && xx<img->x; xx++ ) {
				i = img->index(0, xx, yy, zz, nr);
				j = index(0, xx, yy, zz, nn);
				for ( cc=0; cc<c; cc++, i++, j++ ) set(j, (*img)[i]);
			}
		}
	}
	
	return 0;
}

/**
@brief 	Finds the highest value in a kernel.
@param 	idx		index in multi-image.
@param 	ksize	kernel edge half size
@return long 	index of highest value.
**/
long		Bimage::kernel_min(long idx, long ksize)
{
	long			xx, yy, zz, nn(idx/(x*y*z));
	long			i, ihi(idx);
	double			min_val(1e30);
	
	//calculate the bounds on the kernel
	Vector3<long>	lo = kernel_low(idx, ksize);
	Vector3<long>	hi = kernel_high(idx, ksize);
	
	// Loop through the kernel and check neighbors
	for ( zz=lo[2]; zz<=hi[2]; zz++ ) {
		for ( yy=lo[1]; yy<=hi[1]; yy++ ) {
			for ( xx=lo[0]; xx<=hi[0]; xx++ ) {
				i = index(0,xx,yy,zz,nn);
				if ( min_val > (*this)[i] ) {
					min_val = (*this)[i];
					ihi = i;
				}
			}
		}
	}
	
	return ihi;
}

/**
@brief 	Finds the highest value in a kernel.
@param 	idx		index in multi-image.
@param 	ksize	kernel edge half size
@return long 		index of highest value.
**/
long		Bimage::kernel_max(long idx, long ksize)
{
	long			xx, yy, zz, nn(idx/(x*y*z));
	long			i, ihi(idx);
	double			max_val(-1e30);
	
	//calculate the bounds on the kernel
	Vector3<long>	lo = kernel_low(idx, ksize);
	Vector3<long>	hi = kernel_high(idx, ksize);
	
	// Loop through the kernel and check neighbors
	for ( zz=lo[2]; zz<=hi[2]; zz++ ) {
		for ( yy=lo[1]; yy<=hi[1]; yy++ ) {
			for ( xx=lo[0]; xx<=hi[0]; xx++ ) {
				i = index(0,xx,yy,zz,nn);
				if ( max_val < (*this)[i] ) {
					max_val = (*this)[i];
					ihi = i;
				}
			}
		}
	}
	
	return ihi;
}

/**
@brief 	Calculates the average value in a kernel.
@param 	idx		index in multi-image.
@param 	ksize	kernel edge half size.
@param 	tmin	miminum to exclude.
@param 	tmax	maximum to exclude
@return double 	average value.
**/
double		Bimage::kernel_average(long idx, long ksize, double tmin, double tmax)
{
	
	long			xx, yy, zz, nn(idx/(x*y*z));
	long			i, nv(0);
	double			v1, avg_val(0);
	
	//calculate the bounds on the kernel
	Vector3<long>	lo = kernel_low(idx, ksize);
	Vector3<long>	hi = kernel_high(idx, ksize);
	
	// Loop through the kernel
	for ( zz=lo[2]; zz<=hi[2]; zz++ ) {
		for ( yy=lo[1]; yy<=hi[1]; yy++ ) {
			for ( xx=lo[0]; xx<=hi[0]; xx++ ) {
				i = index(0,xx,yy,zz,nn);
				v1 = (*this)[i];
				if ( v1 >= tmin && v1 <= tmax ) {
					avg_val += v1;
					nv++;
				}
			}
		}
	}
	
	if ( nv ) avg_val /= nv;
	else avg_val = (tmin + tmax)/2;
	
	return avg_val;
}

/**
@brief 	Calculates the sum in a kernel.
@param 	idx		index in multi-image.
@param 	ksize	kernel edge half size.
@return double 	sum.
**/
double		Bimage::kernel_sum(long idx, long ksize)
{
	
	long			i, xx, yy, zz, nn(idx/(x*y*z));
	double			sum(0);
	
	//calculate the bounds on the kernel
	Vector3<long>	lo = kernel_low(idx, ksize);
	Vector3<long>	hi = kernel_high(idx, ksize);
	
	// Loop through the kernel
	for ( zz=lo[2]; zz<=hi[2]; zz++ ) {
		for ( yy=lo[1]; yy<=hi[1]; yy++ ) {
			for ( xx=lo[0]; xx<=hi[0]; xx++ ) {
				i = index(0,xx,yy,zz,nn);
				sum += (*this)[i];
			}
		}
	}
	
	return sum;
}

/**
@brief 	Finds the highest value in a kernel.
@param 	idx		index in multi-image.
@param 	ksize	kernel edge half size.
@return double 	average value.
**/
double		Bimage::kernel_neighbor_average(long idx, long ksize)
{
	
	long			xx, yy, zz, nn(idx/(x*y*z));
	long			i, nv(0);
	double			v1, avg_val(0);
	
	//calculate the bounds on the kernel
	Vector3<long>	lo = kernel_low(idx, ksize);
	Vector3<long>	hi = kernel_high(idx, ksize);
	
	// Loop through the kernel
	for ( zz=lo[2]; zz<=hi[2]; zz++ ) {
		for ( yy=lo[1]; yy<=hi[1]; yy++ ) {
			for ( xx=lo[0]; xx<=hi[0]; xx++ ) {
				i = index(xx,yy,zz,nn);
				if ( compoundtype == TComplex ) {
					v1 = complex(i).power();
				} else {
					v1 = (*this)[i];
				}
				if ( i != idx ) {
					avg_val += v1;
					nv++;
				}
			}
		}
	}
	
	if ( nv ) avg_val /= nv;
	
	return avg_val;
}

/**
@brief 	Finds the highest value in a kernel excluding the central voxel.
@param 	idx		index in multi-image.
@param 	ksize	kernel edge half size
@return long 		index of highest value.
**/
long		Bimage::kernel_max_neigbor(long idx, long ksize)
{
	long			xx, yy, zz, nn(idx/(x*y*z));
	long			i, ihi(idx);
	double			max_val(-1e30);
	
	//calculate the bounds on the kernel
	Vector3<long>	lo = kernel_low(idx, ksize);
	Vector3<long>	hi = kernel_high(idx, ksize);
	
	// Loop through the kernel and check neighbors
	for ( zz=lo[2]; zz<=hi[2]; zz++ ) {
		for ( yy=lo[1]; yy<=hi[1]; yy++ ) {
			for ( xx=lo[0]; xx<=hi[0]; xx++ ) {
				i = index(0,xx,yy,zz,nn);
				if ( max_val < (*this)[i] ) {
					if ( i != idx ) {
						max_val = (*this)[i];
						ihi = i;
					}
				}
			}
		}
	}
	
	return ihi;
}

/**
@brief 	Finds the highest value in a kernel with wrapping.
@param 	idx		index in multi-image.
@param 	ksize	kernel edge half size
@return long 	index of highest value.
**/
long		Bimage::kernel_max_wrap(long idx, long ksize)
{
	long			xx, yy, zz, ix, iy, iz, nn(idx/(x*y*z));
	long			i, ihi(idx);
	double			max_val(-1e30);
	
	//calculate the bounds on the kernel
	Vector3<long>	lo(coordinates(idx)-ksize);
	Vector3<long>	hi(coordinates(idx)+ksize);
	
	// Loop through the kernel and check neighbors
	for ( zz=lo[2]; zz<=hi[2]; zz++ ) {
		iz = zz;
		while ( iz < 0 ) iz += z;
		while ( iz >= z ) iz -= z;
		for ( yy=lo[1]; yy<=hi[1]; yy++ ) {
			iy = yy;
			while ( iy < 0 ) iy += y;
			while ( iy >= y ) iy -= y;
			for ( xx=lo[0]; xx<=hi[0]; xx++ ) {
				ix = xx;
				while ( ix < 0 ) ix += x;
				while ( ix >= x ) ix -= x;
				i = index(0,ix,iy,iz,nn);
				if ( max_val < (*this)[i] ) {
					max_val = (*this)[i];
					ihi = i;
				}
			}
		}
	}
	
	return ihi;
}

/**
@brief 	Orders the values in a kernel.
@param 	idx			index in multi-image.
@param 	ksize		kernel edge half size
@return multimap<double,long> 	index of highest value.
**/
multimap<double,long>	Bimage::kernel_order(long idx, long ksize)
{
	long			i, xx, yy, zz, nn(idx/(x*y*z));
	multimap<double,long>	korder;
	
	//calculate the bounds on the kernel
	Vector3<long>	lo = kernel_low(idx, ksize);
	Vector3<long>	hi = kernel_high(idx, ksize);
	
	// Loop through the kernel and check neighbors
	for ( zz=lo[2]; zz<=hi[2]; zz++ ) {
		for ( yy=lo[1]; yy<=hi[1]; yy++ ) {
			for ( xx=lo[0]; xx<=hi[0]; xx++ ) {
				i = index(0,xx,yy,zz,nn);
				korder.emplace((*this)[i],i);
			}
		}
	}
	
	return korder;
}

/**
@brief 	Orders the neigbor values in a kernel.
@param 	idx			index in multi-image.
@param 	ksize		kernel edge half size
@return multimap<double,long> 	index of highest value.
**/
multimap<double,long>	Bimage::kernel_order_neighbors(long idx, long ksize)
{
	long			i, xx, yy, zz, nn(idx/(x*y*z));
	multimap<double,long>	korder;
	
	//calculate the bounds on the kernel
	Vector3<long>	lo = kernel_low(idx, ksize);
	Vector3<long>	hi = kernel_high(idx, ksize);
	
	// Loop through the kernel and check neighbors
	for ( zz=lo[2]; zz<=hi[2]; zz++ ) {
		for ( yy=lo[1]; yy<=hi[1]; yy++ ) {
			for ( xx=lo[0]; xx<=hi[0]; xx++ ) {
				i = index(0,xx,yy,zz,nn);
				if ( i != idx )
					korder.emplace((*this)[i],i);
			}
		}
	}
	
	return korder;
}

/**
@author Daniel Nemecek and Bernard Heymann
@brief	Calculates the density in a sphere around a coordinate in an image.
@param 	nn				sub-image number.
@param 	coord			position of density in map (voxels).
@param 	radius			spherical radius.
@param	sigma			standard deviation result.
@return double			average density.

	The density origin is positioned on the component.

**/
double		Bimage::density(long nn, Vector3<double> coord, double radius, double& sigma)
{
	long			i, xx, yy, zz;
	double			x2, y2, z2, rad2(radius*radius);
	double			density(0), dist2, v, w, sum(0), sum2(0), sumw(0), fac(-9/rad2);
	
	Vector3<double>	lo = coord - radius;
	Vector3<double>	hi = coord + radius;
	lo = lo.max(0);
	if ( hi[0] >= x ) hi[0] = x - 1;
	if ( hi[1] >= y ) hi[1] = y - 1;
	if ( hi[2] >= z ) hi[2] = z - 1;

	for ( zz=(long)lo[2]; zz<=hi[2]; zz++ ) {
		z2 = zz - coord[2];
		z2 *= z2;
		for ( yy=(long)lo[1]; yy<=hi[1]; yy++ ) {
			y2 = yy - coord[1];
			y2 *= y2;
			for ( xx=(long)lo[0]; xx<=hi[0]; xx++ ) {
				x2 = xx - coord[0];
				x2 *= x2;
				dist2 = x2 + y2 + z2;
				if ( dist2 <= rad2 ) {
					i = index(0, xx, yy, zz, nn);
//					if ( (*this)[i] > average() ) {
						w = exp(fac*dist2);
						v = (*this)[i];
						sum += w*v;
						sum2 += w*v*v;
						sumw += w;
//					}
				}
			}
		}
	}
	
	if ( sumw ) {
		density = sum/sumw;
		sigma = sum2/sumw - density*density;
		if ( sigma > 0 ) sigma = sqrt(sigma);
		else sigma = 0;
	}
	
	return density;
}

/**
@brief 	Calculates the relative density in a region defined by a mask.
@param 	*pmask		a 4 level mask.
@return double		relative density.

	The mask is assumed to have 4 levels:
	0		region of no interest
	1		region to estimate relative density
	2		reference region
	3		background region

**/
double		Bimage::relative_density(Bimage* pmask)
{
	int				t;
	long			i, np(0), nr(0), nb(0);
	double			v, dp(0), dr(0), db(0);
	
	for ( i=0; i<datasize; i++ ) {
		t = (int) ((*pmask)[i] + 0.5);
		v = (*this)[i];
		if ( t == 1 ) {
			dp += v;	// Region of interest
			np++;
		} else if ( t == 2 ) {
			dr += v;	// Reference region
			nr++;
		} else if ( t == 3 ) {
			db += v;	// Background region
			nb++;
		}
	}
	
	if ( np ) dp /= np;
	if ( nr ) dr /= nr;
	if ( nb ) db /= nb;
	
	if ( dr - db ) dp = (dp - db)/(dr - db);
	else if ( dr ) dp /= dr;
	
	return dp;
}


/**
@brief Inverts the data in the image.

	The Bit data type is negated.
	An unsigned data type is subtracted from its type maximum.
	A signed data type is negated.

**/
void		Bimage::invert()
{
	long			j, ds;
	
	switch ( datatype ) {
		case Bit: ds = (px/8)*y*z*n; for ( j=0; j<ds; j++ ) d.uc[j] = ~d.uc[j]; break;
		case UCharacter: case UShort: case UInteger: case ULong:
			for ( j=0; j<datasize; j++ ) set(j, dtmax - (*this)[j]);
			break;
		default: 
			for ( j=0; j<datasize; j++ ) set(j, -(*this)[j]);
			break;
	}
	
	statistics();
}

/**
@brief Switches axes of an image.
@param order 	string encoding the reslicing order.

	Reslicing refers to switching the axes of a multidimensional image.
	The equivalent is rotation around angles that are multiples of pi/2.
	The reslicing is encoded as a string of x, y and z, with the order
	and sign indicating the type of reslicing to be done.
	The new data replaces the old data.

**/
void 		Bimage::reslice(Bstring order)
{
	if ( !d.uc ) return;

	long			i, j;
	Vector3<long>	nusize, coor, nucoor;
	Vector3<long>	sign(1,1,1), ind(0,1,2);
	Vector3<double>	u(image->sampling()), un;

	order.lower();
	
	if ( !order.contains("x") || !order.contains("y") || !order.contains("z") ) {
		cerr << "Error: Reslice order specification not correct: " << order << endl << endl;
		return;
	}

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::reslice: size" << tab << size() << endl;
		cout << "DEBUG Bimage::reslice:";
	}

	j = 0;
	for ( i=0; i<order.length(); i++ ) {
		if ( order[i] == 'x' ) {
			ind[j] = 0; nusize[j] = x; un[j] = u[0];
			j++;
		}
		if ( order[i] == 'y' ) {
			ind[j] = 1; nusize[j] = y; un[j] = u[1];
			j++;
		}
		if ( order[i] == 'z' ) {
			ind[j] = 2; nusize[j] = z; un[j] = u[2];
			j++;
		}
		if ( j < 3 ) {
			if ( order[i] == '-' ) sign[j] = -1;
			else sign[j] = 1;
		}
		if ( verbose & VERB_DEBUG )
			cout << tab << order[i] << " " << sign[j];
	}

	if ( verbose & VERB_DEBUG )
		cout << endl << "DEBUG Bimage::reslice: nusize" << tab << nusize << endl;
	
	long			elementsize(c*data_type_size());
	long			ds(x*y*z*n);
	unsigned char*	nudata = new unsigned char[alloc_size()];
	for ( i=0; i<alloc_size(); i++ ) nudata[i] = 0;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::reslice: elementsize=" << elementsize << " datasize=" << ds << endl;
	
	long	   		nn, xx, yy, zz;
	Vector3<double>	ori, nuori;
	
	for ( i=nn=0; nn<n; nn++ ) {
		ori = image[nn].origin();
		for ( j=0; j<3; j++ ) {
			nuori[j] = ori[ind[j]];
			if ( sign[j] < 0 ) nuori[j] = nusize[j] - 1 - nuori[j];
		}
		image[nn].origin(nuori);
		for ( zz=0; zz<z; zz++ ) {
			coor[2] = zz;
			for ( yy=0; yy<y; yy++ ) {
				coor[1] = yy;
				for ( xx=0; xx<x; xx++, i++ ) {
					coor[0] = xx;
					for ( j=0; j<3; j++ ) {
						nucoor[j] = coor[ind[j]];
						if ( sign[j] < 0 ) nucoor[j] = nusize[j] - 1 - nucoor[j];
					}
					j = ((nn*nusize[2] + nucoor[2])*nusize[1] + nucoor[1])*nusize[0] + nucoor[0];
					if ( datatype == Bit ) {
						if ( (*this)[i] ) nudata[j/8] |= 0x80 >> ((i/x)*px + i%x)%8;
					} else {
						memcpy(nudata+j*elementsize, d.uc+i*elementsize, elementsize);
//						for ( cc=0, ic=i*c, jc=j*c; cc<c; cc++, ic++, jc++ ) nudata[jc] = (*p)[ic];
					}
				}
			}
		}
	}
	
	size(nusize);
	sampling(un);
	data_assign(nudata);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::reslice: alloc_size=" << alloc_size() << endl;
	
	return;
}

/**
@brief 	Converts to absolute values.
**/
void		Bimage::absolute()
{
	if ( verbose & VERB_FULL )
		cout << "Converting to absolute values" << endl << endl;
	
	for ( long j=0; j<datasize; j++ )
		set(j, fabs((*this)[j]));
	
	statistics();
}

/**
@brief	Adds a constant value to an image.
@param 	v		constant to be added.
**/
void		Bimage::add(double v)
{
	if ( verbose & VERB_FULL )
		cout << setprecision(6) << "Adding " << v << endl << endl;
	
	double			v1;
	
	for ( long j=0; j<datasize; j++ ) {
		v1 = (*this)[j] + v;
		set(j, v1);
	}
	
	for ( long nn=0; nn<n; nn++ )
		image[nn].background(background(nn) + v);

	statistics();
}

/**
@brief	Adds a constant value to a phase image and wraps as necssary.
@param 	v		constant to be added.
**/
void		Bimage::phase_add(double v)
{
	if ( verbose & VERB_FULL )
		cout << setprecision(6) << "Adding " << v << " to phase image" << endl << endl;
	
	double			v1;
	
	for ( long j=0; j<datasize; j++ ) {
		v1 = angle_set_negPI_to_PI((*this)[j] + v);
		set(j, v1);
	}
	
	for ( long nn=0; nn<n; nn++ )
		image[nn].background(angle_set_negPI_to_PI(background(nn) + v));

	statistics();
}

/**
@brief	Multiplies an image with a constant value.
@param 	v		constant multiplier.
**/
void		Bimage::multiply(double v)
{
	if ( verbose & VERB_FULL )
		cout << setprecision(6) << "Multiplying with " << v << endl << endl;
	
	double			v1;
	
	for ( long j=0; j<datasize; j++ ) {
		v1 = (*this)[j] * v;
		set(j, v1);
	}

	for ( long nn=0; nn<n; nn++ )
		image[nn].background(background(nn) * v);

	statistics();
}

/**
@brief Multiplies a sub-image with a constant value.
@param	nn		sub-image.
@param 	v		constant multiplier.
**/
void		Bimage::multiply(long nn, double v)
{
	if ( verbose & VERB_FULL )
		cout << setprecision(6) << "Multiplying with " << v << endl << endl;
	
	long			i, j, imgsize(c*x*y*z);
	double			v1;
	
	for ( i=0, j=nn*imgsize; i<imgsize; i++, j++ ) {
		v1 = (*this)[j] * v;
		set(j, v1);
	}
	
	image[nn].background(background(nn) * v);
	
	statistics();
}

/**
@brief Calculates the power of an image.
@param 	v		power value.
**/
void		Bimage::power(double v)
{
	if ( verbose & VERB_PROCESS )
		cout << "Raising to power " << v << endl << endl;
	
	double			v1;
	
	for ( long j=0; j<datasize; j++ ) {
		v1 = pow((*this)[j], v);
		set(j, v1);
	}
	
	statistics();
}

/**
@brief Calculates the sine of a phase image.

	The values must be in radians.
	
**/
void		Bimage::sine()
{
	if ( verbose & VERB_PROCESS )
		cout << "Calculating a sine image" << endl << endl;
	
	for ( long j=0; j<datasize; ++j )
		set(j, sinl((*this)[j]));
	
	statistics();
}

/**
@brief Calculates the arcsine of an image.

	The values are truncated to be within [-1,1].
	
**/
void		Bimage::arcsine()
{
	if ( verbose & VERB_PROCESS )
		cout << "Calculating an arcsine image" << endl << endl;
	
	double		v;
	
	for ( long j=0; j<datasize; ++j ) {
		v = (*this)[j];
		if ( v < -1 ) v = -1;
		else if ( v > 1 ) v = 1;
		set(j, asinl(v));
	}
	
	statistics();
}

/**
@brief Calculates the cosine of a phase image.

	The values must be in radians.
	
**/
void		Bimage::cosine()
{
	if ( verbose & VERB_PROCESS )
		cout << "Calculating a cosine image" << endl << endl;
	
	for ( long j=0; j<datasize; ++j )
		set(j, cosl((*this)[j]));
	
	statistics();
}

/**
@brief Calculates the arccosine of an image.

	The values are truncated to be within [-1,1].
	
**/
void		Bimage::arccosine()
{
	if ( verbose & VERB_PROCESS )
		cout << "Calculating an arcsine image" << endl << endl;
	
	double		v;
	
	for ( long j=0; j<datasize; ++j ) {
		v = (*this)[j];
		if ( v < -1 ) v = -1;
		else if ( v > 1 ) v = 1;
		set(j, acosl(v));
	}
	
	statistics();
}

/**
@brief Calculates the tangent of a phase image.

	The values must be in radians.
	
**/
void		Bimage::tangent()
{
	if ( verbose & VERB_PROCESS )
		cout << "Calculating a tangent image" << endl << endl;
	
	for ( long j=0; j<datasize; ++j )
		set(j, tanl((*this)[j]));
	
	statistics();
}

/**
@brief Calculates the arctangent of an image.

	The values are truncated to be within [-1,1].
	
**/
void		Bimage::arctangent()
{
	if ( verbose & VERB_PROCESS )
		cout << "Calculating an arctangent image" << endl << endl;
	
	double		v;
	
	for ( long j=0; j<datasize; ++j ) {
		v = (*this)[j];
		if ( v < -1 ) v = -1;
		else if ( v > 1 ) v = 1;
		set(j, atanl(v));
	}
	
	statistics();
}

/**
@brief 	Adds all sub-images.
**/
void		Bimage::sum_images()
{
	if ( n < 2 ) return;
	
	long			i, j, nn, img_size(c*x*y*z);
	
	float*			nudata = new float[img_size];
	
	for ( j=0; j<img_size; j++ ) nudata[j] = 0;
	
	for ( nn=i=0; nn<n; nn++ )
		for ( j=0; j<img_size; j++, i++ ) nudata[j] += (*this)[i];
	
	n = 1;
	Bsub_image* 	nusub = new Bsub_image[1];
	nusub[0] = image[0];
	delete[] image;
	image = nusub;
	
	data_type(Float);
	data_assign((unsigned char *) nudata);

	statistics();
}

/**
@brief 	Averages all sub-images and optionally calculates standard deviation image.
@param	sd			flag to calculate linked standard deviation image.
@return Bimage*		average image.
**/
Bimage*		Bimage::average_images(bool sd)
{
	if ( n < 2 ) return NULL;
	
	Bimage*			pavg = new Bimage(Float, TSimple, size(), 1);
	pavg->clear();
	
	Bimage*			pstd = NULL;
	if ( sd ) pstd = pavg->next = pavg->copy();
	
	long			i, j, nn, img_size(c*x*y*z);
	double 			v, w(1.0/n);
	
	for ( nn=i=0; nn<n; nn++ )
		for ( j=0; j<img_size; j++, i++ ) {
			v = (*this)[i];
			pavg->add(j, v);
			if ( sd ) pstd->add(j, v*v);
		}
	
	pavg->multiply(w);
	
	if ( sd ) {
		pstd->multiply(w);
		for ( j=0; j<img_size; j++ ) {
			v = (*pstd)[j] - (*pavg)[j]*(*pavg)[j];
			if ( v > 0 ) v = sqrt(v);
			else v = 0;
			pstd->set(j, v);
		}
	}
	
	pavg->statistics();
	if ( sd ) pstd->statistics();
	
	return pavg;
}

/**
@brief 	Calculates a moving sum of the sub-images.
@param	window		number of successive images to sum.
@param	step		intervals between windows.
@param	flag		if set, calculate average.
@return Bimage*		moving sum/average image.

	Each successive subset of sub-images of size set by the window parameter are summed.
	Where the window extends beyond the limits for the number of images,
	the summation is just over the existing images.
	The number of images depends on the step size.
**/
Bimage*		Bimage::moving_sum(long window, long step, int flag)
{
	if ( n < 2 ) return NULL;
	if ( window > n ) window = n;
	if ( step < 1 ) step = 1;
	
	long			nimg = (n+step-1)/step;

	if ( verbose & VERB_PROCESS )
		cout << "Moving sum window " << window << " and step " << step << endl;

	Bimage*			pma = copy_header(nimg);
	pma->data_alloc_and_clear();
	
	long			i, j, k, nn, hw((window-1)/2), w, ws, we, img_size(c*x*y*z);
	double 			v;
	
	for ( nn = 0; nn < nimg; ++nn ) {
		ws = nn*step - hw;
		we = ws + window;
		if ( ws < 0 ) ws = 0;
		if ( we > n ) we = n;
		if ( verbose & VERB_PROCESS )
			cout << "Window " << nn+1 << ": " << ws << " - " << we-1 << endl;
		for ( w = ws; w < we; ++w ) {
			i = w*img_size;
			j = nn*img_size;
			for ( k = 0; k < img_size; ++i, ++j, ++k ) {
				v = (*this)[i];
				pma->add(j, v);
			}
		}
		if ( flag ) {
			w = we - ws;
			for ( j = nn*img_size, k = 0; k < img_size; ++j, ++k )
				pma->set(j, (*pma)[j]/w);
		}
	}
	
	pma->statistics();
	
	return pma;
}

/**
@brief 	Progressive sum of the sub-images.
@return int		0.

	Each sub-image is summed with all previous sub-images.
**/
void		Bimage::progressive_sum()
{
	long		i, j, k, nn, is(x*y*z*c);
	
//	cout << "size=" << is << endl;
//	cout << "c=" << c << " ct=" << compoundtype << endl;
	
	for ( nn=1; nn<n; ++nn )
		for ( i=nn*is, j=i-is, k=0; k<is; ++i, ++j, ++k )
			add(i, (*this)[j]);
}


/**
@brief 	Adds another image to an image.
@param 	*p		image to be added.

	Requirement: The images must have the same size.

**/
void		Bimage::add(Bimage* p)
{
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::add", __FILE__, __LINE__);
		return;
	}
	
	long			j;
	double			v1;
	
	for ( j=0; j<datasize; j++ ) {
		v1 = (*this)[j] + (*p)[j];
		set(j, v1);
	}

	statistics();
}

/**
@brief 	Subtracts another image from an image.
@param 	*p		image to be added.

	Requirement: The images must have the same size.

**/
void		Bimage::subtract(Bimage* p)
{
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::subtract", __FILE__, __LINE__);
		return;
	}
	
	double			v1;
	
	for ( long j=0; j<datasize; j++ ) {
		v1 = (*this)[j] - (*p)[j];
		set(j, v1);
	}
	
	statistics();
}

/**
@brief 	Adds another image to a sub-image.
@param	nn		sub-image.
@param 	*p		image to be added.

	Requirement: The images must have the same size.
	If the FOM block is defined, it is used to accumulate a sum of:
		the FOM block of the added image if it exists,
		or otherwise the square of the data added.

**/
void		Bimage::add(long nn, Bimage* p)
{
	long			i, j, imgsize(c*image_size());
	double			v1;
	
	if ( verbose & VERB_FULL )
		cout << "Adding to image " << nn << endl << endl;
	
	for ( i=0, j=nn*imgsize; i<imgsize; i++, j++ ) {
		v1 = (*this)[j] + (*p)[i];
		set(j, v1);
	}

	statistics();
}

/**
@brief 	Adds another image to an image.
@param 	*p		image to be added.
@param 	scale 	density scale to other image
@param 	shift 	density shift to other image.

	Requirement: The images must have the same size.

**/
void		Bimage::add(Bimage* p, double scale, double shift)
{
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::add", __FILE__, __LINE__);
		return;
	}
	
	if ( verbose & VERB_FULL )
		cout << "Adding images" << endl << endl;
	
	double			v1;
	
	for ( long j=0; j<datasize; j++ ) {
		v1 = (*this)[j] + (*p)[j] * scale + shift;
		set(j, v1);
	}
	
	statistics();
}

/**
@brief 	Adds another image to an image.
@param	nn		sub-image.
@param 	*p		image to be added.
@param 	scale 	density scale to other image
@param 	shift 	density shift to other image.

	Requirement: The images must have the same size.

**/
void		Bimage::add(long nn, Bimage* p, double scale, double shift)
{
	long			i, j, imgsize(c*image_size());
	double			v1;
	
	if ( verbose & VERB_FULL )
		cout << "Adding image " << nn << endl << endl;
	
	for ( i=0, j=nn*imgsize; i<imgsize; i++, j++ ) {
		v1 = (*this)[j] + (*p)[i] * scale + shift;
		set(j, v1);
	}
	
	statistics();
}

/**
@brief 	Multiplies an image with another image.
@param 	*p		image multiplier.
@param 	scale 	density scale to other image
@param 	shift 	density shift to other image.

	Requirement: The images must have the same size.

**/
void		Bimage::multiply(Bimage* p, double scale, double shift)
{
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::multiply", __FILE__, __LINE__);
		return;
	}
	
	if ( verbose & VERB_FULL )
		cout << "Multiplying images" << endl << endl;
	
	double			v1;
	
	for ( long j=0; j<datasize; j++ ) {
		v1 = (*this)[j] * ((*p)[j] * scale + shift);
		set(j, v1);
	}
	
	statistics();
}

/**
@brief 	Multiplies a sub-image with another image.
@param 	nn		sub-image to multiply.
@param 	*p		image multiplier.

	Requirement: The images must have the same size.

**/
void		Bimage::multiply(long nn, Bimage* p)
{
	if ( nn >= n ) return;
	
	if ( verbose & VERB_FULL )
		cout << "Multiplying image " << nn << endl << endl;
	
	long			imgsize(x*y*z*c), i(nn*imgsize), j(0), k(0);
	
	if ( p->n == n ) j = i;
	
	for ( k=0; k<imgsize; ++i, ++j, ++k )
		set(i, (*this)[i] * (*p)[j]);
}

/**
@brief 	Multiplies all sub-images with the first sub-image of the other image.
@param 	*p		image multiplier.

	The other image is multiplied with the first.
	Both images are converted to floating point.

**/
void 		Bimage::multiply(Bimage* p)
{
	if ( verbose &  VERB_PROCESS )
		cout << "Multiplying " << n << " sub-images with " << p->n << " sub-images" << endl << endl;
		
    if ( datatype >= Float || p->datatype >= Float )
    	change_type(Float);
    
#ifdef HAVE_GCD
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
		multiply(nn, p);
	});
#else
#pragma omp parallel for
	for ( long nn=0; nn<n; ++nn )
		multiply(nn, p);
#endif
	
	statistics();
}

/**
@brief 	Divides the first image by the second.
@param 	*p		other image.
@param 	scale 	density scale to apply to other image
@param 	shift 	density shift to apply to other image.

	The image is divided by the other with exceptions:
		image /= other*scale + shift
	Both images are converted to floating point.

**/
void		Bimage::divide(Bimage* p, double scale, double shift)
{
	double			minval = p->min*scale + shift;

	if ( size() != p->size() || n != p->n ) {
		check_if_same_size(p);
		error_show("Bimage::divide", __FILE__, __LINE__);
		return;
	}
	
	if ( minval < 0 ) {
		cerr << "Error: The scaled and shifted second image must be positive! (min = "
			<< minval << ")" << endl << endl;
		return;
	}
	
    long			j;
	double			div;
	Complex<float>	cv;
//	if ( minval <= 0 ) minval = p->standard_deviation()/1000;
	if ( minval < 1e-10 ) minval = scale;
	if ( minval < 1e-10 ) minval = shift;
	if ( minval < 1e-10 ) minval = 1e-10;
    
	if ( compound_type() != TComplex ) {
		for ( j=0; j<datasize; j++ ) {
			div =  (*p)[j] * scale + shift;
			if ( div < minval ) div = minval;
			set(j, (*this)[j]/div);
		}
	} else if ( p->compound_type() != TComplex ) {
		for ( j=0; j<datasize; j++ ) {
			div =  (*p)[j] * scale + shift;
			if ( div < minval ) div = minval;
			set(j, complex(j)/div);
		}
	} else {
		for ( j=0; j<datasize; j++ ) {
			cv = p->complex(j) * scale + shift;
			div = cv.power();
			if ( div < minval ) div = minval;
			cv = cv.conj();
			set(j, (p->complex(j)*cv)/div);
		}
	}
	
	statistics();
}

/**
@brief 	Divides the first image by the first sub-image of the second.
@param 	*p		other image.
@param 	scale 	density scale to apply to other image
@param 	shift 	density shift to apply to other image.

	The image is divided by the other with exceptions:
		image /= other*scale + shift
	Both images are converted to floating point.

**/
void		Bimage::divide_one(Bimage* p, double scale, double shift)
{
	double			minval = p->image->minimum()*scale + shift;

	if ( minval < 0 ) {
		cerr << "Error: The scaled and shifted second image must be positive! (min = "
			<< minval << ")" << endl << endl;
		return;
	}
	
    long			i, j, nn, imgsize(x*y*z*c);
	double			div;
	Complex<float>	cv;
//	if ( minval <= 0 ) minval = p->standard_deviation()/1000;
	if ( minval < 1e-10 ) minval = scale;
	if ( minval < 1e-10 ) minval = shift;
	if ( minval < 1e-10 ) minval = 1e-10;
    
	if ( compound_type() != TComplex ) {
		for ( i=nn=0; nn<n; nn++ ) {
			for ( j=0; j<imgsize; j++, i++ ) {
				div =  (*p)[j] * scale + shift;
				if ( div < minval ) div = minval;
				set(i, (*this)[i]/div);
			}
		}
	} else {
		imgsize /= 2;
		for ( i=nn=0; nn<n; nn++ ) {
			for ( j=0; j<imgsize; j++, i++ ) {
				cv = p->complex(j) * scale + shift;
				div = cv.power();
				if ( div < minval ) div = minval;
				cv = cv.conj();
				set(i, (p->complex(i)*cv)/div);
			}
		}
	}
	
	statistics();
}

/**
@brief 	Calculates the inverse of the image.
@param 	minval 	the minimum absolute value considered not zero.

	The inverse of every pixel is calculated.
	If the minumum value is given as zero, zero pixels are retained.

**/
void		Bimage::inverse(double minval)
{
    long			j;
	double			v;
	double			invminval(0);
	if ( minval ) invminval = 1.0/minval;
    
	for ( j=0; j<datasize; j++ ) {
		v = (*this)[j];
		if ( fabs(v) > minval )
			set(j, 1/v);
		else if ( v > 0 )
			set (j, invminval);
		else
			set(j, -invminval);
	}
	
	statistics();
}

/**
@brief 	Selects the largest of each pixel from two images.
@param 	*p		other image.
**/
void		Bimage::largest(Bimage* p)
{
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::largest", __FILE__, __LINE__);
		return;
	}
	
    if ( verbose & VERB_PROCESS )
		cout << "Setting to the largest" << endl;

	double		v1;
	for ( long j=0; j<datasize; j++ ) {
		v1 = (*p)[j];
    	if ( (*this)[j] < v1 ) set(j, v1);
    }
	
	statistics();
}

/**
@brief 	Selects the smallest of each pixel from two images.
@param 	*p		other image.

**/
void		Bimage::smallest(Bimage* p)
{
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::smallest", __FILE__, __LINE__);
		return;
	}

    if ( verbose & VERB_PROCESS )
		cout << "Setting to the smallest" << endl;
	
	double		v1;
	for ( long j=0; j<datasize; j++ ) {
		v1 = (*p)[j];
    	if ( (*this)[j] > v1 ) set(j, v1);
    }
	
	statistics();
}

/**
@brief 	Calculates the inverse tangent from two images.
@param 	*p		denominator (x) image.

**/
void		Bimage::arctangent(Bimage* p)
{
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::arctangent", __FILE__, __LINE__);
		return;
	}

    if ( verbose & VERB_PROCESS )
		cout << "Calculating the inverse tangent" << endl;
	
	double		v1;
	for ( long j=0; j<datasize; j++ ) {
		v1 = atan2((*this)[j], (*p)[j]);
    	set(j, angle_set_negPI_to_PI(v1));
    }
	
	statistics();
}

/**
@brief 	Adds two images together, adjusting for size difference.
@param 	&p			image to add.
@return Bimage* 	summed image, NULL on error.

	The second image is added to the first:
		image1 = image1 + image2*scale + shift
	The resultant image size is the bigger of the two.

**/
Bimage*		Bimage::operator+(Bimage& p)
{
	if ( n != p.n ) {
		error_show("Error in Bimage::operator+", __FILE__, __LINE__);
		cerr << "Different numbers of sub-images: " << n << " vs " << p.n << endl << endl;
		return this;
	}
	
	if ( c != p.c ) {
		error_show("Error in Bimage::operator+", __FILE__, __LINE__);
		cerr << "Different numbers of channels: " << c << " vs " << p.c << endl << endl;
		return this;
	}
	
	Vector3<long> 	sz(size().max(p.size()));
	
    long			i, j, xx, yy, zz, cc, nn;
    long     	    i1, i2, x1, x2, y1, y2, z1, z2;
	double			v;
	double			fill((avg + p.avg)/2);
	
	Bimage* 		pnu = new Bimage(Float, compound_type(), sz, n);
	pnu->sampling(sampling(0));

    if ( verbose & VERB_PROCESS ) {
		cout << "Summing two images:" << endl;
		cout << "New size:                       " << sz << endl;
	}

	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<pnu->z; zz++ ) {
			z1 = z2 = -1;
			if ( zz < z ) z1 = zz;
			if ( zz < p.z ) z2 = zz;
			for ( yy=0; yy<pnu->y; yy++ ) {
				y1 = y2 = -1;
				if ( z1 > -1 && yy < y ) y1 = yy;
				if ( z2 > -1 && yy < p.y ) y2 = yy;
				for ( xx=0; xx<pnu->x; xx++ ) {
					x1 = x2 = -1;
					if ( y1 > -1 && xx < x ) x1 = xx;
					if ( y2 > -1 && xx < p.x ) x2 = xx;
					for ( cc=0; cc<c; cc++, i++ ) {
						j = 0;
						v = (*pnu)[i];
						if ( x1 > -1 ) {
							i1 = index(cc,x1,y1,z1,nn);
							v += (*this)[i1];
							j++;
						}
						if ( x2 > -1 ) {
							i2 = p.index(cc,x2,y2,z2,nn);
							v += p[i2];
							j++;
						}
						if ( j ) v /= j;
						else v = fill;
						pnu->set(i, v);
					}
				}
			}
		}
		pnu->image[nn].background(background(nn) + p.background(nn));
	}
	
	pnu->statistics();

	return pnu;
}

/**
@brief 	Converts a one-dimensional image into a plot.
@return Bplot* 		plot structure.

Each image generates one plot with each channel converted to a curve.

**/
Bplot* 		Bimage::plot()
{
	if ( y*z > 1 ) {
		cerr << "Error in plot_from_image: Only one-dimensional images can be converted!" << endl;
		bexit(-1);
	}
	
	change_type(Float);
	
	long			i, j, nn, cc, xx;
	long			npage(n);
	long			nrow(x);
	long			ncol(n*(1 + c));
	double			dist, ux;
	Bstring			title(file_name());
	
	Bplot*			plot = new Bplot(npage, nrow, ncol);
	
	plot->title(title);

	for ( nn=0; nn<npage; nn++ ) {
		ux = sampling(nn)[0];
		title = file_name().c_str() + Bstring(nn+1, ": %d");
		plot->page(nn).title(title);
		plot->page(nn).columns(1+c);
		plot->page(nn).column(0).number(nn*(1+c));
		plot->page(nn).column(0).label("x");
		plot->page(nn).column(0).axis(1);
		for ( i=1, j=1+nn*(1+c); i<=c; i++, j++ ) {
			plot->page(nn).column(i).number(j);
			plot->page(nn).column(i).type(2);
			plot->page(nn).column(i).label("Value");
			plot->page(nn).column(i).axis(3);
		}
/*		if ( c == 2 ) {
			plot->page(nn).column(1).color(1,0,0);
			plot->page(nn).column(2).color(0,0,1);
		}
		if ( c > 2 ) {
			plot->page(nn).column(1).color(1,0,0);
			plot->page(nn).column(2).color(0,1,0);
			plot->page(nn).column(3).color(0,0,1);
		}*/
		for ( cc=0; cc<c; ++cc ) {
			dist = cc*2.0/(c-1);
			if ( dist < 1 ) {
				plot->page(nn).column(cc+1).color(1 - dist, dist, 0);
			} else {
				plot->page(nn).column(cc+1).color(0, 2 - dist, dist - 1);
			}
		}
		plot->page(nn).axis(1).min(-image[nn].origin()[0]*ux);
		plot->page(nn).axis(1).max((x-1-image[nn].origin()[0])*ux);
		plot->page(nn).axis(3).min(minimum());
		plot->page(nn).axis(3).max(maximum());
		for ( xx=0, i=nrow*nn*(1+c); xx<nrow; xx++, i++ )
			(*plot)[i] = (xx - image[nn].origin()[0])*ux;
		for ( cc=0; cc<c; cc++ )
			for ( xx=0, j=cc+nn*nrow*c; xx<nrow; xx++, i++, j+=c )
				(*plot)[i] = (*this)[j];
	}
	
	return plot;
}


/**
@brief A vector image is converted to a simple image.

	The new value for each voxel is the length of each vector.

**/
void		Bimage::vector_to_simple()
{
	if ( compoundtype == TSimple ) return;
	
	long			i, j, cc, ds(x*y*z*n);
	double			d;
	float*			nudata = new float[ds];
	
	if ( verbose & VERB_PROCESS )
		cout << "Converting vectors to lengths" << endl << endl;

	for ( i=j=0; i<ds; ++i ) {
		for ( cc=0, d=0; cc<c; ++cc, ++j )
			d += (*this)[j]*(*this)[j];
		nudata[i] = sqrt(d);
	}
	
	data_type(Float);
	compoundtype = TSimple;
	c = 1;
	
	data_assign((unsigned char *) nudata);
	
	statistics();
}


