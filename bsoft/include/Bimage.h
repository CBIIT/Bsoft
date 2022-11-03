/**
@file	Bimage.h
@brief	Header file for image class
@author Bernard Heymann
@date	Created: 19990321
@date 	Modified: 20220805
**/

//#include <time.h>

#include "json.h"
#include "string_util.h"
#include "Bstring.h"
#include "Complex.h"
#include "fft.h"
#include "FSI_Kernel.h"
#include "Matrix.h"
#include "Vector3.h"
#include "View.h"
#include "View2.h"
#include "Color.h"
#include "UnitCell.h"
#include "symmetry.h"
#include "ps_plot.h"
#include "Bsuperpixel.h"
#include "Bgraphseg.h"

#include <fstream>
#include <ctime>

#define NPOLANG	720

#ifndef _fouriertype_
/**
@enum 	FourierType
@brief 	Fourier transform format specifier.

	Transforms are classified according to where the origin is and whether
	all or only the hermitian half is stored in a file.
**/
enum FourierType {
	NoTransform = 0, 	// No transform
	Standard = 1,		// Standard transform: origin = (0,0,0)
	Centered = 2,		// Centered transform: origin = (nx/2,ny/2,nz/2)
	Hermitian = 3,		// Hermitian half: origin = (0,0,0)
	CentHerm = 4		// Centered hermitian: origin = (0,ny/2,nz/2)
} ;

#define _fouriertype_
#endif

#ifndef _Bimage_

union TypePointer {
	unsigned char*	uc;
	signed char*	sc;
	unsigned short*	us;
	short*			ss;
	unsigned int*	ui;
	int*			si;
	unsigned long*	ul;
	long*			sl;
	float*			f;
	double*			d;
} ;


/**
@brief General sub-image parameter class.

	This contains all the information pertinent to a single image in a
	multi-image file.
**/
class Bsub_image {
private:
	double		min, max;		// Image extremes
	double		avg, std;		// Average and standard deviation
	double		bkg;			// Image background
	double		ux, uy, uz; 	// Sampling (in Angstrom/pixel) 
	double		ox, oy, oz; 	// Origin (position of object origin in voxel coordinates)
	double		vx, vy, vz; 	// View orientation unit vector
	double		angle;			// Rotation around view vector in radians
	double		mag;			// Magnification
	double		fom;			// Figure-of-merit associated with the image
	long		sel;			// Selection
	void			check_sampling() {
		if ( !isfinite(ux) ) ux = 1;
		if ( !isfinite(uy) ) uy = 1;
		if ( !isfinite(uz) ) uz = 1;
		if ( ux <= 1e-10 ) ux = 1;
		if ( uy <= 1e-10 ) uy = 1;
		if ( uz <= 1e-10 ) uz = 1;
	}
	void		check_origin() {
		if ( !isfinite(ox) ) ox = 0;
		if ( !isfinite(oy) ) oy = 0;
		if ( !isfinite(oz) ) oz = 0;
	}
public:
	Bsub_image();
//	~Bsub_image();
	double		minimum() { return min; }
	double		maximum() { return max; }
	double		average() { return avg; }
	double		standard_deviation() { return std; }
	double		variance() { return std*std; }
	void		minimum(double d) { min = d; }
	void		maximum(double d) { max = d; }
	void		average(double d) { avg = d; }
	void		standard_deviation(double d) { std = d; }
	double		background() { return bkg; }
	void		background(double d) { bkg = d; }
	void		sampling(double x, double y, double z) { ux=x; uy=y; uz=z; check_sampling(); }
	void		sampling(Vector3<double> vec) { ux=vec[0]; uy=vec[1]; uz=vec[2]; check_sampling(); }
	void		sampling(vector<double> vec) {
		ux=uy=uz=vec[0];
		if ( vec.size() > 1 ) uy=vec[1];
		if ( vec.size() > 2 ) uz=vec[2];
		check_sampling();
	}
	Vector3<double> sampling() { return Vector3<double>(ux, uy, uz); }
	void		origin(double x, double y, double z) { ox=x; oy=y; oz=z; check_origin(); }
	void		origin(Vector3<double> vec) { ox=vec[0]; oy=vec[1]; oz=vec[2]; check_origin(); }
	void		origin(vector<double> vec) {
		ox=oy=oz=vec[0];
		if ( vec.size() > 1 ) oy=vec[1];
		if ( vec.size() > 2 ) oz=vec[2];
		check_origin();
	}
	Vector3<double> origin() { return Vector3<double>(ox, oy, oz); }
	Vector3<double>	real_coordinates(long ix, long iy, long iz) {
		return Vector3<double>(ux*(ix - ox), uy*(iy - oy), uz*(iz - oz));
	}
	Vector3<double>	real_coordinates(Vector3<long> ic) {
		return real_coordinates(ic[0], ic[1], ic[2]);
	}
	Vector3<long>	image_coordinates(double xx, double yy, double zz) {
		return Vector3<long>(xx/ux + ox, yy/uy + oy, zz/uz + oz);
	}
	Vector3<long>	image_coordinates(Vector3<double> rc) {
		return image_coordinates(rc[0], rc[1], rc[2]);
	}
	void		view(double x, double y, double z, double a) { vx=x; vy=y; vz=z; angle=a; }
	void		view(View vw) { vx=vw[0]; vy=vw[1]; vz=vw[2]; angle=vw.angle(); }
	template <typename T>
	void		view(View2<T> vw) { vx=vw[0]; vy=vw[1]; vz=vw[2]; angle=vw.angle(); }
	void		view(vector<double> vw) { vx=vw[0]; vy=vw[1]; vz=vw[2]; angle=vw[3]; }
	View		view() { return View(vx, vy, vz, angle); }
	void		view_angle(double d) { angle = d; }
	double		view_angle() { return angle; }
	double		magnification() { return mag; }
	void		magnification(double d) { mag = d; }
	double		FOM() { return fom; }
	void		FOM(double d) { fom = d; }
	long		select() { return sel; }
	void		select(long i) { sel = i; }
} ;

/**
@brief General image parameter class.

	All floating point coordinates are in angstroms and converted 
	using the voxel units parameters.
**/
class Bimage {
public:
	Bimage*			next;			// Pointer for linked lists
private:
	Bstring			id;				// Internal identifier
	long			c; 				// Number of channels
	long			x, y, z; 		// Dimensions, xyz
	long			n;				// Number of images
	long			px, py, pz; 	// Page dimensions
	long			sn, sz;			// Sub-image number and slice for display
	long			offset; 		// Data offset
	long			datasize;		// Number of data elements c*x*y*z*n
	DataType		datatype;		// Base data type
	double			dtmin;			// Data type minimum
	double			dtmax;			// Data type maximum
	CompoundType	compoundtype;	// Compound data type
	FourierType		fouriertype;  	// Transform type
	double			min, max;		// Limits
	double			avg, std;		// Average and standard deviation
	double			smin, smax; 	// Limits for display
	double			ss;				// Display scale
	UnitCell		ucell;			// Unit cell dimensions (angstrom) and angles (radian)
	TypePointer		d;				// Unioned data pointer
	JSvalue			metadata;		// Miscelaneous parameters/Metadata
public:
	Bsub_image*		image;			// Sub-images
private:
	void			initialize();
	void			initialize(long nc, long nx,
						long ny, long nz, long nn);
	void			initialize(CompoundType ctype, long nx,
						long ny, long nz, long nn);
	void			internal_copy(const Bimage& p);
public:
	// Construction and destruction
	Bimage();
	Bimage(const Bimage& p);
	Bimage(Bstring& fn, int readdata, int img_select);
	Bimage(DataType type, CompoundType ctype, long nx, long ny,
				long nz, long nn);
	Bimage(DataType type, CompoundType ctype, Vector3<long> size, long nn);
	Bimage(DataType type, CompoundType ctype, vector<long> size, long nn);
	Bimage(DataType type, long nc, long nx, long ny,
				long nz, long nn);
	Bimage(DataType type, long nc, Vector3<long> size, long nn);
	Bimage(Matrix& mat, long scale);
	~Bimage();
	CompoundType	guess_compoundtype(long nc);
	bool			check_compoundtype(long nc, CompoundType ct);
	void			check();
	bool			check_if_same_size(Bimage* p);
	bool			check_if_same_image_size(Bimage* p);
	void			check_sampling();
	void			check_resolution(double& resolution);
	bool 			compatible(Bimage* p);
	// Informational strings
	void			identifier(Bstring s) { id = s; }
	Bstring&		identifier() { return id; }
	// Metadata
	JSvalue&		meta_data() { return metadata; }
	JSvalue&		operator[](string tag) {
		return metadata[tag];
	}
	vector<JSvalue*>	query(string& path) {
		return metadata(path);
	}
	long			erase(string tag) { return metadata.erase(tag); }
	void			meta_data_update();
	void			update_from_meta_data();
	void			meta_data_retain_one_image(long img_num);
	void			file_name(string s) { metadata["filename"] = s; }
	string&			file_name() {
		if ( !metadata.exists("filename") ) {
			metadata["filename"] = "?";
			cerr << "Warning: No proper file name found!" << endl;
		}
		return metadata["filename"].value();
	}
	Bimage*			find(Bstring fn) {
		Bstring			b = fn.base();
		for ( Bimage* p = this; p; p = p->next ) {
			Bstring		f = p->file_name();
			f = f.base();
			if ( f.compare(b) == 0 ) return p;
		}
		return NULL;
	}
	void			label(string s) { metadata["label"] = s; }
	string&			label() {
		if ( !metadata.exists("label") ) metadata["label"] = "?";
		return metadata["label"].value();
	}
	// Time parameter management: Seconds since 00:00:00 January 1, 1970, (UTC)
	time_t			get_time() {
		if ( !metadata.exists("time") ) metadata["time"] = time(NULL);
		return metadata["time"].integer();
	}
	tm*				get_localtime() {
		time_t 			time_sec = get_time();
		return localtime(&time_sec);
	}
	void			set_time(time_t t) { metadata["time"] = t; }
	void			set_time(tm* t) { metadata["time"] = mktime(t); }
	// Data allocation and assignment
	long			alloc_size() const {
		return ( datatype == Bit )? (px/8)*y*z*n: c*x*y*z*n*data_type_size();
	}
	long			data_size() { datasize = c*x*y*z*n; return datasize; }
	unsigned char*	data_alloc() { return data_alloc(alloc_size()); }
	unsigned char*	data_alloc_and_clear() { data_alloc(); clear(); return d.uc; }
	unsigned char*	data_alloc(long nbytes);
	unsigned char*	data_alloc_and_clear(long nbytes) { data_alloc(nbytes); clear(); return d.uc; }
	unsigned char*	data_alloc(DataType type, CompoundType ctype, long nx, 
							   long ny, long nz, long nn);
	unsigned char*	data_alloc_and_clear(DataType type, CompoundType ctype, long nx, 
							   long ny, long nz, long nn) {
		data_alloc(type, ctype, nx, ny, nz, nn); 
		clear(); 
		return d.uc; 
	}
	unsigned char*	data_assign(unsigned char* nudata);
	unsigned char*	data_pointer() { return d.uc; }
	unsigned char*	data_pointer(long offset) { return &d.uc[offset*data_type_size()]; }
	void			data_pointer(unsigned char* ptr) { d.uc = ptr; }
	void			data_delete();
	long			data_offset() { return offset; }
	void			data_offset(long doff) { offset = doff; }
	long			image_size() { return x*y*z; }
	// Assignment
	Bimage&			operator=(const Bimage& p);
	// Element manipulations
	double			operator[](long j) const;
	double			get(long nn, long xx, long yy, long zz, long cc=0);
	double			get(long nn, Vector3<double> vox, long cc=0);
	vector<double>	values(long i) {
		long			cc;
		vector<double>	val(c);
		for ( i*=c, cc=0; cc<c; cc++, i++ ) val[cc] = (*this)[i];
		return val;
	}
	void			values(long i, vector<double> val) {
		long			cc;
		for ( i*=c, cc=0; cc<c; cc++, i++ ) val[cc] = (*this)[i];
	}
	vector<double>	values(long nn, Vector3<double> vox);
	Complex<double>	complex(long j);
	Vector3<double>	vector3(long j);
	RGB<double>		rgb(long j);
	RGBA<double>	rgba(long j);
	CMYK<double>	cmyk(long j);
	TypePointer		fill_value(double v);
	void			set(long j, double v);
	void			set(long j, Complex<double> cv);
	void			set(long j, RGB<double> color);
	void			set(long j, RGBA<double> color);
	void			set(long j, CMYK<double> color);
	void			set(long j, Vector3<double> vec);
	void			set(long j, View view);
	void			add(long j, double v) { set(j, (*this)[j] + v); }
	void			add(long j, Complex<double> cv) { set(j, complex(j) + cv); }
	void			add(double xx, double yy, double zz, long nn, double v);
	void			set_max(double xx, double yy, double zz, long nn, double v);
	void			multiply(long j, Complex<double> cv) { set(j, complex(j) * cv); }
	double			average2D(long cc, double xf, double yf, 
							double zf, long nn, double iscale);
	double			average(long cc, double xf, double yf, 
							double zf, long nn, double iscale);
	double			interpolate(long cc, double xx,
							double yy=0, double zz=0, long nn=0, double fill=0) const;
//	double			interpolate(double xx, double yy=0, double zz=0,
//							long nn=0, double fill=0) const {
//		return interpolate(0, xx, yy, zz, nn, fill);
//	}
	double			interpolate(long cc, Vector3<double> vec,
							long nn=0, double fill=0) const {
		return interpolate(cc, vec[0], vec[1], vec[2], nn, fill);
	}
	double			interpolate(Vector3<double> vec,
							long nn=0, double fill=0) const {
		return interpolate(0, vec[0], vec[1], vec[2], nn, fill);
	}
	double			interpolate(vector<double> vec,
							long nn=0, double fill=0) const {
		return interpolate(0, vec[0], vec[1], vec[2], nn, fill);
	}
	double			interpolate_wrap(long cc, double xx,
							double yy=0, double zz=0, long nn=0) const;
	double			interpolate_wrap(double xx, double yy=0, double zz=0,
							long nn=0) const {
		return interpolate_wrap(0, xx, yy, zz, nn);
	}
	double			interpolate_wrap(long cc,
							Vector3<double> vec, long nn=0) const {
		return interpolate_wrap(cc, vec[0], vec[1], vec[2], nn);
	}
	double			interpolate_wrap(Vector3<double> vec, long nn=0) const {
		return interpolate_wrap(0, vec[0], vec[1], vec[2], nn);
	}
	// Type management
	DataType		data_type() { return datatype; }
	void			data_type(DataType dt) { datatype = dt; dtmin = data_type_min(); dtmax = data_type_max(); }
	CompoundType	compound_type() { return compoundtype; }
	long			compound_type(CompoundType ct);
	long			data_type_bits() const;
	long			data_type_size() const;
	long			compound_type_size() { return c; }
	double			data_type_min();
	double			data_type_max();
	Bstring			data_type_string();
	Bstring			compound_type_string();
	void			fix_type();
	void			change_type(char letter);
	void			change_type(char* string);
	void			change_type(DataType nutype);
	Bimage*			split_channels();
	Bimage*			split_channels_to_images();
	void			combine_channels(long nc, CompoundType ct = TSimple);
	// Fourier transform management
	FourierType		fourier_type() { return fouriertype; }
	void			fourier_type(FourierType tf) { fouriertype = tf; }
	void			zero_fourier_origin();
	// Data parameters
	long			channels() { return c; }
	void			channels(long cc) { c = cc; }
	Vector3<long>	size() const { return Vector3<long>(x, y, z); }
	void			size(long nx, long ny, long nz) { x=nx; y=ny; z=nz; }
	void			size(Vector3<long> vec) { x=vec[0]; y=vec[1]; z=vec[2]; data_size(); }
	void			size(vector<long> vec) { x=vec[0]; y=vec[1]; z=vec[2]; data_size(); }
	void			sizeX(long nx) { x = nx; }
	void			sizeY(long ny) { y = ny; }
	void			sizeZ(long nz) { z = nz; }
	long			sizeX() const { return x; }
	long			sizeY() const { return y; }
	long			sizeZ() const { return z; }
	long			index(long nx, long ny) const {
		return ny*x + nx;
	}
	long			index(long nx, long ny, long nz) const {
		return (nz*y + ny)*x + nx;
	}
	long			index(long nx, long ny, long nz, long nn) const {
		return ((nn*z + nz)*y + ny)*x + nx;
	}
	long			index(long nc, long nx, long ny, long nz, long nn) const {
		return (((nn*z + nz)*y + ny)*x + nx)*c + nc;
	}
	long			index(Vector3<long> vox, long nn) const {
		return ((nn*z + vox[2])*y + vox[1])*x + vox[0];
	}
	long			index(vector<long> vox, long nn) const {
		return ((nn*z + vox[2])*y + vox[1])*x + vox[0];
	}
	long			index_wrap(long nx, long ny, long nz) const {
		while ( nx < 0 ) nx += x;
		while ( ny < 0 ) ny += y;
		while ( nz < 0 ) nz += z;
		while ( nx >= x ) nx -= x;
		while ( ny >= y ) ny -= y;
		while ( nz >= z ) nz -= z;
		return (nz*y + ny)*x + nx;
	}
	long			index_wrap(Vector3<long> coor) const {
		return index_wrap(coor[0], coor[1], coor[2]);
	}
	long			index_wrap(vector<long> coor) const {
		return index_wrap(coor[0], coor[1], coor[2]);
	}
	void			coordinates(long i, long &nx, long &ny, long &nz) {
			nx = i%x;
			i = (i - nx)/x;
			ny = i%y;
			i = (i - ny)/y;
			nz = i%z;
	}
	void			coordinates(long i, long &nx, long &ny, long &nz, long &nn) {
			nx = i%x;
			i = (i - nx)/x;
			ny = i%y;
			i = (i - ny)/y;
			nz = i%z;
			i = (i - nz)/z;
			nn = i%n;
	}
	void			coordinates(long i, long &nc, long &nx, long &ny, long &nz, long &nn) {
			nc = i%c;
			i = (i - nc)/c;
			nx = i%x;
			i = (i - nx)/x;
			ny = i%y;
			i = (i - ny)/y;
			nz = i%z;
			i = (i - nz)/z;
			nn = i%n;
	}
	Vector3<long>	coordinates(long i) {
		long	nc, nx, ny, nz, nn;
		coordinates(i, nc, nx, ny, nz, nn);
		return Vector3<long>(nx,ny,nz);
	}
	Vector3<double>	real_coordinates(long i) {
		return image->real_coordinates(coordinates(i));
	}
	template <typename T>
	bool			within_boundaries(Vector3<T> loc) {
		if ( loc[0] >= 0 && loc[0] < x && loc[1] >= 0 && loc[1] < y &&
			loc[2] >= 0 && loc[2] < z ) return 1;
		return 0;
	}
	bool			within_boundaries(long xx, long yy, long zz) {
		if ( xx >= 0 && xx < x && yy >= 0 && yy < y && zz >= 0 && zz < z ) return 1;
		return 0;
	}
	Vector3<long>	kernel_low(long i, long k=1) {
		Vector3<long>		lo(coordinates(i));
		for ( long j=0; j<3; j++ )
			lo[j] = ( lo[j] < k )? 0: lo[j] -= k;
		return lo;
	}
	Vector3<long>	kernel_high(long i, long k=1) {
		Vector3<long>		hi(coordinates(i)+k);
		return hi.min(size()-1);
	}
	Vector3<long>	kernel_low(long i, Vector3<long> k) {
		Vector3<long>		lo(coordinates(i));
		for ( long j=0; j<3; j++ )
			lo[j] = ( lo[j] < k[j] )? 0: lo[j] -= k[j];
		return lo;
	}
	Vector3<long>	kernel_high(long i, Vector3<long> k) {
		Vector3<long>		hi(coordinates(i)+k);
		return hi.min(size()-1);
	}
	long			kernel_min(long idx, long ksize);
	long			kernel_max(long idx, long ksize);
	double			kernel_average(long idx, long ksize, double tmin, double tmax);
	double			kernel_sum(long idx, long ksize);
	double			kernel_neighbor_average(long idx, long ksize);
	long			kernel_max_neigbor(long idx, long ksize);
	long			kernel_max_wrap(long idx, long ksize);
	multimap<double,long>	kernel_order(long idx, long ksize);
	multimap<double,long>	kernel_order_neighbors(long idx, long ksize);
	Vector3<long>	page_size() { return Vector3<long>(px, py, pz); }
	void			page_size(long nx, long ny, long nz) {
						px=nx; py=ny; pz=nz;
						if ( datatype==Bit ) px = 8*((x-1)/8 + 1);
	}
	void			page_size(Vector3<long> vec) { page_size(vec[0], vec[1], vec[2]); }
	void			page_size(vector<long> vec) { page_size(vec[0], vec[1], vec[2]); }
	Vector3<double>	real_size() {
		Vector3<double>	u = image->sampling();
		return Vector3<double>(u[0]*x, u[1]*y, u[2]*z);
	}
	long			images() { return n; }
	void			images(long nn);
	double			voxel_size() {
		Vector3<double>		u(image->sampling());
		double				vs(1);
		if ( x ) vs *= u[0];
		if ( y ) vs *= u[1];
		if ( z ) vs *= u[2];
		return vs;
	}
	void			sampling(long nn, Vector3<double> u) {
		image[nn].sampling(u);
	}
	void			sampling(long nn, double ux, double uy, double uz) {
		image[nn].sampling(ux, uy, uz);
	}
	void			sampling(Vector3<double> u) {
		for ( long nn=0; nn<n; ++nn ) image[nn].sampling(u);
	}
	void			sampling(double ux, double uy, double uz) {
		for ( long nn=0; nn<n; ++nn ) image[nn].sampling(ux, uy, uz);
	}
	Vector3<double>	sampling(long nn) {
		return image[nn].sampling();
	}
	vector<Vector3<double>>	sampling() {
		vector<Vector3<double>>	sam;
		for ( long nn=0; nn<n; ++nn ) sam.push_back(sampling(nn));
		return sam;
	}
	void			sampling(vector<Vector3<double>> sam) {
		if ( n != sam.size() )
			cerr << "Error in Bimage::sampling: number of images are not the same! (" 
				<< n << " != " << sam.size() << ")" << endl;
		for ( long nn=0; nn<n && nn<sam.size(); ++nn ) image[nn].sampling(sam[nn]);
	}
	// Statistics
	double			minimum() { return min; }
	double			maximum() { return max; }
	double			average() { return avg; }
	double			standard_deviation() { return std; }
	double			variance() { return std*std; }
	void			minimum(double d) { min = d; }
	void			maximum(double d) { max = d; }
	void			average(double d) { avg = d; }
	void			standard_deviation(double d) { std = d; }
	double			background(long nn) {
		return image[nn].background();
	}
	void			background(long nn, double bkg) {
		image[nn].background(bkg);
	}
	void			background(double bkg) {
		for ( long nn=0; nn<n; ++nn ) image[nn].background(bkg);
	}
	// Variables for image display
	long			show_image() { return sn; }
	void			show_image(long nn) { sn = nn; }
	long			show_slice() { return sz; }
	void			show_slice(long nz) { sz = nz; }
	double			show_scale() { return ss; }
	void			show_scale(double scale) { ss = scale; }
	double			show_minimum() { return smin; }
	double			show_maximum() { return smax; }
	void			show_minimum(double v) { smin = v; }
	void			show_maximum(double v) { smax = v; }
	// Geometry
	void			origin(double ox, double oy, double oz) {
		for ( long nn=0; nn<n; ++nn )
			image[nn].origin(ox, oy, oz);
	}
	void			origin(vector<double> vec) {
		for ( long nn=0; nn<n; ++nn )
			image[nn].origin(vec);
	}
	void			origin(Vector3<double> vec) {
		for ( long nn=0; nn<n; ++nn )
			image[nn].origin(vec);
	}
	void			origin(long nn, Vector3<double> vec) {
		image[nn].origin(vec);
	}
	void			origin(long nn, vector<double> ori) {
		while ( ori.size() < 3 ) ori.push_back(0);
		image[nn].origin(ori);
	}
	void			origin(long nn, double ox, double oy, double oz) {
		image[nn].origin(ox, oy, oz);
	}
	Vector3<double>	default_origin() {
		Vector3<double>		ori(double(x/2), double(y/2), double(z/2));
		return ori;
	}
	void			view(double vx, double vy, double vz, double va) { for ( long j=0; j<n; j++ ) image[j].view(vx, vy, vz, va); }
	void			view(View vw) { for ( long j=0; j<n; j++ ) image[j].view(vw); }
	long			space_group() {
		if ( !metadata.exists("spacegroup") ) metadata["spacegroup"] = 1;
		return metadata["spacegroup"].integer();
	}
	void			space_group(unsigned int grp) { metadata["spacegroup"] = grp; }
	string&			symmetry() {
		if ( !metadata.exists("pointgroup") ) metadata["pointgroup"] = "C1";
		return metadata["pointgroup"].value();
	}
	void			symmetry(string grp) { metadata["pointgroup"] = grp; }
	void			unit_cell(UnitCell uc) {
		Vector3<double>		u(image->sampling());
		ucell = uc;
		if ( !isfinite(ucell.a()) || ucell.a() < u[0]*x ) ucell.a(u[0]*x);
		if ( !isfinite(ucell.b()) || ucell.b() < u[1]*y ) ucell.b(u[1]*y);
		if ( !isfinite(ucell.c()) || ucell.c() < u[2]*z ) ucell.c(u[2]*z);
	}
	UnitCell		unit_cell() { return ucell; }
	double			maximum_included_radius();
	// Data organization changes
	int				slices_to_images();
	int				images_to_slices();
	int 			channels_to_images();
	int 			images_to_channels(long nc, CompoundType ct);
	long			set_subset_selection(Bstring list);
	long			delete_images(Bstring list, int retain=0);
	long			select_images(Bstring list);
	// Data I/O
	unsigned char*	read_data(ifstream* fimg, int img_select, int sb, int vax, long pad);
	int				write(Bstring& fn);
	int				unpack_transform(unsigned char* data, FourierType tf);
	int				unpack_transform(int img_select, unsigned char* data, FourierType tf);
	int				pack_transform(unsigned char* data, FourierType tf);
	int				pack_transform(int img_select, unsigned char* data, FourierType tf);
	// Data parameter information
	long			statistics();
	long			statistics(long img_num);
	long			statistics(Bimage* pmask, double& regavg, double& regstd);
	double			poisson_statistics_check();
	long 			stats_within_radii(long nn, Vector3<double> loc,
						double rad_min, double rad_max, double& vavg, double& vstd);
	long			stats_in_shape(long nn, int type, Vector3<long> start,
						Vector3<long> end, double& vavg, double& vstd);
	long			stats_in_poly(long nn, int nvert, Vector3<double>* poly,
						double& vavg, double& vstd);
	long			stats_in_mask(long nn, Bimage* pmask);
private:
	int				kernel_sums(long n, long i, long ik, long nk);
	int				line_sums(long n, long ik, long nk);
public:
	int				variance(long kernel_size, int flag=0);
	int				variance(Vector3<long> kernel_size, int flag=0);
private:
	int				kernel_sums(long nn, long i, Bimage* pweight);
public:
	int				variance(Bimage* pweight);
	int				information();
	int				subimage_information();
	int				moments(long max_order);
	int				moments(long max_order, long nn);
	void			get(Bstring tag);
	// New image generation and obtaining data from other images
	Bimage*			copy();
	Bimage*			copy(long nu_nimg);
	Bimage*			copy_header() { return copy_header(n); }
	Bimage*			copy_header(long nu_nimg);
	Bimage*			extract(long nn);
	Bimage*			extract(long n1, long n2);
	Bimage*			extract(long nn, Vector3<long> coords, Vector3<long> size,
						int fill_type=0, double fill=0);
	Bimage*			extract(long nn, Vector3<double> loc, Vector3<long> size,
						Vector3<double> origin) {
		Matrix3		mat(1);
		return extract(nn, loc, size, origin, mat);
	}
	Bimage*			extract(long nn, Vector3<double> loc, Vector3<long> size,
						Vector3<double> origin, Matrix3 mat);
	Bimage*			extract(long nn, Vector3<double> loc, Vector3<long> size,
						Vector3<double> origin, View view) {
		return extract(nn, loc, size, origin, view.matrix());
	}
	Bimage*			extract_wrap(long nn, Vector3<double> loc, Vector3<long> size,
						Vector3<double> origin, Matrix3 mat);
	Bimage*			extract_wrap(long nn, Vector3<long> size, Matrix3 mat) {
		Vector3<double>		v;
		return extract_wrap(nn, v, size, v, mat);
	}
	Bimage*			extract_shell(long nn, double minrad, double maxrad);
	vector<Vector3<long>>	tile_coordinates(Vector3<long>& start,
						Vector3<long>& region, Vector3<long>& tile_size,
						Vector3<long>& step_size, int exceed);
	vector<Vector3<long>>	tile_coordinates(Vector3<long> tile_size, Vector3<long>& step_size);
	Bimage*			extract_tiles(long nn, vector<Vector3<long>>& coords, Vector3<long> tile_size);
	Bimage*			extract_tiles(long nn, Vector3<long> start,
						Vector3<long> region, Vector3<long> tile_size,
						Vector3<long> step_size, int exceed);
	Bimage*			extract_tiles(long nn, Vector3<long> tile_size, double fraction=0.2);
	Bimage*			extract_tile_stack(Vector3<long> coords, Vector3<long> tile_size, int fill_type=0, double fill=0);
	Bimage**		extract_tile_stacks(vector<Vector3<long>>& coords, Vector3<long> tile_size);
	Bimage*			extract_line(long nn, Vector3<double> start, Vector3<double> end, long width);
	Bimage*			extract_tetrahedron(Vector3<double>* tet, int fill_type=0, double fill=0);
	Bimage*			orthogonal_slices(long nn, Vector3<long> voxel,
						Vector3<long> ext_size);
	Bimage*			orthogonal_montage(Vector3<long> voxel, Vector3<long> ext_size, int pad=0, int fill_type=0, double fill=0);
	int				extract_show_chunk(Bimage* pshow, int aflag, long i, long len);
	Bimage*			extract_show(int aflag);
	Bimage*			extract_magnify(long nn, Vector3<long> center,
						Vector3<long> ext_size, double scale);
	Bimage*			extract_slice(long nz);
	Bimage*			extract_filament(long img_num, double width,
						int axis, long nspline, Vector3<double>* spline);
	int				replace(Bimage* img);
	int				replace(long nn, Bimage* img, long nr=0);
	int				replace(long nn, Bimage* img, long nr, double fill);
	// Whole image manipulations
	void			clear() { data_size(); for ( long j=0; j<datasize; j++ ) set(j, 0); }
//	void			clear() { for ( long j=0; j<alloc_size(); ++j ) d.uc[j] = 0; }
	void			fill(double v) {
		data_size(); for ( long j=0; j<datasize; j++ ) set(j, v);
	}
	// Density methods
	double			density(long nn, Vector3<double> coord, double radius, double& sigma);
	double			density(long nn, Vector3<double> coord, double radius) {
		double		sigma(0);
		return density(nn, coord, radius, sigma);
	}
	double			relative_density(Bimage* pmask);
	void			invert();
	void			reslice(const char* order) { Bstring s(order); reslice(s); s=0; }
	void			reslice(Bstring order);
	void			absolute();
	void			add(double v);
	void			phase_add(double v);
	void			multiply(double v);
	void			multiply(long nn, double v);
	void			power(double v);
	void			sine();
	void			arcsine();
	void			cosine();
	void			arccosine();
	void			tangent();
	void			arctangent();
	void			sum_images();
	void			average_images() { double w(1.0/n); sum_images(); multiply(w); }
	Bimage*			average_images(bool sd);
	Bimage*			moving_sum(long window, long step=1, int flag=0);
	void			progressive_sum();
	void			add(Bimage* p);
	void			add(long nn, Bimage* p);
	void			subtract(Bimage* p);
	void			add(Bimage* p, double scale, double shift);
	void			add(long nn, Bimage* p, double scale, double shift);
	void			multiply(Bimage* p, double scale, double shift=0);
	void			multiply(long nn, Bimage* p);
	void 			multiply(Bimage* p);
	void			divide(Bimage* p, double scale=1, double shift=0);
	void			divide_one(Bimage* p, double scale=1, double shift=0);
	void			inverse(double minval=0);
	void			largest(Bimage* p);
	void			smallest(Bimage* p);
	void			arctangent(Bimage* p);
	Bimage*			operator+(Bimage& p);
	Bplot* 			plot();
	void			vector_to_simple();
	void			sum(long m, Bimage** p);
	void			catenate(long m, Bimage** p);
	Bimage* 		blend(Bimage* p, long number);
	int 			place(long nn, Bimage* p, Vector3<double> loc,
						double radius=0, double scale=1, double shift=0, int operation=0);
	int				place_with_addition(Bimage* p, long nn);
	int				place_with_overlap(Bimage* p, long nn);
	int				place_central_part(Bimage* p, long nn);
	int				assemble_tiles(Bimage* pt, int flag=0);
	double			linear_fit(Bimage* p, Bimage* pmask, double max_exclude);
	int 			histomatch(Bimage* p, long bins);
	int				replace_half(Bimage* p);
	// Filter methods
	int				kernel_gaussian(double sigma, double max);
	int				kernel_laplacian_of_gaussian(double sigma, double max);
private:
	int				average_line_sum(long n, long i, long ik, long nk);
	int				average_line_sums(long n, long ik, long nk);
public:
	int				filter_average(long kernel_size);
	int				filter_average(Vector3<long> k);
	Bimage*			gradient();
private:
	Vector3<double>	gradient_voxel(long i);
public:
	Bimage*			gradient3x3();
private:
	int				gaussian_line_sum(long nn, long i, long ik, Bimage* pk);
	int				gaussian_line_sums(long nn, long ik, Bimage* pk);
public:
	int				filter_gaussian(long kernel_size, double sigma=0);
	int				filter_sinc();
	int				convolve_chunk(Bimage* pkernel, float* nudata, long i, long len);
	int				convolve(Bimage* pkernel);
	int 			filter_ortho(int type);
	int				filter_dog(double sigma1, double sigma2);
	int				filter_bilateral_chunk(Bimage* pkernel, double sigma2,
						int kernel_type, float* nudata, long i, long len, int first);
	int 			filter_bilateral(double sigma1, double sigma2,
						int kernel_type, long kernel_radius);
	int 			filter_rolling_ball(long radius, double scale);
	int				filter_rank_chunk(long kernel_size, double rank, float* nudata, long i, long len);
	int				filter_rank(long kernel_size, double rank);
	Bimage*			filter_peak(long kernel_size);
	Bimage*			periodic_averaging(Vector3<double> period);
private:
	double			aniso_voxel(Bimage* pg, long i, long ksize, double w);
public:
	Bimage*			aniso_average(long ksize, double w);
	int				filter_by_difference(Bimage* p);
private:
	Bimage*			structure_tensor_2D(long nn=0);
	Bimage*			structure_tensor(long zs, long zw);
	int				diffusion_tensor_2D(double lambda, double C, double alpha);
	int				diffusion_tensor(double lambda, double C, double alpha);
	int				diffuse_2D(Bimage* pdn, Bimage* pd, double ht, long nn=0);
	int				diffuse(Bimage* pdn, Bimage* pd, double ht, long zs);
	int				nad_chunk_2D(Bimage* pdn, double lambda, double C,
						double alpha, double ht, long nn=0);
	int				nad_chunk(Bimage* pdn, double lambda, double C,
						double alpha, double ht, long zs, long zw);
public:
	Bimage*			nad_2D(double ht, double lambda, double C, double alpha);
	Bimage*			nad(double ht, long zw, double lambda, double C, double alpha);
	// Complex type methods
	void			simple_to_complex();
	void			two_to_complex();
	void			multi_channel_to_complex();
	Bimage*			complex_split();
	void			phase_to_complex();
	void			complex_to_real();
	void			complex_to_imaginary();
	void			complex_to_intensities();
	void			complex_to_amplitudes();
	void			complex_to_signed_amplitudes();
	void			complex_to_phases();
	void 			complex_conjugate();
	double			complex_power();
	double			complex_normalize();
	int 			complex_invert();
	int 			complex_convert(ComplexConversion conv);
	int				phase_shift(Vector3<double> shift);
	int 			phase_shift(long nn, Vector3<double> shift);
	int				phase_shift_to_origin();
	int				phase_shift_to_center();
	int 			complex_multiply(Bimage* p);
	int 			complex_product(Bimage* p);
	int 			complex_conjugate_product(Bimage* p, int norm=0);
	Bimage*			complex_conjugate_product_one2many(Bimage* p);
	int 			complex_apply_mask(Bimage* pmask);
	int 			complex_apply_negative_mask(Bimage* pmask);
	int 			complex_apply_dual_mask(Bimage* pmask);
	int 			complex_bandpass(double hires, double lores);
	Bimage*			pack_two_in_complex(Bimage* p);
	Bimage* 		unpack_combined_transform();
	int 			combined_complex_product();
	int 			combined_complex_product(Bimage* pmask);
//	int 			combined_complex_product(double hires, double lores);
	int 			combined_complex_product(double hires, double lores, Bimage* pmask=NULL);
	int				combined_phase_product(double hires, double lores, Bimage* pmask=NULL);
	int 			combined_complex_product_implicit_mask(double hires, double lores);
	double			merge_amplitudes_and_phases(Bimage* pamp);
	double			merge_amplitudes_and_phases(Bimage* pref, double res_hi, double res_lo);
	Bimage*			intensities_phase_colored(double scale);
	int				phase_colour_wheel();
	// FFT methods
	fft_plan		fft_setup(fft_direction dir, int opt=0);
	int 			fft(fft_direction dir, int norm_flag, ComplexConversion conv);
	int 			fft(fft_direction dir, int norm_flag);
	int 			fft(fft_plan plan, int norm_flag=1);
	int 			fft() { return fft(FFTW_FORWARD, 1); }
	int 			fft_back() {
		fft(FFTW_BACKWARD, 1, Real);
//		fourier_type(NoTransform);
//		complex_to_real();
		return 0;
	}
	int 			fft_back(fft_plan plan, int norm_flag=1) {
		fft(plan, norm_flag);
		fourier_type(NoTransform);
		complex_to_real();
		return 0;
	}
	int				fft(fft_direction dir, Vector3<long> tile_size, int norm_flag=1);
	int				fftz(fft_direction dir, int norm_flag=1);
	int 			fftz() { return fftz(FFTW_FORWARD, 1); }
	Vector3<double>	change_transform_size(Vector3<long> nusize);
	// Power spectrum methods
	int				power_spectrum(int flags=0);
	Bimage*			powerspectrum_tiled(long img_num, Vector3<long> tile_size, int flags=0);
	Bimage*			powerspectrum_tiled_exact(long img_num, Vector3<long> tile_size, int flags=0);
	Bimage*			powerspectrum_tilt_axis(long img_num, Vector3<long> tile_size,
						double tilt_axis, double tilt_offset, int flags=0);
	Bimage*			defocus_scale(long nn, double df, double df2, double iCL2, int fill_type);
	Bimage*			powerspectrum_tilted(long img_num, Vector3<long> tile_size,
						double tilt_axis, double tilt_angle, double defocus, double iCL2, int flags=0);
	Bimage*			powerspectrum_tiled_and_tilted(Vector3<long> tile_size,
						double tilt_axis, double tilt_angle, double tilt_offset, 
						double defocus, double iCL2, int flags=0);
	vector<double>	powerspectrum_isotropy(long n, double& lores, double& hires);
	// Frequency space interpolation and reconstruction
	long			fspace_maximum_radius(double resolution, double sampling_ratio=1);
	int				fspace_background();
	Complex<double>	fspace_interpolate(long img_num, Vector3<double> m, FSI_Kernel* kernel);
	int				fspace_2D_interpolate(Complex<float> cv, Vector3<double> m,
						double part_weight, int interp_type);
	int				fspace_pack_2D(Bimage* p, Matrix3 mat, double hi_res, double lo_res,
						Vector3<double> scale, double ewald_wavelength=0,
						double part_weight=1, int interp_type=0);
	int				fspace_pack_2D(Bimage* p, View asu_view, Bsymmetry& sym, double hi_res,
						double lo_res, Vector3<double> scale, double ewald_wavelength=0,
						double part_weight=1, int interp_type=0);
	long			fspace_pack_2D_into_central_section(Bimage* p,
						long ft_size, double scale, double hi_res, double lo_res, 
						Matrix3 matr, Matrix3 mat);
	int				fspace_pack_3D(Bimage* p, double hi_res=0, double threshold=0);
	long			fspace_reconstruction_add(Bimage* p);
	long			fspace_reconstruction_weigh();
	int 			fspace_reconstruction_stats(double resolution, double sampling_ratio=1);
	long			fspace_reconstruction_snr();
	int				fspace_translate(Vector3<double> shift);
	int				fspace_translate(long nn, Vector3<double> shift);
	int				fspace_resize(double scale, double res_hi, double res_lo);
	Bimage*			fspace_resize(Bimage* pref);
	// Amplitude weighting methods
	int				fspace_amp_one();
	int				fspace_amp_threshold(double threshold);
	int				fspace_sqrt_amp();
	int				fspace_square_amp();
	int				fspace_bandpass(double res_hi, double res_lo=0, double width=0);
	int				fspace_bandpass(double res_hi, double res_lo, double width, fft_plan planf, fft_plan planb);
	int				fspace_frequency_filter(double freq, double sigma);
	int				fspace_frequency_filter(double freq, double sigma, fft_plan planf, fft_plan planb);
	int				fspace_gabor_filter(Vector3<double> freq, double fsigma, double psigma);
	int				fspace_gabor_filter(Vector3<double> freq, double fsigma, double psigma, fft_plan planf, fft_plan planb);
	Bimage* 		fspace_radial_power(double resolution, double sampling_ratio=1);
	vector<double>	fspace_radial(long nn, long maxrad, int flag=0);
	int 			fspace_weigh(vector<double>& curve, int flag=0);
	int 			fspace_scale(vector<double>& scale, Bimage* pmask=NULL);
	int 			fspace_scale(long nn, vector<double>& scale, Bimage* pmask=NULL);
	double			fspace_fit_B_factor(double res_hi=0);
	int 			fspace_weigh_ramp(double resolution, fft_plan planf, fft_plan planb);
	int 			fspace_weigh_ramp(double resolution, double axis, fft_plan planf, fft_plan planb);
	int 			fspace_weigh_B_factor(double B, double resolution=0);
	int 			fspace_butterworth_band(double res_hi, double res_lo, int order=16);
	int 			fspace_weigh_C_curve(double resolution=0);
	int 			fspace_weigh_LoG(double resolution, double sigma);
	int 			fspace_weigh_RPS_curve(Bplot* plot, double resolution=0);
	int 			fspace_weigh_FSC_curve(Bplot* plot, double resolution=0);
	int 			fspace_weigh_gaussian(long nn, Vector3<double> sigma, int dir=0);
	Bimage*			fspace_gradient(Vector3<double> sigma);
	int 			fspace_weigh(Bimage* pref, Bimage* pmask, double resolution=0);
	int 			fspace_weigh_dose(long nn, double dose_per_frame, vector<double> critdose);
	int 			fspace_weigh_dose(double dose_per_frame, int flag=0);
	int				fspace_weigh_accumulated_dose(vector<double> dose);
	int 			fspace_normalize();
	int 			fspace_normalize_radial(Bimage* pmask, double resolution=0, int flag=0);
	int 			fspace_positive();
	// Friedel symmetry
	double			friedel_check();
	double			friedel_difference();
	int 			friedel_apply();
	// Projection methods
	Bimage*			project(char axis, int flags=1);
	Bimage*			rotate_project(Matrix3 mat, Vector3<double> translate,
						double radial_cutoff, int norm_flag=1);
	Bimage* 		project(View* view, int norm_flag=1);
	Bimage*			central_section(Matrix3 mat, double resolution, FSI_Kernel* kernel, double wavelength=0);
	Bimage*     	project(View* view, double resolution, FSI_Kernel* kernel,
						double wavelength=0, bool back=1, ComplexConversion conv = NoConversion);
	int 			back_project(Bimage* p, double resolution, double axis, 
						fft_plan planf, fft_plan planb);
	int				opposite_ewald();
	int				combine_ewald();
	// Resolution measures
	Bimage*			resolution_prepare(Bimage* p);
	Bimage*			resolution_prepare(Bimage* p, fft_plan plan);
	Bplot*			fsc_dpr(double hi_res, double sampling_ratio=1, int flag=0);
	Bplot*			fsc(double hi_res, double sampling_ratio, vector<double>& fsccut);
	Bplot*			fsc(Bimage* p, double hi_res, double sampling_ratio=1);
	Bimage*			fsc_shell(Bimage* p, double hi_res, double* cutoff, 
						int thickness, int step, int minrad, int maxrad, 
						int pad=1, int smooth=0, double fill=0);
private:
	double*			fsc_local_voxel(Bimage* p, double resolution, int size,
						int pad, int taper, double* cutoff, fft_plan plan, long i);
	double			local_filter_voxel(Bimage* resmap, long size,
						fft_plan planf, fft_plan planb, long i);
public:
	Bimage*			fsc_local(Bimage* p, Bimage* pmask, double resolution, double* cutoff,
						int mask_level, int size, int pad, Vector3<long> vedge, 
						int step=1, int taper=1, double fill=0);
	Bimage*			local_filter(Bimage* pmask, int mask_level, Bimage* resmap, 
						int size, Vector3<long> vedge);
	// Phase difference
	Bimage* 		phase_difference(Bimage* p, int type=0, double res_hi=0, double res_lo=0);
	double	 		average_phase_difference(Bimage* p, double res_hi, double res_lo, int weighting=1);
	int				phase_flip(Bimage* pd);
	int				ewald_sphere(double volt, double t);
	// Correlation methods
	double			correlate(Bimage* p);
	double			correlate(Bimage* p, double rmin, double rmax, Bimage* pmask=NULL, int flag=0);
	double			rotate_correlate(Vector3<double> axis, double angle);
	double			R_factor(Bimage* p);
	int 			auto_correlate(double hires, double lores);
	Bimage* 		cross_correlate(Bimage* p, double hires, double lores,
						fft_plan planf, fft_plan planb);
	Bimage* 		cross_correlate(Bimage* p, Bimage* pmask=NULL) {
		return cross_correlate(p, 0, 0, pmask);
	}
	Bimage* 		cross_correlate(Bimage* p, double hires, double lores, Bimage* pmask=NULL);
/*	Bimage* 		cross_correlate(Bimage* p, double hires, double lores,
						fft_plan planf, fft_plan planb) {
		return cross_correlate(p, hires, lores, NULL, planf, planb);
	}*/
	Bimage* 		phase_correlate(Bimage* p, double hires, double lores, Bimage* pmask=NULL);
	double			correlation_coefficient(Vector3<double> shift);
	Vector3<double>	find_shift_in_transform(double shift_limit);
	Bimage* 		cross_correlate_fspace(Bimage* p, double hires, double lores, double shift_limit);
	Bimage* 		cross_correlate(Bimage* p, double hires, double lores,
						Bimage* pmask, fft_plan planf, fft_plan planb);
	Bimage* 		cross_correlate_two_way(Bimage* p, double hires,
						double lores, fft_plan planf, fft_plan planb);
	Bimage* 		cross_correlate_validate(Bimage* p, Bimage* pmask);
	Vector3<double>	rotate_cross_correlate(Bimage* pref, View view,
						double hires, double lores, double search_radius, Bimage* pmask,
						double& cc, fft_plan planf, fft_plan planb);
	double			rotate_cross_correlate_two_way(Bimage* pref,
						double angle, double res_hi, double res_lo, double shift_limit,
						fft_plan planf, fft_plan planb);
	int				find_shift(Bimage* pref, Bimage* pmask, double hires,
						double lores, double radius, double sigma, int refine_flag);
	Vector3<double>	find_shift(Bimage* pref, Bimage* pmask, double hires,
						double lores, double radius, double sigma, int refine_flag, double& cc);
	Vector3<double>	find_shift(Bimage* pref, double hires,
						double lores, double radius, double sigma, int refine_flag,
						fft_plan planf, fft_plan planb);
	Vector3<double>	find_shift(Bimage* pref, Bimage* pmask, double hires,
						double lores, double radius, double sigma, int refine_flag,
						fft_plan planf, fft_plan planb, double& cc);
	Vector3<double>	find_shift(long nn, Bimage* pref, Bimage* pmask,
						double hi_res, double lo_res, double shift_limit, 
						fft_plan planf, fft_plan planb);
	Bimage*			find_template(Bimage* ptemp, Bimage* pmask,
						double hires, double lores, int bin, fft_plan planf, fft_plan planb);
	int				find_center(Bimage* pmask, double hires,
						double lores, double radius, double sigma, int refine_flag);
	Vector3<double>	rotate_find_shift(Matrix3 mat,
						double hires, double lores, double radius, double sigma,
						int refine_flag, fft_plan planf, fft_plan planb, double& cc);
	int				find_peak(double radius=1e30, double sigma=0);
	Vector3<double>	fit_peak();
	int				refine_peak_new();
	int				refine_peak();
	int				refine_peak(long kernel_size);
	Bimage*			find_peaks(long kernelsize);
	Vector3<double>*	find_peaks(double excl_dist, long& ncoor, double& threshold_min, 
						double& threshold_max, double pix_min=2, double pix_max=10);
	double			ccmap_confidence(long nn);
	double			peak_sigma(long nn, Vector3<long> coor, long kernel_size);
	double			search_views(Bimage* ptemp, View* view,
						double hires, double lores, double search_radius,
						Bimage* pmask, View& currview, Vector3<double>& currshift);
	double			search_volume_view(Bimage* ptemp, View view,
						double hires, double lores, Bimage* pmask, 
						double threshold, Bimage* pfit);
	Bimage*			search_volume(Bimage* ptemp, View* view, double alpha, 
						double alpha_step, double hires, double lores, 
						Bimage* pmask, double threshold);
	// Alignment methods
	Vector3<double>	find_shift_in_transform(long nn, Bimage* pref, double shift_limit);
	Bimage*			align_progressive_fast(long nref, double shift_limit);
	Bimage*			align_progressive(long nref, Bimage* pmask, 
						double hi_res, double lo_res, double shift_limit,
						fft_plan planf, fft_plan planb);
	Bimage*			align_local(long nref, Bimage* pmask,
						double hi_res, double lo_res, double shift_limit,
						fft_plan planf, fft_plan planb);
	vector<Vector3<double>>	align(long ref_num, long window, long step, Bimage* pmask,
						double hi_res, double lo_res, double shift_limit, 
						double edge_width, double gauss_width, Vector3<long> bin, int mode=0);
	JSvalue			align_fast(long ref_num, Bimage* pmask,
						double hi_res, double lo_res, double shift_limit, 
						double edge_width, double gauss_width);
	Bimage*			fspace_sum(int shift=0);
	Bimage*			fspace_shift_sum() { return fspace_sum(1); }
	Bimage*			fspace_subset_sums(int subset, int flag=0);
	Bplot*			fspace_ssnr(long nimg, double res_hi, double sampling_ratio);
	Bplot*			fspace_subset_ssnr(int subset, double res_hi, double sampling_ratio, int flag=0);
	double			correlate_annuli(Bimage* polref,
						int ann_min, int ann_max, double ang_min, double ang_max,
						fft_plan planf, fft_plan planb, double& cc_max);
	vector<double>	pps_angular_correlation(Bimage* pref,
						double res_hi, double res_lo, long nang, fft_plan planf);
	double			align2D_pps(Bimage* pref, 
						double res_hi, double res_lo, double shift_limit, double angle_limit, 
						fft_plan planf, fft_plan planb);
	double			align2D(Bimage* pref, double res_polar,
						int ann_min, int ann_max, Bimage* prs_mask, double shift_limit, double angle_limit, 
						fft_plan planf_1D, fft_plan planb_1D, fft_plan planf_2D, fft_plan planb_2D);
	int				align2D(Bimage* pref, int ann_min, int ann_max,
						double res_lo, double res_hi, double shift_limit, double angle_limit);
	// Color type methods
	void			simple_to_rgb();
	void			simple_to_rgba();
	void			color_to_simple();
	void			rgb_to_rgba();
	void			rgba_to_rgb();
	void			rgb_to_cmyk();
	void			cmyk_to_rgb();
	int 			one_color(int color, double cmin, double cmax, int flag=0);
	int 			color_red(double cmin, double cmax, int flag=0);
	int 			color_green(double cmin, double cmax, int flag=0);
	int 			color_blue(double cmin, double cmax, int flag=0);
	int 			pure_color();
	int 			color_combine(Bimage* p);
	Bimage*			color_spectrum(double cmin, double cmax);
	Bimage*			red_white_blue(double red_min, double white_min,
						double white_max, double blue_max);
	// Rescaling methods
	int				rescale(double scale, double shift);
	int				rescale(long nn, double scale, double shift);
	int 			rescale_to_min_max(double numin, double numax);
	int 			rescale_to_min_max(long nn, double numin, double numax);
	int 			rescale_to_avg_std(double nuavg, double nustd);
	int 			rescale_to_avg_std(long nn, double nuavg, double nustd);
	int 			rescale_to_avg_std(double nuavg, double nustd, Bimage* pmask);
	int 			truncate(double minim, double maxim, double setmin, double setmax);
	int 			truncate_to_min_max(double minim, double maxim);
	int 			truncate_to_avg(double minim, double maxim);
	int 			truncate_to_background(double minim, double maxim);
	int 			limit_levels(int nlevels);
	int				normalize(long imgnum, double average, double stdev, int norm_type, long bins);
	int				normalize(double average, double stdev, int norm_type);
	int				normalize_local(long kernel_size);
	int				normalize_local(Vector3<long> kernel);
	void			square();
	void			square_root();
	void			logarithm();
	void			exponential();
	int 			gradient_correction();
	int				quadric_correct(vector<double> param);
	vector<double>	quadric_fit();
	Bimage*			thickness(double reference, double emfp);
	// Background methods
	int				calculate_background(long nn, int flag);
	int				calculate_background(int flag=0);
	int				calculate_background(Bimage* pmask, long nn, int flag=0);
	int				calculate_background(Bimage* pmask, int flag=0);
	int				correct_background(long nn, int flag);
	int				correct_background(int flag=0);
	int				correct_background(Bimage* pmask, int flag=0);
	int				subtract_background();
	int				shift_background(double bkg);
	// Histogram methods
	vector<long>	histogram(long bins, double& scale, double& offset);
	Bplot* 			histogram(long bins);
	Bplot*			histogram_counts(int flags=0);
	Bplot*			percentiles();
	int				histogram_minmax(double& tmin, double& tmax);
	Bplot*			histogram_otsu_variance(long bins);
	double			otsu_threshold(long bins);
	vector<double>	otsu_variance(vector<long> h);
	vector<double>	histogram_multi_thresholds(long bins, long number);
//	vector<double>	histogram_multi_thresholds(vector<long> h, long number);
	vector<double>	histogram_gauss_fit(long bins, long ngauss=1);
	vector<double>	histogram_gauss_fit2(long bins, long ngauss=1);
	Bplot* 			histogram_gauss_plot(long bins, long ngauss=1);
	Bplot* 			histogram_poisson_fit(long bins, int flag=0);
	// Resizing methods
	int				resize(Vector3<long> nusize, Vector3<long> translate,
						int fill_type=0, double fill=0);
	Bimage*			resize_copy(Vector3<long> nusize, Vector3<long> translate,
						int fill_type=0, double fill=0);
	int				resize_wrap(Vector3<long> nusize, Vector3<long> translate);
	Bimage*			resize_wrap_copy(Vector3<long> nusize, Vector3<long> translate);
	int				pad(long sz, int fill_type=0, double fill=0);
	int				pad(Vector3<long> sz, int fill_type=0, double fill=0);
	Bimage*			pad_copy(long sz, int fill_type=0, double fill=0);
	Bimage*			pad_copy(Vector3<long> sz, int fill_type=0, double fill=0);
	int 			shrink_wrap(Vector3<long> nusize, Vector3<long> translate);
	int 			enlarge(Vector3<long> scale);
	// Montaging methods
private:
	int				montage_one(Bimage* pm, long zz, long mi, long cols, long rows, int flipy);
public:
	Bimage*			montage(int first, int cols, int rows, int skip=0, int flipy=0);
	// Editing methods
	int				shape(int type, Vector3<long> rect, Vector3<double> start,
						double width, int fill_type=0, double fill=0, bool wrap=0);
	int				shape(long nn, int type, Vector3<long> rect, Vector3<double> start,
						double width, int fill_type=0, double fill=0, bool wrap=0);
	int				line(Vector3<double> start, Vector3<double> end, double width, int fill_type=0, double fill=0);
	Bimage*			edge_mask(int type, Vector3<long> rect,
						Vector3<double> start, double width);
	int				edge(int type, Vector3<long> rect, Vector3<double> start,
						double width, int fill_type=0, double fill=0);
	int				edge(long nn, int type, Vector3<long> rect, Vector3<double> start,
						double width, int fill_type=0, double fill=0);
	Bimage*			extract_edge_difference();
	int				hanning_taper(double fill=0);
	int				sphere(Vector3<double> center, double radius,
						double width=0, int fill_type=0, double fill=0, bool wrap=0);
	int				cylinder(Vector3<double> center, double radius,
						double height, double width, int fill_type=0, double fill=0, bool wrap=0);
	int				gaussian_sphere(long nn, Vector3<double> center, double sigma, double amp, bool wrap=0);
	int				shell(Vector3<double> center, double minrad,
						double maxrad, double width, int fill_type=0, double fill=0);
	int				shell(long nn, Vector3<double> center, double minrad,
						double maxrad, double width, int fill_type=0, double fill=0);
	int				shell_wrap(Vector3<double> center, double minrad,
						double maxrad, double width, int fill_type, double fill);
	int				shell_wrap(long nn, Vector3<double> center, double minrad,
						double maxrad, double width, int fill_type=0, double fill=0);
	int				bar(Vector3<double> start, Vector3<double> end,
						double width, double edge_width, int fill_type=0, double fill=0);
	int				bar(long nn, Vector3<double> start, Vector3<double> end,
						double width, double edge_width, int fill_type=0, double fill=0);
	int				quadric(double* param);
	int				chirp(double freq_scale, double freq_shift=0);
	int				fill_gaps(long step);
	int				interpolate_gaps(long step);
	// Transformation methods
	int				transform(Vector3<double> scale,
						Vector3<double> origin, Vector3<double> translate,
						Matrix3 mat, int fill_type=0, double fill=0);
	Bimage*			transform(Vector3<long> nusize, Vector3<double> scale,
						Vector3<double> origin, Vector3<double> translate,
						Matrix3 mat, int fill_type=0, double fill=0);
	void			transform_voxel(long i, Bimage* pt, long nn, Vector3<double> oldorigin, 
						Vector3<double> nuorigin, Matrix3 affmat, double fill);
	Bimage*			transform(long nn, Vector3<long> nusize, Vector3<double> scale,
						Vector3<double> origin, Vector3<double> translate,
						Matrix3 mat, int fill_type=0, double fill=0);
	int				rotate();
	int				rotate(double angle);
	int				rotate(Vector3<double> axis, double angle);
	int				rotate(View view);
	int				rotate(Vector3<double> translate, View view);
	int				rotate(Matrix3 mat);
	int				rotate(Vector3<double> translate, Matrix3 mat);
	Bimage*			rotate(Vector3<long> nusize);
	Bimage*			rotate(Vector3<long> nusize, double angle);
	Bimage*			rotate(Vector3<long> nusize, Vector3<double> axis, double angle);
	Bimage*			rotate(Vector3<long> nusize, View view);
	Bimage*			rotate(Vector3<long> nusize, Vector3<double> translate, View view);
	Bimage*			rotate(Vector3<long> nusize, Matrix3 mat);
	int				rotate_and_add(Bimage* p, Vector3<double> origin, View view);
	Bimage*			orient(View* views);
	int				mirror();
	int 			shift(Vector3<double> vec, int fill_type=0, double fill=0);
	int 			shift(long nn, Vector3<double> vec, int fill_type=0, double fill=0);
	int 			shift_wrap(Vector3<double> vec);
	int 			shift_wrap(long nn, Vector3<double> vec);
	int 			center(int fill_type=0, double fill=0);
	int 			center_wrap();
	int 			zero_origin();
	Vector3<long>	reciprocal_half() {
		return Vector3<long>((x+1)/2, (y+1)/2, (z+1)/2);
	}
	Bimage*			scale_to_same_size(Bimage* pref);
	Bimage*			scale_to_reference(Bimage* pref, Bimage* pmask, double scalemin, 
						double scalemax, double step);
	Bimage*			scale_to_reference(Bimage* pref, Bimage* pmask=NULL);
	// Symmetry methods
	double 			symmetrize(Bstring& symmetry_string, int flag) {
		View	ref_view;
		return symmetrize(symmetry_string, ref_view, flag);
	}
	double 			symmetrize(Bsymmetry sym, int flag) {
		View	ref_view;
		return symmetrize(sym, ref_view, flag);
	}
	double 			symmetrize(Bstring& symmetry_string, View ref_view, int flag) {
		Bsymmetry 	sym(symmetry_string);
		return symmetrize(sym, ref_view, flag);
	}
	double 			symmetrize_cyclic(int cyclic, int flag) {
		Bstring		symmetry_string(cyclic, "C%d");
		return symmetrize(symmetry_string, flag);
	}
	double 			symmetrize(Bsymmetry sym, View ref_view, int flag);
	double 			check_point_group(Bstring& check_string);
	double 			find_cyclic_point_group(Bsymmetry& sym,
						int binfac, double hires, double lores);
	double 			find_point_group(Bsymmetry& sym, double angle_step,
						int binfac, double hires, double lores, int flags);
	long			rotate_to_axis(Bsymmetry& sym, long axis, long axis_flag);
	int				change_symmetry(Bsymmetry& symold, Bsymmetry& symnu,
						double radius, double z_slope);
	Matrix3			symmetry_equivalent(Bimage* ptemp, Bimage* pmask, Bsymmetry& sym);
	Matrix3			symmetry_equivalent_cyclic(Bimage* pref, Bimage* pmask, Bsymmetry& sym);
	Bimage*			levelmask_asymmetric_units(Bsymmetry& sym, int index);
	int				replicate_asymmetric_unit(Bsymmetry& sym);
	Bimage*			find_symmetric_view(Bimage* ptemp, Bsymmetry& sym,
						double phi_step, double theta_step, double alpha_step, Vector3<double> shift);
	// Helical methods
	double			helix_interpolate(long i, double helix_rise, double helix_angle,
						int zmin, int zmax, double radius, int norm_flag=1);
	double			dyad_interpolate(long i, int norm_flag=1);
	double			tube_interpolate(long i, int h, int k,
						double latconst, int zmin, int zmax, double radius, int norm_flag=1);
	double 			helix_symmetrize(double helix_rise, double helix_angle,
						int dyad_axis, int zmin, int zmax, double radius, int norm_flag=1);
	Bplot* 			seamed_helix_symmetrize(double helix_rise, double helix_angle, 
						double seam_shift, int dyad_axis, int zmin, int zmax, 
						double radius, int norm_flag);
	Bimage*			symmetrize_cylinder(int flag);
	int				symmetrize_cylinder();
	double 			tube_symmetrize(int h, int k, double latconst,
						int zmin, int zmax, double radius, int norm_flag=1);
	double 			convert_to_helix(double helix_rise, double helix_angle, 
						Vector3<double> offset);
	int				distort_elliptically(double angle, double shift);
	Vector3<double>	test_helix_parameters(double angle, double hires, double lores,
						Vector3<long> mask_size, Vector3<long> mask_start, long max_iter,
						fft_plan planf, fft_plan planb, double& cc);
	Bplot*			find_helix_parameters(double angle_start, double angle_end,
						double angle_step, int bin, double hires, double lores, double radius);
	double			test_helix_parameters(double rise, double angle,
						Vector3<long> mask_size, Vector3<long> mask_start);
	Bplot*			find_helix_parameters(double rise_start, double rise_end,
						double rise_step, double angle_start, double angle_end, double angle_step,
						int bin, double radius);
	int				helix_segment_correlation_one(long i,
						double angle_start, double angle_end, double angle_step,
						int bin, double hires, double lores, double radius,
						fft_plan planf, fft_plan planb, double* cc);
	Bplot*			helix_segment_correlation(int thickness,
						double angle_start, double angle_end, double angle_step,
						int bin, double hires, double lores, double radius);
	int				transform_lines();
	Bimage*			helical_cross_section(double helix_rise, double helix_angle,
						double scale, double hires);
	double			snvariance(double snradius);
	int				extrude_cross_section(long length, double helix_rise,
						double helix_angle, int fill_type, double fill);
	Bplot*			filament_width(long width, long lim_lo, long lim_hi);
	Bimage*			filament_density(double width);
	Bimage*			filament_from_projections(double hi_res, int flag=0);
	// Polar methods
	Bimage* 		radial(long minrad, long maxrad, double rad_step=1, int wrap=0);
	Bimage* 		radial(long minrad, long maxrad, double rad_step, Bimage* pmask, int wrap=0);
	Bimage* 		radial(long minrad, long maxrad, double rad_step,
						double ellipticity, double angle, Bimage* pmask, int wrap=0);
	Bimage*			radial_symmetry_adjusted(double rad_start,
						double rad_end, double rad_step,
						double spherical_fraction, Bsymmetry& sym);
	double*			radial_fit(Bimage* pref);
	Bimage* 		radial_to_full(Vector3<long> nusize, Vector3<double> origin);
	Bimage*			cartesian_to_spherical(long nannuli, long nphi, long ntheta);
	Bimage*			cartesian_to_cylindrical(long nannuli, long nphi, int flag=0);
	Bimage* 		polar_transform(long nangles, long ann_min, long ann_max,
						long dann, long zmin, long zmax, long zinc);
	Bimage*			polar_power_spectrum(double resolution, long num_angle);
	int				line_powerspectra(fft_plan plan);
	int		 		radial_shells();
	int		 		cylindrical_shells();
private:
	int				set_radial_distances(double spherical_fraction, Bsymmetry& sym);
	int				radial_section(Bimage* prad, long nn,
						long zr, double rad_start, double rad_step, double fill=0);
public:
	Bimage*			radial_sections(double rad_start, double rad_end,
						double rad_step, double spherical_fraction,
						Bsymmetry& sym, int fill_type=FILL_USER, double fill=0);
	Bimage*			radial_coverage(double threshold, double rad_step=1);
	// Topological methods
	Bimage*			topograph_to_surface(Bimage* psd, long nz, double density, double resolution);
	Bimage*			surface_to_topograph(double threshold, int dir=0);
	Bimage*			rotate_height(Matrix3 mat, Vector3<double> translate, double threshold=0);
	Bimage* 		height(View* views, double threshold=0);
	// Binning methods
	int 			integer_interpolation(int integer_factor);
	int 			integer_interpolation(int integer_factor, int odd);
private:
	int				bin(long i, Vector3<long> bk, Bimage* pb);
public:
	int 			bin(long b) {
		Vector3<long>		v = {b,b,b};
		return bin(v);
	}
	int 			bin(Vector3<long> bk);
	Bimage*			bin_copy(long b) {
		Vector3<long>		v = {b,b,b};
		return bin_copy(v);
	}
	Bimage*			bin_copy(Vector3<long> bk);
	Bimage*			bin_around_origin(int bin);
	int 			median_bin(int binning);
	// Noise generation methods
	int 			noise_uniform(double rmin, double rmax);
	int				noise_gaussian(double ravg=0, double rstd=1);
	int				noise_poisson(double ravg);
	int				noise_logistical(double ravg, double rstd);
	int				noise_spectral(double alpha);
	int 			noise_uniform_distance(long number);
	// Binary masks
	long			mask(Bimage* pmask, double fill);
	long			to_mask() { return to_mask((max + min)/2); }
	long			to_mask(double threshold);
	Bimage*			mask_by_threshold(double threshold);
	Bimage*			mask_by_thresholds(vector<double> threshold);
	Bimage*			mask_by_conditional_thresholds(vector<double> threshold);
	long			mask_stats();
	long			mask_invert();
	long			mask_combine(Bimage* p, int operation);
	Bimage*			mask_extract(Bimage* pmask);
	int				max_in_kernel(long ksize);
	long			mask_dilate_erode(unsigned char dir);
	long			mask_dilate(long times=1);
	long			mask_erode(long times=1);
	long 			mask_open(int times=1);
	long 			mask_close(int times=1);
	long 			mask_fill(Vector3<long> voxel);
	long			mask_shell(Vector3<double> origin, double rad_min, double rad_max) {
		return shell(origin, rad_min, rad_max, 0.1, FILL_USER, 1.99);
	}
	long			mask_shell_wrap(Vector3<double> origin, double rad_min, double rad_max) {
		return shell_wrap(origin, rad_min, rad_max, 0.1, FILL_USER, 1.99);
	}
	long			mask_plane(Vector3<double> origin, Vector3<double> normal);
	long			mask_rectangle(double length, double width,
						double rect_angle, int wrap);
	long			mask_symmetrize(Bsymmetry& sym);
	Vector3<double>	mask_fspace_resize(Vector3<long> nusize);
	vector<double>	fspace_default_bands(double res_lo, double res_hi);
	int				mask_fspace_banded(vector<double>& band);
	long			mask_missing_wedge(Vector3<double> origin, double tilt_axis,
						double tilt_neg, double tilt_pos, double resolution);
	long			mask_missing_pyramid(Vector3<double> origin, double tilt_axis1,
						double tilt_axis2, double tilt_neg1, double tilt_pos1,
						double tilt_neg2, double tilt_pos2, double resolution);
	long			mask_missing_cone(Vector3<double> origin,
						double mis_ang, double resolution);
	long			mask_missing_find(Vector3<double> ori, double resolution, Bstring& mis_type);
	long			mask_pack_plane(Matrix3 mat, double hi_res, double scale);
	// Ternary masks
	double	 		variance_threshold(double lowvar);
	Bimage* 		variance_mask(long kernel_size, double lowvar=1e-6, int bkg_flag=0);
	// Level masks
	Bimage*			tile_mask(long step);
	long			mask_split(long voxels_per_level);
	long			levelmask_collapse();
	long			levelmask_add(Bimage* pmask, int add_level=1);
	long			levelmask_dilate(int times);
	long			levelmask_dilate();
	long			levelmask_select(Bstring& select_list, int flag=0);
	long			levelmask_combine(Bstring& select_list);
	long			levelmask_select(long nn, Vector3<long> voxel);
	long			levelmask_select(Bimage* pmask);
	long			levelmask_switch(long index1, long index2);
	int		 		level_masked_stats(Bimage* pmask);
	long			levelmask_symmetrize(Bsymmetry& sym);
	Bimage* 		level_mask_extract(Bimage* pmask, int fill_type=0, double fill=0);
	double	 		levelmask_average_region_size();
	long			levelmask_clean();
	long 			levelmask_colorize();
	Bimage*			levelmask_color_by_size();
	int				levelmask_region_size();
	Bplot*			levelmask_size_histogram();
	Matrix			mask_interface_matrix(int img_num);
	long			mask_region_interfaces(int reg_num);
	long			mask_merge_delete(long min_size, long min_if);
	// Soft masks
	long			apply_soft_mask(long nn, Bimage* pmask, int fill_type, double fill);
	// Segmentation
	Bimage* 		regions(double threshold, int sign);
	long			region_assign(Bimage* pmask, long idx,
						long region_number, double threshold, int sign);
	int				region_threshold_series(double threshold_first,
						double threshold_last, double threshold_step);
	long 			region_flood(Bimage* pmask, double threshold_hi,
						double threshold_lo, double threshold_step, int fill_borders);
	long			check_neighbors(long idx);
	Bimage*			track_gradient(double threshold, int flag=0);
	Bimage*			region_peaks(long kernel_size, double threshold, int flood=0, int wrap=0);
	long 			blobs(double threshold, double min_size, double max_size,
						double setvalue, int sign);
	int				filter_extremes();
	int				filter_extremes(int mod_flag);
	int				filter_extremes(double tmin, double tmax, int kernel=3);
	long			replace_maxima(double threshold);
	double 			mass_threshold(long img_num, double mol_weight, double rho);
	double 			mass_at_threshold(long img_num, double threshold, double rho);
	Bimage*			internal_volume(double threshold);
	Bimage*			internal_volume(double threshold, int mask_out_freq);
	Bimage*			kmeans_segment(long nregion=2, long max_iter=10, double ratio=1);
	GSgraph			graph_setup(int connect_type);
	GSgraph			graph_segment(int type=1, int connect_type=0,
						double complexity=0, long min_size=0);
	Bimage*			graph_segments_to_image(GSgraph& g);
	Bimage*			graph_segments_to_mask(GSgraph& g);
	int				superpixels_update(Bimage* pmask, vector<long> vstep,
				double colorweight, vector<Bsuperpixel>& seg);
	vector<Bsuperpixel>	superpixels_from_mask(long cc, long step);
	vector<Bsuperpixel>	superpixels(long step, double colorweight=0.2, long iterations=10, double stop=1);
	vector<Bsuperpixel>	superpixels(long step, double colorweight=0.2, long iterations=10, long bin_levels=1, double stop=1);
	int				impose_superpixels(Bimage* pmask, vector<Bsuperpixel>& seg, int impose);
} ;
#define _Bimage_
#endif

// Function prototypes
