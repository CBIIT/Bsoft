/**
@file	Bimage_fft.cpp
@brief	General FFT for n-dimensional data
@author Bernard Heymann
@date	Created: 19980805
@date	Modified: 20220805

		Implementing the FFTW3 library
**/

#include "Bimage.h"
#include "timer.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Sets up a plan for fast Fourier transforms.
@param 	dir			direction of transformation (FFTW_FORWARD or FFTW_BACKWARD)
@param 	opt			optimization (0=FFTW_ESTIMATE, 1=FFTW_MEASURE, 2=FFTW_PATIENT, 3=FFTW_EXHAUSTIVE).
@return fft_plan 	FFTW plan.

	FFTW library (www.fftw.org).
	The size and direction determines the plan.
	Both FFTW versions 2 and 3 are supported.

**/
fft_plan	Bimage::fft_setup(fft_direction dir, int opt)
{
	return fft_setup_plan(size(), dir, opt);
}

/**
@brief 	Fast Fourier transforms an image with a specified conversion.
@param 	dir			direction of transformation (FFTW_FORWARD or FFTW_BACKWARD)
@param 	norm_flag	normalization: 0=none, 1=sqrtN, 2=N.
@param 	conv		a flag for converting after transform.
@return int 			error code.

**/
int 		Bimage::fft(fft_direction dir, int norm_flag, ComplexConversion conv)
{
	if ( fft(dir, norm_flag) )
		return error_show("Bimage::fft", __FILE__, __LINE__);
/*
	switch ( conv ) {
		case NoConversion: break;
		case Real: complex_to_real(); break;
		case Imaginary: complex_to_imaginary(); break;
		case Amplitude: complex_to_amplitudes(); break;
		case Intensity: complex_to_intensities(); break;
	}
	
	if ( conv ) fourier_type(NoTransform);
*/
	complex_convert(conv);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fft: datatype=" << datatype << " compoundtype=" << compoundtype << endl;
	
	return 0;
}

/**
@brief 	Fast Fourier transforms an image.
@param 	dir			direction of transformation (FFTW_FORWARD or FFTW_BACKWARD)
@param 	norm_flag	normalization: 0=none, 1=sqrtN, 2=N.
@return int 			error code.

	FFTW library (www.fftw.org).
	A multi-image 1D, 2D and 3D data set is transformed forward or backward 
	and optionally rescaled by 1/sqrt(N) or N.  The forward transformation 
	has a negative signed exponent in the kernel and the backward transform 
	a positive sign. The transformation is done in place and the resultant 
	data are returned within the original image structure.
	For both directions the resultant image is complex.

**/
int 		Bimage::fft(fft_direction dir, int norm_flag)
{
	if ( !d.uc )
		return error_show("Error in Bimage::fft: Cannot Fourier transform - the data block is empty!", __FILE__, __LINE__);
	
	if ( c < 3 ) {
		change_type(Float);
		simple_to_complex();
	} else if ( c > 2 ) {
		multi_channel_to_complex();
	}
	
    long		   i, imgsize(x*y*z);

	if ( sizeof(fft_complex) != c*data_type_size() ) {
		error_show("Error in Bimage::fft", __FILE__, __LINE__);
		cerr << "The FFTW complex number size = " << sizeof(fft_complex) << " (should be " << c*data_type_size() << ")!" << endl;
		return -1;
	}

	if ( verbose & VERB_FULL ) {
		if ( dir == FFTW_FORWARD ) cout << "Doing a forward FFT:" << endl;
		else cout << "Doing a backward FFT:" << endl;
#ifdef FFTW_VERSION
		cout << "Using FFTW:                     Version " << FFTW_VERSION << endl;
#endif
    	cout << "Complex size:                   " << sizeof(fft_complex) << endl;
    	cout << "Image size:                     " << size() << endl;
    	cout << "Number of images:               " << n << endl;
   		cout << "Normalization:                  " << norm_flag << endl;
		cout << endl;
	}

	fft_plan		plan = fft_setup(dir, 0);

	if ( n == 1 ) {
		fftw(plan, (Complex<float> *) d.uc);
	} else {
#ifdef HAVE_GCD
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::fft: GCD FFTW3" << endl;
		dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
			Complex<float>* 	data = (Complex<float> *) (data_pointer(2*nn*imgsize));
			fftw(plan, data);
		});
#else
#ifdef HAVE_OMP
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::fft: OpenMP FFTW3" << endl;
#pragma omp parallel for
#endif
		for ( long nn=0; nn<n; nn++ ) {
			Complex<float>* 	data = (Complex<float> *) (data_pointer(2*nn*imgsize));
			fftw(plan, data);
		}
#endif
	}

	fft_destroy_plan(plan);
	
	// Scale data
	double		scale = 1.0/imgsize;
	if ( norm_flag == 1 ) scale = sqrt(scale);
	if ( norm_flag ) {
		for ( i=0, imgsize *= n; i<imgsize; i++ ) set(i, complex(i) * scale);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::fft: scale = " << scale << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fft: FFT done! (" << dir << ")" << endl;
	
	data_type(Float);
	compound_type(TComplex);
	channels(2);
	fourier_type(Standard);
	
	return 0;
}

/**
@brief 	Fast Fourier transforms an image.
@param 	plan		Fourier transform plan.
@param 	norm_flag	normalization: 0=none, 1=sqrtN, 2=N.
@return int 			error code.

	FFTW library (www.fftw.org).
	A multi-image 1D, 2D and 3D data set is transformed forward or backward 
	and optionally rescaled by 1/sqrt(N) or N.  The forward transformation 
	has a negative signed exponent in the kernel and the backward transform 
	a positive sign. The transformation is done in place and the resultant 
	data are returned within the original image structure.
	For both directions the resultant image is complex.
	Requirement: The plan must be derived from the same size image.

**/
int 		Bimage::fft(fft_plan plan, int norm_flag)
{
	if ( !d.uc )
		return error_show("Error in Bimage::fft: Cannot Fourier transform - the data block is empty!", __FILE__, __LINE__);

	if ( c < 3 ) {
		simple_to_complex();
		change_type(Float);
	} else if ( c > 2 ) {
		multi_channel_to_complex();
	}

    long		   	i, nn, imgsize(x*y*z);

	if ( sizeof(fft_complex) != c*data_type_size() ) {
		error_show("Error in Bimage::fft_complex", __FILE__, __LINE__);
		cerr << "The FFTW complex number size = " << sizeof(fft_complex) << " (should be " << c*data_type_size() << ")!" << endl;
		return -1;
	}

	if ( verbose & VERB_FULL ) {
//		if ( dir == FFTW_FORWARD ) cout << "Doing a forward FFT:" << endl;
//		else cout << "Doing a backward FFT:" << endl;
#ifdef FFTW_VERSION
		cout << "Using FFTW:                     Version " << FFTW_VERSION << endl;
#endif
    	cout << "Complex size:                   " << sizeof(fft_complex) << endl;
    	cout << "Image size:                     " << size() << endl;
    	cout << "Number of images:               " << n << endl;
   		cout << "Normalization:                  " << norm_flag << endl;
		cout << endl;
	}

    // Specify FFTW arrays
	Complex<float>*		data = (Complex<float> *) data_pointer();

	for ( nn=0; nn<n; nn++, data += imgsize ) {
//		cout << nn << tab << data << endl;
		fftw(plan, data);
	}

	// Scale data
	double		scale = 1.0L/imgsize;
	if ( norm_flag == 1 ) scale = sqrt(scale);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fft: scale = " << scale << endl;
	if ( norm_flag )
		for ( i=0, imgsize *= n; i<imgsize; i++ ) set(i, complex(i) * scale);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fft: FFT done!" << endl << endl;
	
	data_type(Float);
	compound_type(TComplex);
	channels(2);
	fourier_type(Standard);
	
	return 0;
}

/**
@brief 	Tiled Fourier transform.
@param 	dir			direction of transformation (FFTW_FORWARD or FFTW_BACKWARD)
@param	tile_size	tile size.
@param 	norm_flag	normalization: 0=none, 1=sqrtN, 2=N.
@return int 			error code.

	The 2d/3d image is tiled, each tile transformed and phase shifted.
	The tiles are then averaged and returned as a new image.

**/
int			Bimage::fft(fft_direction dir, Vector3<long> tile_size, int norm_flag)
{
	Vector3<long> 	start, ext_size(size()), step_size;

	if ( verbose & VERB_PROCESS )
		cout << "Tile size:                      " << tile_size << endl << endl;
	
	Bimage* 		pex = extract_tiles(0, start, ext_size, tile_size, step_size, 0);
	
	if ( verbose & VERB_FULL )
		cout << "FFT" << endl;
	
	if ( pex->fft(dir, norm_flag) )
		return -1;
	
	if ( verbose & VERB_FULL ) {
		cout << "Compound type: " << pex->compound_type() << endl;
		cout << "Shift phases" << endl;
	}
	
	for ( long nn=0; nn<pex->images(); nn++ )
		pex->phase_shift(nn, -pex->image[nn].origin());
	
	if ( verbose & VERB_FULL ) {
		cout << "#\tOrigin" << endl;
		for ( long nn=0; nn<pex->images(); nn++ )
			cout << nn << tab << pex->image[nn].origin() << endl;
	}
	
	pex->statistics();
	pex->information();

	cout << "Summing images" << endl;
	pex->sum_images();
	pex->information();

	size(pex->size());
	page_size(size());
	compound_type(TComplex);
	fourier_type(Standard);
	origin(size()/2);

    data_assign(pex->data_pointer());
//	pex->data_pointer(NULL);
	pex->d.uc = NULL;
//	cout << pex->data_pointer() << endl;
	delete pex;

	return 0;
}

/**
@brief 	Fourier transform only of the z columns.
@param 	dir			direction of transformation (FFTW_FORWARD or FFTW_BACKWARD)
@param 	norm_flag	normalization: 0=none, 1=sqrtZ, 2=Z.
@return int 			error code.

	Only the z columns are transformed.

**/
int			Bimage::fftz(fft_direction dir, int norm_flag)
{
	simple_to_complex();
	
//	information();
	
	if ( verbose )
		cout << "Fourier transforming z columns" << endl << endl;
	
	Complex<float>*	data = (Complex<float> *) data_pointer();

	int 			rank(1); 			// 1d transforms
	int 			zcol[] = {int(z)}; 	// 1d transforms of length z
	int 			ntr(x*y);			// number of transforms
	int 			idist(1), odist(1);
	int 			istride(x*y), ostride(x*y); // distance between two elements in the same column
	int 			*inembed = zcol, *onembed = zcol;
	
	fft_complex*	in = (fft_complex *) data;
	fft_complex*	out = in;
	fft_plan		plan = fftwf_plan_many_dft(rank, zcol, ntr,
						in, inembed, istride, idist,
						out, onembed, ostride, odist, dir, FFTW_ESTIMATE);

	long			i, nn, imgsize(x*y*z);
	
	for ( nn=0; nn<n; nn++, data += imgsize ) {
//		cout << nn << tab << data << endl;
		fftw(plan, data);
	}

	fft_destroy_plan(plan);

	// Scale data
	double		scale = 1.0L/z;
	if ( norm_flag == 1 ) scale = sqrt(scale);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fftz: scale = " << scale << endl;
	if ( norm_flag )
		for ( i=0, imgsize *= n; i<imgsize; i++ ) set(i, complex(i) * scale);

	return 0;
}

/**
@brief 	Resizes a "standard" transform.
@param 	nusize			new image size.
@return Vector3<double>	scale.

	A standard transform is resized by inserting or removing
	rows or columns in the middle of the data set.

**/
Vector3<double>	Bimage::change_transform_size(Vector3<long> nusize)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::change_transform_size: size=" << nusize << endl;

	if ( z < 2 ) nusize[2] = 1;
	Vector3<double>	scale(nusize[0]*1.0L/x, nusize[1]*1.0L/y, nusize[2]*1.0L/z);

	if ( nusize.volume() < 1 ) return scale;
	if ( nusize == size() ) return scale;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::change_transform_size: scale=" << scale << endl;

	if ( compoundtype != TComplex ) {
		cerr << "Error: Bimage::change_transform_size: Only complex numbers supported!" << endl;
		return scale;
    }
	
    long	   		i, j, nn, xx, yy, zz, iy, iz;
    long     	    xold, yold, zold, zerox, zeroy, zeroz;
	Vector3<long>	hold, h, d;
    long	   		ds(n*nusize.volume());
	Complex<float>	cv;
    Complex<float>*	nudata = new Complex<float>[ds];
	float*			fom = NULL;
	float*			nufom = NULL;
	if ( next ) {
		fom = (float *) next->data_pointer();
		nufom = new float[ds];
		for ( i=0; i<ds; i++ ) {
			nudata[i] = cv;
			nufom[i] = 0;
		}
	}
	
	d = nusize - size();
	hold = (size() - 1)/2;
	h = (nusize - 1)/2;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::change_transform_size: datasize=" << ds << endl;
	
	if ( verbose & VERB_FULL )
	    cout << "Changing transform to size:     " << nusize << " voxels" << endl << endl;

    for ( j=nn=0; nn<n; nn++ ) {
	    for ( zz=0; zz<nusize[2]; zz++ ) {
    		zold = zz;
			if ( zz > h[2] ) zold -= d[2];
			zeroz = 0;
			if ( zz > hold[2] && zz <= d[2] + hold[2] ) zeroz += 1;
			iz = (nn*z + zold)*y;
    		for ( yy=0; yy<nusize[1]; yy++ ) {
			 	yold = yy;
				if ( yy > h[1] ) yold -= d[1];
				zeroy = zeroz;
				if ( yy > hold[1] && yy <= d[1] + hold[1] ) zeroy += 1;
				iy = (iz + yold)*x;
				for ( xx=0; xx<nusize[0]; xx++, j++ ) {
					xold = xx;
					if ( xx > h[0] ) xold -= d[0];
					zerox = zeroy;
					if ( xx > hold[0] && xx <= d[0] + hold[0] ) zerox += 1;
					if ( zerox < 1 ) {
						i = iy + xold;
						nudata[j] = complex(i);
						if ( fom ) nufom[j] = fom[i];
					}
				}
	    	}
		}
		image[nn].origin(scale*image[nn].origin());
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::change_transform_size: origin[" << n << "]=" << image[nn].origin() << endl;
	
    }
    
	sampling(image->sampling()/scale);
	size(nusize[0], nusize[1], nusize[2]);
	page_size(size());

    data_assign((unsigned char *) nudata);
	if ( next ) {
		next->size(nusize);
    	next->data_assign((unsigned char *) nufom);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::change_transform_size: sampling=" << sampling(0) << endl;
	
	return scale;
}

