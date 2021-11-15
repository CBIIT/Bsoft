/**
@file	img_combine.cpp
@brief	Functions to combine two images in various ways
@author Bernard Heymann
@date	Created: 19990219
@date	Modified: 20190208
**/

#include "Bimage.h"
#include "img_combine.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Adds multiple images together with optional variance or standard deviation.
@param 	*file_list 	list of file names.
@param 	flags		flags to modify summation.
@return Bimage* 	resultant image (floating point).

	Images are read from a number files and added.
	All the images must be the same size, but could have different numbers of sub-images.
	The flags that can set are:
		1	calculate the average in stead of the sum.
		2	calculate the variance as FOM
		4	calculate the standard deviation as FOM (supercedes the variance)
	All images are converted to floating point.

**/
Bimage* 	img_add(Bstring* file_list, int flags)
{
	if ( !file_list ) return NULL;
	
	int				calcavg = flags & 1;
	int 			calcfom = (flags & 2) >> 1;
	long			nfiles = count_list((char *) file_list);
	
	long			i, j, n(0), nf, nimg(0);
	double			v, va, w;
	Bstring*		filename;
	Bimage*			p = NULL;

	Bimage*	 		psum = read_img(*file_list, 1, 0);
	psum->clear();
	psum->change_type(Float);
	
	Bimage*			pfom = psum->next = psum->copy();
	
	long			imgsize = psum->sizeX()*psum->sizeY()*psum->sizeZ()*psum->channels();

	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) )
		cout << endl << "Adding " << nfiles << " files together:" << endl;	
	
	for ( nf=0, filename = file_list; filename && nf<nfiles; filename = filename->next, nf++ ) {
		p = read_img(*filename, 1, -1);
		if ( p != NULL ) {
			if ( verbose & VERB_LABEL )
				cout << "Adding file " << nf << " with images " << p->images() << endl;
			for ( n=j=0; n<p->images(); n++ ) {
				for ( i=n*imgsize; i<(n+1)*imgsize; i++, j++ ) {
					v = (*p)[j];
					psum->add(j, v);
					if ( calcfom ) pfom->add(j, v*v);
				}
			}
			nimg += p->images();
			delete p;
		}
	}

	w = 1.0/nimg;
	
	for ( i=n=0; n<psum->images(); n++ ) {
		for ( j=0; j<imgsize; j++, i++ ) {
			va = (*psum)[i] * w;
			if ( calcavg ) psum->set(i, va);
			else psum->set(i, (*psum)[i]);
			v = (*pfom)[i];
			if ( calcfom > 0 ) {
				v = v * w - va*va;
				if ( v < 0 ) v = 0;
			}
			if ( calcfom > 1 ) v = sqrt(v);
			pfom->set(i, v);
		}
	}
	
	psum->statistics();
	pfom->statistics();
	
	return psum;
}


/**
@brief 	Sets up a list of images for concatenation or summation.
@param 	*file_list		list of file names.
@param 	&nimg			number of concatenated images.
@param 	cat				flag to indicate concatenation.
@return Bimage*			new image into which to write data.

	The images can have different numbers of sub-images, sizes and data types.

**/
Bimage*	 	img_setup_combined(Bstring* file_list, long& nimg, int cat)
{
	if ( !file_list ) return NULL;
	
	long			nf(0), nn(0);
	Vector3<double>	sam(1,1,1);
	Bstring*		filename;
	Bimage*			p;
	Bimage*			pc = new Bimage;
	
	nimg = 0;
	
	for ( filename = file_list; filename; filename = filename->next, nf++ ) {
		p = read_img(*filename, 0, -1);
		if ( p != NULL ) {
			if ( pc->channels() < p->channels() ) pc->channels(p->channels());
			if ( pc->data_type() < p->data_type() ) pc->data_type(p->data_type());
			if ( pc->compound_type() < p->compound_type() ) pc->compound_type(p->compound_type());
			pc->size(pc->size().max(p->size()));
			nimg += p->images();
			if ( nn < p->images() ) nn = p->images();
			if ( nf == 0 ) sam = p->image->sampling();
			delete p;
		}
	}

	if ( cat ) pc->images(nimg);
	else pc->images(nn);
	pc->origin(pc->default_origin());
	pc->sampling(sam);

	pc->data_alloc_and_clear();
	
	pc->information();
	
	return pc;
}



/**
@brief 	Catenates a list of images into a multi-image structure.
@param 	*file_list		list of file names.
@param 	&rawstring		format for re-interpretation of file.
@param 	nudatatype		new data type (default from first image).
@param 	nusize			new size (default from images).
@param 	setZslices		flag to create 2D images from slices.
@param 	fill_type		fill type for expanding images.
@param 	fill			fill value for expanding images.
@param 	newavg			new average to set each individual image.
@param 	newstd			new standard deviation to set each individual image.
@return Bimage*			catenated image.

	The images can have different numbers of sub-images, sizes and data types.

**/
Bimage*		img_catenate(Bstring* file_list, Bstring& rawstring, DataType nudatatype, 
				Vector3<long> nusize, int setZslices, int fill_type, double fill,
				double newavg, double newstd)
{
	
	// Set up image structures
	Bimage*			p = NULL;
	Bimage*			pcat = NULL;
	Bstring			filename, *name;
	int				err(0);
	long			i, j, k;

	if ( rawstring.length() > 0 && rawstring[0] != '#' ) rawstring = "#" + rawstring;
	
	// First read the file headers to identify problems
	for ( i=j=0, name = file_list; name; name = name->next ) {
		filename = *name + rawstring;
		p = read_img(filename, 0, -1);
		if ( p ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG img_catenate: " << p->file_name() << " " << 
						p->sizeX() << " " << p->sizeY() << " " << p->sizeZ() << " " << p->images() << " " << p->channels() << endl;
			if ( i == 0 ) {
				pcat = p;
				if ( fill_type == FILL_AVERAGE ) fill = p->average();
				if ( fill_type == FILL_BACKGROUND ) fill = p->background(long(0));
				if ( p->fourier_type() ) p->fourier_type(Standard);
			} else {
				if ( setZslices && p->images() > 1 ) {
					cerr << p->file_name() << " must be a single image file to pack slices" << endl;
					err += -1;
				}
//				err += img_compatibility(pcat, p);
				if ( !p->compatible(pcat) ) err += -1;
				if ( pcat->sizeX() < p->sizeX() ) pcat->sizeX(p->sizeX());
				if ( pcat->sizeY() < p->sizeY() ) pcat->sizeY(p->sizeY());
				if ( pcat->sizeZ() < p->sizeZ() ) pcat->sizeZ(p->sizeZ());
				if ( fill_type == FILL_AVERAGE ) fill += p->average();
				if ( fill_type == FILL_BACKGROUND ) fill += p->background(long(0));
				delete p;
			}
			if ( setZslices ) i += p->sizeZ();
			else i += p->images();
			j++;
		} else err += -1;
	}
	
	if ( err < 0 || j == 0 ) {
		cerr << "Files not concatenated! Number of errors: " << -err << endl;
		bexit(err);
	}
	
	// Set up for catenation
	if ( nusize.volume() > 0 ) pcat->size(nusize);
	nusize = {pcat->sizeX(), pcat->sizeY(), pcat->sizeZ()};
	if ( fill_type != FILL_USER ) fill /= j;
	if ( setZslices ) {
		pcat->sizeZ(i);
		pcat->images(1);
	} else {
		pcat->images(i);
	}		
	if ( nudatatype > Unknown_Type ) pcat->data_type(nudatatype);
	pcat->data_alloc();
	
	if ( verbose & VERB_FULL )
		cout << "Catenated size: n=" << pcat->images() << " size=" << pcat->size()
			<< " c=" << pcat->channels() << endl << endl;
	
	// Read the images including their data
	long			imgsize = (long) pcat->size().volume()*pcat->channels();
	Vector3<long>	translate;
	
	if ( setZslices ) imgsize = (long) pcat->sizeX()*pcat->sizeY()*pcat->channels();
	
	if ( verbose )
		cout << "File\tImages\tStart" << endl;
	for ( i=0, name = file_list; name; name = name->next ) {
		filename = *name + rawstring;
		p = read_img(filename, 1, -1);
		if ( p != NULL ) {
			if ( verbose )
				cout << p->file_name() << tab << p->images()
					<< tab << i << endl;
			if ( setZslices ) nusize[2] = p->sizeZ();
			if ( newstd ) p->rescale_to_avg_std(newavg, newstd);
			translate = (nusize - p->size())/2;
			p->resize(nusize, translate, fill_type, fill);
			p->change_type(nudatatype);
			for ( j=i*imgsize, k=0; k<imgsize*p->images(); j++, k++ ) pcat->set(j, (*p)[k]);
			if ( setZslices ) i += p->sizeZ();
			else for ( j=0; j<p->images(); j++, i++ ) pcat->image[i] = p->image[j];
			delete p;
		}
	}
	
	if ( verbose )
		cout << "Number of images:               " << pcat->images() << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_catenate: final_size=" << i*imgsize << endl;
	
	return pcat;
}

/**
@brief 	Adds multiple images together with given weights.
@param 	*file_list 	list of file names.
@param 	weight		list of weights (can be NULL).
@param 	newavg		new average for rescaling.
@param 	newstd		new standard deviation for rescaling.
@param 	flags		flags to modify summation.
@return Bimage* 	resultant image (floating point).

	Images are read from a number files and added to each other, using
	the given weights to determine each contribution.
	The images are rescaled to a new average and standard deviation before 
	weighted addition. If the given standard deviation is zero or less,
	this step is omitted.
	The weighed average is calculated and returned as a new image.
	The flags that can set are:
		1	calculate the average in stead of the sum.
		2	calculate the variance as FOM
		4	calculate the standard deviation as FOM (supercedes the variance)
		8	center each image before summation
	All images are converted to floating point.

**/
Bimage* 	img_add_weighed(Bstring* file_list, vector<double> weight,
					double newavg, double newstd, int flags)
{
	if ( !file_list ) return NULL;
	
	int					calcavg = flags & 1;
	int 				calcfom = (flags & 2) >> 1;
	int					center = flags & 8;
	
	long				i, j, n(0), nf, nimg(0), c(0);
	double				v;
	Bstring*			filename;
	Bimage*				p = NULL;

	Bimage*	 			psum = img_setup_combined(file_list, nimg, 0);
	long				nfiles = count_list((char *) file_list);

	double				weightsum(0);
	if ( weight.size() < nfiles )
		for ( i=weight.size(); i<nfiles; i++ ) weight.push_back(1);

	for ( i=0; i<nfiles; i++ ) weightsum += weight[i];;

	double				degrees_of_freedom = weightsum*(1 - 1.0L/nfiles);
	
	if ( verbose & VERB_PROCESS ) {
		cout << endl << "Adding " << nimg << " images together:" << endl;
		cout << "New image size:                 " << psum->size() << endl;
		cout << "Degrees of freedom:             " << degrees_of_freedom << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << endl << "Adding " << nimg << " images together" << endl << endl;
	
	long			imgsize = psum->sizeX()*psum->sizeY()*psum->sizeZ();
	long			datasize = imgsize*psum->images();
	
	float*			fom = NULL;
	if ( calcfom ) {
		psum->next = new Bimage(Float, TSimple, psum->size(), psum->images());
		psum->next->sampling(psum->image->sampling());
		psum->next->origin(psum->image->origin());
		fom = (float *) psum->next->data_pointer();
	}

	vector<double>	imgweight(psum->images());
	
	for ( nf=0, filename = file_list; filename && nf<nfiles; filename = filename->next, nf++ ) {
		p = read_img(*filename, 1, -1);
		if ( p != NULL ) {
			if ( newstd > 0 ) p->rescale_to_avg_std(newavg, newstd);
			if ( center ) p->center_wrap();
			if ( verbose & VERB_LABEL )
				cout << "Adding image " << nf << " with weight " << weight[nf] << endl << endl;
			for ( n=j=0; n<p->images(); n++ ) {
				imgweight[n] += weight[nf];
				for ( i=n*imgsize; i<(n+1)*imgsize; i++ ) {
					for ( c=0; c<psum->channels(); c++, j++ ) {
						v = (*p)[j];
						psum->add(j, weight[nf]*v);
						if ( calcfom ) fom[i] += weight[nf]*v*v;
					}
				}
			}
			delete p;
		}
	}
	
	if ( calcfom ) {
		for ( i=n=0; n<psum->images(); n++ )
			for ( j=0; j<imgsize; j++, i++ )
				fom[i] = (fom[i] - (*psum)[i]*(*psum)[i]/imgweight[n])/degrees_of_freedom;
		if ( calcfom > 1 ) {
			for ( i=0; i<datasize; i++ ) {
				if ( fom[i] > 0 ) fom[i] = sqrt(fom[i]);
				else fom[i] = 0;
			}
		}
	}
	
	if ( calcavg )
		for ( i=n=0; n<psum->images(); n++ )
			for ( j=0; j<imgsize; j++, i++ )
				psum->set(i, (*psum)[i] / imgweight[n]);

	psum->statistics();
	
	return psum;
}

