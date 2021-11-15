/**
@file	mg_ctf_sim.cpp
@brief	Functions to simulate tilted micrographs
@author Bernard Heymann
@date	Created: 20150224
@date	Modified: 20200525
**/

#include "mg_ctf.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen



int			img_ttf_simulate_line(Bimage* p, Bimage* pn, long img_num, long i, 
				char ax, CTFparam& ctf, int action, double wiener,
				double scale, fft_plan plan)
{
	double		ddef = ((double)i - p->image[img_num].origin()[0])*scale;
	if ( ax == 'x' )
		ddef = ((double)i - p->image[img_num].origin()[1])*scale;
	
	CTFparam 	ctf_copy = ctf;
	ctf_copy.defocus_average(ctf.defocus_average()+ddef);

	Bimage*		pctf = img_ctf_calculate(ctf_copy, action, wiener, pn->size(), pn->sampling(0), 1e6, 0);
	
	Bimage*		pnc = pn->extract(img_num);
	
	pnc->complex_multiply(pctf);
	
	delete pctf;

	pnc->fft_back(plan);
	
	long		xx, yy, j;
	if ( ax == 'x' )
		for ( xx=0, i*=p->sizeX(), j=i+img_num*p->image_size(); xx<p->sizeX(); xx++, i++, j++ )
			p->set(j, (*pnc)[i]);
	else
		for ( yy=0, j=i+img_num*p->image_size(); yy<p->sizeY(); yy++, i+=p->sizeX(), j+=p->sizeX() )
			p->set(j, (*pnc)[i]);
	
	delete pnc;
	
	return 0;
}

long		img_ttf_simulate_one(Bimage* p, Bimage* pn, long img_num,
				CTFparam& ctf, int action, double wiener, double tilt,
				double axis, fft_plan plan)

{
	char			ax = (fabs(axis)<M_PI/4)? 'x': 'y';
	long			nn = (ax=='x')? pn->sizeY(): pn->sizeX();
	double			scale(pn->sampling(0)[1]*tan(tilt));

	if ( verbose )
		cout << "Image " << img_num << ": " << tilt*180.0/M_PI << endl;
	
#ifdef HAVE_GCD
	__block	long		n(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(nn, dispatch_get_global_queue(0, 0), ^(size_t i){
		img_ttf_simulate_line(p, pn, img_num, i, ax, ctf, action, wiener, scale, plan);
		dispatch_sync(myq, ^{
			n++;
			if ( verbose )
				cout << "Line " << n << "\r" << flush;
		});
	});
#else
	long			n(0);
#pragma omp parallel for
	for ( long i=0; i<nn; i++ ) {
		img_ttf_simulate_line(p, pn, img_num, i, ax, ctf, action, wiener, scale, plan);
	#pragma omp critical
		{
			n++;
			if ( verbose )
				cout << "Line " << n << "\r" << flush;
		}
	}
#endif
	
	if ( verbose )
		cout << "Lines calculated [" << img_num << "]:            " << n << endl;

	return n;
}

/**
@brief 	Calculates an exact tilted CTF imposed image.
@param 	*pn			image structure.
@param 	&ctf		CTF parameter structure.
@param 	action		0=apply, 1=flip.
@param 	wiener		Wiener factor.
@param 	tilt		tilt angle in radians.
@param	tilt_inc	tilt increment per sub-image in radians.
@param	axis		tilt axis direction in radians.
@return Bimage*		new image.

	A CTF function is calculated for each line parallel to the tilt axis,
	applied to the whole image, and the transformed line written into a new image.

**/
Bimage*		img_ttf_simulate(Bimage* pn, CTFparam& ctf, int action,
				double wiener, double tilt, double tilt_inc, double axis)

{
	if ( pn->sampling(0).volume() < 0.1 ) pn->sampling(Vector3<double>(1,1,1));

	pn->fft();
	
	Bimage*			p = new Bimage(Float, TSimple, pn->size(), pn->images());
	p->sampling(pn->sampling(0));
	p->origin(p->size()/2);

	if ( verbose ) {
		ctf.show();
		cout << "Action:                         " << action << endl;
		cout << "Size:                           " << p->size() << endl;
		cout << "Images:                         " << p->images() << endl;
		cout << "Sampling:                       " << p->sampling(0) << endl;
		cout << "Tilt angle and increment:       " << tilt*180/M_PI << " " << tilt_inc*180/M_PI << endl;
		cout << "Tilt axis:                      " << axis*180/M_PI << endl << endl;
	}
	
	fft_plan		plan = p->fft_setup(FFTW_BACKWARD, 0);

	for ( long nn=0; nn<pn->images(); ++nn, tilt += tilt_inc )
		img_ttf_simulate_one(p, pn, nn, ctf, action, wiener, tilt, axis, plan);
		
	fft_destroy_plan(plan);

	return p;
}

