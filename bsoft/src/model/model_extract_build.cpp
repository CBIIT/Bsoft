/**
@file	model_extract_build.cpp
@brief	Functions to extract subvolumes and build new maps.
@author Bernard Heymann
@date	Created: 20060411
@date	Modified: 20150915
**/

#include "Bimage.h"
#include "model_extract_build.h"
#include "model_util.h"
#include "file_util.h"
#include "linked_list.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/*
	The first axis is along the link.
	The second axis is the average of the component normals, usually directed
	away from the molecule group center.
	The third axes is generated as the cross product of the first two. i.e.,
	perpendicular to both.
	A new second axis is generated as the cross product between the other two,
	ensuring that all the axes are orthogonal.
*/
Matrix3		get_link_matrix(Blink* link)
{
	Matrix3			mat;
	Vector3<double>	u, v, w;
	
	u = link->comp[1]->location() - link->comp[0]->location();
	w = link->comp[0]->view().vector3();
	w += link->comp[1]->view().vector3();
	v = w.cross(u);
	w = u.cross(v);
	u.normalize();
	v.normalize();
	w.normalize();
//	for ( int i=0; i<3; i++ ) {
//		mat[i] = u[i];
//		mat[i+3] = v[i];
//		mat[i+6] = w[i];
//	}
	mat[0] = u;
	mat[1] = v;
	mat[2] = w;
	mat = mat.transpose();
//	cout << mat << endl;

	return mat;
}

Vector3<double>	component_find_shift(Bcomponent* comp, Bimage* p, Bimage* ptemp,
					Bimage* pmask, Bimage* pfsmask, View2<float> view,
					double hires, double lores, double shift_limit,
					int shift_flag, fft_plan planf, fft_plan planb)
{
	long			n, imax(0);
	Vector3<long>	size(ptemp->size());
	Vector3<double>	origin(ptemp->image->origin());
	Vector3<double>	scale(ptemp->sampling(0));
	Vector3<double>	loc = comp->location()/scale + p->image->origin();

	Matrix3			mat = view.matrix();
	Bimage*			pcomp = p->extract(0, loc, size, origin, mat);
	Bimage*			pone = NULL;
	Bimage*			pmaskrot = NULL;

	double			cc, ccmax(-1);
	Vector3<double>	shift, best_shift;
	
	if ( pmask ) {
		pcomp->rescale_to_avg_std(0, 1, pmask);
		pcomp->multiply(pmask);
	}
	
	if ( pfsmask ) pmaskrot = pfsmask->extract_wrap(0, size, mat);
	
	for ( n=0; n<ptemp->images(); n++ ) {
		pone = ptemp->extract(n);
		pone->origin(origin);
		shift = pcomp->find_shift(pone, pmaskrot, hires, lores, shift_limit, 0, 1, planf, planb, cc);
		if ( ccmax < cc ) {
			ccmax = cc;
			imax = n;
			if ( shift_flag == 1 ) shift[0] = shift[1] = 0;	// Shift only along the normal
			else if ( shift_flag == 2 ) shift[2] = 0;		// Shift only lateral
			best_shift = mat * shift;						// Shift in voxels
		}
		delete pone;
	}
	
	delete pcomp;
	delete pmaskrot;

	comp->select(imax + 1);
	comp->FOM(ccmax);
	
	return best_shift;
}

Vector3<double>	component_find_shift2(Bcomponent* comp, Bimage* p, Bimage* ptemp,
					Bimage* pmask, Bimage* pfsmask, View2<float> view,
					double hires, double lores, double shift_limit,
					int shift_flag, fft_plan planf, fft_plan planb)
{
	long			n, imax(0);
	Vector3<long>	size(ptemp->size());
	Vector3<double>	origin(ptemp->image->origin());
	Vector3<double>	scale(ptemp->sampling(0));
	Vector3<double>	loc = comp->location()/scale + p->image->origin();

	Matrix3			mat = view.matrix();
	mat = mat.transpose();

	Bimage*			pcomp = p->extract(0, loc, size, origin);
	Bimage*			pone = NULL;
	Bimage*			pmaskrot = NULL;

	double			cc, ccmax(-1);
	Vector3<double>	shift, best_shift;
	
	if ( pmask ) {
		pmaskrot = pmask->rotate(size, mat);
		pcomp->rescale_to_avg_std(0, 1, pmaskrot);
		pcomp->multiply(pmaskrot);
		delete pmaskrot;
	}
	
	for ( n=0; n<ptemp->images(); n++ ) {
		pone = ptemp->extract(n, origin, size, origin, mat);
		shift = pcomp->find_shift(pone, pfsmask, hires, lores, shift_limit, 0, 1, planf, planb, cc);
		if ( ccmax < cc ) {
			ccmax = cc;
			imax = n;
			if ( shift_flag == 1 ) shift[0] = shift[1] = 0;	// Shift only along the normal
			else if ( shift_flag == 2 ) shift[2] = 0;		// Shift only lateral
			best_shift = shift;								// Shift in voxels
		}
		delete pone;
	}
	
	delete pcomp;

	comp->select(imax + 1);
	comp->FOM(ccmax);
	
	return best_shift;
}

int		component_refine_random(Bcomponent* comp, Bimage* p, Bimage* ptemp, Bimage* pmask,
				Bimage* pfsmask, int max_iter, double hires, double lores, double max_shift,
				double max_view_angle, double max_rot_angle, int shift_flag, fft_plan planf, fft_plan planb)
{
	if ( verbose & VERB_FULL )
		cout << "Refining component " << comp->identifier() << " by random view variations" << endl;
	
	int				imax(0), iter, itm;
	double			shift_limit(max_shift/ptemp->sampling(0)[0]);
	double			cc, ccmax(-1), a, s, irm(1.0/get_rand_max());
	View2<float>	view(comp->view());
	View2<float>	comp_view(comp->view());
	Vector3<double>	scale(ptemp->sampling(0));
	Vector3<double>	r, comp_shift;
	Vector3<double>	shift;
	
	for ( iter=itm=0; iter<=max_iter; iter++ ) {
		if ( iter > 0 ) {
			r = vector3_random_gaussian(0, max_view_angle);
			a = max_rot_angle*(random()*irm*2 - 1);
			view[0] += r[0];
			view[1] += r[1];
			view[2] += r[2];
			view[3] += a;
			view.normalize();
		}
		shift = component_find_shift(comp, p, ptemp, pmask, pfsmask, view, hires, lores,
				shift_limit, shift_flag, planf, planb);
		cc = comp->FOM();
		if ( ccmax < cc ) {
			ccmax = cc;
			imax = comp->select() - 1;
			comp_shift = shift;						// Shift in voxels
			comp_view = view;
			itm = iter;
		}
	}
	
	comp_shift *= scale;		// Shift in angstroms
	s = comp_shift.length();
	if ( s > max_shift ) {
		comp_shift *= max_shift/s;
		s = comp_shift.length();
	}
	comp->shift(comp_shift);
	comp->force(comp_shift);
	comp->view(comp_view);
	comp->select(imax + 1);
	comp->FOM(ccmax);

	return itm;
}

int		component_refine(Bcomponent* comp, Bimage* p, Bimage* ptemp, Bimage* pmask,
				Bimage* pfsmask, double viewmax, double rotmax,
				double viewstep, double rotstep,
				double hires, double lores, double accuracy, double max_shift,
				int shift_flag, fft_plan planf, fft_plan planb)
{
	if ( verbose & VERB_FULL )
		cout << "Refining component " << comp->identifier() << " by contracting grid search" << endl;
	
	if ( viewmax < 0 ) viewmax = 0;
	if ( rotmax < 0 ) rotmax = 0;
	if ( viewstep < 0 ) viewstep = 0;
	if ( rotstep < 0 ) rotstep = 0;
		
	int				imax(0), itm;
	double			shift_limit(max_shift/ptemp->sampling(0)[0]);
	double			cc, ccp(-2), ccmax(-1), s;
	Vector3<double>	scale(ptemp->sampling(0));
	Vector3<double>	comp_shift;
	Vector3<double>	shift;
	
	View2<float>	comp_view(comp->view());
	list<View2<float>>	views;
//	View*			v;
	
	double			tol(1e-6);
	if ( viewstep && accuracy > viewstep ) accuracy = viewstep;
	if ( rotstep && accuracy > rotstep ) accuracy = rotstep;

	int				notdone(1);
	for ( itm=0; notdone; itm++ ) {
		views = views_within_limits2(comp_view, viewstep, viewstep, rotstep, viewmax, rotmax);
//		for ( v = views; v; v = v->next ) {
		for ( auto v = views.begin(); v != views.end(); ++v ) {
			shift = component_find_shift2(comp, p, ptemp, pmask, pfsmask, *v, hires, lores,
					shift_limit, shift_flag, planf, planb);
			cc = comp->FOM();
			if ( ccmax < cc ) {
				ccmax = cc;
				imax = comp->select() - 1;
				comp_shift = shift;						// Shift in voxels
				comp_view = *v;
			}
		}
//		kill_list((char *) views, sizeof(View));
		if ( fabs(ccp - ccmax) < tol ) {
			viewstep /= 2;			// Contract around best point
			rotstep /= 2;
		}
		ccp = ccmax;
		if ( viewstep < tol && rotstep < tol ) notdone = 0;
		if ( viewstep <= accuracy && rotstep <= accuracy ) notdone = 0;
	}

	comp_shift *= scale;		// Shift in angstroms
	s = comp_shift.length();
	if ( s > max_shift ) {
		comp_shift *= max_shift/s;
		s = comp_shift.length();
	}
	comp->shift(comp_shift);
	comp->force(comp_shift);
	comp->view(comp_view);
	comp->select(imax + 1);
	comp->FOM(ccmax);

	return itm;
}

/**
@brief 	Refines component views and positions by cross-correlation.
@param 	*model			model.
@param 	*ct_names		list of names associated with template sub-images.
@param 	*ptemp			density template.
@param 	*pmask			real space mask.
@param 	*pfsmask		cross-correlation mask.
@param 	max_iter		maximum number of iterations, 0 means only positional refinement.
@param 	viewstep		first view direction angular step size (radians).
@param 	rotstep			rotation around view angular step size (radians).
@param 	hires			high resolution limit for cross-correlation.
@param 	lores			low resolution limit for cross-correlation.
@param 	accuracy		angular accuracy (radians).
@param 	max_shift		maximum shift in coordinates (angstrom).
@param 	max_view_angle	maximum angular change in view vector (radians).
@param 	max_rot_angle	maximum angular change in rotation around view vector (radians).
@param 	shift_flag		flag to shift only along the normal (1) or perpendicular to it (2).
@return int				0, <0 on error.

	The density origin is positioned on the component.
	The component views must already be set.
	The number of component type names should be equal to the number of 
	sub-images in the template.

**/
int			model_refine_components(Bmodel* model, Bstring* ct_names, Bimage* ptemp,
				Bimage* pmask, Bimage* pfsmask, int max_iter, double viewstep, double rotstep,
				double hires, double lores, double accuracy, double max_shift,
				double max_view_angle, double max_rot_angle, int shift_flag)
{
	if ( max_iter < 0 ) max_iter = 0;
	if ( max_shift <= 0 ) max_shift = ptemp->sizeX()*ptemp->sampling(0)[0]/10.0;
//	if ( max_view_angle <= 0 ) max_view_angle = M_PI;
//	if ( max_rot_angle <= 0 ) max_rot_angle = M_PI;

//	double			alpha_step = M_PI/60.0;		// 3 degree starting step size
//	if ( alpha_step > max_view_angle ) alpha_step = max_view_angle;
	
//	double			alpha_step = hires/ptemp->real_size()[0];		// resolution-linked starting step size
//	if ( alpha_step > max_rot_angle ) alpha_step = max_rot_angle;
	
	long			i, nname;
	Bstring			base("C"), id, tempname(ptemp->file_name());
	Bstring*		name;
	Bcomptype*		ct;
	Bcomptype*		ct2;
	
	if ( !ct_names ) {
		for ( ct = model->type; ct && !ct->select(); ct = ct->next ) ;
		if ( ct ) base = ct->identifier();
		if ( ptemp->images() < 2 ) string_add(&ct_names, base);
		else for ( i=1; i<=ptemp->images(); i++ ) {
			id = base + Bstring(i, "%03d");
			string_add(&ct_names, id);
		}
	}
	
	for ( nname=0, name = ct_names; name; name = name->next ) nname++;
	
	if ( nname != ptemp->images() ) {
		cerr << "Error: The number of component type names must equal" << endl;
		cerr << "	the number of template sub-images (" << ptemp->images() << ")" << endl << endl;
		return -1;
	}

	// If a type is selected, it will be replaced
	for ( ct = model->type; ct && !ct->select(); ct = ct->next ) ;
	if ( ct ) {
		for ( ct2 = model->type; ct2; ct2 = ct2->next ) if ( ct2->select() ) ct2->select(-1);
		if ( ptemp->images() == 1 ) {
			ct->identifier() = ct_names->str();
			ct->select(1);
		} else {
			for ( i=0, name = ct_names; name; name = name->next, i++ )
				ct = model->add_type(name->str(), tempname.str(), i);
			for ( ct = ct2 = model->type; ct; ) {
				if ( ct->select() < 0 ) {
					if ( ct == model->type ) {
						model->type = ct2 = ct->next;
						delete ct;
						ct = ct2;
					} else {
						ct2->next = ct->next;
						delete ct;
						ct = ct2->next;
					}
				} else {
					ct2 = ct;
					ct = ct->next;
				}
			}
		}
	}

	random_seed();
	
	Bimage*			p = read_img(model->mapfile(), 1, model->image_number());
	if ( p == NULL ) {
		cerr << "Error: No model map file read! (" << model->mapfile() << ")" << endl;
		bexit(-1);
	}
	
	p->change_type(Float);
	ptemp->change_type(Float);
	p->calculate_background();
	
	if ( pmask ) {
		ptemp->rescale_to_avg_std(0, 1, pmask);
		ptemp->multiply(pmask);
	}
	
	if ( hires > lores ) swap(hires, lores);
	
	double			radius = ptemp->sizeX()/2;

	View2<float>	ref_view;
	list<View2<float>>	views = views_within_limits2(ref_view, viewstep, viewstep, rotstep, max_view_angle, max_rot_angle);
	long			nv = views.size();
//	long			nv = count_list((char *) views);
//	kill_list((char *) views, sizeof(View));

	if ( verbose ) {
		cout << "Refining component positions and views:" << endl;
		cout << "Component types:                ";
		for ( name = ct_names; name; name = name->next ) cout << *name << " ";
		cout << endl;
		cout << "Reference map:                  " << ptemp->file_name() << endl;
		if ( pmask ) cout << "Mask:                           " << pmask->file_name() << endl;
		if ( pfsmask ) cout << "Fourier space mask:             " << pfsmask->file_name() << endl;
		cout << "Iterations:                     " << max_iter << endl;
		cout << "Resolution limits:              " << hires << " " << lores << " A" << endl;
		cout << "Maximum shift:                  " << max_shift << " A" << endl;
		cout << "Maximum view rotation:          " << max_view_angle*180.0/M_PI << " degrees" << endl;
		cout << "Maximum angle rotation:         " << max_rot_angle*180.0/M_PI << " degrees" << endl;
		cout << "Step sizes:                     " << viewstep*180.0/M_PI << "  "
			<< rotstep*180.0/M_PI << " degrees" << endl;
		cout << "Accuracy:                       " << accuracy*180.0/M_PI << " degrees" << endl;
		cout << "Number of views:                " << nv << endl;
		cout << "Density radius:                 " << radius << endl << endl;
	}
	
	long			n, ncomp;
	double			s, a, d, d1, d2, ccavg(0), d1avg(0), d2avg(0);
	Vector3<double>	map_origin(p->image->origin());
	Vector3<double>	scale(ptemp->sampling(0));
	Vector3<double>	loc, comp_shift, shift_avg;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_refine_components: map_origin=" << map_origin << " scale=" << scale << endl;

	Bcomponent*		comp;
	Blink*			link;
	int*			num = new int[ptemp->images()];
	for ( n=0; n<ptemp->images(); n++ ) num[n] = 0;

	for ( ncomp=0, comp = model->comp; comp; comp = comp->next )
		if ( comp->select() ) ncomp++;

	Bcomponent**	comparr = new Bcomponent*[ncomp];

	for ( i=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		comparr[i++] = comp;
		loc = comp->location()/scale + map_origin;
		comp->density(p->density(0, loc, radius));
	}
	
	fft_plan		planf = ptemp->fft_setup(FFTW_FORWARD, 1);
	fft_plan		planb = ptemp->fft_setup(FFTW_BACKWARD, 1);
	
	if ( verbose )
		cout << "Comp\tType\tdx\tdy\tdz\tCC" << endl;
#ifdef HAVE_GCD
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(ncomp, dispatch_get_global_queue(0, 0), ^(size_t i){
		if ( max_iter )
			component_refine_random(comparr[i], p, ptemp, pmask, pfsmask, max_iter, hires, lores,
					max_shift, max_view_angle, max_rot_angle, shift_flag, planf, planb);
		else
			component_refine(comparr[i], p, ptemp, pmask, pfsmask, max_view_angle,
					max_rot_angle, viewstep, rotstep, hires, lores, accuracy,
					max_shift, shift_flag, planf, planb);
		if ( verbose )
			dispatch_sync(myq, ^{
				cout << comparr[i]->identifier() << tab << comparr[i]->select() << tab
					<< setprecision(2) << comparr[i]->force()
					<< tab << setprecision(4) << comparr[i]->FOM() << endl;
			});
	});
#else
#pragma omp parallel for
	for ( i=0; i<ncomp; i++ ) {
		if ( max_iter )
			component_refine_random(comparr[i], p, ptemp, pmask, pfsmask, max_iter, hires, lores,
					max_shift, max_view_angle, max_rot_angle, shift_flag, planf, planb);
		else
			component_refine(comparr[i], p, ptemp, pmask, pfsmask, max_view_angle,
					max_rot_angle, viewstep, rotstep, hires, lores, accuracy,
					max_shift, shift_flag, planf, planb);
	#pragma omp critical
		if ( verbose )
				cout << comparr[i]->identifier() << tab << comparr[i]->select() << tab
					<< setprecision(2) << comparr[i]->force()
					<< tab << setprecision(4) << comparr[i]->FOM() << endl;
	}
#endif

	if ( verbose )
		cout << "Comp\tType\tdx\tdy\tdz\tShift\tCC\tDens1\tDens2\tdDens" << endl;
	for ( i=0; i<ncomp; i++ ) {
		comp = comparr[i];
		s = comp->force().length();
		num[comp->select()-1]++;
		loc = comp->location()/scale + map_origin;
		shift_avg += comp->force();
		ccavg += comp->FOM();
		d1 = comp->density();
//		d2 = img_density_at_coordinates(p, 0, loc, radius);
		d2 = p->density(0, loc, radius);
		d1avg += d1;
		d2avg += d2;
		comp->density(d2);
		for ( n=1, name = ct_names; n<comp->select() && n<ptemp->images(); n++, name = name->next ) ;
//		if ( name ) comp->type = model_add_type_by_id(model, *name);
//		if ( name ) comp->type(model->add_type(*name));
		if ( verbose )
			cout << comp->identifier() << tab << comp->select() << tab
				<< setprecision(2) << comp->force()
				<< tab << s << tab
				<< setprecision(4) << comp->FOM() << tab << d1 << tab << d2 << tab << d2 - d1 << endl;
	}

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	shift_avg /= ncomp;
	s = shift_avg.length();
	ccavg /= ncomp;
	d1avg /= ncomp;
	d2avg /= ncomp;
	
	if ( verbose ) {
		cout << "Average\t\t" << shift_avg[0] << tab << shift_avg[1] << tab << 
			shift_avg[2] << tab << s << tab << ccavg << tab <<
			d1avg << tab << d2avg << tab << d2avg - d1avg << endl;
		cout << "Component type distribution:\n#\tType\tCount" << endl;
		for ( n=0, name = ct_names; n<ptemp->images(); n++, name = name->next )
			 cout << n+1 << tab << *name << tab << num[n] << endl;
		cout << endl;
		if ( model->link ) {
			for ( n=0, a=s=0, link = model->link; link; link = link->next, n++ ) {
				d = link->comp[0]->location().distance(link->comp[1]->location());
				a += d;
				s += d*d;
			}
			if ( n ) {
				a /= n;
				s = s/n - a*a;
				if ( s > 0 ) s = sqrt(s);
				else s = 0;
			}
			cout << "Average link length:            " << a << " A (" << s << ")" << endl << endl;
		}
	}
	
	delete[] num;
	delete p;

	return 0;
}

/**
@brief 	Refines link positions by cross-correlation.
@param 	*model			model.
@param 	*ptemp			density template.
@param 	*pmask			real space mask.
@param 	*pfsmask		cross-correlation mask.
@param 	hires			high resolution limit for cross-correlation.
@param 	lores			low resolution limit for cross-correlation.
@param 	max_shift		maximum shift in coordinates (angstrom).
@param 	shift_flag		flag to shift only along the normal (1) or perpendicular to it (2).
@param 	bias			bias to apply to first correlation coefficient.
@return Bimage*			density around the component.

	The density origin is positioned on the link center.
	The component views must already be set.

**/
int			model_refine_link_positions(Bmodel* model, Bimage* ptemp, Bimage* pmask, 
				Bimage* pfsmask, double hires, double lores, double max_shift, int shift_flag, double bias)
{
	if ( max_shift <= 0 ) max_shift = 1e37;
	
	random_seed();
	
	Bimage*			p = read_img(model->mapfile(), 1, model->image_number());
	if ( p == NULL ) {
		cerr << "Error: No model map file read! (" << model->mapfile() << ")" << endl;
		bexit(-1);
	}
	
	p->change_type(Float) ;
	ptemp->change_type(Float) ;

	if ( pmask ) ptemp->multiply(pmask);
	
	if ( hires > lores ) swap(hires, lores);
	
	double			radius = ptemp->sizeX()/2;

	if ( verbose ) {
		cout << "Refining link positions:" << endl;
		cout << "Resolution limits:              " << hires << " " << lores << " A" << endl;
		cout << "Maximum shift:                  " << max_shift << " A" << endl;
		cout << "Density radius:                 " << radius << endl;
		cout << "Bias:                           " << bias << endl << endl;
	}
	
	long			n, nlink, imax(0);
	Vector3<long>	size(ptemp->size());
	Vector3<double>	origin(ptemp->image->origin());
	Vector3<double>	map_origin(p->image->origin());
	Vector3<double>	scale(ptemp->sampling(0));
	Vector3<double>	link_shift, shift_avg, link_center;
	Vector3<double>	shift;
	Bimage*			plink = NULL;
	Bimage*			pone = NULL;
	Bimage*			pmaskrot = NULL;

	Blink*			link = NULL;
	Matrix3			mat;
	double			ccmax, s, d1, d2, ccavg(0), d1avg(0), d2avg(0);
	int*			num = new int[ptemp->images()];
	double*			cc = new double[ptemp->images()];
	for ( n=0; n<ptemp->images(); n++ ) cc[n] = num[n] = 0;

	fft_plan		planf = ptemp->fft_setup(FFTW_FORWARD, 1);
	fft_plan		planb = ptemp->fft_setup(FFTW_BACKWARD, 1);
		
	if ( verbose ) {
		cout << "Bond\tdx\tdy\tdz\tShift\tDens1\tDens2\tdDens";
		for ( n=0; n<ptemp->images(); n++ ) cout << "\tCC" << n+1;
		cout << endl;
	}
	for ( nlink=0, link = model->link; link; link = link->next, nlink++ ) if ( link->select() ) {
		link_center = (link->comp[0]->location() + link->comp[1]->location())*0.5;
		d1 = p->density(0, link_center/scale + map_origin, radius);
		mat = get_link_matrix(link);
		plink = p->extract(0, map_origin + link_center/p->sampling(0)[0], size, origin, mat);
		if ( pmask ) plink->multiply(pmask);
		if ( pfsmask ) pmaskrot = pfsmask->extract_wrap(0, size, mat);
		ccmax = -1;
		for ( n=0; n<ptemp->images(); n++ ) {
			pone = ptemp->extract(n);
			pone->origin(origin);
			shift = plink->find_shift(pone, pmaskrot, hires, lores, ptemp->sizeX()/10.0, 0, 1, planf, planb, cc[n]);
			if ( n==0 ) cc[n] *= bias;
			if ( ccmax < cc[n] ) {
				ccmax = cc[n];
				imax = n;
				if ( shift_flag == 1 ) shift[0] = shift[1] = 0;	// Shift only along the normal
				else if ( shift_flag == 2 ) shift[2] = 0;		// Shift only lateral
				link_shift = scale * (mat * shift);
			}
			delete pone;
		}
		delete plink;
		delete pmaskrot;
		num[imax]++;
		s = link_shift.length();
		if ( s > max_shift ) {
			link_shift *= max_shift/s;
			s = max_shift;
		}
		link->comp[0]->shift(link_shift);
		link->comp[1]->shift(link_shift);
		link_center = (link->comp[0]->location() + link->comp[1]->location())*0.5;
		link->select(imax);
		link->FOM(ccmax);
		shift_avg += link_shift;
		ccavg += ccmax;//		d2 = img_density_at_coordinates(p, 0, link_center/scale + map_origin, radius);
		d2 = p->density(0, link_center/scale + map_origin, radius);
		d1avg += d1;
		d2avg += d2;
		if ( verbose ) {
			cout << nlink+1 << tab << link_shift[0] << tab << link_shift[1] << tab << 
				link_shift[2] << tab << s << tab << d1 << tab << d2 << tab << d2 - d1;
			for ( n=0; n<ptemp->images(); n++ ) cout << tab << cc[n];
			cout << endl;
		}
	}

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	shift_avg /= nlink;
	s = shift_avg.length();
	ccavg /= nlink;
	d1avg /= nlink;
	d2avg /= nlink;

	delete[] cc;
	delete[] num;
	delete p;

	if ( verbose ) {
		cout << "Average\t" << shift_avg[0] << tab << shift_avg[1] << tab << 
			shift_avg[2] << tab << s << tab << d1avg << tab << d2avg << tab << 
			d2avg - d1avg << tab << ccavg << endl << endl;
		cout << "Distribution:";
		for ( n=0; n<ptemp->images(); n++ ) cout << tab << num[n];
		cout << endl << endl;
	}
	
	return 0;
}

/**
@brief 	Averages the density associated with each component type in a model.
@param 	*model		model.
@param 	size		size of component density to extract.
@param 	origin		origin of new component density image.
@param	npt			number per type.
@return Bimage*		average density around a component.

	The densities associated with each vertex type is extracted and averaged.
	The extracted density origin is placed on the component coordinates.
	The component views must already be set.
	Only the first model in the list is processed.
	The map file for the model must exist.

**/
Bimage*		model_average_component_density(Bmodel* model, Vector3<long> size, Vector3<double> origin, int npt)
{
	if ( origin[0] == 0 && origin[1] == 0 && origin[2] == 0 ) 
		origin = size/2;

	if ( npt < 1 ) npt = 1;
	
	long			i, j, n, navg, nc(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomptype*		ct = NULL;
	Bimage*			p = NULL;
	Bstring*		comptypelist = NULL;
	Bstring*		comptype = NULL;
	Bstring*		temptype = NULL;

	// Compile a list of component types
	for ( navg=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( ct = mp->type; ct; ct = ct->next ) if ( ct->select() ) {
			for ( temptype = comptypelist; temptype && *temptype != ct->identifier(); temptype = temptype->next ) ;
			if ( !temptype ) {
				comptype = string_add(&comptype, ct->identifier().c_str());
				if ( !comptypelist ) comptypelist = comptype;
				navg++;
			}
		}
	}
	
	// Associate a consistent image number to each component type
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( ct = mp->type; ct; ct = ct->next ) if ( ct->select() ) {
			for ( n=0, temptype = comptypelist; temptype && *temptype != ct->identifier(); temptype = temptype->next ) n++;
			ct->image_number(n);
		}
	}

	string_kill(comptypelist);
	
	if ( verbose ) {
		cout << "Averaging component densities:" << endl;
		cout << "Number of selected types:       " << navg << endl;
		cout << "Number of averages per type:    " << npt << endl;
		cout << "Size:                           " << size << endl;
		cout << "Origin:                         " << origin << endl;
	}
	
	navg *= npt;
	
	int*			ncomp = new int[navg];
	for ( n=0; n<navg; n++ ) ncomp[n] = 0;
	
	Bimage*			pcomp = new Bimage(Float, TSimple, size, navg);
	pcomp->origin(origin);
	
	long			imgsize = pcomp->sizeX()*pcomp->sizeY()*pcomp->sizeZ();
	double			radius = origin[0];
	Vector3<double>	imgori, loc;
	Matrix3			mat(1);
	Bimage*			pone = NULL;
	pcomp->next = new Bimage(Float, TSimple, pcomp->size(), pcomp->images());
	float*			fom = (float *) pcomp->next->data_pointer();
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		p = read_img(mp->mapfile(), 1, mp->image_number());
		if ( p == NULL ) {
			cerr << "Error: No model map file read! (" << mp->mapfile() << ")" << endl;
			bexit(-1);
		}
		pcomp->sampling(p->sampling(0));
		p->change_type(Float) ;
		p->calculate_background();
		imgori = p->image->origin();
		if ( verbose & VERB_PROCESS ) {
			cout << "Model:                          " << mp->identifier() << endl;
			cout << "Map:                            " << mp->mapfile() << " (" << mp->image_number() << ")" << endl;
			cout << "Comp\tAvg#\tx\ty\tz\tvx\tvy\tvz\tva" << endl;
		}
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			ct = comp->type();
			n = ct->image_number() * npt + nc%npt;
			ncomp[n]++;
			if ( verbose & VERB_PROCESS )
				cout << comp->identifier() << tab << n+1 << tab << comp->location() << tab << comp->view() << endl;
			mat = comp->view().matrix();
			loc = imgori + comp->location()/p->sampling(0);
			comp->density(p->density(0, loc, radius));
			pone = p->extract(0, loc, size, origin, View(comp->view()[0],comp->view()[1],comp->view()[2],comp->view()[3]));
			for ( i=n*imgsize, j=0; j<imgsize; i++, j++ ) {
				pcomp->add(i, (*pone)[j]);
				fom[i] += (*pone)[j]*(*pone)[j];
			}
			delete pone;
			nc++;
		}
		delete p;
	}

	for ( i=n=0; n<navg; n++ ) if ( ncomp[n] ) {
		pcomp->multiply(n, 1.0L/ncomp[n]);
		for ( j=0; j<imgsize; i++, j++ ) {
			fom[i] = fom[i]/ncomp[n] - (*pcomp)[i]*(*pcomp)[i];
			if ( fom[i] > 0 ) fom[i] = sqrt(fom[i]);
			else fom[i] = 0;
		}
	}
	
	if ( verbose ) {
		cout << "Component types:\n#\tType\tCount" << endl;
		for ( n=0, ct = model->type; ct; ct = ct->next ) if ( ct->select() ) {
			for ( i=0; i<npt; i++, n++ )
				cout << n+1 << tab << ct->identifier() << tab << ncomp[n] << endl;
		}
		cout << endl;
	}
	
	delete[] ncomp;

	return pcomp;
}

/**
@brief 	Extracts all densities associated with components in a model.
@param 	*model		model.
@param 	size		size of component density to extract.
@param 	origin		origin of new component density image.
@return Bimage*		all densities around components.

	The densities associated with each component is extracted.
	The extracted density origin is placed on the component coordinates.
	The component views must already be set.
	Only the first model in the list is processed.
	The map file for the model must exist.

**/
Bimage*		model_extract_component_densities(Bmodel* model, Vector3<long> size, Vector3<double> origin)
{
	Bimage*			p = read_img(model->mapfile(), 1, model->image_number());
	if ( p == NULL ) {
		cerr << "Error: No model map file read! (" << model->mapfile() << ")" << endl;
		bexit(-1);
	}
	
	p->change_type(Float) ;
	p->calculate_background();

	long			ncomp(0);

	if ( origin[0] == 0 && origin[1] == 0 && origin[2] == 0 ) 
		origin = Vector3<double>(size/2);

	Bcomponent**	carr = component_get_array(model, ncomp);
	
	if ( verbose ) {
		cout << "Extracting component densities:" << endl;
		cout << "Number of selected components:  " << ncomp << endl;
		cout << "Size:                           " << size << endl;
		cout << "Origin:                         " << origin << endl;
	}
	
	double			radius = origin[0];

	Bimage*			pcomp = new Bimage(Float, p->compound_type(), size, ncomp);
	pcomp->origin(origin);
	pcomp->sampling(p->sampling(0));
	
	Vector3<double>	imgori(p->image->origin());

	if ( verbose & VERB_PROCESS ) {
		cout << "Model:                          " << model->identifier() << endl;
		cout << "Map:                            " << model->mapfile() << " (" << model->image_number() << ")" << endl;
		cout << "Map origin:                     " << p->image->origin() << endl;
		cout << "Comp\tx\ty\tz\tvx\tvy\tvz\tva" << endl;
	}
	
#ifdef HAVE_GCD
	dispatch_apply(ncomp, dispatch_get_global_queue(0, 0), ^(size_t i){
		Bcomponent*		comp = carr[i];
		Vector3<double>	loc(comp->location()/p->sampling(0) + imgori);
		comp->density(p->density(0, loc, radius));
		Bimage*			pone = p->extract(0, loc, size, origin, View(comp->view()[0],comp->view()[1],comp->view()[2],comp->view()[3]));
		pcomp->replace(i, pone);
		delete pone;
		if ( verbose & VERB_PROCESS )
			cout << comp->identifier() << tab << loc << tab << comp->view() << endl;
	});
#else
#pragma omp parallel for
	for ( long i=0; i<ncomp; i++ ) {
		Bcomponent*		comp = carr[i];
		Vector3<double>	loc(comp->location()/p->sampling(0) + imgori);
		comp->density(p->density(0, loc, radius));
		Bimage*			pone = p->extract(0, loc, size, origin, View(comp->view()[0],comp->view()[1],comp->view()[2],comp->view()[3]));
		pcomp->replace(i, pone);
		delete pone;
		if ( verbose & VERB_PROCESS )
			cout << comp->identifier() << tab << loc << tab << comp->view() << endl;
	}
#endif

	delete p;
	delete[] carr;

	return pcomp;
}

/**
@brief 	Extracts a density associated with each link in a model.
@param 	*model		model.
@param 	size		size of link density to extract.
@param 	origin		origin of new link density image.
@return Bimage*		new image with the density around the link.

	The link density origin is positioned on the center of the link.
	The component views must already be set.

**/
Bimage*		model_average_link_density(Bmodel* model, Vector3<long> size, Vector3<double> origin)
{
	if ( origin[0] == 0 && origin[1] == 0 && origin[2] == 0 ) 
		origin = Vector3<double>(size[0]/2, size[1]/2, size[2]/2);

	long	i, j, n, navg;
	Bmodel*			mp;
	Blink*			link;
	Matrix3			mat;
	Bimage*			p = NULL;

	for ( navg=0, link = model->link; link; link = link->next )
		if ( link->select() > navg ) navg = link->select();

	if ( verbose ) {
		cout << "Averaging link density" << endl << endl;
		cout << "Number of averages:             " << navg << endl;
		cout << "Size:                           " << size << endl;
		cout << "Origin:                         " << origin << endl;
	}
	
	int*			nlink = new int[navg+1];
	for ( n=0; n<=navg; n++ ) nlink[n] = 0;
	
	Bimage*			plink = new Bimage(Float, TSimple, size, navg);
	plink->origin(origin);
	
	long			imgsize = plink->sizeX()*plink->sizeY()*plink->sizeZ();
	Vector3<double>	imgori, loc;
	Bimage*			pone = NULL;
	plink->next = new Bimage(Float, TSimple, plink->size(), plink->images());
	float*			fom = (float *) plink->next->data_pointer();
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		p = read_img(mp->mapfile(), 1, mp->image_number());
		if ( p == NULL ) {
			cerr << "Error: No model map file read! (" << mp->mapfile() << ")" << endl;
			bexit(-1);
		}
		plink->sampling(p->sampling(0));
		p->change_type(Float) ;
		imgori = p->image->origin();
		if ( verbose & VERB_PROCESS )
			cout << "Comp1\tComp2\tAvg#" << endl;
		for ( link = mp->link; link; link = link->next ) if ( link->select() ) {
			n = (long) link->select();
			nlink[n]++;
			if ( n ) {
				if ( verbose & VERB_PROCESS )
					cout << link->comp[0]->identifier() << tab << link->comp[1]->identifier() << tab << n << endl;
				mat = get_link_matrix(link);
				loc = imgori + (link->comp[1]->location() + link->comp[0]->location())/(2*p->sampling(0)[0]);
//				pone = img_extract_density(p, size, origin, loc, mat);
				pone = p->extract(0, loc, size, origin, mat);
				for ( i=(n-1)*imgsize, j=0; j<imgsize; i++, j++ ) {
					plink->add(i, (*pone)[j]);
					fom[i] += (*pone)[j]*(*pone)[j];
				}
				delete pone;
			}
		}
		delete p;
	}

	for ( i=n=0; n<navg; n++ ) if ( nlink[n] ) {
		plink->multiply(n, 1.0L/nlink[n]);
		for ( j=0; j<imgsize; i++, j++ ) {
			fom[i] = fom[i]/nlink[n] - (*plink)[i]*(*plink)[i];
			if ( fom[i] > 0 ) fom[i] = sqrt(fom[i]);
			else fom[i] = 0;
		}
	}
	
	delete[] nlink;
	
	return plink;
}

/**
@brief 	Builds a new map from a densities of components in a model.
@param 	*model		model.
@param 	size		size of new map.
@param 	origin		origin of new map with respect to the model.
@param 	flags		flags to weigh by contributions (1) and build separate maps (2).
@return Bimage*		new map.

	The number of new maps depends on the number of selected component types
	and access to their density maps.
	The component views must already be set.
	The sampling must be the same for all component type maps.

**/
Bimage*		model_build_from_component_density(Bmodel* model, Vector3<long> size, 
				Vector3<double> origin, int flags)
{
	long			i, ntype(1), nmap(1), nb;
	Bcomptype*		ct = NULL;
	Bcomponent*		comp;
	Bimage*			p = NULL;
	Vector3<double>	sampling;
	
	// Set up for the number of component types and their image numbers
	for ( ntype = 0, ct = model->type; ct; ct = ct->next )
		if ( ct->select() ) ntype++;
	for ( ct = model->type; ct; ct = ct->next )
		if ( ct->select() ) break;
	if ( flags & 2 ) nmap = ntype;

	if ( !ct ) {
		cerr << "Error: No component types defined for model " << model->identifier() << endl;
		bexit(-1);
	}
	
	if ( file_type(ct->file_name().c_str() ) == Image ) {
		p = read_img(ct->file_name(), 0, ct->image_number());
		if ( !p ) {
			cerr << "Error: No component map file read! (" << ct->file_name()  << ")" << endl;
			bexit(-1);
		}
		sampling = p->sampling(0);
		delete p;
		p = NULL;
	}

	if ( origin[0] == 0 && origin[1] == 0 && origin[2] == 0 ) 
		origin = size/2;

	Bimage*			pmap = new Bimage(Float, TSimple, size, nmap);
	pmap->origin(origin);
	pmap->sampling(sampling);
	
	if ( verbose ) {
		cout << "Building new maps from component densities:" << endl;
		cout << "Component types selected:       " << ntype << endl;
		cout << "Number of maps to build:        " << pmap->images() << endl;
		cout << "Size:                           " << pmap->size() << endl;
		cout << "Origin:                         " << origin << endl;
	}
	
	float*			fom = NULL;
	if ( flags & 1 ) {
		pmap->next = new Bimage(Float, TSimple, pmap->size(), pmap->images());
		fom = (float *) pmap->next->data_pointer();
	}
	
	long			x, y, z, n, ncomp, nlink;
 	double			radius, linklength, weight;
	Vector3<double>	v, loc, ori_comp, ori_map(pmap->image->origin());
	Matrix3			mat;
	
	for ( n = 0, ct = model->type; ct; ct = ct->next )
				if ( ct->select() && file_type(ct->file_name().c_str() ) == Image ) {
		p = read_img(ct->file_name() , 1, ct->image_number());
		if ( !p ) {
			cerr << "Error: No component map file read! (" << ct->file_name()  << ")" << endl;
			bexit(-1);
		}
		p->change_type(Float) ;
		if ( verbose & VERB_PROCESS )
			cout << "Comp\tType\tFile\tImage" << endl;
		for ( ncomp=nb=0, comp = model->comp; comp; comp = comp->next, ncomp++ ) 
				if ( comp->select() && comp->type() == ct ) {
			if ( verbose & VERB_PROCESS )
				cout << comp->identifier() << tab << n+1 << tab << ct->file_name()  << tab << ct->image_number() << endl;
			ori_comp = ori_map + comp->location()/p->sampling(0)[0];
			mat = comp->view().matrix();
			for ( i=0, nlink=0, linklength=0; comp->link[i]; i++ ) {
				linklength += comp->location().distance(comp->link[i]->location());
				nlink++;
			}
			if ( nlink ) linklength /= 2*nlink;
			for ( z = 1; z < p->sizeZ()-1; z++ ) {
				v[2] = z - p->image->origin()[2];
				for ( y = 1; y < p->sizeY()-1; y++ ) {
					v[1] = y - p->image->origin()[1];
					for ( x = 1; x < p->sizeX()-1; x++ ) {
						v[0] = x - p->image->origin()[0];
						loc = (mat * v) + ori_comp;
						if ( loc[0]>=0 && loc[0]<=pmap->sizeX()-2 && loc[1]>=0 && loc[1]<=pmap->sizeY()-2 && loc[2]>=0 && loc[2]<=pmap->sizeZ()-2 ) {
							i = p->index(0, x, y, z, 0);
							weight = 1;
							if ( nlink ) {
								radius = v.length() - linklength;
								if ( radius > 0 ) weight = 1 - 0.002*radius*radius;
							}
							if ( weight > 0 ) {
								if ( flags & 1 )
									pmap->add(loc[0], loc[1], loc[2], n, weight*(*p)[i]);
								else
									pmap->set_max(loc[0], loc[1], loc[2], n, weight*(*p)[i]);
								if ( fom ) fom[pmap->index(loc, n)] += weight;
							}
						}
					}
				}
			}
			nb++;
		}
		delete p;
		if ( flags & 2 ) n++;
		if ( verbose )
			cout << "Built " << nb << " components of type " << ct->identifier() << " into map" << endl;
	}

	if ( flags & 1 ) for ( i=0; i<pmap->data_size(); i++ )
		if ( fom[i] ) pmap->set(i, (*pmap)[i]/fom[i]);
	
	return pmap;
}

/**
@brief 	Builds a new map from a density of a link in a model.
@param 	*model			model.
@param 	&linkmap		link map filename.
@param 	size			size of new map.
@param 	origin			origin of new map with respect to the model.
@param 	link_select		link selection number to build (first = 1).
@param 	flags			flags to weigh by contributions (1) and build separate maps (2).
@return Bimage*			new map.

	The component views must already be set.

**/
Bimage*		model_build_from_link_density(Bmodel* model, Bstring& linkmap, 
				Vector3<long> size, Vector3<double> origin, int link_select, int flags)
{
	Bimage*			p = read_img(linkmap, 1, -1);
	if ( p == NULL ) {
		cerr << "Error: No link map file read! (" << model->type->file_name()  << ")" << endl;
		bexit(-1);
	}
	
	p->change_type(Float) ;

	if ( origin[0] == 0 && origin[1] == 0 && origin[2] == 0 ) 
		origin = Vector3<double>(size[0]/2, size[1]/2, size[2]/2);

	long	i;
	
	Bimage*			pmap = new Bimage(Float, p->compound_type(), size[0], size[1], size[2], p->images());
	pmap->origin(origin);
	pmap->sampling(p->sampling(0));
	
	if ( verbose ) {
		cout << "Building a new map from a link density" << endl << endl;
		cout << "Size:                           " << pmap->size() << endl;
		cout << "Origin:                         " << origin << endl;
	}
	
	pmap->next = new Bimage(Float, TSimple, pmap->size(), pmap->images());
	float*			fom = (float *) pmap->next->data_pointer();
	
	long			x, y, z, n, xmin, xmax, nlink;
 	double			linklength, fx, weight;
	Vector3<double>	v, loc, ori_link;
	Matrix3			mat;
	Blink*			link;

	if ( verbose & VERB_PROCESS )
		cout << "Comp1\tComp2\tMap" << endl;
	for ( nlink=0, link = model->link; link; link = link->next, nlink++ ) if ( link->select() ) {
		n = (long) link->select() - 1;
		if ( n >= p->images() ) n = p->images() - 1;
		if ( link_select < 1 || link_select == n+1 ) {
			if ( verbose & VERB_PROCESS )
				cout << link->comp[0]->identifier() << tab << link->comp[1]->identifier() << tab << n+1 << endl;
			linklength = link->comp[1]->location().distance(link->comp[0]->location());
			ori_link = origin + (link->comp[1]->location() + link->comp[0]->location())/(2*p->sampling(0)[0]);
			xmin = (long) (p->image->origin()[0] - linklength/(2*p->sampling(0)[0]) + 0.5);
			xmax = (long) (p->image->origin()[0] + linklength/(2*p->sampling(0)[0]) + 0.5);
			mat = get_link_matrix(link);
			for ( z = 1; z < p->sizeZ()-1; z++ ) {
				v[2] = z - p->image->origin()[2];
				for ( y = 1; y < p->sizeY()-1; y++ ) {
					v[1] = y - p->image->origin()[1];
					for ( x = 1; x < p->sizeX()-1; x++ ) {
						v[0] = x - p->image->origin()[0];
						loc = (mat * v) + ori_link;
						if ( loc[0]>=0 && loc[0]<=pmap->sizeX()-2 && loc[1]>=0 && loc[1]<=pmap->sizeY()-2 && loc[2]>=0 && loc[2]<=pmap->sizeZ()-2 ) {
							i = p->index(0, x, y, z, n);
							weight = 1;
							fx = 0;
							if ( x < xmin ) fx = xmin - x;
							else if ( x > xmax ) fx = x - xmax;
							if ( fx ) weight = 1 - 0.002*fx*fx;
							if ( weight > 0 )
								pmap->add(loc[0], loc[1], loc[2], 0, weight*(*p)[i]);
						}
					}
				}
			}
		}
	}

	if ( flags & 1 ) for ( i=0; i<pmap->data_size(); i++ )
		if ( fom[i] ) pmap->set(i, (*pmap)[i]/fom[i]);
	
	delete p;

	return pmap;
}

