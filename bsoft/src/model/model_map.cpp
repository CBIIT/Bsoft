/**
@file	model_map.cpp
@brief	Function to generate a map from a model.
@author Bernard Heymann
@date	Created: 20081112
@date	Modified: 20220419
**/

#include "rwmodel.h"
#include "rwmodel_param.h"
#include "model_links.h"
#include "model_transform.h"
#include "model_extract_build.h"
#include "model_mol.h"
#include "model_util.h"
#include "rwimg.h"
#include "img_combine.h"
#include "scatter.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Creates one or more models from 3D images.
@param 	*plist			list of images.
@return Bmodel*			new model.

	The model ID's are derived from the image filenames.

**/
Bmodel*		model_from_images(Bimage* plist)
{
	long			n;
	Bimage*			p = plist;
	Bstring			base, id, filename(p->file_name());
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	
	base = filename.base().alnum();
	
	for ( n=0; n<p->images(); n++ ) {
		id = base;
		if ( p->images() > 1 ) id += Bstring(n+1, "_%04d");
//		mp = model_add(&mp, id);
//		if ( !model ) model = mp;
		if ( mp ) mp = mp->add(id.str());
		else model = mp = new Bmodel(id.str());
		mp->mapfile() = p->file_name();
		mp->image_number(n);
	}
	
	return model;
}

/**
@brief 	Creates one or more models from graph segments.
@param 	*p				image.
@param 	&gs				graph segments.
@return Bmodel*			new model.

	The regions are assumed to be number consecutively.

**/
Bmodel*		model_from_graph_segments(Bimage* p, GSgraph& gs)
{
	long			i, j, k(0), nr(gs.voxel_count());
	double			maxlen(1.1*p->image->sampling()[0]), f(100.0/p->image_size());
	Bstring			id;
	vector<Bmodel*>	marr(nr);
	vector<Bcomponent*>	carr(nr);
	
	if ( verbose )
		cout << "Converting " << nr << " regions to models" << endl;
	
	Bmodel*			model = new Bmodel(1);
	Bmodel*			mp = model;
	for ( i=1; i<nr; ++i ) mp = mp->add(i+1);
	
	for ( i=0, mp = model; mp && i<nr; ++i, mp = mp->next ) {
		marr[i] = mp;
		carr[i] = mp->comp;
	}
	
	for ( i=0; i<p->image_size(); ++i ) {
		k = gs.region_find(i);
		j = gs.voxel(i).joins() - 1;
//		cout << i << tab << j << tab << gs.voxel(k).voxels() << endl;
		if ( gs.voxel(k).voxels() < 10000 ) {
			if ( carr[j] ) {
				id = carr[j]->identifier();
				k = id.integer() + 1;
				carr[j] = carr[j]->add(k);
			} else {
				carr[j] = marr[j]->add_component(1);
				k = 1;
			}
//			cout << i << tab << j << tab << k << endl;
			carr[j]->location(p->real_coordinates(i));
			carr[j]->density((*p)[i]);
		}
		if ( ( i%1000 == 0 ) && verbose )
			cerr << "Complete:                       " << setprecision(3)
							<< f*i << " %    \r" << flush;
	}
	
	if ( verbose )
		cout << "Model\tSize\tLinks" << endl;
	for ( mp = model; mp; mp = mp->next ) {
		mp->mapfile(p->file_name());
		model_link_list_generate(mp, maxlen);
		cout << mp->identifier() << tab << mp->component_count() << tab
			<< mp->link_count() << endl;
	}

	return model;
}

/**
@brief 	Calculates a map from model component coordinates.  
@param 	*model		model.
@param 	ori			origin within new map.
@param 	size		size of new map.
@param 	sam			voxel size of new map.
@param 	sigma		gaussian sphere width.
@return Bimage*		new map.
**/
Bimage*		img_from_model(Bmodel* model, Vector3<double> ori,
				Vector3<long> size, Vector3<double> sam, double sigma)
{
	if ( !model ) {
		cerr << "Error: No model selected!" << endl << endl;
		return NULL;
	}

	if ( sam.volume() < 1e-6 )
		sam[0] = sam[1] = sam[2] = 1;
	else
		sam = sam.max(0.01);
	
	if ( sigma < sam[0] ) sigma = sam[0];
	
	long			i, n;
	long			x, y, z, xlo, ylo, zlo, xhi, yhi, zhi;
	double			x2, y2, z2;
	double			pix_sigma = sigma/sam[0];
	double			pix_rad = 3*pix_sigma;
	double			nhis2 = -0.5/(pix_sigma*pix_sigma);
	Vector3<double>	min, max, c;
	Bmodel*			mp;
	Bcomponent*		comp;
	
	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		mp->image_number(n++);
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			min = min.min(comp->location());
			max = max.max(comp->location());
		}
	}
	
	if ( size[0] < 1 ) size[0] = (long) ((max[0] - min[0] + 3*sigma)/sam[0]);
	if ( size[1] < 1 ) size[1] = (long) ((max[1] - min[1] + 3*sigma)/sam[1]);
	if ( size[2] < 1 ) size[2] = (long) ((max[2] - min[2] + 3*sigma)/sam[2]);
	size = size.min(1024);
	size = size.max(1);

	Bimage* 	p = new Bimage(Float, TSimple, size, n);
	p->file_name(model->identifier());
	p->origin(ori);
	p->sampling(sam);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating density by simple gaussian expansion:" << endl;
	    cout << "Models:                         " << p->images() << endl;
	    cout << "Dimensions:                     " << p->size() << " voxels" << endl;
    	cout << "Origin:                         " << p->image->origin()<< endl;
		cout << "Sampling:                       " << p->sampling(0) << " A/voxel" << endl;
    	cout << "Sigma:                          " << sigma << " A (" << pix_sigma << " pixels)" << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << endl << "Calculating density" << endl << endl;

	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			c = comp->location()/sam + ori;
			xlo = (long) (c[0] - pix_rad + 0.5);
			ylo = (long) (c[1] - pix_rad + 0.5);
			zlo = (long) (c[2] - pix_rad + 0.5);
			if ( xlo < 0 ) xlo = 0;
			if ( ylo < 0 ) ylo = 0;
			if ( zlo < 0 ) zlo = 0;
			xhi = (long) (c[0] + pix_rad + 0.5);
			yhi = (long) (c[1] + pix_rad + 0.5);
			zhi = (long) (c[2] + pix_rad + 0.5);
			if ( xhi >= p->sizeX() ) xhi = p->sizeX() - 1;
			if ( yhi >= p->sizeY() ) yhi = p->sizeY() - 1;
			if ( zhi >= p->sizeZ() ) zhi = p->sizeZ() - 1;
			for ( z = zlo; z <= zhi; z++ ) {
				z2 = z - c[2];
				z2 *= z2;
				for ( y = ylo; y <= yhi; y++ ) {
					y2 = y - c[1];
					y2 *= y2;
					for ( x = xlo; x <= xhi; x++ ) {
						x2 = x - c[0];
						x2 *= x2;
						i = p->index(0, x, y, z, n);
						p->add(i, exp(nhis2*(x2 + y2 + z2)));
					}
				}
			}
		}
		n++;
	}

	return p;
}

/**
@brief 	Concatenates all model maps into one multi-image file.  
@param 	*model			model.
@param 	&filename		new map file name.
@return int				0.

	The model map file name must be set and point to a valid file.

**/
int			model_catenate_maps(Bmodel* model, Bstring& filename)
{
	long			i;
	Bmodel*			mp;
	Bstring*		file_list = NULL;
	Vector3<long>	nusize;
	Bstring			rawstring;
	
	for ( i=0, mp = model; mp; mp = mp->next, i++ ) {
		string_add(&file_list, mp->mapfile().c_str());
		mp->mapfile() = filename.str();
		mp->image_number(i);
	}

	Bimage*		pc = img_catenate(file_list, rawstring, Unknown_Type, nusize, 0, 0, 0, 0, 0);
	
	write_img(filename, pc, 0);
	
	delete pc;
	
	return 0;
}

/**
@brief 	Fits a shell model as a rigid body to a map.  
@param 	*model			model.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	neg				flag to set contrast negative.
@return int				0.

	The model map file name must be set and point to a valid file.

**/
int			model_shell_fit(Bmodel* model, double hires, double lores, int neg)
{
	if ( lores < 1 ) lores = 1000;
	if ( hires > lores ) swap(hires, lores);
	
	Bimage*			p = NULL;
	 
	if ( ( p = read_img(model->mapfile(), 1, model->image_number()) ) == NULL )
		error_show("Error in model_shell_fit: No map read!", __FILE__, __LINE__);
	
	double			sigma = hires;
	Vector3<double>	origin(p->image->origin());
	Vector3<long>	size = {p->sizeX(), p->sizeY(), p->sizeZ()};
	Vector3<double>	sam(p->sampling(0));
	
	Bimage*			pmod = img_from_model(model, origin, size, sam, sigma);

	if ( neg ) pmod->invert();
	
	if ( verbose )
		cout << "Fitting the model to the map" << endl;
	
	double			cc(0), radius = p->sizeX()/4;

	Vector3<double>	shift = p->find_shift(pmod, NULL, hires, lores, radius, 0, 1, cc);
	model->FOM(cc);

	if ( verbose ) {
		cout << "Fit CC:                         " << cc << endl;
		cout << "Shifting the model by:          " << shift << endl << endl;
	}
	
	model_shift(model, shift*pmod->sampling(0));
	
	delete p;
	delete pmod;
	
	return 0;
}

/**
@brief 	Calculates a radial profile from a shell model.  
@param 	*model			model.
@return int				0.

	The model map file name must be set and point to a valid file.

**/
int			model_shell_radial_profile(Bmodel* model)
{
	if ( model->mapfile().length() < 1 )
		return error_show("Error in model_shell_radial_profile: No map specified!", __FILE__, __LINE__);
	
	Bimage*			p = NULL;
	 
	if ( ( p = read_img(model->mapfile(), 1, model->image_number()) ) == NULL )
		error_show("Error in model_shell_radial_profile: No map read!", __FILE__, __LINE__);
	
	p->change_type(Float);
	
	int				ncomp, ntype, n;
	long	i, j;
	double			radstep = p->sampling(0)[0];
	int				maxrad = p->sizeX()/2;
	double			d, davg, v;
	Vector3<double>	nvec, vec, origin(p->image->origin());
	
	Bcomponent*		comp;
	Bcomptype*		ct = NULL;
	
	for ( ntype=0, ct = model->type; ct; ct = ct->next ) ntype++;
	
	if ( ntype < 1 ) ntype = 1;
	
	float*			r = new float[ntype*maxrad];
	float*			r2 = new float[ntype*maxrad];
	float*			w = new float[ntype*maxrad];
	for ( i=0; i<ntype*maxrad; i++ ) r[i] = r2[i] = w[i] = 0;
	
	for ( ncomp=0, davg=0, comp = model->comp; comp; comp = comp->next, ncomp++ )
		davg += comp->location().length()/radstep;

	if ( ncomp < 1 )
		error_show("Error in model_shell_radial_profile: No components!", __FILE__, __LINE__);
	
	davg /= ncomp;

	if ( verbose ) {
		cout << "Calculating a radial profile based on a shell model:" << endl;
		cout << "Radial step size:               " << radstep << endl;
		cout << "Average component distance:     " << davg << endl << endl;
	}
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		ct = comp->type();
		if ( !ct ) n = 0;
		nvec = comp->location()/radstep;
		d = nvec.normalize() - davg;
		for ( j=0; j<maxrad; j++ ) {
			i = n*maxrad + j;
			vec = nvec*(d+j) + origin;
			if ( vec[0]>=0 && vec[0]<=p->sizeX()-2 && vec[1]>=0 && vec[1]<=p->sizeY()-2 && vec[2]>=0 && vec[2]<=p->sizeZ()-2 ) {
//				v = p->interpolate(vec[0], vec[1], vec[2], 0, p->image->background());
				v = p->interpolate(vec, 0, p->background(long(0)));
				r[i] += v;
				r2[i] += v*v;
				w[i] += 1;
			}
		}
	}

	for ( i=0; i<ntype*maxrad; i++ ) if ( w[i] ) {
		r[i] /= w[i];
		r2[i] = r2[i]/w[i] - r[i]*r[i];
		if ( r2[i] > 0 ) r2[i] = sqrt(r2[i]);
		else r2[i] = 0;
	}
	if ( verbose ) {
		for ( n=0, ct = model->type; ct && n<ntype; ct = ct->next, n++ ) {
			if ( ct ) cout << "Component type:                 " << ct->identifier() << endl;
			cout << "#\tRadius\tDensity\tSTDev\tWeight" << endl;
			for ( j=0, i=n*maxrad; j<maxrad; j++, i++ )
				cout << j << tab << j*radstep << tab << r[i] << tab << r2[i] << tab << w[i] << endl;
			cout << endl;
		}
	}
	
	delete[] r;
	delete[] r2;
	delete[] w;
	delete p;
	
	return 0;
}



/**
@brief 	Calculates a radial profile from a shell model.  
@param 	*model		model.
@param 	size		size of component density to extract.
@param 	origin		origin for extracted densities.
@param 	ft_size		Fourier transform size.
@param 	ann_min		minimum annulus for rotational alignment.
@param 	ann_max		maximum annulus for rotational alignment.
@param 	hires		high resolution limit for cross-correlation.
@param 	lores		low resolution limit for cross-correlation.
@return Bimage*		0.

	The model map file name must be set and point to a valid file.

**/
Bimage*		model_shell_power_spectrum(Bmodel* model, Vector3<long> size,
				Vector3<double> origin, int ft_size, int ann_min, int ann_max,
				double hires, double lores)
{
	if ( ft_size < size[0] ) ft_size = findNextPowerOf(size[0], 2);
	if ( ft_size < 128 ) ft_size = 128;

//	Vector3<double>	origin(size[0]/2, size[1]/2, size[2]/2);

	Bimage*			p = model_extract_component_densities(model, size, origin);
	
	if ( ann_min < 1 ) ann_min = 1;
	if ( ann_max < 1 ) ann_max = p->sizeX()/2;
	if ( ann_max < ann_min ) swap(ann_min, ann_max);
	
	if ( lores < 1 ) lores = p->sampling(0)[0]*p->sizeX()/2;
	if ( hires < 1 ) hires = 2*p->sampling(0)[0];
	if ( lores < hires ) swap(hires, lores);
	
	p->project('z', 1);
	size[2] = 1;
	origin[2] = 0;
	
	Vector3<double>	start(1,1,0);
	Vector3<int>	esize(size[0] - 2, size[1] - 2, 1);

	p->edge(1, esize, start, 1, FILL_BACKGROUND, 0);

	int				i, j, k, n, jb, nmax(6);

	int 			nannuli = p->sizeX()/2;			// Sampling same as image
	int 			nangles(NPOLANG); 				// 0.5 degree step size
	double			cc1, ccbest(0);					// Best correlation coefficient
	double			pcc(0);							// Polar power spectrum correlation coefficient
	Vector3<double>	scale(1,1,1), translate, shift1;
	Vector3<double>	shift;
	Vector3<double>	axis(0,0,1);
	Matrix3			mat;
	View			view;
	Bcomponent*		comp;
	Bimage*			pone;
	Bimage*			pol;
	Bimage*			polref;
	Bimage*			pt;
	
	Bcomptype*		ct[6];
	Bstring			ctstr[6];
	ctstr[0] = "MON";
	ctstr[1] = "DI";
	ctstr[2] = "TRI";
	ctstr[3] = "TET";
	ctstr[4] = "PEN";
	ctstr[5] = "HEX";
	
	comp_type_list_kill(model->type);
	model->type = NULL;

	Bimage**		pref = new Bimage*[nmax];
	
	for ( j=0; j<nmax; j++ ) {
		pref[j] = new Bimage(Float, TSimple, size, 1);
		pref[j]->origin(origin);
		pref[j]->sampling(p->sampling(0));
		ct[j] = model->add_type(ctstr[j]);
	}

	int*			num = new int[nmax];
	float*			cc = new float[nmax];
	float*			ccs = new float[nmax];
	float*			ccr = new float[nmax];
	Vector3<float>*	ori = new Vector3<float>[nmax];
	for ( i=0; i<nmax; i++ ) cc[i] = ccs[i] = ccr[i] = num[i] = 0;
	
	cc[0] = 1;

    // Specify FFTW arrays
	fft_plan		planf = fft_setup_plan(nangles, 1, 1, FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(nangles, 1, 1, FFTW_BACKWARD, 1);
	
	if ( verbose )
		cout << "Comp\tOrder\tAngle\tCC" << endl;
	for ( i=0, comp = model->comp; i<p->images() && comp; i++, comp = comp->next ) {
		pone = p->extract(i);
		for ( j=1; j<nmax; j++ ) {
			n = j + 1;
			cc[j] = ccs[j] = 0;
			ori[j] = 0;
			for ( k=1; k<n; k++ ) {
				mat = Matrix3(axis, TWOPI*k*1.0L/n);
				pt = pone->transform(size, scale, origin, translate, mat, FILL_BACKGROUND, 0);
				shift = pone->find_shift(pt, NULL, hires, lores, pone->sizeX()/4, 0, 1, cc1);
				delete pt;
				cc[j] += cc1;
				ccs[j] += cc1*cc1;
				mat = Matrix3(1) - mat;
				mat = mat.singular_value_decomposition();
				ori[j] += mat * shift + origin;
			}
			cc[j] /= j;
			ccs[j] = ccs[j]/j - cc[j]*cc[j];
			if ( ccs[j] > 0 ) ccs[j] = sqrt(ccs[j]);
			else ccs[j] = 0;
			cc[j] += ccs[j];
			ori[j] /= j;
		}
		ccr[1] = 4*cc[1]/(cc[2]+cc[3]+cc[4]+cc[5]);
		ccr[2] = 4*cc[2]/(cc[1]+cc[3]+cc[4]+cc[5]);
		ccr[3] = 3*(cc[1]+cc[3])/(2*(cc[2]+cc[4]+cc[5]));
		ccr[4] = 4*cc[4]/(cc[1]+cc[2]+cc[3]+cc[5]);
		ccr[5] = 2*(cc[1]+cc[2]+cc[5])/(3*(cc[3]+cc[4]));
		for ( j=1, jb=0, ccbest = 1.1; j<nmax; j++ ) {
			if ( ccbest < ccr[j] ) {
				ccbest = ccr[j];
				jb = j; 
			}
			cout << cc[j] << tab << ccr[j] << tab;
		}
		if ( !pref[jb] ) {
			view[3] = 0;
		} else {
//			pone->image->origin(ori[jb]);
			pone->origin(ori[jb]);
//			polref = img_polar_transform(pref[jb], nannuli, nangles, ann_min, ann_max);
//			pol = img_polar_transform(pone, nannuli, nangles, ann_min, ann_max);
			polref = pref[jb]->cartesian_to_cylindrical(nannuli, nangles, 1);
			pol = pone->cartesian_to_cylindrical(nannuli, nangles, 1);
//			view[3] = comp->view()[3] = img_correlate_annuli(polref, pol, 1, ann_min, ann_max, planf, planb, pcc);
			view[3] = comp->view()[3] = polref->correlate_annuli(pol, ann_min, ann_max,
				0, TWOPI, planf, planb, pcc);
			delete pol;
			delete polref;
		}
		pref[jb]->rotate_and_add(pone, ori[jb], view);
		delete pone;
		comp->FOM(cc[jb]);
		comp->type(ct[jb]);
		switch ( jb+1 ) {
			case 2: comp->color(RGBA<float>(1,1,0,1)); break;
			case 3: comp->color(RGBA<float>(1,0.5,0,1)); break;
			case 4: comp->color(RGBA<float>(0.3,1,0.3,1)); break;
			case 5: comp->color(RGBA<float>(0,0,1,1)); break;
			case 6: comp->color(RGBA<float>(1,0,0,1)); break;
			default: comp->color(RGBA<float>(1,1,1,1));
		}
		num[jb]++;
		if ( verbose )
			cout << comp->identifier() << tab << jb+1 << tab << comp->view().angle()*180.0/M_PI << tab << pcc << endl;
	}
	
    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	delete p;

	if ( verbose ) {
		cout << "Order\tCount" << endl;
		for ( j=0; j<nmax; j++ )
			cout << j+1 << tab << num[j] << endl;
	}
	
	for ( j=0; j<nmax; j++ ) if ( num[j] ) pref[j]->multiply(1.0L/num[j]);

	delete[] num;
	delete[] cc;
	delete[] ccs;
	delete[] ccr;
	delete[] ori;

	Bimage*			pfin = new Bimage();
	pfin->catenate(nmax, pref);

	for ( j=0; j<nmax; j++ ) delete pref[j];
	delete[] pref;

	pfin->pad(ft_size, FILL_BACKGROUND, pfin->background(long(0)));

	pfin->fft();

	pfin->complex_to_intensities();
	
	for ( i=0; i<pfin->images(); i++ )
		pfin->image[i].origin(0,0,0);

	pfin->center_wrap();

	return pfin;
}

int			comp_set_fom_sym(Bcomponent* comp, Bimage* p, long minorder, long maxorder)
{
	if ( minorder < 1 ) minorder = 1;
	
	long		i, j, k, is(0), is2, xx, nyz(p->sizeY()*p->sizeZ()), nf, nb;
	double		v, sf, sb, mx(0), mx2(0);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG comp_set_fom_sym: comp_id=" << comp->identifier() << endl;
	
	Bimage*		psym = new Bimage(Float, TSimple, maxorder+1, p->sizeY(), p->sizeZ(), 1);
	
	for ( j=0; j<nyz; j++ ) {
		for ( k=minorder; k<=maxorder; k++ ) {
			i = j*p->sizeX();
			nf = nb = 0;
			sf = sb = 0;
			for ( xx=minorder, i+=xx; xx<p->sizeX()/2; xx++, i++ ) {
				v = (*p)[i] * sqrt(1.0*xx);
				if ( k%xx == 0 ) {
					sf += v;
					nf++;
				} else {
					sb += v;
					nb++;
				}
			}
			sf /= nf;
			sb /= nb;
			if ( sb < 0.001 ) sb = 0.001;
			if ( sf < 0.001 ) sf = 0.001;
//			psym->set(j*psym->sizeX()+k, 1 - sb/sf);
//			psym->set(j*psym->sizeX()+k, sf - sb);
			psym->set(j*psym->sizeX()+k, sf/sb - 1);
//			cout << k << tab << nf << "," << sf << tab << nb << "," << sb << endl;
		}
		mx = mx2 = 0;
		is = is2 = 0;
		for ( k=minorder; k<=maxorder; k++ ) {	// find the highest value
			v = (*psym)[j*psym->sizeX()+k];
			if ( mx < v ) {
				mx = v;
				is = k;
			}
		}
		for ( k=minorder; k<=maxorder; k++ ) if ( k != is ) {	// find the second highest value
			v = (*psym)[j*psym->sizeX()+k];
			if ( mx2 < v ) {
				mx2 = v;
				is2 = k;
			}
		}
		if ( is2%is == 0 ) {	// if the highest is a factor of the second, find another second highest
			mx2 = 0;
			for ( k=minorder; k<=maxorder; k++ ) if ( k != is && k != is2 ) {
				v = (*psym)[j*psym->sizeX()+k];
				if ( mx2 < v ) {
					mx2 = v;
				}
			}
			is = is2;
		}
	}
	
	delete psym;
	
	comp->FOM(1 - mx2/mx);
	comp->select(is);

	return 0;
}

/**
@brief 	Determines the cyclic symmetry of components.
@param 	*model			model parameters.
@param 	nangles			number of angles.
@param 	ann_min			minimum annulus.
@param 	ann_max			maximum annulus.
@param 	ann_width		annular width.
@param 	zmin			minimum z limit.
@param 	zmax			maximum z limit.
@param 	zinc			z increment.
@param	minorder		minimum cyclic order to search for.
@param	maxorder		maximum cyclic order to search for.
@return int				0.
**/
int			model_component_symmetry(Bmodel* model, long nangles,
				long ann_min, long ann_max, long ann_width, 
				long zmin, long zmax, long zinc, long minorder, long maxorder)
{
	if ( ann_width < 1 ) ann_width = ann_max - ann_min;
	
	long			i, nsym[maxorder+1];
	double			fomsym[maxorder+1], sep(0), sep_avg(0);
	Vector3<long>	size;
	Vector3<double>	origin;
	Matrix3			mat(1);
	Bstring			id;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bimage*			p = NULL;
	Bimage*			pex = NULL;
	Bimage*			pps = NULL;
	
	for ( i=0; i<=maxorder; i++ ) fomsym[i] = nsym[i] = 0;

	fft_plan		plan = fft_setup_plan(nangles, 1, 1, FFTW_FORWARD, 1);
	
	p = read_img(model->mapfile(), 0, 0);

	if ( zmax <= 0 ) {
		zmin = p->sizeZ() / 4;
		zmax = p->sizeZ() - zmin - 1;
		zinc = zmax - zmin + 1;
	}
	
	delete p;
	
	if ( verbose ) {
		cout << "Determining cyclic symmetry:" << endl;
		cout << "Maximum order:                " << nangles << endl;
		cout << "Number of angles:             " << maxorder << endl;
		cout << "Annuli:                       " << ann_min << " - " << ann_max << " @ " << ann_width << " pixels" << endl;
		cout << "Z limits:                     " << zmin << " - " << zmax << " @ " << zinc << " pixels" << endl;
	}

	if ( verbose )
		cout << "Model\tComponent\tOrder\tPower\tSeparation" << endl;
	for ( i=0, mp = model; mp; mp = mp->next ) {
		p = read_img(model->mapfile(), 1, mp->image_number());
		for ( comp = mp->comp; comp; comp = comp->next, i++ ) {
			size = Vector3<long>(4*comp->radius(), 4*comp->radius(), 4*comp->radius());
			size = size.min(p->size());
			origin = size/2;
			pex = p->extract(0, p->image->origin() + comp->location()/p->sampling(0), size, origin, mat);
//			pex->rotate(pex->size()/2, comp->view());
			if ( ann_max > 0 ) {
				if ( ann_width < 1 ) ann_width = ann_max - ann_min;
				if ( ann_width < 1 ) ann_width = 1;
			} else if ( ann_width > 0 ) {
				ann_min = comp->radius()/p->sampling(0)[0] - ann_width/2;
				ann_max = ann_min + ann_width - 1;
			} else {
				ann_min = comp->radius()/(2*p->sampling(0)[0]);
				ann_max = 3*ann_min;
				ann_width = ann_max - ann_min + 1;
			}
			pps = pex->polar_transform(nangles, ann_min, ann_max, ann_width, zmin, zmax, zinc);
			delete pex;
			pps->line_powerspectra(plan);
			comp_set_fom_sym(comp, pps, minorder, maxorder);
			delete pps;
			sep = TWOPI*comp->radius()/comp->select();
			sep_avg += sep;
			if ( verbose )
				cout << mp->identifier() << tab << comp->identifier() << tab << comp->select() << tab <<
					setprecision(4) << comp->FOM() << tab << sep << endl;
			nsym[comp->select()]++;
			fomsym[comp->select()] += comp->FOM();
		}
		delete p;
	}
	
	if ( i ) sep_avg /= i;
	
    fft_destroy_plan(plan);

	for ( i=minorder; i<=maxorder; i++ ) if ( nsym[i]) {
		fomsym[i] /= nsym[i];
	}
	
	if ( verbose ) {
		cout << "Order\tCount\tFOM" << endl;
		for ( i=minorder; i<=maxorder; i++ )
			cout << i << tab << nsym[i] << tab << fomsym[i] << endl;
		cout << "Separation average:             " << sep_avg << endl << endl;
	}

	return 0;
}

/*
	The new structure factors are added to the input image.
*/
int			scatter_one(Bcomponent* comp, Bimage* p,
				vector<double>& scurve, CTFparam cp, double ds, double scut, int flag)
{
	bool			ab(flag&1), ew(flag&2);
	long			i, t, xx, yy;
	double			ids(1/ds), s, s2, sc2(scut*scut), sx, sy, sx2, sy2;
	double			f1, f, phi, pilz(0);
	Vector3<long>	h((p->size()+1)/2);
	
//	if ( ab ) cp.defocus_average(cp.defocus_average() + comp->location()[2]);

	if ( ew ) pilz = M_PI*cp.lambda()*comp->location()[2];	// Ewald sphere phase shift prefactor

	for ( i=yy=0; yy<p->sizeY(); ++yy ) {
		sy = ( yy < h[1] )? yy: yy - p->sizeY();
		sy /= p->real_size()[1];
		sy2 = sy*sy;
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			sx = ( xx < h[0] )? xx: xx - p->sizeX();
			sx /= p->real_size()[0];
			sx2 = sx*sx;
			s2 = sx2 + sy2;
			if ( s2 <= sc2 ) {
				s = ids*sqrt(s2);		// Sampling relative to the scattering curve
				t = long(s);
				f1 = s - t;				// Fraction for interpolation
				f = comp->density()*((1-f1)*scurve[t] + f1*scurve[t+1]);		// Amplitude
				phi = MIN2PI*(sx*comp->location()[0] + sy*comp->location()[1]);	// Phase shift
				if ( pilz ) phi -= pilz*s2;					// Ewald sphere phase shift
//				if ( ab ) f *= cp.calculate(s2, atan2(sy,sx));		// Aberrations (CTF)
				if ( ab ) phi += M_PI_2 + cp.delta_phi(s2, atan2(sy,sx));		// Aberrations (CTF)
//				p->add(i, Complex<float>(f*cosl(phi), f*sinl(phi)));
				p->add(i, complex_polar(f, phi));
			}
		}
	}

	return 1;
}

/*
	The input image is not modified
*/
Bimage*		scatter_one_img(Bcomponent* comp, Bimage* p,
				vector<double>& scurve, CTFparam cp, double ds, double scut, int flag)
{
	Bimage*			pone = new Bimage(Float, TComplex, p->size(), 1);

	scatter_one(comp, pone, scurve, cp, ds, scut, flag);

	return pone;
}


Bimage*		scatter_to_img(const vector<Bcomponent*>& carr, long ib, long ie, Bimage* p,
				map<string, vector<double>>& scat, CTFparam& cp, double ds, double scut, int flag)
{
	long			ic;
	Bcomponent*		comp;
	string			cel;

	Bimage*			pone = new Bimage(Float, TComplex, p->size(), 1);

	for ( ic=ib; ic<ie; ++ic ) {
		comp = carr[ic];
		
		cel = component_element(comp);
		
		if ( scat.find(cel) == scat.end() ) {
			cerr << "Warning: Scattering curve for element " << cel << " not found!" << endl;
		} else {
			vector<double>&	scurve = scat.at(cel);
			
			scatter_one(comp, pone, scurve, cp, ds, scut, flag);
		}
	}

	return pone;
}

Bimage*		scatter_slice_to_img(const vector<Bcomponent*>& carr, Bimage* p,
				map<string, vector<double>>& scat, CTFparam& cp, double ds, double scut, int flag)
{
	Bimage*			pone = new Bimage(Float, TComplex, p->size(), 1);

	if ( verbose & VERB_DEBUG ) {
		cout << "Scattering curves: " << scat.size() << endl;
		for ( auto it = scat.begin(); it != scat.end(); ++it )
			cout << it->first << endl;
	}

	for ( auto comp: carr ) {
		string		cel = component_element(comp);
		
		if ( scat.find(cel) == scat.end() ) {
			cerr << "Warning: Scattering curve for element " << cel << " not found!" << endl;
		} else {
			vector<double>&	scurve = scat.at(cel);
			
			scatter_one(comp, pone, scurve, cp, ds, scut, flag);
		}
	}

	return pone;
}



int			img_add(Bimage* ps, vector<Complex<float>>& v)
{
	for ( long i=0; i<ps->data_size(); ++i )
		ps->add(i, v[i]);

	return 0;
}

int			img_add_fast_old(Bimage* ps, Bimage* p)
{
	for ( long i=0; i<ps->data_size(); ++i )
		ps->add(i, (*p)[i]);

	return 0;
}

int			img_add_fast(Bimage* ps, Bimage* p)
{
	float*		fds = (float *) ps->data_pointer();
	float*		fd = (float *) p->data_pointer();
	
	for ( long i=0; i<ps->data_size(); ++i, ++fd, ++fds )
		*fds += *fd;

	return 0;
}

int			img_electron_scattering(Bmodel* model, Bimage* p,
				CTFparam& cp, double dose, double stdev, map<string, vector<double>>& scat, double ds, double scut, int flag)
{
	if ( dose && stdev ) model_random_displace_number(model, dose, stdev);

	vector<Bcomponent*>	carr = models_get_component_array(model);
	long				ncomp(carr.size());

#ifdef HAVE_GCD
	__block long		nd(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(ncomp, dispatch_get_global_queue(0, 0), ^(size_t i){
		string				cel = component_element(carr[i]);
		if ( scat.find(cel) == scat.end() ) {
			cerr << "Warning: Scattering curve for element " << cel << " not found!" << endl;
		} else {
			vector<double>	scurve = scat.at(cel);
//			__block vector<Complex<float>>	v = scatter_one(carr[i], p, scurve, cp, ds, scut, flag);
			Bimage*			pone = scatter_one_img(carr[i], p, scurve, cp, ds, scut, flag);
			dispatch_sync(myq, ^{
//				img_add(p, v);
				img_add_fast(p, pone);
				nd++;
				cout << " " << nd << "/" << ncomp << "\r" << flush;
			});
		}
	});
#else
	long			nd(0);
#pragma omp parallel for
	for ( long i=0; i<ncomp; ++i ) {
		string				cel = component_element(carr[i]);
		if ( scat.find(cel) == scat.end() ) {
			cerr << "Warning: Scattering curve for element " << cel << " not found!" << endl;
		} else {
			vector<double>	scurve = scat.at(cel);
//			vector<Complex<float>>	v = scatter_one(carr[i], p, scurve, cp, ds, scut, flag);
			Bimage*			pone = scatter_one_img(carr[i], p, scurve, cp, ds, scut, flag);
#pragma omp critical
			{
//				img_add(p, v);
				img_add_fast(p, pone);
				nd++;
				cout << " " << nd << "/" << ncomp << "\r" << flush;
			}
			delete pone;
		}
	}
#endif

	return 0;
}

int			img_electron_scattering_chunks(Bmodel* model, Bimage* p,
				CTFparam& cp, double dose, double stdev, map<string, vector<double>>& scat, double ds, double scut, int flag)
{
	if ( dose && stdev ) model_random_displace_number(model, dose, stdev);

	vector<Bcomponent*>	carr = models_get_component_array(model);
	long				ncomp(carr.size());
	long				nc(100), n(ncomp/nc+1);
	
	cout << ncomp << tab << nc << tab << n << endl;

#ifdef HAVE_GCD
	__block long		nd(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(nc, dispatch_get_global_queue(0, 0), ^(size_t i){
		long			j = ((i+1)*n<ncomp)? (i+1)*n: ncomp;
		Bimage*			pone = scatter_to_img(carr, i*n, j, p, scat, cp, ds, scut, flag);
		dispatch_sync(myq, ^{
			img_add_fast(p, pone);
			nd += j-i*n;
			cout << " " << nd << "/" << ncomp << "\r" << flush;
		});
	});
#else
	long				nd(0);
#pragma omp parallel for
	for ( long i=0; i<ncomp; i+=n ) {
		long			j = (i+n<ncomp)? i+n: ncomp;
		Bimage*			pone = scatter_to_img(carr, i, j, p, scat, cp, ds, scut, flag);
#pragma omp critical
		{
			img_add_fast(p, pone);
			nd += j-i;
			cout << " " << nd << "/" << ncomp << "\r" << flush;
		}
		delete pone;
	}
#endif

	return 0;
}

/*
	Takes the input number of sub-images as the number of slices
*/
int			img_electron_scattering_slices(Bmodel* model, Bimage* p, CTFparam& cp, 
				map<string, vector<double>>& scat, double ds, double scut, int flag)
{
	vector<Vector3<double>>	bounds = models_calculate_bounds(model);
	double				thickness(p->image->sampling()[2]);
	double				bottom(bounds[0][2]);
	double				top(p->images()*thickness + bottom);
	vector<vector<Bcomponent*>>	comp_slice = model_split_into_slices(model, bottom, top, thickness);

#ifdef HAVE_GCD
	__block long		nd(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(p->images(), dispatch_get_global_queue(0, 0), ^(size_t i){
		Bimage*			pone = scatter_slice_to_img(comp_slice[i], p, scat, cp, ds, scut, flag);
		dispatch_sync(myq, ^{
			p->replace(i, pone);
			nd++;
			cout << " " << nd << "/" << p->images() << "\r" << flush;
		});
		delete pone;
	});
#else
	long				nd(0);
#pragma omp parallel for
	for ( long i=0; i<p->images(); ++i ) {
		Bimage*			pone = scatter_slice_to_img(comp_slice[i], p, scat, cp, ds, scut, flag);
#pragma omp critical
		{
			p->replace(i, pone);
			nd++;
			cout << " " << nd << "/" << p->images() << "\r" << flush;
		}
		delete pone;
	}
#endif

	return 0;
}



int			img_electron_scattering(Bmodel* model, Bimage* p,
				CTFparam& cp, double dose, double stdev, Bstring& atompropfile, int flag)
{
	if ( dose ) (*p)["dose"] = dose;	// Dose per frame

	double			scut = cp.frequency_cutoff();
	
//	if ( scut > 0.5/p->sampling(0)[0] ) scut = 0.5/p->sampling(0)[0];

	double			ds(0.01), smax(5);
	
	map<string,Bcomptype>	atompar = read_atom_properties(atompropfile);
	
	JSvalue			el = model_elements(model, atompar);
	
	if ( verbose ) {
		cout << "Calculating structure factors:" << endl;
		cout << "Atomic properties file:         " << atompropfile << endl;
		cout << "Size:                           " << p->size() << tab << p->images() << endl;
		cout << "Dose:                           " << dose << " e/frame (" << stdev << ")" << endl;
		cout << "Frequency cutoff:               " << scut << " (" << 1/scut << " A)" << endl;
		cout << "Aberration flag:                " << (flag&1) << endl;
		cout << "Ewald sphere flag:              " << (flag&2) << endl;
		cp.show();
		cout << endl;
	}
	
	map<string, vector<double>>	scat = calculate_scattering_curves(el, atompar, ds, smax);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "Scattering curves: " << scat.size() << endl;
		for ( auto it = scat.begin(); it != scat.end(); ++it )
			cout << it->first << endl;
	}
	
	Bimage*			p1;
	double			lambda = cp.lambda();
	double			b2 = beta2(cp.volt());
	double			scale = lambda/sqrt(1 - b2);
	scale *= 1/p->real_size().volume();
	
	if ( dose ) {
		for ( long nn=0; nn<p->images(); ++nn ) {
			if ( verbose )
				cout << "Calculating image " << nn+1 << endl;
			p1 = p->extract(nn);
			img_electron_scattering(model, p1, cp, dose, stdev, scat, ds, scut, flag);
			p->replace(nn, p1);
		}
		if ( verbose )
			cout << endl;
		p->multiply(scale);
	} else {
		img_electron_scattering_slices(model, p, cp, scat, ds, scut, flag);
	}

	return 0;
}

