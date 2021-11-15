/**
@file	mg_helix.cpp
@brief	Functions to process helical data
@author	Bernard Heymann
@date	Created: 20061110
@date	Modified: 20170623
**/

#include "Bimage.h"
#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_helix.h"
#include "mg_extract.h"
#include "spline.h"
#include "linked_list.h"
#include "Complex.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates a filament profile and node profile.
@param 	*fnode		node list.
@param 	*p			image with filament.
@param	img_num		sub-image number.
@param 	id			identifier of selected node.
@param 	width		width of box around filament.
@param 	&n			number of elements in the profile.
@return double*		two profile arrays, each length n.

	A filament image is extracted from the micrograph along a spline through
	the nodes. The overall profile and the profile around a specified node
	are calculated and returned in a single array.
	It only works for 2D images.

**/
double*		filament_profile(Bfilnode* fnode, Bimage* p, long img_num, int id, double width, long& n)
{
	long				nspline(0);

	Vector3<double>*	spline = vector3_spline_from_nodes(fnode, nspline);

	Vector3<long>		size((long)width, nspline, 1);

	n = size[0];

	long				i, x, y, z;
	double				len(0);
	Vector3<double>		loc(fnode->loc);
	Bfilnode* 			fn;

	for ( fn = fnode; fn && fn->id != id; fn = fn->next ) {
		if ( fn != fnode ) len += (fn->loc - loc).length();
		loc = fn->loc;
	}
	if ( fn ) len += (fn->loc - loc).length();
	else {
		cerr << "Error: Filament node with ID " << id << " not found!" << endl;
		return NULL;
	}
	
	long				ns = (long) (len - width);
	long				ne = (long) (len + width);
	if ( ns < 0 ) ns = 0;
	if ( ne > nspline ) ne = nspline;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG filament_profile: id=" << id << " fn->id=" << fn->id <<
			" fn->loc=" << fn->loc << endl;
		cout << "DEBUG filament_profile: len=" << len << " ns=" << ns <<
			" ne=" << ne << " nspline=" << nspline << endl;
	}

	Bimage*				pfil = p->extract_filament(img_num, width, 2, nspline, spline);

	delete[] spline;
	
	double*				profile = new double[2*n];
	
	for ( i=0; i<2*n; i++ ) profile[i] = 0;
	
	for ( i=z=0; z<pfil->sizeZ(); z++ ) {
		for ( y=0; y<pfil->sizeY(); y++ ) {
			for ( x=0; x<pfil->sizeX(); x++, i++ ) {
				profile[x] += (*pfil)[i];
				if ( y >= ns && y < ne )
					profile[x+n] += (*pfil)[i];
			}
		}
	}

	delete pfil;
	
	for ( i=0; i<n; i++ ) profile[i] /= nspline;
	for ( ne -= ns; i<2*n; i++ ) profile[i] /= ne;
	
	return profile;
}

double		nodes_center(Bfilnode* fnode, Bimage* p, long img_num, int width)
{
	long				nspline(0);

	Vector3<double>*	spline = vector3_spline_from_nodes(fnode, nspline);

	Vector3<long>		size(width, nspline, 1);

	Bimage*				pfil = p->extract_filament(img_num, width, 2, nspline, spline);

	delete[] spline;
	
	long		i, x, y, n(0);
	long				ns, ne;
	Complex<float>*		profile = new Complex<float>[pfil->sizeX()];
	double				len(0), max, pwr, off, cc(0);
	Vector3<double>		loc(fnode->loc), vec;
	Bfilnode* 			fn;
	Matrix3				mat(0,1,0,-1,0,0,0,0,1);

	fft_plan			planf = fft_setup_plan(pfil->sizeX(), 1, 1, FFTW_FORWARD, 1);
	fft_plan			planb = fft_setup_plan(pfil->sizeX(), 1, 1, FFTW_BACKWARD, 1);
	
	if ( verbose )
		cout << "Node\tOffset\tCC" << endl;
	for ( fn = fnode; fn; fn = fn->next ) {
		len += fn->loc.distance(loc);
		ns = (long) (len - width);
		ne = (long) (len + width);
		if ( ns < 0 ) ns = 0;
		if ( ne > nspline ) ne = nspline;
		for ( i=0; i<pfil->sizeX(); i++ ) profile[i] = 0;
		for ( y=ns; y<ne; y++ ) {
			i = y*pfil->sizeX();
			for ( x=0; x<pfil->sizeX(); x++, i++ ) {
				profile[x][0] += (*pfil)[i];
			}
		}
		fftw(planf, profile);
		for ( pwr=0, i=0; i<pfil->sizeX(); i++ ) {
			pwr += profile[i].power();
			profile[i] *= profile[i];
		}
		fftw(planb, profile);
		for ( max=-1000, off=0, i=0; i<pfil->sizeX(); i++ ) if ( max < profile[i].real() ) {
			max = profile[i].real();
			off = i;
		}
//		max /= sqrt(2*pwr);
		max /= pwr;
		cc += max;
		n++;
		if ( off > pfil->sizeX()/2 ) off -= (long)pfil->sizeX();
		off /= 2;
		if ( verbose )
//			cout << fn->id << tab << m - (long)pfil->sizeX()/2 << endl;
			cout << fn->id << tab << off << tab << max << endl;
		if ( fn->next ) vec = fn->next->loc;
		else vec = fn->loc;
		vec -= loc;
		vec = mat * vec;
		vec.normalize();
		fn->loc += vec * off;
		loc = fn->loc;
	}

	delete pfil;
	delete[] profile;
	fft_destroy_plan(planf);
	fft_destroy_plan(planb);
	
	if ( n ) cc /= n;
	
	return cc;
}

/**
@brief 	Centers filament nodes in a micrograph.
@param 	*fillist		list of filaments.
@param 	*p				image with filament.
@param	img_num			sub-image number.
@param 	filament_width	width of box around filament.
@return double			average correlation coefficient.

	The profile for each node is calculated and centered by cross-correlation
	with the mirrored profile.

**/
double		filaments_center(Bfilament* fillist, Bimage* p, long img_num, int filament_width)
{
	long				n(0);
	double				cc(0), avg(0);
	Bfilament* 			fil;
				
	for ( fil = fillist; fil; fil = fil->next ) {
		if ( verbose )
			cout << "Filament: " << fil->id << endl;
		cc = nodes_center(fil->node, p, img_num, filament_width);
		avg += cc;
		n++;
		if ( verbose )
			cout << "Average correlation coefficient: " << cc << endl;
	}
	
	if ( n ) avg /= n;

	return avg;
}

/**
@brief 	Centers filament nodes in all micrographs in a project.
@param 	*project	project object.
@param 	filament_width	width of box around filament.
@return double				average correlation coefficient.

	The profile for each node is calculated and centered by cross-correlation
	with the mirrored profile.

**/
double		project_center_filaments(Bproject* project, int filament_width)
{
	long				n(0);
	double				cc(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bimage*				p;
	
	if ( project->select < 1 ) {
		if ( verbose )
			cout << "Centering filaments from micrographs" << endl << endl;
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( verbose )
					cout << "Centering filaments from micrograph: " << mg->id << endl << endl;
				p = read_img(mg->fmg, 1, mg->img_num);
				if ( p ) {
					cc += filaments_center(mg->fil, p, 0, filament_width);
					n++;
					delete p;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( verbose )
				cout << "Centering filaments from reconstruction: " << rec->id << endl << endl;
			p = read_img(rec->frec, 1, 0);
			if ( p ) {
				cc += filaments_center(rec->fil, p, 0, filament_width);
				n++;
				delete p;
			}
		}
	}
	
	if ( n ) cc /= n;

	return cc;
}

/**
@brief 	Calculates and average power spectrum from filaments.
@param 	*project		project parameter structure.
@param 	pad					additional padding before transformation.
@param 	rotated				flag to indicate if particles are already rotated.
@param 	&path			path to write power spectra.
@return int						0.

	Particles derived from picked filaments are extracted and transformed
	to orient the helical axis along the x-axis. These images are then
	Fourier transformed and their power spectra averaged.
	One average per micrograph is calculated.
	Requirements: The filaments must be picked and converted to particle locations.

**/
int			project_filament_powerspectrum(Bproject* project, int pad, int rotated, Bstring& path)
{
	if ( !project || !project->field || !project->field->mg ) return -1;
	
	if ( pad < 0 ) pad = 0;
//	if ( pad > 1000 ) pad = 1000;
	
	if ( path.length() && path[-1] != '/') path += "/";
	
	long				i;
	Bfield*				field = project->field;
	Bmicrograph*		mg = field->mg;
	Bparticle*			part = part_find_first(project);
	Bstring				insert("_filps.");
	Bstring				ext("map");
	
	if ( !part ) {
		cerr << "Error in project_filament_powerspectrum: No particle!" << endl;
		return -1;
	}
	
	Bimage*				p = particle_read_img(part, 0);
	if ( !p ) {
		cerr << "Error in project_filament_powerspectrum: No image!" << endl;
		return -1;
	}
	
	Vector3<long>		size(p->size());
	if ( size[0] > 1 ) size[0] += pad;
	if ( size[1] > 1 ) size[1] += pad;
	if ( size[2] > 1 ) size[2] += pad;
	
	Bimage*				pps = new Bimage(Float, TSimple, size, 1);
	pps->sampling(p->sampling(0));
	
	double				angle;
	Vector3<double>		translate, origin;
	Vector3<double>		scale(1,1,1);
	Vector3<double>		axis(0,0,1);
	Matrix3				mat;
	
	delete p;
	
	if ( verbose ) {
		cout << "Calculating an average filament power spectrum" << endl;
		cout << "Padding:                        " << pad << endl;
		cout << "Power spectrum size:            " << size << endl;
		cout << "Powerspectrum file path:        " << path << endl << endl;
	}

	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select && mg->part ) {
			pps->clear();
			for ( i=0, part = mg->part; part; part = part->next ) if ( part->sel ) {
				p = particle_read_img(part, 1);
				if ( p ) {
					p->change_type(Float);
					if ( rotated < 1 ) {
//						origin = p->image->origin();
						angle = part->view.angle() - atan2(part->view[1], part->view[0]) - M_PI_2;
//						mat = Matrix3(axis, angle);
//						img_transform(p, p->size(), scale, origin, translate, mat, FILL_BACKGROUND, 0);
						p->rotate(angle);
					}
					if ( pad )
						p->pad(size, FILL_BACKGROUND, 0);
					p->power_spectrum(4);
					pps->add(p);
//					pps->image->origin(p->image->origin());
					pps->origin(p->image->origin());
					delete p;
					i++;
				}
			}
			if ( i ) {
				mg->fps = mg->fmg.pre_rev('.') + insert + ext;
				if ( path.length() ) {
					if ( mg->fps.contains("/") ) mg->fps = mg->fps.post_rev('/');
					mg->fps = path + mg->fps;
				}
				pps->rescale(1.0L/i, 0);
				write_img(mg->fps, pps, 0);
			}
		}
	}

	delete pps;
	
	return 0;
}

/**
@brief 	Extracts filament images from micrographs defined in a project and estimates their density per length.
@param 	*project		micrograph project.
@param 	filament_width	extracted filament width.
@return Bimage*			image with all densities.
**/
Bimage*		project_filament_density(Bproject* project, int filament_width)
{
	long				nfil(0), nspline, width(3*filament_width);
	Bfield*				field;
	Bmicrograph*		mg;
//	Breconstruction*	rec;
	Bfilament*			fil;
	Bimage*				p, *pfil, *pd;
	Vector3<double>*	spline;
	
	Bimage*				plist = NULL;
	Bimage*				pi = NULL;
	
	if ( project->select < 1 ) {
		if ( verbose )
			cout << "Extracting filaments from micrographs" << endl << endl;
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) if ( mg->fmg.length() && mg->fil ) {
				p = read_img(mg->fmg, 1, mg->img_num);
				if ( p ) {
					if ( verbose )
						cout << "Micrograph " << mg->id << endl;
					for ( fil = mg->fil; fil; fil=fil->next ) {
						if ( mg->filament_width < width )
							mg->filament_width = width;
						spline = vector3_spline_from_nodes(fil->node, nspline);
//						pfil = img_extract_filament(p, 0, mg->filament_width, 'y', nspline, spline);
						pfil = p->extract_filament(0, mg->filament_width, 2, nspline, spline);
						pd = pfil->filament_density(filament_width);
						vector<double>	gauss = pd->histogram_gauss_fit(100, 1);
						if ( verbose )
							cout << fil->id << tab << gauss[1] << tab << gauss[2] << endl;
						delete[] spline;
						delete pfil;
						nfil++;
						if ( pi ) { pi->next = pd; pi = pd; }
						else { plist = pi = pd; }
					}
					delete p;
				} else {
					cerr << "Error: Micrograph " << mg->fmg << " not read!" << endl;
				}
			}
		}
	}
	
	long				i, x, nv(0);
	for ( pi = plist; pi; pi = pi->next ) nv += pi->sizeX();
	
	pd = new Bimage(Float, TSimple, nv, 1, 1, 1);
	pd->sampling(plist->sampling(0));
//	cout << "all values = " << nv << endl;
	
	for ( i=0, pi = plist; pi; pi = pi->next )
		for ( x=0; x<pi->sizeX(); x++, i++ ) pd->set(i, (*pi)[x]);
//	cout << "transferred = " << i << endl;

	pd->statistics();
	
	for ( pi = plist; pi; ) {
		plist = pi->next;
		delete pi;
		pi = plist;
	}
	
	return pd;
}


/**
@brief 	Generates layer lines given the unit cell vectors.
@param 	*mg			micrograph.
@param 	rad_lim		layer line radial limit.
@return int			number of layer lines generated, <0 on error.

	The structure factor location is given by:
		x = uh + vk
	where u and v are the unit cell vectors,
	and h and k are the associated Miller indices.

**/
int			mg_generate_layer_lines(Bmicrograph* mg, int rad_lim)
{
	int					n, ind_lim;
	double				r;
	Blayerline*			line;
	
	for ( n=0, r=0, line = mg->layer; line; line = line->next )
		if ( line->number ) {
			r += line->distance/line->number;
			n++;
		}
	
	if ( n < 1 ) return -1;
	
	r /= n;		// distance between subsequent lines
	
	ind_lim = (int) (rad_lim/r);
	
	kill_list((char *) mg->layer, sizeof(Blayerline));
	mg->layer = line = NULL;

	for ( n=-ind_lim; n<=ind_lim; n++ ) {
		line = (Blayerline *) add_item((char **) &line, sizeof(Blayerline));
		if ( !mg->layer ) mg->layer = line;
		line->number = n;
		line->order = abs(n);
		line->distance = n*r;
		line->amp = 1;
		line->fom = 1;
		line->sel = 1;
	}
	
	return n;
}

/**
@brief 	Masks the image using the list of layer lines.
@param 	*p			complex image.
@param 	*layer_line	layer line list.
@param 	helix_axis	helix axis angle.
@param 	width		width of mask for a line.
@return int			number of lines or error code.
**/
int			img_mask_layer_lines(Bimage* p, Blayerline* layer_line, float helix_axis, float width)
{
	if ( p->compound_type() != TComplex ) return -1;
	
	p->change_type(Float);
	
	long				i, n, datasize = (long) p->size().volume();
	int					ind_lim(0);
	double				rad(0), hlen;
	Complex<double>		cv;
	Vector3<double>		start, end;
	Blayerline*			line;

	Bimage*				pmask = new Bimage(Float, TSimple, p->size(), 1);

	for ( line = layer_line; line; line = line->next ) {
		if ( ind_lim < abs(line->number) ) {
			ind_lim = abs(line->number);
			rad = fabs(1.1*line->distance);
		}
	}
	
	double				cosax = cos(helix_axis);
	double				sinax = sin(helix_axis);
	
	for ( n=0, line = layer_line; line; line = line->next, n++ ) {
		hlen = sqrt(rad*rad - line->distance*line->distance);
		start[0] = p->image->origin()[0] + line->distance*cosax - hlen*sinax;
		start[1] = p->image->origin()[1] + line->distance*sinax + hlen*cosax;
		end[0] = p->image->origin()[0] + line->distance*cosax + hlen*sinax;
		end[1] = p->image->origin()[1] + line->distance*sinax - hlen*cosax;
		pmask->bar(start, end, width, 0, FILL_USER, 1);
	}

	for ( i=0; i<datasize; i++ )
		if ( (*pmask)[i] < 1 ) p->set(i, cv);

//	write_img("m.map", pmask);
	
	delete pmask;
	
    return n;
}

/**
@brief 	Extracts one layer line from an image.
@param 	*p				power spectrum.
@param 	*line		layer line.
@param 	helix_axis		helix axis angle.
@param 	length				length of line.
@return double*					array with line values.
**/
double*		img_extract_layer_line(Bimage* p, Blayerline* line, float helix_axis, int length)
{
	long			i, hlen(length/2);
	double			cosax = cos(helix_axis);
	double			sinax = sin(helix_axis);
	Vector3<double>	c, v;
	double*			d = new double[length];
	
	c[0] = p->image->origin()[0] + line->distance*cosax - hlen*sinax;
	c[1] = p->image->origin()[1] + line->distance*sinax + hlen*cosax;
	v[0] = length*sinax;
	v[1] = -length*cosax;

	v.normalize();
	
	for ( i=0; i<length; i++, c+=v )
//		d[i] = p->interpolate(c[0], c[1], c[2], 0);
		d[i] = p->average(0, c[0], c[1], c[2], 0, 2);

	return d;
}

/**
@brief 	Extracts and prints layer lines from an image with corresponding Bessel functions.
@param 	*mg			micrograph.
@param 	length				length of line.
@param 	show				show: 1=extracted layer lines, 2=Bessel functions, 3=both
@return int						0.
**/
int			mg_extract_show_layer_lines(Bmicrograph* mg, int length, int show)
{
	if ( show < 1 || show > 3 ) return 0;
	
	Bimage*			p = read_img(mg->fps, 1, 0);
	if ( !p ) return -1;
	
	int				i, n, nll, k, hlen(length/2);
	Blayerline*		line;
	double			s, j, xscale = 1/(p->sizeX()*p->sampling(0)[0]), scale = TWOPI*mg->helix_radius;
	double*			d;

	for ( nll=0, line = mg->layer; line; line = line->next, nll++ ) ;
	
	double*			dm = new double[nll*length];

	for ( n=0, line = mg->layer; line; line = line->next, n++ ) {
		d = img_extract_layer_line(p, line, mg->helix_axis, length);
		for ( i=0, k=n; i<length; i++, k+=nll ) {
			s = (i - hlen)*xscale;
			dm[k] = d[i];
		}
		delete[] d;
	}

	cout << "Micrograph: " << mg->id << endl;
	if ( show == 1 ) cout << "Layer lines:" << endl;
	else if ( show == 2 ) cout << "Bessel functions associated with layer lines:" << endl;
	else cout << "Layer lines and associated Bessel functions:" << endl;
	
	cout << setprecision(4) << "s";
	for ( line = mg->layer; line; line = line->next ) {
		if ( show & 1 ) cout << "\ta(" << line->number << ")";
		if ( show & 2 ) cout << "\tj(" << line->order << ")";
	}
	cout << endl;
	
	for ( i=k=0; i<length; i++ ) {
		s = (i - hlen)*xscale;
		cout << s;
		for ( line = mg->layer; line; line = line->next, k++ ) {
			j = jn(line->order, s*scale);
			if ( show & 1 ) cout << tab << dm[k];
			if ( show & 2 ) cout << tab << j*j;
		}
		cout << endl;
	}
	
	delete[] dm;
	delete p;
	
	return n;
}


