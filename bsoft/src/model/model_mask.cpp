/**
@file	model_mask.cpp
@brief	Generating a mask from an atomic model
@author Bernard Heymann
@date	Created: 20060301
@date	Modified: 20200329
**/

#include "rwmodel.h"
#include "model_mask.h"
#include "model_shell.h"
#include "model_util.h"
#include "matrix_linear.h"
#include "Vector3.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Calculates a mask of the components and links of a model structure.
@param 	*model		model structure.
@param 	size		size of the new mask.
@param 	ori			origin of the new mask.
@param 	sam			pixel size of the new mask.
@param 	edge		edge width in angstrom.
@return Bimage*		new mask.

	Each component is used to generate a sphere and each link a bar.
	Only the first model in the linked list is used.

**/
Bimage*		model_create_mask(Bmodel* model, Vector3<long> size,
				Vector3<double> ori, Vector3<double> sam, double edge)
{
	if ( sam[0] <= 0) sam = sam.max(1);

	long			i, nv;
	double			width = edge/sam[0];
	Vector3<double>	loc, start, end;
	Bcomponent*		comp;
	Blink*			link;
	
	Bimage*			p = new Bimage(Float, TSimple, size, 1);
	if ( ori[0] == 0 && ori[1] == 0 && ori[2] == 0 )
		p->origin(p->default_origin());
	else
		p->origin(ori);
	ori = p->image->origin();
	p->sampling(sam);
	
	long			datasize = p->data_size();
	
	if ( verbose ) {
		cout << "Generating a model mask:" << endl;
		cout << "Size:                           " << p->size() << endl;
		cout << "Origin:                         " << p->image->origin() << endl;
		cout << "Sampling:                       " << p->sampling(0) << " A/voxel" << endl;
	}

	for ( comp = model->comp; comp; comp = comp->next ) {
		loc[0] = (long) (comp->location()[0]/p->sampling(0)[0] + p->image->origin()[0] + 0.5);
		loc[1] = (long) (comp->location()[1]/p->sampling(0)[1] + p->image->origin()[1] + 0.5);
		loc[2] = (long) (comp->location()[2]/p->sampling(0)[2] + p->image->origin()[2] + 0.5);
		p->sphere(loc, comp->radius()/p->sampling(0)[0], width, FILL_USER, 1);
	}
	
	for ( link = model->link; link; link = link->next ) {
		start[0] = (long) (link->comp[0]->location()[0]/p->sampling(0)[0] + p->image->origin()[0] + 0.5);
		start[1] = (long) (link->comp[0]->location()[1]/p->sampling(0)[1] + p->image->origin()[1] + 0.5);
		start[2] = (long) (link->comp[0]->location()[2]/p->sampling(0)[2] + p->image->origin()[2] + 0.5);
		end[0] = (long) (link->comp[1]->location()[0]/p->sampling(0)[0] + p->image->origin()[0] + 0.5);
		end[1] = (long) (link->comp[1]->location()[1]/p->sampling(0)[1] + p->image->origin()[1] + 0.5);
		end[2] = (long) (link->comp[1]->location()[2]/p->sampling(0)[2] + p->image->origin()[2] + 0.5);
		p->bar(start, end, 2*link->radius()/p->sampling(0)[0], width, FILL_USER, 1);
	}

	// Calculate the mask volume
	for ( i=nv=0; i<datasize; i++ ) if ( (*p)[i] > 0.5 ) nv++;
	
	if ( verbose )
		cout << "Mask volume:                    " << nv * sam.volume() << " A3" << endl << endl;
	
	return p;
}

/**
@brief 	Calculates a mask based on the periphery of a model structure.
@param 	*model		model structure.
@param 	size		size of the new mask.
@param 	ori			origin of the new mask.
@param 	sam			pixel size of the new mask.
@param 	curv_flag	flag to indicate curved surface.
@param 	fast		flag to use fast algorithm.
@return Bimage*		new mask.

	Each point in the new image is tested for inclusion in the mask,
	by calculating whether it falls inside the closest 3 vertices.
	Only the first model in the linked list is used.

**/
Bimage*		model_create_hull_mask(Bmodel* model, Vector3<long> size,
				Vector3<double> ori, Vector3<double> sam, int curv_flag, int fast)
{
	if ( sam[0] <= 0) sam = sam.max(1);

	long			i, nv, x, y, z;
	double			d, min = 1e30, max = 0;
	Vector3<double>	vec;
	Bcomponent*		comp;
	
	Bimage*			p = new Bimage(Float, TSimple, size, 1);
	if ( ori[0] == 0 && ori[1] == 0 && ori[2] == 0 )
		p->origin(p->default_origin());
	else
		p->origin(ori);
	ori = p->image->origin();
	p->sampling(sam);
	
	if ( verbose ) {
		cout << "Generating a hull mask:" << endl;
		cout << "Size:                           " << p->size() << endl;
		cout << "Origin:                         " << p->image->origin() << endl;
		cout << "Sampling:                       " << p->sampling(0) << " A/voxel" << endl;
		cout << "Surface:                        ";
		if ( curv_flag ) cout << "Curved" << endl;
		else cout << "Planar" << endl;
	}

	// Calculate the center of the vertices
	Vector3<double>	center = model_center_of_mass(model);
	center /= sam;

//	cout << "center=" << center << endl;

	// Calculate the minimum and maximum distances from the center
	for ( comp = model->comp; comp; comp = comp->next ) {
		comp->force((comp->location()/p->sampling(0)) - center);
		d = comp->force().length();
		if ( min > d ) min = d;
		if ( max < d ) max = d;
	}
	
	if ( curv_flag < 1 ) min = 0;

	center += ori;
	
//	cout << "center=" << center << endl;
	
	// Calculate the mask
	for ( i=z=nv=0; z<p->sizeZ(); z++ ) {
		if ( verbose & VERB_FULL )
			cout << "Calculating slice " << z << endl;
		vec[2] = z - center[2];
		for ( y=0; y<p->sizeY(); y++ ) {
			vec[1] = y - center[1];
			for ( x=0; x<p->sizeX(); x++, i++ ) {
				vec[0] = x - center[0];
				d = vec.length();
				if ( d < min ) p->set(i,1);
				else if ( d > max ) p->set(i,0);
				else if ( model_inside_outside(vec, model, curv_flag, fast) > 0 ) p->set(i,1);
				if ( (*p)[i] ) nv++;
			}
		}
	}
	
	if ( verbose )
		cout << "Mask volume:                    " << nv * sam.volume() << " A3" << endl << endl;
	
	return p;
}

/**
@brief 	Calculates a shell mask covering a model structure.
@param 	*model		model structure.
@param 	size		size of the new mask.
@param 	ori			origin of the new mask.
@param 	sam			pixel size of the new mask.
@param 	shell_width	width of shell mask.
@param 	curv_flag	flag to indicate curved surface.
@param 	fast		flag to use fast algorithm.
@return Bimage*		new mask.

	Each point in the new image is tested for inclusion in the mask,
	by calculating whether it falls inside the closest 3 vertices.
	Only the first model in the linked list is used.

**/
Bimage*		model_create_shell_mask(Bmodel* model, Vector3<long> size,
				Vector3<double> ori, Vector3<double> sam, double shell_width, int curv_flag, int fast)
{
	if ( sam[0] <= 0) sam = sam.max(1);

	long			i, x, y, z, nv=0, ni=0, no=0;
	double			d, min = 1e30, max = 0, half_width = shell_width/(2*sam[0]);
	if ( half_width < 0.5 ) half_width = 0.5;
	shell_width = 2*sam[0]*half_width;
	
	Vector3<double>	vec;
	Bcomponent*		comp;
	
	Bimage*			p = new Bimage(Float, TSimple, size, 1);
	if ( ori[0] == 0 && ori[1] == 0 && ori[2] == 0 )
		p->origin(p->default_origin());
	else
		p->origin(ori);
	ori = p->image->origin();
	p->sampling(sam);
	
	if ( verbose ) {
		cout << "Generating a shell mask:" << endl;
		cout << "Size:                           " << p->size() << endl;
		cout << "Origin:                         " << p->image->origin() << endl;
		cout << "Sampling:                       " << p->sampling(0) << " A/voxel" << endl;
		cout << "Shell width:                    " << shell_width << " A" << endl;
		if ( curv_flag ) cout << "Curved" << endl;
		else cout << "Planar" << endl;
	}

	// Calculate the center of the vertices
	Vector3<double>	center = model_center_of_mass(model);
	center /= sam;
	
	// Calculate the minimum and maximum distances from the center
	for ( comp = model->comp; comp; comp = comp->next ) {
		comp->force((comp->location()/p->sampling(0)) - center);
		d = comp->force().length();
		if ( min > d ) min = d;
		if ( max < d ) max = d;
	}

	min -= half_width;
	max += half_width;

	if ( curv_flag < 1 ) min = 0;

	center += ori;
	
	// Calculate the mask
	for ( i=z=0; z<p->sizeZ(); z++ ) {
		if ( verbose & VERB_FULL )
			cout << "Calculating slice " << z << endl;
		vec[2] = z - center[2];
		for ( y=0; y<p->sizeY(); y++ ) {
			vec[1] = y - center[1];
			for ( x=0; x<p->sizeX(); x++, i++ ) {
				vec[0] = x - center[0];
				d = vec.length();
//				data[i] = 0;
				if ( d < min ) ni++;
				else if ( d > max ) no++;
				else {
					d = model_inside_outside(vec, model, curv_flag, fast);
					if ( fabs(d) <= half_width ) p->set(i,1);
					else if ( d > 0 ) ni++;
					else no++;
				}
				if ( (*p)[i] ) nv++;
			}
		}
	}
	
	if ( verbose ) {
		cout << "Mask volume:                    " << nv * sam.volume() << " A3" << endl;
		cout << "Inner volume:                   " << ni * sam.volume() << " A3" << endl;
		cout << "Outer volume:                   " << no * sam.volume() << " A3" << endl << endl;
	}
	
	return p;
}

/**
@brief 	Calculates a mask with one level per model.
@param 	*model		model structure.
@param 	size		size of the new mask.
@param 	ori			origin of the new mask.
@param 	sam			pixel size of the new mask.
@return Bimage*		new mask.

	A level is defined as those voxels closest to the vertices of one model.

**/
Bimage*		model_create_level_mask(Bmodel* model, Vector3<long> size,
				Vector3<double> ori, Vector3<double> sam)
{
	if ( sam[0] <= 0) sam = sam.max(1);

	Vector3<double>	loc;
	Bmodel*			mp;
	Bcomponent*		comp;
	
	Bimage*			p = new Bimage(Integer, TSimple, size, 1);
	if ( ori[0] == 0 && ori[1] == 0 && ori[2] == 0 )
		p->origin(p->default_origin());
	else
		p->origin(ori);
	ori = p->image->origin();
	p->sampling(sam);
	
	if ( verbose ) {
		cout << "Generating a level mask:" << endl;
		cout << "Size:                           " << p->size() << endl;
		cout << "Origin:                         " << p->image->origin() << endl;
		cout << "Sampling:                       " << p->sampling(0) << " A/voxel" << endl;
	}

	int				level;
	
	// Calculate the levels
	for ( level=1, mp = model; mp; mp = mp->next, level++ ) {
		for ( comp = mp->comp; comp; comp = comp->next ) {
			loc[0] = (long) (comp->location()[0]/p->sampling(0)[0] + p->image->origin()[0] + 0.5);
			loc[1] = (long) (comp->location()[1]/p->sampling(0)[1] + p->image->origin()[1] + 0.5);
			loc[2] = (long) (comp->location()[2]/p->sampling(0)[2] + p->image->origin()[2] + 0.5);
			p->sphere(loc, comp->radius()/p->sampling(0)[0], 0, FILL_USER, level);
		}
	}
	
	return p;
}

/**
@brief 	Extracts segmentss from a multi-level mask around points defined by a model.  
@param 	*p				multi-level mask.
@param 	*model			model marking regions in mask.
@param 	multi_level		flag to retain multiple levels.
@return Bimage*			new mask from marked segments.

	Each segment should only be marked by a maximum of one component.

**/
Bimage*		img_extract_segments_using_model(Bimage* p, Bmodel* model, int multi_level)
{
	if ( !p || !model ) return NULL;
	if ( !p->data_pointer() ) return NULL;
	if ( !model->comp ) return NULL;
	
	Bmodel*				mp;
	Bcomponent*			comp;
	int					ncomp = model->component_count();
	
	long       i, j, k, x, y, z, n;
	
	if ( verbose ) {
		cout << "Extracting segments from a multi-level mask:" << endl;
		cout << "Model:                          " << model->identifier() << endl;
		cout << "Number of components:           " << ncomp << endl;
		if ( multi_level ) cout << "Multi-level" << endl;
		else cout << "Single level" << endl;
		cout << endl;
	}

	Bimage*     		pseg = new Bimage(Integer, TSimple, p->sizeX(), p->sizeY(), p->sizeZ(), p->images());
	pseg->sampling(p->sampling(0));
	for ( i=0; i<p->images(); i++ ) pseg->image[i] = p->image[i];
	
	int					level;
	long				max = (long) (p->maximum()+1);
	long				datasize = p->data_size();
	int*				lut = new int[max];
	for ( i=0; i<max; i++ ) lut[i] = 0;

	for ( k=1, mp = model; mp; mp = mp->next ) {
		for ( comp = model->comp; comp; comp = comp->next, k++ ) {
			n = mp->image_number();
			x = (long) (comp->location()[0]/p->sampling(0)[0] + p->image[n].origin()[0] + 0.5);
			y = (long) (comp->location()[1]/p->sampling(0)[1] + p->image[n].origin()[1] + 0.5);
			z = (long) (comp->location()[2]/p->sampling(0)[2] + p->image[n].origin()[2] + 0.5);
			if ( x < p->sizeX() && y < p->sizeY() && z < p->sizeZ() ) {
				i = p->index(0, x, y, z, n);
				level = (int) ((*p)[i] + 0.5);
				if ( level > 0 ) lut[level] = k;
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG img_extract_segments_using_model: " << x << tab << y << tab << z << tab << n << tab << level << endl;
			}
		}
    }

	for ( i=j=0; i<datasize; i++ ) {
		level = (int) ((*p)[i] + 0.5);
		if ( level > 0 && lut[level] ) {
			if ( multi_level ) pseg->set(i, lut[level]);
			else pseg->set(i, 1);
		}
	}

	delete[] lut;
	
	pseg->statistics();
	
	return pseg;
}

/**
@brief 	Adds to the mask the component locations.  
@param 	*p				mask.
@param 	*model			model.
@return int				0, <0 on error.

	At each component location, the voxel is added to the mask.

**/
int			img_add_model_to_mask(Bimage* p, Bmodel* model)
{
	if ( !p || !model ) return -1;
	if ( !p->data_pointer() ) return -1;
	if ( !model->comp ) return -1;
	
	p->change_type(Integer);
	
	Bmodel*				mp;
	Bcomponent*			comp;
	
	long       			i, x, y, z, n;
	int					new_value, count = (int) p->maximum();
	
	if ( verbose )
		cout << "Adding model " << model->identifier() << " to mask " << p->file_name() << endl << endl;

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
			n = mp->image_number();
			x = (long) (comp->location()[0]/p->sampling(0)[0] + p->image[n].origin()[0] + 0.5);
			y = (long) (comp->location()[1]/p->sampling(0)[1] + p->image[n].origin()[1] + 0.5);
			z = (long) (comp->location()[2]/p->sampling(0)[2] + p->image[n].origin()[2] + 0.5);
			if ( x < p->sizeX() && y < p->sizeY() && z < p->sizeZ() ) {
				i = p->index(0, x, y, z, n);
				if ( (*p)[i] < 1 ) {
					new_value = p->check_neighbors(i);
					if ( new_value > 0 ) p->set(i, new_value);
					else p->set(i, count++);
				}
			}
		}
    }
	
	p->maximum(count);
	
	return 0;
}

/**
@brief 	Generates components at the center of mass of each mask level.
@param 	*p				multi-level mask.
@return Bmodel*			new model.

	The mask is expected to be of integer data type.

**/
Bmodel*		model_from_multilevel_mask(Bimage* p)
{
	p->change_type(Integer);
	
	if ( verbose ) {
		cout << "Generating a model from a multi-level mask:" << endl;
		cout << "Mask file:                      " << p->file_name() << endl;
		cout << "Number of levels:               " << p->maximum() << endl << endl;
	}
	
	long			i, j, x, y, z, n;
	long			nlev = (long)p->maximum() + 1;
	double			u = p->sampling(0).volume();
	Bmodel*			model = NULL;
	Bmodel*			m = NULL;
	Bstring			id, comptype("SEG");
	Bcomponent*		comp = NULL;
	Bcomptype*		ct;
	Vector3<double>	coor;
	
	int*			num = new int[nlev];
	Vector3<float>*	loc = new Vector3<float>[nlev];
	
	int*			mask = (int *) p->data_pointer();
	
	for ( n=i=0; n<p->images(); n++ ) {
		id = Bstring(n, "Segments_%d");
//		m = model_add(&m, id);
//		if ( !model ) model = m;
		if ( m ) m = m->add(id.str());
		else model = m = new Bmodel(id.str());
		m->mapfile(p->file_name());
		m->image_number(n);
		ct = m->add_type(comptype);
		for ( j=0; j<=p->maximum(); j++ ) {
			num[j] = 0;
			loc[j] = 0;
		}
		for ( z=0; z<p->sizeZ(); z++ ) {
			coor[2] = (z - p->image[n].origin()[2])*p->sampling(0)[2];
			for ( y=0; y<p->sizeY(); y++ ) {
				coor[1] = (y - p->image[n].origin()[1])*p->sampling(0)[1];
				for ( x=0; x<p->sizeX(); x++, i++ ) if ( mask[i] > 0 ) {
					coor[0] = (x - p->image[n].origin()[0])*p->sampling(0)[0];
					num[mask[i]]++;
					loc[mask[i]] += coor;
				}
			}
		}
		comp = NULL;
		for ( j=1; j<=p->maximum(); j++ ) if ( num[j] ) {
//			id = Bstring(j, "%d");
//			comp = component_add(&comp, id);
//			if ( !m->comp ) m->comp = comp;
			if ( comp ) comp = comp->add(j);
			else m->comp = comp = new Bcomponent(j);
			comp->type(ct);
			comp->location(loc[j]/num[j]);
			comp->radius(pow(3*num[j]*u/(4*M_PI), 1.0/3.0));
			comp->view(View2<float>(comp->location()));
			comp->color()[0] = 0;
		}
	}

	delete[] num;
	delete[] loc;
	
	return model;
}

/**
@brief 	Calculates a set of 2D masks for the projections of a model structure.
@param 	*model		model structure.
@param 	size		size of the new mask.
@param 	ori			origin of the new mask.
@param 	sam			pixel size of the new mask.
@param 	dang		angular step size.
@param 	sym			point group symmetry.
@return Bimage*		new mask.

	Each component is used to set  the corresponding projected pixel in each 2D image.
	The projection directions are calculated within the asymmetric unit.

**/
Bimage*		model_create_projected_mask(Bmodel* model,
				Vector3<long> size, Vector3<double> ori,
				Vector3<double> sam, double dang, Bsymmetry& sym)
{
	if ( sam.volume() <= 0) sam = sam.max(1);

	View*			v;
	View*			views = asymmetric_unit_views(sym, dang, dang, 1);
	
	long			i, j, n, nm;
	long			nv = count_list((char *)views);
	double			f, fa(0);

	Bmodel*			mp = model;
	Bcomponent*		comp;
	Matrix3			mat;
	Vector3<double>	loc;

	size[2] = 1;
	Bimage*			p = new Bimage(UCharacter, TSimple, size, nv);
	if ( ori[0] == 0 && ori[1] == 0 && ori[2] == 0 )
		p->origin(p->default_origin());
	else
		p->origin(ori);
	ori = p->image->origin();
	sam[2] = 1;
	p->sampling(sam);

	if ( verbose ) {
		cout << "Generating projection masks:" << endl;
		cout << "Size:                           " << p->size() << endl;
		cout << "Origin:                         " << p->image->origin() << endl;
		cout << "Sampling:                       " << p->sampling(0) << " A/voxel" << endl;
		cout << "Angular step size:              " << dang*180.0/M_PI << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Projections:                    " << nv << endl;
	}

	if ( verbose & VERB_PROCESS )
		cout << "Mask\tView\t\t\t\t\tCount\tFraction" << endl;
	for ( v = views, n=0; v; v = v->next, ++n ) {
		mat = v->matrix();
		for ( mp = model; mp; mp = mp->next ) {
			for ( comp = mp->comp; comp; comp = comp->next ) {
				loc = (mat * comp->location())/p->image->sampling() + ori + 0.5;
				loc[2] = 0;
				if ( !p->within_boundaries(loc) )
					cerr << "Out of bounds! " << loc << endl;
				i = p->index(loc, n);
				p->set(i, 1);
			}
		}
		for ( i=n*p->image_size(), j=0, nm=0; j<p->image_size(); ++i, ++j )
			nm += (*p)[i];
		f = nm*1.0/p->image_size();
		fa += f;
		if ( verbose & VERB_PROCESS )
			cout << n+1 << setprecision(4) << tab << *v << tab
				<< nm << tab << f << endl;
	}
	
	fa /= n;
	
	// Calculate the mask area
	for ( i=nv=0; i<p->data_size(); i++ ) if ( (*p)[i] > 0.5 ) nv++;
	
	if ( verbose ) {
		cout << "Mask area:                      " <<
			nv * sam.volume()/n << " A2" << endl;
		cout << "Average mask fraction:          " <<
			fa << endl << endl;
	}
	
	return p;
}
