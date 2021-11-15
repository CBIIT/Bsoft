/**
@file	bafm.cpp 
@brief	Simulation of AFM experiments
@author Bernard Heymann
@date	Created: 19990124
@date	Modified: 20170612
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
Bimage* 	afm_generate_tip(Vector3<long> size, Vector3<double> sampling,
				double radius, double angle, double resolution);
Bplot* 		afm_simulate(Bimage* ptip, Bimage* p, double spring,
				double modulus, double thick);
Bimage* 	afm_force_to_height(Bimage* pf, double force_step);

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bafm [options] input.map height.img",
"------------------------------------------",
"Simulates AFM imaging.",
"An optional force map is produced (where every 2D slice is a force profile ",
"		at a cantilever height).",
"A set of height images is written.",
" ",
"Actions:",
"-invert                  Invert inside & outside (default not).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5            Sampling (default 1 A/pixel).",
" ",
"Tip generation:",
"-tipradius 5.3           Generate tip: radius (angstrom).",
"-tipsize 10,12,5         Tip size (x,y,z pixels).",
"-tipangle 55.7           Tip side angle (degrees).",
" ",
"Simulation parameters:",
"-spring 0.1              Spring constant (default 0.05 N/m).",
"-modulus 2e8             Bulk modulus (default 1e9 N/m2).",
"-thickness 55            Thickness (default 50 A).",
"-forcestep 15            Force intervals for height maps (default 10 pN).",
" ",
"Input:",
"-Tip tip.map             Read tip from file.",
" ",
"Output:",
"-Postscript file.ps      Postscript force curve.",
"-Force force.map         Force map file name.",
" ",
NULL
};

int 	 	main(int argc, char **argv)
{
    /* Initialize variables */
    Vector3<long>	tipsize = {1,1,1};		// Tip dimensions
	double			tipradius(0); 			// Tip radius also flag to generate tip
	double			tipangle(45);			// Tip side angle
	double			resmin(1 + sqrt(5.0));	// Resolution for tip generation
    Vector3<double>	sam;			// Unit dimension scaling (Angstrom/pixel)
    int 			setinvert(0); 			// Don't invert inside-outside
    double			spring(0.05); 			// Spring constant (typically 0.01-0.1 N/m)
    double			modulus(1e9); 			// Bulk modulus (typically ~10 pN/A2 ???)
	double			thick(50);	 			// Thickness of sample in Angstrom
	double			forcestep(10);			// Force intervals for height maps
    Bstring			tipfile;				// Tip file name
	Bstring			psfile;					// Postscript file name
	Bstring			forcefile;				// Force map file name
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "tipradius" )
			if ( ( tipradius = curropt->value.real() ) < 1 )
				cerr << "-tipradius: A radius must be specified!" << endl;
		if ( curropt->tag == "tipsize" )
			tipsize = curropt->size();
		if ( curropt->tag == "tipangle" )
			if ( ( tipangle = curropt->value.real() ) < 1 )
				cerr << "-tipangle: An angle must be specified!" << endl;
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "spring" )
			if ( ( spring = curropt->value.real() ) < 0.00000001 )
				cerr << "-spring: The spring constant must be specified!" << endl;
		if ( curropt->tag == "modulus" )
			if ( ( modulus = curropt->value.real() ) < 0.001 )
				cerr << "-modulus: The bulk modulus must be specified!" << endl;
		if ( curropt->tag == "thickness" )
			if ( ( thick = curropt->value.real() ) < 1 )
				cerr << "-thickness: The thickness must be specified!" << endl;
		if ( curropt->tag == "forcestep" )
			if ( ( forcestep = curropt->value.real() ) < 0.001 )
				cerr << "-forcestep: The force step must be specified!" << endl;
		if ( curropt->tag == "Tip" )
			tipfile = curropt->filename();
 		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
		if ( curropt->tag == "Force" )
			forcefile = curropt->filename();
    }
	option_kill(option);
    
	double		ti = timer_start();
	
    // Read the input file if a file name is given
	Bimage*		p = read_img(argv[optind++], 1, 0);
	Bimage* 	ptip = NULL;
        
    if ( !p ) {
		cerr << "Error: No input file given!" << endl;
		bexit(-1);
	}

	p->change_type(Float);
	
	p->truncate_to_min_max(0, p->maximum());
	
	if ( sam.volume() > 0 ) p->sampling(sam);
	else sam = p->sampling(0);
	
	// Tip density generation
	if ( tipradius > 0 ) {
		tipangle *= M_PI/180.0;
		ptip = afm_generate_tip(tipsize, sam, tipradius, tipangle, resmin);
		write_img("tip.map", ptip, 0);
//		bexit(0);
	} else {
		if ( tipfile.length() ) {
			ptip = read_img(tipfile, 1, 0);
			ptip->change_type(Float);
		}
	}
	
	// Invert if requested
	if ( setinvert ) p->invert();
	
	// Simulate an AFM scan if a tip is given
	Bplot*		plot = NULL;
	if ( ptip )
		plot = afm_simulate(ptip, p, spring, modulus, thick);
	
	delete ptip;

	if ( psfile.length() && plot ) ps_plot(psfile, plot);
	
	delete plot;
	
	// Copy the command line to the label field of the image structure
	
    // Write an output force map file if a file name is given
	if ( forcefile.length() )
		write_img(forcefile, p, 0);

	// Calculate height images
	Bimage* ph = afm_force_to_height(p, forcestep);
	
	delete p;
	
    // Write an output file if a file name is given
	if ( argc > optind && strspn(argv[optind],"-") != 1 )
		write_img(argv[optind], ph, 0);

	delete ph;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/**
@brief 	Generates a rounded conical AFM tip.
@param 	size		size of the tip density
@param 	sampling	sampling/spacing (in angstrom/pixel)
@param 	radius		radius of the tip curvature (in angstrom)
@param 	angle		angle of the tip side (in radians)
@param 	resolution	the resolution affects the tip surface softness
@return Bimage* 	tip 3D map.

	A 2D image of a conical tip is generated, with a 45 degree angle and
	a rounding tip with the desired radius and softness.
	The 2D image is then converted to a 3D surface using the function 
	img_to_surface and returned.  The density of beta-silicon nitride
	of 1.925 Da/A3 is used.

**/
Bimage* 	afm_generate_tip(Vector3<long> size, Vector3<double> sampling,
				double radius, double angle, double resolution)
{
	long			i, xx, yy;
	double			dx, dy, dz, dis, dis2;
	
	if ( size[0] < 4*radius/sampling[0] ) size[0] = (long) (4*radius/sampling[0]);
	if ( size[1] < 4*radius/sampling[1] ) size[1] = (long) (4*radius/sampling[1]);
	if ( size[2] < 4*radius/sampling[2] ) size[2] = (long) (4*radius/sampling[2]);
	
    if ( verbose & VERB_PROCESS ) {
		cout << "Generating a 2D tip image" << endl;
    	cout << "Tip dimensions:                 " << size << endl;
    	cout << "Tip radius:                     " << radius << " A" << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Generating a 2D tip image" << endl;
	
	// Generate a 2D tip height image to be converted to 3D later
    Bimage*		p = new Bimage(Float, TSimple, size[0], size[1], 1, 1);
	p->sampling(sampling);
	p->origin(p->default_origin());
	
    radius /= sampling[0]; 					// Convert the radius to pixel units
	double		radcut(radius/sqrt(2.0));		// Edge between sphere and cone
//	double		radcut(radius*cos(angle));		// Edge between sphere and cone
    double		rad2(radius*radius);
	double		max(sqrt(1.0*p->sizeX()*p->sizeX() + p->sizeY()*p->sizeY())/2);
//	double		ta(1/tan(angle));
   
    // Generate 45 degree cone with spherical tip with input radius
	for ( i=yy=0; yy<p->sizeY(); yy++ ) {
		dy = yy - p->image->origin()[1];
		for ( xx=0; xx<p->sizeX(); xx++, i++ ) {
			dx = xx - p->image->origin()[0];
    		dis2 = dx*dx + dy*dy;
    		dis = sqrt(dis2);
    		if ( dis < radcut )
    			dz = max - (2*radcut - sqrt(rad2 - dis2));
//    			dz = max - radius/sin(angle) + sqrt(rad2 - dis2);
    		else
    			dz = max - dis;
//   			dz = max - dis*ta;
//			if ( dz < 0 ) dz = 0;
			p->set(i, dz);
    	}
    }
	
	p->statistics();
	
//	p->show_maximum((sqrt(2.0) - 1)*radius*p->sampling(0)[2]);
//	p->show_maximum(cos(angle)*radius*p->sampling(0)[2]);
//	p->show_maximum(sin(angle)*p->maximum()*p->sampling(0)[2]);
//	p->show_maximum((1-cos(angle))*p->maximum()*p->sampling(0)[2]);
	p->show_maximum(radius*p->sampling(0)[2]);
	p->show_minimum(p->maximum()*p->sampling(0)[2]/sin(angle));
    
	Bimage*		ptip = p->topograph_to_surface(NULL, size[2], 1.925, resolution);

	ptip->image->origin(ptip->size()/2);
	
	delete p;
	
    return ptip;
}

int 	afm_simulate_at_xy(Bimage* ptip, Bimage* p, long ii,
			float* force, double spring, double modulus, double thick)
{
	long			i, ti, di, xx, yy, zz, tx, ty, tz, hiz;
	long			dx, dy, dz, zshift;
	long			ss(p->sizeX()*p->sizeY()), tss(ptip->sizeX()*ptip->sizeY());
	double			fraction, overlap_volume(1);
	double			scale(p->sampling(0).volume()/(p->maximum()*ptip->maximum()));	// Adjustment for correct units
	
	// Tip displacement = proportion * volume overlap
	double			proportion(scale*modulus/(spring*thick));

	double*			shift = new double[p->sizeZ()];
	
	yy = ii/p->sizeX();
	xx = ii - yy*p->sizeX();
	for ( zz=0; zz<p->sizeZ(); zz++ ) {
		shift[zz] = 0;
		if ( overlap_volume > 1e-30 ) {
			overlap_volume = 0;
			i = zz*ss + ii;
			hiz = ptip->sizeZ();
			if ( hiz > p->sizeZ() - zz ) hiz = p->sizeZ() - zz;
			for ( ti=ty=0; ty<ptip->sizeY(); ty++ ) {	// Calculate overlap volume
				dy = yy + ty - (long)ptip->image->origin()[1];
				while ( dy < 0 ) dy += p->sizeY();		// Wrap y
				while ( dy >= p->sizeY() ) dy -= p->sizeY();
				for ( tx=0; tx<ptip->sizeX(); tx++, ti++ ) {
					dx = xx + tx - (long)ptip->image->origin()[0];
					while ( dx < 0 ) dx += p->sizeX();		// Wrap x
					while ( dx >= p->sizeX() ) dx -= p->sizeX();
					di = dy*p->sizeX()+dx;
					for ( tz=0; tz<hiz; tz++ ) {
						dz = zz + tz;
						overlap_volume += (*p)[dz*ss+di]*(*ptip)[tz*tss+ti];
					}
				}
			}
		}
		shift[zz] = overlap_volume*proportion;
//		cout << "%ld\t%g\n", z, shift[z]);
	}
	
	for ( zz=0; zz<p->sizeZ(); zz++ ) {	// Calculate zshift where the forces are equal
		zshift = 0;
		while ( zshift*p->sampling(0)[2] < shift[zz+zshift] &&
						zshift+zz < p->sizeZ() ) zshift++;
		if ( zshift > 0 ) zshift--;
		fraction = (zshift*p->sampling(0)[2] - shift[zz+zshift])/
						(zshift*p->sampling(0)[2] - shift[zz+zshift] -
						(zshift+1)*p->sampling(0)[2] + shift[zz+zshift+1]);
		if ( !isfinite(fraction) )
			cout << "Problem:\t" << xx << tab << yy << tab << zz << tab <<
					fraction << tab << zshift*p->sampling(0)[2] << tab <<
					shift[zz+zshift] << tab << shift[zz+zshift+1] << endl;
		i = zz*ss + ii;
		force[i] = spring*(zshift + fraction)*p->sampling(0)[2];
	}
//	cout << "%ld\t%ld\t%ld\n", x, y, zshift);

	delete[] shift;
	
	return 0;
}

/**
@brief 	Simulates and AFM experiment.
@param 	*ptip		AFM tip density map - converted to force map.
@param 	*p			specimen density map.
@param 	spring		AFM cantilever spring constant (N/m).
@param 	modulus		bulk modules (N/m2).
@param 	thick 		sample thickness (angstrom).
@return Bplot* 		force curve.

	The elastic force on an AFM tip is calculated as:
		F = kt*dz = kb*dV/d
	where	kt is the cantilever spring constant, typically 0.01 - 0.1 N/m.
			dz is the tip displacement due to specimen interaction.
			kb is the bulk modulus of the specimen, typically 1e9 N/m2.
			dV is the volume of the specimen displaced by the tip. This
				volume is here approximated as the density overlap volume.
			d is the specimen thickness in angstrom.
	The tip is positioned at each point (x,y,z) in the density map and the
	overlap volume calculated (i.e., a type of convolution).  This tip
	position corresponds to a zero displacement.  The correct tip position
	is found by shifting the tip upwards (in the z-direction) until the 
	displacement and elastic forces are approximately equal.  The tip 
	displacement is refined by interpolation and the force calculated.  
	The force map is returned in place of the original density map.
	Forces are calculated in piconewton. 
	The origin of the tip density is taken as nx/2, ny/2, 0.

**/
Bplot*		afm_simulate(Bimage* ptip, Bimage* p, double spring, double modulus, double thick)
{
	long			i, ii, zz;
	long			ss(p->sizeX()*p->sizeY());
	double			avg_force;
	
    if ( verbose & VERB_PROCESS ) {
		cout << "Simulating an AFM experiment" << endl;
    	cout << "Spring constant:                " << spring << " N/m" << endl;
		cout << "Bulk modulus:                   " << modulus << " N/m2" << endl;
    	cout << "Specimen thickness:             " << thick << " A" << endl;
    	cout << "Density scale:                  " << p->sampling(0).volume()/(p->maximum()*ptip->maximum()) << " A/pixel" << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Simulating an AFM experiment" << endl << endl;
	
	// Conversions to Angstrom and pN
	spring *= 100;			// 100 pN/A = N/m
	modulus *= 1e-8;		// 0.01 aN/A2 = N/m2

	float*			force = new float[p->sizeX()*p->sizeY()*p->sizeZ()];
	
#ifdef HAVE_GCD
	dispatch_apply(ss, dispatch_get_global_queue(0, 0), ^(size_t i){
	 	afm_simulate_at_xy(ptip, p, i, force, spring, modulus, thick);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<ss; i++ )
	 	afm_simulate_at_xy(ptip, p, i, force, spring, modulus, thick);
#endif

	p->maximum(0);
	for ( i=0; i<ss*p->sizeZ(); i++ )
		if ( p->maximum() < force[i] ) p->maximum(force[i]);
	
	Bstring			title("AFM simulation");
	int				ncol(3), nv(p->sizeZ());
	Bplot*			plot = new Bplot(1, nv, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Distance");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(2);
	plot->page(0).column(1).label("Deflection");
	plot->page(0).column(1).axis(4);
	plot->page(0).column(2).type(2);
	plot->page(0).column(2).label("Force");
	plot->page(0).column(2).axis(3);
//	plot->page(0).axis(1).min(0);
//	plot->page(0).axis(1).max(360);
//	plot->page(0).axis(1).inc(30);
	
	for ( zz=0; zz<nv; zz++ ) {
		for ( i=0, ii=zz*ss, avg_force=0; i<ss; i++, ii++ )
			avg_force += force[ii];
		avg_force /= ss;
		(*plot)[zz] = zz*p->sampling(0)[2];
		(*plot)[zz+nv] = avg_force/spring;
		(*plot)[zz+2*nv] = avg_force;
	}

	Bstring		txt;
	txt = "Cantilever spring constant:   " + Bstring(spring/100, "%lg N/m");
	plot->page(0).add_text(txt);
	txt = "Bulk modulus:                 " + Bstring(modulus*1e8, "%lg N/m2");
	plot->page(0).add_text(txt);
	txt = "Sample thickness:             " + Bstring(thick, "%lg A");
	plot->page(0).add_text(txt);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Distance(A)\tDeflection(A)\tForce(pN)" << endl;
		for ( zz=0; zz<nv; zz++ )
			cout << setw(12) << (*plot)[zz] << tab << setw(12)
				<< (*plot)[zz+nv] << tab << (*plot)[zz+2*nv] << endl;
		cout << endl;
	}
	
	p->data_assign((unsigned char *) force);

    return plot;
}

/**
@brief 	Converts an AFM force map into a set of 2D height images at different constant force values.
@param 	*pf			force map.
@param 	force_step	force intervals for height images
@return *ph			height images.
**/
Bimage* 	afm_force_to_height(Bimage* pf, double force_step)
{
	long			i, j, xx, yy, zz, nn;
	long			ss(pf->sizeX()*pf->sizeY());
	double			force, denom, diff, scale, min, max;
	
	nn = (long) (pf->maximum()/force_step); 	// A height image every 10 pN
	Bimage*			ph = pf->copy_header(nn);
	ph->sizeZ(1);
	ph->data_type(Float);
	ph->data_alloc();
	
    if ( verbose & VERB_PROCESS ) {
		cout << "Converting an AFM force map into " << ph->images() << " height images" << endl;
		cout << "Force(pN)\tHeight_min(A)\tHeight_range(A)" << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Converting an AFM force map into " << ph->images() << " height images" << endl << endl;
	
	for ( nn=0; nn<ph->images(); nn++ ) {
		force = force_step*(nn+1);
		min = 1e30;
		max = -1e30;
   		for ( i=nn*ss, yy=0; yy<ph->sizeY(); yy++ ) {
			for ( xx=0; xx<ph->sizeX(); xx++, i++ ) {
 				zz = 1;
				j = yy*pf->sizeX() + xx;
				while ( zz < pf->sizeZ() && (*pf)[zz*ss + j] > force )
					zz++;
				j += zz*ss;
				denom = (*pf)[j-ss] - (*pf)[j];
				diff = force - (*pf)[j];
//				cout << zz << tab << denom << tab << diff << endl;
//				if ( diff > 0 && denom > 0 )
//					ph->set(i, (*ph)[i] - ph->sampling(0)[2]*diff/denom);
				ph->set(i, ph->sampling(0)[2]*zz);
				if ( denom != 0 ) {
					diff /= denom;
					if ( fabs(diff) < 1 )
						ph->set(i, (*ph)[i] - ph->sampling(0)[2]*diff);
				}
				if ( min > (*ph)[i] ) min = (*ph)[i];
				if ( max < (*ph)[i] ) max = (*ph)[i];
			}
		}
		scale = max - min;
    	if ( verbose & VERB_PROCESS )
			cout << setw(12) << force << tab << setw(12) << min << tab << scale << endl;
		for ( i=0, j=nn*ss; i<ss; i++, j++ )
				ph->set(j, (*ph)[j] - min);
	}
	cout << endl;
	
	ph->statistics();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG afm_force_to_height: Done!" << endl << endl;
	
    return ph;
}
