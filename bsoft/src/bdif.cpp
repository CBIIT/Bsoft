/**
@file	bdif.cpp
@brief	Analyze diffraction patterns
@author Bernard Heymann
@date	Created: 20050217
@date	Modified: 20220524
**/

#include "rwimg.h"
#include "matrix_linear.h"
#include "simplex.h"
#include "math_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
int		img_fix_powerspectrum_cross(Bimage* p, int flag, double fixtol);
int		img_orthogonal_radial_profiles(Bimage* p, double angle);
int		img_analyze_reflection(Bimage* p, double ref_res, double threshold, long kedge, int sym);

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bdif [options] [input.img] output.img",
"--------------------------------------------",
"Analyzes diffraction patterns.",
" ",
"Actions:",
"-fixcross x,2.5          Fix high values in the cross of a powerspectrum (x/y/xy, default tolerance 2).",
"-angle 36.8              Calculate orthogonal radial profiles rotated with this angle.",
"-reflection 2.5,15       Analyze reflections: spatial frequency (angstrom) and threshold (default 2).",
"-findorigin              Find the origin of the diffraction pattern.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-origin 60.4,-32.2,44    Origin of image (voxels).",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-resolution 4.5,130      Resolution range for correlation (default 0 - 1e6 angstrom).",
" ",
NULL
};

int 	main(int argc, char* argv[])
{
    // Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double> sam;						// Units for the three axes (A/pixel)
	int				fixcross(0);				// Flag to fix the cross in a power spectrum
	double			fixtol(2);					// Tolerance to fix the cross in a power spectrum
	double			angle(0);					// Angle for directional averages
	int				profiles(0);				// Flag to do orthogonal profiles
	double			ref_res(0);					// Reflection resolution
	double			threshold(2);				// Threshold for selecting reflections
	int				kedge(5);					// Kernel edge size
	int				sym(6);						// Cyclic symmetry order
	int				find_origin(0);				// Flag to find diffraction pattern origin
	int				set_origin(0);
	Vector3<double>	origin;						// Origin in file
	double			hires(0), lores(0);			// Limiting resolution range (hires must be > 0 to be set)

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "fixcross" ) {
			if ( curropt->value.contains("x") ) fixcross |= 1;
			if ( curropt->value.contains("y") ) fixcross |= 2;
			if ( curropt->value.contains(",") ) fixtol = curropt->value.post(',').real();
		}
		if ( curropt->tag == "angle" ) {
        	angle = curropt->value.real() * M_PI/180.0;
        	profiles = 1;
		}
		if ( curropt->tag == "reflection" )
			if ( curropt->values(ref_res, threshold, kedge, sym) < 1 )
				cerr << "-reflection: A reflection resolution must be specified!" << endl;
		if ( curropt->tag == "findorigin" )
			find_origin = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
    }
	option_kill(option);
        
	double		ti = timer_start();
	
    Bimage* 	p = read_img(argv[optind++], 1, -1);
	Bimage*		pmask = NULL;
	
	if ( !p ) {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	if ( nudatatype == Unknown_Type ) nudatatype = p->data_type(); // Preserve the old type
	
	if ( nudatatype >= Float ) p->change_type(nudatatype);		// Convert for efficiency
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}
		
	if ( sam.volume() ) p->sampling(sam);

	if ( hires < p->sampling(0)[0] ) hires = p->sampling(0)[0]*10;
	if ( lores < hires ) lores = hires*100;
	
	if ( fixcross ) img_fix_powerspectrum_cross(p, fixcross, fixtol);
	
	double			radius(p->sizeX()/4.0);
	if ( find_origin ) p->find_center(pmask, hires, lores, radius, 0, 1);

	if ( profiles ) img_orthogonal_radial_profiles(p, angle);
	
	if ( ref_res ) img_analyze_reflection(p, ref_res, threshold, kedge, sym);
	
    // Write an output file if a file name is given
    if ( optind < argc ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}
	
	delete p;
	delete pmask;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

int		img_fix_powerspectrum_cross(Bimage* p, int flag, double fixtol)
{
	long		i, xx, yy, x(p->sizeX()), y(p->sizeY()), hx(x/2), hy(y/2), npx(0);
	double		a, r(fixtol/2);
	
	if ( verbose )
		cout << "Fixing the powerspectrum cross with tolerance " << fixtol << endl;
	
	if ( flag & 1 ) {
		if ( verbose ) cout << "Direction x" << endl;
		for ( i=hy*x, xx=0; xx<x; ++xx, ++i ) {
			a = (*p)[i-x] + (*p)[i+x];
			if ( (*p)[i] > r*a ) {
				p->set(i, a/2);
				npx++;
			}
		}
	}
	
	if ( flag & 2 ) {
		if ( verbose ) cout << "Direction y" << endl;
		for ( i=hx, yy=0; yy<y; ++yy, i+=x ) {
			a = (*p)[i-1] + (*p)[i+1];
			if ( (*p)[i] > r*a ) {
				p->set(i, a/2);
				npx++;
			}
		}
	}
	
	if ( verbose ) cout << "Pixels fixed:                " << npx << endl << endl;
	
	return 0;
}

int		img_orthogonal_radial_profiles(Bimage* p, double angle)
{
	long			i, len(p->sizeX()/2);
	double			s;
	Vector3<double>	end(cos(angle)*len, sin(angle)*len, 0);
	Bimage*			prad = p->radial(0, len, 1);
	Bimage*			pdir1 = p->extract_line(0, p->image->origin(), end, 5);;
	end = Vector3<double>(sin(angle)*len, cos(angle)*len, 0);
	Bimage*			pdir2 = p->extract_line(0, p->image->origin(), end, 5);;
	
	cout << "Radius\ts\tr\tAvg\tDir1\tDir2" << endl;
	for ( i=1; i<prad->sizeX(); i++ ) {
		s = i/p->real_size()[0];
		cout << i << tab << s << tab << 1/s << tab << (*prad)[i] << tab << (*pdir1)[i] << tab << (*pdir2)[i] << endl;
	}
	cout << endl;
	
	delete prad;
	delete pdir1;
	delete pdir2;

	return 0;
}


Vector3<double>	sinc_fit(vector<Vector3<double>> coor, vector<double> value, double& R)
{
	R = 0;
	
	long			i;
	double			k(2), d, R1;
	Vector3<double>	ori, bori, off;
	
	for ( ori[1]=-k; ori[1]<=k; ori[1]+=0.1 ) {
		for ( ori[0]=-k; ori[0]<=k; ori[0]+=0.1 ) {
			for ( i=0, R1=0; i<coor.size(); ++i ) {
				off = coor[i] - ori;
				d = value[i]*sinc(off[0])*sinc(off[1]);
				R1 += d;
			}
			R1 = sqrt(R1/coor.size());
//			cout << ori << tab << R << endl;
			if ( R < R1 ) {
				R = R1;
				bori = ori;
			}
		}
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Best origin:            " << bori << tab << R << endl;
	}
	
	return bori;
}

Vector3<double>	img_sinc_fit(Bimage* p, Vector3<double> loc, double bkg, double peak, double& R)
{
	long					i, xx, yy, k(2);
//	double					sum(0);
	vector<Vector3<double>>	coor;
	vector<double>			value;
	
	for ( yy=loc[1]-k; yy<=loc[1]+k; ++yy ) {
		for ( xx=loc[0]-k; xx<=loc[0]+k; ++xx ) {
			coor.push_back(Vector3<double>(xx,yy,0)-loc);
			i = p->index(xx, yy);
			value.push_back(((*p)[i]-bkg)/peak);
//			sum += value.back();
		}
	}
	
	loc += sinc_fit(coor, value, R);
	
	return loc;
}

double		ellipse_R(Bsimplex& simp)
{
	long			i;
	double			R(0), df;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); i++ ) {
		df = simp.parameter(0)*x[i]*x[i] + simp.parameter(1)*x[i]*f[i] + simp.parameter(2)*f[i]*f[i] - simp.constant(0);
		R += df*df;
	}
	
	R = sqrt(R/i);
			
	return R;
}

double		fit_ellipse_iter(vector<Vector3<double>> v)
{
	long			i;
	double			d2, davg(0), dstd(0);
	vector<double>	vx, vy;
	
	for ( auto p: v ) {
		d2 = p.length2();
		davg += d2;
		dstd += d2*d2;
		vx.push_back(p[0]);
		vy.push_back(p[1]);
	}
	davg /= v.size();
	dstd = sqrt(dstd/v.size() - davg*davg);
	
//	cout << davg << tab << dstd << endl;
	
	Bsimplex		simp(1, 3, 1, vx.size(), vx, vy);
	
	simp.constant(0, davg);
	simp.parameter(0, 1);
	simp.parameter(1, 0);
	simp.parameter(2, 1);
	simp.limits(0, 0.9, 1.1);
	simp.limits(1, -0.1, 0.1);
	simp.limits(2, 0.9, 1.1);
	
	double			R = simp.run(1000, 1e-6, ellipse_R, 100);
	R /= dstd;
	
	vector<double>	c(3,0);
	for ( i=0; i<3; ++i ) c[i] = simp.parameter(i);
	
	double			a = atan((c[2]-c[0]-sqrt((c[0]-c[2])*(c[0]-c[2])+c[1]*c[1]))/c[1]);
	
	if ( verbose ) {
		cout << "Constant: " << simp.constant(0) << endl;
		cout << "Coefficients:" << endl;
		for ( i=0; i<3; ++i ) cout << c[i] << endl;
		cout << "Angle: " << a*180.0/M_PI << endl;
		cout << c[1]/(2*sin(a)*cos(a)) << endl;
		cout << "R: " << R << endl << endl;
	}
	
	return R;
}

int			img_analyze_reflection(Bimage* p, double ref_res, double threshold, long kedge, int sym)
{
	if ( threshold < 0.1 ) threshold = 20;
	if ( kedge < 1 ) kedge = 1;
	if ( sym < 1 ) sym = 1;
	
	long				i, j, na;
	long				ks(kedge/2);				// Summation kernel half edge and volume
	double				a, sf;
	double				s(1/ref_res);				// Spatial frequency
	Vector3<double>		k(p->real_size()*s);		// Frequency space pixel distance
	k[2] = 0;
	double				rd(1+5/k.length());			// Relative displacement for reference kernel
	double				angsam(ks);					// Angular sampling in pixels
	double				da(angsam/k.length());		// Angle increment
	double				pmax, ravg, bkg(0), tR;
	Vector3<double>		tloc, rloc, sloc;
	vector<Vector3<double>> loc;
	vector<double>		R, Rfit;
	
	if ( verbose ) {
		cout << "Analyzing a reflection:" << endl;
		cout << "Spatial frequency:              " << s << " 1/A" << endl;
		cout << "Threshold:                      " << threshold << " e/px" << endl;
		cout << "Kernel edge size:               " << 2*ks+1 << endl;
		cout << "Cyclic symmetry:                " << sym << endl;
		cout << "Frequency space pixels:         " << k << endl;
		cout << "Frequency space pixel size:     " << 1/p->real_size()[0] << " 1/Å" << endl;
		cout << "Angular increment:              " << da*180.0/M_PI << " degrees" << endl << endl;
	}

	for ( a=0, j=na=0; a<TWOPI; a+=da, ++na ) {
//		tloc = Vector3<double>(p->image->origin()[0] + k[0]*cos(a), p->image->origin()[1] + k[1]*sin(a), 0);
		tloc = Vector3<double>(k[0]*cos(a), k[1]*sin(a), 0);
		rloc = tloc * rd;
		tloc += p->image->origin();
		rloc += p->image->origin();
//		cout << na << tab << tloc << tab << rloc << endl;
		i = p->kernel_max(p->index(tloc,0), ks);
		pmax = (*p)[i];
		ravg = p->kernel_average(p->index(rloc,0), ks, 0, 1e30);
		bkg += ravg;
		if ( pmax > threshold ) {
			tloc = p->coordinates(i);
			tloc = img_sinc_fit(p, tloc, ravg, pmax - ravg, tR);
			if ( loc.size() && tloc.distance(loc.back()) < 2 ) {
				if ( pmax > R.back() ) {
					loc.back() = tloc;
					R.back() = pmax;
					Rfit.back() = tR;
				}
			} else {
				loc.push_back(tloc);
				R.push_back(pmax);
				Rfit.push_back(tR);
			}
		}
	}
	
	bkg /= na;
	
	if ( loc.size() < 1 ) {
		cerr << "Error: No reflections found, check the threshold." << endl;
		return -1;
	}
	
	if ( loc[0].distance(loc.back()) < 2 ) {
		if ( R[0] < R.back() ) {
			loc[0] = loc.back();
			R[0] = R.back();
			Rfit[0] = Rfit.back();
		}
		loc.pop_back();
		R.pop_back();
		Rfit.pop_back();
	}
	
	if ( verbose )
		cout << "#\tx\ty\tz\ta\tR\tRfit" << endl;
	for ( i=0; i<loc.size(); ++i ) {
		sloc = loc[i] - p->image->origin();
		if ( verbose )
			cout << setprecision(2) << i+1 << tab << loc[i] << tab <<
				atan2(sloc[1], sloc[0])*180.0/M_PI << tab << setprecision(3) << R[i] << tab << Rfit[i] << endl;
		R[i] = p->kernel_sum(p->index(loc[i],0), 1) - 9*bkg;
		loc[i] = sloc;
	}
	
	if ( verbose ) {
		cout << "Number of reflections found:    " << loc.size() << endl;
		cout << "Background intensity:           " << bkg << " e/px" << endl << endl;
	}
	
	long				n, ns(0);
	double				res_avg(0), res_std(0);
	Matrix3				mat(Vector3<double>(0,0,1), TWOPI/sym);
	vector<int>			sel(loc.size(), 0);
	
	if ( verbose )
		cout << "x\ty\tz\tAngle\tPeak\t#\tResolution" << endl;
	for ( i=0; i<loc.size(); ++i ) if ( !sel[i] ) {
		ns++;
		sel[i] = ns;
		rloc = sloc = loc[i];
		sloc /= p->real_size();
		sf = sloc.length();
		res_avg += 1/sf;
		res_std += 1/(sf*sf);
		cout << setprecision(2) << loc[i] << tab << atan2(sloc[1], sloc[0])*180.0/M_PI << tab <<
			R[i] << tab << sel[i] << tab << setprecision(4) << 1/sf << endl;
		for ( n=1; n<sym; ++n ) {
			rloc = mat*rloc;
			for ( j=0; j<loc.size(); ++j ) if ( !sel[j] ) {
				sloc = loc[j];
				if ( sloc.distance(rloc) < ks ) {
					sel[j] = sel[i];
					sloc /= p->real_size();
					sf = sloc.length();
					res_avg += 1/sf;
					res_std += 1/(sf*sf);
					if ( verbose )
						cout << setprecision(2) << loc[j] << tab << atan2(sloc[1], sloc[0])*180.0/M_PI << tab <<
							R[j] << tab << sel[j] << tab << setprecision(4) << 1/sf << endl;
				}
			}
		}
		if ( verbose )
			cout << endl;
	}
	
	res_avg /= loc.size();
	res_std = sqrt(res_std/loc.size() - res_avg*res_avg);

	if ( verbose ) {
		cout << "Resolution average:             " << res_avg << tab << "(" << res_std << ")" << endl;
		cout << "Relative magnification:         " << res_avg/ref_res << endl << endl;
	}

	fit_ellipse_iter(loc);

	return 0;
}
