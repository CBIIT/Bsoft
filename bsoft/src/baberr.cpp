/**
@file	baberr.cpp
@brief	Program to generate images from aberration weights and analyze them
@author Bernard Heymann
@date	Created: 20220106
@date	Modified: 20221212
**/

#include "rwimg.h"
#include "rwPostScript.h"
#include "rwmg.h"
#include "mg_ctf.h"
#include "mg_ctf_fit.h"
#include "simplex.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			weights_add(map<pair<long,long>,double>& weights, vector<double> v, int flag);
int			weights_add(map<pair<long,long>,double>& weights, string s, int flag);
map<pair<long,long>,double>	convert_zernike_weights(map<pair<long,long>,double>& weights);
map<pair<long,long>,double>	convert_to_zernike_weights(map<pair<long,long>,double>& weights);
vector<Vector3<double>>	shifts_from_file(Bstring& shiftfile);
int			aberration_show(map<string,CTFparam>& cpa);
int			img_ewald_sphere(Bimage* p, CTFparam& cp, double t);

// Usage assistance
const char* use[] = {
" ",
"Usage: baberr [options] input.star",
"----------------------------------",
"Generates or analyzes aberration images.",
" ",
"Actions:",
"-create 200,200,1        Create a new image of this size (pixels/voxels).",
//"                         Note: No input image is read and the first file name is taken as the output image.",
"-basis 2                 Generate basis images (0=full, 1=odd, 2=even).",
"-zernike                 Treat weights as zernike polynomials (only with the -even and -odd options).",
"-fit 2,10000             Fit phase image; flag: 0=full, 1=odd, 2=even; fitting iterations.",
"-halfshift               Shift output image to center the origin.",
"-compare another.star    Compare aberration parameters with this parameter file.",
"-show                    Show aberration weights and statistics.",
"-ewald 125               Create an Ewald sphere weight map for the given thickness (Å).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-resolution 27.5,125.3   Resolution limits for aberration fitting (default 0.1,1e6).",
" ",
#include "use_ctf.inc"
"-Defocus 1.2             Defocus average (default 2 um).",
"-Astigmatism 0.3,-34     Set defocus deviation and astigmatism angle (um, degrees).",
" ",
"Parameters for creating an aberration image:",
"-weight 2,1,0.46         Indices n, m and weight (can be used repeatedly).",
"-odd [0.5,-0.3,...]      Odd weights (comma-separated list).",
"-even [4,3.2,...]        Even weights (comma-separated list).",
"-fromoptics 1            Create from optics group parameters, flag: 0=full, 1=odd, 2=even (default 0).",
" ",
"Input:",
"-phasedifference p.grd   Particle phase difference image.",
"-shifts shifts.txt       Image shifts.",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-aberration ab.grd       Created or fitted aberration image.",
"-difference dif.grd      Fit difference image.",
"-json weights.json       Output CTF parameters with aberration weights.",
"-Postscript img.ps       Output Postscript file with aberration images.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	Vector3<long>	nusize;					// New image size
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;					// Units for the three axes (A/pixel)
	int				basis(-1);				// Flag to generate basis images
	int				zflag(0);				// Flag to indicate zernike polynomial weights
	int				from_optics(-1);		// Flag to generate aberration images from optics groups
	int				fit_flag(-1);			// Fit type: 0=full, 1=odd, 2=even
	double			hi_res(0), lo_res(0);	// Fit resolution limits
	long			fit_iter(0);			// Fit iterations
	int				halfshift(0);			// Output image shift
	long			m(0), n(0);				// Weight indices
	double			v(0);					// Weight value
	map<pair<long,long>,double>	weights;	// Aberration weights
	string			wodd, weven;			// Strings for weight lists
    map<string,CTFparam>	cpa;			// CTF parameters
    double			ew_thick(0);			// Ewald sphere thickness
    int				show(0);				// Flag to show aberration weights
	Bstring			compfile;				// Input comparison parameter file name
	Bstring			phifile;				// Input correlation image file name
	Bstring			shiftfile;				// Input image shifts file name
	Bstring			outfile;				// Output parameter file name
	Bstring			abfile;					// Output aberration image file name
	Bstring			diffile;				// Output fit difference image file name
	Bstring			jsfile;					// Output fit difference image file name
	Bstring			psfile;					// Output fit difference image file name
	int				write_flags(0);			// Flags to pass to the parameter file writing function

//	double			v;
	JSvalue			jsctf(JSobject);
	double			def_avg(0);				// In angstrom
	double			def_dev(0);				// In angstrom
	double			ast_angle(0);	 		// Used to limit astigmatism

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "create" )
			nusize = curropt->size();
		if ( curropt->tag == "basis" )
			if ( ( basis = curropt->value.integer() ) < 0 )
				cerr << "-basis: A flag value must be specified!" << endl;
		if ( curropt->tag == "fromoptics" )
			if ( ( from_optics = curropt->value.integer() ) > 2 )
				cerr << "-fromoptics: A value of 0, 1 or 2 must be specified!" << endl;
		if ( curropt->tag == "zernike" )
			zflag = 2;
		if ( curropt->tag == "fit" )
    	    if ( curropt->values(fit_flag, fit_iter) < 1 )
				cerr << "-fit: A fitting flag value must be specified!" << endl;
		if ( curropt->tag == "halfshift" ) halfshift = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "resolution" )
			if ( curropt->values(hi_res, lo_res) < 1 )
				cerr << "-resolution: At least one resolution limit must be specified!" << endl;
		if ( curropt->tag == "weight" ) {
    	    if ( curropt->values(n, m, v) < 3 )
				cerr << "-weight: Three values must be specified!" << endl;
			else
				weights[{n,m}] = v;
		}
		if ( curropt->tag == "odd" ) wodd = curropt->value.str();
		if ( curropt->tag == "even" ) weven = curropt->value.str();
		if ( curropt->tag == "ewald" )
 			if ( ( ew_thick = curropt->value.real() ) < 1 )
				cerr << "-ewald: A thickness must be specified!" << endl;
#include "ctf.inc"
		if ( curropt->tag == "Defocus" )
			def_avg = curropt->real_units();
		if ( curropt->tag == "Astigmatism" ) {
			if ( curropt->values(def_dev, ast_angle) < 1 )
				cerr << "-Astigmatism: A defocus value must be specified!" << endl;
			else {
				if ( def_dev < 1e3 ) def_dev *= 1e4;			// Assume um
				ast_angle *= M_PI/180.0;						// Assume degrees
			}
		}
 		if ( curropt->tag == "show" ) show = 1;
		if ( curropt->tag == "compare" )
			compfile = curropt->filename();
		if ( curropt->tag == "phasedifference" )
			phifile = curropt->filename();
		if ( curropt->tag == "shifts" )
			shiftfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "aberration" )
			abfile = curropt->filename();
		if ( curropt->tag == "difference" )
			diffile = curropt->filename();
		if ( curropt->tag == "json" )
			jsfile = curropt->filename();
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
	}
	option_kill(option);
	
	double			ti = timer_start();

	CTFparam		cp = ctf_from_json(jsctf);
	cp.aberration_weights(weights);
	if ( def_avg ) cp.defocus_average(def_avg);
	if ( def_dev ) cp.astigmatism(def_dev, ast_angle);
	
//	cout << jsctf << endl;
	cp.show();

	if ( nudatatype == Unknown_Type ) nudatatype = Float;

	Bproject*		project = NULL;
    Bimage* 		p = NULL;
    
	if ( nusize.volume() > 0 ) {
		long		n(1);
		if ( basis >= 0 ) {
			if ( basis == 0 ) n = 15;
			if ( basis == 1 ) n = 6;
			if ( basis == 2 ) n = 9;
		}
		p = new Bimage(nudatatype, TSimple, nusize, n);
		if ( sam.volume() > 0 ) p->sampling(sam);
		if ( basis >= 0 ) {
			img_aberration_basis(p, basis);
		} else if ( ew_thick > 0 ) {
			if ( cp.volt() < 1 ) {
				cerr << "Error: The acceleration voltage must be specified!" << endl;
				bexit(-1);
			}
//			p->ewald_sphere(volt, ew_thick);
			img_ewald_sphere(p, cp, ew_thick);
		} else {
			if ( wodd.length() || weven.length() ) {
				if ( wodd.length() ) weights_add(weights, wodd, zflag|1);
				if ( weven.length() ) weights_add(weights, weven, zflag);
			} else {
				weights = cp.aberration_weights();
			}
			img_create_aberration(p, weights, 0);
		}
	} else if ( phifile.length() ) {
		p = read_img(phifile, 1, -1);
		if ( p == NULL ) bexit(-1);
		if ( sam.volume() > 0 ) p->sampling(sam);
	}
	
	if ( argc > optind ) {
		Bstring				filename = argv[optind++];
		project = read_project(filename);
//		if ( volt ) project_set_volts(project, volt);
		project_update_ctf(project, jsctf);
		cpa = project_ctf_optics_groups(project);
		if ( p && p->images() != cpa.size() ) {
			cerr << "Error: Mismatch in the number of optics groups and correlation images!" << endl;
			bexit(-1);
		}
		if ( show )
			aberration_show(cpa);
//			for ( auto cp: cpa ) cp.second.show_aberration();
	}
	
	if ( from_optics>=0 && project ) {
		if ( p ) delete p;
		Bmicrograph*	mg = project->field->mg;
		nusize = mg->box_size;
		p = new Bimage(nudatatype, TSimple, nusize, cpa.size());
		p->sampling(mg->part->pixel_size);
		if ( sam.volume() > 0 ) p->sampling(sam);
		img_create_aberration(p, cpa, from_optics);
	}
	
	if ( compfile.length() ) {
		Bproject*	project2 = read_project(compfile);
		project_aberration_compare(project, project2);
		project_kill(project);
	}

	vector<map<pair<long,long>,double>> wa;
	if ( weights.find({0,0}) == weights.end() )
		weights[{0,0}] = cp.aberration_weight(0,0);
	if ( weights.find({2,0}) == weights.end() )
		weights[{2,0}] = cp.aberration_weight(2,0);
	if ( weights.find({4,0}) == weights.end() )
		weights[{4,0}] = cp.aberration_weight(4,0);
	if ( weights.size() ) wa.push_back(weights);
	cout << "Checking weights:" << endl;
	for ( auto w: weights ) cout << w.first.first << tab << w.first.second << tab << w.second << endl;

	Bimage*		pd = NULL;
	if ( fit_flag >= 0 && p ) {
		if ( fit_flag > 3 )
			img_ctf_fit_prepare(p, 10);
		wa = img_aberration_phase_fit(p, lo_res, hi_res, wa, fit_flag, fit_iter);
		if ( abfile.length() ) {
			Bimage*		pf = p->copy();
			img_create_aberration(pf, wa, fit_flag);
			if ( diffile.length() ) {
				p->subtract(pf);
				pd = p;
			} else {
				delete p;
			}
			p = pf;
		}
		if ( cpa.size() < 1 ) {
			for ( long i=1; i<=wa.size(); ++i ) {
				cpa[to_string(i+50)]=CTFparam(to_string(i+50));
			}
		}
		if ( wa.size() == cpa.size() ) {
			if ( verbose )
				cout << "Transferring aberration weights:" << endl;
			long		i(0);
			for ( auto& cp: cpa ) {
//				if ( volt > 0 ) cp.second.volt(volt);
				ctf_update_from_json(cp.second, jsctf);
				cp.second.aberration_weights(wa[i++]);
				if ( verbose )
					cout << i << tab << cp.second.aberration_weight_string() << endl;
			}
		} else {
			cerr << "Error: The number of weights records (" << wa.size() <<
				") does not agree with the number of optics groups (" << cpa.size() << ")!" << endl;
		}
		if ( project )
			project_update_ctf_aberration(project, cpa, 1);
		if ( jsfile.length() ) {
			JSvalue		js(JSarray);
			for ( auto cp: cpa )
				js.push_back(ctf_to_json(cp.second));
			js.write(jsfile.str());
		}
	}
	
	if ( halfshift ) {
		p->shift_wrap(-p->size()/2);
		p->origin(p->size()/2);
		if ( pd ) {
			pd->shift_wrap(-p->size()/2);
			pd->origin(p->size()/2);
		}
	}
	
	if ( outfile.length() && project )
		write_project(outfile, project, write_flags);
	
	if ( p && abfile.length() ) {
		if ( phifile.length() ) img_create_aberration(p, wa, 4);
		p->change_type(nudatatype);
		write_img(abfile, p, 0);
	}

	if ( pd && diffile.length() ) {
		pd->change_type(nudatatype);
		write_img(diffile, pd, 0);
	}
	
	if ( shiftfile.length() && psfile.length() ) {
		vector<Vector3<double>>	shifts = shifts_from_file(shiftfile);
		p->file_name(psfile.str());
		writePostScriptImage(p, shifts);
	}

	delete p;
	delete pd;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/*
	flag	0=even, 1=odd
*/
int			weights_add(map<pair<long,long>,double>& weights, vector<double> v, int flag)
{
	// Converting Zernike weights
	if ( flag == 2 ) {
		v[0] += -v[2] + v[6];
		v[1] -= 3*v[5];
		v[2] = 2*v[2] - 6*v[6];
		v[3] -= 3*v[7];
		v[5] *= 4;
		v[6] *= 6;
		v[7] *= 4;
	}
	if ( flag == 3 ) {
		v[0] -= v[3];
		v[1] -= v[4];
		v[3] *= 3;
		v[4] *= 3;
	}

	if ( verbose & VERB_FULL ) {
		cout << "Adding weights" << endl;
		cout << "#\tn\tm\tw" << endl;
	}
	
	long					i, n(flag&1), m, t(2);
	if ( v.size() > 9 ) t = 1;
		
	for ( i=0; n<=10 && i<v.size(); n+=t ) {
		for ( m=-n; m<=n; m+=2 ) {
//			if ((n - m) % 2 == 0) {
				weights[{n,m}] = v[i++];
				if ( verbose & VERB_FULL )
					cout << i << "\t" << n << "\t" << m << "\t" << weights[{n,m}] << endl;
//			}
		}
	}
	
	return 0;
}

int			weights_add(map<pair<long,long>,double>& weights, string s, int flag)
{
	if ( s[0] == '[' ) s = s.substr(1);
	
	vector<double>	v = parse_real_vector(s);

	return weights_add(weights, v, flag);
}

map<pair<long,long>,double>	convert_zernike_weights(map<pair<long,long>,double>& weights)
{
	map<pair<long,long>,double>	w = weights;
	
	w[{0,0}] += -w[{2,0}] + w[{4,0}];
	w[{1,-1}] -= w[{3,-1}];
	w[{1,1}] -= w[{3,1}];
	w[{2,-2}] -= 3*w[{4,-2}];
	w[{2,0}] = 2*w[{2,0}] - 6*w[{4,0}];
	w[{2,2}] -= 3*w[{4,2}];
	w[{3,-1}] *= 3;
	w[{3,1}] *= 3;
	w[{4,-2}] *= 4;
	w[{4,0}] *= 6;
	w[{4,2}] *= 4;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Simple weights:" << endl;
		cout << "n\tm\tw" << endl;
		for ( auto w1: w )
			cout << w1.first.first << tab << w1.first.second << tab << w1.second << endl;
	}
	
	return w;
}

map<pair<long,long>,double>	convert_to_zernike_weights(map<pair<long,long>,double>& weights)
{
	map<pair<long,long>,double>	w = weights;
	
	w[{0,0}] += w[{2,0}]/2 + w[{4,0}]/3;
	w[{1,-1}] += w[{3,-1}]/3;
	w[{1,1}] += w[{3,1}]/3;
	w[{2,-2}] += 0.75*w[{4,-2}];
	w[{2,0}] = (w[{2,0}] + w[{4,0}])/2;
	w[{2,2}] += 0.75*w[{4,2}];
	w[{3,-1}] /= 3;
	w[{3,1}] /= 3;
	w[{4,-2}] /= 4;
	w[{4,0}] /= 6;
	w[{4,2}] /= 4;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Zernike weights:" << endl;
		cout << "n\tm\tw" << endl;
		for ( auto w1: w )
			cout << w1.first.first << tab << w1.first.second << tab << w1.second << endl;
	}
	
	return w;
}

vector<Vector3<double>>	shifts_from_file(Bstring& shiftfile)
{
	Vector3<double>				v;
	vector<Vector3<double>>		t;
	
    ifstream		ft(shiftfile.str());
    if ( ft.fail() ) return t;
    
    cout << "Reading " << shiftfile << endl;

 	while ( !ft.eof() ) {
 		ft >> v[0] >> v[1];
 		if ( !ft.eof() ) {
 			t.push_back(v);
 			cout << v << endl;
		}
	}
	
	cout << "Shifts read: " << t.size() << endl;
	
	ft.close();
	
	return t;
}

int			aberration_show(map<string,CTFparam>& cpa)
{
	long						n(cpa.size());
	map<pair<long,long>,double>	wab, wavg, wstd;
	
	for ( auto cp: cpa ) {
		wab = cp.second.aberration_weights();
		for ( auto w: wab ) {
			wavg[w.first] += w.second;
			wstd[w.first] += w.second*w.second;
		}
		cp.second.show_aberration();
	}

	cout << "Statistics:" << endl;
	cout << "n\tm\tavg\tstd" << endl;
	for ( auto& w: wavg ) {
		w.second /= n;
		wstd[w.first] = sqrt(wstd[w.first]/n - w.second*w.second);
		cout << w.first.first << tab << w.first.second << tab << w.second << tab << wstd[w.first] << endl;
	}

	
	return 0;
}

/**
@brief 	Imposes an Ewald sphere weighting on a 3D frequency space volume.
@param 	*p				frequency space 3D volume.
@param	cp				wavelength.
@param	t				thickness.
@return int				0.
**/
int			img_ewald_sphere(Bimage* p, CTFparam& cp, double t)
{
	double			wl = cp.lambda();
	
	long			i, nn, xx, yy, zz;
	double			u, v, w, wew, s2, a, theta, phi, amp;
	Vector3<long>	h(p->size()/2);
	Vector3<double>	fspace_scale(1/p->real_size());
	
	if ( verbose ) {
		cout << "Generating an Ewald sphere weight map:" << endl;
		cout << "Thickness:                     " << t << " A" << endl;
		cout << "Frequency sampling:            " << fspace_scale << " 1/A" << endl;
		cp.show();
	}
	
	for ( i=nn=0; nn<p->images(); nn++ ) {
		for ( zz=0; zz<p->sizeZ(); ++zz ) {
			w = (zz < h[2])? zz: zz-p->sizeZ();
			w *= fspace_scale[2];
			w = fabs(w);
			for ( yy=0; yy<p->sizeY(); ++yy ) {
				v = (yy < h[1])? yy: yy-p->sizeY();
				v *= fspace_scale[1];
				for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
					u = (xx < h[0])? xx: xx-p->sizeX();
					u *= fspace_scale[0];
					s2 = u*u + v*v;
					a = atan2(v, u);
					theta = atan2(sqrt(s2), 1.0/wl-w);
					wew = (1.0/wl)*(1.0-cos(theta));
					phi = M_PI*t*(w-wew);
					if ( phi ) amp = sin(phi)/phi;
					else amp = 1;
					amp *= -cp.calculate(s2, a);
					p->set(i, amp*amp);
				}
			}
		}
	}
	
	return 0;
}
