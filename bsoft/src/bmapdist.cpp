/**
@file	bmapdist.cpp
@brief	Calculates a distance criterion between a set of 2D images or 3D maps
@author Bernard Heymann
@date	Created: 20010308
@date	Modified: 20210508
**/

#include "rwimg.h"
#include "rwmg.h"
#include "mg_img_proc.h"
#include "Matrix.h"
#include "cluster.h"
#include "linked_list.h"
#include "file_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
Matrix		img_map_distance(Bstring* file_list, int type, double* cutoff, 
				double hires, double lores, Bimage* pmask);
int			project_set_selection(Bproject* project, long n, vector<long> idx);
int			project_rec_select(Bproject* project, vector<long> idx);

Bstring*	project_get_filenames(Bproject* project, int flag)
{
	Bstring*			file_list = NULL;
	Bparticle*			part = project->class_avg;
	Breconstruction*	rec;
	
	if ( flag ) {
		if ( part )
//			for ( part = project->class_avg; part; part = part->next )
				if ( part->fpart.length() > 3 ) string_add(&file_list, part->fpart);
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			if ( rec->frec.length() > 3 ) string_add(&file_list, rec->frec);
	}
	
	return file_list;
}


// Usage assistance
const char* use[] = {
" ",
"Usage: bmapdist [options] input.star/input.img [input.star/input.img]",
"---------------------------------------------------------------------",
"Calculates distance criterions between a set of 2D images or 3D maps.",
"CC: Correlation coefficient.",
"R: R factor (Euclidean distance).",
"FSC: Resolution based on Fourier shell correlation.",
"DPR: Resolution based on differential phase residual.",
"Clustering: The preference value should be set to the approximate average in the",
"	similarity matrix. The damping factor should be ~0.5, or between 0.1 and 0.9.",
" ",
"Actions:",
"-type fsc                Output type of criterion: cc (default), r, fsc, dpr.",
"-preference -10          Homogenous preference value for clustering.",
"-TwoD                    Compare 2D images rather than 3D maps.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-resolution 15.3,100     Resolution limits (angstrom).",
"-cutoff 0.5,45           Resolution cutoffs for FSC and DPR.",
"-lambda 0.6              Damping factor for clustering (default 0.5).",
" ",
"Input:",
"-mask mask.tif           Real space mask (same size as input maps).",
"-Matrix mat.dat          Input comparison matrix file.",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-matrix mat.dat          Output comparison matrix file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	int				type(0);				// Default type is CC
	int				flag(0);				// Flag to indicate 2D images
	double			hires(0.01);			// Default resolution is all-inclusive
	double			lores(1e37);			// Default resolution is all-inclusive
	double			cutoff[4] = {0.5,45,1,0.5};	// Cutoff values
	double			pref(-1e37);			// Preference value for clustering
	double			lambda(0.5);			// Damping factor for clustering
	Bstring			real_mask_file;			// Real space mask file
	Bstring			outfile;				// Output parameter file
	Bstring			matin;					// Input comparison matrix file
	Bstring			matout;					// Output comparison matrix file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "type" ) {
			if ( curropt->value[0] == 'r' ) type = 1;
			if ( curropt->value[0] == 'f' ) type = 2;
			if ( curropt->value[0] == 'd' ) type = 3;
		}
		if ( curropt->tag == "preference" )
			if ( ( pref = curropt->value.real() ) < -1e36 )
				cerr << "-preference: A preference value must be specified!" << endl;
		if ( curropt->tag == "TwoD" ) flag = 1;
		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "cutoff" )
			if ( curropt->values(cutoff[0], cutoff[1], cutoff[2], cutoff[3]) < 1 )
				cerr << "-cutoff: A cuttoff must be specified!" << endl;
		if ( curropt->tag == "lambda" )
			if ( ( lambda = curropt->value.real() ) < 0.001 )
				cerr << "-lambda: A damping factor must be specified!" << endl;
		if ( curropt->tag == "mask" )
			real_mask_file = curropt->filename();
		if ( curropt->tag == "Matrix" )
			matin = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "matrix" )
			matout = curropt->filename();
    }
	option_kill(option);
	
 	double		ti = timer_start();
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No image or parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = NULL;

	if ( file_type(*file_list) == Micrograph ) {
		project = read_project(file_list);
		string_kill(file_list);
		file_list = project_get_filenames(project, flag);
	}

	if ( !file_list ) {
		cerr << "Error: No image file names found!" << endl;
		bexit(-1);
	}

	Bimage*			pmask = NULL;
	if ( real_mask_file.length() )
		pmask = read_img(real_mask_file, 1, 0);
	
	Matrix			matrix;
	
	if ( matin.length() )
		matrix = Matrix(matin);
	else
		matrix = img_map_distance(file_list, type, cutoff, hires, lores, pmask);

	if ( matout.length() )
		matrix.write(matout);

	long			i, ncluster(0);
	vector<long>	idx;
	if ( matrix.rows() ) {
		if ( pref > -1e30 ) {
			for ( i=0; i<matrix.rows(); i++ )
				matrix[i][i] = pref;
			idx = affin_prop_clustering(matrix, 500, 50, lambda, ncluster);
			if ( project ) {
				if ( flag )
					project_set_selection(project, ncluster, idx);
				else
					project_rec_select(project, idx);
			}
		}
	

	}
	
	// Write an output parameter file if a name is given
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}

	project_kill(project);
	delete pmask;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

Matrix		img_map_distance(Bstring* file_list, int type, double* cutoff, 
				double hires, double lores, Bimage* pmask)
{
	if ( type < 0 || type > 3 ) type = 0;
	
	if ( hires < 0.01 ) hires = 0.01;
	
	// Main loop over all reference maps
	long				i, j, ii, jj, n1, n2, nmap(0);
	Bstring*			fn1;
	Bstring*			fn2;
	Bimage*				p1 = NULL;
	Bimage*				p2 = NULL;

	for ( fn1 = file_list; fn1; fn1 = fn1->next ) {
		p1 = read_img(*fn1, 0, -1);
		nmap += p1->images();
		delete p1;
	}
	
	if ( verbose ) {
		cout << "Calculating the similarities and differences between " << nmap << " images:" << endl;
		cout << "Resolution limit:               " << hires << endl;
		if ( pmask ) 
			cout << "Mask:                           " << pmask->file_name() << endl;
	}
	
	Matrix				matrix;

	if ( nmap < 2 ) return matrix;
	
	Bimage*				pc = NULL;
	vector<double>		CC(nmap*nmap);
	vector<double>		R(nmap*nmap);
	vector<double>		FSC(nmap*nmap);
	vector<double>		DPR(nmap*nmap);
	Bplot*				plot;
	
	cout << "Map1\tMap2\tCC\tR\tFSC\tDPR" << endl;
	for ( i=0, fn1 = file_list; fn1; fn1 = fn1->next ) {
		p1 = read_img(*fn1, 0, -1);
		if ( !p1 ) bexit(-1);
		n1 = p1->images();
		delete p1;
		for ( ii=0; ii<n1; ++ii, ++i ) {
			p1 = read_img(*fn1, 1, ii);
			p1->change_type(Float);
			p1->rescale_to_avg_std(0, 1);
			if ( pmask ) p1->multiply(pmask);
			for ( j=0, fn2 = file_list; j<i && fn2; fn2 = fn2->next ) {
				p2 = read_img(*fn2, 0, -1);
				if ( !p2 ) bexit(-1);
				n2 = p2->images();
				delete p2;
				for ( jj=0; j<i && jj<n2; ++jj, ++j ) {
					p2 = read_img(*fn2, 1, jj);
					p2->change_type(Float);
					p2->rescale_to_avg_std(0, 1);
					if ( pmask ) p2->multiply(pmask);
					CC[i*nmap+j] = CC[j*nmap+i] = p1->correlate(p2);
					R[i*nmap+j] = R[j*nmap+i] = p1->R_factor(p2);
					pc = p1->resolution_prepare(p2);
					delete p2;
					plot = pc->fsc_dpr(hires, 1);
					delete pc;
					FSC[i*nmap+j] = FSC[j*nmap+i] = 1/plot->cut(1, cutoff[0], -1);
					DPR[i*nmap+j] = DPR[j*nmap+i] = 1/plot->cut(2, cutoff[1], 1);
					delete plot;
					cout << *fn1 << "(" << i << ")" << tab << *fn2 << "(" << j << ")" << tab <<
							CC[i*nmap+j] << tab << R[i*nmap+j] << tab << 
							FSC[i*nmap+j] << tab << DPR[i*nmap+j] << endl;
				}
			}
			delete p1;
		}
	}
	
	// All CC values
	cout << endl << "Correlation coefficients:" << endl;
	for ( i=0; i<nmap; i++ ) {
		CC[i*nmap+i] = 1;
		cout << tab << i+1;
	}
	cout << endl;
	for ( i=0; i<nmap; i++ ) {
		cout << i+1;
		for ( j=0; j<nmap; j++ )
			cout << tab << CC[i*nmap+j];
		cout << endl;
	}
	cout << endl;
	
	// The negative logarithm of the CC values
	cout << "-log(Correlation coefficients):" << endl;
	for ( i=0; i<nmap; i++ )
		cout << tab << i+1;
	cout << endl;
	for ( i=0; i<nmap; i++ ) {
		cout << i+1;
		for ( j=0; j<nmap; j++ )
			cout << tab << -log(CC[i*nmap+j]);
		cout << endl;
	}
	cout << endl;
		
	// All R values
	cout << "R factors:" << endl;
	for ( i=0; i<nmap; i++ ) {
		R[i*nmap+i] = 0;
		cout << tab << i+1;
	}
	cout << endl;
	for ( i=0; i<nmap; i++ ) {
		cout << i+1;
		for ( j=0; j<nmap; j++ )
			cout << tab << R[i*nmap+j];
		cout << endl;
	}
	cout << endl;
	
	// The FSC values
	cout << "Resolution based on Fourier shell correlation:" << endl;
	cout << "Cutoff:                         " << cutoff[0] << endl;
	for ( i=0; i<nmap; i++ ) {
		cout << tab << i+1;
		FSC[i*nmap+i] = hires;
	}
	cout << endl;
	for ( i=0; i<nmap; i++ ) {
		cout << i+1;
		for ( j=0; j<nmap; j++ )
			cout << tab << FSC[i*nmap+j];
		cout << endl;
	}
	cout << endl;
		
	// The DPR values
	cout << "Resolution based on differential phase residuals:" << endl;
	cout << "Cutoff:                         " << cutoff[1] << endl;
	for ( i=0; i<nmap; i++ ) {
		cout << tab << i+1;
		DPR[i*nmap+i] = hires;
	}
	cout << endl;
	for ( i=0; i<nmap; i++ ) {
		cout << i+1;
		for ( j=0; j<nmap; j++ )
			cout << tab << DPR[i*nmap+j];
		cout << endl;
	}
	cout << endl;

	matrix = Matrix(nmap, nmap);
	
	Bstring			label;
	for ( i=0; i<nmap; i++ ) {
		switch ( type ) {
			case 0:
				for ( j=0; j<nmap; j++ ) matrix[i][j] = CC[i*nmap+j];
				break;
			case 1:
				for ( j=0; j<nmap; j++ ) matrix[i][j] = R[i*nmap+j];
				break;
			case 2:
				for ( j=0; j<nmap; j++ ) matrix[i][j] = FSC[i*nmap+j];
				break;
			case 3:
				for ( j=0; j<nmap; j++ ) matrix[i][j] = DPR[i*nmap+j];
				break;
		}
	}
		
	return matrix;
}

int			project_set_selection(Bproject* project, long n, vector<long> idx)
{
	long			i, j;
	Bimage*			p = NULL;

	if ( project->class_avg && project->class_avg->fpart.length() )
		p = read_img(project->class_avg->fpart, 1, -1);
	else {
		cerr << "Error: No class averages!" << endl;
		return -1;
	}
	
	Bstring			filename(project->class_avg->fpart);
	filename = filename.pre_rev('.') + "_cluster." + filename.post_rev('.');
	
	if ( verbose ) {
		cout << "Replacing class averages:" << endl;
		cout << "Old number:                      " << p->images() << endl;
		cout << "New number:                      " << n << endl;
		cout << "New file name:                   " << filename << endl;
	}

	Bfield*			field = project->field;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( part = mg->part; part; part = part->next )
				if ( part->sel > 0 )
					part->sel = idx[part->sel-1] + 1;
	
	Bimage*			pnu = p->copy_header(n);
	pnu->data_alloc_and_clear();
	Bimage*			p1;
	Bparticle*		partlist = NULL;
	Bparticle*		nupart = NULL;
	
	for ( i=j=0, part = project->class_avg; i < p->images() && part; ++i, part = part->next ) {
		if ( i == idx[i] ) {
			nupart = particle_copy(&partlist, part);
			nupart->fpart = filename;
			nupart->id = ++j;
			nupart->group = i;
			nupart->sel = 0;
		}
	}
	
	for ( i=0, part = project->class_avg; i < p->images(); ++i, part = part->next ) {
		for ( j=0, nupart = partlist; nupart; ++j, nupart = nupart->next )
			if ( nupart->group == idx[i] ) break;
		nupart->sel += part->sel;
		p1 = p->extract(i);
		pnu->add(j, p1);
		delete p1;
	}
	
	pnu->statistics();
	
	for ( j=0; j<pnu->images(); ++j )
		pnu->rescale_to_avg_std(j, 0, 1);
	
	particle_kill(project->class_avg);
	project->class_avg = partlist;
	
	write_img(filename, pnu, 0);
	
	delete pnu;
	delete p;
	
	return 0;
}


int			project_rec_select(Bproject* project, vector<long> idx)
{
	long				i;
	Breconstruction*	rec;
	
	for ( i=0, rec = project->rec; rec; rec = rec->next, i++ ) {
		rec->select = idx[i] + 1;
	}
	
	return 0;
}

