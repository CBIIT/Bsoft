/**
@file	btile.cpp
@brief	Program to split an image into overlapping tiles
@author Bernard Heymann
@date	Created: 20040712
@date	Modified: 20160601
**/

#include "rwmg.h"
#include "mg_processing.h"
#include "rwimg.h"
#include "file_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			img_tile_with_overlap(Bimage* p, Bstring outputfile, Vector3<int> size, 
				Vector3<int> overlap, int digits, int first_number, DataType nudatatype);

// Usage assistance
const char* use[] = {
" ",
"Usage: btile [options] input.img/input.star output.img",
"------------------------------------------------------",
"Splits an image into overlapping tiles.",
"During tiling, the original size and tile origins are written into a text",
"   file with the base of the output file name and the extension \".tiles\".",
"The program bpatch can reassemble the tiles using this text file.",
" ",
"Action:",
"-split                   Split into multiple images.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-size 120,150,70         Size of tiles (default 512,512,512).",
"-overlap 13,15,8         Overlap of tiles (default 0,0,0).",
"-first 5                 Number given to the first file (default 0).",
"-digits 3                Number of digits inserted before the last period in the output file name (default 3).",
" ",
"Parameters for reconstructions:",
"-recid rec_1             Reconstruction ID to use.",
" ",
"Output:",
"-output output.star      Output parameter file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	int				split(0);					// Flag to split into individual images
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double> sampling;					// Units for the three axes (A/pixel)
	Vector3<long>	size(512,512,512);			// Tile size
	Vector3<long>	overlap;					// Tile overlap
	int 			first_number(0);			// Number given to first file
	int 			digits(3);					// File number size
	Bstring			recid;						// Reconstruction ID to tile
	Bstring			outfile;					// Output micrograph parameter file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "split" ) split = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sampling = curropt->scale();
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "overlap" )
			overlap = curropt->size();
		if ( curropt->tag == "first" ) {
			if ( ( first_number = curropt->value.integer() ) < 0 )
				cerr << "-first: A number must be specified!" << endl;
			if ( first_number < 0 ) first_number = 0;
		}
		if ( curropt->tag == "digits" )
			if ( ( digits = curropt->value.integer() ) < 1 )
				cerr << "-digits: A number of digits must be specified!" << endl;
		if ( curropt->tag == "recid" ) {
			recid = curropt->value;
			if ( recid.length() < 1 )
				cerr << "-recid: A reconstruction ID must be specified!" << endl;
		}
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bstring				filename = argv[optind++];
	Bproject*			project = NULL;
	Breconstruction*	rec;
	
	if ( file_type(filename) == Micrograph ) {
		project = read_project(filename);
		project->select = 1;
		for ( rec = project->rec; rec && rec->id != recid; rec = rec->next ) ;
		if ( !rec ) {
			cerr << "Error: Reconstruction " << recid << " not found!" << endl;
			bexit(-1);
		}
		filename = rec->frec;
		if ( sampling[0] > 0 ) rec->voxel_size = sampling;
	}

	Bimage*		p = read_img(filename, 1, -1);
	if ( p == NULL ) bexit(-1);
    
	if ( optind >= argc ) bexit(0);
	
	if ( sampling.volume() > 0 ) p->sampling(sampling);

	
	if ( argc > optind ) {
		Bstring					filename(argv[optind]);
		if ( p->images() > 1 ) {
			Bstring					fn1;
			char					format[16];
			snprintf(format, 16, "%%0%dd.", digits);
			Vector3<long>			step(size-overlap);
			vector<Vector3<long>>	coor = p->tile_coordinates(size, step);
			Bimage**	parr = p->extract_tile_stacks(coor, size);
			for ( long n=0; n<coor.size(); ++n ) {
				fn1 = filename.pre_rev('.') + Bstring(n+first_number, format) + filename.post_rev('.');
				parr[n]->change_type(nudatatype);
				write_img(fn1, parr[n], 0);
			}
		} else if ( split ) {
			img_tile_with_overlap(p, argv[optind], size, overlap, digits, first_number, nudatatype);
		} else {
			Vector3<long>	start, step(size-overlap);
			Bimage*			ptile = p->extract_tiles(0, start, p->size(), size, step, 1);
			ptile->change_type(nudatatype);
			write_img(argv[optind], ptile, 0);
			delete ptile;
		}
	}
	
    delete p;
	
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/**
@brief 	Writes overlapping tiles into new images.
@param 	*p				image = tile.
@param 	outputfile		output file name.
@param 	size			tile size.
@param 	overlap			tile overlap.
@param 	digits			number of digits in tile file number.
@param 	first_number	first tile file number.
@param 	nudatatype		tile datatype.
@return int				number of tiles.

	Each tile is copied into a new image and written with a numbered file name.
	The tiles at the ends of the block may not be the same size as the
	other tiles, being truncated at the end of the input data.
	A text file with the extension ".tiles" is written with the positions
	of the tiles, and it is needed to patch the image together again.

**/
int			img_tile_with_overlap(Bimage* p, Bstring outputfile, Vector3<int> size, 
				Vector3<int> overlap, int digits, int first_number, DataType nudatatype)
{
	Bstring			filename;
	char			format[16];
    Bimage*			pone = NULL;
	Vector3<double>	start;
	
	if ( digits < 2 ) digits = 2;
	
	snprintf(format, 16, "%%0%dd.", digits);
	
	if ( size[0] < 1 ) size[0] = 512;
	if ( size[1] < 1 ) size[1] = 512;
	if ( size[2] < 1 ) size[2] = 512;
	if ( size[0] > p->sizeX() ) size[0] = p->sizeX();
	if ( size[1] > p->sizeY() ) size[1] = p->sizeY();
	if ( size[2] > p->sizeZ() ) size[2] = p->sizeZ();
	if ( overlap[0] < 0 ) overlap[0] = 0;
	if ( overlap[1] < 0 ) overlap[1] = 0;
	if ( overlap[2] < 0 ) overlap[2] = 0;
	if ( overlap[0] > p->sizeX() ) overlap[0] = p->sizeX()/2;
	if ( overlap[1] > p->sizeY() ) overlap[1] = p->sizeY()/2;
	if ( overlap[2] > p->sizeZ() ) overlap[2] = p->sizeZ()/2;
	
	int				n, i = first_number;
	Vector3<long>	step = size - overlap;
	Vector3<long>	xsize;
	
	filename = outputfile.pre_rev('.') + ".tiles";
	n = (p->sizeX()/step[0] + 1)*(p->sizeY()/step[1] + 1)*(p->sizeZ()/step[2] + 1);
	
	if ( verbose ) {
		cout << "Splitting image " << p->file_name() << " into tiles:" << endl;
		cout << "Tile size:                      " << size << endl;
		cout << "Overlap:                        " << overlap << endl;
		cout << "Number of tiles:                " << n << endl << endl;
	}

	ofstream	fd(filename.c_str());
	if ( fd.fail() ) return 0;
	
	fd << p->sizeX() << " " << p->sizeY() << " " << p->sizeZ() << endl;
	fd << overlap[0] << " " << overlap[1] << " " << overlap[2] << endl;

	for ( n=0, start[2]=0; start[2]<p->sizeZ(); start[2]+=step[2] ) {
		xsize[2] = p->sizeZ() - (long)start[2];
		if ( xsize[2] > size[2] ) xsize[2] = size[2];
		for ( start[1]=0; start[1]<p->sizeY(); start[1]+=step[1] ) {
			xsize[1] = p->sizeY() - (long)start[1];
			if ( xsize[1] > size[1] ) xsize[1] = size[1];
			for ( start[0]=0; start[0]<p->sizeX(); start[0]+=step[0], i++ ) {
				xsize[0] = p->sizeX() - (long)start[0];
				if ( xsize[0] > size[0] ) xsize[0] = size[0];
				pone = p->extract(0, start, xsize);
				pone->label(p->label());
				pone->change_type(nudatatype);
				filename = outputfile.pre_rev('.') + Bstring(i, format) + outputfile.post_rev('.');
				if ( verbose & VERB_LABEL )
					cout << "File " << i << ": " << filename << endl;
				write_img(filename, pone, 0);
				delete pone;
				fd << start[0] << " " << start[1] << " " << start[2] << endl;
				n++;
			}
		}
	}
	
	fd.close();
	
    if ( verbose & VERB_LABEL )
		cout << "Number of tiles:                " << n << endl << endl;
    
	return n;
}
