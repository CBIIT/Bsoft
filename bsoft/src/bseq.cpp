/**
@file	bseq.cpp
@brief	A program to manipulate DNA and protein sequences
@author Bernard Heymann
@date	Created: 20000808
@date	Modified: 20190613
**/

#include "rwmolecule.h"
#include "rwgencode.h"
#include "seq_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bseq [options] input.seq output.seq",
"------------------------------------------",
"Manipulates DNA and protein sequences.",
" ",
"Actions:",
"-show                    Show all DNA and protein sequences after other operations.",
"-Mass                    Show the calculated mass for each sequence.",
"-complement              Complement all DNA sequences before other operations.",
"-translate 1             Translate all DNA sequences with this frame (0,1,2).",
"-finddna ggattcga        DNA sequence to look for in DNA sequences.",
"-findprotein fghwereaas  Protein sequence to look for in protein sequences.",
"-findcoding aklsdrtv     Coding sequence to look for in DNA sequences.",
"                         Can be a file name with a protein sequence.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-nogaps                  Strip gaps from input sequences.",
"-length 1100,1500        Sequence length range to select for.",
"-number 10,14            Number of residues on either side to include in output.",
"-threshold 50            Threshold percentage to report hits (default only best hit).",
" ",
"Input:",
"-geneticcode file.star   Genetic code file (default = internal values).",
"-elements prop.star      Calculate elemental composition.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	int 			setshow(0);				// Show flag
	int 			setmass(0);				// Show flag
	int 			setfind(0);				// Search flag
	int 			setcomplement(0);		// Complement flag
	int 			settranslate(0);		// Translation flag
	int 			frame(0);				// Frame for translation
	int 			side1(0), side2(0);  	// Additional residues
	double			threshold(0);			// Only best match
	int 			seqlenmin(0);			// Minimum sequence length
	int 			seqlenmax(100000);		// Minimum sequence length
	Bstring			atom_select("all");		// Selection
	
	// Genetic code parameter file and template sequence file
	Bstring			gcfile;
	Bstring			seqfile;				// Sequence to search for in a file
	Bstring			rpfile;					// Residue properties file
	Bstring			paramfile;
    
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "nogaps" )
			atom_select = "nogap";
		if ( curropt->tag == "show" )
			setshow = 1;
		if ( curropt->tag == "Mass" )
			setmass = 1;
		if ( curropt->tag == "complement" )
			setcomplement = 1;
		if ( curropt->tag == "translate" ) {
        	 if ( ( frame = curropt->value.integer() ) < 0 )
				cerr << "-translate: A frame must be specified!" << endl;
			while ( frame < 0 ) frame += 3;
			while ( frame > 2 ) frame -= 3;
			settranslate = 1;
		}
		if ( curropt->tag == "length" ) {
        	if ( curropt->values(seqlenmin, seqlenmax) < 1 )
				cerr << "-length: Minimum and maximum lengths must be specified!" << endl;
			else {
				if ( seqlenmin < 0 ) seqlenmin = 0;
				if ( seqlenmax < seqlenmin ) {
					i = seqlenmin;
					seqlenmin = seqlenmax;
					seqlenmax = i;
				}
				if ( seqlenmax == seqlenmin ) seqlenmax = seqlenmin + 1;
			}
		}
		if ( curropt->tag == "finddna" ) {
			seqfile = curropt->filename();
			setfind = 1;
		}
		if ( curropt->tag == "findprotein" ) {
			seqfile = curropt->filename();
			setfind = 2;
		}
		if ( curropt->tag == "findcoding" ) {
			seqfile = curropt->filename();
			setfind = 3;
		}
		if ( curropt->tag == "Sequence" ) {
        	seqfile = curropt->filename();
			if ( seqfile.length() ) setfind = 3;
		}
		if ( curropt->tag == "number" )
			if ( curropt->values(side1, side2) < 1 )
				cerr << "-number: At least one number must be specified!" << endl;
		if ( curropt->tag == "threshold" )
        	if ( ( threshold = curropt->value.real() ) < 1e-30 )
				cerr << "-threshold: A threshold must be specified!" << endl;
		if ( curropt->tag == "geneticcode" )
            gcfile = curropt->filename();
		if ( curropt->tag == "elements" )
            rpfile = curropt->filename();
    }
	option_kill(option);
    
	double		ti = timer_start();
	
	// Read the sequence request file if given
	Bmolgroup*		seq_request = NULL;
	Bstring			seq;		// Sequence to search for
	if ( seqfile.length() ) {
		if ( seqfile.contains(".") ) {
			seq_request = read_molecule(seqfile, atom_select, paramfile);
			if ( seq_request->mol->naseq.length() ) {
				seq = seq_request->mol->naseq;
				seq_request->mol->naseq = NULL;
			} else {
				seq = seq_request->mol->seq;
				seq_request->mol->seq = NULL;
			}
			molgroup_kill(seq_request);
		} else {
			seq = seqfile;
		}
		// Convert threshold to number of residues
		threshold = floor(threshold*seq.length()/100.0);
	}
	
	if ( optind >= argc ) {
		cerr << "Error: No input file given!" << endl;
		bexit(-1);
	}
	
	if ( verbose && seq.length() ) cout << "Search string: " << seq << endl;
	
    // Read the sequence file
	seqfile = argv[optind++];
	Bmolgroup*		molgroup = read_molecule(seqfile, atom_select, paramfile);
	if ( !molgroup ) {
		cerr << "Error: No input file given!" << endl;
		bexit(-1);
	}
	
	if ( setcomplement )
		seq_complement_all(molgroup);
	
	if ( settranslate )
		seq_translate_all(molgroup, frame, gcfile);
		
	switch ( setfind ) {
		case 1: seq_find_dna(molgroup, seq); break;
		case 2: seq_find_protein(molgroup, seq); break;
		case 3: seq_find_protein_in_dna(molgroup, seq, seqlenmin, seqlenmax, side1, side2, threshold, gcfile); break;
		default: break;
	}

	if ( setshow )
		seq_show(molgroup);
	
	if ( setmass )
		seq_mass(molgroup);
	
	if ( rpfile.length() )
		seq_elements(molgroup, rpfile);
		
	molecule_update_comment(molgroup, argc, argv);
	
	if ( optind < argc )
		write_molecule(argv[optind], molgroup);
	
	molgroup_kill(molgroup);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

