/**
@file	ps_sequence.cpp
@brief	Postscript output for sequence analysis functions.
@author Bernard Heymann 
@date	Created: 20010515
@date	Modified: 20210426
**/
 
#include "ps_plot.h" 
#include "Complex.h"
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototype
int 		ps_seq_representation(ofstream* fps, vector<double>& nseq);
int 		ps_seq_periodicity(ofstream* fps, vector<Complex<float>>& per);

/**
@brief 	Generates a postscript plot from an information content analysis of a protein sequence alignment.
@param 	&filename	postscript file name.
@param 	&title 		title.
@param 	nres		number of residue types.
@param 	&info 		moving average of information content.
@param 	&nseq 		number of sequences represented at each position.
@param 	&per		periodicity analysis.
@param 	&freq 		frequence of each residue type at each position.
@param 	&pattern	residue types at each position, in frequency order.
@return int 			0.
**/
int 		ps_seq_info(Bstring& filename, Bstring& title, int nres,
				vector<double>& info, vector<double>& nseq,
				vector<Complex<float>>& per, vector<double>& freq, string& pattern)
{
	long 		i, j, n, page(1), length(info.size());
	long 		left(70), bottom(500), width(500), height(200);
	
	long 		npages = 2 + (length-1)/500;
	
	ofstream*	fps = ps_open_and_init(filename, title, npages, 600, 800);
		
	// The first page contains the information content and representation plots
	*fps << "%%Page: Info 1" << endl;
	*fps << "/Helvetica findfont 14 scalefont setfont" << endl;
	*fps << "40 740 moveto (" << filename << ": " << title << ") show" << endl;
	*fps << "/Left " << left << " def\n/Bottom " << bottom << " def" << endl;
	*fps << "/Height " << height << " def\n/Width " << width << " def" << endl;
	
	*fps << "/Data [\n%Pos Info nSeq" << endl;
	for ( i=0; i<length; i++ )
		*fps << i+1 << " " << info[i] << " " << nseq[i] << endl;
	*fps << "] def" << endl;

	*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
	*fps << "/x_Min 0 def\n/y_Min 0 def\n/x_Max " << length << " def\n/y_Max 4.3 def" << endl;
	*fps << "/x_Scale Width x_Max x_Min sub div def" << endl;
	*fps << "/y_Scale Height y_Max y_Min sub div def" << endl;
	*fps << "/Frame { newpath 0 0 moveto 0 Height lineto Width Height lineto Width 0 lineto closepath } def" << endl;

	*fps << left-45 << " " << bottom+40 << " moveto gsave 90 rotate (Information content) show grestore" << endl;

	ps_scale(fps, left, 200, left+width, 200, 0, length, 10*(int)(length/100 + 1),
			10, M_PI_2, 0, 12, 0);
	
	ps_scale(fps, left, bottom, left, bottom+height, 0, 4.3, 1.0,
			10, 0, 1, 12, 0);
	
	*fps << "gsave" << endl;
	*fps << "	Left Bottom translate Frame stroke Frame clip" << endl;
	*fps << "	/x Data 0 get x_Min sub x_Scale mul def" << endl;
	*fps << "	/y Data 1 get y_Min sub y_Scale mul def" << endl;
	*fps << "	x y moveto" << endl;
	*fps << "	3 3 Data length 3 sub {" << endl;
	*fps << "		/x_Index exch def" << endl;
	*fps << "		/y_Index x_Index 1 add def" << endl;
	*fps << "		/x Data x_Index get x_Min sub x_Scale mul def" << endl;
	*fps << "		/y Data y_Index get y_Min sub y_Scale mul def" << endl;
	*fps << "		x y lineto" << endl;
	*fps << "	} for stroke" << endl;
	*fps << "grestore" << endl;
	
	ps_seq_periodicity(fps, per);

	// The residue representation plot
	ps_seq_representation(fps, nseq);
	
	*fps << "showpage" << endl;
	
	// The following pages contain the sequence logo
	int 	pagestart = 0;
	for ( page=2; page<=npages; page++ ) {
		*fps << "%Page: Logo " << page << endl;
		*fps << "/Helvetica findfont 16 scalefont setfont" << endl;
		*fps << "40 740 moveto (Sequence Logo: " << filename << ") show" << endl;
		*fps << "/TopLine 640 def\n/LineHeight 60 def" << endl;
		*fps << "/FontSize 14 def" << endl;
		*fps << "/FontHeight FontSize 0.6 mul def" << endl;
		*fps << "/CharWidth FontSize 0.6 mul def" << endl;
		*fps << "/Courier-Bold findfont FontSize scalefont setfont" << endl;
		*fps << "/Letters (ABCDEFGHIJKLMNOPQRSTUVWXYZ) def" << endl;
		*fps << "/Red   [0.5 0.0 1.0 1.0 1.0 0.5 0.8 0.0 0.5 0.9 0.0 0.5 1.0 0.0 0.9 " << endl;
		*fps << "		0.8 0.0 0.0 0.8 0.8 0.9 0.5 0.5 0.1 0.5 0.0 ] def " << endl;
		*fps << "/Green [0.5 1.0 1.0 0.0 0.0 0.5 0.5 0.0 0.5 0.9 0.0 0.5 1.0 1.0 0.9 " << endl;
		*fps << "		0.5 1.0 0.0 0.5 0.5 0.9 0.5 0.5 0.1 0.5 1.0 ] def " << endl;
		*fps << "/Blue  [0.5 0.0 0.0 0.0 0.0 0.5 0.8 1.0 0.5 0.9 1.0 0.5 0.0 0.0 0.9 " << endl;
		*fps << "		0.8 0.0 1.0 0.8 0.8 0.9 0.5 0.5 0.1 0.5 0.0 ] def " << endl;
		*fps << "/SeqPattern [" << endl;
		for ( i=pagestart; i<pagestart+500 && i<length; i++ ) {
			*fps << "(";
			n = 0;
			for ( j=nres-1; j>=0; j-- ) {
				if ( pattern[nres*i+j] >= 'A' && freq[nres*i+j] > 0.1 ) {
//					*fps << "%c", pattern[nres*i+j]);
					*fps << pattern[nres*i+j];
					n++;
				}
			}
			if ( n < 1 ) *fps << "X";
			*fps << ")" << endl;
		}
		*fps << "] def" << endl;
		*fps << "/Height [" << endl;
		for ( i=pagestart; i<pagestart+500 && i<length; i++ ) {
			n = 0;
			for ( j=nres-1; j>=0; j-- ) {
				if ( pattern[nres*i+j] >= 'A' && freq[nres*i+j] > 0.1 ) {
					*fps << freq[nres*i+j] << " ";
					n++;
				}
			}
			if ( n < 1 ) *fps << "0.05 ";
			*fps << endl;
		}
		*fps << "] def" << endl;
		*fps << "/HeightIndex 0 def" << endl;
		*fps << "/LineIndex 0 def" << endl;
		*fps << "/theLine TopLine def" << endl;
		*fps << "gsave" << endl;
		*fps << "	40 theLine 15 sub moveto (" << pagestart << ") show" << endl;
		*fps << "	40 theLine moveto" << endl;
		*fps << "	0 1 SeqPattern length 1 sub { %% Get a pattern at each position" << endl;
		*fps << "		/PatternIndex exch def" << endl;
		*fps << "		LineIndex PatternIndex 50 div truncate lt {" << endl;
		*fps << "			/theLine theLine LineHeight sub def" << endl;
		*fps << "			/LineIndex PatternIndex 50 div truncate def" << endl;
		*fps << "			40 theLine 15 sub moveto" << endl;
		*fps << "			PatternIndex " << pagestart << " add (xxxx) cvs show" << endl;
		*fps << "			40 theLine moveto" << endl;
		*fps << "		} if" << endl;
		*fps << "		/LineIndex PatternIndex 50 div def" << endl;
		*fps << "		/aPattern {SeqPattern PatternIndex get} def" << endl;
		*fps << "		gsave" << endl;
		*fps << "			0 1 aPattern length 1 sub { %% Get the character and its colour and height" << endl;
		*fps << "				/CharIndex exch def" << endl;
		*fps << "				/aChar aPattern CharIndex 1 getinterval def" << endl;
		*fps << "				/CharHeight Height HeightIndex get def" << endl;
		*fps << "				0 1 Letters length 1 sub {" << endl;
		*fps << "					/anIndex exch def" << endl;
		*fps << "					aChar Letters anIndex 1 getinterval eq { /LetterIndex anIndex def } if" << endl;
		*fps << "				} for" << endl;
		*fps << "				gsave" << endl;
		*fps << "					Red LetterIndex get Green LetterIndex get Blue LetterIndex get setrgbcolor" << endl;
		*fps << "					1 CharHeight scale" << endl;
		*fps << "					aChar show" << endl;
		*fps << "				grestore" << endl;
		*fps << "				0 CharHeight FontHeight mul rmoveto" << endl;
		*fps << "				/HeightIndex HeightIndex 1 add def" << endl;
		*fps << "			} for" << endl;
		*fps << "		grestore" << endl;
		*fps << "		CharWidth 0 rlineto" << endl;
		*fps << "	} for stroke" << endl;
		*fps << "grestore" << endl;
		*fps << "showpage" << endl;
		pagestart += 500;
	}
	
	ps_close(fps);
	
	return 0;
}

/**
@brief 	Generates a postscript plot from a hydrophobicity analysis of a protein sequence alignment.
@param 	&filename	postscript file name.
@param 	&title 		title.
@param 	*Hphob 		moving average of average hydrophobicity at each position.
@param 	*HPseg		assignment of hydrophobic segments.
@param 	*nseq 		number of sequences represented at each position.
@param 	*per		periodicity analysis.
@return int 			0.
**/
int 		ps_seq_hydrophob(Bstring& filename, Bstring& title,
				vector<double>& Hphob, vector<int>& HPseg,
				vector<double>& nseq, vector<Complex<float>>& per)
{
	long 		i, length(Hphob.size());
	long 		left(70), bottom(500), width(500), height(200);
	double		y_min, y_max;
	
	ofstream*		fps = ps_open_and_init(filename, title, 1, 600, 800);
		
	// The first page contains the hydrophobicity and representation plots
	*fps << "%%Page: Hphob 1" << endl;
	*fps << "/Helvetica findfont 16 scalefont setfont" << endl;
	*fps << "40 740 moveto (" << filename << ": " << title << ") show" << endl;
	*fps << "/Left " << left << " def\n/Bottom " << bottom << " def" << endl;
	*fps << "/Height " << height << " def\n/Width " << width << " def" << endl;
	
	y_min = y_max = 0;
	*fps << "/Data [" << endl << "%Pos Hphob HPseg" << endl;
	for ( i=0; i<length; i++ ) {
		*fps << i+1 << " " << Hphob[i] << " " << HPseg[i] << endl;
		if ( y_min > Hphob[i] ) y_min = Hphob[i];
		if ( y_max < Hphob[i] ) y_max = Hphob[i];
	}
	*fps << "] def" << endl;

	y_min = 0.01*floor(100*y_min);
	y_max = 0.01*floor(100*y_max + 1);
	
	*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
	*fps << "/x_Min 0 def\n/y_Min " << y_min << " def" << endl;
	*fps << "/x_Max " << length << " def\n/y_Max " << y_max << " def" << endl;
	*fps << "/x_Scale Width x_Max x_Min sub div def" << endl;
	*fps << "/y_Scale Height y_Max y_Min sub div def" << endl;
	
	*fps << "/Frame { newpath 0 0 moveto 0 Height lineto Width Height lineto Width 0 lineto closepath } def" << endl;

	*fps << left-45 << " " << bottom+30 << " moveto gsave 90 rotate (Hydrophobicity (kcal/mol)) show grestore" << endl;

	ps_scale(fps, left, 200, left+width, 200, 0, length, 10*(int)(length/100),
			10, M_PI_2, 0, 12, 0);
	
	ps_scale(fps, left, bottom, left, bottom+height, y_min, y_max, (y_max - y_min)/5,
			10, 0, 1, 12, 0);
	
	*fps << "gsave" << endl;
	*fps << "	Left Bottom translate Frame stroke Frame clip stroke" << endl;
	*fps << "	20 setlinewidth 0.5 setgray" << endl;
	*fps << "	/x Data 0 get x_Min sub x_Scale mul def" << endl;
	*fps << "	/y 1 y_Min sub y_Scale mul def" << endl;
	*fps << "	x y moveto" << endl;
	*fps << "	3 3 Data length 3 sub {" << endl;
	*fps << "		/x_Index exch def" << endl;
	*fps << "		/x Data x_Index get x_Min sub x_Scale mul def" << endl;
	*fps << "		Data x_Index 2 add get 0 gt {" << endl;
	*fps << "			x y lineto" << endl;
	*fps << "		}{" << endl;
	*fps << "			x y moveto" << endl;
	*fps << "		} ifelse" << endl;
	*fps << "	} for stroke" << endl;
	*fps << "grestore" << endl;
	
	*fps << "gsave" << endl;
	*fps << "	Left Bottom translate Frame stroke Frame clip" << endl;
	*fps << "	/x Data 0 get x_Min sub x_Scale mul def" << endl;
	*fps << "	/y Data 1 get y_Min sub y_Scale mul def" << endl;
	*fps << "	x y moveto" << endl;
	*fps << "	3 3 Data length 3 sub {" << endl;
	*fps << "		/x_Index exch def" << endl;
	*fps << "		/y_Index x_Index 1 add def" << endl;
	*fps << "		/x Data x_Index get x_Min sub x_Scale mul def" << endl;
	*fps << "		/y Data y_Index get y_Min sub y_Scale mul def" << endl;
	*fps << "		x y lineto" << endl;
	*fps << "	} for stroke" << endl;
	*fps << "grestore" << endl;
	
	ps_seq_periodicity(fps, per);
	
	// The residue representation plot
	ps_seq_representation(fps, nseq);
	
	*fps << "showpage" << endl;
	
	ps_close(fps);
	
	return 0;
}

// The residue representation plot
int 		ps_seq_representation(ofstream* fps, vector<double>& nseq)
{
	int 		i, length(nseq.size()), bar_height(30);
	double		y_scale(bar_height);
	
	*fps << "/Data [" << endl << "%Pos nSeq" << endl;
	for ( i=0; i<length; i++ ) {
		*fps << i+1 << " " << nseq[i] << endl;
		if ( nseq[i]*y_scale > bar_height ) y_scale = bar_height/nseq[i];
	}
	*fps << "] def" << endl;
	*fps << "/Left 70 def\n/Bottom 70 def" << endl;
	*fps << "/Width 500 def" << endl;
	*fps << "/x_Min 1 def" << endl;
	*fps << "/x_Max " << length << " def" << endl;
	*fps << "/x_Scale Width x_Max x_Min sub div def" << endl;
	*fps << "/y_Scale " << y_scale << " def" << endl;
	*fps << "gsave" << endl;
	*fps << "	Left Bottom translate" << endl;
	*fps << "	/x Data 0 get x_Min sub x_Scale mul def" << endl;
	*fps << "	/y " << bar_height << " Data 1 get y_Scale mul add def" << endl;
	*fps << "	x y moveto" << endl;
	*fps << "	2 2 Data length 2 sub {" << endl;
	*fps << "		/x_Index exch def" << endl;
	*fps << "		/y_Index x_Index 1 add def" << endl;
	*fps << "		/x Data x_Index get x_Min sub x_Scale mul def" << endl;
	*fps << "		/y " << bar_height << " Data y_Index get y_Scale mul add def" << endl;
	*fps << "		x y lineto" << endl;
	*fps << "	} for stroke" << endl;
	*fps << "	/x Data 0 get x_Min sub x_Scale mul def" << endl;
	*fps << "	/y " << bar_height << " Data 1 get y_Scale mul sub def" << endl;
	*fps << "	x y moveto" << endl;
	*fps << "	2 2 Data length 2 sub {" << endl;
	*fps << "		/x_Index exch def" << endl;
	*fps << "		/y_Index x_Index 1 add def" << endl;
	*fps << "		/x Data x_Index get x_Min sub x_Scale mul def" << endl;
	*fps << "		/y " << bar_height << " Data y_Index get y_Scale mul sub def" << endl;
	*fps << "		x y lineto" << endl;
	*fps << "	} for stroke" << endl;
	*fps << "grestore" << endl;
	
	return 0;
}

// The periodicity plot
int 		ps_seq_periodicity(ofstream* fps, vector<Complex<float>>& per)
{
	int 		i, length(per.size());
	int 		left(70), bottom(300), width(500), height(200);
	double		y_max;
	
	*fps << "/Per [" << endl << "%Pos Amp Phi" << endl;
	y_max = 0;
	for ( i=0; i<length; i++ ) {
		*fps << i+1 << " " << per[i].amp() << " " << per[i].phi()*180/M_PI << endl;
		if ( y_max < per[i].amp() ) y_max = per[i].amp();
	}
	*fps << "] def" << endl;
	
	y_max = 0.1*floor(10*y_max + 1);
	
	*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
	*fps << "/Left " << left << " def\n/Bottom " << bottom << " def" << endl;
	*fps << "/Height " << height << " def\n/Width " << width << " def" << endl;
	*fps << "/x_Min 1 def\n/y_Min 0 def" << endl;
	*fps << "/x_Max " << length << " def\n/y_Max " << y_max << " def" << endl;
	*fps << "/x_Scale Width x_Max x_Min sub div def" << endl;
	*fps << "/y_Scale Height y_Max y_Min sub div def" << endl;
	*fps << left-45 << " " << bottom+70 << " moveto gsave 90 rotate (Amplitude) show grestore" << endl;
	
	ps_scale(fps, left, bottom, left, bottom+height, 0, y_max, y_max/4.5,
			10, 0, 1, 12, 0);

	*fps << "/Frame { newpath 0 0 moveto 0 Height lineto Width Height lineto Width 0 lineto closepath } def" << endl;
	*fps << "gsave" << endl;
	*fps << "	Left Bottom translate Frame stroke Frame clip" << endl;
	*fps << "	/x Per 0 get x_Min sub x_Scale mul def" << endl;
	*fps << "	/y Per 1 get y_Min sub y_Scale mul def" << endl;
	*fps << "	x y moveto" << endl;
	*fps << "	3 3 Per length 3 sub {" << endl;
	*fps << "		/x_Index exch def" << endl;
	*fps << "		/y_Index x_Index 1 add def" << endl;
	*fps << "		/x Per x_Index get x_Min sub x_Scale mul def" << endl;
	*fps << "		/y Per y_Index get y_Min sub y_Scale mul def" << endl;
	*fps << "		x y lineto" << endl;
	*fps << "	} for stroke" << endl;
	*fps << "grestore" << endl;
	
	bottom = 200;
	height = 100;

	*fps << "/Left " << left << " def\n/Bottom " << bottom << " def" << endl;
	*fps << "/Height " << height << " def\n/Width " << width << " def" << endl;
	*fps << left-40 << " " << bottom << " moveto (-180) show" << endl;
	*fps << left-40 << " " << bottom+height/2-8 << " moveto (   0) show" << endl;
	*fps << left-40 << " " << bottom+height-16 << " moveto ( 180) show" << endl;
	*fps << left-45 << " " << bottom << " moveto gsave 90 rotate (Phase (degrees)) show grestore" << endl;

	*fps << left+190 << " " << bottom-40 << " moveto (Alignment position) show" << endl;
	*fps << "/y_Min -180 def\n/y_Max 180 def" << endl;
	*fps << "/y_Scale Height y_Max y_Min sub div def" << endl;
	*fps << "/Frame { newpath 0 0 moveto 0 Height lineto Width Height lineto Width 0 lineto closepath } def" << endl;
	*fps << "gsave" << endl;
	*fps << "	Left Bottom translate Frame stroke Frame clip" << endl;
	*fps << "	/x Per 0 get x_Min sub x_Scale mul def" << endl;
	*fps << "	/y Per 2 get y_Min sub y_Scale mul def" << endl;
	*fps << "	x y moveto" << endl;
	*fps << "	/last_x x def" << endl;
	*fps << "	/last_y y def" << endl;
	*fps << "	3 3 Per length 3 sub {" << endl;
	*fps << "		/x_Index exch def" << endl;
	*fps << "		/y_Index x_Index 2 add def" << endl;
	*fps << "		/x Per x_Index get x_Min sub x_Scale mul def" << endl;
	*fps << "		/y Per y_Index get y_Min sub y_Scale mul def" << endl;
	*fps << "		/y_diff y last_y sub def" << endl;
	*fps << "		y_diff 0 lt {" << endl;
	*fps << "			y_diff Height 2 div neg gt { x y lineto }" << endl;
	*fps << "				{ x y Height add lineto last_x last_y Height sub moveto x y lineto" << endl;
	*fps << "			} ifelse" << endl;
	*fps << "		}{" << endl;
	*fps << "			y_diff Height 2 div lt { x y lineto }" << endl;
	*fps << "				{ x y Height sub lineto last_x last_y Height add moveto x y lineto" << endl;
	*fps << "			} ifelse" << endl;
	*fps << "		} ifelse" << endl;
	*fps << "		x y lineto" << endl;
	*fps << "		/last_x x def" << endl;
	*fps << "		/last_y y def" << endl;
	*fps << "	} for stroke" << endl;
	*fps << "grestore" << endl;
		
	return 0;
}

