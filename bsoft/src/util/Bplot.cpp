/**
@file	Bplot.cpp
@brief	Postscript output functions.
@author Bernard Heymann
@date	Created: 20010515
@date	Modified: 20160318
**/
 
#include "Bplot.h"
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

ostream& operator<<(ostream& output, Bplot* plot) {
	long			i, j, k;
//	output.setf(ios::fixed, ios::floatfield);
	output << plot->page(0).column(0).label();
	for ( j=1; j<plot->columns(); j++ ) cout << tab << plot->page(0).column(j).label();
	output << endl;
	for ( i=0; i<plot->rows(); i++ ) {
		output << (*plot)[i];
		for ( j=1, k=plot->rows()+i; j<plot->columns(); j++, k+=plot->rows() )
			output << tab << (*plot)[k];
		output << endl;
	}
	return output;
}

Bplot::Bplot(Bstring& filename)
{
	if ( filename.length() < 1 ) return;
	
	if ( verbose & VERB_PROCESS )
		cout << "# Reading plot file:       " << filename << endl;

	ifstream			fps(filename.c_str());
	if ( fps.fail() ) return;

	char				aline[1024];
	long				i(0), j;
	int					ncol(0), nrow(0);
	Bstring				sline;

	// Title
	while ( fps.getline(aline, 1024) && strncmp(aline, "%%Title:", 8) ) ;
	
	Bstring				plot_title(aline);
	plot_title = plot_title.substr(8,plot_title.length());

	// Start of data
	while ( fps.getline(aline, 1024) && strncmp(aline, "/Data [", 7) ) ;

	// Lables line
	fps.getline(aline, 1024);
	sline = aline;
	Bstring*			labels = sline.split();
	Bstring*			label;
	
	// Determine number of rows and columns
	while ( fps.getline(aline, 1024) && aline[0] != ']' ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bplot::Bplot: aline = \"" << aline << "\"" << endl;
		sline = aline;
		sline = sline.substr(0, sline.length() - 1);
		vector<double>	temp = sline.split_into_doubles(0);
		if ( ncol < temp.size() ) ncol = temp.size();
		nrow++;
	}
	
	if ( verbose ) {
		cout << "Plot file:                      " << filename << endl;
		cout << "Columns:                        " << ncol << endl;
		cout << "Rows:                           " << nrow << endl;
		cout << "Labels:                         ";
	}

	initialize(1, nrow, ncol);
	title(plot_title);
	page(0).title(plot_title);
	page(0).columns(ncol);
	for ( i=0, label = labels; i<ncol; i++, label = label->next ) {
		page(0).column(i).number(i);
		page(0).column(i).label(*label);
		if ( verbose )
			cout << " " << *label;
	}

	if ( verbose )
		cout << endl << endl;

	string_kill(labels);
	
	fps.seekg(0);
	
	// Start of data
	while ( fps.getline(aline, 1024) && strncmp(aline, "/Data [", 7) ) ;

	// Lables line
	fps.getline(aline, 1024);

	// Read the data
	i = 0;
	while ( fps.getline(aline, 1024) && aline[0] != ']' ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bplot::Bplot: aline = \"" << aline << "\"" << endl;
		sline = aline;
		vector<double>	temp = sline.split_into_doubles(0);
		for ( j=0; j<ncol; j++ ) data[j*nrow+i] = temp[j];
		i++;
	}
}

void		Bplot::resolution_display(double* fsccut, double* dprcut)
{
	long		i;
	double		res(999), FSCsum(0), w, IVR;
	double		*s, *FSC, *DPR;
	Bstring		txt, label(page(0).axis(3).label());
	
	s = data;
	FSC = s + nr;
	DPR = FSC + nr;
	
	cout << "Resolution curves:" << endl;
	cout << "s(1/A)\tRes(A)\t" << label;
	if ( nc > 2 ) cout << "\tDPR";
	cout << endl;
	
	for ( i=0; i<nr; i++ ) {
		if ( i ) res = 1/s[i];
		w = M_PI*2*i;
		if ( label == "FSC" ) w *= 2.0*i;
		FSCsum += FSC[i]*w;
		cout << showpoint << fixed << setprecision(5) << setw(7) << s[i] << tab
					<< setprecision(3) << setw(7) << res << tab
					<< setw(7) << FSC[i];
		if ( nc > 2 ) cout << tab << setw(7) << DPR[i];
		cout << endl;
	}
	
	if ( label == "FSC" ) IVR = pow(4.0*M_PI/(3.0*FSCsum), 1.0/3.0);
	else IVR = sqrt(M_PI/FSCsum);
	IVR /= s[1];

	if ( fsccut || dprcut )
		cout << endl << "Resolution estimates:" << endl;
	
	if ( fsccut ) for ( i=0; i<4; i++ ) if ( fsccut[i] ) {
		txt = label + Bstring(fsccut[i], "(%g): ") + Bstring(1/cut(1,fsccut[i],-1), "%g A");
		page(0).add_text(txt);
		cout << txt << endl;
	}

	if ( dprcut ) for ( i=0; i<4; i++ ) if ( dprcut[i] ) {
		txt = Bstring(dprcut[i], "DPR(%g): ") + Bstring(1/cut(2,dprcut[i],1), "%g A");
		page(0).add_text(txt);
		cout << txt << endl;
	}

	txt = Bstring(IVR, "Information volume radius: %g A");
	page(0).add_text(txt);
	cout << txt << endl << endl;
}


