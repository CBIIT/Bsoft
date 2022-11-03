/**
@file	Bplot.h
@brief	Header file for postscript output functions.
@author Bernard Heymann 
@date	Created: 20010515
@date	Modified: 20221101
**/
 
#include "Bstring.h"
#include "utilities.h"
#include <fstream>

#ifndef _Bplot_
#define _Bplot_
class Baxis {
private:
	long		id;			// Axis: 1=x, 2=x2, 3=y, 4=y2, 5=v
	Bstring		lab;		// Label
	double		mn;			// Minimum
	double		mx;			// Maximum
	double		ic;			// Increment
	double		cl[3];		// Color
	long		flg;		// Flags: 1=inverse, 2=logarithm
public:
	Baxis() { id=0; lab=0; mn=0; mx=0; ic=0; flg=0; }
	long		which() { return id; }
	void		which(long i) { id = i; }
	Bstring&	label() { return lab; }
	void		label(const Bstring& s) { lab = s; }
	double		min() { return mn; }
	void		min(double v) { mn = v; }
	double		max() { return mx; }
	void		max(double v) { mx = v; }
	double		inc() { return ic; }
	void		inc(double v) { ic = v; }
	double*		color() { return cl; }
	void		color(double r, double g, double b) { cl[0]=r; cl[1]=g; cl[2]=b; }
	long		flags() { return flg; }
	void		flags(long v) { flg = v; }
} ;

class Bcolumn {
private:
	long		num;		// Column index
	Bstring		lab;		// Label
	long		tp;			// Type of plot: 0=points/dots, 1=bars, 2=line, 3=error, 4=variable dots
	long		ax;			// Axis: 0=noaxis, 1=x, 2=x2, 3=y, 4=y2, 5=v
	double		mn;			// Minimum
	double		mx;			// Maximum
	double		es;			// Size of dot or line
	double		cl[3];		// Color
	void		initialize() {
		num=0; lab=0;
		tp=ax=0;
		mn=mx=es=0;
		cl[0]=cl[1]=cl[2]=0;
	}
public:
	Bcolumn() { initialize(); }
	Bcolumn(Bstring& s) { initialize(); lab=s; }
	Bcolumn(long i, Bstring& s) { initialize(); num=i; lab=s; }
	long		number() { return num; }
	void		number(long i) { num = i; }
	Bstring&	label() { return lab; }
	void		label(const Bstring& s) { lab = s; }
	void		label(const char* s) { lab = s; }
	long		type() { return tp; }
	void		type(long v) { tp = v; }
	long		axis() { return ax; }
	void		axis(long v) { ax = v; }
	double		min() { return mn; }
	void		min(double v) { mn = v; }
	double		max() { return mx; }
	void		max(double v) { mx = v; }
	double		element_size() { return es; }
	void		element_size(double v) { es = v; }
	double*		color() { return cl; }
	void		color(double r, double g, double b) { cl[0]=r; cl[1]=g; cl[2]=b; }
	double		red() { return cl[0]; }
	double		green() { return cl[1]; }
	double		blue() { return cl[2]; }
} ;

class Bpage {
private:
	Bstring		tit;		// Page title
	long		num;		// Page number
	long		nc;			// Number of columns
//	Bcolumn*	col;		// Columns
	vector<Bcolumn>	col;	// Columns
	Baxis		ax[5];		// Axes: x, x2, y, y2, v
	Bstring*	txt;		// Set of strings to be plotted at the bottom
	void		initialize() {
		tit = 0; txt = NULL; num = 0; nc = 0;
//		col = NULL;
		axes();
	}
public:
	Bpage() { initialize(); }
	Bpage(long n, long ncol) { initialize(); num = n; columns(ncol); }
//	~Bpage() { tit = 0; delete[] col; col = NULL; }
	~Bpage() { tit = 0; col.clear(); }
	Bstring&	title() { return tit; }
	void		title(const Bstring& t) { tit = t; }
	void		title(const char* t) { tit = t; }
	void		number(long n) { num = n; }
	long		number() { return num; }
	void		columns(long ncol) {
		nc = ncol;
//		for ( long i=0; i<nc; ++i ) col.push_back(Bcolumn());
//		col = new Bcolumn[nc];
		col.resize(nc);
	}
	long		columns() { return nc; }
	void		add_columns(long ncol) {
		nc += ncol;
		col.resize(nc);
	}
	Bcolumn&	column(long i) {
		if ( i<nc ) return col[i];
		cerr << "Error: The column requested, " << i << " for page " << num <<
			" is greater than the number of columns, " << nc << endl << endl;
		exit(-1);
	}
	void		axes() { for ( long i=0; i<5; i++ ) ax[i].which(i+1); }
	Baxis&		axis(long id) {
		if ( id > 0 && id < 6 ) return ax[id-1];
		else return ax[0];
	}
	Bstring*	text() { return txt; }
	void		add_text(Bstring& s) { string_add(&txt, s); }
	void		limits() {
//		if ( !col ) return;
		if ( col.size() < 1 ) return;
		long		c, i;
		for ( c=0; c<nc; c++ ) {
			i = col[c].axis();
			if ( i > 0 && i < 6 ) {
				i--;
//				cout << "col[" << c << "]: min=" << col[c].min() << " max=" << col[c].max() << endl;
				if ( ax[i].min() == ax[i].max() ) {
					ax[i].min(col[c].min());
					ax[i].max(col[c].max());
				}
//				cout << "col[" << c << "] axis[" << i+1 << "]: min=" << ax[i].min() << " max=" << ax[i].max() << endl;
			}
		}
		for ( i=0; i<5; i++ ) {
			if ( ax[i].min() == ax[i].max() ) ax[i].max(ax[i].min() + 1);
			if ( ax[i].inc() == 0 ) ax[i].inc((ax[i].max() - ax[i].min())/10);
//			cout << "axis[ " << i+1 << "]: min=" << ax[i].min() << " max=" << ax[i].max() << endl;
		}
	}
} ;

class Bplot {
private:
	Bstring		tit;		// Plot title
	long		pw, ph;		// Page size
	long		gw,	gh;		// Plot size
	long		gl, gb;		// Plot left and bottom
	long		np;			// Number of pages
	long		nr;			// Number of points (rows)
	long		nc;			// Number of columns (x is first)
	Bpage*		pg;			// Array of page information objects
	double*		data;		// Plot data
	void	initialize(long npage=0, long nrow=0, long ncol=0) {
		tit = 0;
		pw = 600; ph = 800;
		gw = 450; gh = 450;
		gl = 70; gb = 200;
		np = npage;
		nr = nrow;
		nc = ncol;
		pg = NULL;
		data = NULL;
		if ( np ) pages(np);
		if ( nc && nr ) {
			data = new double[nc*nr];
			for ( long i=0; i<nr*nc; i++ ) data[i] = 0;
		}
	}
public:
	Bplot() { initialize(); }
	Bplot(long npage, long nrow, long ncol) {
		initialize(npage, nrow, ncol);
	}
	Bplot(Bstring& filename);
	Bplot(Bstring& filename, long skip, int type);
	~Bplot() {
		tit = 0;
		delete[] pg; pg = NULL;
		delete[] data; data = NULL;
	}
	void		reset(long npage, long nrow, long ncol) {
		if ( pg ) delete[] pg;
		if ( data ) delete[] data;
		initialize(npage, nrow, ncol);
	}
	double&	operator[](long i) {
		while ( i < 0 ) i += nr*nc;
		while ( i >= nr*nc ) i -= nr*nc;
		return data[i];
	}
	Bstring&	title() { return tit; }
	void		title(const Bstring& t) { tit = t; }
	void		title(const char* t) { tit = t; }
	void		page_size(long w, long h) { pw = w; ph = h; }
	long		page_width() { return pw; }
	long		page_height() { return ph; }
	void		size(long w, long h) { gw = w; gh = h; }
	long		width() { return gw; }
	long		height() { return gh; }
	void		origin(long l, long b) { gl = l; gb = b; }
	long		left() { return gl; }
	long		bottom() { return gb; }
	void		pages(long n) {
		np = n;
		pg = new Bpage[np];
		for ( long i=0; i<np; i++ ) pg[i].number(i+1);
	}
	long		pages() { return np; }
	void		rows(long n) { nr = n; }
	long		rows() { return nr; }
	void		columns(long n) { nc = n; }
	long		columns() { return nc; }
	Bpage&		page(long i) { // cout << "Getting page " << i << endl;
		if ( pg && i >=0 && i < np ) return pg[i];
		else {
			cerr << "Error: The plot page " << i << " cannot be found!" << endl;
			return pg[0];
		}
	}
	double*		data_pointer() { return data; }
	void		add_columns(int num) {
		int			old_nc = nc;
		double*		old_data = data;
		nc += num;
		data = new double[nc*nr];
		for ( long i=0; i<nr*old_nc; ++i ) data[i] = old_data[i];
	}
	void		limits() {
		if ( !pg ) return;
		long		i, j, n, c;
		double		min, max;
		Bcolumn		col;
		for ( n=0; n<np; n++ ) {
			for ( c=0; c<pg[n].columns(); c++ ) {
				col = pg[n].column(c);
				if ( col.min() == col.max() ) {
					min = 1e100; max = -1e100;
					for ( i=col.number()*nr, j=0; j<nr; i++, j++ ) {
						if ( min > data[i] ) min = data[i];
						if ( max < data[i] ) max = data[i];
					}
					pg[n].column(c).min(min);
					pg[n].column(c).max(max);
//					cout << "col[" << c << "]: min=" << min << " max=" << max << endl;
				}
			}
		}
	}
	void		add(Bplot* p, long start=0) {
		for ( long i=start*nr; i<nr*nc; i++ ) data[i] += p->data[i];
	}
	void		multiply(double d) {
		for ( long i=0; i<nr*nc; i++ ) data[i] *= d;
	}
	double		cut(long icol, double threshold, int dir=1) {
		long			i, i1, j, j1;
		double			val((data[0]+data[1])/2);
		for ( i=0, i1=1, j=icol*nr, j1=j+1; i1<nr; i++, i1++, j++, j1++ ) {
			if ( dir < 0 ) {	// decreasing threshold
				if ( data[j] >= threshold && data[j1] < threshold )
					val = data[i] + (data[i1] - data[i])*(threshold - data[j])/(data[j1] - data[j]);
			} else {			// increasing threshold
				if ( data[j] <= threshold && data[j1] > threshold )
					val = data[i] + (data[i1] - data[i])*(threshold - data[j])/(data[j1] - data[j]);
			}
		}
		if ( val < data[1] ) {
			if ( dir < 0 ) {
				if ( data[j] > threshold ) val = data[i];
			} else {
				if ( data[j] < threshold ) val = data[i];
			}
		}
		return val;
	}
	void		write_tsv(Bstring& filename) {
		long			i, j, k;
		ofstream		ftsv(filename.c_str());
		if ( ftsv.fail() ) return;
//		ftsv << page(0).column(0).label();
//		for ( j=1; j<nc; j++ ) ftsv << tab << page(0).column(j).label();
//		ftsv << endl;
		for ( i=0; i<nr; i++ ) {
			ftsv << data[i];
			for ( j=1, k=nr+i; j<nc; j++, k+=nr ) ftsv << tab << data[k];
			ftsv << endl;
		}
		ftsv.close();
	}
	void		resolution_display(vector<double>& fsccut, vector<double>& dprcut);
} ;

ostream& operator<<(ostream& output, Bplot* plot);

#endif

