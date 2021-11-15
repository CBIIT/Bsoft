/**
@file	Bsuperpixel.h
@brief	Super pixel segmentation.
@author	Bernard Heymann
@date	Created: 20160320
@date	Modified: 20200502
**/

#include "utilities.h"

#ifndef _Bsuperpixel_

#define	NNEIGHBOR	64

class Bsuperpixel {
private:
	long			idx;	// index
	long			cnt;	// pixel count
	vector<double>	coor;	// coordinate average
	vector<double>	chan;	// channel averages
	vector<double>	var;	// channel variances
	double			wt;		// weight
	vector<long>	nb;		// neigbor list
public:
	Bsuperpixel() : idx(0), cnt(0), wt(0) {
		nb.resize(NNEIGHBOR, -1);
		for ( long i=0; i<3; ++i ) coor.push_back(0);
	}
	void			index(long i) { idx = i; }
	long			index() { return idx; }
	void			count(long i) { cnt = i; }
	long			count() { return cnt; }
	void			coordinates(vector<double>& c) { coor = c; }
	vector<double>&	coordinates() { return coor; }
	void			channels(long nc) {
		chan.resize(nc, 0);
		var.resize(nc, 0);
	}
	void			channels(vector<double>& v) {
		if ( v.size() != chan.size() ) {
			cerr << "Error: vector != channels" << endl;
			return;
		}
		for ( size_t i=0; i<v.size(); ++i ) chan[i] = v[i];
	}
	void			variances(vector<double>& v) {
		if ( v.size() != var.size() ) {
			cerr << "Error: vector != channels" << endl;
			return;
		}
		for ( size_t i=0; i<v.size(); ++i ) var[i] = v[i];
	}
	vector<double>&	channels() { return chan; }
	vector<double>&	variances() { return var; }
	double			channel(long i) { return chan[i]; }
	double			variance(long i) { return var[i]; }
	void			weight(double w) { wt = w; }
	double			weight() { return wt; }
	void			neighbor(long i, long j) { nb[i] = j; }
	long			neighbor(long i) { return nb[i]; }
	void			neighbors(vector<long>& v) {
		for ( size_t i=0; i<NNEIGHBOR && i<v.size(); ++i ) nb[i] = v[i];
	}
	void			add_neighbor(long i) {
		long		j;
		for ( j=0; j<NNEIGHBOR && nb[j] >= 0; ++j ) ;
		if ( j < NNEIGHBOR ) nb[j] = i;
		else {
			cerr << "Too many neighbors for segment " << idx << endl;
			exit(-1);
		}
	}
	void			clear_neighbors() {
//		for ( long i=0; i<NNEIGHBOR; i++ ) nb[i] = -1;
		fill(nb.begin(), nb.end(), -1);
	}
	void			clear() {
		cnt = 0;
		fill(coor.begin(), coor.end(), 0);
		fill(chan.begin(), chan.end(), 0);
		fill(var.begin(), var.end(), 0);
		wt = 0;
	}
	void			add(long cx, long cy, long cz, vector<double> v) {
		cnt++;
		coor[0] += cx;
		coor[1] += cy;
		coor[2] += cz;
		for ( size_t i=0; i<chan.size() && i<v.size(); ++i ) {
			chan[i] += v[i];
			var[i] += v[i]*v[i];
		}
	}
	void			add(vector<long> c, vector<double> v) {
		add(c[0], c[1], c[2], v);
	}
	void			average() {
		if ( cnt ) {
			coor /= cnt;
			for ( size_t i=0; i<chan.size(); ++i ) {
				chan[i] /= cnt;
				var[i] /= cnt;
				var[i] -= chan[i] * chan[i];
			}
		}
	}
	double			difference2(vector<double> v) {
		double 			d(0);
		for ( size_t i=0; i<chan.size() && i<v.size(); ++i )
			d += (chan[i] - v[i])*(chan[i] - v[i]);
		return d;
	}
	void			scale(double s) {
		coor *= s;
		cnt *= s;
	}
	void			show() {
		cout << idx << tab << cnt;
		for ( auto const& i : coor ) cout << tab << i;
		for ( auto const& i : chan ) cout << tab << i;
		for ( auto const& i : var ) cout << tab << i;
		cout << endl;
	}
};

#define _Bsuperpixel_
#endif
