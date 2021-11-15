/**
@file	Bgraphseg.h
@brief	Graph segmentation classes.
@author	Bernard Heymann
@date	Created: 20110318
@date	Modified: 20201228
**/

#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

#ifndef _GraphSeg_

/**
@brief Graph segmentation edge.
**/
class GSedge {
private:
	long		u, v;	// GSedge ends in the 1D data array
	float		w;		// GSedge weight: Neighbor value difference
public:
	GSedge() : u(0), v(0), w(0) {}
	GSedge(long start, long end, double weight) : u(start), v(end), w(weight) {}
	long		operator[](int i) { return (i%2)? u: v; }
	void		weight(double val) { w = val; }
	double		weight() const { return w; }
} ;

//bool sort_edge(const GSedge &e1, const GSedge &e2);

/**
@brief Graph segmentation region.
**/
class GSvoxel {
private:
	long		p;		// GSvoxel index
	int			r;		// Number of region joins
	int			s;		// Number of voxels in the region
	float		a;		// Average value for the region
public:
	GSvoxel() : p(0), r(0), s(0), a(0) {}
	GSvoxel(long index, long nj, long nv, double avg) : p(index), r(nj), s(nv), a(avg) {}
	void		index(long i) { p = i; }
	long		index() { return p; }
	void		joins(int i) { r = i; }
	int			joins() { return r; }
	void		increment_joins() { r++; }
	void		voxels(int i) { s = i; }
	int			voxels() { return s; }
	void		increment_voxels() { s++; }
	void		add_to_voxels(int i) { s += i; }
	void		average(double d) { a = d; }
	double		average() { return a; }
	void		add_to_average(double d) { a += d; }
	void		divide_average_by_voxels() { if ( s ) a /= s; }
	void		clear() { s=0; a=0; }
} ;

/**
@brief Graph segmentation container.
**/
class GSgraph {
private:
	vector<GSvoxel>		r;
	vector<GSedge>		e;
public:
	vector<GSvoxel>&	voxels() { return r; }
	vector<GSedge>&		edges() { return e; }
	GSvoxel&	voxel(long i) { return r[i]; }
	GSvoxel&	add_voxel(long index, long nj, long nv, double avg) {
		r.push_back(GSvoxel(index, nj, nv, avg));
		return r.back();
	}
	long		voxel_count() {
		long		i, j;
		for ( i=j=0; i<r.size(); ++i )
			if ( r[i].index() == i ) j++;
		return j;
	}
	GSedge&		edge(long i) { return e[i]; }
	GSedge&		add_edge(long u, long v, double w) {
		e.push_back(GSedge(u, v, w));
		return e.back();
	}
	long		edge_count() { return e.size(); }
	void		edge_sort() {
		sort (e.begin(), e.end(), [](const GSedge &e1, const GSedge &e2) { return e1.weight() < e2.weight(); });
	}
	long		region_find(long i) {
		long			j(i);
		while ( j != r[j].index() ) j = r[j].index();
		r[i].index(j);
		return j;
	}
	int			region_join(long i, long j) {
		if ( i == j ) return 0;
		int				s(r[i].voxels() + r[j].voxels());
		double			a((r[i].voxels()*r[i].average() + r[j].voxels()*r[j].average())/s);
		if ( r[i].joins() > r[j].joins() ) {
			r[j].index(i);
			r[i].add_to_voxels(r[j].voxels());
			r[i].average(a);
		} else {
			r[i].index(j);
			r[j].add_to_voxels(r[i].voxels());
			r[j].average(a);
			if ( r[i].joins() == r[j].joins() ) r[j].increment_joins();
		}
		return 1;
	}
	long		region_merging(double threshold) {
		long			i, u, v, nrc(r.size());
		GSedge			edge;
		for ( i=0; i<r.size(); ++i ) r[i].average(threshold);
		for ( i=0; i<e.size(); ++i ) {
			edge = e[i];
			u = region_find(edge[0]);
			v = region_find(edge[1]);
			if ( u != v ) {
				if ( edge.weight() <= r[u].average() && edge.weight() <= r[v].average() ) {
					nrc -= region_join(u, v);
					u = region_find(u);
					r[u].average(edge.weight() + threshold/r[u].voxels());
				}
			}
		}
		return nrc;
	}
	long		statistical_region_merging(double threshold) {
		// Find segment memberships
		long			i, u, v, nrc(r.size());
		double			fac(threshold*threshold);
		double			d2, ln1, ln2, s1, s2;
		double			lnd(2*log(6*nrc));
		for ( i=0; i<e.size(); ++i ) {
			u = region_find(e[i][0]);
			v = region_find(e[i][1]);
			if ( u != v ) {
				d2 = r[u].average() - r[v].average(); d2 *= d2;
				s1 = r[u].voxels(); s2 = r[v].voxels();
//				ln1 = lnd + log(1.0 + s1) * ((g < s1)? g: s1);
//				ln2 = lnd + log(1.0 + s2) * ((g < s2)? g: s2);
				ln1 = lnd + log(1.0 + s1);
				ln2 = lnd + log(1.0 + s2);
				if ( d2 < fac * (ln1/s1 + ln2/s2) )
					nrc -= region_join(u, v);
			}
		}
		return nrc;
	}
	long		region_merge_small(long nrc, long min_size)	{
		if ( min_size < 1 ) return nrc;
		long			i, u, v;
		for ( i=0; i<e.size(); ++i ) {
			u = region_find(e[i][0]);
			v = region_find(e[i][1]);
			if ( ( u != v ) && ( r[u].voxels() < min_size || r[v].voxels() < min_size ) ) {
				nrc -= region_join(u, v);
			}
		}
		return nrc;
	}
	long		region_number() {
		long		i, j;
		for ( i=j=0; i<r.size(); ++i )
			if ( r[i].index() == i ) r[i].joins(++j);
		for ( i=0; i<r.size(); ++i ) {
			j = region_find(i);
			r[i].joins(r[j].joins());
			r[i].voxels(r[j].voxels());
		}
		return j;
	}
	void		show_regions() {
		long		i;
		cout << "Region\tSize\tAverage" << endl;
		for ( i=0; i<r.size(); ++i )
			if ( r[i].index() == i )
				cout << r[i].joins() << tab << r[i].voxels() << tab << r[i].average() << endl;
	}
} ;

#define _GraphSeg_
#endif
