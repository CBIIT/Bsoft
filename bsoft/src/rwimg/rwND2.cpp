/**
@file	rwND2.cpp
@brief	Functions for reading Nikon ND2 files
@author Bernard Heymann
@date	Created: 20210627
@date 	Modified: 20210628
**/

#include "rwND2.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Fixed values for the ND2 format
int			signature(180276954);
string		mastertag("ND2 CHUNK MAP SIGNATURE 0000001!");

class ChunkHeader {
public:
	int			sig;
	int			off;
	long		len;
	string		tag;
	ChunkHeader() {}
	ChunkHeader(int s, int o, long l) : sig(s), off(o), len(l) { }
};

class ChunkEntry {
public:
	string		tag;
	long		pos;
	long		len;
	ChunkEntry() {}
	ChunkEntry(string t, long p, long l) : tag(t), pos(p), len(l) {}
};

class ImageAttribute {
public:
	bool		type;	// 0=number, 1=string
	string		tag;
	double		val;
	string		sval;
	ImageAttribute() : type(0), val(0) {}
};

string			read_tag(ifstream& fimg)
{
	char			c;
	string 			tag;
	
	while ( fimg.get(c) ) {
		tag.push_back(c);
		if ( c == '!' ) break;
	}
	
	return tag;
}

string			read_param(ifstream& fimg)
{
	char			c;
	long			i(0);
	string			param;
	
	while ( fimg.get(c) ) {
		if ( i%2 == 0 ) {
//			if ( isalnum(c) )
			if ( c )
				param.push_back(c);
			else break;
		}
		++i;
	}
	
	fimg.get(c);
	
	return param;
}

ChunkHeader		read_chunk_header(ifstream& fimg)
{
	ChunkHeader		chunk;
	
	fimg.read((char *)&chunk.sig, sizeof(int));
	if ( chunk.sig != signature ) {
		cerr << "Signature does not agree!" << endl;
		exit(-1);
	}
	
	fimg.read((char *)&chunk.off, sizeof(int));
	fimg.read((char *)&chunk.len, sizeof(long));
	
	chunk.tag = read_tag(fimg);
	
	if ( verbose & VERB_DEBUG_ND2 ) {
		cout << chunk.sig << "\t" << chunk.off << "\t" << chunk.len << endl;
		cout << chunk.tag << "\t" << chunk.tag.length() << endl;
	}
	
	return chunk;
}

ChunkEntry		read_chunk_entry(ifstream& fimg)
{
	ChunkEntry		entry;

	entry.tag = read_tag(fimg);

	fimg.read((char *)&entry.pos, sizeof(long));
	fimg.read((char *)&entry.len, sizeof(long));

	return entry;
}

ImageAttribute	read_attribute(ifstream& fimg)
{
	char			c(0);
	short			si(0);
	int				ii(0);
	unsigned int	ui(0);
	long			ll(0);
	unsigned long	ul(0);
	double			d;
	ImageAttribute	attr;
	
//	while ( si == 0 )
		fimg.read((char *)&si, sizeof(short));
	
	attr.tag = read_param(fimg);
	
	switch ( si%256 ) {
		case 1:
			fimg.get(c);
			attr.val = c;
			break;
		case 2:
			fimg.read((char *)&ii, sizeof(int));
			attr.val = ii;
			break;
		case 3:
			fimg.read((char *)&ui, sizeof(int));
			attr.val = ui;
			break;
		case 4:
			fimg.read((char *)&ll, sizeof(long));
			attr.val = ll;
			break;
		case 5:
			fimg.read((char *)&ul, sizeof(long));
			attr.val = ul;
			break;
		case 6:
			fimg.read((char *)&d, sizeof(double));
			attr.val = d;
			break;
		case 7:		// Pointer
			fimg.read((char *)&ll, sizeof(long));
			attr.val = ll;
			break;
		case 8:
			attr.type = 1;
			attr.sval = read_param(fimg);
			break;
		case 11:
			for ( int i=0; i<3; ++i ) fimg.read((char *)&ii, sizeof(int));
			attr.val = ii;
			break;
	}

	if ( verbose & VERB_DEBUG_ND2 ) {
		cout << attr.tag << "\t" << attr.val << endl;
		cout << fimg.tellg() << endl;
	}
	
	return attr;
}

int 	readND2(Bimage* p, int readdata, int img_select)
{
	string		infile;
	int			vi;
	long		chunkmap, len;
	vector<ChunkEntry>	entry_list;
	vector<ChunkHeader>	chunk_list;
	
    ifstream	fimg(p->file_name());
 	if ( fimg.fail() ) return -1;
   
	fimg.read((char *)&vi, sizeof(int));
	if ( verbose & VERB_DEBUG_ND2 )
		cout << "Signature: " << vi << endl;
    
    fimg.seekg(-8, fimg.end);
    
 	fimg.read((char *)&chunkmap, sizeof(long));
	if ( verbose & VERB_DEBUG_ND2 )
		cout << "Chunk map: " << chunkmap << endl;
    
	fimg.seekg(chunkmap, fimg.beg);

	chunk_list.push_back(read_chunk_header(fimg));

    while ( !fimg.eof() ) {
		ChunkEntry	entry = read_chunk_entry(fimg);

		if ( entry.tag == mastertag ) break;

		entry_list.push_back(entry);
	}
	
	if ( verbose & VERB_DEBUG_ND2 )
		for ( auto it: entry_list )
			cout << it.tag << "\t" << it.pos << "\t" << it.len << endl;

	long			i, j, start, line_size(0);
	double			px(1);
	string			param;
	vector<ImageAttribute>	attr;
	map<int,long>	img_pos;
	
	for ( auto it: entry_list ) {
		if ( it.tag == "ImageAttributesLV!" ) {
			fimg.seekg(it.pos, fimg.beg);
			chunk_list.push_back(read_chunk_header(fimg));
			ChunkHeader&  header = chunk_list.back();
			start = it.pos + 18 + header.off;
			fimg.seekg(start, fimg.beg);
			param = read_param(fimg);
			fimg.read((char *)&vi, sizeof(int));
			if ( verbose & VERB_DEBUG_ND2 ) {
				cout << start << endl;
				cout << param << endl;
				cout << fimg.tellg() << endl;
				cout << "Attributes = " << vi << endl;
			}
			fimg.read((char *)&len, sizeof(long));
			for ( i=0; i<vi; ++i )
				attr.push_back(read_attribute(fimg));
		}
		if ( it.tag.find("ImageDataSeq") != string::npos ) {
			vi = stoi(it.tag.substr(it.tag.find("|")+1));
			fimg.seekg(it.pos, fimg.beg);
			chunk_list.push_back(read_chunk_header(fimg));
			ChunkHeader&  header = chunk_list.back();
			img_pos[vi] = it.pos + 18 + header.off;
		}
		if ( it.tag.find("ImageCalibration") != string::npos ) {
			fimg.seekg(it.pos, fimg.beg);
			chunk_list.push_back(read_chunk_header(fimg));
			ChunkHeader&  header = chunk_list.back();
			start = it.pos + 18 + header.off;
			fimg.seekg(start, fimg.beg);
			param = read_param(fimg);
			fimg.read((char *)&vi, sizeof(int));
			if ( verbose & VERB_DEBUG_ND2 ) {
				cout << start << endl;
				cout << param << endl;
				cout << fimg.tellg() << endl;
				cout << "Calibration = " << vi << endl;
			}
			fimg.read((char *)&len, sizeof(long));
			for ( i=0; i<vi; ++i )
				attr.push_back(read_attribute(fimg));
		}
	}
	
	if ( verbose & VERB_DEBUG_ND2 )
		for ( auto it: img_pos )
			cout << it.first << "\t" << it.second << endl;

	for ( auto it: attr ) {
		if ( it.tag == "uiWidth" ) p->sizeX(it.val);
		if ( it.tag == "uiWidthBytes" ) line_size = it.val;
		if ( it.tag == "uiHeight" ) p->sizeY(it.val);
		if ( it.tag == "uiComp" ) p->channels(it.val);
		if ( it.tag == "uiBpcInMemory" ) {
			switch (int(it.val)) {
				case 8: p->data_type(UCharacter); break;
				case 16: p->data_type(UShort); break;
				case 32: p->data_type(UInteger); break;
				default: p->data_type(UCharacter);
			}
		}
		if ( it.tag == "uiSequenceCount" ) p->images(it.val);
		if ( it.tag == "dCalibration" ) px = 1.0e4*it.val;
	}
	
	p->sampling(px, px, 1);
	p->sizeZ(1);
	if ( p->channels() == 1 ) p->compound_type(TSimple);
	else p->compound_type(TMulti);
	p->data_offset(img_pos[0]);

	if ( readdata ) {
		p->data_alloc();
		unsigned char*	data = p->data_pointer();
		// Should this be line_size?
		long			readsize(p->sizeX()*p->channels()*p->data_type_size());
		for ( auto it: img_pos ) {
			for ( i=0, j=it.second; i<p->sizeY(); ++i, j+=line_size ) {
				fimg.seekg(j, fimg.beg);
				fimg.read((char *)data, readsize);
				data += readsize;
			}
		}
	}
	
	fimg.close();
    
	return 0;
}

