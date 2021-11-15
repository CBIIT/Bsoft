/**
@file	rwEER.cpp
@brief	Reading EER files
@author Bernard Heymann
@date	Created: 20210427
@date 	Modified: 20210731
	
**/

#include "rwEER.h"

#ifdef HAVE_XML
#include "rwxml.h"
#endif

#include "ElectronCountedFramesDecompressor.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief	Reading an electron event representation file format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@param 	img_select	image selection in multi-image file (-1 = all images).
@param 	supres		super-resolution level.
@return	int			error code (<0 means failure).
A format for recording individual electron events on a direct detector.
	Uses the TIFF library with custom compression (a type of run-length encoding).
**/
int			readEER(Bimage* p, int readdata, int img_select, int supres)
{
	if ( supres == 0 || supres == -1 ) supres = 1;
	
	unsigned			x, y, z, f;
	ElectronCountedFramesDecompressor ecfd(p->file_name());
	
	JSvalue				mdata = json_from_xml(ecfd.GetAcquisitionMetadata().c_str());
	if ( verbose & VERB_DEBUG_EER )
		cout << mdata << endl;
	double				dose = mdata["totalDose"].real();
	double				exptime = mdata["exposureTime"].real();

	(*p)["dose"] = dose;
	(*p)["exposure"] = exptime;

	ecfd.getSize(x,y,z);
	f = ecfd.getNFrames();
	
    if (supres > 0) { x *= supres; y *= supres; }
    if (supres < 0) { x /= (-supres); y /= (-supres); }
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readEER: size: " << x << " " << y << " " << z << " " << f << endl;
//		cout << ecfd.GetAcquisitionMetadata() << endl;
	}

	p->size(x,y,1);
	p->channels(1);
	p->data_type(Float);
	
	long				nimg(f);		// Read all frames
	if ( img_select >= 0 ) nimg = 1;	// Read one frame
	else img_select = 0;
	p->images(nimg);
	
//	p->information();
	p->meta_data_update();

	if ( !readdata ) return 0;
	
	p->data_alloc();

	long				i, j, k(0), n, ne, nc(0), ncp(0);
	
	if ( img_select >= 0 ) i = img_select;
	for ( i=img_select; i<img_select+nimg; ++i ) {
		vector<unsigned char>	d(x*y, 0);
		ecfd.decompressImage_AddTo(&d[0], supres, i);
		if ( nimg == 1 ) k = 0;		// Select one frames
		for ( j=n=ne=0; j<x*y; ++j, ++k ) {
			if ( d[j] ) {
				ne += d[j];
				p->add(k, d[j]);
//				cout << ne << "\t" << n << endl;
				n = 0;
			} else if ( n==127 ) {
//				cout << ne << "\t" << n << endl;
				n = 1;
			} else {
				n++;
			}
		}
		nc = ecfd.getNElectronsCounted();
		dose = (nc-ncp)*1.0L/(x*y);
		if ( nimg > 1 ) {
			p->image[i].FOM(nc-ncp);
			(*p)["image"][i]["electrons"] = nc-ncp;
			(*p)["image"][i]["dose"] = dose;
			(*p)["image"][i]["exposure"] = exptime/f;
		}
//		cout << "Electrons counted: " << ne << "\t" << nc-ncp << "\t" << nc-ncp-ne << endl;
		ncp = nc;
	}
	
	dose = nc*1.0L/(p->images()*x*y);
	(*p)["electrons"] = nc;
	(*p)["dose"] = dose;
	(*p)["exposure"] = exptime;
	if ( nimg == 1 ) {
		p->image->FOM(nc);
		(*p)["image"][0]["electrons"] = nc;
		(*p)["image"][0]["dose"] = dose;
		(*p)["image"][0]["exposure"] = exptime;
	}

	if ( verbose ) {
		cout << "Electrons counted:              " << nc << endl;
		cout << "Electrons per pixel:            " << dose << endl << endl;
	}
	
//	for ( i=0; i<p->images(); ++i )
//		cout << i << tab << p->image[i].FOM() << endl;
	
	return 0;
}

vector<ElectronPos>	read_eer_positions(string filename)
{
	unsigned			x, y, z, f;
	ElectronCountedFramesDecompressor ecfd(filename);
	
	ecfd.getSize(x,y,z);
	f = ecfd.getNFrames();
	
	if ( verbose ) {
		cout << "Size: " << x << " " << y << " " << z << " " << f << endl;
//		cout << ecfd.GetAcquisitionMetadata() << endl;
	}

	JSvalue				mdata = json_from_xml(ecfd.GetAcquisitionMetadata().c_str());

	long				i, j, ne, nest, nep(0);
	ElectronPos			p0(0,0);
	vector<ElectronPos>	p;
	
	for ( i=0; i<f; ++i ) {
		nest = ecfd.nElectronFractionUpperLimit(i,i+1);
//		cout << "Estimated number of electrons: " << nest << endl;
		vector<ElectronPos>	p1(nest, p0);
		ecfd.decompressCoordinateList((ElectronPos*)p1.data(), i);
		ne = ecfd.getNElectronsCounted();
		j = ne - nep;
//		cout << "Electrons counted: " << j << endl;
//		cout << p1[0].x << " " << p1[0].y << endl;
//		cout << pL1[j-1].x << " " << p1[j-1].y << endl;
		nep = ne;
		p.insert(p.end(), p1.begin(), p1.begin()+j);
	}

	if ( verbose ) {
		cout << "Electrons counted:              " << ne << endl;
		cout << "Electrons per pixel:            " << ne/(1.0*x*y*f) << endl;
	}
	
	return p;
}

int			write_electron_list(string filename, vector<ElectronPos>& p)
{
	if ( verbose )
		cout << "Writing " << p.size() << " positions to " << filename << endl;
	
	ofstream		fpos(filename);
	if ( fpos.fail() )  return -1;
	
	fpos.write((char *)p.data(), p.size()*sizeof(ElectronPos));
	
	fpos.close();
	
	return 0;
}

