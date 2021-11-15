/**
@file	rwimg.cpp
@brief	Library for 2D and 3D image I/O
@author Bernard Heymann
@date	Created: 19990321
@date 	Modified: 20210628
**/

#include "rwimg.h"
#include "file_util.h"
#include "UnitCell.h"
#include "Complex.h"
#include "Vector3.h"
#include "utilities.h"

#include "rwRAW.h"
#include "rwASCII.h"
#include "rwBCR.h"
#include "rwBIORAD.h"
#include "rwBRIX.h"
#include "rwBrookhavenSTEM.h"
#include "rwCCP4.h"
#include "rwDI.h"
#include "rwDM.h"
#include "rwDSN6.h"
#include "rwDX.h"
#include "rwEER.h"
#include "rwEM.h"
#include "rwGOODFORD.h"
#include "rwGRD.h"
#include "rwHKL.h"
#include "rwIMAGIC.h"
#include "rwIP.h"
#include "rwJPEG.h"
#include "rwkernel.h"
#include "rwMFF.h"
#include "rwMIFF.h"
#include "rwMRC.h"
#include "rwND2.h"
#include "rwPIC.h"
#include "rwPIF.h"
#include "rwPNG.h"
#include "rwPNM.h"
#include "rwPostScript.h"
#include "rwSER.h"
#include "rwSitus.h"
#include "rwSPE.h"
#include "rwSPIDER.h"
#include "rwSUPRIM.h"
#include "rwTGA.h"
#include "rwTIFF.h"
#include "rwXPLOR.h"

// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern string	command;		// Command line

const char* month[] = {"Jan","Feb","Mar","Apr","May","Jun",
                    "Jul","Aug","Sep","Oct","Nov","Dec"} ;

// Internal function prototypes

/**
@brief	General driver function to read multiple image formats
@param 	filename	file name (plus any tags for the RAW format).
@param 	readdata	flag to activate reading of image data.
@param 	img_select	image selection in multi-image file (-1 = all images).
@return	Bimage*		the image structure, NULL if reading failed.
	This is the only image reading function that should be called from programs.
	A Bimage structure is initialized with default values.
	The file format is deduced from the file name extension.
	Every file format has its own funtion to read the file.
	The selection argument is used to read only one image from a multi-image
	file if it is greater than -1. This selection must be handled within 
	each format to ensure the correct allocation of the sub-image structure.
	If the requested selection is equal or larger than the number of
	images, the selection is set to the last image.
**/
Bimage* 	read_img(char* filename, int readdata, int img_select)
{
	string		fn(filename);
	return read_img(fn, readdata, img_select);
}

Bimage* 	read_img(Bstring filename, int readdata, int img_select)
{
	return read_img(filename.str(), readdata, img_select);
}

Bimage* 	read_img(string filename, int readdata, int img_select)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_img: Filename=" << filename << " dataflag=" << readdata << endl;
	
	if ( !filename.length() ) {
		cerr << "Error: No input image file given!" << endl;
		return NULL;
	}
	
	int 		err = 0;

//	Bstring		path;
	string		path;
//	Bstring		ext = filename.extension();	
	string		ext = extension(filename);	
//	Bstring		clean_filename(filename.pre('#'));
//	string		clean_filename(filename.pre('#').str());
	string		clean_filename(filename.substr(0, filename.find("#")));
	
//	if ( clean_filename.contains("@") ) {
//		img_select = clean_filename.post('@').integer();
//		clean_filename = clean_filename.pre('@');
//	}
	
	if ( clean_filename.find("@") != string::npos ) {
		img_select = stoi(clean_filename.substr(clean_filename.find("@")+1));
		clean_filename = clean_filename.substr(0, clean_filename.find("@"));
	}
	
//	if ( filename.contains("://") )
	if ( filename.find("://") != string::npos )
		clean_filename = find_file(clean_filename, path);
	else
//		clean_filename = clean_filename.pre(':');
		clean_filename = clean_filename.substr(0, clean_filename.find(":"));

//	if ( access(clean_filename.c_str(), R_OK) && ext.contains("krn") )
//		clean_filename = parameter_file_path(clean_filename);
	
	if ( access(clean_filename.c_str(), R_OK) && ext.find("krn") != string::npos )
		clean_filename = parameter_file_path(clean_filename);
	
	if ( access(clean_filename.c_str(), R_OK) ) {
		error_show(clean_filename.c_str(), __FILE__, __LINE__);
		return NULL;
	}

	Bimage* 	p = new Bimage;
	
//	if ( filename.contains("#") )
	if ( filename.find("#") != string::npos )
		p->file_name(filename);
	else
		p->file_name(clean_filename);
	
//	cerr << "file before reading: " << p->file_name() << endl;
//	cout << (*p)["filename"] << endl;
	
//	bexit(-1);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG read_img: Transferred filename=" << p->file_name() << endl;
		cout << "DEBUG read_img: Extension=" << ext << endl;
		cout << "DEBUG read_img: img_select=" << img_select << " readdata=" << readdata << endl;
	}

	if ( filename.find("#") != string::npos || ext.length() < 1 )
		err = readRAW(p, img_select);
	else if ( ext.find("raw") != string::npos )
		err = readRAW(p, img_select);
	else if ( ext.find("asc") != string::npos || ext.find( "txt") != string::npos )
		err = readASCII(p, readdata);
	else if ( ext.find("bcr") != string::npos )
		err = readBCR(p, readdata);
	else if ( ext.find("pic") != string::npos )
		err = readBIORAD(p, readdata);
	else if ( ext.find("brx" ) != string::npos || ext.find("brix" ) != string::npos )
		err = readBRIX(p, readdata);
	else if ( ext.find("dat" ) != string::npos )
		err = readBrookhavenSTEM(p, readdata);
	else if ( ext.find("ccp") != string::npos || ext.find("map") != string::npos )
		err = readCCP4(p, readdata);
	else if ( ext.find("di") != string::npos )
		err = readDI(p, readdata, img_select);
	else if ( ext.find("dm") != string::npos )
		err = readDM(p, readdata, img_select);
	else if ( ext.find("omap" ) != string::npos || ext.find("dsn6") != string::npos || ext.find("dn6")  != string::npos )
		err = readDSN6(p, readdata);
	else if ( ext.find("dx") != string::npos )
		err = readDX(p, readdata);
	else if ( ext.find("eer") != string::npos )
		err = readEER(p, readdata, img_select, -readdata);
	else if ( ext.find("em") != string::npos )
		err = readEM(p, readdata);
	else if ( ext.find("pot") != string::npos )
		err = readGOODFORD(p, readdata);
	else if ( ext.find("grd") != string::npos )
		err = readGRD(p, readdata, img_select);
	else if ( ext.find("hkl") != string::npos )
		err = readHKL(p, readdata);
	else if ( ext.find("img") != string::npos || ext.find("hed") != string::npos )
		err = readIMAGIC(p, readdata, img_select);
	else if ( ext.find("ip") != string::npos )
		err = readIP(p, readdata);
	else if ( ext.find("jpg") != string::npos || ext.find("jpeg") != string::npos )
		err = readJPEG(p, readdata);
	else if ( ext.find("krn") != string::npos )
		err = readKernel(p);
	else if ( ext.find("mif") != string::npos )
		err = readMIFF(p, readdata, img_select);
	else if ( ext.find("mff") != string::npos )
		err = readMFF(p, readdata);
	else if ( ext.find("mrc") != string::npos || ext == "stk" )
		err = readMRC(p, readdata, img_select);
	else if ( ext == "st" || ext.find("ali") != string::npos || ext.find("rec") != string::npos )	// IMOD extensions
		err = readMRC(p, readdata, img_select);
	else if ( ext.find("nd2") != string::npos )
		err = readND2(p, readdata, img_select);
	else if ( ext.find("pif") != string::npos || ext.find("sf") != string::npos )
		err = readPIF(p, readdata, img_select);
	else if ( ext.find("bp") != string::npos || ext.find("bq") != string::npos )
		err = readPIC(p, readdata);
	else if ( ext.find("png") != string::npos )
		err = readPNG(p, readdata);
	else if ( ext.find("pbm") != string::npos || ext.find("pgm") != string::npos || ext.find("ppm") != string::npos )
		err = readPNM(p, readdata, img_select);
	else if ( ext.find("ser") != string::npos )
		err = readSER(p, readdata);
	else if ( ext.find("sit") != string::npos )
		err = readSitus(p, readdata);
	else if ( ext.find("spe") != string::npos )
		err = readSPE(p, readdata, img_select);
	else if ( ext.find("spi") != string::npos )
		err = readSPIDER(p, readdata, img_select);
	else if ( ext.find("spm") != string::npos || ext.find("sup") != string::npos || ext == "f" )
		err = readSUPRIM(p, readdata);
	else if ( ext.find("tga") != string::npos )
		err = readTGA(p, readdata);
//	else if ( ext.find("tif") != string::npos || ext.find("eer") != string::npos )
	else if ( ext.find("tif") != string::npos )
		err = readTIFF(p, readdata, img_select);
	else if ( ext.find("xpl") != string::npos || ext.find("cns") != string::npos || ext.find("rfl") != string::npos )
		err = readXPLOR(p, readdata);
	else {
		cerr << "Error: File format with extension \"" << ext << "\" not supported!" << endl;
		err = -1;
	}

//	cerr << "file after reading: " << p->file_name() << endl;

	if ( err < 0 ) {
		error_show(filename, __FILE__, __LINE__);
		delete p;
		return NULL;
	}
	
//	write_img("test.pif", p);
	
	if ( p->fourier_type() > Standard ) img_convert_fourier(p, Standard);

	p->check();
	
	p->show_slice((long)p->image->origin()[2]);
//	if ( p->show_slice() < 0 ) p->show_slice(0);
	if ( p->show_slice() >= p->size()[2] ) p->show_slice(p->size()[2] - 1);
		
	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) ) {
		if ( img_select < 0 )
			cout << "Reading file:                   " << p->file_name() << endl;
		else
			cout << "Reading image " << img_select << " from file:      " << p->file_name() << endl;
	}
	
	if ( verbose & VERB_PROCESS ) {
		p->information();
		if ( !p->data_pointer() ) cout << "Data not read!" << endl << endl;
	}
	
	if ( verbose & VERB_STATS )
		p->subimage_information();
	
	return p;
}

/**
@brief	General driver function to write multiple image formats
@param 	filename		file name (plus any tags for the RAW format).
@param	*p				the image structure.
@param	compression		compression type: 0=none, 5=LZW(Tiff)
@return	int				error code (<0 means failure).
	This is the only image writing function that should be called
	from programs.
	The file format is deduced from the file name extension.
	Every file format has its own funtion to write the file.
**/
int 		write_img(const char* filename, Bimage* p, int compression)
{
	string		thefile(filename);
//	cout << "char*: length = " << thefile.length() << endl;
	return write_img(thefile, p, compression);
}

int 		write_img(Bstring filename, Bimage* p, int compression)
{
	return write_img(filename.str(), p, compression);
}

int 		write_img(string filename, Bimage* p, int compression)
{
	if ( !p ) return -1;
	
	if ( !filename.length() ) {
		cerr << "Error: No output image file name given!" << endl;
		return -1;
	}
	
	int				err(0), flags(compression);

//	Bstring			ext = filename.extension();
	string			ext = extension(filename);
//	cout << "extension = " << ext << endl;
		
	if ( ext.length() < 1 ) {
		cerr << "Error: No extension found in the file name " << filename << endl << endl;
		error_show(filename, __FILE__, __LINE__);
		return -2;
	}
	
	if ( p->compound_type() == TRGB || p->compound_type() == TRGBA ) {
		if ( ext.find("tif") == string::npos && ext.find("jp") == string::npos 
				&& ext.find("mif") == string::npos && ext.find("png") == string::npos 
				&& ext.find("pbm") == string::npos && ext.find("pgm") == string::npos 
				&& ext.find("ppm") == string::npos && ext.find("ps") == string::npos 
				&& ext.find("tga") == string::npos && ext.find("grd") == string::npos ) {
			cerr << "Warning: Only TIFF, JPEG, PNG, PBM, PGM, PPM, MIFF, TGA, GRD or postscript files are acceptable as RGB output formats!" << endl;
		}
	}
	
//	filename = filename.pre(':');
	
//	filename = filename.pre('#');

	string::size_type	pos = filename.find(":");
	if ( pos != string::npos ) filename = filename.substr(pos);
	
	pos = filename.find("#");
	if ( pos != string::npos ) filename = filename.substr(pos);

	p->file_name(filename);

	p->statistics();
	
    p->set_time(time(NULL));
	
	p->check();

	p->label(command);

	p->meta_data_update();
	(*p)["compression"] = compression;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG write_img: Image setup done for file " << p->file_name() << " (extension " << ext << ")" << endl;
		cout << p->meta_data() << endl;
	}

	if ( ext.find("raw") != string::npos )
		err = writeRAW(p);
	else if ( ext.find("asc") != string::npos || ext.find( "txt") != string::npos )
		err = writeASCII(p);
	else if ( ext.find("bcr") != string::npos )
		err = writeBCR(p);
	else if ( ext.find("pic") != string::npos )
		err = writeBIORAD(p);
	else if ( ext.find("brx" ) != string::npos || ext.find("brix" ) != string::npos )
		err = writeBRIX(p);
	else if ( ext.find("dat" ) != string::npos )
		err = writeBrookhavenSTEM(p);
	else if ( ( ext.find("ccp") != string::npos || ext.find("map") != string::npos ) )
		err = writeCCP4(p);
	else if ( ext.find("omap" ) != string::npos || ext.find("dsn6") != string::npos || ext.find("dn6") != string::npos )
		err = writeDSN6(p);
	else if ( ext.find("em") != string::npos )
		err = writeEM(p);
	else if ( ext.find("dm") != string::npos )
		err = writeDM(p);
	else if ( ext.find("dx") != string::npos )
		err = writeDX(p);
	else if ( ext.find("pot") != string::npos )
		err = writeGOODFORD(p);
	else if ( ext.find("grd") != string::npos )
		err = writeGRD(p, flags);
	else if ( ext.find("hkl") != string::npos )
		err = writeHKL(p);
	else if ( ext.find("img") != string::npos || ext.find("hed") != string::npos )
		err = writeIMAGIC(p);
	else if ( ext.find("ip") != string::npos )
		err = writeIP(p);
	else if ( ext.find("jpg") != string::npos || ext.find("jpeg") != string::npos )
		err = writeJPEG(p, flags);
	else if ( ext.find("krn") != string::npos )
		err = writeKernel(p);
	else if ( ext.find("mff") != string::npos )
		err = writeMFF(p);
	else if ( ext.find("mif") != string::npos )
		err = writeMIFF(p);
	else if ( ext.find("mrc") != string::npos || ext == "stk" )
		err = writeMRC(p);
	else if ( ext == "st" || ext.find("ali") != string::npos || ext.find("rec") != string::npos )	// IMOD extensions
		err = writeMRC(p);
	else if ( ext.find("pif") != string::npos || ext.find("sf") != string::npos )
		err = writePIF(p);
	else if ( ext.find("bp") != string::npos || ext.find("bq") != string::npos )
		err = writePIC(p);
	else if ( ext.find("png") != string::npos )
		err = writePNG(p);
	else if ( ext.find("pbm") != string::npos || ext.find("pgm") != string::npos || ext.find("ppm") != string::npos )
		err = writePNM(p);
	else if ( ext.find("ps") != string::npos )
		err = writePostScriptImage(p);
	else if ( ext.find("ser") != string::npos )
		err = writeSER(p);
	else if ( ext.find("sit") != string::npos )
		err = writeSitus(p);
	else if ( ext.find("spe") != string::npos )
		err = writeSPE(p);
	else if ( ext.find("spi") != string::npos )
		err = writeSPIDER(p);
	else if ( ext.find("spm") != string::npos || ext.find("sup") != string::npos || ext == "f" )
		err = writeSUPRIM(p);
	else if ( ext.find("tga") != string::npos )
		err = writeTGA(p);
	else if ( ext.find("tif") != string::npos ) {
		if ( p->sizeZ() > 1 ) flags |= 1;
		err = writeTIFF(p, flags);
	} else if ( ext.find("xpl") != string::npos || ext.find("cns") != string::npos || ext.find("rfl") != string::npos )
		err = writeXPLOR(p);
	else {
		cout << "File format with extension \"" << ext << "\" not supported!" << endl << endl;
		err = -1;
	}
	
	// Deallocate string
//	ext = 0;
	
	if ( err < 0 ) {
		error_show(filename, __FILE__, __LINE__);
		return err;
	}
	
	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) )
		cout << "Writing file:                   " << p->file_name() << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_img: p->data_pointer()=" << static_cast<void*>(p->data_pointer()) << endl;

	if ( verbose & VERB_PROCESS )
		p->information();
	
	if ( verbose & VERB_STATS )
		p->subimage_information();
	
	return err;
}

int			img_convert_fourier_one(Bimage* p, Complex<float>* nufdata, int n,
				long hermx, long friedel, long oldx, long nux,
				long xo1, long xo2, long yo1, long yo2,
				long zo1, long zo2, int z_even, int y_even)
{
    Complex<short>*		sdata = (Complex<short> *) p->data_pointer();
    Complex<int>*		idata = (Complex<int> *) p->data_pointer();
    Complex<float>*		fdata = (Complex<float> *) p->data_pointer();
    Complex<double>*	ddata = (Complex<double> *) p->data_pointer();
	
	long				x, y, z, nz, i, j, izi, izf, iyi, iyf, jz, jy;
	long				xf, yf, zf, xi, yi, zi;
	int					phisgn = 1;
    Complex<double>*	nuddata = (Complex<double> *) nufdata;
    Complex<int>*		nuidata = (Complex<int> *) nufdata;
    Complex<short>*		nusdata = (Complex<short> *) nufdata;
	
	nz = n*p->sizeZ();
	for ( z=0; z<p->sizeZ(); z++ ) {
		zi = z + zo1 - zo2;
		if ( zi < 0 ) zi += p->sizeZ();
		if ( zi >= (int)p->sizeZ() ) zi -= p->sizeZ();
		if ( z_even )
			zf = (zi)? p->sizeZ()-zi: zi;
		else
			zf = 2*zo1 - zi;
		if ( zf < 0 ) zf += p->sizeZ();
		izi = (nz + zi)*p->sizeY();
		izf = (nz + zf)*p->sizeY();
		jz = (nz + z)*p->sizeY();
		for ( y=0; y<p->sizeY(); y++ ) {
			yi = y + yo1 - yo2;
			if ( yi < 0 ) yi += p->sizeY();
			if ( yi >= (int)p->sizeY() ) yi -= p->sizeY();
			if ( y_even )
				yf = (yi)? p->sizeY()-yi: yi;
			else
				yf = 2*yo1 - yi;
			if ( yf < 0 ) yf += p->sizeY();
			iyi = (izi + yi)*oldx;
			iyf = (izf + yf)*oldx;
			jy = (jz + y)*nux;
			for ( x=0; x<nux; x++ ) { //cout << "%d %d %d\n", x, y, z);
				xi = x + xo1 - xo2;
				if ( xi < 0 ) xi += p->sizeX();
				if ( xi >= (int)p->sizeX() ) xi -= p->sizeX();
				xf = (xi)? p->sizeX()-xi: xi;
				if ( !friedel || ( xi < (int)hermx ) ) {
					i = iyi + xi;
					phisgn = 1;
				} else {
					i = iyf + xf;
					phisgn = -1;
				}
				j = jy + x;
				if ( p->data_type() == Short ) {
					if ( phisgn > 0 ) nusdata[j] = sdata[i];
					else nusdata[j] = sdata[i].conj();
				} else if ( p->data_type() == Integer ) {
					if ( phisgn > 0 ) nuidata[j] = idata[i];
					else nuidata[j] = idata[i].conj();
				} else if ( p->data_type() == Float ) {
					if ( phisgn > 0 ) nufdata[j] = fdata[i];
					else nufdata[j] = fdata[i].conj();
				} else if ( p->data_type() == Double ) {
					if ( phisgn > 0 ) nuddata[j] = ddata[i];
					else nuddata[j] = ddata[i].conj();
				}
			}
		}
	}

	return 0;
}

/**
@brief	Converts Fourier transform types.
@param 	*p			the image.
@param	nutransform	new transform type.
@return	int			error code (<0 means failure).
Fourier transform classification:
	0=NoTransform:	No transform: Just a complex data set
	1=Standard:		Standard transform with origin = (0,0,0)
						(Suprim)
	2=Centered:		Centered transform with origin = (nx/2,ny/2,nz/2)
						(Imagic)
	3=Hermitian:	Hermitian transform with origin = (0,0,0) and size (nx/2+1,ny,nz)
						(Spider, EM)
	4=CentHerm:		Centered hermitian transform with origin = (0,ny/2,nz/2)
					and size (nx/2+1,ny,nz)
						(MRC)
	Assumption: 	The correct dimensions for a standard transform is stored
					in the x, y, and z fields.
**/
int			img_convert_fourier(Bimage* p, FourierType nutransform)
{
	if ( !p ) return -1;
	
	if ( !p->data_pointer() ) return -1;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG img_convert_fourier: old transform type = " << p->fourier_type() << endl;
		cout << "DEBUG img_convert_fourier: new transform type = " << nutransform << endl;
		cout << "DEBUG img_convert_fourier: input size = " << p->size() << endl;
	}
	
	if ( p->compound_type() != TComplex ) {
		p->fourier_type(NoTransform);
		return 0;
	}
	
	if ( p->fourier_type() == NoTransform ) p->fourier_type(Standard);
	
	if ( nutransform == p->fourier_type() ) return 0;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Converting Fourier transform:   ";
		switch ( p->fourier_type() ) {
			case Standard: cout << "Standard -> "; break;
			case Centered: cout << "Centered -> "; break;
			case Hermitian: cout << "Hermitian -> "; break;
			case CentHerm: cout << "Centered hermitian -> "; break;
			default: break;
		}
		switch ( nutransform ) {
			case Standard: cout << "Standard"; break;
			case Centered: cout << "Centered"; break;
			case Hermitian: cout << "Hermitian"; break;
			case CentHerm: cout << "Centered hermitian"; break;
			default: break;
		}
		cout << endl << endl;
	}
		
	long	hermx = p->sizeX()/2 + 1;
	long	friedel = 0;
	long	oldx = p->sizeX(), nux = p->sizeX();
	if ( p->fourier_type() == Hermitian || p->fourier_type() == CentHerm ) {
		oldx = hermx;
		if ( nutransform == Standard || nutransform == Centered )
			friedel = 1;		// Flag => fill in friedel pairs
	}
	if ( nutransform == Hermitian || nutransform == CentHerm )
		nux = hermx;
	
	long	xo1, yo1, zo1, xo2, yo2, zo2;	// Origins
	switch ( p->fourier_type() ) {
		case Centered:
			xo1 = p->sizeX()/2;
			yo1 = p->sizeY()/2;
			zo1 = p->sizeZ()/2;
			break;
		case CentHerm:
			xo1 = 0;
			yo1 = p->sizeY()/2;
			zo1 = p->sizeZ()/2;
			break;
		default:
			xo1 = yo1 = zo1 = 0;
			break;
	}
	
	switch ( nutransform ) {
		case Centered:
			xo2 = p->sizeX()/2;
			yo2 = p->sizeY()/2;
			zo2 = p->sizeZ()/2;
			break;
		case CentHerm:
			xo2 = 0;
			yo2 = p->sizeY()/2;
			zo2 = p->sizeZ()/2;
			break;
		default:
			xo2 = yo2 = zo2 = 0;
			break;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG img_convert_fourier: hermx=" << hermx << " oldx=" << oldx << " nux=" << nux << endl;
		cout << "DEBUG img_convert_fourier: input origin={" << xo1 << "," << yo1 << "," << zo1 << endl;
		cout << "DEBUG img_convert_fourier: output origin={" << xo2 << "," << yo2 << "," << zo2 << endl;
	}
	
	long				datasize = p->alloc_size();
    Complex<float>*		nudata = new Complex<float>[datasize];
	
	int					z_even = 0, y_even = 0;
	if ( p->sizeZ()%2 == 0 ) z_even = 1;
	if ( p->sizeY()%2 == 0 ) y_even = 1;

#ifdef HAVE_GCD
	dispatch_apply(p->images(), dispatch_get_global_queue(0, 0), ^(size_t n){
		img_convert_fourier_one(p, nudata, n, hermx, friedel, oldx, nux,
				xo1, xo2, yo1, yo2, zo1, zo2, z_even, y_even);
	});
#else
#pragma omp parallel for
	for ( long n=0; n<p->images(); n++ )
		img_convert_fourier_one(p, nudata, n, hermx, friedel, oldx, nux,
				xo1, xo2, yo1, yo2, zo1, zo2, z_even, y_even);
#endif

	p->data_assign((unsigned char *) nudata);
	
	p->fourier_type(nutransform);
	
	return 0;
}

