/**
@file	rwDM.h
@brief	Header file for reading and writing DM files
@author Bernard Heymann
@date	Created: 20020619
@date	Modified: 20150916

	Format: 2D electron microscope CCD image file format for the Digital Micrograph package
**/

#include "rwimg.h"

enum DMDataType {
	NULL_DATA, SIGNED_INT16_DATA, REAL4_DATA, COMPLEX8_DATA, OBSELETE_DATA,	// 4
	PACKED_DATA, UNSIGNED_INT8_DATA, SIGNED_INT32_DATA, RGB_DATA,			// 8
	SIGNED_INT8_DATA, UNSIGNED_INT16_DATA, UNSIGNED_INT32_DATA, REAL8_DATA,	// 12
	COMPLEX16_DATA, BINARY_DATA, RGBA_FLOAT32_DATA, RGB_UINT16_DATA ,		// 16
	RGB_FLOAT64_DATA, RGBA_FLOAT64_DATA, RGBA_UINT16_DATA,					// 19
	RGB_UINT8_DATA , RGBA_UINT8_DATA, LAST_DATA, OS_RGBA_UINT8_DATA			// 23
} ;

struct DMhead {             // file header for DM Fixed Format Specification
	unsigned int endian;	// if not 0x0000FFFF then byte swap (including  next 5 parameters).
	int xSize;
	int ySize;
	int zSize;
	int depth; 				// data type size in bytes - may be determined by data atype below
	DMDataType type; 		// An enumerated value
} ;
	
struct DMMachead {			// file header for old DM Macintosh Format Specification
	short width;
	short height;
	short bytes_per_pixel;
	short type;
} ;

/* The header information is stored in four short integers (16 bits each). 
Each integer is stored with the high byte first (Motorola or 'big-endian' format). 
The image data is stored row by row as a continuous stream of data (i.e. no line feeds,
carriage returns, or end of file control characters are present in the file). 
The amount of storage space that each pixel in the image occupies depends on the 
data type. The header and the image data are stored in the binary format 
(not in the ASCII format).

The format for the small header image file is described below.

Field			 Position  Length
Width of image   0  	   2 bytes (integer)
Height of image  2  	   2 bytes (integer)
Bytes per pixel  4  	   2 bytes (integer)
Encoded data type6  	   2 bytes (integer)
Image data  	 8  	   Width * height * bytes per pixel

The encoding of the data type (at position 6) is shown below.

Value Data Type
1	  2 byte signed integer
2	  Floating point IEEE standard 754
3	  8 byte complex floating point (real, imaginary)
5	  Packed complex
6	  1 byte unsigned integer
7	  4 byte signed integer
8	  4 byte RGB
9	  1 byte signed integer
10    2 byte unsigned integer
11    4 byte unsigned integer
12    8 byte double
13    16 byte complex double (real, imaginary)
14    1 byte binary data
*/

struct DM3head {		// All 3 fields must be big-endian
	int version;		// 2/3/4
	int file_length;	// Actually file length - 16 = size of root tag directory
	int endianness;		// 0=big, 1=little
	char sorted;		// Flag for sorted tag (1)
	char open;			// Flag for open tag (1)
	int	ntag;			// Number of tags
} ;

struct DM4head {		// All 3 fields must be big-endian
	int version;		// 2/3/4
	unsigned long file_length;	// Actually file length - 16 = size of root tag directory
	int endianness;		// 0=big, 1=little
	char sorted;		// Flag for sorted tag (1)
	char open;			// Flag for open tag (1)
	int	ntag;			// Number of tags
} ;

// I/O prototypes
int 		readDM(Bimage* p, int readdata, int img_select);
int 		writeDM(Bimage* p);
