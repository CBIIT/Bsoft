/**
@file	utilities.h 
@brief	Header file for general utilities 
@author Bernard Heymann
@date	Created: 19990722
@date	Modified: 20210721
**/

#define BVERSION "2.1.3-20211110"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h> 
#include <pwd.h>
//#include <time.h>
#include <math.h>

#include <climits>
#include <cfloat>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <list>
#include <tuple>
#include <regex>
#include <ctime>

using namespace std;

#ifdef SUNOS
#include <ieeefp.h>
#endif

#ifdef HAVE_OMP
#include <omp.h>
#endif

#ifdef HAVE_GCD
#include <dispatch/dispatch.h>
#endif

#include "Bstring.h"

// Constants
#define MAXLINELEN			1024		// Maximum length of line for input
#define SMALLFLOAT			1e-30		// Threshold for considering a number to be zero
#define MIN_SIZE_FOR_THREADING 100000
#define TRIGPRECISION		1e-30

// Physical constants
#define PLANCK				6.626070E-34	// Js
#define	ECHARGE				1.602177E-19	// C
#define EMASS				9.109384E-31	// kg
#define LIGHTSPEED			299792458		// m/s
#define	AVOGADRO			6.02214076e23	// Entities per mole

/************************************************************************
@Constant: M_PI
@Description:
	Mathematical PI and variants.
*************************************************************************/
#ifndef M_PI 
#define M_PI	3.14159265358979323846264338327950288
#endif 
#ifndef M_PI_2 
#define M_PI_2	1.57079632679489661923132169163975144
#endif 
#ifndef M_PI_4 
#define M_PI_4	0.785398163397448309615660845819875721
#endif 
#ifndef TWOPI 
#define TWOPI	6.28318530717958647692528676655900576
#endif 
#ifndef MIN2PI 
#define MIN2PI	-6.28318530717958647692528676655900576
#endif
//end

/************************************************************************
@Constant: GOLDEN
@Description:
	Golden number: g = (sqrt(5) + 1)/2.
*************************************************************************/
#ifndef GOLDEN
#define GOLDEN	1.61803398874989484820458683436563811772
#endif
//end

/************************************************************************
@Constant: RHO
@Description:
	Protein density.
*************************************************************************/
#ifndef RHO 
#define RHO 		0.81    	// Protein density in Da/A3
#endif 
//end

/************************************************************************
@Constants: VERB_
@Description:
	Verbosity designations - controlling output to stdout.
@Features:
	The notion is that any program should not automatically generate output,
	so that it could be used within a script. Increasing levels of verbosity
	then results in increasing output.
*************************************************************************/
#define VERB_NONE		0
#define VERB_RESULT		1
#define VERB_LABEL		2
#define VERB_PROCESS	4
#define VERB_STATS		8
#define VERB_FULL		16
#define VERB_TIME		32
#define VERB_MEMORY 	64
#define VERB_DEBUG		128
#define VERB_DEBUG_STAR	256
#define VERB_DEBUG_XML	512
#define VERB_DEBUG_DM	1024
#define VERB_DEBUG_ND2	2048
#define VERB_DEBUG_EER	4096
//end

#ifndef _systype_ 
/************************************************************************
@Object: enum SysType
@Description:
	System type enumeration.
@Features:
	BigIEEE must be the lowest number type 
	LittleIEEE must be the little-endian type with the lowest number 
*************************************************************************/
enum SysType { 
	Unknown_System = 0,	// Indeterminate 
	BigIEEE = 1,		// Big-endian IEEE (unix: MIPS, Motorola PPC) 
	BigOther = 5,		// Big-endian systems other than IEEE 
	LittleIEEE = 6, 	// Little-endian IEEE (unix, Intel, VMS Alpha) 
	LittleAlpha = 7,	// Little-endian VMS non-IEEE (Alpha) 
	LittleVAX = 8,		// Little-endian VMS (VAX) 
	LittleOther = 10	// Little-endian other systems 
} ; 
#define _systype_ 
#endif 

#ifndef _datatype_
/**
@enum	DataType
@brief 	Base data type specifier.

	This determines what simple data type is used in an image.
**/
enum DataType {
	Unknown_Type = 0,	// Undefined data type
	Bit = 1,			// Bit/binary type
	UCharacter = 2,		// Unsigned character or byte type
	SCharacter = 3,		// Signed character
	UShort = 4,			// Unsigned integer (2-byte)
	Short = 5,			// Signed integer (2-byte)
	UInteger = 6,		// Unsigned integer (4-byte)
	Integer = 7,		// Signed integer (4-byte)
	ULong = 8,			// Unsigned integer (4 or 8 byte, depending on system)
	Long = 9,			// Signed integer (4 or 8 byte, depending on system)
	Float = 10,			// Floating point (4-byte)
	Double = 11,		// Double precision floating point (8-byte)
} ;

/**
@enum 	CompoundType
@brief 	Compound data type specifier.

	This determines what compound data type is used in an image.
**/
enum CompoundType {
	TSimple = 0,		// Single value data type
	TComplex = 1,		// 2-value complex type
	TVector2 = 2,		// 2-value vector type
	TVector3 = 3,		// 3-value vector type
	TView = 4,			// 4-value view type
	TRGB = 10,			// Red-green-blue interleaved
	TRGBA = 11,			// Red-green-blue-alpha interleaved
	TCMYK = 12,			// Cyan-magenta-yellow-black interleaved
	TMulti = 99			// Arbitrary number of channels
} ;

#define _datatype_
#endif 


/************************************************************************
@Constants: FILL_
@Description:
	Fill value designations - controlling the type of filling new pixels.
@Features:
	The fill value for pixels/voxels when new ones are created or old
	ones replaced.
*************************************************************************/
#define FILL_USER		0
#define FILL_AVERAGE	1
#define FILL_BACKGROUND 2
#define FILL_MIN		3
#define FILL_MAX		4
//end




// Function prototypes 
ostream 	&tab(ostream &out);
int			cmd_line(int argc, char **argv);
Bstring 	command_line();
Bstring 	command_line_time();
string 		command_line_time2();
string		get_user_name();
void		usage(const char** use, int all);
SysType 	systype(int show); 
size_t		system_processors();
long		system_memory();
int			memory_check(long mem_required);
int			bexit(int value);
DataType	getdatatype(char letter);
CompoundType	getcompoundtype(int ch, string sct);
int			select_numbers(Bstring& string, int n, int* numsel);
vector<int>	select_numbers(Bstring& string, int n);
int 		findNextPowerOf(int startNumber, int powerOf);
int			get_integer(char* ptr, int len);
float		get_float(char* ptr, int len);
void		swapbytes(unsigned char* v, size_t n);
void		swapbytes(size_t size, unsigned char* v, size_t n);
void		vax2ieee(unsigned char* v, int sb);
double		angle_set_negPI_to_PI(double angle);
float		remove_negative_zeros(float value0, float threshold);
size_t		get_chunk_size(size_t datasize);
size_t		get_chunk_size(size_t datasize, size_t c);
//int			error_show(const char* message, const char* file, int line);
int			error_show(string message, string file, int line);

