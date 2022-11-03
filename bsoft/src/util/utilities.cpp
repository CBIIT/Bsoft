/**
@file	utilities.cpp
@brief	Library functions useful in all the package
@author Bernard Heymann
@date	Created: 19990722
@date	Modified: 20220722
**/

#include <errno.h>
#include <sys/sysctl.h>
#include <stdio.h>
#include "utilities.h"

#include <iostream>

#define MAX_GET_LEN	32

// Definition of global variables 
int 	verbose(0);		// Default: No output to shell
string	command;		// Command line
int		thread_limit(1000000);	// Default: Effectively no thread limit

ostream 	&tab(ostream &out)
{
    return out << '\t';
}

/**
@brief 	Reconstructs the command line as a global string.
@param 	argc		the number of command line arguments.
@param 	**argv 		the command line arguments.
@return int			0.

	Concatenates the command line arguments into one string.

**/
int			cmd_line(int argc, char **argv)
{ 
	int				i;
	
	command.clear();
	
	for ( i=0; i<argc; i++ )
		command = command + argv[i] + " ";

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG: cmd_line: " << command << endl;

	return 0;
} 

/**
@brief 	Returns the command line.
@return Bstring 		new string.

	This is designed to pack the command line into a string.

**/
Bstring 	command_line()
{
	return Bstring(command);
}

/**
@brief 	Returns the command line and time in a string.
@return Bstring 		new string.

	This is designed to pack the command line into a string followed by
	a second string for the time.

**/
Bstring 	command_line_time()
{
	time_t		ti = time(NULL);
	
	Bstring		clstring = "# " + command;
	
	clstring = clstring + "\n# " +  asctime(localtime(&ti)) + "\n";
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG command_line_time: " << clstring << endl;
	 
	return clstring;
}

string 		command_line_time2()
{
	time_t		ti = time(NULL);
	
	string		clstring = "# " + command;
	
	clstring = clstring + "\n# " +  asctime(localtime(&ti)) + "\n";
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG command_line_time: " << clstring << endl;
	 
	return clstring;
}

/**
@brief 	Returns the user name.
@return string 		user name.

	It uses getpwuid_r (thread safe) to find the user in the password file.

**/
string		get_user_name()
{
	int 		bufsize;

	if ( (bufsize = sysconf(_SC_GETPW_R_SIZE_MAX)) == -1 )
		bufsize = 16384;

	char 		buffer[bufsize];
	passwd 		pw, *result = NULL;

	string		user_name("unknown");

	if ( getpwuid_r(getuid(), &pw, buffer, bufsize, &result) == 0 && result)
		user_name = pw.pw_name;

	return user_name;
}

/**
@brief 	Prints usage information.
@param 	**use		the string array.
@param 	all			flag to output all usage information.

	The usage information must be written into an array of srings, with each
	string a line and following a specific convention for the Bsoft package.
	The first line with non-space characters must start with "Usage:" 
	followed by the command-line syntax. The next lines should describe the
	program. The options are indicated by lines strating with "-".
	The options are categorized as "Actions", "Parameters", "Input",
	and "Output". This constitutes the brief form.
	An additional section can be added as "Examples", that is shown
	only when the "all" argument is set.

**/
void		usage(const char** use, int all)
{
    int     i;

	if ( all )
		for ( i = 0; use[i] != NULL; i++ )
			cout << use[i] << endl;
	else
		for ( i = 0; use[i] && strstr(use[i], "Example") == 0; i++ )
			cout << use[i] << endl;
	
	cout << endl << "Bsoft version " << BVERSION << " (" << 8*sizeof(long) << " bit)" << endl << endl;
	
    bexit(0);
}

/**
@brief 	Finds the system type - mostly just for byte order.
@param 	show		a flag to indicate if the result should be shown.
@return SysType		an enumerated type.

	Test the byte order of an arbitrary byte sequence by interpreting it as
	an integer or a floating point number.

**/
SysType 	systype(int show)
{ 
	SysType 	type = Unknown_System; 
	char*		test = new char[8]; 
	int*		itest = (int *) test; 
	float*		ftest = (float *) test; 
	 
	memcpy(test, "jbh     ", 8); 
	 
	if ( itest[0] == 1784834080 && fabs(ftest[0] - 6.84272e+25) < 1e+21 ) 
		type = BigIEEE; 
	 
	if ( itest[0] == 543711850 && fabs(ftest[0] - 1.96837e-19) < 1e-23 ) 
		type = LittleIEEE; 
	 
	if ( ( ( verbose > VERB_PROCESS ) && show ) || ( verbose & VERB_DEBUG ) ) { 
		switch ( type ) { 
			case BigIEEE: 
				cout << "System type:                    Big-endian IEEE"; 
				break; 
			case LittleIEEE: 
				cout << "System type:                    Little-endian IEEE"; 
				break; 
			default: 
				cout << "System type:                    Unknown"; 
				break; 
		} 
		cout << " (" << 8*sizeof(long) << ")" << endl << endl; 
	} 
 
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG systype: " << test << " " << itest[0] << " " << ftest[0] << " " << type << endl; 
	 
	delete[] test;
	
	return type;
} 
 
/**
@brief 	Returns the number of processors.
@return long		number of processors.
**/
size_t		system_processors()
{
	size_t		np(1);
	
#ifdef HAVE_OMP
	np = omp_get_num_procs();
#endif
#ifdef HAVE_GCD
	np = sysconf( _SC_NPROCESSORS_ONLN );
#endif

	return np;
}

/**
@brief 	Returns system memory size.
@return long		memory size.
**/
long		system_memory()
{
#if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))

	int 		mib[2];
	mib[0] = CTL_HW;
#if defined(HW_MEMSIZE)
	mib[1] = HW_MEMSIZE;            // OSX
#elif defined(HW_PHYSMEM64)
	mib[1] = HW_PHYSMEM64;          // NetBSD, OpenBSD
#endif
	u_int 		namelen = sizeof(mib) / sizeof(mib[0]);
	long	 	size(0);
	size_t 		len = sizeof(size);

    if (sysctl(mib, namelen, &size, &len, NULL, 0) == 0)
		return size;

#elif defined(_SC_PHYS_PAGES)

#if defined(_SC_PAGESIZE)	// FreeBSD, Linux, OpenBSD, and Solaris
	return (long)sysconf( _SC_PHYS_PAGES ) *
		(long)sysconf( _SC_PAGESIZE );

#elif defined(_SC_PAGE_SIZE)
	return (long)sysconf( _SC_PHYS_PAGES ) *
		(long)sysconf( _SC_PAGE_SIZE );
#endif

#elif defined(_SC_AIX_REALMEM)	// AIX
	return (long)sysconf( _SC_AIX_REALMEM ) * (long)1024L;

#endif

	perror("system_memory: sysctl error");

    return 0;
}

/**
@brief 	Checking if there is enough memory, abort if not.
@param 	mem_required	memory required.
@return int				1 if enough.
**/
int			memory_check(long mem_required)
{
	long		mem_avail = system_memory();
	
	if ( verbose ) {
		cout << "Memory available:               " << mem_avail << endl;
		cout << "Memory required:                " << mem_required << " ("
			<< setprecision(2) << mem_required*100.0/mem_avail << " %)" << endl << endl;
	}
	
	if ( mem_avail < mem_required ) {
		cerr << "Error: Not enough memory!" << endl;
		cerr << "Memory available:               " << mem_avail << endl;
		cerr << "Memory required:                " << mem_required << endl << endl;
		bexit(-1);
	}
	
	return 1;
}

/**
@brief 	Exit function for cleanup and debugging.
@param 	value		exit value.
@return int			given value.
**/
int			bexit(int value)
{
	exit(value);
	return value;
}

/**
@brief 	Get the data type indicated by a single letter code.
@param 	letter 		letter indicating data type.
@return DataType 		data type.

	This function is used in optional command-line arguments to indicate 
	a new data type for an image.

**/
DataType	getdatatype(char letter)
{
	DataType		type;
	
	switch ( letter ) {
		case '1': type = Bit; break;
		case 'b': type = UCharacter; break;
		case 'c': type = SCharacter; break;
		case 'u': type = UShort; break;
		case 's': type = Short; break;
		case 'j': type = UInteger; break;
		case 'i': type = Integer; break;
		case 'k': type = ULong; break;
		case 'l': type = Long; break;
		case 'f': type = Float; break;
		case 'd': type = Double; break;
		case 'S': type = Short; break;
		case 'I': type = Integer; break;
		case 'F': type = Float; break;
		case 'D': type = Double; break;
		default: type = Unknown_Type; break;
	}
	
	return type;
}

/**
@brief 	Get the compound type indicated by a single letter code.
@param 	ch 				number of channels.
@param 	sct 				string indicating compound type.
@return CompoundType 		compound type.

	This function is used in optional command-line arguments to indicate
	a new compound type for an image.

**/
CompoundType	getcompoundtype(int ch, string sct)
{
	CompoundType		type(TSimple);
	
	transform(sct.begin(), sct.end(), sct.begin(), ::tolower);
	
	switch ( ch ) {
		case 1: type = TSimple; break;
		case 2:
			type = TComplex;		// Most common
			if ( sct[0] == 'v' ) type = TVector2;
			break;
		case 3:
			type = TRGB;			// Most common
			if ( sct[0] == 'v' ) type = TVector3;
			break;
		case 4: type = TView; break;
		default:
			type = TSimple;			// Most common
			if ( sct[0] == 'm' ) type = TMulti;
	}

	return type;
}


/**
@brief 	Converts a string with a selection specification into an integer array.
@param 	&string		string.
@param 	n			length of integer array.
@param 	*numsel		pre-allocated integer array.
@return int			number of levels.

	The integer array must be allocated to a length that would accommodate
	the highest number in the selection.
	If the string length is zero, all elements are selected.
	Multiple subsets are separated by colons

**/
int			select_numbers(Bstring& string, int n, int* numsel)
{
	vector<int>		ns = select_numbers(string, n);
	for ( long i=0; i<n; ++i ) numsel[i] = ns[i];
	return 0;
}

vector<int>	select_numbers(Bstring& string, int n)
{
	vector<int>		numsel(n, 0);

	int			i, j, k;
	
	if ( string.length() < 1 || string == "all" ) {
		for ( i=0; i<n; i++ ) numsel[i] = 1;
		return numsel;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG select_numbers: string=" << string << " n=" << n << endl;
	
	Bstring*	s, *s2, *sc;
	Bstring*	strcol = string.split(":");
	Bstring*	strarr = NULL;
	
	for ( k=1, sc = strcol; sc; sc = sc->next, k++ ) {
		strarr = sc->split(",");

		for ( s = strarr; s; s = s->next ) {
			s2 = s->split("-");
			i = j = s2->integer();
			if ( s2->next ) j = s2->next->integer();
			for ( ; i <= j && i < n; i++ ) numsel[i] = k;
			string_kill(s2);
		}
		
		string_kill(strarr);
	}

	string_kill(strcol);

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG select_numbers: ";
		for ( i=0; i<n; i++ ) cout << numsel[i];
		cout << endl;
	}
	
	return numsel;
}

/**
@author Dan Krainak
@brief 	Finds the next greatest number that is a power of a given number.
@param	startNumber	number to begin from.
@param	powerOf		power of this number is the number returned.
@return int 		the next greatest power (i.e. 128) or 0 on error.

	Loop through the powerOf variable, multiplying it each successive
	iteration until it is greater than the starting number.
	Eg., the next greatest power of 2 starting at 100 is 128.

**/
int 		findNextPowerOf(int startNumber, int powerOf)
{
	int	rVal = powerOf;
	
	if (powerOf <= 0 ) {
		rVal = 0; 
	} else if (powerOf == 1) {
		rVal = startNumber;
	} else {
		while (rVal < startNumber) {
			rVal *= powerOf;
		}
	} 

	return rVal;
}

/**
@brief 	Converts a defined length string into an integer.
@param 	*ptr		pointer to the string.
@param 	len			length to be scanned.
@return int 		the integer.

	The string is copied, 0-terminated, and scanned for an integer.

**/
int			get_integer(char* ptr, int len)
{
	if ( len > MAX_GET_LEN ) {
		cerr << "Error in get_integer: scanning length too large!" << endl;
		return 0;
	}
	char	temp[MAX_GET_LEN];
	int		d(0);
	strncpy(temp, ptr, len);
	temp[len] = 0;
	sscanf(temp, "%d", &d);
	return d;
}

/**
@brief 	Converts a defined length string into a floating point number.
@param 	*ptr		pointer to the string.
@param 	len			length to be scanned.
@return float 		the floating point number.

	The string is copied, 0-terminated, and scanned for a floating point number.

**/
float		get_float(char* ptr, int len)
{
	if ( len > MAX_GET_LEN ) {
		cerr << "Error in get_integer: scanning length too large!" << endl;
		return 0;
	}
	char	temp[MAX_GET_LEN];
	float	f(0);
	strncpy(temp, ptr, len);
	temp[len] = 0;
	sscanf(temp, "%f", &f);
	return f;
}

/**
@brief 	Swaps bytes.
@param	*v 			a pointer to the bytes.
@param 	n			number of bytes to swap.

	Byte swapping is done in place. 

**/
void		swapbytes(unsigned char* v, size_t n)
{
	unsigned char	t;
	size_t	i;
	
	for ( i=0, n--; i<n; i++, n-- ) {
		t = v[i];
		v[i] = v[n];
		v[n] = t;
	}
}

/**
@brief 	Swaps bytes.
@param 	size		size of the block to be swapped.
@param 	*v 			a pointer to the bytes.
@param 	n			number of bytes to swap.

	Byte swapping is done in place. 

**/
void		swapbytes(size_t size, unsigned char* v, size_t n)
{
	if ( n < 2 ) return;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG swapbytes: size=" << size << " n=" << n << endl;
	
	size_t	i;
	
	for ( i=0; i<size; i+=n, v+=n ) swapbytes(v, n);
}

/**
@brief 	Converts VAX floating point format to IEEE floating point format.
@param 	*v 			four-byte array holding the floating point value
@param 	sb			flag to swap bytes before conversion

	Swap bytes prior to conversion if the swap flag is set.
	Handle special cases of zero, infinity, NaN or normalized values
	Otherwise assign the new byte values
Reference:  Derived from CCP4 code

**/
void		vax2ieee(unsigned char* v, int sb) 
{ 
    unsigned char    t[4]; 
    int				i, shift, expnt; 
 
    if ( sb ) { 
		t[0] = v[2];  					// Swap byte order 
		t[1] = v[3]; 
		t[2] = v[0]; 
		t[3] = v[1]; 
    } else {
	    t[0] = v[1];  					// Swap for VAX format 
	    t[1] = v[0]; 
	    t[2] = v[3]; 
	    t[3] = v[2];
	}
	for ( i=0; i<4; i++ ) v[i] = t[i];

    expnt = (v[0] << 1) | (v[1] >> 7);	// Get exponent 
    if ( expnt == 0 ) 
    	if ( v[0] == 0 )				// Special case zero 
    	    memset(t, '\0', 4); 
    	else {			    			// Special case infinity or NaN 
	    	t[0] = 0xff; 
	    	t[1] = v[1] | 0x80; 
    	} 
    else if ( expnt > 2 )		    	// Normalized value 
    	t[0] = v[0] - 0x01;				// Subtract 2 from the exponent 
    else {								// Denormalized value 
    	shift = 3 - expnt; 
    	t[0] = v[0] & 0x80;	    		// Keep sign, expnt = 0  
    	t[1] = ( (v[1] & 0x7f) >> shift ) | ( 0x10 << expnt ); 
		t[2] = ( v[1] << ( 8-shift ) ) | ( v[2] >> shift ); 
		t[3] = ( v[2] << ( 8-shift ) ) | ( v[3] >> shift ); 
    } 
	for ( i=0; i<4; i++ ) v[i] = t[i];
}

/**
@brief 	Returns an angle between -PI and PI.
@param 	angle		input angle.
@return double		angle between -PI and PI.

**/
double		angle_set_negPI_to_PI(double angle)
{
	if ( !isfinite(angle) ) return 0;

	long 		n(angle/M_PI);
	
	angle -= (n+n%2)*M_PI;

	return angle;
}

/**
@brief 	Returns an angle between 0 and 2*PI.
@param 	angle		input angle.
@return double		angle between 0 and 2*PI.

**/
double		angle_set_0_to_TWOPI(double angle)
{
	if ( !isfinite(angle) ) return 0;

	long 		n(angle/TWOPI);
	
	angle -= n*TWOPI;
	if ( angle < 0 ) angle += TWOPI;

	return angle;
}

/**
@author   David Belnap
@brief 	Prevent a negative sign from being placed in front of zero value in
	a text file.
@param	value0      input value to be tested
@param	threshold   a small negative number
@return float   	"Corrected" or input value

	This function is intended to be used when obvious zero values are set
	to a very small negative number.  Input value is set to zero if
	              value0 > threshold  and  value0 < 0
	If so, the value is reset to zero.  Otherwise, the input value is 
	returned.
Reference:  Derived from CCP4 code

**/
float   remove_negative_zeros(float value0, float threshold)
{
	float  value;
	
	if ( (value0 > threshold) && (value0 < 0) )
		value = 0;
	else
		value = value0;

	return value;
}

/**
@brief 	Returns the chunk size per thread.
@param 	datasize	size of data to be divided into chunks.
@return long		chunk size.

	If multiple processors are used, the chunk size is set to the data size
	divided by the number of processors.
	Otherwise, the chunk size is equal to the data size.
Reference:  Derived from CCP4 code

**/
size_t		get_chunk_size(size_t datasize)
{
	size_t			np = system_processors();
	size_t			chunk_size = datasize/np;
	
	if ( chunk_size < datasize && chunk_size < MIN_SIZE_FOR_THREADING ) {
		chunk_size = MIN_SIZE_FOR_THREADING;
		np = (datasize - 1)/chunk_size + 1;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of processors:           " << np << endl << endl;
	
	if ( verbose & VERB_FULL )
		cout << "Chunk size:                     " << chunk_size << endl << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_chunk_size: chunk_size = " << chunk_size << endl;

	return chunk_size;
}

/**
@brief 	Returns the chunk size per thread.
@param 	datasize	size of data to be divided into chunks.
@param	c			channel size to align chunk to.
@return long		chunk size.

	If multiple processors are used, the chunk size is set to the data size
	divided by the number of processors.
	Otherwise, the chunk size is equal to the data size.
Reference:  Derived from CCP4 code

**/
size_t		get_chunk_size(size_t datasize, size_t c)
{
	size_t			np = system_processors();
	size_t			chunk_size(datasize/np);
	
	while ( chunk_size%c ) chunk_size++;
	
	if ( chunk_size < datasize && chunk_size < MIN_SIZE_FOR_THREADING ) {
		chunk_size = MIN_SIZE_FOR_THREADING;
		np = (datasize - 1)/chunk_size + 1;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of processors:           " << np << endl << endl;
	
	if ( verbose & VERB_FULL )
		cout << "Chunk size:                     " << chunk_size << endl << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_chunk_size: chunk_size = " << chunk_size << endl;

	return chunk_size;
}


/**
@brief 	Displays the error with file and line reference.
@param 	*message	a string to be included.
@param 	*file		the file name (should be __FILE__).
@param 	line		the line number (should be __LINE__).
@return int			error number.

	The function uses perror() to display a message containing the source
	file and line number where it originated.
Reference:  Derived from CCP4 code

**/
int			error_show(const char* message, const char* file, int line)
{
	char		string[1024];
	
	snprintf(string, 1024, "[%s:%d] %s", file, line, message);
	
	perror(string);
	
	if ( errno ) cerr << "errno " << errno << ": ";
	
	switch ( errno ) {
		case ENOENT:
			cerr << "Make sure the file exists and the path to the file is correct" << endl;
			break;
		case EACCES:
			cerr << "Make sure you have permission to write to this directory" << endl;
			break;
		case ENOMEM:
			cerr << "The memory available on this computer may be insufficient" << endl;
			break;
		default:
			cerr << endl;
	}
	
	return errno;
}

/**
@brief 	Displays the error with file and line reference.
@param 	*message	a string to be included.
@param 	*file		the file name (should be __FILE__).
@param 	line		the line number (should be __LINE__).
@return int			error number.

	The function uses perror() to display a message containing the source
	file and line number where it originated.
Reference:  Derived from CCP4 code

**/
int			error_show(string message, string file, int line)
{
//	char		string[1024];
	
//	snprintf(string, 1024, "[%s:%d] %s", file, line, message);
	
	string		s = "[" + file + ":" + to_string(line) + "] " + message;
	
	perror(s.c_str());
//	cout << s << endl;
	
	if ( errno ) cerr << "errno " << errno << ": ";
	
	switch ( errno ) {
		case ENOENT:
			cerr << "Make sure the file exists and the path to the file is correct" << endl;
			break;
		case EACCES:
			cerr << "Make sure you have permission to write to this directory" << endl;
			break;
		case ENOMEM:
			cerr << "The memory available on this computer may be insufficient" << endl;
			break;
		default:
			cerr << endl;
	}
	
	return errno;
}

