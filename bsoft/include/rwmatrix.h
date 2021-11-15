/**
@file	rwmatrix.h 
@brief	Header file for reading (and writing) matrices.
@author Bernard Heymann 
@date	Created: 20010723
@date	Modified: 20120211
**/
 
#include "Bmatrix.h"

// Function prototypes
Bmatrix		read_matrix(Bstring& filename);
int 		write_matrix(Bstring& filename, Bmatrix matrix);

