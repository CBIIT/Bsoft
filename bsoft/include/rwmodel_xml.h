/**
@file	rwmodel_xml.h
@brief	Read and write models in XML format
@author Bernard Heymann
@date	Created: 20081029
@date	Modified: 20081029
**/

#include "rwmodel.h"

// Function prototypes
Bmodel*		read_model_xml(Bstring* file_list);
int			write_model_xml(Bstring& filename, Bmodel* model);

