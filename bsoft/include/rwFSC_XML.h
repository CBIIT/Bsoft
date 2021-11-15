/**
@file	rwFSC_XML.h
@brief	Library routines to read and write FSC curves
@author Bernard Heymann
@date	Created: 20121216
@date	Modified: 20140423
**/

#include "ps_plot.h"

/// Function prototypes
Bplot*		xml_read_fsc(Bstring& filename, Bstring& xsdfile);
Bplot*		xml_read_fsc(Bstring& filename);
void		xml_write_fsc(Bstring& filename, Bplot* plot);

