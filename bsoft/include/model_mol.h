/**
@file	model_mol.h
@brief	Library routines for processing molecular models
@author Bernard Heymann
@date	Created: 20220215
@date	Modified: 20220324
**/

#include "rwmodel.h"
#include "json.h"

// Function prototypes
string		component_element(Bcomponent* comp);
string		component_element(Bcomponent* comp, map<string,Bcomptype>& atompar);
JSvalue		model_elements(Bmodel* model, map<string,Bcomptype>& atompar);
