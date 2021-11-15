/**
@file	model_select.h
@brief	Header file for reading and writing atomic model files
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20210318
**/

#include "rwmodel.h"
#include "Bimage.h"
#include "Bstring.h"

/* Function prototypes */
long 		model_select(Bmodel* model, Bstring& comptype);
long		model_select(Bmodel* model, long number);
long		model_select_all(Bmodel* model);
long		model_select_unknowns(Bmodel* model);
long		model_reset_selection(Bmodel* model);
long		model_unset_selection(Bmodel* model);
long		model_invert_selection(Bmodel* model);
long		model_select_sets(Bmodel* model, int size, int flag);
long		model_select_number_of_components(Bmodel* model, int ncomp_min, int ncomp_max);
long		model_select_closed(Bmodel* model, int closure_rule, int val_order);
long		model_select_fullerene(Bmodel* model);
long		model_select_non_fullerene(Bmodel* model);
long		model_select_valence(Bmodel* model, int valence);
long		model_select_polygons(Bmodel* model, int order);
long 		model_select_first(Bmodel* model, int first);
long 		model_select_within_shell(Bmodel* model, Vector3<double> center, double minrad, double maxrad);
long 		model_select_in_mask(Bmodel* model, Bimage* pmask);
long 		model_delete(Bmodel** model);
long 		model_delete_comp_type(Bmodel* model, Bstring& comptype);
long 		model_delete_non_selected(Bmodel** model);
long		model_selection_stats(Bmodel* model);
long		model_type_from_selection(Bmodel* model, Bstring* comp_type, Bstring& filename);
long		model_fom_deselect(Bmodel* model, double fom_cutoff);
long		model_fom_max_fraction_deselect(Bmodel* model, double fom_fraction);
//long		models_radius_deselect(Bmodel* model, double minrad, double maxrad);
long		model_fom_histogram(Bmodel* model, double fom_step);
long		model_fom_ranking(Bmodel* model, int nrank);
long		model_delete_overlapped_components(Bmodel** model, double distance);
long		model_average_overlapped_components(Bmodel* model, double distance);
long		model_prune_simple(Bmodel* model, double mindist);
long		model_prune_fom(Bmodel* model, double distance);
long		model_prune_similar(Bmodel* model);
long		model_prune_fit(Bmodel* model, double distance);
long		model_prune_large(Bmodel* model, double sampling);
long		model_find_overlap(Bmodel* model, Bstring& reffile, double distance);


