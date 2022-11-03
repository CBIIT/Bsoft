/**
@file	model_map.h
@brief	Header file for model to map conversions.
@author Bernard Heymann
@date	Created: 20081112
@date	Modified: 20220414
**/

#include "rwmodel.h"
#include "rwimg.h"
#include "ctf.h"
#include "Bgraphseg.h"

/* Function prototypes */
Bmodel*		model_from_images(Bimage* plist);
Bmodel*		model_from_graph_segments(Bimage* p, GSgraph& gs);
Bimage*		img_from_model(Bmodel* model, Vector3<double> ori,
				Vector3<long> size, Vector3<double> sam, double sigma);
int			model_catenate_maps(Bmodel* model, Bstring& filename);
int			model_shell_fit(Bmodel* model, double hires, double lores, int neg);
int			model_shell_radial_profile(Bmodel* model);
Bimage*		model_shell_power_spectrum(Bmodel* model, Vector3<long> size,
				Vector3<double> origin, int ft_size, int ann_min, int ann_max,
				double hires, double lores);
int			model_component_symmetry(Bmodel* model, long nangles,
				long ann_min, long ann_max, long ann_width, 
				long zmin, long zmax, long zinc, long minorder, long maxorder);
int			img_electron_scattering(Bmodel* model, Bimage* p,
				CTFparam& cp, double dose, double stdev, Bstring& atompropfile, int flag);

