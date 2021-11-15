/**
@file	mg_ctf.h
@brief	Functions for CTF (contrast transfer function) processing
@author Bernard Heymann
@date	Created: 19970715
@date	Modified: 20210817
**/

#include "mg_processing.h"
#include "rwimg.h"

// Function prototypes
Bimage*		img_ctf_calculate(CTFparam cp, int action, double wiener, Vector3<long> size, 
				Vector3<double> sam, double lores=0, double hires=0);
Bimage*		img_wave_aberration(CTFparam cp, Vector3<long> size, Vector3<double> sam);
int 		img_ctf_apply(Bimage* p, CTFparam em_ctf, int action, double wiener,
				double lores, double hires);
int 		img_ctf_apply(Bimage* p, CTFparam em_ctf, int action, double wiener,
				double lores, double hires, fft_plan planf_2D, fft_plan planb_2D);
int			img_ttf_apply(Bimage* p, CTFparam ctf, int action, double wiener,
				Vector3<long> tile_size, double tilt, double axis, double res_lo, double res_hi);
int			img_ctf_apply_to_proj(Bimage* proj, CTFparam em_ctf, double defocus, double res_lo, double res_hi, fft_plan planf_2D, fft_plan planb_2D);
int 		project_ctf_prepare(Bproject* project, int action, double lores,
				double hires, Vector3<long> tile_size, 
				double def_start, double def_end, double def_inc,
				Bstring& path, Bstring& newname, int flags);
int 		project_ctf(Bproject* project, int action, double lores,
				double hires, Vector3<long> tile_size, double wiener, 
				DataType datatype, Bstring& partpath, Bstring& newname, int flags);
int			project_powerspectrum_isotropy(Bproject* project, double lores, double hires);
JSvalue		project_defocus_range(Bproject* project);
int 		project_ctf_average(Bproject* project, Bstring& psname);
Bimage*		project_powerspectrum_average(Bproject* project, double deftarget);
int			project_merge_CTF_parameters(Bproject* project, Bproject* ctfproject);
int			project_CTF_to_part(Bproject* project);
int			project_set_defocus(Bproject* project, double def_avg,
					double def_dev, double ast_angle);
int			project_set_astigmatism(Bproject* project,
					double def_dev, double ast_angle);
int			project_update_ctf(Bproject* project, JSvalue& jsctf);
int			project_set_volts(Bproject* project, double volts);
int			project_set_Cs(Bproject* project, double Cs);
int			project_set_amp_shift(Bproject* project, double amp_shift);
int			project_set_focal_length(Bproject* project, double focal_length);
int			project_set_aperture(Bproject* project, double aperture);
int			project_set_slit_width(Bproject* project, double slit);
int			project_set_alpha(Bproject* project, double alpha);
int			project_set_envelope_type(Bproject* project, int type);
int			project_set_envelope(Bproject* project, int type, double* coeff);
int			project_set_coherence_envelope(Bproject* project);
int			project_set_baseline_type(Bproject* project, int type);
int			project_set_baseline(Bproject* project, int type, double* coeff);
int			project_update_first_zero(Bproject* project);
int			project_plot_ctf(Bproject* project, Bstring& filename);

