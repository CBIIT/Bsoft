/**
@file	mg_ctf.h
@brief	Functions for CTF (contrast transfer function) processing
@author 	Bernard Heymann
@date	Created: 19970715
@date	Modified: 20230127
**/

#include "mg_processing.h"
#include "rwimg.h"

// Function prototypes
Bimage*		img_ctf_calculate(CTFparam& cp, bool flip, double wiener,
				Vector3<long> size, Vector3<double> sam, double lores, double hires);
Bimage*		img_ctf_calculate(CTFparam cp, int action, double wiener, Vector3<long> size,
				Vector3<double> sam, double lores=0, double hires=0);
Bimage*		img_ctf_gradient(CTFparam& cp, double def_min, double def_max, double def_inc,
				Vector3<long> size, Vector3<double> sam, double lores, double hires);
double		aberration(long n, long m, double s, double p);
vector<double>	aberration_terms(map<pair<long,long>,double>& wa, double u, double v);
vector<double>	aberration_terms(long nt, double u, double v);
vector<double>	aberration_even_terms(long nt, double u, double v);
vector<double>	aberration_odd_terms(long nt, double u, double v);
map<pair<long,long>,double>	aberration_weights(vector<double> v, int flag);
string		aberration_weight_string(map<pair<long,long>,double>& weights);
int			img_aberration_basis(Bimage* p, int flag);
int			img_create_aberration(Bimage* p, map<pair<long,long>,double>& weights, int flag);
int			img_create_aberration(Bimage* p, vector<map<pair<long,long>,double>>& weights, int flag);
int			img_create_aberration(Bimage* p, CTFparam cp, int flag);
int			img_create_aberration(Bimage* p, map<string,CTFparam>& cpa, int flag);
Bimage*		img_wave_aberration(CTFparam cp, Vector3<long> size, Vector3<double> sam);
int 		img_ctf_apply(Bimage* p, CTFparam em_ctf, int action, double wiener,
				double lores, double hires, bool invert);
int 		img_ctf_apply(Bimage* p, CTFparam em_ctf, int action, double wiener,
				double lores, double hires, bool invert, fft_plan planf_2D, fft_plan planb_2D);
int 		img_ctf_apply_complex(Bimage* p, CTFparam& cp, bool flip,
				double wiener, double lores, double hires);
int			img_apply_phase_aberration(Bimage* p, CTFparam em_ctf);
int			img_ttf_apply(Bimage* p, CTFparam ctf, int action, double wiener,
				Vector3<long> tile_size, double tilt, double axis, double res_lo, double res_hi, int invert);
int			img_ctf_apply_to_proj(Bimage* proj, CTFparam em_ctf, double defocus,
				double res_lo, double res_hi, bool invert, fft_plan planf_2D, fft_plan planb_2D);
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
int			project_set_envelope(Bproject* project, int type, vector<double>& coeff);
int			project_set_coherence_envelope(Bproject* project);
int			project_set_baseline_type(Bproject* project, int type);
int			project_set_baseline(Bproject* project, int type, vector<double>& coeff);
map<string,CTFparam>	project_ctf_optics_groups(Bproject* project);
long		project_update_ctf_aberration(Bproject* project, map<string,CTFparam>& cpa, int flag);
//int			project_convert_CTF_to_aberration_weights(Bproject* project);
int			project_aberration_compare(Bproject* project, Bproject* project2);
int			project_delete_aberration(Bproject* project, int which);
int			project_update_first_zero(Bproject* project);
int			project_plot_ctf(Bproject* project, Bstring& filename);
long		project_ctf_statistics(Bproject* project);

