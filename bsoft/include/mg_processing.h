/**
@file	mg_processing.h
@brief	Header file for micrograph processing
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20210728
**/

#include "ctf.h"
#include "symmetry.h"
#include "marker.h"
#include "View.h"
#include "Euler.h"
#include "Bstring.h"

#define	NFOM 10

#define	MODE_PPS		0
#define	MODE_SCC		1
#define	MODE_CCC		2
#define	MODE			3
#define	FULL_ASU		8
#define	MULTI_FILE		16
#define	APPLY_CTF		32
#define	PART_LOG		64
#define	WRITE_PPX		128
#define	CHECK_PPX		256


#ifndef _ProjectParamStructs_
/************************************************************************
@Object: enum FOMType
@Description:
	Particle figure-of-merit type.
@Features:
	A variety of FOM's that can be associated with a comparison of images.
*************************************************************************/
enum FOMType {
	NoFOM,
	FOM,
	FOM_CC,
	FOM_CV,
	FOM_SNR,
	FOM_CC_AVG,
	FOM_CC_STD,
	FOM_HAND_A,
	FOM_HAND_B,
	FOM_PFT_CC,
	FOM_PFT_PRJ,
	FOM_PFT_CMP,
	FOM_RFACTORAB,
	COVERAGE,
	DENSITY,
	FOM1,
	FOM2,
	FOM3,
	FOM4,
	FOM5,
	FOM6,
	FOM7,
	FOM8,
	FOM9,
	FOMlast
} ;

class Bmicrograph;
class Breconstruction;

/************************************************************************
@Object: class Bframe
@Description:
	Movie frame parameter structure.
@Features:
	Shift for each frame.
*************************************************************************/
class Bframe {
private:
	void	initialize() {
		next = NULL;
		id = 0;
		fom = 1;
		sel = 1;
	}
public:
	Bframe*			next;			// Next frame in list
	int				id;				// Frame number in file (starts at 1)
	Vector3<double>	shift;			// Shift relative to reference frame
	double			fom;			// Figure-of-merit
	long			sel;			// Selection flag
	Bframe() { initialize(); }
} ;

/************************************************************************
@Object: class Bparticle
@Description:
	Single particle image parameter structure.
@Features:
	Orientation parameters for each particle.
*************************************************************************/
class Bparticle {
private:
	void	initialize() {
		next = NULL;
		fpart = 0;
		id = 0;
		group = 1;
		mag = 1;
		def = 0;
		dev = 0;
		ast = 0;
		view = View(0,0,1,0);
		for ( int i=0; i<NFOM; ++i ) fom[i] = 0;
		fom[0] = 1;
		sel = 1;
		mg = NULL;
		rec = NULL;
	}
public:
	Bparticle*		next;			// Next particle in list
	Bstring			fpart;			// Particle image file name
	int				id;				// Particle number in file (starts at 1)
	int				group;			// Group membership (such as filaments or crystals)
	double			mag;			// Magnification
	double			def;			// Defocus average (angstrom)
	double			dev;			// Defocus deviation (angstrom)
	double			ast;			// Astigmatism angle
	Vector3<double>	pixel_size;		// Particle pixel size - usually same as micrograph
	Vector3<double>	loc;			// Coordinates in the micrograph or tomogram
	Vector3<double>	ori;			// Origin of particle in voxel units
	View			view; 			// View: 3-value vector and angle (radians)
	double			fom[NFOM];		// Figure-of-merit, types defined in project structure
	long			sel;			// Selection flag
	Bmicrograph*	mg;
	Breconstruction*	rec;
	Bparticle() { initialize(); }
} ;

/************************************************************************
@Object: class Bfilnode
@Description:
	Filament node parameter structure.
@Features:
	Parameters defining a filament node.
*************************************************************************/
class Bfilnode {
private:
	void	initialize() {
		next = NULL;
		id = 0;
	}
public:
	Bfilnode*		next;			// Next node in list
	int				id;				// Filament node number (starts at 1)
	Vector3<double>	loc;			// Coordinates in the micrograph or tomogram
	Bfilnode() { initialize(); }
} ;

/************************************************************************
@Object: class Bfilament
@Description:
	Filament parameter structure.
@Features:
	Parameters defining a filament.
*************************************************************************/
class Bfilament {
private:
	void	initialize() {
		next = NULL;
		ffil = 0;
		id = 0;
		node = NULL;
	}
public:
	Bfilament*		next;			// Next filament in list
	Bstring			ffil;			// Filament image file name
	int				id;				// Filament number (starts at 1)
	Bfilnode*		node;			// Linked list of nodes
	Bfilament() { initialize(); }
} ;

/************************************************************************
@Object: class Bbadarea
@Description:
	Structure to mark a bad area.
@Features:
	3D coordinates for the center of the bad area.
*************************************************************************/
class Bbadarea {
private:
	void	initialize() {
		next = NULL;
		id = 0;
	}
public:
	Bbadarea*		next;			// Next bad area in list
	int				id;				// Bad area number
	Vector3<double>	loc;			// Coordinates in the micrograph or tomogram
	Bbadarea() { initialize(); }
} ;

/************************************************************************
@Object: class Bstrucfac
@Description:
	Structure factor structure.
@Features:
	3D Miller indices and structure factor locations.
*************************************************************************/
class Bstrucfac {
private:
	void	initialize() {
		next = NULL;
		amp = 0;
		sigamp = 0;
		phi = 0;
		sigphi = 0;
		fom = 0;
		sel = 1;
	}
public:
	Bstrucfac*		next;			// Next structure factor in list
	Vector3<double>	loc;			// Coordinates in the transform
	Vector3<int>	index;			// Miller indices
	double			amp;			// Amplitude
	double			sigamp;			// Amplitude deviation
	double			phi;			// Phase
	double			sigphi;			// Phase deviation
	double			fom;			// Figure-of-merit
	long			sel;			// Selection flag
	Bstrucfac() { initialize(); }
} ;

/************************************************************************
@Object: class Blayerline
@Description:
	Layer line structure.
@Features:
	3D Miller indices and structure factor locations.
*************************************************************************/
class Blayerline {
private:
	void	initialize() {
		next = NULL;
		number = 0;
		order = 0;
		distance = 0;
		freq = 0;
		amp = 0;
		fom = 0;
		sel = 1;
	}
public:
	Blayerline*		next;			// Next layer line in list
	int				number;			// Layer line number
	int				order;			// Bessel order
	double			distance;		// Distance along helical axis from origin
	double			freq;			// Spatial frequency
	double			amp;			// Amplitude
	double			fom;			// Figure-of-merit
	long			sel;			// Selection flag
	Blayerline() { initialize(); }
} ;

/************************************************************************
@Object: class Bmicrograph
@Description:
	General micrograph parameter structure.
@Features:
	This contains all parameters associated with a micrograph:
		The micrograph number and file name.
		The associated STAR data block number.
		The picked particle file name and parameters.
		The filament file name and parameters.
		Defocus parameters for the micrograph.
		Orientation parameters for each particle.
		Filament node locations.
		Marker locations.
*************************************************************************/
class Bmicrograph {
private:
	void	initialize() {
		next = NULL;
		id = 0;
		select = 1;
		block = 0;
		fmg = 0;
		fframe = 0;
		fpart = 0;
		ffil = 0;
		fft = 0;
		fps = 0;
		img_num = 0;
		magnification = 0;
		sampling = 0;
		pixel_size = Vector3<double>(1,1,1);
		dose = 0;
		intensity = 0;
		wri = 0;
		tilt_axis = 0;
		tilt_angle = 0;
		level_angle = 0;
		rot_angle = 0;
		scale = Vector3<double>(1,1,1);
		matrix = Matrix3(1);
		helix_axis = 0;
		helix_rise = 0;
		helix_angle = 0;
		helix_radius = 0;
		filament_width = 0;
		fil_node_radius = 0;
		bad_radius = 0;
		sf_radius = 0;
		mark_radius = 0;
		fom = 0;
		frame = NULL;
		ctf = NULL;
		part = NULL;
		fil = NULL;
		bad = NULL;
		mark = NULL;
		sf = NULL;
		layer = NULL;
	}
public:
	Bmicrograph*	next;			// Next micrograph in list
	Bstring			id;				// Micrograph identifier
	long			select;			// Selection flag
	int				block;			// STAR data block number
	Bstring			fmg;			// Micrograph image file
	Bstring			fframe;			// Micrograph image frames file
	Bstring			fpart;			// Image file with picked particles
	Bstring			ffil;			// Image file with extracted filaments
	Bstring			fft;			// Image file with Fourier transform
	Bstring			fps;			// Image file with power spectrum
	int 			img_num;		// Image number in file
	double			magnification;	// Microscope magnification
	double			sampling;		// Scanner sampling (angstrom)
	Vector3<double>	frame_pixel_size; 	// Nominal frame pixel size
	Vector3<double>	pixel_size; 	// Nominal micrograph pixel size
	double			exposure; 		// Acquisition time (seconds)
	double			dose; 			// Electron dose (electrons/angstrom^2)
	double			intensity;		// Micrograph average
	double			wri;			// Water ring index
	double			tilt_axis;		// Tilt axis angle, origin at x-axis (radians)
	double			tilt_angle; 	// Tilt angle, right-handed around tilt axis (radians)
	double			level_angle; 	// Level angle, deviation of tilt axis from xy plane (radians)
	double			rot_angle;		// In-plane rotation angle of micrograph or specimen (radians)
	Vector3<double>	origin;			// Origin of micrograph (usually the center)
	Vector3<double>	scale;			// Scale with respect to field-of-view
	Matrix3			matrix;			// Affine matrix for non-rigid transformations
	Vector3<double>	hvec;			// Vector for first Miller index
	Vector3<double>	kvec;			// Vector for second Miller index
	Vector3<double>	lvec;			// Vector for third Miller index
	double			helix_axis;		// Rotation angle defining the helical axis
	double			helix_rise;		// Helical rise per subunit
	double			helix_angle;	// Helical rotation angle per subunit
	double			helix_radius;	// Helical radius
	Vector3<long>	box_size;		// Particle box size
	double			filament_width;	// Filament width
	double			fil_node_radius; // Filament node radius
	double			bad_radius;		// Radius of bad area around coordinates
	double			sf_radius;		// Structure factor radius
	double			mark_radius;	// Radius of marker
	double			fom;			// Figure-of-merit for the micrograph
	Bframe*			frame;			// Movie frame parameters
	CTFparam*		ctf;			// Contrast transfer function parameters
	Bparticle*		part;			// First particle in linked list
	Bfilament*		fil;			// First filament in linked list
	Bbadarea*		bad;			// First bad area in linked list
	Bmarker*		mark;			// First marker in linked list
	Bstrucfac*		sf;				// First structure factor in linked list
	Blayerline*		layer;			// First layer line in linked list
	Bmicrograph() { initialize(); }
} ;

/************************************************************************
@Object: class Bfield
@Description:
	Field-of-view parameter structure.
@Features:
	One to many micrographs can be taken of a field-of-view, including
	focal and tomographing series.
	The orientation parameters, origin and matrix, are intended to deal 
	with tomographic dual tilt series.
*************************************************************************/
class Bfield {
private:
	void	initialize() {
		next = NULL;
		id = 0;
		matrix = Matrix3(1);
		select = 1;
		fom = 0;
		mg = NULL;
	}
public:
	Bfield* 		next;			// Next field-of-view in list
	Bstring			id;				// Field identifier
	Vector3<double>	origin;			// Origin of reference micrograph in the field (usually the center)
	Matrix3			matrix;			// Affine matrix for non-rigid transformations
	long			select;			// Selection flag
	double			fom;			// Figure-of-merit for the field
	Bmicrograph*	mg;				// First micrograph in list
	Bfield() { initialize(); }
} ;

/************************************************************************
@Object: class Breconstruction
@Description:
	Reconstruction parameter structure.
@Features:
	This contains all parameters associated with a reconstruction:
		The reconstruction number and file name.
		The associated STAR data block number.
		The picked particle file name and parameters.
		Orientation parameters for each particle.
		Filament node locations.
		Marker locations.
*************************************************************************/
class Breconstruction {
private:
	void	initialize() {
		next = NULL;
		id = 0;
		select = 1;
		block = 0;
		type = 0;
		symmetry = 0;
		frec = 0;
		fpart = 0;
		ffil = 0;
		fft = 0;
		fps = 0;
		voxel_size = Vector3<double>(1,1,1);
		scale = Vector3<double>(1,1,1);
		view = View(0,0,1,0);
		filament_width = 0;
		fil_node_radius = 0;
		bad_radius = 0;
		sf_radius = 0;
		mark_radius = 0;
		fom = 0;
		ctf = NULL;
		part = NULL;
		fil = NULL;
		bad = NULL;
		mark = NULL;
		sf = NULL;
		model = NULL;
	}
public:
	Breconstruction* 	next;			// Next reconstruction in list
	Bstring				id;				// Reconstruction identifier
	long				select;			// Selection flag
	int					block;			// STAR data block number
	int					type;			// Reconstruction type: 1-99 = spa; 100-199 = fil/hel; 200-299 = xtal; 300+ = tomo
	Bstring				symmetry;		// Symmetry string
	Bstring				frec;			// Reconstruction image file name
	Bstring				fpart;			// Particle file name
	Bstring				ffil;			// Filament file name
	Bstring				fft;			// Fourier transform file name
	Bstring				fps;			// Power spectrum file name
	Vector3<double>		voxel_size; 	// Nominal reconstruction voxel size
	Vector3<double>		origin;			// Image origin (voxels)
	Vector3<double>		scale;			// Scale
	View				view; 			// View: 3-value vector and angle (radians)
	Vector3<long>		box_size;		// Particle box size
	double				filament_width;	// Filament width
	double				fil_node_radius; // Filament node radius
	double				bad_radius;		// Radius of bad area around coordinates
	double				sf_radius;		// Structure factor radius
	double				mark_radius;	// Radius of marker
	double				fom;			// Figure-of-merit for the reconstruction
	CTFparam*			ctf;			// Contrast transfer function parameters
	Bparticle*			part;			// First particle in linked list
	Bfilament*			fil;			// First filament in linked list
	Bbadarea*			bad;			// First bad area in linked list
	Bmarker*			mark;			// First marker in linked list
	Bstrucfac*			sf;				// First structure factor in linked list
	Bstring*			model;			// First atomic model name in linked list
	Breconstruction() { initialize(); }
} ;

/************************************************************************
@Object: class Bproject
@Description:
	Project parameter structure.
@Features:
	A project is defined as a set of micrographs taken of a single 
	specimen and the associated reconstructions.
	A specimen is understood to be the result of a particular
	preparation method.
	This structure is the top level in a hierarchy describing
	micrographs and relationships between them.
	Orientation parameters can be given as view vector and rotation angle
	around the view vector (default).
	Alternatively, orientation can be described by Euler angles.
	Only one orientation representation is allowed to avoid ambiguity,
	and the choice is indicated by the "Euler" flag.
	Older software relying on the negative of the third Euler angle, -psi
	or omega, can make use of the omega flag.
	Multiple FOM's are allowed and the tags provide program-specific
	descriptions of the FOM's in the particle structure
*************************************************************************/
class Bproject {
private:
	void	initialize() {
		next = NULL;
		split = 0;
		filename = 0;
		comment = 0;
		reference = NULL;
		select = 0;
		euler_flag = 0;
		omega_flag = 0;
		field = NULL;
		rec = NULL;
		class_avg = NULL;
		for ( int f=0; f<NFOM; ++f ) fom_tag[f] = NoFOM;
		fom_tag[0] = FOM;
	}
public:
	Bproject*			next;			// Pointer to the next structure in the list
	Bstring				filename;		// Parameter file name
	Bstring				comment;		// Variable length comment field
	int					split;			// Split code, 0=no split, 1-4=number of digits in file name, 9=use mg id's
	Bstring*			reference;		// Linked list of reference map file names
	int					select;			// Selects micrographs (0) or reconstructions (1)
	int					euler_flag; 	// Set if Euler angles are used
	int					omega_flag; 	// Set if omega is used in stead of psi
	FOMType				fom_tag[NFOM];	// STAR tags for up to NFOM FOM values
	Bfield*				field;			// First field-of-view in list
	Breconstruction*	rec;			// First reconstruction in list
	Bparticle*			class_avg;		// 2D or 3D particle classes
	Bproject() { initialize(); }
} ;

#define _ProjectParamStructs_
#endif

// Function prototypes
Bproject*	project_create(int nmg, int nrec);
int			project_equal_mg_part_files(Bproject* project);
Bfield* 	field_add(Bfield** field, Bstring& field_id);
Bmicrograph*	micrograph_add(Bmicrograph** mg, Bstring& mg_id);
Bframe* 	frame_add(Bframe** frame, int pid);
Bparticle* 	particle_add(Bparticle** part, int pid);
int			mg_part_links(Bmicrograph* mg);
int			rec_part_links(Breconstruction* rec);
Bparticle*	part_find_first(Bproject* project);
Bfilament* 	filament_add(Bfilament** fil, int pid);
Bfilnode* 	filament_node_add(Bfilnode** fnode, int pid);
double 		filament_length(Bfilament* fil);
Breconstruction*	reconstruction_add(Breconstruction** rec, Bstring& rec_id);
int 		project_kill(Bproject* project);
int 		field_kill(Bfield* field);
int 		micrograph_kill(Bmicrograph* mg);
int 		particle_kill(Bparticle* part);
int 		filament_kill(Bfilament* fil);
int 		reconstruction_kill(Breconstruction* rec);
Bparticle*	particle_copy(Bparticle** partlist, Bparticle* part);
Bparticle*	particle_copy(Bparticle* partlist);
Bmicrograph*	micrograph_copy(Bmicrograph* mg);
int			project_update(Bproject* project, Bproject* proj_new, int fom_index);
int			micrograph_update(Bmicrograph* mg, Bmicrograph* nu_mg, int fom_index, int flags);
int			reconstruction_update(Breconstruction* rec, Breconstruction* nu_rec, int fom_index, int flags);
int			particle_update(Bparticle* part, Bparticle* nu_part);
int			particle_update(Bparticle** pnt_part, Bparticle* new_part, int fom_index);
int			project_merge_part_parameters(Bproject* project, Bproject* partproject);
long		project_set_part_links(Bproject* project);
long		micrograph_set_part_links(Bmicrograph* mg);
long		reconstruction_set_part_links(Breconstruction* rec);
long		project_divide(Bproject* project, long n);
Bmicrograph**	project_micrograph_array(Bproject* project, long &nmg);
Bparticle**	project_mg_particle_array(Bproject* project, int part_select, long &npart);
Bparticle**	project_rec_particle_array(Bproject* project, int part_select, long &npart);
Bparticle**	particle_array(Bparticle* partlist, int part_select, long &npart);
int			project_revert_filenames(Bproject* project, Bproject* project_old, int flag);
int			project_set_field_id(Bproject* project, int nseries, Bstring& field_id);
long		project_show_hierarchy(Bproject* project);
long		project_show_class_averages(Bproject* project);
long		project_dump(Bproject* project, Bstring& filename);
int			project_set_micrograph_path(Bproject* project, Bstring& path);
int			project_set_frame_path(Bproject* project, Bstring& path);
int			project_set_powerspectrum_path(Bproject* project, Bstring& path);
int			project_set_particle_path(Bproject* project, Bstring& path);
int			project_set_filament_path(Bproject* project, Bstring& path);
int			project_set_magnification(Bproject* project, double mag);
int			project_set_scan_sampling(Bproject* project, double sampling);
int			project_set_mg_pixel_size(Bproject* project, Vector3<double> pixel_size);
int			project_set_frame_pixel_size(Bproject* project, Vector3<double> pixel_size);
int			project_set_rec_voxel_size(Bproject* project, Vector3<double> pixel_size);
int			project_set_part_pixel_size(Bproject* project, Vector3<double> pixel_size);
int			project_set_tilt(Bproject* project, double tilt_axis, double tilt_angle);
int			project_sort_by_tilt(Bproject* project);
int			project_set_exposure(Bproject* project, double exposure);
int			project_set_dose(Bproject* project, double dose);
int			project_set_dose(Bproject* project, JSvalue& dose_frac);
int			project_set_micrograph_origins(Bproject* project, Vector3<double> origin);
int			project_add_origins_to_coords(Bproject* project);
int			project_flip_origins(Bproject* project, int flip);
long		project_renumber_particles(Bproject* project);
int			project_set_particle_box_size(Bproject* project, Vector3<long> box_size);
int			project_set_particle_box_size(Bproject* project, long box_size);
int			project_set_particle_origins(Bproject* project, Vector3<double> origin);
int			project_set_particle_asu_views(Bproject* project, Bstring& symmetry_asu);
int			project_set_particle_asu_views(Bproject* project, Bsymmetry& sym);
int			project_rotate_particle_views(Bproject* project, View view);
int			project_apply_map_magnifications(Bproject* project, int mag_num, float* mag);
long		project_reset(Bproject* project, Bstring& reset);
View*		views_from_project(Bproject* project, int selection);
Bstring		get_fom_tag(FOMType fom_type);



