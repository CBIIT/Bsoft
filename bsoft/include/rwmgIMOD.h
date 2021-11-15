/**
@file	rwmgIMOD.h
@brief	Converts between IMOD files and a micrograph parameter file
@author	Bernard Heymann
@date	Created: 20070501
@date	Modified: 20170120
**/

#include "mg_processing.h"

#define	IMODSIZE	232
#define	OBJTSIZE	176
#define	CONTSIZE	16

struct IMOD {
	char name[128];		// Name of model file.
   	int   xmax;         //  Maximum values for x,y and z. These are usually
   	int   ymax;         //  set to the image dimensions.
   	int   zmax;

   	int  objsize;      //  Number of objects in model.

   	uint flags;        //  Model flags  (IMODF_FLAG...)
                      //  Bit 12 on : multiple clip planes are possible
                      //  Bit 13 on : mat1 and mat3 are stored as bytes
                      //  Bit 14 on : otrans has image origin values
                      //  Bit 15 on : current tilt angles are stored correctly
                      //  Bit 16 on : model last viewed on Y/Z flipped image
   int   drawmode;     //  1 to draw model, -1 not to
   int   mousemode;    //  Mouse editing mode. (0=movie,1=model)

   int   blacklevel;   //  Contrast adjustment. (0-256) Default = 0.
   int   whitelevel;   //  Contrast adjustment. (0-256) Default = 255.

   float xoffset;      //  Offsets used for display & conversion,
   float yoffset;      //  should be set to 0.
   float zoffset;      //  (unused, superceded by MINX information)

   float xscale;       //  Scaling for pixel size. xscale and yscale should
   float yscale;       //  be 1.0 and zscale should be adjusted to account
   float zscale;       //  for the image thickness.

   int   object;       //  Current object, contour and point, used for
   int   contour;      //  model editing.
   int   point;

   int   res;          //  Minimum number of pixels between points when adding
                      //  points by dragging mouse with middle button down.
   int   thresh;       //  Threshold level for auto contour generator.

   float pixsize;      //  Size of one pixel, in the units given next
   int   units;        //  0 = pixels, 3 = km, 1 = m, -2 = cm, -3 = mm,
                      //  -6 = microns, -9 = nm, -10 = Angstroms, -12 = pm

   int   csum;         //  Checksum storage. Used for autosave only.

   float alpha;        //  Angles used for orientation of the model with
   float beta;         //  image. Should be set to 0.0
   float gamma;        //  (unused, superceded by MINX information)
};

struct OBJT {
   char  name[64];      //  Name of Object. Null terminated string.
   uint  extra[16];     //       Extra data for future use, indexed by IOBJ_EX_ defines
                      // 0 (PNT_LIMIT): Limit on number of points per contour

   int   contsize;   // Number of Contours in object.
   uint  flags;      // bit flags for object (IMOD_OBJFLAG...).
                        // Bit 1 on : Off - turns off display
                        // Bit 2 on : Draw using depth cue
                        // Bit 3 on : Open - contours are not closed
                        // Bit 4 on : Wild - contours not constrained in Z
                        // Bit 5 on : Inside-out  - light inside surface
                        // Bit 6 on : Use fill color for spheres
                        // Bit 7 on : Draw spheres on central section only
                        // Bit 8 on : Fill - draw filled polygons
                        // Bit 9 on : Scattered - contours have scattered points
                        // Bit 10 on : Mesh - draw mesh data in 3D
                        // Bit 11 on : Lines - do not draw contours in 3D
                        // Bit 12 on : Use stored values to modify drawing
                        // Bit 13 on : Keep contours planar in open object
                        // Bit 14 on : Fill color - use instead of regular color
                        // Bit 15 on : Anti-aliasing on
                        // Bit 16 on : Normals have magnitudes
                        // Bit 17 on : Draw stored values with false color
                        // Bit 18 on : Time - contours have time indexes
                        // Bit 19 on : Light both sides of surface
                        // Bit 20 on : Draw current contour thicker in model view
                        // Bit 21 on : Draw extra object in model view
                        // Bit 22 on : Allow extra object editing in model view
                        // Bit 23 on : Do not draw spheres in model view
                        // Bit 24 on : Draw extra object only in model view
                        // Bit 31 reserved for temporary use
  
                        // In 3dmodv, Draw Data Type - Contour = Mesh flag off,
                        //                            Mesh    = Mesh flag on
                        // Drawing Style - Points = Lines flag on, Fill flag off
                        //                 Lines  = Lines flag off, Fill flag off
                        //                 Fill   = Lines flag on, Fill flag on
                        //          Fill Outline  = Lines flag off, Fill flag on

    int   axis;       // Z = 0, X = 1, Y = 2. (unused)
    int   drawmode;   // Tells type of scattered points to draw (unused)

    float red;        // Color values, range is (0.0 - 1.0)
    float green;
    float blue;

    int   pdrawsize;    // Default radius in pixels of scattered points in 3D.

    unsigned char symbol;       // Point Draw symbol in 2D, default is 1.
                       //    0 = circle, 1 = none, 2 = square, 3 = triangle,
	                     // 4 = star.
    unsigned char symsize;      // Size of 2D symbol; radius in pixels.
    unsigned char linewidth2;   // Linewidth in 2-D view.
    unsigned char linewidth;    // Linewidth in 3-D view.
    unsigned char linesty;      // Line draw style, 0 = solid; 1 = dashed (unused).
    unsigned char symflags;     // Bit 0 on : Fill the symbol.
                       // Bit 1 on : Draw beginning and end symbols.
    unsigned char sympad;       // Set to 0, for future use.
    unsigned char trans;        // Transparency, range is (0 to 100), maps to
                        //   (1.0 to 0.0) in call to glColor4f or glMaterialfv

    int   meshsize;     // Number of meshes in object.
    int   surfsize;     // Max surfaces in object.
};

struct CONT {
      int          psize;   // Number of points in contour.
      uint         flags;   // Bit 3 on : open, do not connect endpoints
                           // Bit 4 on : wild, not in one Z plane
                           // Bit 5 on : draw stippled lines
                           // Bit 6 on : draw if  mouse in window, ignore Z
                           // Bit 7 on : draw on all planes, regardless of Z
                           // Bit 8 on : draw only in model mode
                           // Bit 17 on : contour has pairs of scan lines
                           // Bits 18, 19, 20, 31: reserved for temporary use
      int          time;    // The time index used by this contour.
      int          surf;    // The surface index used by this contour.
//	float*    pt;     // Array of point triplets.
};

// Function prototypes
int			read_imod_tlt(Bstring& filename, Bproject* project);
int			read_project_imod(Bstring* file_list, Bproject* project, int flag);
int			write_project_imod(Bstring& imodfile, Bproject* project);

