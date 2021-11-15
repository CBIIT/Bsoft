/**
@file	radon.cpp
@brief	Program to do radon transforms.
@author	Salvatore Lanzavecchia, Francesca Cantele and Pier Luigi Bellon
         Dip. Chimica Strutturale e Stereochimica Inorganica, Via Venezian 21, 20133 Milano, Italy
@author	Bernard Heymann
         Rm 1515, 50 South Dr., NIH, Bethesda, MD, 20892, USA

@date	Created: 2003 07 04
@date	Modified: 20160728 (BH)
**/

#include "rwimg.h"
#include "img_radon.h"
#include "utilities.h"
#include "options.h"

// Declaration of global variables
extern int verbose;		// Level of output to the screen

// Usage assistence
const char *use[] = {
  " ",
  "Usage: radon [options] input.img output.img",
  "-------------------------------------------",
  "Program to do radon transforms.",
  "Full transform: 2PI*2PI",
  "Quarter transform: PI*PI",
  " ",
  "Types of forward radon transformation:",
  "1 From 3D structure to single axis proj.",
  "2 From single axis proj. to quarter radon transform",
  "3 From 3D structure to full radon transform",
  "4 From 3D structure to quarter radon transform",
  " ",
  "Types of backward radon transformation:",
  "1 From single axis proj. to 3D structure",
  "2 From quarter radon transform to single axis proj.",
  "3 From full radon transform to 3D structure",
  "4 From quarter radon transform to 3D structure",
  " ",
  "POCS specifications:",
  "POCS filter works only on cube quarter radon transform,",
  "that means theta = 2 * side of the cubic 3D structure map.",
  "Side must be equal to a power of 2.",
  "If radon transform has holes, POCS filter needs mask file.",
  " ",
  "Actions:",
  "-forward 4               Forward radon transform (default not)",
  "-pocs 5,5,30,30,0        POCS density modification on radon transform (default not)",
  "                         Parameters (in order): outer cycles, inner cycles, ",
  "                         3D limit, projection limit, support for finite extension",
  "-backward 4              Backward radon transform (default not)",
  " ",
  "Parameters:",
  "-verbose 7               Verbosity of output",
  "-datatype u              Force writing of a new data type (default output floating point)",
  "-theta 64                Number of angles (default 2*ncol)",
  "-kernel 15,4             Kernel size and power (default 11,2)",
  "-nopadding               Turn off padding (default on)",
  " ",
  "Input:",
  "-mask mask.map           Binary mask for POCS filter",
  " ",
  NULL
};

int
main (int argc, char **argv)
{
  // Initialize variables
  DataType newdatatype = Unknown_Type;	// Conversion to new data type

  // Radon transformation types
  int forward(0);		// Forward radon transform
  int backward(0);		// Backward radon transform
  int n_theta(0);		// Number of angular samples per 2*PI
  int nkernel(11);		// Kernel width
  int kernel_power(2);		// Kernel power
  int padd(1);			// Padding??

  // POCS parameters: Both numbers of cycles must be specified to turn this on
  int n_cyc_out(0);		// cycles of swopping between r,phi and r,theta
  int n_cyc_in(0);		// pocs cycles within r,phi and r,theta planes
  double rad_3D(0);		// maximum radius in 3D
  double rad_plane(0);		// maximum radius on single axis proj
  int support(0);		// if (1) the condition of finite extension is imposed also in direct space. (use with care)
  Bstring maskfile;		// File name for mask file
  Bimage *pmask = NULL;		// Binary mask for POCS filter

  int optind;
   Boption *option = get_option_list (use, argc, argv, optind);
   Boption *curropt;
  for (curropt = option; curropt; curropt = curropt->next)
    {
		if ( curropt->tag == "datatype" )
			newdatatype = curropt->datatype();
      if ( curropt->tag == "forward" )
	{
	  if ( ( forward = curropt->value.integer() ) < 1)
	    cerr << "-forward: A type of forward transformation must be specified!" << endl;
	  else if (forward > 4)
	    forward = 4;
	}
      if ( curropt->tag == "pocs" )
	if (sscanf (curropt->value.c_str(), "%d,%d,%lf,%lf,%d", &n_cyc_out, &n_cyc_in, &rad_3D, &rad_plane, &support) < 2)
	  cerr << "-pocs: At least the outer and inner cycles must be specified!" << endl;
      if ( curropt->tag == "backward" )
	{
	  if ( ( backward = curropt->value.integer() ) < 1)
	    cerr << "-backward: A type of backward transformation must be specified!" << endl;
	  else if (backward > 4)
	    backward = 4;
	}
      if ( curropt->tag == "theta" )
	if ( ( n_theta = curropt->value.integer() ) < 1)
	  cerr << "-theta: The number of angles must be specified!" << endl;
      if ( curropt->tag == "kernel" )
	if ( curropt->values(nkernel, kernel_power) < 1)
	  cerr << "-kernel: At least the kernel size must be specified!" << endl;
      if ( curropt->tag == "nopadding" )
	padd = 0;
      if ( curropt->tag == "mask" )
	{
	  maskfile = curropt->filename();
	  if (maskfile.length())
	    {
	      pmask = read_img (maskfile, 1, -1);
	      pmask->change_type(Float);
	    }
	}
    }
  option_kill (option);

  Bimage *p = read_img (argv[optind++], 1, -1);

  if(backward!=3)
    {
    p->calculate_background();
    img_resize_to_next_power2(p, FILL_BACKGROUND, p->background(long(0)));
    }

  if (!p)
    exit (-1);

  if (n_theta < 1)
    n_theta = p->sizeX() * 2;

  if (rad_3D < 1)
    rad_3D = p->sizeX() / 2 - 1;

  if (rad_plane < 1)
    rad_plane = p->sizeX() / 2 - 1;

  Bimage *prad = NULL;

  if (forward)
    {
      prad = img_radon_transform (p, forward, nkernel, kernel_power, padd, n_theta);
      delete p;
      p = prad;
    }

  if (n_cyc_out && n_cyc_in) {
	if ( n_cyc_in > 1 && !pmask ) {
		cout << "Error: To use POCS effectively, a mask file must be given!" << endl;
		exit(-1);
	}
    img_radon_pocs_filter (p, n_cyc_out, n_cyc_in, rad_3D, rad_plane, support, pmask);
  }
	
  if (backward)
    {
      prad = img_radon_inverse_transform (p, backward, nkernel, kernel_power, padd);
      delete p;
      p = prad;
    }

  if (p && optind < argc)
    {
      p->change_type(newdatatype);
      write_img (argv[optind], p, 0);
    }

  delete p;
  delete pmask;

  return 0;
}
