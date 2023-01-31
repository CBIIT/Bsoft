		if ( curropt->tag == "Volt" ) {
			if ( ( v = curropt->real_units() ) < 1 )
				cerr << "-Volt: A voltage must be specified!" << endl;
			else
				jsctf["Volt"] = v;
		}
		if ( curropt->tag == "Amplitude" ) {
			v = curropt->value.real();
			if ( v < -1 || v > 1 )
				cerr << "-Amplitude: A fraction must be specified!" << endl;
			else
				jsctf["Amplitude_phase"] = asin(v);
		}
		if ( curropt->tag == "Cs" ) {
			v = curropt->real_units();
			jsctf["Cs"] = v;
		}
		if ( curropt->tag == "Cc" ) {
			if ( ( v = curropt->real_units() ) < 0.1 )
				cerr << "-Cc: A Cc value must be specified!" << endl;
			else
				jsctf["Cc"] = v;
		}
		if ( curropt->tag == "alpha" ) {
			if ( ( v = curropt->value.real() ) < 1e-10 )
				cerr << "-alpha: A value must be specified!" << endl;
			else {
				if ( v > 0.01 ) v /= 1000;	// Assume milliradians
				jsctf["Beam_convergence"] = v;
			}
		}
		if ( curropt->tag == "energyspread" ) {
			if ( ( v = curropt->value.real() ) < 1e-10 )
				cerr << "-energyspread: A value must be specified!" << endl;
			else
				jsctf["Energy_spread"] = v;
		}
		if ( curropt->tag == "aperture" ) {
			if ( ( v = curropt->real_units() ) < 1 )
				cerr << "-aperture: An aperture diameter must be specified!" << endl;
			else
				jsctf["Objective_aperture"] = v;
		}
		if ( curropt->tag == "focallength" ) {
			if ( ( v = curropt->real_units() ) < 1 )
				cerr << "-focallength: A focal length must be specified!" << endl;
			else
				jsctf["Focal_length"] = v;
		}
		if ( curropt->tag == "slitwidth" ) {
			if ( ( v = curropt->value.real() ) < 1 )
				cerr << "-slitwidth: A slit width must be specified!" << endl;
			else
				jsctf["Slit_width"] = v;
		}
