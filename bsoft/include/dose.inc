		if ( curropt->tag == "dose" ) {
			dose_frac["dose_per_frame"] = curropt->value.real();
    	    if ( dose_frac["dose_per_frame"].real() < 0 )
                cerr << "-dose: The dose per frame must be specified!" << endl;
		}
		if ( curropt->tag == "exposure" ) {
			dose_frac["exposure"] = curropt->value.real();
    	    if ( dose_frac["exposure"].real() < 0 )
				cerr << "-exposure: An exposure time in seconds must be specified." << endl;
		}
		if ( curropt->tag == "rate" ) {
			dose_frac["frame_rate"] = curropt->value.real();
    	    if ( dose_frac["frame_rate"].real() < 0 )
				cerr << "-rate: The number of frames per second must be specified." << endl;
		}
