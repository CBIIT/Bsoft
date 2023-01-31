/**
@file	ps_plot.cpp
@brief	Postscript output functions.
@author Bernard Heymann
@date	Created: 20010515
@date	Modified: 20200312
**/
 
#include "ps_plot.h" 
#include "utilities.h" 

#include <iostream>
#include <time.h>

// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern string	command;		// Command line


/**
@brief 	Opens and initializes a postscript file.
@param 	&filename	output postscript file name.
@param 	&title		title.
@param 	npages		number of pages.
@param 	width		page width.
@param 	height		page height.
@return ofstream*	postscript file descriptor.

	Opens a stream and writes the postscript header line, the creator and
	creation date lines, and a document data specifier.

**/
ofstream*	ps_open_and_init(Bstring filename, Bstring title, int npages, 
					int width, int height)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_open_and_init: title=" << title << " npages=" 
			<< npages << " width=" << width << " height=" << height << endl;

	ofstream*	fps = new ofstream(filename.c_str());
	
	time_t		t = time(NULL);

	*fps << "%!PS-Adobe-2.0" << endl;
	*fps << "%%Creator: Bsoft" << endl;
	*fps << "%%CreationDate: " << asctime(localtime(&t));
	*fps << "%%Title: " << title << endl;
	*fps << "%%DocumentData: Clean7Bit" << endl;
	*fps << "%%LanguageLevel: 2" << endl;
	*fps << "%%Origin: 0 0" << endl;
	*fps << "%%Pages: " << npages << endl;
	*fps << "%%BoundingBox: 0 0 " << width << " " << height << endl;
	*fps << "%%EndComments" << endl;
	*fps << "%%BeginSetup" << endl;
	*fps << "%%EndSetup" << endl;
	*fps << "%%EndProlog" << endl;
	
	return fps;
}

/**
@brief 	Opens and initializes a postscript file.
@param 	&filename	output postscript file name.
@param 	*plot		plot information.
@return ofstream* 	postscript file stream.

	Opens a file and writes the postscript header line, the creator and 
	creation date lines, and a document data specifier.

**/
ofstream*	ps_open_and_init(Bstring filename, Bplot* plot)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_open_and_init: title=" << plot->title() << " npages=" 
			<< plot->pages() << " width=" << plot->page_width() << " height=" << plot->page_height() << endl;

	ofstream*	fps = new ofstream(filename.c_str());
	
	time_t		t = time(NULL);

	*fps << "%!PS-Adobe-2.0" << endl;
	*fps << "%%Creator: Bsoft" << endl;
	*fps << "%%CreationDate: " << asctime(localtime(&t));
	*fps << "%%Title: " << command << endl;
	*fps << "%%DocumentData: Clean7Bit" << endl;
	*fps << "%%LanguageLevel: 2" << endl;
	*fps << "%%Origin: 0 0" << endl;
	*fps << "%%Pages: " << plot->pages() << endl;
	*fps << "%%BoundingBox: 0 0 " << plot->page_width() << " " << plot->page_height() << endl;
	*fps << "%%EndComments" << endl;
	*fps << "%%BeginSetup" << endl;
	*fps << "%%EndSetup" << endl;
	*fps << "%%EndProlog" << endl;
	
	return fps;
}

/**
@brief 	Closes a postscript file.
@param 	*fps		file handle.
@return int			0.

	Closes a file and writes the postscript trailer line.

**/
int			ps_close(ofstream* fps)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_close" << endl;
		
	*fps << "%%Trailer" << endl;
	*fps << "%%EOF" << endl << endl;

	fps->close();
	
	delete fps;
	
	return 0;
}

/**
@brief 	Generates postscript graphs.
@param 	&filename	output postscript file name.
@param 	*plot		plot information.
@return int 		0, <0 on error.

	Any number of functions of a single independant variable (x) can be plot.

**/
int			ps_plot(Bstring filename, Bplot* plot)
{
	if ( verbose )
		cout << "Plotting a graph to " << filename << endl << endl;
	
	ofstream*	fps = ps_open_and_init(filename, plot);
	
	for ( int i=1; i<=plot->pages(); i++ ) ps_graph(fps, plot, i);
	
	ps_close(fps);

	return 0;
}

/**
@brief 	Reads a postscript file with a data table.
@param 	&filename	input postscript file name.
@return Bplot* 		new plot.

**/
Bplot*		ps_read(Bstring& filename)
{
	if ( filename.length() < 1 ) return NULL;
	
	if ( verbose & VERB_PROCESS )
		cout << "# Reading Postscript file:       " << filename << endl;

	Bplot*				plot = NULL;
	
	ifstream			fps(filename.c_str());
	if ( fps.fail() ) return plot;

	char				aline[1024];
	long				i(0), j;
	int					ncol(0), nrow(0);
	Bstring				sline;

	// Title
	while ( fps.getline(aline, 1024) && strncmp(aline, "%%Title:", 8) ) ;
	
	Bstring				title(aline);
	title = title.substr(8,title.length());

	// Start of data
	while ( fps.getline(aline, 1024) && strncmp(aline, "/Data [", 7) ) ;

	// Lables line
	fps.getline(aline, 1024);
	sline = aline;
	Bstring*			labels = sline.split();
	Bstring*			label;
	
	// Determine number of rows and columns
	while ( fps.getline(aline, 1024) && aline[0] != ']' ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG ps_read: aline = \"" << aline << "\"" << endl;
		sline = aline;
		sline = sline.substr(0, sline.length() - 1);
		vector<double>	temp = sline.split_into_doubles(0);
		if ( ncol < temp.size() ) ncol = temp.size();
		nrow++;
	}
	
	if ( verbose ) {
		cout << "Plot file:                      " << filename << endl;
		cout << "Columns:                        " << ncol << endl;
		cout << "Rows:                           " << nrow << endl;
		cout << "Labels:                        ";
	}

	plot = new Bplot(1, nrow, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0, label = labels; i<ncol; i++, label = label->next ) {
		plot->page(0).column(i).number(i);
		plot->page(0).column(i).label(*label);
		if ( verbose )
			cout << " " << *label;
	}

	if ( verbose )
		cout << endl << endl;

	string_kill(labels);
	
	fps.seekg(0);
	
	// Start of data
	while ( fps.getline(aline, 1024) && strncmp(aline, "/Data [", 7) ) ;

	// Lables line
	fps.getline(aline, 1024);

	// Read the data
	i = 0;
	while ( fps.getline(aline, 1024) && aline[0] != ']' ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG ps_read: aline = \"" << aline << "\"" << endl;
		sline = aline;
		vector<double>	temp = sline.split_into_doubles(0);
		for ( j=0; j<ncol; j++ ) (*plot)[j*nrow+i] = temp[j];
		i++;
	}

	return plot;
}

/**
@brief 	Generates postscript graphs.
@param 	*fps		output postscript file descriptor.
@param 	*plot		plot information.
@param 	page_number	page number (starts at 1).
@return int 		0, <0 on error.

	Any number of functions of a single independant variable (x) can be plot.

**/
int 		ps_graph(ofstream* fps, Bplot* plot, int page_number)
{
	if ( page_number < 1 ) {
		cerr << "Error in ps_graph: A page number less than 1 is not valid!" << endl;
		bexit(-1);
	}
	
	plot->limits();
	
	Bpage		page = plot->page(page_number-1);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_graph: starting page " << page.number() << endl;
	
	if ( page.columns() < 2 ) {
		cerr << "Error in ps_graph: At least two columns must be specified for page " << page.number() << endl;
		return -1;
	}

	page.limits();
	
	long 		i, j, left = plot->left(), bottom = plot->bottom(), width = plot->width(), height = plot->height();
	long		c, ncol(0), nrow(0), x_inverse, y_inverse, y2_inverse;
	long		max_row(32767);
	Bstring		x_label, y_label, y2_label, *txt;
	double		x_min, x_max, x_inc;
	double		y_min, y_max, y_inc;
	double		y2_min, y2_max, y2_inc;
	double		v_min, v_max;
	double		y_red(0), y_green(0), y_blue(0), y2_red(0), y2_green(0), y2_blue(0);
	Bcolumn		col;
	
	nrow = plot->rows();
	ncol = page.columns();
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG ps_graph: page index = " << page_number-1 << " columns = " << ncol << endl;
		cout << "DEBUG ps_graph: indices =";
		for ( i=0; i<ncol; i++ ) cout << " " << page.column(i).number();
		cout << endl;
		cout << "DEBUG ps_graph: axes =";
		for ( i=0; i<ncol; i++ ) cout << " " << page.column(i).axis();
		cout << endl;
	}

	x_min = page.axis(1).min();
	x_max = page.axis(1).max();
	x_inc = page.axis(1).inc();
	x_inverse = page.axis(1).flags() & 1;
	y_min = page.axis(3).min();
	y_max = page.axis(3).max();
	y_inc = page.axis(3).inc();
	y_inverse = page.axis(3).flags() & 1;
	y2_min = page.axis(4).min();
	y2_max = page.axis(4).max();
	y2_inc = page.axis(4).inc();
	y2_inverse = page.axis(4).flags() & 1;
	v_min = page.axis(5).min();
	v_max = page.axis(5).max();
	
//	cout << "v_min = " << v_min << " v_max = " << v_max << endl;

	x_label = page.axis(1).label();
	y_label = page.axis(3).label();
	y2_label = page.axis(4).label();

	for ( i=0; i<ncol; i++ ) {
		col = page.column(i);
		if ( col.axis() == 1 ) {
			if ( x_label.length() < 1 ) x_label = col.label();
		} else if ( col.axis() == 3 ) {
			if ( y_label.length() < 1 ) y_label = col.label();
			y_red = col.red();
			y_green = col.green();
			y_blue = col.blue();
		} else if ( col.axis() == 4 ) {
			if ( y2_label.length() < 1 ) y2_label = col.label();
			y2_red = col.red();
			y2_green = col.green();
			y2_blue = col.blue();
		}
	}
	
	long		hw = (int) (width*1.0/nrow);	// Histogram bar width
	if ( hw < 1 ) hw = 1;
	long		hhw = hw/2;
	
	long		x_digits = (long) (3 - log(x_max - x_min)/log(10));
	long		y_digits = (long) (3 - log(y_max - y_min)/log(10));
	long		y2_digits = (long) (3 - log(y2_max - y2_min)/log(10));
	
	*fps << "%%Page: " << page.number() << " " << page.number() << endl;
	*fps << "/Helvetica findfont 14 scalefont setfont" << endl;
	*fps << "40 740 moveto (" << page.title() << ") show" << endl;
	*fps << "/Left " << left << " def" << endl;
	*fps << "/Bottom " << bottom << " def" << endl;
	*fps << "/Height " << height << " def" << endl;
	*fps << "/Width " << width << " def" << endl;
	*fps << setprecision(4) << "/x_Min " << x_min << " def" << endl;
	*fps << "/y_Min " << y_min << " def" << endl;
	*fps << "/y2_Min " << y2_min << " def" << endl;
	*fps << "/v_Min " << v_min << " def" << endl;
	*fps << "/x_Max " << x_max << " def" << endl;
	*fps << "/y_Max " << y_max << " def" << endl;
	*fps << "/y2_Max " << y2_max << " def" << endl;
	*fps << "/v_Max " << v_max << " def" << endl;
	*fps << "/x_Scale Width x_Max x_Min sub div def" << endl;
	*fps << "/y_Scale Height y_Max y_Min sub div def" << endl;
	*fps << "/y2_Scale Height y2_Max y2_Min sub div def" << endl;
//	*fps << "/v_Scale Height v_Max v_Min sub div def" << endl;
	*fps << "/v_Scale 1 v_Max v_Min sub div def" << endl;
	*fps << "/Ncol " << ncol << " def" << endl;
	*fps << "/Frame { newpath 0 0 moveto 0 Height lineto Width Height lineto Width 0 lineto closepath } def" << endl;
	
	*fps << "/Data [" << endl;
	*fps << "%";
	for ( j=0; j<ncol; j++ ) {
		col = page.column(j);
		*fps << col.label() << tab;
	}
	*fps << endl;
	*fps << setprecision(8);
	for ( i=0; i<nrow && i<max_row; i++ ) {
		for ( j=0; j<ncol; j++ ) {
			c = page.column(j).number();
//			cout << c << ":" << i << " " << (*plot)[c*nrow+i] << endl;
			*fps << (*plot)[c*nrow+i] << tab;
		}
		*fps << endl;
	}
	*fps << "] def" << endl;

	*fps << "/Helvetica findfont 12 scalefont setfont" << endl;
	for ( i=140, txt = page.text(); txt; txt = txt->next, i-=18 )
		*fps << "40 " << i << " moveto (" << *txt << ") show" << endl;
	
	*fps << "/Helvetica findfont 14 scalefont setfont" << endl;
	*fps << "/DeviceRGB setcolorspace" << endl;
	*fps << "0 0 0 setcolor" << endl;

	if ( x_label.length() )
		*fps << left+(width-6*x_label.length())/2 << " " << bottom-40 <<
			" moveto (" << x_label << ") show" << endl;

	if ( y_label.length() ) {
		*fps << y_red << " " << y_green << " " << y_blue << " setcolor" << endl;
		*fps << left-40 << " " << bottom+(height-6*y_label.length())/2 <<
			" moveto gsave 90 rotate (" << y_label << ") show grestore" << endl;
	}
	
	if ( y2_label.length() ) {
		*fps << y2_red << " " << y2_green << " " << y2_blue << " setcolor" << endl;
		*fps << left+width+35 << " " << bottom+(height-6*y2_label.length())/2 <<
			" moveto gsave -90 rotate (" << y2_label << ") show grestore" << endl;
	}
	
	*fps << "0 0 0 setcolor" << endl;
	
	ps_scale(fps, left, bottom, left+width, bottom, x_min, x_max, x_inc,
			10, M_PI_2, x_digits, 10, x_inverse);
	
	ps_scale(fps, left, bottom, left, bottom+height, y_min, y_max, y_inc,
			10, 0, y_digits, 10, y_inverse);

	if ( y2_label.length() )
		ps_scale(fps, left+width, bottom, left+width, bottom+height, y2_min, y2_max, y2_inc,
			10, M_PI, y2_digits, 10, y2_inverse);

	*fps << "gsave" << endl;
	*fps << "	/DeviceRGB setcolorspace" << endl;
	*fps << "	Left Bottom translate Frame stroke Frame clip" << endl;
	*fps << "	/y 0 y_Min sub y_Scale mul def" << endl;
	*fps << "	0 y moveto Width y lineto stroke" << endl;
	*fps << "	/x 0 x_Min sub x_Scale mul def" << endl;
	*fps << "	x 0 moveto x Height lineto stroke" << endl;
	
	for ( i=0; i<ncol; i++ ) {
		col = page.column(i);
		if ( col.axis() == 3 || col.axis() == 4 ) {
			*fps << "	" << col.red() << " " << col.green() << " " << col.blue() << " setcolor" << endl;
			*fps << "	/x Data 0 get x_Min sub x_Scale mul def" << endl;
			*fps << "	/y Data " << i << " get y_Min sub y_Scale mul def" << endl;
			*fps << "	x y moveto" << endl;
			*fps << "	Ncol Ncol Data length Ncol sub {" << endl;
			*fps << "		/x_Index exch def" << endl;
			*fps << "		/y_Index x_Index " << i << " add def" << endl;
			if ( col.type() == 3 ) *fps << "		/py_Index y_Index 1 sub def" << endl;
			if ( col.type() == 4 ) *fps << "		/v_Index y_Index 1 add def" << endl;
			*fps << "		/x Data x_Index get x_Min sub x_Scale mul def" << endl;
			if ( col.axis() == 3 ) {
				if ( col.type() == 3 ) *fps << "		/y Data y_Index get y_Scale mul def" << endl;
				else *fps << "		/y Data y_Index get y_Min sub y_Scale mul def" << endl;
			} else *fps << "		/y Data y_Index get y2_Min sub y2_Scale mul def" << endl;
			if ( col.type() == 3 ) *fps << "		/py Data py_Index get y_Min sub y_Scale mul def" << endl;
			*fps << "		/w " << hw << " def" << endl;
			*fps << "		/xm x " << hhw << " sub def" << endl;
//			if ( col.type() == 4 ) *fps << "		/v Data v_Index get v_Min sub v_Scale mul def" << endl;
			if ( col.type() == 4 ) {
				*fps << "		/v v_Max Data v_Index get sub v_Scale mul def" << endl;
				*fps << "		v setgray" << endl;
			}
			switch ( col.type() ) {
				case 1: *fps << "		xm 0 w y rectstroke" << endl; break;	// bar
				case 2: *fps << "		x y lineto" << endl; break;				// line
				case 3:															// std line
					*fps << "		x py y sub moveto 0 y 2 mul rlineto" << endl;
					*fps << "		x 3 sub py y sub moveto 6 0 rlineto" << endl;
					*fps << "		x 3 sub py y add moveto 6 0 rlineto" << endl;
					break;
				case 4: *fps << "		x y 5 0 360 arc fill" << endl; break;	// variable circle
				default: *fps << "		x y " << col.element_size() << " 0 360 arc fill" << endl; // circle
			}
			*fps << "	} for stroke" << endl;
		}
	}
	*fps << "grestore" << endl;
	*fps << "showpage" << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ps_graph: finished page " << page.number() << endl;
		
	return nrow;
}

/**
@brief 	Generates a scale for a postscript graph.
@param 	*fps		postscript file descriptor.
@param 	x1			starting x coordinate.
@param 	y1			starting y coordinate.
@param 	x2			ending x coordinate.
@param 	y2			ending y coordinate.
@param 	min			minimum value.
@param 	max			maximum value.
@param 	increment 	value increment.
@param 	tick_length	length of ticks.
@param 	tick_angle	angle of ticks relative to horizontal axis.
@param 	digits		number of decimal digits for labels.
@param 	fontsize	font size for labels.
@param 	inverse		flag for inverse values.
@return int 		0.
**/
int 		ps_scale(ofstream* fps, double x1, double y1, double x2, double y2,
					double min, double max, double increment,
					double tick_length, double tick_angle,
					int digits, int fontsize, int inverse)
{
	ios_base::fmtflags	fl = fps->flags();
	
	long 		i, j;
	double		value, showval;
	double		x_tick = tick_length*cos(tick_angle);
	double		y_tick = tick_length*sin(tick_angle);
	double		x_step = (x2-x1)*increment/(max - min);
	double		y_step = (y2-y1)*increment/(max - min);
	
	if ( digits < 0 ) digits = 0;
	if ( fontsize < 2 ) fontsize = 10;
	
	long 		integer_size = (int) (log(fabs(max))/log(10.0) + 1);
	if ( fabs(max) < fabs(min) ) integer_size = (int) (log(fabs(min))/log(10.0) + 1);
	if ( inverse ) {
		if ( min <= 0 ) integer_size = (int) (log(fabs(10/max))/log(10.0) + 1);
		else integer_size = (int) (log(fabs(1/min))/log(10.0) + 1);
	}
	if ( integer_size < 1 ) integer_size = 1;
	integer_size++;		// Add one for the sign
	
	long 		number_size = integer_size + 1 + digits;
	double		x_offset = x1 - number_size*fontsize/5;
	double		y_offset = y1 - fontsize/3;
	if ( x2 - x1 > y2 - y1 )
		y_offset -= fontsize;
	else
		x_offset -= number_size*fontsize/3;
	
	if ( tick_angle > 0.75*M_PI && tick_angle < 1.25*M_PI )
		x_offset = x1 + fontsize/5;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG ps_scale: min=" << min << " max=" << max << " increment=" << increment << endl;
		cout << "DEBUG ps_scale: number_size=" << number_size << " digits=" << digits << endl;
	}
	
	// Ticks and labels
	*fps << "/Helvetica findfont " << fontsize << " scalefont setfont " << endl;
	for ( i=0, value=min; value<max+0.01*increment; value+=increment, i++ ) {
		*fps << fixed << setprecision(2) <<
				x1+i*x_step << " " << y1+i*y_step << " moveto " <<
				x1+i*x_step+x_tick << " " << y1+i*y_step+y_tick << " lineto stroke" << endl;
		showval = value;
		if ( inverse && value ) { if ( value ) showval = 1/value; else showval = 999; }
//		cout << showval << endl;
		*fps << x_offset + i*x_step << " " << y_offset + i*y_step << " moveto (" 
			<< fixed << setprecision(digits) << right << setw(number_size) << showval << ") (";
		for ( j=0; j<integer_size; j++ ) *fps << "x";
		*fps << ".";
		for ( j=0; j<digits; j++ ) *fps << "x";
		*fps << ") cvs show" << endl;
	}
	
//	fps->unsetf(fps->flags());
	fps->flags(fl);

	return 0;
}

/**
@brief 	Generates 1-3 line plot.
@param 	nrow		number of points.
@param 	*c0			independent variable.
@param 	*c1			curve 1.
@param 	*c2			curve 2.
@param 	*c3			curve 3.
@return Bplot* 		plot structure.
**/
Bplot*		plot_curve(long nrow, double* c0, double* c1, double* c2, double* c3)
{
	int			ncol(4);
	Bstring		title("Curve");
	Bplot*		plot = new Bplot(1, nrow, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( int i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("x");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(2);
	plot->page(0).column(1).label("y");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(2).type(2);
	plot->page(0).column(2).label("f");
	plot->page(0).column(2).axis(3);
	plot->page(0).column(3).type(2);
	plot->page(0).column(3).label("c");
	plot->page(0).column(3).axis(3);
	plot->page(0).column(1).color(1,0,0);
	plot->page(0).column(2).color(0,0,1);
	plot->page(0).column(3).color(0,1,0);
//	plot->page(0).axis(1).min(0);
//	plot->page(0).axis(1).max(nrow);
//	plot->page(0).axis(3).min(prad->minimum());
//	plot->page(0).axis(3).max(prad->maximum());

	for ( int i=0, j=nrow, k=j+nrow, l=k+nrow; i<nrow; i++, j++, k++, l++ ) {
		if ( c0 ) (*plot)[i] = c0[i];
		if ( c1 ) (*plot)[j] = c1[i];
		if ( c2 ) (*plot)[k] = c2[i];
		if ( c3 ) (*plot)[l] = c3[i];
	}
	
	return plot;
}

/**
@brief 	Generates an arrow and line.
**/
int			ps_define_arrowline(ofstream* fps)
{
	*fps << "/ArrowHeadSize " << 30 << " def" << endl;
	*fps << "/ahead {" << endl;
	*fps << "    1 index 4 index sub" << endl;
	*fps << "    1 index 4 index sub" << endl;
	*fps << "    exch atan" << endl;
	*fps << "    ArrowHeadSize -.8 mul" << endl;
	*fps << "    dup" << endl;
	*fps << "    2 index cos mul 4 index add" << endl;
	*fps << "    exch" << endl;
	*fps << "    2 index sin mul 3 index add" << endl;
	*fps << "    5 2 roll" << endl;
	*fps << "    gsave" << endl;
	*fps << "        3 1 roll" << endl;
	*fps << "        translate" << endl;
	*fps << "        rotate" << endl;
	*fps << "        newpath" << endl;
	*fps << "        0 0 moveto" << endl;
	*fps << "        ArrowHeadSize dup neg exch .25 mul" << endl;
	*fps << "        2 copy lineto" << endl;
	*fps << "        ArrowHeadSize -.8 mul 0" << endl;
	*fps << "        2 copy" << endl;
	*fps << "        6 4 roll" << endl;
	*fps << "        neg curveto" << endl;
	*fps << "        closepath fill" << endl;
	*fps << "    grestore" << endl;
	*fps << "} bind def" << endl;
	*fps << "/arrowline {" << endl;
	*fps << "    ahead" << endl;
	*fps << "    moveto" << endl;
	*fps << "    lineto" << endl;
	*fps << "    stroke" << endl;
	*fps << "} bind def" << endl;
	
	return 0;
}
