##
# @file		bshow_diffraction.tcl
#
# @brief	Procedures to analyze diffraction patterns
#
# @author	Bernard Heymann
# @date		Created: 20050203
# @date		Modified: 20050204

## @brief Dialog box to analyze diffraction patterns
#

proc Diffraction { } {
	global helv12 theimg
	global gorix goriy gwidth gheight
	global diffplot

	if ![Bimage exists $theimg] {
		tk_dialog .dialog "No image in memory!" \
				{Please read in an appropriate image first.} \
				info 0 {OK}
		return
	}
	
	if { [info exists diffplot] } { unset diffplot }
	if ![winfo exists .wdif] {
		toplevel .wdif
		set gwidth 500
		set gheight 350
		set gorix 100
		set goriy 400
		set xmin 0
		set xmax [expr [Bimage get $theimg gwidth]/2]
		set ymin 0
		set ymax 10
		set xscale 1.0
#		for { set i 0 } { $i < $xmax } { incr i } {
#			set diffplot(x,$i) [expr $i*$xscale + $xmin]
#			set diffplot(y,$i) [lindex $raddiff $i]
#		}

		menu .wdif.menuBar -tearoff 0
		.wdif.menuBar add cascade -menu .wdif.menuBar.dif -label "Diffraction" -underline 0
		menu .wdif.menuBar.dif -tearoff 0
		.wdif configure -menu .wdif.menuBar
		.wdif.menuBar.dif add command -label "Help" \
				-command { showHelp bshow_diffraction.hlp } -underline 0
		.whdif.menuBar.dif add command -label "Reset plot" -underline 0 \
				-command { drawGraph .wdif diffplot "line" "Frequency" "Intensity" } -accelerator "Ctrl-r"
		.wdif.menuBar.dif add command -label "Write data" \
				-command { writePlotData } -underline 0
		.wdif.menuBar.dif add command -label "Close" \
				-command "destroy .wdif" -underline 0
		
		set g .wdif.plot
		canvas $g -relief raised -width [expr $gorix + $gwidth + 50] \
				-height [expr $goriy + 100]
		pack $g -side top -fill x

		frame .wdif.point -bg green
		label .wdif.point.xtag -font $helv12 -text "x" -bg green
		label .wdif.point.x -text "0" -relief sunken -bd 1 -width 8 \
				-font $helv12 -anchor w
		label .wdif.point.ytag -font $helv12 -text "y" -bg green
		label .wdif.point.y -text "0" -relief sunken -bd 1 -width 8 \
				-font $helv12 -anchor w
		pack .wdif.point.xtag .wdif.point.x .wdif.point.ytag .wdif.point.y \
				-side left -pady 10 -padx 10

		frame .wdif.lim -bg cyan
		label .wdif.lim.x_tag -font $helv12 -text "X limits" -bg cyan
		entry .wdif.lim.x_min -width 6 -validate key -vcmd { string is double %P }
		entry .wdif.lim.x_max -width 6 -validate key -vcmd { string is double %P }
		label .wdif.lim.y_tag -font $helv12 -text "Y limits" -bg cyan
		entry .wdif.lim.y_min -width 6 -validate key -vcmd { string is double %P }
		entry .wdif.lim.y_max -width 6 -validate key -vcmd { string is double %P }
		pack .wdif.lim.x_tag .wdif.lim.x_min .wdif.lim.x_max .wdif.lim.y_tag \
				.wdif.lim.y_min .wdif.lim.y_max -side left -pady 10 -padx 10
		pack .wdif.point .wdif.lim -side left		

		updateLimits .wdif $xmin $xmax $ymin $ymax
		plotGraph .wdif $g diffplot "line"
		drawAxes .wdif "Frequency" "Intensity"
	} else {
		wm deiconify .wdif
		raise .wdif
    }
    wm title .wdif "Diffraction"
    wm iconname .wdif "Diffraction"
	$g bind point <Any-Enter> { showPoint .wdif diffplot %x %y }
	$g bind point <Any-Leave> { showPoint .wdif diffplot %x %y }
	bind $g <Button-1> {
		set x1 [.wdif.plot canvasx %x]
		set y1 [.wdif.plot canvasy %y]
	}
	bind $g <B1-Motion> { drawRectangle .wdif.plot $x1 $y1 [.wdif.plot canvasx %x] [.wdif.plot canvasy %y] }
	bind $g <B1-ButtonRelease> { drawSubGraph .wdif .wdif.plot diffplot $x1 $y1 \
		[.wdif.plot canvasx %x] [.wdif.plot canvasy %y] "bar" "Frequency" "Intensity" }
	bind .wdif <Return> { plotGraph .wdif .wdif.plot diffplot "bar" "Frequency" "Intensity" }
	bind .wdif <Control-w> {
		if { [info exists diffplot] } { unset diffplot }
		destroy .wdif
	}
}

