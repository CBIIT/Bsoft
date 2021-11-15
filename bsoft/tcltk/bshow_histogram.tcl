##
# @file		bshow_histogram.tcl
#
# @brief	Procedures to draw a histogram in Bshow.
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20130424

## @brief Draws a histogram of the image data.
#

proc Histogram { } {
	global helv12 theimg
	global gorix goriy gwidth gheight
	global histplot
	set w .whis
	if { [info exists histplot] } { unset histplot }
	if ![winfo exists $w] {
		toplevel $w
		set width 500
		set height 350
		set orix 100
		set oriy 400
		set xmin [Bimage get $theimg min]
		set xmax [Bimage get $theimg max]
		set datatype [Bimage get $theimg datatype]
		set bins [expr int($xmax - $xmin + 1)]
		if { $bins < 10 } { set bins 10 }
		if { $bins > 500 } { set bins 500 }
		if { [string first "float" $datatype] > -1 } { set bins 250 }
		set bar_width [expr $width/(2*$bins)]
		set hist [Bimage histogram $theimg $bins]
		
		set ymin 0.0
		set ymax [lindex $hist 0]
		for { set i 0 } { $i < $bins } { incr i } {
			if { $ymax < [lindex $hist $i] } { set ymax [lindex $hist $i] }
		}
		if { $xmin == $xmax } { set xmax [expr $xmin + 1] }
		if { $ymin == $ymax } { set ymax [expr $ymin + 1] }
#		puts "$datatype $xmin $xmax $ymin $ymax $bins"

		set xscale [expr ($xmax - $xmin)/$bins]
		for { set i 0 } { $i < $bins } { incr i } {
			set histplot(x,$i) [expr $i*$xscale + $xmin]
			set histplot(y,$i) [lindex $hist $i]
		}

		menu $w.menuBar -tearoff 0
		$w.menuBar add cascade -menu $w.menuBar.his -label "Histogram" -underline 0
		menu $w.menuBar.his -tearoff 0
		$w configure -menu $w.menuBar
		$w.menuBar.his add command -label "Help" \
				-command { showHelp bshow_graph.hlp } -underline 0
		$w.menuBar.his add command -label "Reset plot" -underline 0 \
				-command { drawGraph $w $w.plot histplot "bar" } -accelerator "Ctrl-r"
		bind . <Control-r> { drawGraph $w $w.plot histplot "bar" }
		
		set h $w.plot
		canvas $h -relief raised -width [expr $orix + $width + 50] \
				-height [expr $oriy + 100]
		pack $h -side top -fill x

		frame $w.point -bg green
		label $w.point.xtag -font $helv12 -text "x" -bg green
		label $w.point.x -text "0" -relief sunken -bd 1 -width 8 \
				-font $helv12 -anchor w
		label $w.point.ytag -font $helv12 -text "y" -bg green
		label $w.point.y -text "0" -relief sunken -bd 1 -width 8 \
				-font $helv12 -anchor w
		pack $w.point.xtag $w.point.x $w.point.ytag $w.point.y \
				-side left -pady 10 -padx 10

		frame $w.lim -bg cyan
		label $w.lim.x_tag -font $helv12 -text "X limits" -bg cyan
		entry $w.lim.x_min -width 6 -validate key -vcmd { string is double %P }
		entry $w.lim.x_max -width 6 -validate key -vcmd { string is double %P }
		label $w.lim.y_tag -font $helv12 -text "Y limits" -bg cyan
		entry $w.lim.y_min -width 6 -validate key -vcmd { string is double %P }
		entry $w.lim.y_max -width 6 -validate key -vcmd { string is double %P }
		pack $w.lim.x_tag $w.lim.x_min $w.lim.x_max $w.lim.y_tag \
				$w.lim.y_min $w.lim.y_max -side left -pady 10 -padx 10
		pack $w.point $w.lim -side left		

		updateLimits $w $xmin $xmax $ymin $ymax
		plotGraph $w $h histplot "bar" "Bin" "Count"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Histogram"
    wm iconname $w "Histogram"
	$h bind point <Any-Enter> { showPoint $w histplot %x %y }
	$h bind point <Any-Leave> { showPoint $w histplot %x %y }
	bind $h <Button-1> {
		set x1 [$w.plot canvasx %x]
		set y1 [$w.plot canvasy %y]
	}
	bind $h <B1-Motion> { drawRectangle $w.plot $x1 $y1 [$w.plot canvasx %x] [$w.plot canvasy %y] }
	bind $h <B1-ButtonRelease> { drawSubGraph $w $w.plot histplot $x1 $y1 \
		[$w.plot canvasx %x] [$w.plot canvasy %y] "bar" }
	bind $w <Return> { plotGraph $w $w.plot histplot "bar" "Bin" "Count" }
	bind $w <Control-w> {
		if { [info exists histplot] } { unset histplot }
		destroy $w
	}
}

