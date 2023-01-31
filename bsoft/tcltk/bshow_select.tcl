##
# @file		bshow_select.tcl
#
# @brief	Procedures to draw selections in Bshow.
#
# @author	Bernard Heymann
# @date		Created: 20060811
# @date		Modified: 20130725

set seltype "rectangle"
set calc_stats 0

## @brief Handles selections from the image data.
#

proc Selection { } {
	global helv12 theimg
	global selectiontype
	
	set c [getImageCanvas $theimg]
	
	if { [info exists selpoly] } { unset selpoly }
	if ![winfo exists .wsel] {
		toplevel .wsel
		
#		menu .wsel.menuBar -tearoff 0
#		.wsel.menuBar add cascade -menu .wsel.menuBar.his -label "Selection" -underline 0
#		menu .wsel.menuBar.his -tearoff 0
#		.wsel configure -menu .wsel.menuBar
#		.wsel.menuBar.his add command -label "Help" \
#				-command { showHelp bshow_select.hlp } -underline 0
#		.wsel.menuBar.his add command -label "Clear" \
#				-command {
#					if { [info exists selpoly] } { unset selpoly }
#					$c delete selection
#				} -underline 0
#		.wsel.menuBar.his add command -label "Close" -underline 0 \
#				-command { 
#					if { [info exists selpoly] } { unset selpoly }
#					destroy .wsel
#				}
		
		frame .wsel.selectiontype
		label .wsel.selectiontype.tag -text "Type" -width 20 -anchor w
		tk_optionMenu .wsel.selectiontype.menu selectiontype "Square" "Rectangle" \
			"Cube" "Circle" "Ellipse" "Sphere"
	
		frame .wsel.size
		label .wsel.size.tag -font $helv12 -text "Size" -width 5 -anchor w
		label .wsel.size.x -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label .wsel.size.y -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w
		label .wsel.size.z -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame .wsel.num
		label .wsel.num.tag -font $helv12 -text "Number" -width 25 -anchor w
		label .wsel.num.value -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame .wsel.sum
		label .wsel.sum.tag -font $helv12 -text "Sum" -width 25 -anchor w
		label .wsel.sum.value -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame .wsel.avg
		label .wsel.avg.tag -font $helv12 -text "Average" -width 25 -anchor w
		label .wsel.avg.value -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame .wsel.std
		label .wsel.std.tag -font $helv12 -text "Standard deviation" -width 25 -anchor w
		label .wsel.std.value -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame .wsel.min
		label .wsel.min.tag -font $helv12 -text "Minimum" -width 25 -anchor w
		label .wsel.min.value -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w

		frame .wsel.max
		label .wsel.max.tag -font $helv12 -text "Maximum" -width 25 -anchor w
		label .wsel.max.value -width 10 -text "" -relief sunken -bd 1 \
				-font $helv12 -anchor w

		checkbutton .wsel.calc -text "Calculate" -variable calc_stats \
				-command { updateSelectionStats }

		frame .wsel.buttons
		button .wsel.buttons.clear -text Clear \
				-command {
					if { [info exists selpoly] } { unset selpoly }
					"$c delete selection"
				} -underline 0
		button .wsel.buttons.close -text Close \
				-command { 
					if { [info exists selpoly] } { unset selpoly }
					destroy .wsel
				}
		
		pack .wsel.selectiontype.tag .wsel.selectiontype.menu -side left -pady 5 -padx 10
		pack .wsel.size.tag .wsel.size.x .wsel.size.y .wsel.size.z -side left -pady 5 -padx 5
		pack .wsel.num.tag .wsel.num.value -side left -pady 5 -padx 20
		pack .wsel.sum.tag .wsel.sum.value -side left -pady 5 -padx 20
		pack .wsel.avg.tag .wsel.avg.value -side left -pady 5 -padx 20
		pack .wsel.std.tag .wsel.std.value -side left -pady 5 -padx 20
		pack .wsel.min.tag .wsel.min.value -side left -pady 5 -padx 20
		pack .wsel.max.tag .wsel.max.value -side left -pady 5 -padx 20
		pack .wsel.buttons.clear .wsel.buttons.close -side left -expand 1
		pack .wsel.selectiontype .wsel.size .wsel.num .wsel.sum \
			.wsel.avg .wsel.std .wsel.min .wsel.max \
			.wsel.calc .wsel.buttons -side top -pady 2

	} else {
		wm deiconify .wsel
		raise .wsel
    }
    wm title .wsel "Selection"
    wm iconname .wsel "Selection"
	.menuBar.window entryconfigure "Selection" -state normal
	bind $c <1> "mousePressed %W %x %y"
	bind $c <B1-Motion> "mouseMoved %W %x %y 1"
	bind $c <B1-ButtonRelease> "mouseReleased %W"
	bind .wsel <Control-w> {
		if { [info exists selpoly] } { unset selpoly }
		destroy .wsel
	}
}

## @brief Updates selection object.
#

proc updateSelection { } {
	global tool theimg
#	global px py
#	global mx my
	global selectiontype
	global sx1 sx2 sy1 sy2 sz1 sz2

	set c [getImageCanvas $theimg]
	$c delete selection
	
	if { [winfo exists .wsel] && $tool == "select" && $sx1 != $sx2 && $sy1 != $sy2 } {
		set wc [getControlWindow $theimg]
		set slice_num [$wc.slice.scale get]
		if { $sz1 <= $slice_num && $sz2 >= $slice_num } {
			set scale [$wc.scale.scale get]
			set h [lindex [$c cget -scrollregion] 3]
			if { $selectiontype == "Square" || $selectiontype == "Rectangle" || $selectiontype == "Cube" } {
				set cx1 [expr int($sx1 * $scale)]
				set cy1 [expr int($h - $sy1 * $scale) ]
				set cx2 [expr int($sx2 * $scale)]
				set cy2 [expr int($h - $sy2 * $scale) ]
				$c create rectangle $cx1 $cy1 $cx2 $cy2 -tags selection -outline red
			}
			if { $selectiontype == "Circle" || $selectiontype == "Ellipse" || $selectiontype == "Sphere" } {
				set t 0.0
				set z [expr ($sz1 + $sz2)*0.5 - $slice_num]
				if { $z != 0 } {
					set r [expr ($sz1 - $sz2)*0.5]
					set t [expr $z/$r]
				}
				set t [expr sqrt(1.0 - $t*$t)]
				set tn [expr 0.5*(1.0 - $t)]
				set tp [expr 0.5*(1.0 + $t)]
				set cx1 [expr $sx1 * $scale]
				set cy1 [expr $h - $sy1 * $scale ]
				set cx2 [expr $sx2 * $scale]
				set cy2 [expr $h - $sy2 * $scale ]
				set ccx1 [expr int($tp*$cx1 + $tn*$cx2)]
				set ccx2 [expr int($tn*$cx1 + $tp*$cx2)]
				set ccy1 [expr int($tp*$cy1 + $tn*$cy2)]
				set ccy2 [expr int($tn*$cy1 + $tp*$cy2)]
				$c create oval $ccx1 $ccy1 $ccx2 $ccy2 -tags selection -outline red
			}
		}
	}
}

## @brief Updates selection object statistics.
#

proc updateSelectionStats { } {
	global theimg
	global selectiontype calc_stats
	global sx1 sx2 sy1 sy2 sz1 sz2
	if { ! $calc_stats } { return }
#	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set slice_num [$wc.slice.scale get]
	set vlist "$sx1 $sy1 $sz1 $sx2 $sy2 $sz2"
	if { $selectiontype == "Square" || $selectiontype == "Rectangle" || $selectiontype == "Cube" } {
		set stats [Bimage stats $theimg $img_num 1 $vlist]
	}
	if { $selectiontype == "Circle" || $selectiontype == "Ellipse" || $selectiontype == "Sphere" } {
		set stats [Bimage stats $theimg $img_num 2 $vlist]
	}
	if { $selectiontype == "Polygon" } {
		set poly ""
		append poly "$sx1 $sy1 $slice_num "
		append poly "$sx1 $sy2 $slice_num "
		append poly "$sx2 $sy2 $slice_num "
		append poly "$sx2 $sy1 $slice_num "
#		puts $poly
		set stats [Bimage stats $theimg $img_num 3 $poly]
	}
	.wsel.num.value config -text [lindex $stats 0]
	.wsel.min.value config -text [lindex $stats 1]
	.wsel.max.value config -text [lindex $stats 2]
	.wsel.avg.value config -text [lindex $stats 3]
	.wsel.std.value config -text [lindex $stats 4]
	.wsel.sum.value config -text [lindex $stats 5]
#	puts "Expected size [expr (abs($sx1 - $sx2) + 1)*(abs($sy1 - $sy2) + 1)]"
}

proc Measure {} {
	global theimg
	set c [getImageCanvas $theimg]
	bind $c <1> "mousePressed %W %x %y"
}


