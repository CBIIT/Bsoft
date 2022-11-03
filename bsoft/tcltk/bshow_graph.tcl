##
# @file		bshow_graph.tcl
#
# @brief	Procedures to draw graphs in Bshow.
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20210110

set gwidth 500
set gheight 350
set gorix 100
set goriy 400
set color { blue green red yellow black SkyBlue2 }
set ngauss 0

## @brief Draws a histogram of the image data.
#

proc Histogram { } {
	global helv12 theimg filename
	global gorix goriy gwidth gheight
	global histplot histfit
	global color ngauss
	if { [info exists histplot] } { unset histplot }
	if { [info exists histfit] } { unset histfit }
	set w .whis
	if ![winfo exists $w] {
		toplevel $w
		
		set xmin [Bimage get $theimg min]
		set xmax [Bimage get $theimg max]
		set datatype [Bimage get $theimg datatype]
		set bins [expr int($xmax - $xmin + 1)]
		if { $bins < 10 } { set bins 10 }
		if { $bins > 500 } { set bins 500 }
		if { [string first "float" $datatype] > -1 } { set bins 250 }
#		set hist [Bimage histogram $theimg $bins]
		
#		set xscale [expr ($xmax - $xmin)/$bins]
#		for { set i 0 } { $i < $bins } { incr i } {
#			set histplot(x,$i) [expr $i*$xscale + $xmin]
#			set histplot(y,$i) [lindex $hist $i]
#		}

		menu $w.menuBar -tearoff 0
		$w.menuBar add cascade -menu $w.menuBar.his -label "Histogram" -underline 0
		menu $w.menuBar.his -tearoff 0
		$w configure -menu $w.menuBar
		$w.menuBar.his add command -label "Help" \
				-command { showHelp bshow_graph.hlp } -underline 0
		$w.menuBar.his add command -label "Reset plot" -underline 0 \
				-command { drawHistogram .whis $ngauss 1 } -accelerator "Ctrl-r"
		$w.menuBar.his add command -label "Write data" \
				-command { writePlotData histplot } -underline 0
		$w.menuBar.his add command -label "Close" -underline 0 \
				-command { 
					if { [info exists histplot] } { unset histplot }
					if { [info exists histfit] } { unset histfit }
					destroy $w
				}
		
		set h $w.plot
		canvas $h -relief raised -width [expr $gorix + $gwidth + 50] \
				-height [expr $goriy + 100] -bg white
		pack $h -side top -fill x

		frame $w.bins -bg SkyBlue2
		label $w.bins.tag -font $helv12 -text "Bins" -bg SkyBlue2
		entry $w.bins.entry -width 10 -validate key -vcmd { string is double %P }
		label $w.bins.tag2 -font $helv12 -text "Gaussians" -bg SkyBlue2
		$w.bins.entry insert 0 $bins
		tk_optionMenu $w.bins.ng ngauss 0 1 2 3 4 5
		button $w.bins.update -text "Update" -relief raised \
				-command { drawHistogram .whis $ngauss 1 }
		pack $w.bins.tag $w.bins.entry $w.bins.tag2 $w.bins.ng \
				$w.bins.update -side left -pady 10 -padx 10
		pack $w.bins -side bottom
		
		setupGraphVariables $w
		
		getLimits $w histplot
		drawHistogram $w $ngauss 1
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Histogram: $filename"
    wm iconname $w "Histogram"
#	puts "Plotting histogram (bshow_graph)"
	$h bind point <Any-Enter> { showPoint .whis histplot %x %y }
	$h bind point <Any-Leave> { showPoint .whis histplot %x %y }
	bind $h <Button-1> {
		set x1 [.whis.plot canvasx %x]
		set y1 [.whis.plot canvasy %y]
	}
	bind $h <B1-Motion> { drawRectangle .whis.plot $x1 $y1 [.whis.plot canvasx %x] [.whis.plot canvasy %y] }
	bind $h <B1-ButtonRelease> {
		setSubGraphLimits .whis $x1 $y1 \
			[.whis.plot canvasx %x] [.whis.plot canvasy %y]
		drawHistogram .whis $ngauss 0
	}
#	drawSubGraph .whis histplot $x1 $y1 \
#		[.whis.plot canvasx %x] [.whis.plot canvasy %y] "bar" "Bin" "Count"
	bind $w <Control-r> { drawHistogram .whis $ngauss 1 }
	bind $w <Return> { drawHistogram .whis $ngauss 0 }
	bind $w <Control-w> {
		if { [info exists histplot] } { unset histplot }
		if { [info exists histfit] } { unset histfit }
		destroy .whis
	}
}

## @brief Calculates the histogram
#
# @param	ngauss		Number of gaussians to fit.

proc calcHistogram { ngauss } {
	global histplot histfit gprm theimg
	global color
	if { [info exists histplot] } { unset histplot }
	if { [info exists histfit] } { unset histfit }
	set xmin [Bimage get $theimg min]
	set xmax [Bimage get $theimg max]
	set bins [.whis.bins.entry get]
	if { $bins < 2 } {
		set bins 250
		.whis.bins.entry delete 0 end
		.whis.bins.entry insert 0 $bins
	}
	set hist [Bimage histogram $theimg $bins]
	if { $ngauss > 0 } {
		set gprm [Bimage histofit $theimg $bins $ngauss]
	}
	set n [expr [llength $hist] / $bins]
	set color { LightBlue }
	if { $n == 3 } {
		set color { red green blue }
	}
#	puts "Number of channels = $n"
	set xscale [expr ($xmax - $xmin)/$bins]
	for { set i 0 } { $i < $bins } { incr i } {
		set histplot(x,$i) [expr $i*$xscale + $xmin]
		set histplot(y,$i) [lindex $hist $i]
	}
	if { $ngauss > 0 } {
		for { set i 0 } { $i < $bins } { incr i } {
			set x [expr $i*$xscale + $xmin]
			set v [expr ([lindex $gprm 1] - $x)/[lindex $gprm 2]]
			set histfit(x,$i) $x
			set histfit(y,$i) [expr [lindex $gprm 0] * exp(-0.5 * $v * $v)]
		}
	}
	for { set j 1 } { $j < $n } { incr j } {
		for { set i 0 } { $i < $bins } { incr i } {
			append histplot(y,$i) " " [lindex $hist [expr $j*$bins + $i]]
		}
	}
	if { $ngauss > 1 } {
		for { set j 1 } { $j < $ngauss } { incr j } {
			set a [lindex $gprm [expr 3*$j]]
			set d [lindex $gprm [expr 3*$j+1]]
			set s [lindex $gprm [expr 3*$j+2]]
			for { set i 0 } { $i < $bins } { incr i } {
				set x [expr $i*$xscale + $xmin]
				set v [expr ($d - $x)/$s]
				append histfit(y,$i) " " [expr $a * exp(-0.5 * $v * $v)]
			}
		}
	}
}

## @brief Draws a histogram and fit if indicated.
#
# @param	w			Window object.
# @param	ngauss		Number of gaussians to fit.
# @param	lim			Flag to recalculate limits.

proc drawHistogram { w ngauss lim } {
	global histplot histfit gprm
	global color
	calcHistogram $ngauss
	$w.plot delete all
	if { $lim } {
		getLimits $w histplot
	}
	drawAxes $w "Bin" "Count"
	plotGraph $w histplot "bar"
	if { $ngauss > 0 } {
		set color ""
		for { set j 0 } { $j < $ngauss } { incr j } {
			lappend color "red"
		}
		plotGraph $w histfit "line"
		for { set j 0 } { $j < $ngauss } { incr j } {
			set a [lindex $gprm [expr 3*$j]]
			set d [lindex $gprm [expr 3*$j+1]]
			set s [lindex $gprm [expr 3*$j+2]]
			$w.plot create text 120 [expr 20*($j+1)] \
				-text [format "y=%10.6f exp(-0.5((x-%10.6f)/%10.6f)^2)" $a $d $s] \
				-anchor w -fill black -tags rps
		}
		set R [lindex $gprm [expr 3*$ngauss]]
		$w.plot create text 120 [expr 20*($ngauss+1)] -text [format "R = %10.6f" $R] -anchor w -fill black -tags rps
	}
}

## @brief Dialog box to analyze diffraction patterns
#

proc Diffraction { } {
	global theimg
	global helv12
	global gorix goriy gwidth gheight
	global diffplot
	set c [getImageCanvas $theimg]
	if { [info exists diffplot] } { unset diffplot }

	set w .wdif
	if ![winfo exists $w] {
		toplevel $w

		menu $w.menuBar -tearoff 0
		$w.menuBar add cascade -menu $w.menuBar.dif -label "Diffraction" -underline 0
		menu $w.menuBar.dif -tearoff 0
		$w configure -menu $w.menuBar
		$w.menuBar.dif add command -label "Help" \
				-command { showHelp bshow_diffraction.hlp } -underline 0
		$w.menuBar.dif add command -label "Reset plot" -underline 0 \
				-command { drawGraph $w diffplot "line" "Frequency" "Intensity" } -accelerator "Ctrl-r"
		$w.menuBar.dif add command -label "Write data" \
				-command { writePlotData diffplot } -underline 0
		$w.menuBar.dif add command -label "Close" \
				-command "destroy $w" -underline 0
		
		set g $w.plot
		canvas $g -relief raised -width [expr $gorix + $gwidth + 50] \
				-height [expr $goriy + 100]
		pack $g -side top -fill x
		
		frame $w.origin
		label $w.origin.xtag -font $helv12 -text "ox"
		entry $w.origin.x -width 10 -validate key -vcmd { string is double %P }
		label $w.origin.ytag -font $helv12 -text "oy"
		entry $w.origin.y -width 10 -validate key -vcmd { string is double %P }
		label $w.origin.ztag -font $helv12 -text "oz"
		entry $w.origin.z -width 10 -validate key -vcmd { string is double %P }
		label $w.origin.angletag -font $helv12 -text "angle"
		entry $w.origin.angle -width 10 -validate key -vcmd { string is double %P }
		button $w.origin.find -text "Find origin" -relief raised \
				-command { calcRadialAverage .wdif diffplot 1 }
		pack $w.origin $w.origin.xtag $w.origin.x $w.origin.ytag \
				$w.origin.y $w.origin.ztag $w.origin.z \
				$w.origin.angletag $w.origin.angle \
				$w.origin.find -side left -pady 4 -padx 4
		
		frame $w.upd_mask
		button $w.upd_mask.update -text "Update" -relief raised \
				-command { calcRadialAverage .wdif diffplot 0 }
#		label $w.upd_mask.mask_label -text "Mask"
		checkbutton $w.upd_mask.mask_check -text "Mask" -variable use_mask
		entry $w.upd_mask.mask_entry -width 50
		button $w.upd_mask.mask_button -text "Browse ..." -relief raised \
				-command {
					set maskfile [tk_getOpenFile -filetypes $filetypes -parent $w]
					$w.upd_mask.mask_entry delete 0 end
					$w.upd_mask.mask_entry insert 0 $maskfile
				}
		pack $w.upd_mask.update $w.upd_mask.mask_check \
				$w.upd_mask.mask_entry $w.upd_mask.mask_button -side left -pady 4 -padx 4

		pack $w.origin $w.upd_mask -side bottom -fill x
		
		set ori [Bimage get $theimg origin]
		$w.origin.x insert 0 [lindex $ori 0]
		$w.origin.y insert 0 [lindex $ori 1]
		$w.origin.z insert 0 [lindex $ori 2]
		$w.origin.angle insert 0 0

		setupGraphVariables $w

		calcRadialAverage $w diffplot 0
	} else {
		wm deiconify $w
		raise $w
		set g $w.plot
    }
    wm title $w "Diffraction"
    wm iconname $w "Diffraction"
	bind $c <1> "setDiffOrigin .wdif %x %y"
	$g bind point <Any-Enter> { showPoint .wdif diffplot %x %y }
	$g bind point <Any-Leave> { showPoint .wdif diffplot %x %y }
	bind $g <Button-1> {
		set x1 [.wdif.plot canvasx %x]
		set y1 [.wdif.plot canvasy %y]
	}
	bind $g <B1-Motion> { drawRectangle .wdif.plot $x1 $y1 [.wdif.plot canvasx %x] [.wdif.plot canvasy %y] }
	bind $g <B1-ButtonRelease> { drawSubGraph .wdif diffplot $x1 $y1 \
		[.wdif.plot canvasx %x] [.wdif.plot canvasy %y] "line" "Frequency" "Intensity" }
	bind $w <Control-r> { drawGraph .wdif diffplot "line" "Frequency" "Intensity" }
	bind $w <Return> {
		plotGraph .wdif diffplot "line"
		drawAxes "Frequency" "Intensity"
	}
	bind $w <Control-w> {
		if { [info exists diffplot] } { unset diffplot }
		destroy .wdif
	}
}

## @brief Sets the origin from canvas coordinates
#
# @param	w			Plot window object.
# @param	x
# @param	y			Position in canvas coordinates.

proc setDiffOrigin { w x y } {
	global theimg
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	if ![winfo exists $w] {
		$c delete center
		return
	}
	set scale [$wc.scale.scale get]
	set slice_num [$wc.slice.scale get]
	set xr [expr [$c canvasx $x] / $scale ]
	set yr [yFlip [expr [$c canvasy $y] / $scale] ]
	updateOrigin $w $xr $yr $slice_num
}

## @brief Updates the origin
#
# @param	w			Plot window object.
# @param	x
# @param	y
# @param	z			Origin.

proc updateOrigin { w x y z } {
	global PI theimg
	if ![winfo exists $w] {
		$c delete center
		return
	}
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	$w.origin.x delete 0 end
	$w.origin.x insert 0 $x
	$w.origin.y delete 0 end
	$w.origin.y insert 0 $y
	$w.origin.z delete 0 end
	$w.origin.z insert 0 $z
	set a [expr [$w.origin.angle get] * $PI/180.0]
	set scale [$wc.scale.scale get]
	set height [Bimage get $theimg height]
	set rad [expr 10*$scale]
	set len [expr $height*$scale/4]
	set cx [expr $x*$scale]
	set cy [expr ($height - $y - 1)*$scale]
	set cxmin [expr $cx - $rad]
	set cxmax [expr $cx + $rad]
	set cymin [expr $cy - $rad]
	set cymax [expr $cy + $rad]
	set ca [expr $len * cos($a)]
	set sa [expr $len * sin($a)]
	$c delete center
	$c create oval $cxmin $cymin $cxmax $cymax -outline yellow -width 1 -tags center
	$c create line $cx $cy [expr $cx + $ca] [expr $cy - $sa] -width 1 -fill red -tags center
	$c create line $cx $cy [expr $cx - $sa] [expr $cy - $ca] -width 1 -fill blue -tags center
}

## @brief Calculates the radial average of the first image
#
# @param	w			Plot window object.
# @param	plotdata	Array of data.
# @param	findori		Flag to find diffraction origin

proc calcRadialAverage { w plotdata findori } {
	global PI theimg
	global color
	global use_mask
	set c [getImageCanvas $theimg]
	upvar $plotdata diffplot
	set ox [$w.origin.x get]
	set oy [$w.origin.y get]
	set oz [$w.origin.z get]
	set a [expr [$w.origin.angle get] * $PI/180.0]
	set maskname [$w.upd_mask.mask_entry get]
	if { $findori > 0 } {
		if { $use_mask } {
			Bimage diffori $theimg $ox $oy $oz $maskname
		} else {
			Bimage diffori $theimg $ox $oy $oz
		}
		set ori [Bimage get $theimg origin]
		set ox [lindex $ori 0]
		set oy [lindex $ori 1]
		set oz [lindex $ori 2]
	}
	if { $use_mask } {
		set raddiff [Bimage radial $theimg $ox $oy $oz $a $maskname]
	} else {
		set raddiff [Bimage radial $theimg $ox $oy $oz $a]
	}
	set n [expr [llength $raddiff] / 3]
	set color { black red blue }
#	puts $raddiff
#	puts "Radial plot length = $n"
	set ori [Bimage get $theimg origin]
	updateOrigin $w [lindex $ori 0] [lindex $ori 1] [lindex $ori 2]
#	set xscale 1.0
	for { set i 0 } { $i < $n } { incr i } {
#		set diffplot(x,$i) [expr $i*$xscale + $xmin]
		set diffplot(x,$i) $i
		set diffplot(y,$i) [lindex $raddiff $i]
	}
	for { set j 1 } { $j < 3 } { incr j } {
		for { set i 0 } { $i < $n } { incr i } {
			append diffplot(y,$i) " " [lindex $raddiff [expr $j*$n + $i]]
		}
	}
#	for { set i 0 } { $i < $n } { incr i } {
#		puts "$i [lindex $raddiff $i] [lindex $raddiff [expr $n + $i]] [lindex $raddiff [expr 2*$n + $i]]"
#	}
	drawGraph $w diffplot "line" "Frequency" "Intensity"
}

## @brief Sets up the point position and limits variables in the window
#
# @param	w				Window object.

proc setupGraphVariables { w } {
	global helv12
	frame $w.point -bg green
	label $w.point.xtag -font $helv12 -text "x" -bg green
	label $w.point.x -text "0" -relief sunken -bd 1 -width 10 \
			-font $helv12 -anchor w
	label $w.point.ytag -font $helv12 -text "y" -bg green
	label $w.point.y -text "0" -relief sunken -bd 1 -width 10 \
			-font $helv12 -anchor w
	pack $w.point.xtag $w.point.x $w.point.ytag $w.point.y \
			-side left -pady 10 -padx 6

	frame $w.lim -bg cyan
	label $w.lim.x_tag -font $helv12 -text "X limits" -bg cyan
	entry $w.lim.x_min -width 8 -validate key -vcmd { string is double %P }
	entry $w.lim.x_max -width 8 -validate key -vcmd { string is double %P }
	label $w.lim.y_tag -font $helv12 -text "Y limits" -bg cyan
	entry $w.lim.y_min -width 8 -validate key -vcmd { string is double %P }
	entry $w.lim.y_max -width 8 -validate key -vcmd { string is double %P }
	pack $w.lim.x_tag $w.lim.x_min $w.lim.x_max $w.lim.y_tag \
			$w.lim.y_min $w.lim.y_max -side left -pady 6 -padx 6 -fill x

	pack $w.point $w.lim -side left -fill x
}

## @brief Plots a graph within specific limits.
#
# @param	w			Window object.
# @param	plotdata	Array of data.
# @param	type		type of graph (oval, bar, line)

proc plotGraph { w plotdata type } {
	upvar $plotdata plot
	global gorix goriy gwidth gheight
	global helv12 helv14
	global color
	set gc $w.plot
#	set color { blue green red yellow black SkyBlue2 }
#	set color { red green blue yellow black SkyBlue2 }
#	$gc delete all
	# Default type is oval
	set typenum 0
	if { [string first "bar" $type] > -1 } { set typenum 1 }
	if { [string first "line" $type] > -1 } { set typenum 2 }
	set number [expr int( [array size plot] / 2 ) ]
	if { $number < 2 } { return }
	set nsets [llength $plot(y,0)]
	set xmin [$w.lim.x_min get]
	set xmax [$w.lim.x_max get]
	set ymin [$w.lim.y_min get]
	set ymax [$w.lim.y_max get]
	set xrange [expr $xmax - $xmin]
	set yrange [expr $ymax - $ymin]
	set xscale [expr $gwidth * 1.0 / $xrange]
	set yscale [expr $gheight * 1.0 / $yrange]
	set bar_width [expr ($plot(x,1) - $plot(x,0))*$xscale/$nsets]
	set bwh [expr $bar_width /2.0]
#	puts "$gorix $goriy $gwidth $gheight"
#	puts "$xmin $xmax $ymin $ymax $xrange $yrange $xscale $yscale"
	# Plot the data
	for { set j 0 } { $j < $nsets } { incr j } {
		set theline ""
		set offset [expr $bar_width * ($j - ($nsets - 1)/2.0)]
		for { set i 0 } { $i < $number } { incr i } {
			set vx $plot(x,$i)
			if { $vx >= $xmin && $vx <= $xmax } {
				set x [expr $gorix + ($vx - $xmin)*$xscale]
				set vy [lindex $plot(y,$i) $j]
				if { $typenum == 1 } {
					if { $vy > $ymax } { set vy $ymax }
				}
				set y [expr $goriy - ($vy - $ymin)*$yscale]
				if { $typenum == 2 } {
					append theline " " $x " " $y
				} else { if { $vy >= $ymin && $vy <= $ymax } {
#					puts "$x $y"
					if { $typenum == 0 } {
						set item [$gc create oval [expr $x-3] [expr $y-3] \
							[expr $x+3] [expr $y+3] -width 1 -outline {} \
							-fill [lindex $color $j]]
					} else {
						set x [expr $x + $offset]
						set item [$gc create rectangle [expr $x-$bwh] $goriy \
							[expr $x+$bwh] $y -width 1 -outline {} \
							-fill [lindex $color $j]]
					}
					$gc addtag point withtag $item
				} }
			}
		}
		if { $typenum == 2 } {
			$gc create line $theline -fill [lindex $color $j] -smooth yes -width 1 -tags point
		}
#		puts "Plotting line $j with color [lindex $color $j]"
	}
}

## @brief Draws the axes for a graph.
#
# @param	w			Window object.
# @param	xlab
# @param	ylab		axis labels
proc drawAxes { w xlab ylab } {
	global gorix goriy gwidth gheight
	global helv12 helv14
	set gc $w.plot
	set xmin [$w.lim.x_min get]
	set xmax [$w.lim.x_max get]
	set ymin [$w.lim.y_min get]
	set ymax [$w.lim.y_max get]
	set xrange [expr $xmax - $xmin]
	set yrange [expr $ymax - $ymin]
	set ex 0
	set e1 0
	if { $xmax != 0 } { set ex [expr floor(log10(abs($xmax)))] }
	if { $xmin != 0 } { set e1 [expr floor(log10(abs($xmin)))] }
	if { abs($e1) > abs($ex) } { set ex $e1 }
	set xrange [expr $xrange / pow(10,$ex)]
	set xmin [expr $xmin / pow(10,$ex)]
	set ey 0
	set e1 0
	if { $ymax != 0 } { set ey [expr floor(log10(abs($ymax)))] }
	if { $ymin != 0 } { set e1 [expr floor(log10(abs($ymin)))] }
	if { abs($e1) > abs($ey) } { set ey $e1 }
	set yrange [expr $yrange / pow(10,$ey)]
	set ymin [expr $ymin / pow(10,$ey)]
	# Draw the axes and labels
	$gc create line $gorix $goriy [expr $gorix + $gwidth] $goriy -width 2
	$gc create line $gorix $goriy $gorix [expr $goriy - $gheight] -width 2
	for { set i 0 } { $i <= 10 } { incr i } {
		set x [expr $gorix + $i*$gwidth/10]
		$gc create line $x $goriy $x [expr $goriy - 5] -width 1
		set t [format "%6.3g" [expr $xrange*$i/10.0 + $xmin]]
		$gc create text $x [expr $goriy + 4] -text $t -anchor n -font $helv12
	}
	set ex [expr int($ex)]
	set ey [expr int($ey)]
	set x [expr $gorix + $gwidth/2]
	set y [expr $goriy + 40]
	$gc create text $x $y -text $xlab -anchor e -font $helv14
	$gc create text $x $y -text " (10   )" -anchor w -font $helv14
	$gc create text [expr $x + 25] [expr $y - 5] -text "$ex" -anchor w -font $helv12
	for { set i 0 } { $i <= 10 } { incr i } {
		set y [expr $goriy - $i*$gheight/10]
		$gc create line $gorix $y [expr $gorix + 5] $y -width 1
		set t [format "%6.3g" [expr $yrange*$i/10.0 + $ymin]]
		$gc create text [expr $gorix - 4] $y -text $t -anchor e -font $helv12
	}
	set x [expr $gorix - 70]
	set y [expr $goriy - $gheight/2]
	$gc create text $x $y -text $ylab -anchor c -font $helv12
	$gc create text $x [expr $y + 20] -text "(10   )" -anchor c -font $helv14
	$gc create text [expr $x + 5] [expr $y + 15] -text "$ey" -anchor w -font $helv12
}

## @brief Gets the limits from the input data.
#
# @param	w				Window object.
# @param	plotdata		Array of data.

proc getLimits { w plotdata } {
	upvar $plotdata plot
#	puts "$plot(x,0) $plot(y,0)"
	set number [expr int( [array size plot] / 2 ) ]
	if { $number < 2 } { return }
	set xmin $plot(x,0)
	set xmax $plot(x,0)
	set ymin [lindex $plot(y,0) 0]
	set ymax [lindex $plot(y,0) 0]
	for { set i 0 } { $i < $number } { incr i } {
		if { $xmin > $plot(x,$i) } { set xmin $plot(x,$i) }
		if { $xmax < $plot(x,$i) } { set xmax $plot(x,$i) }
		for { set j 0 } { $j < [llength $plot(y,$i)] } { incr j } {
			if { $ymin > [lindex $plot(y,$i) $j] } { set ymin [lindex $plot(y,$i) $j] }
			if { $ymax < [lindex $plot(y,$i) $j] } { set ymax [lindex $plot(y,$i) $j] }
		}
	}
	if { $xmin == $xmax } { set xmax [expr $xmin + 1] }
	if { $ymin == $ymax } { set ymax [expr $ymin + 1] }
	updateLimits $w $xmin $xmax $ymin $ymax
}

## @brief Updates the entry boxes for the plot limits.
#
# @param	w			Window object.
# @param	xmin
# @param	xmax		x range.
# @param	ymin
# @param	ymax		y range.

proc updateLimits { w xmin xmax ymin ymax } {
	$w.lim.x_min delete 0 end
	$w.lim.x_min insert 0 $xmin
	$w.lim.x_max delete 0 end
	$w.lim.x_max insert 0 $xmax
	$w.lim.y_min delete 0 end
	$w.lim.y_min insert 0 $ymin
	$w.lim.y_max delete 0 end
	$w.lim.y_max insert 0 $ymax
}

## @brief Sets the subgraph limits.
#
# @param	w			Window object.
# @param	x1
# @param	y1
# @param	x2
# @param	y2			Canvas coordinates of selected region

proc setSubGraphLimits { w x1 y1 x2 y2 } {
	global gorix goriy gwidth gheight
	if { $x1 == $x2 } { return }
	if { $y1 == $y2 } { return }
	if { $x2 < $x1 } {
		set t $x1
		set x1 $x2
		set x2 $t
	}
	if { $y2 > $y1 } {	# The y-scale is upside down
		set t $y1
		set y1 $y2
		set y2 $t
	}
	set xmin [$w.lim.x_min get]
	set xmax [$w.lim.x_max get]
	set ymin [$w.lim.y_min get]
	set ymax [$w.lim.y_max get]
#	if { $y1 > $goriy } { set y1 $goriy }
	set xrange [expr ($xmax - $xmin) * 1.0 / $gwidth]
	set yrange [expr ($ymax - $ymin) * 1.0 / $gheight]
	set xmax [expr $xmin + ($x2 - $gorix) * $xrange]
	set xmin [expr $xmin + ($x1 - $gorix) * $xrange]
	set ymax [expr $ymin + ($goriy - $y2) * $yrange]
	set ymin [expr $ymin + ($goriy - $y1) * $yrange]
	if { $ymin < 0 } { set ymin 0 }
	updateLimits $w $xmin $xmax $ymin $ymax
}


## @brief Draws the graph when new data is read in.
#
# @param	w			Window object.
# @param	plotdata	Array of data.
# @param	type		Plot type.
# @param	xlab
# @param	ylab		axis labels

proc drawGraph { w plotdata type xlab ylab } {
	upvar $plotdata plot
	$w.plot delete all
	getLimits $w plot
	plotGraph $w plot $type
	drawAxes $w $xlab $ylab
}

## @brief Draws the graph when new data is read in.
#
# @param	w			Window object.
# @param	plotdata	Array of data.
# @param	x1
# @param	y1
# @param	x2
# @param	y2			Canvas coordinates of selected region
# @param	type		Plot type.
# @param	xlab
# @param	ylab		axis labels

proc drawSubGraph { w plotdata x1 y1 x2 y2 type xlab ylab } {
	upvar $plotdata plot
	global gorix goriy gwidth gheight
	if { $x1 == $x2 } { return }
	if { $y1 == $y2 } { return }
	if { $x2 < $x1 } {
		set t $x1
		set x1 $x2
		set x2 $t
	}
	if { $y2 > $y1 } {	# The y-scale is upside down
		set t $y1
		set y1 $y2
		set y2 $t
	}
	set xmin [$w.lim.x_min get]
	set xmax [$w.lim.x_max get]
	set ymin [$w.lim.y_min get]
	set ymax [$w.lim.y_max get]
#	if { $y1 > $goriy } { set y1 $goriy }
	set xrange [expr ($xmax - $xmin) * 1.0 / $gwidth]
	set yrange [expr ($ymax - $ymin) * 1.0 / $gheight]
	set xmax [expr $xmin + ($x2 - $gorix) * $xrange]
	set xmin [expr $xmin + ($x1 - $gorix) * $xrange]
	set ymax [expr $ymin + ($goriy - $y2) * $yrange]
	set ymin [expr $ymin + ($goriy - $y1) * $yrange]
	if { $ymin < 0 } { set ymin 0 }
	updateLimits $w $xmin $xmax $ymin $ymax
	plotGraph $w plot $type
	drawAxes $w $xlab $ylab
}

## @brief Gets the limits from the input data.
#
# @param	w			Window object.
# @param	plotdata	Array of data.
# @param	x
# @param	y			Point coordinates in the window.

proc showPoint { w plotdata x y } {
	upvar $plotdata plot
	global gorix goriy gwidth gheight
	set number [expr int( [array size plot] / 2 ) ]
	if { $number < 2 } { return }
	set xmin [$w.lim.x_min get]
	set xmax [$w.lim.x_max get]
	set ymin [$w.lim.y_min get]
	set ymax [$w.lim.y_max get]
	set x [expr ($x - $gorix)*($xmax - $xmin)/$gwidth + $xmin]
	set y [expr ($goriy - $y)*($ymax - $ymin)/$gheight + $ymin]
	set dxmin 1e30
	for { set i 0 } { $i < $number } { incr i } {
		set dx [expr abs($plot(x,$i) - $x)]
		if { $dxmin > $dx } {
			set dxmin $dx
			set vx $plot(x,$i)
			set dymin 1e30
			for { set j 0 } { $j < [llength $plot(y,$i)] } { incr j } {
				set dy [expr abs([lindex $plot(y,$i) $j] - $y)]
				if { $dymin > $dy } {
					set dymin $dy
					set vy [lindex $plot(y,$i) $j]
				}
			}
		}
	}
#	puts "$vx $vy"
	$w.point.x config -text "$vx"
	$w.point.y config -text "$vy"
}

## @brief Draws a rectangle on the canvas.
#
# @param	c		Canvas object.
# @param	x1
# @param	y1
# @param	x2
# @param	y2		Canvas coordinates of rectangle.

proc drawRectangle { c x1 y1 x2 y2 } {
	$c delete rect
	$c create rectangle $x1 $y1 $x2 $y2 -tags rect -outline red
}

## @brief Writes the plot data to a text file.
#
# @param	plotdata		Array of data.

proc writePlotData { plotdata } {
	upvar $plotdata plot
	set datafile [tk_getSaveFile -defaultextension .txt]
	if { [string length $datafile] < 1 } { return }
	set fd [open $datafile w]
	set number [expr int( [array size plot] / 2 ) ]
	for { set i 0 } { $i < $number } { incr i } {
		puts $fd "$plot(x,$i) $plot(y,$i)"
	}
}
