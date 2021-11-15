##
# @file		bshow_scales.tcl
#
# @brief	Procedures to manage scales for image display
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20130820

set window_scale 1.0
set std_scale 5.0

## @brief Configures scales for the image window
#
# These scales are used as stores for the slice and image numbers
# as well as the scaling of the image itself
#
# @param	w 	 	 	Control window to put scales in.

proc setupScales { w } {
	global helv12 window_scale std_scale aniproc
	
	frame $w.image
	label $w.image.label -font $helv12 -text "Image" -width 6 -anchor w
	scale $w.image.scale -orient horizontal -from 0 -to 0 \
		-command { setImage }
	pack $w.image.label -side left
	pack $w.image.scale -side bottom -fill x
	$w.image.scale set 0

	frame $w.slice
	label $w.slice.label -font $helv12 -text "Slice" -width 6 -anchor w
	scale $w.slice.scale -orient horizontal -from 0 -to 0 \
		-command { setSlice }
	pack $w.slice.label -side left
	pack $w.slice.scale -side bottom -fill x
	$w.slice.scale set 0

	frame $w.scale
	label $w.scale.label -font $helv12 -text "Scale" -width 6 -anchor w
	pack $w.scale.label -side left -pady 2
	scale $w.scale.scale -orient horizontal -from 0.1 -to 10 \
		-command { setScale } -resolution 0.1
	pack $w.scale.label -side left
	pack $w.scale.scale -side bottom -fill x
	$w.scale.scale set $window_scale

	frame $w.min
	label $w.min.label -font $helv12 -text "Min" -width 6 -anchor w
	scale $w.min.scale -orient horizontal -from 0 -to 0 \
		-command { setMin }
	pack $w.min.label -side left
	pack $w.min.scale -side bottom -fill x
	
	frame $w.max
	label $w.max.label -font $helv12 -text "Max" -width 6 -anchor w
	scale $w.max.scale -orient horizontal -from 0 -to 0 \
		-command { setMax }
	pack $w.max.label -side left
	pack $w.max.scale -side bottom -fill x

	setupEntry $w.step "Range step" double 0.01 "Step size for minimum and maximum scales"
	bind $w.step.e <Return> reconfigureScales

	frame $w.as
	checkbutton $w.as.autoscale -text "Autoscale" -command { autoScale }
	checkbutton $w.as.imageonly -text "Image only" -command { autoScale }
	entry $w.as.scale -width 6 -validate key -vcmd { string is double %P }
	$w.as.scale insert 0 $std_scale
	pack $w.as.autoscale $w.as.scale $w.as.imageonly -side left -pady 2 -padx 2 -anchor w

	frame $w.ani
	button $w.ani.animate_button -text Animate \
		-relief raised -command do_animate_slices
	label $w.ani.animate_tag -font $helv12 -text "Speed"
	entry $w.ani.animate -width 6 -validate key -vcmd { string is int %P }
	$w.ani.animate insert 0 15
	bind $w.ani.animate <Return> do_animate_slices
	pack $w.ani.animate_button $w.ani.animate_tag $w.ani.animate -side left -pady 2 -padx 2 -anchor w
	set aniproc 0		

	bind . <Left> { .slice.scale set [expr [.slice.scale get] - 1] }
	bind . <Right> { .slice.scale set [expr [.slice.scale get] + 1] }
	bind . <Up> { .image.scale set [expr [.image.scale get] + 1] }
	bind . <Down> { .image.scale set [expr [.image.scale get] - 1] }
	bind $w.as.scale <Return> { autoScale }
}

## @brief Resets the scales when the image changes

proc reconfigureScales { } {
	global theimg
	global window_scale montover
	if ![Bimage exists $theimg] { return }
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	if { [info exists montover(active)] && $montover(active) } {
		set width [expr $montover(ncol) * ($montover(nx) + 2*$montover(pad))]
		set height [expr $montover(nrow) * ($montover(ny) + 2*$montover(pad))]
		set nimages 1
		set nslices 1
	} else {
		set width [Bimage get $theimg width]
		set height [Bimage get $theimg height]
		set nimages [Bimage get $theimg nimages]
		set nslices [Bimage get $theimg nslices]
	}
	set width [expr $width * $window_scale]
	set height [expr $height * $window_scale]
	$wc.image.scale config -to [expr $nimages - 1] -tickinterval [expr $nimages - 1]
	$wc.slice.scale config -to [expr $nslices - 1] -tickinterval [expr $nslices - 1]
	reconfigureMinMax
	set size "0 0 $width $height"
	$c configure -scrollregion $size
}

## @brief Resets the min and max scales when the range step size changes

proc reconfigureMinMax { } {
	global theimg autoscale
	set wc [getControlWindow $theimg]
	set range_step [$wc.step.e get]
	set min [Bimage get $theimg min]
	set max [Bimage get $theimg max]
	Bimage set $theimg show_min $min
	Bimage set $theimg show_max $max
	$wc.min.scale config -from $min -to $max -resolution $range_step
	$wc.max.scale config -from $min -to $max -resolution $range_step
	if $autoscale {
		autoScale
	} else {
		$wc.min.scale set $min
		$wc.max.scale set $max
	}
}

## @brief Resets the range step size to the default value

proc reconfigureRangeStep { } {
	global theimg imgtype
	if ![Bimage exists $theimg] { return }
	set wc [getControlWindow $theimg]
	set min [Bimage get $theimg min]
	set max [Bimage get $theimg max]
	set step 1000
	if { $imgtype == "ps" || $imgtype == "ft" } {
		set step 100000
	}
	set range_step [expr ($max - $min) / $step]
	$wc.step.e delete 0 end
	$wc.step.e insert 0 $range_step
}

## @brief A dialog to type in a scale factor and new min and max values.
#

proc setScaleDialog { } {
	global helv12
	set wd .scaleDialog
	catch {destroy $wd}
	toplevel $wd
	wm title $wd "Set display parameters"
	wm iconname $wd "Set display"

	frame $wd.panel
	pack $wd.panel -side top -fill both
	
	frame $wd.panel.image
	label $wd.panel.image.tag -text "Image" -width 6
	entry $wd.panel.image.entry -width 10 -validate key -vcmd { string is integer %P }
	$wd.panel.image.entry insert 0 [.image.scale get]
	pack $wd.panel.image.tag $wd.panel.image.entry \
		-side left -fill both -expand true

	frame $wd.panel.slice
	label $wd.panel.slice.tag -text "Slice" -width 6
	entry $wd.panel.slice.entry -width 10 -validate key -vcmd { string is integer %P }
	$wd.panel.slice.entry insert 0 [.slice.scale get]
	pack $wd.panel.slice.tag $wd.panel.slice.entry \
		-side left -fill both -expand true

	frame $wd.panel.scale
	label $wd.panel.scale.tag -text "Scale" -width 6
	entry $wd.panel.scale.entry -width 10 -validate key -vcmd { string is double %P }
	$wd.panel.scale.entry insert 0 [.scale.scale get]
	pack $wd.panel.scale.tag $wd.panel.scale.entry \
		-side left -fill both -expand true

	frame $wd.panel.min
	label $wd.panel.min.tag -text "Min" -width 6
	entry $wd.panel.min.entry -width 10 -validate key -vcmd { string is double %P }
	$wd.panel.min.entry insert 0 [.min.scale get]
	pack $wd.panel.min.tag $wd.panel.min.entry \
		-side left -fill both -expand true

	frame $wd.panel.max
	label $wd.panel.max.tag -text "Max" -width 6
	entry $wd.panel.max.entry -width 10 -validate key -vcmd { string is double %P }
	$wd.panel.max.entry insert 0 [.max.scale get]
	pack $wd.panel.max.tag $wd.panel.max.entry \
		-side left -fill both -expand true

	pack $wd.panel.image $wd.panel.slice $wd.panel.scale $wd.panel.min $wd.panel.max -side top -padx 5 -pady 5
	
	frame $wd.buttons
	pack $wd.buttons -side bottom -fill x -pady 2m
	button $wd.buttons.ok -text OK -command "setScaleOK $wd"
	button $wd.buttons.cancel -text Cancel -command "destroy $wd"
	pack $wd.buttons.ok $wd.buttons.cancel -side left -expand 1
	
	bind $wd <Return> "setScaleOK $wd"
}

## @brief Sets the scales given in the dialog.
#
# @param	wd 	 	 	Dialog window.

proc setScaleOK { wd } {
	global theimg
	set wc [getControlWindow $theimg]
	set n [$wd.panel.image.entry get]
	set z [$wd.panel.slice.entry get]
	set scale [$wd.panel.scale.entry get]
	set min [$wd.panel.min.entry get]
	set max [$wd.panel.max.entry get]
	destroy $wd
	$wc.image.scale set $n
	$wc.slice.scale set $z
	$wc.scale.scale set $scale
	$wc.min.scale set $min
	$wc.max.scale set $max
}

## @brief Sets the scale for display.
#
# @param	scale 	 	Display scale stored in the scale scale

proc setScale { scale } {
	global theimg
	global window_width window_height window_scale
	global img_load_phase montover
#	if { $window_scale == $scale } { return }
	set window_scale $scale
	if ![Bimage exists theimg] { return }
	set c [getImageCanvas $theimg]
	if { [info exists montover(active)] && $montover(active) } {
		set width [expr $montover(ncol) * ($montover(nx) + 2*$montover(pad))]
		set height [expr $montover(nrow) * ($montover(ny) + 2*$montover(pad))]
	} else {
		set width [Bimage get $theimg width]
		set height [Bimage get $theimg height]
	}
	set width [expr $width * $window_scale]
	set height [expr $height * $window_scale]
	set size "0 0 $width $height"
#	$c configure -scrollregion $size
	$c configure -width $window_width -height $window_height -scrollregion $size
#	puts "set scale to $scale"
	if { !$img_load_phase } { Update 1 }
}

## @brief Loads the requested slice in the current file into a photo image for display.
#
# @param	slice_num 	Current slice number stored in the slice scale

proc setSlice { slice_num } {
	global theimg img_load_phase
	if ![Bimage exists $theimg] { return }
	set wc [getControlWindow $theimg]
	set nslices [Bimage get $theimg nslices]
#	set img_num [$wc.image.scale get]
	if { $nslices < 1 } { return }
#	if { $img_num < 0 } { set img_num 0 }
	if { $slice_num >= $nslices } { set slice_num 0 }
	if { [winfo exists .wmag] } {
		.wmag.q.slice set [$wc.slice.scale get]
	}
#	puts "set slice to $slice_num"
	if { !$img_load_phase } { Update 1 }
}

## @brief Loads the requested image in the current file into a photo image for display.
#
# @param	img_num 	Current image number stored in the image scale

proc setImage { img_num } {
	global theimg img_load_phase imageonly
	if ![Bimage exists theimg] { return }
	set nimages [Bimage get $theimg nimages]
	if { $nimages < 1 } { return }
	if { $img_num < 0 } { set img_num 0 }
	if { $img_num >= $nimages } { set img_num 0 }
#	if $imageonly { reconfigureMinMax }
	if $imageonly { autoScale }
#	puts "set image to $img_num"
	if { !$img_load_phase } { Update 1 }
}

## @brief Sets the minimum display value.
#
# @param	smin 			New display minimum.

proc setMin { smin } {
	global theimg mode
	global img_load_phase autoscale
	if ![Bimage exists $theimg] { return }
	set wc [getControlWindow $theimg]
	set step [$wc.step.e get]
	if { [expr 2*abs([Bimage get $theimg show_min] - $smin)] < $step } { return }
	if { [$wc.min.scale cget -from] || [$wc.min.scale cget -to] } {
		set autoscale 0
		set scale [$wc.scale.scale get]
		set slice_num [$wc.slice.scale get]
		set img_num [$wc.image.scale get]
#		puts "set show_min to $smin"
		Bimage set $theimg show_min $smin
		if { !$img_load_phase } {
			Bimage show $theimg $img_num $slice_num $scale $mode
		}
	}
}

## @brief Sets the maximum display value.
#
# @param	smax 			New display maximum.

proc setMax { smax } {
	global theimg mode
	global img_load_phase autoscale
	if ![Bimage exists $theimg] { return }
	set wc [getControlWindow $theimg]
	set step [$wc.step.e get]
	if { [expr 2*abs([Bimage get $theimg show_max] - $smax)] < $step } { return }
	if { [$wc.max.scale cget -from] || [$wc.max.scale cget -to] } {
		set autoscale 0
		set scale [$wc.scale.scale get]
		set slice_num [$wc.slice.scale get]
		set img_num [$wc.image.scale get]
#		puts "set show_max to $smax"
		Bimage set $theimg show_max $smax
		if { !$img_load_phase } {
			Bimage show $theimg $img_num $slice_num $scale $mode
		}
	}
}

## @brief Sets reasonable minimum and maximum display values.
#

proc autoScale {} {
	global theimg autoscale std_scale imageonly
	if ![Bimage exists $theimg] { return }
	set wc [getControlWindow $theimg]
	
#	set autoscale [expr 1 - $autoscale]
	
	if $imageonly {
		set img_num [$wc.image.scale get]
		set min [Bimage get $theimg min $img_num]
		set max [Bimage get $theimg max $img_num]
		set avg [Bimage get $theimg average $img_num]
		set std [Bimage get $theimg standard_deviation $img_num]
	} else {
		set min [Bimage get $theimg min]
		set max [Bimage get $theimg max]
		set avg [Bimage get $theimg average]
		set std [Bimage get $theimg standard_deviation]
	}
	
	set std_scale [$wc.as.scale get]
		
	set smin [expr $avg - $std_scale*$std]
	set smax [expr $avg + $std_scale*$std]
	if { $smin >= $smax } {
		set smin [expr $avg - 1]
		set smax [expr $avg + 1]
	}
	if { $smin < $min } {
		set smin $min
	}
	if { $smax > $max } {
		set smax $max
	}

	$wc.min.scale set $smin
	$wc.max.scale set $smax

	update
	
	set autoscale 1

#	puts "autoscale = $autoscale"
}

## @brief Switches the minimum and maximum display values.
#

proc invertContrast {} {
	global theimg
	if ![Bimage exists $theimg] { return }
	set wc [getControlWindow $theimg]
	set smin [Bimage get $theimg show_min]
	set smax [Bimage get $theimg show_max]
	$wc.min.scale set $smax
	$wc.max.scale set $smin
}

## @brief Pages through slices.
#

proc do_animate_slices {} {
	global theimg aniproc
	set wc [getControlWindow $theimg]
	if { [Bimage get $theimg nimages] > [Bimage get $theimg nslices] } {
		set scale_to_move images
	} else {
		set scale_to_move slices
	}
	if { [$wc.ani.animate get] == "" || [$wc.ani.animate get] == "0" } {
		$wc.ani.animate delete 0 end
		$wc.ani.animate insert 0 15
	}
	if { [$wc.ani.animate get] > 200 } {
		$wc.ani.animate delete 0 end
		$wc.ani.animate insert 0 200
	}
	if { $aniproc == 0} { 
		$wc.ani.animate_button configure -text "Stop"
		set aniproc 1
		$wc.ani.animate configure -state disabled
	} else {
		$wc.ani.animate_button configure -text "Animate"
		set aniproc 0
		$wc.ani.animate configure -state normal
	}
	if { $scale_to_move == "slices" } {
		set nslices1 [expr [Bimage get $theimg nslices] - 1]
		set slice [$wc.slice.scale get]
		set dir 1
		while { $aniproc } {
			after [expr 1000 / [$wc.ani.animate get]]
			if { $slice <= 0 } { set dir 1 }
			if { $slice >= $nslices1 } { set dir -1 }
			set slice [expr $slice + $dir]
			$wc.slice.scale set $slice
			update
		}
	} elseif { $scale_to_move == "images" } {
		set nimages1 [expr [Bimage get $theimg nimages] - 1]
		set currimage [$wc.image.scale get]
		set dir 1
		while { $aniproc } {
			after [expr 1000 / [$wc.ani.animate get]]
			if { $currimage <= 0 } { set dir 1 }
			if { $currimage >= $nimages1 } { set dir -1 }
			set currimage [expr $currimage + $dir]
			$wc.image.scale set $currimage
			update
		}
	}
}

