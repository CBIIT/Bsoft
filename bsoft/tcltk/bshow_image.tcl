##
# @file		bshow_image.tcl
#
# @brief	Procedures for manipulating images in Bshow.
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20191001

set show_mont 1
set show_mont_labels 1
set show_mont_select 1
set	avg_check 0
set kernel 1
set lg 0
set avgtileps 0
set avgframeps 0
set tiled 0
set tilex 512
set tiley 512
set tilez 512
set defocus 2.0

## @brief Menu for manipulating image files
#

proc setupImageMenu {} {
	.menuBar add cascade -menu .menuBar.image -label "Image" -underline 0
	menu .menuBar.image -tearoff 0
	.menuBar.image add command -label "Set scales" -underline 0 \
		-command { setScaleDialog } -accelerator "Ctrl-z"
	.menuBar.image add command -label "Switch slices and images" -underline 0 \
		-command { switchImage }
	.menuBar.image add command -label "Modify image" -underline 0 \
		-command { modifyImage }

	menu .menuBar.image.dt -tearoff 0
	.menuBar.image add cascade -menu .menuBar.image.dt -label "Change data type" -underline 0
	.menuBar.image.dt add command -label "Unsigned char (byte)" -underline 0 \
		-command { changeDataType byte } -accelerator "Ctrl-z"
	.menuBar.image.dt add command -label "Signed char" -underline 0 \
		-command { changeDataType char } -accelerator "Ctrl-z"
	.menuBar.image.dt add command -label "Unsigned short" -underline 0 \
		-command { changeDataType "unsigned short" } -accelerator "Ctrl-z"
	.menuBar.image.dt add command -label "Signed short" -underline 0 \
		-command { changeDataType short } -accelerator "Ctrl-z"
	.menuBar.image.dt add command -label "Unsigned int" -underline 0 \
		-command { changeDataType "unsigned int" } -accelerator "Ctrl-z"
	.menuBar.image.dt add command -label "Signed int" -underline 0 \
		-command { changeDataType int } -accelerator "Ctrl-z"
	.menuBar.image.dt add command -label "Unsigned long" -underline 0 \
		-command { changeDataType "unsigned long" } -accelerator "Ctrl-z"
	.menuBar.image.dt add command -label "Signed long" -underline 0 \
		-command { changeDataType long } -accelerator "Ctrl-z"
	.menuBar.image.dt add command -label "Float" -underline 0 \
		-command { changeDataType float } -accelerator "Ctrl-z"
	.menuBar.image.dt add command -label "Double" -underline 0 \
		-command { changeDataType double } -accelerator "Ctrl-z"

	.menuBar.image add command -label "Fix type sign" -underline 0 \
		-command { fixType }
	.menuBar.image add command -label "Set origin" -underline 0 \
		-command { setOrigin }
	.menuBar.image add command -label "Center origin" -underline 0 \
		-command { centerOrigin }
	.menuBar.image add command -label "Shift origin to zero" -underline 0 \
		-command { zeroOrigin }
	.menuBar.image add command -label "Scale bar" -underline 0 \
		-command { scaleBar }
	.menuBar.image add command -label "Crop" -underline 0 \
		-command { cropImage }
	.menuBar.image add command -label "Pad" -underline 0 \
		-command { Pad }
	.menuBar.image add command -label "Bin" -underline 0 \
		-command { Bin }
	.menuBar.image add command -label "Montage" -underline 0 \
		-command { Montage }
	.menuBar.image add command -label "Invert" -underline 0 \
		-command { invertImage }
	.menuBar.image add command -label "Histogram" -underline 0 \
		-command { Histogram }
	.menuBar.image add command -label "Filter" -underline 0 \
		-command { Filter }
	.menuBar.image add command -label "Fourier transform" -underline 0 \
		-command { FourierTransform -1 }
	.menuBar.image add command -label "Fourier back transform" -underline 0 \
		-command { FourierTransform 1 }
	.menuBar.image add command -label "Power spectrum" -underline 0 \
		-command { PowerSpectrum }
	.menuBar.image add command -label "Diffraction" -underline 0 \
		-command { Diffraction }
		

}

## @brief Switches the images and slices.
#

proc switchImage { } {
	global theimg
	Bimage switch $theimg
	reconfigureScales
	Update 1
}

## @brief Sets the image type.
#

proc setImageType { } {
	global imgtype
	set w .wtyp
	if ![winfo exists $w] {
		toplevel $w
		
		radiobutton $w.notdef -text "Not defined" -variable imgtype -value ""
		radiobutton $w.mg -text "Micrograph(s)" -variable imgtype -value "mg"
		radiobutton $w.frame -text "Frame(s)" -variable imgtype -value "frame"
		radiobutton $w.part -text "Particle(s)" -variable imgtype -value "part"
		radiobutton $w.fil -text "Filament(s)" -variable imgtype -value "fil"
		radiobutton $w.ft -text "Fourier transform(s)" -variable imgtype -value "ft"
		radiobutton $w.ps -text "Power spectrum(s)" -variable imgtype -value "ps"
		radiobutton $w.rec -text "Reconstruction(s)" -variable imgtype -value "rec"
		radiobutton $w.rec -text "Mask(s)" -variable imgtype -value "mask"
		
		frame $w.buttons
		button $w.buttons.ok -text OK -command "destroy $w"

		pack $w.notdef $w.mg $w.frame $w.part $w.fil $w.ft \
			$w.ps $w.rec $w.buttons -side top -fill x -pady 5		
		pack $w.buttons.ok -side left -expand 1
	
		bind $w <Return> "destroy $w"
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Image type"
    wm iconname $w "Image type"
	tkwait window $w
}


## @brief Switchs the sign convention for the data type.
#

proc fixType { } {
	global theimg
	Bimage fixtype $theimg
	reconfigureScales
	Update 1
}

## @brief Sets the image origin.
#

proc setOrigin { } {
	global theimg
	global helv12
	
	set ori [Bimage get $theimg origin]
	
	set w .ori
	catch {destroy $w}
	toplevel $w
	wm title $w "Set origin"
	wm iconname $w "Origin"

	label $w.msg -font $helv12 -wraplength 4i -justify left \
		-text {Selecting OK sets the image origin}
	pack $w.msg -side top -padx 5 -pady 5

	setupVectorEntry $w.coor "Coordinates" double $ori "The image origin coordinates"
	pack $w.coor -side top -fill x -pady 2m

	frame $w.buttons
	button $w.buttons.center -text Center -command "setImageOrigin $w 1"
	button $w.buttons.ok -text OK -command "setImageOrigin $w 0"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
	
	pack $w.buttons -side bottom -fill x -pady 2m
	pack $w.buttons.center $w.buttons.ok $w.buttons.cancel -side left -expand 1
	
	bind $w <Return> "setImageOrigin $w 0"
}

## @brief Sets the image origin.
#
# @param	w				window.
# @param	t				0=coordinates, 1=center.

proc setImageOrigin { w t } {
	global theimg
	if { $t } {
		set ox [expr [Bimage get $theimg width] / 2]
		set oy [expr [Bimage get $theimg height] / 2]
		set oz [expr [Bimage get $theimg nslices] / 2]
	} else {
		set ox [$w.coor.x get]
		set oy [$w.coor.y get]
		set oz [$w.coor.z get]
	}
	Bimage set $theimg origin $ox $oy $oz
	destroy $w
	reconfigureScales
	Update 1
}

## @brief Gets the x origin in canvas coordinates.
#
# @param	c				Canvas object.

proc getOriginX { c } {
	global theimg
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set origin [Bimage get $theimg origin]
	return [expr $scale * [lindex $origin 0]]
}

## @brief Gets the y origin in canvas coordinates.
#
# @param	c				Canvas object.

proc getOriginY { c } {
	global theimg
	set wc [getControlWindow $theimg]
	set h [lindex [$c cget -scrollregion] 3]
	set scale [$wc.scale.scale get]
	set origin [Bimage get $theimg origin]
	return [expr $h - $scale * ([lindex $origin 1] + 1)]
}

## @brief Centers the image origin with wrapping.
#

proc centerOrigin { } {
	global theimg
	global project_item
	if ![Bimage exists $theimg] { return }
	Bimage center $theimg
	set origin [Bimage get $theimg origin]
	if [Bmg exists] {
		Bmg set $project_item origin [lindex $origin 0] [lindex $origin 1] [lindex $origin 2]
	}
	Update 1
}

## @brief Shifts the image origin to zero  {0,0,0} with wrapping.
#

proc zeroOrigin { } {
	global theimg
	global project_item
	if ![Bimage exists $theimg] { return }
	Bimage zero_origin $theimg
	set origin [Bimage get $theimg origin]
	if [Bmg exists] {
		Bmg set $project_item origin [lindex $origin 0] [lindex $origin 1] [lindex $origin 2]
	}
	Update 1
}

## @brief Extracts part of an image.
#

proc extractImage { } {
	global theimg
	global filename
	global imgdir
	global filetypes
	global tool
	global sx1 sx2 sy1 sy2 sz1 sz2
	if ![Bimage exists $theimg] { return }
	set wc [getControlWindow $theimg]
	set x [Bimage get $theimg width]
	set y [Bimage get $theimg height]
	set z [$wc.slice.scale get]
	set vlist "0 0 $z $x $y $z"
	if { [winfo exists .wsel] && $tool == "select" } {
		if { $sx1 != $sx2 || $sy1 != $sy2 } {
			set vlist "$sx1 $sy1 $sz1 $sx2 $sy2 $sz2"
		}
	}
#	puts $vlist
	set img_num [$wc.image.scale get]
	set currdir $imgdir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
	set extfile [tk_getSaveFile -filetypes $filetypes \
		-initialfile [file tail $filename] -initialdir $currdir -defaultextension .mrc]
	if { [string length $extfile] < 1 } { return }
	Bimage extract $theimg $extfile $img_num $vlist
}

## @brief Crops an image to a smaller size.
#

proc cropImage { } {
	global theimg
	global tool
	global sx1 sx2 sy1 sy2 sz1 sz2
	if ![Bimage exists $theimg] { return }
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set x [Bimage get $theimg width]
	set y [Bimage get $theimg height]
	set z [$wc.slice.scale get]
	set vlist "0 0 $z $x $y $z"
	if { [winfo exists .wsel] && $tool == "select" } {
		set vlist "$sx1 $sy1 $sz1 $sx2 $sy2 $sz2"
		set sx1 0
		set sx2 0
		set sy1 0
		set sy2 0
		set sz1 0
		set sz2 0
	}
#	puts $vlist
	Bimage crop $theimg $img_num $vlist
	reconfigureScales
	Update 1
}

## @brief Pads an image to a larger size.
#

proc Pad { } {
	global theimg
	global helv12
	global datatype
	set bkg [Bimage get $theimg background]
	set w .pad
	catch {destroy $w}
	toplevel $w
	wm title $w "Pad image"
	wm iconname $w "Pad"

	label $w.msg -font $helv12 -wraplength 4i -justify left \
		-text {Setting these options modify the internal image on selecting OK}
	pack $w.msg -side top -padx 5 -pady 5

#	set filltype 2

	frame $w.fill
	label $w.fill.tag -font $helv12 -text "Fill type" -width 12
	radiobutton $w.fill.avg -text "Average" -variable filltype -value 1
	radiobutton $w.fill.back -text "Background" -variable filltype -value 2
	radiobutton $w.fill.val -text "Value" -variable filltype -value 0
	entry $w.fill.entry -width 6 -validate key -vcmd { string is double %P }
	$w.fill.entry insert 0 $bkg
	$w.fill.back select
	pack $w.fill.tag $w.fill.avg $w.fill.back $w.fill.val \
		$w.fill.entry -side left -fill both -expand true
	
	set start "0 0 0"
	set end "0 0 0"

	setupVectorEntry $w.start "Translation" int $start "Translation in pixels (can be negative)"
	setupVectorEntry $w.end "Size" int $end "Final size in pixels"
	pack $w.start $w.end $w.fill -side top -padx 5 -pady 5
	
	frame $w.buttons
	button $w.buttons.ok -text OK -command "padImage $w"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
	
	pack $w.buttons -side bottom -fill x -pady 2m
	pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
	
	bind $w <Return> "padImage $w"
}

## @brief Pads the image.
#
# @param	w 			the window to destroy.

proc padImage { w } {
	global theimg
	global filltype
	set xm [$w.start.x get]
	set xp [$w.end.x get]
	set ym [$w.start.y get]
	set yp [$w.end.y get]
	set zm [$w.start.z get]
	set zp [$w.end.z get]
	set fill [$w.fill.entry get]
	destroy $w
	Bimage pad $theimg $xm $ym $zm $xp $yp $zp $filltype $fill
	reconfigureScales
	Update 1
}

## @brief Invokes a dialog box with modifyable image parameters.
#

proc modifyImage { } {
	global theimg
	global helv12
	global datatype
	set datatype [Bimage get $theimg datatype]
	set origin [Bimage get $theimg origin]
	set nimages [Bimage get $theimg nimages]
	set nslices [Bimage get $theimg nslices]
	set w .mod
	catch {destroy $w}
	toplevel $w
	wm title $w "Modify image"
	wm iconname $w "Modify"

	label $w.msg -font $helv12 -wraplength 4i -justify left \
		-text {Setting these options modify the internal image on selecting OK}

	frame $w.panel
	frame $w.panel.datatype
	label $w.panel.datatype.tag -text "Data type"
	tk_optionMenu $w.panel.datatype.menu datatype "unknown" "unsigned char (byte)" \
		"signed char" "unsigned short" "short" "unsigned int" "int" \
		"unsigned long" "long" "float" "double"
	
	frame $w.panel.origin
	label $w.panel.origin.tag -font $helv12 -text "Origin" 
	entry $w.panel.origin.entryx -width 6 -validate key -vcmd { string is double %P }
	entry $w.panel.origin.entryy -width 6 -validate key -vcmd { string is double %P }
	entry $w.panel.origin.entryz -width 6 -validate key -vcmd { string is double %P }
	$w.panel.origin.entryx insert 0 [lindex $origin 0]
	$w.panel.origin.entryy insert 0 [lindex $origin 1]
	$w.panel.origin.entryz insert 0 [lindex $origin 2]
	
	frame $w.panel.reslice
	label $w.panel.reslice.tag -font $helv12 -text "Reslice" 
	entry $w.panel.reslice.entry -width 10
	$w.panel.reslice.entry insert 0 "xyz"

	frame $w.panel.rotate
	label $w.panel.rotate.tagpsi -font $helv12 -text "Rotate: psi" 
	label $w.panel.rotate.tagtheta -font $helv12 -text "theta" 
	label $w.panel.rotate.tagphi -font $helv12 -text "phi" 
	entry $w.panel.rotate.entrypsi -width 4 -validate key -vcmd { string is double %P }
	entry $w.panel.rotate.entrytheta -width 4 -validate key -vcmd { string is double %P }
	entry $w.panel.rotate.entryphi -width 4 -validate key -vcmd { string is double %P }
	$w.panel.rotate.entrypsi insert 0 "0"
	$w.panel.rotate.entrytheta insert 0 "0"
	$w.panel.rotate.entryphi insert 0 "0"

	frame $w.panel.truncate
	label $w.panel.truncate.tagmin -font $helv12 -text "Truncate: min" 
	label $w.panel.truncate.tagmax -font $helv12 -text "max" 
	entry $w.panel.truncate.entrymin -width 6 -validate key -vcmd { string is double %P }
	entry $w.panel.truncate.entrymax -width 6 -validate key -vcmd { string is double %P }
	$w.panel.truncate.entrymin insert 0 [Bimage get $theimg min]
	$w.panel.truncate.entrymax insert 0 [Bimage get $theimg max]

	frame $w.panel.rescale
	label $w.panel.rescale.tagavg -font $helv12 -text "Rescale: avg" 
	label $w.panel.rescale.tagstd -font $helv12 -text "std" 
	entry $w.panel.rescale.entryavg -width 6 -validate key -vcmd { string is double %P }
	entry $w.panel.rescale.entrystd -width 6 -validate key -vcmd { string is double %P }
	$w.panel.rescale.entryavg insert 0 [Bimage get $theimg average]
	$w.panel.rescale.entrystd insert 0 [Bimage get $theimg standard_deviation]

	frame $w.buttons
	button $w.buttons.ok -text OK -command "modifyImage_and_DestroyWindow $w"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
	
	pack $w.msg -side top -padx 5 -pady 5
	pack $w.panel -side top -fill both
	pack $w.panel.datatype $w.panel.origin $w.panel.reslice $w.panel.rotate \
		 $w.panel.truncate $w.panel.rescale -side top -fill both
	
	pack $w.panel.datatype.tag  $w.panel.datatype.menu \
			-side left -padx 5 -pady 5
	pack $w.panel.origin.tag $w.panel.origin.entryx $w.panel.origin.entryy \
			$w.panel.origin.entryz -side left -padx 5 -pady 5
	pack $w.panel.reslice.tag $w.panel.reslice.entry -side left -padx 5 -pady 5
	pack $w.panel.rotate.tagpsi $w.panel.rotate.entrypsi $w.panel.rotate.tagtheta \
		$w.panel.rotate.entrytheta $w.panel.rotate.tagphi $w.panel.rotate.entryphi \
		-side left -padx 5 -pady 5
	pack $w.panel.truncate.tagmin $w.panel.truncate.entrymin \
		$w.panel.truncate.tagmax $w.panel.truncate.entrymax -side left -padx 5 -pady 5
	pack $w.panel.rescale.tagavg $w.panel.rescale.entryavg \
		$w.panel.rescale.tagstd $w.panel.rescale.entrystd -side left -padx 5 -pady 5
	pack $w.buttons -side bottom -fill x -pady 2m
	pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
	
	bind $w <Return> "modifyImage_and_DestroyWindow $w"
}

## @brief Modifies the image.
#
# @param	w 			the window to destroy.

proc modifyImage_and_DestroyWindow { w } {
	global theimg
	global datatype project_item
	set do_trunc 0
	set do_rescale 0
	set ox [$w.panel.origin.entryx get]
	set oy [$w.panel.origin.entryy get]
	set oz [$w.panel.origin.entryz get]
	set order [$w.panel.reslice.entry get]
	set psi [$w.panel.rotate.entrypsi get]
	set theta [$w.panel.rotate.entrytheta get]
	set phi [$w.panel.rotate.entryphi get]
	set min [$w.panel.truncate.entrymin get]
	set max [$w.panel.truncate.entrymax get]
	if { [expr abs($min - [Bimage get $theimg min])] > 0.001 \
			|| [expr abs($max - [Bimage get $theimg max])] > 0.001 } {
		set do_trunc 1
	}
	set avg [$w.panel.rescale.entryavg get]
	set std [$w.panel.rescale.entrystd get]
	if { [expr abs($avg - [Bimage get $theimg average])] > 0.001 \
			|| [expr abs($std - [Bimage get $theimg standard_deviation])] > 0.001 } {
		set do_rescale 1
	}
	destroy $w
	Bimage set $theimg origin $ox $oy $oz
	if [Bmg exists] { Bmg set $project_item origin $ox $oy $oz }
	Bimage reslice $theimg $order
	Bimage rotate $theimg $psi $theta $phi
	if { $do_trunc } {
		Bimage truncate $theimg $min $max
	}
	if { $do_rescale } {
		Bimage rescale $theimg $avg $std
	}
	set olddatatype [Bimage get $theimg datatype]
	if { $datatype != $olddatatype } {
		Bimage set $theimg datatype $datatype
	}
	Update 1
	reconfigureScales
}

## @brief This chnages the data type.
#

proc changeDataType { datatype } {
	global theimg
	set olddatatype [Bimage get $theimg datatype]
	if { $datatype != $olddatatype } {
		Bimage set $theimg datatype $datatype
	}
	Update 1
	reconfigureScales
}

## @brief This inverts the image.
#

proc invertImage {} {
	global theimg
	Bimage invert $theimg
	Update 1
	reconfigureScales
}

## @brief This invokes a dialog box with binning parameters.
#

proc Bin {} {
	global theimg
	global helv12

	set w .wbin
	catch {destroy $w}
	toplevel $w
	wm title $w "Bin image"
	wm iconname $w "Bin"

	set bin "1 1 1"

	label $w.msg -justify left \
		-text {Setting these options modify the internal image on selecting OK}
	pack $w.msg -side top -padx 5 -pady 5

	setupVectorEntry $w.bin "Binning" int $bin "The level of binning"

	frame $w.buttons
	button $w.buttons.ok -text OK -command "binImage $w"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
	pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
	
	pack $w.msg $w.bin $w.buttons -side top -pady 2m
	
	bind $w <Return> "binImage $w"
}

## @brief This inverts the image.
#

proc binImage { w } {
	global theimg
	set bx [$w.bin.x get]
	set by [$w.bin.y get]
	set bz [$w.bin.z get]
	destroy $w
	Bimage bin $theimg $bx $by $bz
	Update 1
	updatePixelsize
	reconfigureScales
}

## @brief This invokes a dialog box with filter parameters.
#

proc Filter {} {
	global theimg filter_type
	global helv12 kernel avg_check

	set w .wfilt
	catch {destroy $w}
	toplevel $w
	wm title $w "Filter image"
	wm iconname $w "Filter"

	label $w.msg -justify left \
		-text {Setting these options modify the internal image on selecting OK}
	pack $w.msg -side top -padx 5 -pady 5

	frame $w.type   
	label $w.type.l -font $helv12 -text "Type" -bg orange
	tk_optionMenu $w.type.opt filter_type "Averaging" "Variance" "Normalizing" "Gaussian"
	pack $w.type.l $w.type.opt -side left -ipadx 2

	setupEntry $w.kernel "Kernel" integer $kernel \
		"The size of the filter kernel\nThe Gaussian sigma is 1/6 of the kernel size"

	frame $w.buttons
	button $w.buttons.pref -text "Set preference" -command { set avg_check 1 }
	button $w.buttons.ok -text Filter -command "filterImage $w"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
	pack $w.buttons.pref $w.buttons.ok $w.buttons.cancel -side left -expand 1
	
	pack $w.msg $w.type $w.kernel $w.buttons -side top -pady 2m
	
	bind $w <Return> "filterImage $w"
}

## @brief This filters the image.
#

proc filterImage { w } {
	global theimg filter_type kernel
	set kernel [$w.kernel.e get]
	destroy $w
	if { [string compare $filter_type "Averaging"] == 0 } {
		Bimage average $theimg $kernel $kernel $kernel
	} elseif { [string compare $filter_type "Variance"] == 0 } {
		Bimage variance $theimg $kernel $kernel $kernel
	} elseif { [string compare $filter_type "Normalizing"] == 0 } {
		Bimage normalize $theimg $kernel $kernel $kernel
	} elseif { [string compare $filter_type "Gaussian"] == 0 } {
		set sigma [expr $kernel/6.0]
#		puts "Gaussian filter: $sigma"
		Bimage gaussian $theimg $kernel $sigma
	}
	Update 1
	reconfigureScales
}

## @brief This invokes a dialog box with montage parameters.
#

proc Montage {} {
	global theimg
	global helv12 font_size selected_color
	global flipx flipy
	set bkg [Bimage get $theimg background]
	set w .wmont
	catch {destroy $w}
	toplevel $w
	wm title $w "Montage image"
	wm iconname $w "Montage"

	label $w.msg -font $helv12 -wraplength 4i -justify left \
		-text {Setting these options modify the internal image on selecting OK}
	pack $w.msg -side top -padx 5 -pady 5

	setupEntry $w.col "Columns" integer 0 "Number of columns"
	setupEntry $w.row "Rows" integer 0 "Number of rows"
	setupEntry $w.first "First panel" integer 0 "First image or slice to start from"
	setupEntry $w.skip "Skip" integer 0 "Number of images or slices to skip between panels"
	setupEntry $w.pad "Padding" integer 0 "Amount of padding between panels"
	setupEntry $w.fill "Fill value" double $bkg "User-defined fill value between panels"
	frame $w.flip
	checkbutton $w.flip.x -text "Flip X" -variable flipx
	checkbutton $w.flip.y -text "Y" -variable flipy
	pack $w.flip.x $w.flip.y -side left
	pack $w.col $w.row $w.first $w.skip $w.pad \
		 $w.fill $w.flip -side top -padx 5 -pady 5

	frame $w.switches
	label $w.switches.tag -text "Show:"
	checkbutton $w.switches.boxes -text "Boxes" -variable show_mont \
				-command { drawMontageOverlay }
	checkbutton $w.switches.labels -text "Labels" -variable show_mont_labels \
				-command { drawMontageOverlay }
	checkbutton $w.switches.select -text "Selection" -variable show_mont_select \
				-command { drawMontageOverlay }
	pack $w.switches.tag $w.switches.boxes $w.switches.labels $w.switches.select \
			-side left -ipadx 5 -ipady 5 -padx 1m

	frame $w.font
	label $w.font.tag -font $helv12 -text "Font size:"
	tk_optionMenu $w.font.size font_size "6" "9" "12" "18" "24" "32" "40"
	canvas $w.font.color -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $selected_color
	pack $w.font.tag $w.font.size $w.font.color -side left -ipadx 5 -ipady 5 -padx 1m

	frame $w.buttons
	button $w.buttons.ok -text Montage -command "montageImage $w"
	button $w.buttons.cancel -text Unset -command "unsetMontage"	
	button $w.buttons.save -text Save -command "saveMontage"
	button $w.buttons.del -text "Delete selected" -command "deleteMontageSelection"
	pack $w.buttons.ok $w.buttons.cancel $w.buttons.save $w.buttons.del -side left -expand 1

	pack $w.switches $w.font $w.buttons -side top -fill x -pady 2m
	
	bind $w.font.color <1> { selectColor .wmont.font.color }
	bind $w <Return> "montageImage $w"
}

## @brief Montages the image.
#

proc montageImage { w } {
	global theimg
	global flipx flipy montover font_size selected_color
	set first [$w.first.e get]
	set ncol [$w.col.e get]
	set nrow [$w.row.e get]
	set skip [$w.skip.e get]
	set pad [$w.pad.e get]
	set fill [$w.fill.e get]
	set flip [expr $flipx + 2 * $flipy]
	set nx [Bimage get $theimg width]
	set ny [Bimage get $theimg height]
	Bimage set $theimg montage $ncol $nrow $first $skip $pad $fill $flip
	set montover(active) 1
	set montover(nx) $nx
	set montover(ny) $ny
	set montover(ncol) $ncol
	set montover(nrow) $nrow
	set montover(first) $first
	set montover(skip) $skip
	set montover(pad) $pad
	set montover(flip) $flip
	set montover(fontsize) $font_size
	set montover(fontcolor) $selected_color
	reconfigureScales
	Update 1
}

proc unsetMontage { } {
	global theimg montover
	set c [getImageCanvas $theimg]
	$c delete mont
	set montover(active) 0
	Bimage set $theimg montage 0 0
	reconfigureScales
	Update 1
}

proc drawMontageOverlay { } {
	global theimg
	global flipx flipy montover show_mont show_mont_labels show_mont_select
	set c [getImageCanvas $theimg]
	$c delete mont
	if { $montover(active) < 1 } { return }
	if { $montover(nx) < 1 } { return }
	if { $montover(ny) < 1 } { return }
#	puts Overlay
	set nimg [Bimage get $theimg nimages]
	set scale [.scale.scale get]
	set nx [expr $scale * $montover(nx)]
	set ny [expr $scale * $montover(ny)]
	set first $montover(first)
	set skip $montover(skip)
	set ncol $montover(ncol)
	set nrow $montover(nrow)
	set pad [expr $scale * $montover(pad)]
	if { $flipy } {
		set ymin [expr $pad + $ny]
	} else {
		set ymin [expr $pad * (2 * $nrow - 1) + $ny * $nrow]
	}
	set i $first
	set show_mont_font [list "Helvetica" $montover(fontsize) "bold"]
#	puts $show_mont_font
	for { set row 0 } { $row < $nrow } { incr row 1 } {
		set ymax [expr $ymin - $ny + 1]
		if { $flipx } {
			set xmin [expr $nx * ($ncol - 1) + $pad * (2 * $ncol - 1)]
		} else {
			set xmin $pad
		}
		for { set col 0 } { $col < $ncol } { incr col 1 } {
			set xmax [expr $xmin + $nx - 1]
#			puts "$xmin $ymin $xmax $ymax"
			set xt [expr $xmin + 5]
			set yt [expr $ymin - 3]
			set sel [array get montover "$col,$row" ]
			if [llength $sel] {
				set color red
			} else {
				set color yellow
			}
			if $show_mont {
				$c create rect $xmin $ymin $xmax $ymax -outline $color \
					-width 1 -tags mont
			}
			if { $i < $nimg } {
				if $show_mont_labels {
					$c create text $xt $yt -text $i -anchor sw -fill $color \
						-font $show_mont_font -tags mont
				}
				set xt [expr $xmax - 5]
				if $show_mont_select {
					set img_sel [Bimage get $theimg select $i]
					$c create text $xt $yt -text $img_sel -anchor se \
						-fill $montover(fontcolor) -font $show_mont_font -tags mont
				}
			}
			if { $flipx } {
				set xmin [expr $xmin - $nx - 2 * $pad]
			} else {
				set xmin [expr $xmin + $nx + 2 * $pad]
			}
			set i [expr $i + $skip + 1]
		}
		if { $flipy } {
			set ymin [expr $ymin + $ny + 2 * $pad]
		} else {
			set ymin [expr $ymin - $ny - 2 * $pad]
		}
	}
	bind $c <1> "montageSelect %x %y"
	bind . <Key-o> { saveMontageSelection }
}

proc montageSelect { x y } {
	global flipx flipy montover
	set scale [.scale.scale get]
	set nx [expr $scale * $montover(nx)]
	set ny [expr $scale * $montover(ny)]
	set ncol $montover(ncol)
	set nrow $montover(nrow)
	set pad [expr $scale * $montover(pad)]
	if { $flipx } {
		set col [expr $ncol - int($x/($nx + 2 * $pad)) - 1]
	} else {
		set col [expr int($x/($nx + 2 * $pad))]
	}
	if { $flipy } {
		set row [expr int($y/($ny + 2 * $pad))]
	} else {
		set row [expr $nrow - int($y/($ny + 2 * $pad)) - 1]
	}
	set sel [array get montover "$col,$row" ]
	if [llength $sel] {
		unset montover($col,$row)
	} else {
		set montover($col,$row) 1
	}
#	puts "$col $row"
#	foreach key [array names montover "*,*"] {
#		puts "$key $montover($key)"
 #	}
	drawMontageOverlay
}

proc saveMontage { } {
	global theimg montover
	global imgdir filename
	global filetypes
	set currdir $imgdir
	if { [string length $currdir] < 2 } { set currdir [pwd] }
	set new_filename [tk_getSaveFile -filetypes $filetypes \
			-initialfile [file tail $filename] -initialdir $currdir \
			-defaultextension .mrc]
	if { [string length $new_filename] < 1 } { return }
	set new_filename [relativePath [pwd] $new_filename]
	Bimage save_montage $theimg $new_filename
}

proc saveMontageSelection { } {
	global montover
	set w .wtxt
	set wt [openTextWindow $w]
	set first $montover(first)
	set skip $montover(skip)
	set ncol $montover(ncol)
	foreach key [array names montover "*,*"] {
#		puts "$key $montover($key)"
		set coor [split $key ,]
		set col [lindex $coor 0]
		set row [lindex $coor 1]
		set i [expr $first + $ncol * $row + $col]
#		$wt insert end "$col\t$row\t$i\n"
		$wt insert end "$i\n"
 	}
	saveAsTextFile $wt
}

proc deleteMontageSelection { } {
	global montover theimg
	set first $montover(first)
	set skip $montover(skip)
	set ncol $montover(ncol)
	set delist ""
	foreach key [array names montover "*,*"] {
#		puts "$key $montover($key)"
		set coor [split $key ,]
		set col [lindex $coor 0]
		set row [lindex $coor 1]
		set i [expr $first + $ncol * $row + $col]
		if [string length $delist] {
			append delist ",$i"
		} else {
			append delist $i
		}
		unset montover($key)
 	}
#	puts $delist
	Bimage delete $theimg $delist
	reconfigureScales
	Update 1
}

## @brief Calculates a Fourier transform of the image data.
#

proc FourierTransform { dir } {
	global theimg
	global imgtype
	Bimage fft $theimg $dir
	if { $dir < 0 } {
		set imgtype "ft"
	} else {
		set imgtype ""
	}
	reconfigureScales
	.scale.scale set 1.0
}

## @brief Calculates a power spectrum of the image data.
#

proc PowerSpectrum { } {
	global theimg project_item
	global helv12
	global lg avgtileps avgframeps tiled
	global tilex tiley tilez
	global tilt_axis tilt_angle tilt_offset defocus

	set nx [Bimage get $theimg width]
	set ny [Bimage get $theimg height]
	set nz [Bimage get $theimg nslices]
	if { $tilex > $nx } { set tilex $nx }
	if { $tiley > $ny } { set tiley $ny }
	if { $tilez > $nz } { set tilez $nz }
	
	set tilesize "$tilex $tiley $tilez"
	set Cs 2.0
	set volt 300

	if [Bmg exists] {
		set defocus [ expr 1e-4 * [Bmg get $project_item defocus] ]
		set tilt_axis [Bmg get $project_item axis]
		set tilt_angle [Bmg get $project_item tilt]
		set Cs [expr 1e-7 * [Bmg get $project_item Cs] ]
		set volt [expr 1e-3 * [Bmg get $project_item volt] ]
	}
	
#	puts "$Cs $volt"
	
	set w .wps
	catch {destroy $w}
	toplevel $w
	wm title $w "Power spectrum"
	wm iconname $w "PS"

	label $w.msg -font $helv12 -wraplength 4i -justify left \
		-text {Select parameters to calculate a power spectrum}
	pack $w.msg -side top -padx 5 -pady 5

	checkbutton $w.logbutton -text "Logarithm" -variable lg	
	checkbutton $w.avgtiles -text "Average tile power spectra" -variable avgtileps
	checkbutton $w.avgframes -text "Average frame power spectra" -variable avgframeps
	
	setupCheckSizeEntry $w.tiles "Tiles" int $tilesize tiled "Size of tiles to divide the micrograph"

	labelframe $w.tilt -text "Tilt parameters"
	setupEntry $w.tilt.axis "Tilt axis angle" double $tilt_axis "Micrograph tilt axis angle conter-clockwise from the x-axis in degrees"
	setupEntry $w.tilt.angle "Tilt angle" double $tilt_angle "Micrograph tilt angle in degrees"
	setupEntry $w.tilt.offset "Offset from origin" double $tilt_offset "Tilt axis offset from the micrograph origin"
	setupEntry $w.tilt.defocus "Defocus" double $defocus "Estimate/guess of defocus"
	setupEntry $w.tilt.cs "Spherical aberration" double $Cs "Cs in mm"
	setupEntry $w.tilt.volt "Acceleration voltage" double $volt "Acceleration voltage in kV"
	pack $w.tilt.axis $w.tilt.angle $w.tilt.offset $w.tilt.defocus \
		$w.tilt.cs $w.tilt.volt -side top -expand 1
	
	frame $w.buttons
	button $w.buttons.ok -text OK -command "calcPowerSpec $w"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
	pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
		
	pack $w.logbutton $w.avgframes $w.avgtiles $w.tiles $w.tilt $w.buttons -side top

	bind $w <Return> "calcPowerSpec $w"
}

## @brief Calculates a power spectrum of the image data.
#
# @param	w 			the window to destroy.

proc calcPowerSpec {w} {
	global theimg
	global lg avgtileps avgframeps tiled
	global tilex tiley tilez
	global tilt_axis tilt_angle tilt_offset defocus
	global imgtype
	global filename project_item
	global window_width window_scale

	set tilex [$w.tiles.x get]
	set tiley [$w.tiles.y get]
	set tilez [$w.tiles.z get]
	set tilt_axis [$w.tilt.axis.e get]
	set tilt_angle [$w.tilt.angle.e get]
	set tilt_offset [$w.tilt.offset.e get]
	set defocus [$w.tilt.defocus.e get]
	set cs [$w.tilt.cs.e get]
	set volt [$w.tilt.volt.e get]
	if { $tilt_angle != 0 } { set tiled 1 }
	destroy $w
#	puts "Powerspectrum: $lg $tiled"
	set flags [expr 1 + 2*$avgtileps + 4 + 8*$lg + 16*$avgframeps]
	if $tiled {
		Bimage powerspec $theimg $flags $tilex $tiley $tilez $tilt_axis $tilt_angle $tilt_offset $defocus $cs $volt
	} else {
		Bimage powerspec $theimg $flags
	}
	set imgtype "ps"
	if [Bmg exists] {
		setMgParam
		Bmg set $project_item filename $filename ps
		Bmg set $project_item axis [expr $tilt_axis]
		Bmg set $project_item tilt [expr $tilt_angle]
	}
	reconfigureRangeStep
	reconfigureScales
	set scale 1.0
	set width [Bimage get $theimg width]
	if { $width > $window_width } { set scale [expr 0.1*((10*$window_width)/$width) ] }
	if { $scale < 0.1 } { set scale 0.1 }
	set window_scale $scale
	.scale.scale set $scale
	set ps [Bimage get $theimg pixel_size 0]
#	puts "Power spectrum pixel size: $ps"
}

