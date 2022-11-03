##
# @file		bshow_misc.tcl
#
# @brief	Miscellaneous procedures for Bshow and other scripts.
#
# @author	Bernard Heymann
# @date		Created: 20010516
# @date		Modified: 20220728


set PI 3.141592654
set HPI [expr $PI/2]

set px 1
set py 1
set mx 0
set my 0
set m0x 0
set m0y 0
set sx1 0
set sx2 0
set sy1 0
set sy2 0
set sz1 0
set sz2 0

set textfile ""
set cursorName "crosshair"
set trace_line 0

set rem_scale 1
set ori_color "blue"
set scalebar_color "white"
set scalebar_width 0
set scalebar_height 0
set scalebar_x -10
set scalebar_y -10
#set show_scalebar 1
set box_color "yellow"
set box_select_color "green"
set bad_color "red"
set bad_select_color "purple"
set filament_color "yellow"
set filament_first_color "blue"
set filament_select_color "green"
set spot_color "yellow"
set spot_select_color "green"
set ll_color "yellow"
set ll_select_color "green"
set marker_color "yellow"
set marker2_color "red"
set marker_select_color "green"
set error_color "orange"
set error_select_color "blue"

proc getControlWindow { theimg } {
	set wc ".wcontrol"
	if [string equal "theimg" $theimg] {
		set wc ""
	}
	return $wc
}

proc getImageCanvas { theimg } {
	set c ""
	if { [string length $theimg] < 1 } {
		puts "Error: image not selected!"
		return $c
	}
	set c ".w$theimg"
	if [string equal "theimg" $theimg] {
		set c ""
	}
	append c ".frame.c"
	return $c
}

## @brief Creates the right-side controls for image display
#
# @param	w 	 	 	Control window.

proc setupControls { w } {
	global theimg helv12
	global scalebar_width

	frame $w.img
	label $w.img.tag -text "Image type"
	tk_optionMenu $w.img.menu imgtype "" "mg" "frame" "part" "fil" "rec" "ft" "ps" "mod" "mask"
	$w.img.menu configure -width 10
#	trace add variable imgbigtype write imageType

	frame $w.mode
	checkbutton $w.mode.avg -text "Averaging mode" -variable mode \
	    -command { Update 1 }
	button $w.mode.invert -text "Invert" -command { invertContrast }

	frame $w.pixel
	label $w.pixel.xtag -font $helv12 -text "x"
	label $w.pixel.x -text "0" -relief sunken -bd 1 -width 5 \
		-font $helv12 -anchor w
	label $w.pixel.ytag -font $helv12 -text "y"
	label $w.pixel.y -text "0" -relief sunken -bd 1 -width 5 \
		-font $helv12 -anchor w
	label $w.pixel.vtag -font $helv12 -text "value"
	label $w.pixel.value -text "0" -relief sunken -bd 1 -width 30 \
		-font $helv12 -anchor w

	set pixel_size "1 1 1"
	setupVectorEntry $w.pixel_size "Pixel size" double $pixel_size "Image pixel size in Å/pixel"
	bind $w.pixel_size.x <Return> { setPixelsize $theimg }
	bind $w.pixel_size.y <Return> { setPixelsize $theimg }	
	bind $w.pixel_size.z <Return> { setPixelsize $theimg }
	
	frame $w.origin
	label $w.origin.tag -font $helv12 -text "Origin"
	label $w.origin.x -text "0" -relief sunken -bd 1 -width 6 \
		-font $helv12 -anchor w
	label $w.origin.y -text "0" -relief sunken -bd 1 -width 6 \
		-font $helv12 -anchor w
	label $w.origin.z -text "0" -relief sunken -bd 1 -width 6 \
		-font $helv12 -anchor w
	checkbutton $w.origin.show -text "Show" -variable show_origin -command { drawOrigin }
	$w.origin.show toggle

	frame $w.scalebar
	label $w.scalebar.tag -font $helv12 -text "Scale bar"
	entry $w.scalebar.width -validate key -invcmd bell -width 6 \
		-vcmd [list checkType double %P]
	$w.scalebar.width insert 0 $scalebar_width
	checkbutton $w.scalebar.show -text "Show" -variable show_scalebar -command { drawScaleBar }
	$w.scalebar.show toggle

	frame $w.distance
	label $w.distance.tag -font $helv12 -text "Distance" -width 10
	label $w.distance.value -text "" -relief sunken -bd 1 -width 16 \
		-font $helv12 -anchor w
	label $w.distance.unit -text "Å" -width 3 -font $helv12 -anchor w

	setupScales $w

	pack $w.img $w.mode -side top -pady 2
	pack $w.img.tag $w.img.menu -side left -pady 2 -anchor w
	pack $w.mode.avg $w.mode.invert -side left -pady 2 -anchor w

	pack $w.pixel.xtag $w.pixel.x $w.pixel.ytag $w.pixel.y $w.pixel.vtag \
		$w.pixel.value -side left -pady 10 -padx 1

	pack $w.origin.tag $w.origin.x $w.origin.y $w.origin.z $w.origin.show -side left -pady 6 -padx 2
	pack $w.scalebar.tag $w.scalebar.width $w.scalebar.show -side left -pady 6 -padx 2
	pack $w.distance.tag $w.distance.value $w.distance.unit -side left -pady 6

	pack $w.pixel $w.pixel_size $w.origin $w.scalebar $w.distance $w.image $w.slice $w.scale \
		$w.min $w.max $w.step $w.as $w.ani -side top -fill x -pady 2
}

## @brief Flips the display y axis.
#
# @param	y				Canvas y value.

proc yFlip { y } {
	global theimg
	set h [Bimage get $theimg height]
	return [expr $h - $y - 1]
}

## @brief Returns the image coordinates given canvas coordinates.
#
# @param	x
# @param	y				Canvas coordinates.

proc imageCoordinatesFromCanvas { x y } {
	global theimg
	global px py
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set z [$wc.slice.scale get]
	set px [expr floor( [$c canvasx $x] / $scale ) ]
	set py [expr floor( [$c canvasy $y] / $scale ) ]
	set pyi [yFlip $py]
	return "$px $pyi $z"
}

## @brief Selects a color.
#
# @param	w			parent window

proc selectColor { w } {
	global selected_color
    set initialColor $selected_color
    set color [tk_chooseColor -title "Choose a color" -parent $w \
		-initialcolor $selected_color]
    if { [string compare $color ""] } {
		set selected_color $color
		$w configure -background $color
    }
	Update 0
	return $selected_color
}

if {[tk windowingsystem] eq "aqua"} {
	proc ::tk::mac::ShowPreferences { } { setPreferences }
}

proc setPreferences { } {
	global helv12 rem_scale avg_check kernel
	global ori_color scalebar_color
	global box_color box_select_color bad_color bad_select_color
	global filament_color filament_first_color filament_select_color
	global ll_color ll_select_color
	global spot_color spot_select_color
	global marker_color marker_select_color error_color error_select_color

	set w .wprf

	if ![winfo exists .wprf] {
		toplevel $w

		checkbutton $w.rem_scale -text "Remember scale" -variable rem_scale
		
		setupCheckEntry $w.avg "Kernel" list kernel avg_check "The size of the filter kernel"

		labelframe $w.col -text "Colors"
		label $w.col.lbl -text "Draw  Select  First" -anchor e

		frame $w.col.origins
		label $w.col.origins.label -text "Origins" -width 15 -anchor e
		canvas $w.col.origins.well -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $ori_color

		frame $w.col.scalebar
		label $w.col.scalebar.label -text "Scale bar" -width 15 -anchor e
		canvas $w.col.scalebar.well -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $scalebar_color

		frame $w.col.box
		label $w.col.box.label -text "Boxes" -width 15 -anchor e
		canvas $w.col.box.well -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $box_color
		canvas $w.col.box.well2 -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $box_select_color

		frame $w.col.bad
		label $w.col.bad.label -text "Bad areas" -width 15 -anchor e
		canvas $w.col.bad.well -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $bad_color
		canvas $w.col.bad.well2 -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $bad_select_color

		frame $w.col.filament
		label $w.col.filament.label -text "Filaments" -width 15 -anchor e
		canvas $w.col.filament.well -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $filament_color
		canvas $w.col.filament.well2 -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $filament_select_color
		canvas $w.col.filament.well3 -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $filament_first_color

		frame $w.col.spot
		label $w.col.spot.label -text "Diffraction spots" -width 15 -anchor e
		canvas $w.col.spot.well -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $spot_color
		canvas $w.col.spot.well2 -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $spot_select_color

		frame $w.col.ll
		label $w.col.ll.label -text "Layer lines" -width 15 -anchor e
		canvas $w.col.ll.well -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $ll_color
		canvas $w.col.ll.well2 -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $ll_select_color

		frame $w.col.marker
		label $w.col.marker.label -text "Markers" -width 15 -anchor e
		canvas $w.col.marker.well -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $marker_color
		canvas $w.col.marker.well2 -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $marker_select_color

		frame $w.col.error
		label $w.col.error.label -text "Marker errors" -width 15 -anchor e
		canvas $w.col.error.well -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $error_color
		canvas $w.col.error.well2 -width 30 -height 30 -relief sunken \
				-borderwidth 2 -background $error_select_color

		button $w.update -text "Update" -relief raised -command { Update 0 }

		pack $w.col.origins.label $w.col.origins.well -side left
		pack $w.col.scalebar.label $w.col.scalebar.well -side left
		pack $w.col.box.label $w.col.box.well $w.col.box.well2 -side left
		pack $w.col.bad.label $w.col.bad.well $w.col.bad.well2 -side left
		pack $w.col.filament.label $w.col.filament.well $w.col.filament.well2 \
			$w.col.filament.well3 -side left
		pack $w.col.spot.label $w.col.spot.well $w.col.spot.well2 -side left
		pack $w.col.ll.label $w.col.ll.well $w.col.ll.well2 -side left
		pack $w.col.marker.label $w.col.marker.well $w.col.marker.well2 -side left
		pack $w.col.error.label $w.col.error.well $w.col.error.well2 -side left
		pack $w.col.lbl $w.col.origins $w.col.scalebar $w.col.box $w.col.bad $w.col.filament \
			$w.col.spot $w.col.ll $w.col.marker $w.col.error -side top -fill x -padx 5
		pack $w.rem_scale $w.avg $w.col $w.update -side top -pady 5 -padx 5
	} else {
		wm deiconify $w
		raise $w
    }
    wm title $w "Preferences"
    wm iconname $w "Preferences"
	bind $w.col.origins.well <1> { set ori_color [selectColor .wprf.col.origins.well] }
	bind $w.col.scalebar.well <1> { set scalebar_color [selectColor .wprf.col.scalebar.well] }
	bind $w.col.box.well <1> { set box_color [selectColor .wprf.col.box.well] }
	bind $w.col.box.well2 <1> { set box_select_color [selectColor .wprf.col.box.well2] }
	bind $w.col.bad.well <1> { set bad_color [selectColor .wprf.col.bad.well] }
	bind $w.col.bad.well2 <1> { set bad_select_color [selectColor .wprf.col.bad.well2] }
	bind $w.col.filament.well <1> { set filament_color [selectColor .wprf.col.filament.well] }
	bind $w.col.filament.well2 <1> { set filament_select_color \
		[selectColor .wprf.col.filament.well2] }
	bind $w.col.filament.well3 <1> { set filament_first_color \
		[selectColor .wprf.col.filament.well3] }
	bind $w.col.spot.well <1> { set spot_color [selectColor .wprf.col.spot.well] }
	bind $w.col.spot.well2 <1> { set spot_select_color [selectColor .wprf.col.spot.well2] }
	bind $w.col.ll.well <1> { set ll_color [selectColor .wprf.col.ll.well] }
	bind $w.col.ll.well2 <1> { set ll_select_color [selectColor .wprf.col.ll.well2] }
	bind $w.col.marker.well <1> { set marker_color [selectColor .wprf.col.marker.well] }
	bind $w.col.marker.well2 <1> { set marker_select_color [selectColor .wprf.col.marker.well2] }
	bind $w.col.error.well <1> { set error_color [selectColor .wprf.col.error.well] }
	bind $w.col.error.well2 <1> { set error_select_color [selectColor .wprf.col.error.well2] }
}

## @brief Saves the settings file .bshow
#

proc writeSettingsFile { } {
	global theimg project_root
	global window_width window_height window_scale
	global autoscale std_scale imageonly
	global rem_scale img_filter
	global lg avgtileps avgframeps tiled
	global tilex tiley tilez defocus
 	global tilt_axis tilt_angle tilt_offset show_tilt_axis
 	global direct_beam
	global volt Cs amp_fac
	global focal_length obj_apert slit
	global ctf_win_x ctf_win_y
	global ctf_scale_x ctf_scale_y
	global ctf_hires ctf_lores
    global def_start def_end def_inc
	global cont_update show_rings base_type env_type sub_base
	global box_size bad_radius
	global filament_width filament_node_radius boxing_interval helix_rise helix_angle
	global marker_radius
	global shift_limit thickness
	global spot_radius spot_resolution
	global ll_length ll_width ll_resolution
	global comp_radius link_radius selected_color line_width
	global gwidth gheight gorix goriy color
	global magnification magsize
	global ori_color 
	global scalebar_color scalebar_width scalebar_height scalebar_x scalebar_y
	global box_color box_select_color bad_color bad_select_color
	global filament_color filament_first_color filament_select_color
	global ll_color ll_select_color
	global spot_color spot_select_color
	global marker_color marker_select_color error_color error_select_color

	set c [getImageCanvas $theimg]
	set window_width [expr [winfo width $c] - 10]
	set window_height [expr [winfo height $c] - 10]
#	puts "Settings: $window_width $window_height $window_scale $autoscale"

	if { [catch {open .bshow w} fset] } { return }
	
	puts $fset "# .bshow: Bshow parameter file"
	puts $fset "# ----------------------------"
	puts $fset "# [clock format [clock seconds] -format {%Y-%m-%d %H:%M:%S}]"
	puts $fset "# This must be valid Tcl/Tk code."
	puts $fset ""
	puts $fset "# File access:"
	puts $fset "set img_filter $img_filter"
	puts $fset ""
	puts $fset "# Image window parameters:"
	puts $fset "set window_width $window_width"
	puts $fset "set window_height $window_height"
	puts $fset "set window_scale $window_scale"
	puts $fset "set rem_scale $rem_scale"
	puts $fset "set autoscale $autoscale"
	puts $fset "set std_scale $std_scale"
	puts $fset "set imageonly $imageonly"
	puts $fset "set ori_color $ori_color"
	puts $fset "set scalebar_color $scalebar_color"
	puts $fset "set scalebar_width $scalebar_width"
	puts $fset "set scalebar_height $scalebar_height"
	puts $fset "set scalebar_x $scalebar_x"
	puts $fset "set scalebar_y $scalebar_y"
	puts $fset ""
	puts $fset "# Micrograph parameters:"
	puts $fset "set project_root $project_root"
	puts $fset "set tilt_axis $tilt_axis"
	puts $fset "set tilt_angle $tilt_angle"
	puts $fset "set tilt_offset $tilt_offset"
	puts $fset "set show_tilt_axis $show_tilt_axis"
	puts $fset "set direct_beam $direct_beam"
	puts $fset ""
	puts $fset "# Powerspectrum parameters:"
	puts $fset "set lg $lg"
	puts $fset "set avgtileps $avgtileps"
	puts $fset "set avgframeps $avgframeps"
	puts $fset "set tiled $tiled"
	puts $fset "set tilex $tilex"
	puts $fset "set tiley $tiley"
	puts $fset "set tilez $tilez"
	puts $fset "set defocus $defocus"
	puts $fset ""
	puts $fset "# CTF parameters:"
	puts $fset "set volt $volt"
	puts $fset "set Cs $Cs"
	puts $fset "set amp_fac $amp_fac"
	puts $fset "set focal_length $focal_length"
	puts $fset "set obj_apert $obj_apert"
	puts $fset "set slit $slit"
	puts $fset "set ctf_hires $ctf_hires"
	puts $fset "set ctf_lores $ctf_lores"
	puts $fset "set ctf_win_x $ctf_win_x"
	puts $fset "set ctf_win_y $ctf_win_y"
	puts $fset "set ctf_scale_x $ctf_scale_x"
	puts $fset "set ctf_scale_y $ctf_scale_y"
	puts $fset "set def_start $def_start"
	puts $fset "set def_end $def_end"
	puts $fset "set def_inc $def_inc"
	puts $fset "set cont_update $cont_update"
	puts $fset "set show_rings $show_rings"
	puts $fset "set base_type $base_type"
	puts $fset "set env_type $env_type"
	puts $fset "set sub_base $sub_base"
	puts $fset ""
	puts $fset "# Boxing parameters:"
	puts $fset "set box_size \"$box_size\""
	puts $fset "set bad_radius $bad_radius"
	puts $fset "set box_color $box_color"
	puts $fset "set box_select_color $box_select_color"
	puts $fset "set bad_color $bad_color"
	puts $fset "set bad_select_color $bad_select_color"
	puts $fset ""
	puts $fset "# Filament parameters:"
	puts $fset "set filament_width $filament_width"
	puts $fset "set filament_node_radius $filament_node_radius"
	puts $fset "set boxing_interval $boxing_interval"
	puts $fset "set helix_rise $helix_rise"
	puts $fset "set helix_angle $helix_angle"
	puts $fset "set filament_color $filament_color"
	puts $fset "set filament_first_color $filament_first_color"
	puts $fset "set filament_select_color $filament_select_color"
	puts $fset ""
	puts $fset "# Tomography parameters:"
	puts $fset "set marker_radius $marker_radius"
	puts $fset "set marker_color $marker_color"
	puts $fset "set marker_select_color $marker_select_color"
	puts $fset "set error_color $error_color"
	puts $fset "set error_select_color $error_select_color"
	puts $fset "set shift_limit $shift_limit"
	puts $fset "set thickness $thickness"
	puts $fset ""
	puts $fset "# Crystallography parameters:"
	puts $fset "set spot_radius $spot_radius"
	puts $fset "set spot_resolution $spot_resolution"
	puts $fset "set spot_color $spot_color"
	puts $fset "set spot_select_color $spot_select_color"
	puts $fset ""
	puts $fset "# Helix parameters:"
	puts $fset "set ll_length $ll_length"
	puts $fset "set ll_width $ll_width"
	puts $fset "set ll_resolution $ll_resolution"
	puts $fset "set ll_color $ll_color"
	puts $fset "set ll_select_color $ll_select_color"
	puts $fset ""
	puts $fset "# Model parameters:"
	puts $fset "set comp_radius $comp_radius"
	puts $fset "set link_radius $link_radius"
	puts $fset "set selected_color $selected_color"
	puts $fset "set line_width $line_width"
	puts $fset ""
	puts $fset "# Graph parameters:"
	puts $fset "set gwidth $gwidth"
	puts $fset "set gheight $gheight"
	puts $fset "set gorix $gorix"
	puts $fset "set goriy $goriy"
	puts $fset "set color { blue green red yellow black SkyBlue2 }"
	puts $fset ""
	puts $fset "# Magnification parameters:"
	puts $fset "set magnification $magnification"
	puts $fset "set magsize $magsize"
	puts $fset ""
	close $fset
}

## @brief getFileName
# Edits a file name entry
#
# @param w 			the entry window

proc getFileName { w } {
	global imgdir
	set name [tk_getOpenFile -initialdir $imgdir]
	if { $name > " " } {
		$w.e delete 0 end
		$w.e insert 0 $name
	}
}

proc checkType {type val} {
	if { [expr {[string match {[.-+]} $val] || [string is $type $val]}] } {
		return true
	} else {
		return false
	}
}

proc setupEntry { w name type var desc } {
	frame $w
	label $w.l -text $name -width 20 -anchor w
	entry $w.e -validate key -invcmd bell -width 10 \
		-vcmd [list checkType $type %P]
	$w.e insert 0 $var
	button $w.b -text "?" -command "tk_messageBox -message \"$desc\" -type ok"
	pack $w.l $w.e $w.b -side left -padx 1m
}

proc setupResolutionEntry { w } {
	global hires lores
	labelframe $w.res -text "Resolution limits"
	setupEntry $w.res.hi "High" double $hires "High resolution limit in angstrom"
	setupEntry $w.res.lo "Low" double $lores "Low resolution limit in angstrom"
	pack $w.res.hi $w.res.lo -side top -fill x
}

proc setupCheckEntry { w name type var cvar desc } {
	frame $w
	checkbutton $w.c -text $name -variable $cvar
	entry $w.e -validate key -invcmd bell -width 10 \
		-vcmd [list checkType $type %P] -textvariable $var
#	$w.e insert 0 $var
	button $w.b -text "?" -command "tk_messageBox -message \"$desc\" -type ok"
	pack $w.c $w.e $w.b -side left -padx 1m
}

proc setupScrollEntry { w name var desc } {
	frame $w
	label $w.l -text $name -anchor w
	entry $w.e -validate key -invcmd bell -width 40 \
		 -xscrollcommand "$w.s set"
	$w.e insert 0 $var
	scrollbar $w.s -relief sunken -orient horiz -command "$w.e xview"
	button $w.b -text "Browse ..." -command "getFileName $w" -width 10
	button $w.i -text "?" -command "tk_messageBox -message \"$desc\" -type ok"
	pack $w.s $w.e -side bottom -fill x -expand 1 -padx 1m 
	pack $w.l $w.b $w.i -side left -fill x -padx 1m
}

proc setupVectorEntry { w name type var desc } {
	frame $w
	label $w.l -text $name -anchor w
	entry $w.x -width 6 -validate key -vcmd [list checkType $type %P]
	entry $w.y -width 6 -validate key -vcmd [list checkType $type %P]
	entry $w.z -width 6 -validate key -vcmd [list checkType $type %P]
	$w.x insert 0 [lindex $var 0]
	$w.y insert 0 [lindex $var 1]
	$w.z insert 0 [lindex $var 2]
	button $w.b -text "?" -command "tk_messageBox -message \"$desc\" -type ok"
	pack $w.l $w.x $w.y $w.z $w.b -side left -padx 1m
}

proc setupCheckSizeEntry { w name type var cvar desc } {
	frame $w
	checkbutton $w.c -text $name -variable $cvar -command [list updateCheckSizeEntry $w $var]
	entry $w.x -width 6 -validate key -vcmd [list checkType $type %P]
	entry $w.y -width 6 -validate key -vcmd [list checkType $type %P]
	entry $w.z -width 6 -validate key -vcmd [list checkType $type %P]
	$w.x insert 0 [lindex $var 0]
	$w.y insert 0 [lindex $var 1]
	$w.z insert 0 [lindex $var 2]
	button $w.b -text "?" -command "tk_messageBox -message \"$desc\" -type ok"
	pack $w.c $w.x $w.y $w.z $w.b -side left -padx 1m
}

proc updateCheckSizeEntry { w var } {
	set ax [$w.x get]
	set ay [$w.y get]
	set az [$w.z get]
	set kernel [list $ax $ay $az]
	return $kernel
}

proc setupFileEntry { w name file desc } {
	frame $w
	label $w.l -text $name -anchor w
	entry $w.e -width 40
	$w.e insert 0 $file
	button $w.i -text "?" -command "tk_messageBox -message \"$desc\" -type ok"
	pack $w.l $w.e $w.i -side left -pady 5
}

proc setupNewFileEntry { w name newfile desc } {
	frame $w
	label $w.l -text $name -anchor w
	entry $w.e -width 40
	$w.e insert 0 $newfile
	button $w.b -text "Load" \
		-relief raised -command " loadParam \[$w.e get] "
	button $w.i -text "?" -command "tk_messageBox -message \"$desc\" -type ok"
	pack $w.l $w.e $w.b $w.i -side left -pady 5
}

proc setupPostscriptEntry { w name psfile desc } {
	frame $w
	label $w.l -text $name -anchor w
	entry $w.e -width 40
	$w.e insert 0 $psfile
	button $w.b -text "Open" \
		-relief raised -command "openPostscript \[$w.e get] "
	button $w.i -text "?" -command "tk_messageBox -message \"$desc\" -type ok"
	pack $w.l $w.e $w.b $w.i -side left -pady 5
}

proc openPostscript { file } {
    global tcl_platform
    if { $tcl_platform(os) == "Darwin" } {
		exec open $file
	} elseif { $tcl_platform(os) == "Linux" } {
		exec evince $file
	}
}

## @brief Opens a text documentation file
#
# @param	helpfile			File name.

proc showHelp { helpfile } {
	global bshow_script
	if ![winfo exists .help] {
		toplevel .help
		frame .help.frame
		pack  .help.frame -expand yes -fill both -padx 1 -pady 1
		text .help.text -height 40 -wrap word -background white\
			-xscrollcommand ".help.xscroll set" \
			-yscrollcommand ".help.yscroll set" \
			-setgrid 1 -highlightthickness 0 -pady 2 -padx 3
		scrollbar .help.xscroll -command ".help.text xview" \
			-highlightthickness 0 -orient horizontal
		scrollbar .help.yscroll -command ".help.text yview" \
			-highlightthickness 0 -orient vertical

		grid .help.text -in .help.frame -padx 1 -pady 1 \
			-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
		grid .help.yscroll -in .help.frame -padx 1 -pady 1 \
			-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
		grid rowconfig    .help.frame 0 -weight 1 -minsize 0
		grid columnconfig .help.frame 0 -weight 1 -minsize 0
	} else {
		wm deiconify .help
		raise .help
    }
    wm title .help $helpfile
    wm iconname .help $helpfile
	readTextFile .help "$bshow_script/$helpfile"
}

## @brief Opens a URL in a web browser.
#
# @param	url				URL.

proc showURL { url } {
    global tcl_platform
	if [catch {
		switch $tcl_platform(os) {
			Darwin {
				exec open $url
			}
			Linux  -
			FreeBSD -
			OpenBSD -
			HP-UX -
			IRIX -
			AIX -
			SunOS {
				foreach pot_exe {firefox mozilla netscape iexplorer opera lynx
                       w3m links epiphany galeon konqueror mosaic amaya
                       browsex elinks} {
					set executable [auto_execok $pot_exe]
					break
				}
				if [string length $executable] {
					exec $executable $url &
				}
			}
			{Windows 95} -
			{Windows 98} -
			{Windows NT} {
				exec [auto_execok start] $url
			}
		}
	} err ] {
		tk_messageBox -icon error -message "Error: '$err' displaying '$url'"
	}
}

proc toolTip { w t } {
# 	bind $w <Any-Enter> "after 1000 showTip %X %Y $t; after 3000 hideTip"
#	bind $w <Any-Leave> "hideTip"
#	bind $w <ButtonPress-3> "set btn3 1; showTip %X %Y $t"
#	bind $w <ButtonRelease-3> "set btn3 0; hideTip"
	bind $w <ButtonPress-3> "showTip %X %Y $t"
	bind $w <ButtonRelease-3> "hideTip"
}

proc showTip { x y args } {
	wm geometry .tooltip +$x+$y
	.tooltip.message configure -text $args
	raise .tooltip
}

proc hideTip { } {
	lower .tooltip
	wm geometry .tooltip +0+0
	.tooltip.message configure -text ""
}

toplevel .tooltip
wm geometry .tooltip 1x1+0+0
wm overrideredirect .tooltip 1
label .tooltip.message -bd 1 -fg black -bg lightyellow -font fixed -text ""
pack .tooltip.message
hideTip

## @brief Updates the pixel size from parameter information
#

proc updatePixelsize { } {
	global filename project_item theimg imgtype

	set ps [Bimage get $theimg pixel_size]
	if [Bmg exists] {
		set n 1
#		puts "in updatePixelsize: Item: $project_item"
		if [regexp "Field" $project_item] {
			set n [Bmg get $project_item number_of_mg]
			if { $imgtype == "mg" || $imgtype == "frame" } {
				set ps [Bmg get $project_item pixel_size]
			} else {
				set ps [Bimage get $theimg pixel_size]
			}
		} elseif [regexp "Micrograph" $project_item] {
			if { $imgtype == "mg" || $imgtype == "frame" } {
				set ps [Bmg get $project_item pixel_size]
			} elseif { $imgtype == "part" } {
#				set n [Bmg get $project_item number_of_part]
				set ps [Bmg get $project_item part_pixel_size]
			} else {
				set ps [Bimage get $theimg pixel_size]
			}
		} elseif [regexp "Reconstruction" $project_item] {
			if { $imgtype == "rec" } {
				set ps [Bmg get $project_item pixel_size]
			} elseif { $imgtype == "part" } {
#				set n [Bmg get $project_item number_of_part]
				set ps [Bmg get $project_item part_pixel_size $filename]
			} else {
				set ps [Bimage get $theimg pixel_size]
			}
		} elseif [regexp "Class" $project_item] {
			if { $imgtype == "part" } {
				set ps [Bmg get $project_item part_pixel_size $filename]
			}
		}
#		puts "updatePixelsize: $ps"
		set i 0
		foreach {x y z} $ps {
			Bimage set $theimg pixel_size $i $x $y $z
			incr i
		}
	}
	updatePixelSizeEntry
}

## @brief Updates the current pixel size entry displayed
#

proc updatePixelSizeEntry { } {
	global theimg
	global scalebar_width scalebar_height
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set ps [Bimage get $theimg pixel_size $img_num]
#	puts "updatePixelSizeEntry: $ps"
	if [llength $ps] {
		$wc.pixel_size.x delete 0 end
		$wc.pixel_size.x insert 0 [lindex $ps 0]
		$wc.pixel_size.y delete 0 end
		$wc.pixel_size.y insert 0 [lindex $ps 1]
		$wc.pixel_size.z delete 0 end
		$wc.pixel_size.z insert 0 [lindex $ps 2]
		if { $scalebar_width < 2 } {
			set scalebar_width [expr 10*[lindex $ps 0]]
			set scalebar_height [expr $scalebar_width/5]
			$wc.scalebar.width config -text $scalebar_width
		}
	}
}

## @brief Returns the current pixel size entry
#

proc	getPixelSizeEntry { } {
	global project_item theimg
	set wc [getControlWindow $theimg]
	set psx [$wc.pixel_size.x get]
	set psy [$wc.pixel_size.y get]
	set psz [$wc.pixel_size.z get]
	set ps "$psx $psy $psz"
#	puts "getPixelSizeEntry: $ps"
	return $ps
}

## @brief Sets the pixel size of the image and micrograph
#
# theimg		Image identifier.

proc setPixelsize { theimg } {
	global project_item
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set ps [getPixelSizeEntry]
#	puts "setPixelsize: $ps"
	Bimage set $theimg pixel_size $img_num [lindex $ps 0] [lindex $ps 1] [lindex $ps 2]
	if [Bmg exists] {
		Bmg set $project_item pixel_size [lindex $ps 0] [lindex $ps 1] [lindex $ps 2]
		if [winfo exists .wctf] { updateCTF }
	}
	set ps [Bimage get $theimg pixel_size $img_num]
#	puts "setPixelsize: $ps"
}

## @brief Sets the pixel size of the image and micrograph
#
# theimg		Image identifier.

proc setAllPixelsizes { theimg } {
	global project_item
	set wc [getControlWindow $theimg]
	set ps [getPixelSizeEntry]
#	puts "setPixelsize: $ps"
	Bimage set $theimg pixel_size -1 [lindex $ps 0] [lindex $ps 1] [lindex $ps 2]
	if [Bmg exists] {
		Bmg set all pixel_size [lindex $ps 0] [lindex $ps 1] [lindex $ps 2]
		if [winfo exists .wctf] { updateCTF }
	}
	set ps [Bimage get $theimg pixel_size $img_num]
#	puts "setPixelsize: $ps"
}

## @brief Updates the image origins from parameters
#

proc updateOrigins { } {
	global theimg filename
	global imgtype project_item
	if [Bmg exists] {
		set n 1
#		puts "in updateOrigins: Item: $project_item"
		if [regexp "Field" $project_item] {
			set n [Bmg get $project_item number_of_mg]
			if { $imgtype == "mg" || $imgtype == "frame" } {
				set origin [Bmg get $project_item origin]
			} else {
				set origin [Bimage get $theimg origin]
			}
		} elseif [regexp "Micrograph" $project_item] {
			if { $imgtype == "mg" || $imgtype == "frame" } {
				set origin [Bmg get $project_item origin]
			} elseif { $imgtype == "part" } {
#				set n [Bmg get $project_item number_of_part]
				set origin [Bmg get $project_item part_origin $filename]
			} else {
				set origin [Bimage get $theimg origin]
			}
		} elseif [regexp "Reconstruction" $project_item] {
			if { $imgtype == "rec" } {
				set origin [Bmg get $project_item origin]
			} elseif { $imgtype == "part" } {
#				set n [Bmg get $project_item number_of_part]
				set origin [Bmg get $project_item part_origin $filename]
			} else {
				set origin [Bimage get $theimg origin]
			}
		} elseif [regexp "Class" $project_item] {
			if { $imgtype == "part" } {
				set origin [Bmg get $project_item part_origin $filename]
			}
		}
#		puts "Origin: $origin"
		set i 0
		foreach {x y z} $origin {
			Bimage set $theimg origin $x $y $z $i
			incr i
		}
	}
	drawOrigin
}

## @brief Draws the current image origin
#

proc drawOrigin { } {
	global theimg
	global show_origin ori_color
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	$c delete origin
	set img_num [$wc.image.scale get]
	set slice_num [$wc.slice.scale get]
	set scale [$wc.scale.scale get]
	set origin [Bimage get $theimg origin $img_num]
	set ox [lindex $origin 0]
	set oy [lindex $origin 1]
	set oz [lindex $origin 2]
	if { $show_origin } {
		if { [expr abs($oz - $slice_num)] < 1 } {
			set h [lindex [$c cget -scrollregion] 3]
			$c create oval [expr int(($ox - 4.5)*$scale)] \
				[expr int($h - ($oy - 4.5)*$scale)] \
				[expr int(($ox + 5.5)*$scale)] \
				[expr int($h - ($oy + 5.5)*$scale)] \
				-tags origin -outline $ori_color
		}
	}
	$wc.origin.x config -text $ox
	$wc.origin.y config -text $oy
	$wc.origin.z config -text $oz
}

proc getPixelCoordinates { c x y } {
	global theimg
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set h [lindex [$c cget -scrollregion] 3]
	set cx [$c canvasx $x]
	set cy [$c canvasy $y]
	set px [expr int($cx / $scale) ]
	set py [expr int(($h - $cy) / $scale) ]
	set p [list $px $py]
	return $p
}

## @brief Update or grabs items closeby whenever the mouse is pressed
#
# @param	c			Canvas object.
# @param	x
# @param	y			Coordinates on the canvas.

proc mousePressed { c x y } {
	global tool theimg
	global px py m0x m0y mx my
	set m0x $mx
	set m0y $my
	set mx $x
	set my $y
    set p [getPixelCoordinates $c $x $y ]
    set px [lindex $p 0]
    set py [lindex $p 1]
	showCoords
	updateMagnify $px $py default mouse
}

## @brief Updates or moves grabbed items whenever the mouse is moved.
#
# @param	c			Canvas object.
# @param	x
# @param	y			Coordinates on the canvas.
# @param	b			button pressed.

proc mouseMoved { c x y b } {
	global filename theimg
	global tool
	global px py
	global mx my
	global selectiontype trace_line
	if { [string length $filename] < 1 } { return }
	set wc [getControlWindow $theimg]
	set scale [$wc.scale.scale get]
	set cx [$c canvasx $x]
	set cy [$c canvasy $y]
	set cmx [$c canvasx $mx]
	set cmy [$c canvasy $my]
    set p1 [getPixelCoordinates $c $mx $my ]
    set p1x [lindex $p1 0]
    set p1y [lindex $p1 1]
    set p [getPixelCoordinates $c $x $y ]
    set px [lindex $p 0]
    set py [lindex $p 1]
	$c delete measure
	set w .wsel
	if { [winfo exists $w] && $b == 1 && $tool == "select" } {
		$c delete selection
		set xsize [expr abs($px - $p1x) + 1]
		set ysize [expr abs($py - $p1y) + 1]
		set zsize 1
		if { $selectiontype == "Square" || $selectiontype == "Cube" \
				|| $selectiontype == "Circle" || $selectiontype == "Sphere" } {
			if { $cy > $cmy } {
				set cy [expr $cmy + abs($cx - $cmx) ]
			} else {
				set cy [expr $cmy - abs($cx - $cmx) ]
			}
			set ysize $xsize
			if { $selectiontype == "Cube" || $selectiontype == "Sphere" } {
				set zsize $xsize
			}
		}
		if { $selectiontype == "Square" || $selectiontype == "Rectangle" || $selectiontype == "Cube" } {
			$c create rectangle $cmx $cmy $cx $cy -tags selection -outline red
		}
		if { $selectiontype == "Circle" || $selectiontype == "Ellipse" || $selectiontype == "Sphere" } {
			$c create oval $cmx $cmy $cx $cy -tags selection -outline red
		}
		$w.size.x config -text $xsize
		$w.size.y config -text $ysize
		$w.size.z config -text $zsize
	}
	if { $tool == "measure" || $trace_line } {
		$c create line $cmx $cmy $cx $cy -fill yellow -width 1 -tags measure
	}
	showCoords
	updateMagnify $px $py default mouse
}

## @brief Updates grabbed or nearby items whenever the mouse is released.
#
# @param	c			Canvas object.

proc mouseReleased { c } {
	global tool theimg
	global px py
	global mx my
	global m0x m0y
	global selectiontype calc_stats trace_line line_width
	global sx1 sx2 sy1 sy2 sz1 sz2
	set wc [getControlWindow $theimg]
	set slice_num [$wc.slice.scale get]
	set n [$wc.image.scale get]
	set scale [$wc.scale.scale get]
#	puts "m=$mx,$my p=$px,$py"
    set p0 [getPixelCoordinates $c $m0x $m0y ]
    set p1 [getPixelCoordinates $c $mx $my ]
	set sx0 [lindex $p0 0]
	set sx1 [lindex $p1 0]
	set sx2 $px
	set sy0 [lindex $p0 1]
	set sy1 [lindex $p1 1]
	set sy2 $py
	set sz0 $slice_num
	set sz1 $slice_num
	set sz2 $slice_num
	if { [winfo exists .wsel] && $tool == "select" } {
#		puts "Statistics for $mx $my $px $py"
		if { $selectiontype == "Square" || $selectiontype == "Cube" \
				|| $selectiontype == "Circle" || $selectiontype == "Sphere" } {
			if { $sy2 > $sy1 } {
				set sy2 [expr $sy1 + abs($sx2 - $sx1) ]
			} else {
				set sy2 [expr $sy1 - abs($sx2 - $sx1) ]
			}
			if { $selectiontype == "Cube" || $selectiontype == "Sphere" } {
				set dz [expr abs($sx1 - $sx2)]
				set zrad [expr $dz/2]
				set sz1 [expr $slice_num - $zrad]
				set sz2 [expr $slice_num + $zrad]
				if { [expr 2*$zrad] != $dz } { set sz1 [expr $sz1 - 1] }
#				puts "$sx1 $sx2 $zrad $sz1 $sz2"
			}
		}
		if { $calc_stats } { updateSelectionStats }
	} elseif { $tool == "measure" } {
		set ps [getPixelSizeEntry]
		set dx [expr ($sx2 - $sx0) * [lindex $ps 0]]
		set dy [expr ($sy2 - $sy0) * [lindex $ps 1]]
		set d [expr sqrt($dx*$dx + $dy*$dy)]
#		puts "Distance = $d angstrom"
#		puts [expr 57.296 * atan2($dy, $dx)]
		$wc.distance.value config -text "$d"
	} elseif { $tool == "voxels" } {
		set w .wvox
		if ![winfo exists $w] { Voxels }
		if { $trace_line } {
			set line_width [$w.type.width.e get]
			set line [Bimage get $theimg line $sx0 $sy0 $sz0 $sx2 $sy2 $sz2 $n $line_width]
			foreach {n x y z v} $line {
				$w.text insert end "$n $x $y $z $v\n"
			}
			set trace_line 0
		} else {
			set value [Bimage get $theimg pixel $sx2 $sy2 $sz2 $n]
#			puts "$sx2 $sy2 $value"
			$w.text insert end "$n $sx2 $sy2 $sz2 $value\n"
		}
		$w.text see end
	}
}

## @brief Displays updated coordinates whenever the mouse is moved.
#

proc showCoords { } {
	global filename
	global theimg imgtype
	global tool
	global px py
	if { [string length $filename] < 1 } { return }
	if ![Bimage exists $theimg] { return }
#	puts "px = $px  py = $py"
	set wc [getControlWindow $theimg]
	set img_num [$wc.image.scale get]
	set slice_num [$wc.slice.scale get]
	set value [Bimage get $theimg pixel $px $py $slice_num $img_num]
	if { $value eq " " } { return }
	$wc.pixel.x config -text "$px"
	$wc.pixel.y config -text "$py"
	$wc.pixel.value config -text "$value"
	if { $tool != "measure" } {
		set ps [getPixelSizeEntry]
#		puts "pixel size = $ps"
		if { $imgtype == "ft" || $imgtype == "ps" } {
			set width [Bimage get $theimg width]
			set height [Bimage get $theimg height]
			set depth [Bimage get $theimg nslices]
			set dx [expr ($px*1.0/$width - 0.5) / [lindex $ps 0]]
			set dy [expr ($py*1.0/$height - 0.5) / [lindex $ps 1]]
			if { $slice_num > 1 } {
				set dz [expr ($slice_num*1.0/$depth - 0.5) / [lindex $ps 2]]
			} else {
				set dz 0
			}
#			puts "$dx $dy $dz"
			set d [expr sqrt($dx*$dx + $dy*$dy + $dz*$dz)]
			if { $d < 0.0001 } { set d 0.0001 }
			set d [expr 1/$d]
			set u "1/A"
		} else {
			set origin [Bimage get $theimg origin]
			set dx [expr [lindex $ps 0]*($px - [lindex $origin 0])]
			set dy [expr [lindex $ps 1]*($py - [lindex $origin 1])]
			set dz [expr [lindex $ps 2]*($slice_num - [lindex $origin 2])]
			set d [expr sqrt($dx*$dx + $dy*$dy + $dz*$dz)]
			set u "A"
		}
		$wc.distance.value config -text "$d"
		$wc.distance.unit config -text "$u"
	}
}

## @brief Sets and updates the cutoff FOM.
#
# @param	cutoff		FOM cutoff value.

proc setFOMcutoff { cutoff } {
	global fom_cutoff
	set fom_cutoff $cutoff
	Update 0
}


## @brief Selects a cursor for the canvas window.
#

proc selectCursor { } {
	global theimg
    global tcl_platform
    global cursorName

	set c [getImageCanvas $theimg]

    set pad 4; # extra space around buttons
    set cols 6; # number of button columns

    set cursor_list {
		X_cursor arrow based_arrow_down based_arrow_up bottom_left_corner
		bottom_right_corner bottom_side bottom_tee box_spiral center_ptr circle
		cross cross_reverse crosshair diamond_cross dot dotbox double_arrow
		draft_large draft_small draped_box exchange fleur hand1 hand2
		icon iron_cross left_ptr left_side left_tee ll_angle lr_angle
		plus right_ptr right_side
		right_tee sb_down_arrow sb_h_double_arrow sb_left_arrow
		sb_right_arrow sb_up_arrow sb_v_double_arrow sizing
		target tcross top_left_arrow top_left_corner top_right_corner top_side
		top_tee ul_angle ur_angle
    }

    if { $tcl_platform(platform) == "windows" } {
        lappend cursor_list no starting size_ne_sw size_ns size_nw_se size_we uparrow wait
    }

    if { $tcl_platform(os) == "Darwin" } {
        lappend cursor_list text cross-hair
    }
	
	set w .wcurs
	toplevel $w
    wm title $w "Select a cursor"
    wm iconname $w "Cursors"

    foreach cursor $cursor_list {
        set cb [button $w.b_$cursor -text $cursor -cursor $cursor -command "configCursor_destroyWindow $c $w $cursor"]
        lappend cblist $cb
        if { [llength $cblist] >= $cols } {
            # place whole row of buttons
            eval grid $cblist -ipadx $pad -ipady $pad -sticky nswe
            set cblist {}
        }
    }
    if { [llength $cblist] > 0 } {
        # place rest of buttons
        eval grid $cblist -ipadx $pad -ipady $pad -sticky nswe
    }

    eval grid [button $w.cancel -text Cancel -relief flat -overrelief raised -command "destroy $w"] \
		-ipadx $pad -ipady $pad -sticky nswe

}

proc configCursor_destroyWindow { c w name } {
    configCursor $c $name
	destroy $w
}

proc configCursor { c name } {
    global cursorName tcl_platform
    set cursorName $name
#	$c config -cursor $cursorName
	switch -- $tcl_platform(platform) {
		unix {
#			$c config -cursor "$cursorName #ffffff #ffff00"
			$c config -cursor $cursorName
		}
		windows {
			$c config -cursor $cursorName
		}
	}
}

if {[tk windowingsystem] eq "aqua"} {
	proc ::tk::mac::ShowHelp { } {
		global Bsoft
		showURL "file://$Bsoft/doc/bshow/bshow.html"
	}
}

proc canvas2Photo { canvId } {
    # The following line grabs the contents of the canvas canvId into photo image ph.
    set retVal [catch {image create photo -format window -data $canvId} ph]
    if { $retVal != 0 } {
       puts "\n\tFATAL ERROR: Cannot create photo from canvas window"
       exit 1
    }
    return $ph
}

proc capture {W format file} {
	if {[catch {package require Img}]} {
		# Error thrown - package not found.
		# Set up a dummy interface to work around the absence
		puts "package Img not found!"
		return
	}
#	raise $W
#	update idletasks
	set image [image create photo -format window -data $W]
	$image write -format $format $file
	puts "capture -> '$file' ([file size $file] bytes)"
	image delete $image
 }

## @brief Draws a scale bar
#

proc drawScaleBar { } {
	global theimg
	global show_scalebar scalebar_color
	global scalebar_width scalebar_height
	global scalebar_x scalebar_y
	set c [getImageCanvas $theimg]
	set wc [getControlWindow $theimg]
	$c delete scalebar
	if { $show_scalebar } {
		set scale [$wc.scale.scale get]
		set ps [$wc.pixel_size.x get]
		set scalebar_width [$wc.scalebar.width get]
		if { $scalebar_width < 2 } {
			set scalebar_width [expr 10*$ps]
			set scalebar_height [expr $scalebar_width/5]
			$wc.scalebar.width config -text $scalebar_width
		}
		set w [expr $scale*$scalebar_width/$ps]
		set h [expr $scale*$scalebar_height/$ps]
		set sr [$c cget -scrollregion]
		set x1 [expr [lindex $sr 2] + $scalebar_x]
		set y1 [expr [lindex $sr 3] + $scalebar_y]
		set x [expr $x1 - $w]
		set y [expr $y1 - $h]
		$c create rect $x $y $x1 $y1 \
			-outline $scalebar_color -fill $scalebar_color -tags scalebar
	}
}

## @brief Sets scale bar parameters.
#

proc scaleBar { } {
	global theimg
	global helv12
	global show_scalebar scalebar_color
	global scalebar_width scalebar_height
	global scalebar_x scalebar_y
	
	set w .sb
	catch {destroy $w}
	toplevel $w
	wm title $w "Set scale bar"
	wm iconname $w "Scale bar"

#	label $w.msg -font $helv12 -wraplength 4i -justify left \
#		-text {Selecting OK sets the scale bar parameters}
#	pack $w.msg -side top -padx 5 -pady 5

	setupEntry $w.width "Width" double $scalebar_width "The width of the scale bar in angstroms"
	setupEntry $w.height "Height" double $scalebar_height "The height of the scale bar in angstroms"
	setupEntry $w.x "X offset" integer $scalebar_x "The offset of the scale bar from the right edge"
	setupEntry $w.y "Y offset" integer $scalebar_y "The offset of the scale bar from the bottom edge"

	frame $w.buttons
	button $w.buttons.ok -text Set -command "setScaleBar $w"
	button $w.buttons.cancel -text Cancel -command "destroy $w"
	
	pack $w.buttons.ok $w.buttons.cancel -side left -expand 1
	pack $w.width $w.height $w.x $w.y $w.buttons -side top -fill x -pady 2m -padx 2m
	
	bind $w <Return> "setScaleBar $w"
}

## @brief Sets scale bar parameters.
#

proc setScaleBar { w } {
	global theimg
	global show_scalebar scalebar_color
	global scalebar_width scalebar_height
	global scalebar_x scalebar_y
	set wc [getControlWindow $theimg]
	set scalebar_width [$w.width.e get]
	set scalebar_height [$w.height.e get]
	set scalebar_x [$w.x.e get]
	set scalebar_y [$w.y.e get]
	set show_scalebar 1
	$wc.scalebar.width config -text $scalebar_width
	drawScaleBar
	destroy $w
}

